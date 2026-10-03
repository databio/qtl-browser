"""Refget anchor: name the reference sequences a release sits on, and check the alleles against them.

`refget_store` opens (or builds) a refgetstore and writes `_tables/reference.json`: the seqcol
collection digest and each chromosome's sequence digest. `variants_refcheck` then looks up every
allele-bearing variant's reference bases and writes `_checks/refcheck/summary.json`. The manifest
copies both into its `reference` block. None of this touches the pack format.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .common import ROOT, Config, die, log, variants_path, write_parquet

# match classes, in the order `classify` codes them
MATCH = ("a1", "a2", "both", "neither", "unchecked")
A1, A2, BOTH, NEITHER, UNCHECKED = range(len(MATCH))
_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def reference_json(cfg: Config) -> Path:
    return cfg.tables / "reference.json"


def refcheck_summary(cfg: Config) -> Path:
    return cfg.derived / "_checks" / "refcheck" / "summary.json"


def refcheck_table(cfg: Config, chrom: str) -> Path:
    """Per-chromosome classification, one row per variant. The summary counts these; this is where
    you go to look at an individual call."""
    return cfg.tables / "refcheck" / f"{chrom}.parquet"


def open_store(cfg: Config, build: bool = False):
    """The refgetstore named by config: remote when `store_url` is set, else the local one. With
    `build`, a missing local store is created from the configured FASTA."""
    from gtars.refget import RefgetStore
    ref = cfg["reference"]
    store_path = ROOT / ref["store"]
    if ref.get("store_url"):
        store = RefgetStore.open_remote(str(store_path), ref["store_url"])
    elif (store_path / "rgstore.json").exists():
        store = RefgetStore.open_local(str(store_path))
    elif build:
        fasta = cfg.raw / ref["fasta"]
        if not fasta.exists():
            die(f"refget_store: {fasta} is missing; download the reference FASTA first")
        log(f"refget_store: building {store_path} from {fasta.name} (digests every sequence once)")
        store = RefgetStore.on_disk(str(store_path))
        store.add_sequence_collection_from_fasta(str(fasta))
    else:
        die(f"no refgetstore at {store_path}; run `build --step refget_store`")
    store.set_quiet(True)
    return store


def _collection_digest(cfg: Config, store) -> str:
    want = cfg["reference"].get("collection")
    if want:
        if store.get_collection_metadata(want) is None:
            die(f"refget_store: collection {want} (config reference.collection) is not in the store")
        return want
    found = [m.digest for m in store.list_collections(page_size=1000)["results"]]
    if len(found) != 1:
        die(f"refget_store: the store holds {len(found)} collections; set reference.collection in config.yaml")
    log(f"refget_store: collection {found[0]}; paste it into config.yaml as reference.collection")
    return found[0]


def refget_store(cfg: Config) -> None:
    ref = cfg["reference"]
    store = open_store(cfg, build=True)
    digest = _collection_digest(cfg, store)
    meta = store.get_collection_metadata(digest)
    # get_collection returns stubs: names, lengths and digests without sequence bytes
    by_name = {r.metadata.name: r.metadata for r in store.get_collection(digest)}
    missing = [c for c in ref["chromosomes"] if c not in by_name]
    if missing:
        die(f"refget_store: collection {digest} has no sequence named {missing}")
    out = {
        "assembly": ref["assembly"], "fasta": ref["fasta"], "store_url": ref.get("store_url"),
        "refget": {"collection": digest, "names_digest": meta.names_digest,
                   "lengths_digest": meta.lengths_digest, "sequences_digest": meta.sequences_digest},
        "sequences": {c: {"digest": by_name[c].sha512t24u, "length": by_name[c].length, "md5": by_name[c].md5}
                      for c in ref["chromosomes"]},
    }
    reference_json(cfg).parent.mkdir(parents=True, exist_ok=True)
    reference_json(cfg).write_text(json.dumps(out, indent=2))
    log(f"refget_store: {len(out['sequences'])} sequences -> {reference_json(cfg).relative_to(cfg.derived)}")


def _matches(seq: np.ndarray, p: int, allele: str) -> bool:
    """VCF-style: an allele matches when the reference reads the same from `p` (1-based) on."""
    a = allele.upper().encode()
    return p - 1 + len(a) <= len(seq) and seq[p - 1:p - 1 + len(a)].tobytes() == a


def classify(seq: np.ndarray, pos: np.ndarray, a1: list, a2: list) -> np.ndarray:
    """Which allele the reference carries, coded as indexes into MATCH. `seq` is the upper-case
    chromosome as uint8; `pos` is 1-based; a null allele makes the variant `unchecked`. SNPs compare
    one base, vectorized. Other alleles compare by prefix, so a deletion's shorter allele matches
    too; when both match with different lengths, the longer one is the reference."""
    n = len(pos)
    out = np.full(n, UNCHECKED, dtype=np.int8)
    has = np.array([x is not None and y is not None for x, y in zip(a1, a2)], dtype=bool)
    snp = has & np.array([x is not None and y is not None and len(x) == 1 and len(y) == 1
                          for x, y in zip(a1, a2)], dtype=bool)
    inside = (pos >= 1) & (pos <= len(seq))
    idx = np.flatnonzero(snp & inside)
    if len(idx):
        ref = seq[pos[idx] - 1]
        b1 = np.frombuffer("".join(a1[i] for i in idx).upper().encode(), dtype=np.uint8)
        b2 = np.frombuffer("".join(a2[i] for i in idx).upper().encode(), dtype=np.uint8)
        m1, m2 = ref == b1, ref == b2
        out[idx] = np.where(m1 & m2, BOTH, np.where(m1, A1, np.where(m2, A2, NEITHER)))
    out[snp & ~inside] = NEITHER
    for i in np.flatnonzero(has & ~snp):
        p, x, y = int(pos[i]), a1[i], a2[i]
        m1, m2 = _matches(seq, p, x), _matches(seq, p, y)
        if m1 and m2 and len(x) != len(y):
            m1, m2 = len(x) > len(y), len(y) > len(x)
        out[i] = BOTH if m1 and m2 else A1 if m1 else A2 if m2 else NEITHER
    return out


def _revcomp(a: str) -> str:
    return a.upper().encode().translate(_COMP)[::-1].decode()


def neither_detail(seq: np.ndarray, pos: np.ndarray, a1: list, a2: list, codes: np.ndarray) -> dict:
    """For `neither` variants: how many match after reverse complement (a strand flip) and how
    many sit on an N in the reference."""
    flip = on_n = 0
    for i in np.flatnonzero(codes == NEITHER):
        p = int(pos[i])
        if 1 <= p <= len(seq) and seq[p - 1] == ord("N"):
            on_n += 1
        elif _matches(seq, p, _revcomp(a1[i])) or _matches(seq, p, _revcomp(a2[i])):
            flip += 1
    return {"strand_flip": flip, "on_n": on_n}


def load_sequence(store, digest: str, length: int) -> np.ndarray:
    """One whole chromosome as upper-case uint8 (chr1 is 249 MB; load one at a time)."""
    return np.frombuffer(store.get_substring(digest, 0, length).upper().encode("ascii"), dtype=np.uint8)


def match_fraction(tally: dict) -> float:
    n = sum(tally[k] for k in ("a1", "a2", "both", "neither"))
    return (tally["a1"] + tally["a2"] + tally["both"]) / n if n else 0.0


def check_threshold(cis: dict, minimum: float) -> float:
    frac = match_fraction(cis)
    if frac < minimum:
        raise ValueError(f"variants_refcheck: cis alleles match the reference for {frac:.5f} of variants, "
                         f"below reference.allele_check.min_match_fraction {minimum}")
    return frac


def variants_refcheck(cfg: Config) -> None:
    if not reference_json(cfg).exists():
        die("variants_refcheck: _tables/reference.json is missing; run `build --step refget_store`")
    ref = json.loads(reference_json(cfg).read_text())
    store = open_store(cfg)
    table = pq.read_table(variants_path(cfg), columns=["chr", "position", "A1", "A2", "in_cis"])
    blank = {k: 0 for k in MATCH} | {"strand_flip": 0, "on_n": 0}
    total = {"cis": dict(blank), "trans_only": dict(blank)}
    per_chrom, examples = {}, {}
    for chrom in sorted(set(table["chr"].to_pylist())):
        if chrom not in ref["sequences"]:
            die(f"variants_refcheck: variants on {chrom}, which reference.chromosomes does not name")
        t = table.filter(pc.equal(table["chr"], chrom))
        s = ref["sequences"][chrom]
        seq = load_sequence(store, s["digest"], s["length"])
        pos = t["position"].to_numpy().astype(np.int64)
        a1, a2 = t["A1"].to_pylist(), t["A2"].to_pylist()
        codes = classify(seq, pos, a1, a2)
        in_cis = t["in_cis"].to_numpy(zero_copy_only=False)
        per_chrom[chrom] = {}
        for label, mask in (("cis", in_cis), ("trans_only", ~in_cis)):
            tally = {k: int(np.sum(codes[mask] == j)) for j, k in enumerate(MATCH)}
            tally |= neither_detail(seq, pos[mask], [x for x, m in zip(a1, mask) if m],
                                    [y for y, m in zip(a2, mask) if m], codes[mask])
            per_chrom[chrom][label] = tally
            for k, v in tally.items():
                total[label][k] += v
        bad = np.flatnonzero(codes == NEITHER)[:20]
        examples[chrom] = [{"position": int(pos[i]), "A1": a1[i], "A2": a2[i], "in_cis": bool(in_cis[i]),
                            "ref": chr(seq[pos[i] - 1]) if 1 <= pos[i] <= len(seq) else None} for i in bad]
        # one row per variant, so a single call can be looked at without re-running the step
        write_parquet(pa.table({
            "chr": pa.array([chrom] * len(pos), pa.string()),
            "position": pa.array(pos, pa.int32()),
            "A1": pa.array(a1, pa.string()), "A2": pa.array(a2, pa.string()),
            "in_cis": pa.array(in_cis, pa.bool_()),
            "match": pa.DictionaryArray.from_arrays(pa.array(codes, pa.int32()), pa.array(MATCH, pa.string())),
        }), refcheck_table(cfg, chrom), 200_000, stats_columns=["chr", "position"])
        c = per_chrom[chrom]["cis"]
        log(f"variants_refcheck: {chrom:5s} cis a1={c['a1']:,} a2={c['a2']:,} both={c['both']:,} neither={c['neither']:,}")
        del seq
    minimum = float(cfg["reference"]["allele_check"]["min_match_fraction"])
    out = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "collection": ref["refget"]["collection"],
           "min_match_fraction": minimum, "cis_match_fraction": match_fraction(total["cis"]),
           "total": total, "per_chrom": per_chrom, "neither_examples": examples}
    refcheck_summary(cfg).parent.mkdir(parents=True, exist_ok=True)
    refcheck_summary(cfg).write_text(json.dumps(out, indent=2))
    log(f"variants_refcheck: per-variant calls in {refcheck_table(cfg, '<chr>').parent.relative_to(cfg.derived)}/")
    c = total["cis"]
    log(f"variants_refcheck: A1 is the reference for {c['a1']:,} cis variants, A2 for {c['a2']:,}; "
        f"both {c['both']:,}, neither {c['neither']:,} (strand flip {c['strand_flip']:,}, on N {c['on_n']:,})")
    frac = check_threshold(c, minimum)
    log(f"variants_refcheck: cis match fraction {frac:.5f} (minimum {minimum})")


def manifest_block(cfg: Config) -> dict:
    """The manifest's `reference` block. Data about the release, not part of the pack format."""
    for path, step in ((reference_json(cfg), "refget_store"), (refcheck_summary(cfg), "variants_refcheck")):
        if not path.exists():
            raise FileNotFoundError(f"manifest: {path.relative_to(cfg.derived)} is missing; run `build --step {step}`")
    ref = json.loads(reference_json(cfg).read_text())
    rc = json.loads(refcheck_summary(cfg).read_text())
    if rc["collection"] != ref["refget"]["collection"]:
        raise ValueError("manifest: refcheck ran against another collection; re-run `build --step variants_refcheck --force`")
    cis, tr = rc["total"]["cis"], rc["total"]["trans_only"]
    return {
        "assembly": ref["assembly"],
        "refget": {"collection": ref["refget"]["collection"], "store_url": ref["store_url"]},
        "sequences": {c: {"digest": s["digest"], "length": s["length"]} for c, s in ref["sequences"].items()},
        "allele_check": {
            "cis_with_alleles": sum(cis[k] for k in ("a1", "a2", "both", "neither")),
            "a1_is_reference": cis["a1"], "a2_is_reference": cis["a2"], "both": cis["both"], "neither": cis["neither"],
            "trans_only_with_alleles": sum(tr[k] for k in ("a1", "a2", "both", "neither")),
            "trans_only_neither": tr["neither"],
            "min_match_fraction": rc["min_match_fraction"],
            "source": f"{refcheck_summary(cfg).relative_to(cfg.derived)} (generated {rc['generated']})",
        },
    }


def validate(cfg: Config, man: dict, check) -> None:
    """The manifest's reference block agrees with the store and the refcheck summary, and names
    every chromosome a pack does."""
    import re
    ref = man.get("reference")
    if not ref or not ref.get("refget", {}).get("collection") or not ref.get("sequences"):
        check(False, "manifest.reference has refget.collection and sequences (run `build --step manifest`)")
        return
    seqs = ref["sequences"]
    digest_ok = re.compile(r"^[A-Za-z0-9_-]{32}$")
    bad = [c for c, s in seqs.items() if not digest_ok.match(s.get("digest", ""))
           or not isinstance(s.get("length"), int) or s["length"] <= 0]
    bad += [] if digest_ok.match(ref["refget"]["collection"]) else ["collection"]
    check(not bad, f"manifest.reference digests are 32-character sha512t24u and lengths positive ({bad[:3]})")

    files = (man.get("packs") or {}).get("files") or {}
    named = {c for kind in ("variants", "eqtl", "sqtl", "trans", "gwas", "hits") for c in files.get(kind, {})}
    check(named <= set(seqs), f"every pack chromosome is in manifest.reference.sequences ({sorted(named - set(seqs))})")

    store = open_store(cfg)
    coll = ref["refget"]["collection"]
    if store.get_collection_metadata(coll) is None:
        check(False, f"the refgetstore holds collection {coll}")
    else:
        by_name = {r.metadata.name: r.metadata for r in store.get_collection(coll)}
        wrong = [c for c, s in seqs.items() if c not in by_name
                 or by_name[c].sha512t24u != s["digest"] or by_name[c].length != s["length"]]
        check(not wrong, f"manifest.reference sequences match the store's collection {coll} ({wrong[:3]})")

    rc_path = refcheck_summary(cfg)
    if not rc_path.exists():
        check(False, "_checks/refcheck/summary.json exists (run `build --step variants_refcheck`)")
        return
    rc = json.loads(rc_path.read_text())
    minimum = float(cfg["reference"]["allele_check"]["min_match_fraction"])
    frac = match_fraction(rc["total"]["cis"])
    check(frac >= minimum, f"cis alleles match the reference for {frac:.5f} of variants (minimum {minimum})")
    ac, cis = ref.get("allele_check", {}), rc["total"]["cis"]
    same = (ac.get("a1_is_reference"), ac.get("a2_is_reference"), ac.get("both"), ac.get("neither")) == \
           (cis["a1"], cis["a2"], cis["both"], cis["neither"])
    check(same, "manifest.reference.allele_check equals _checks/refcheck/summary.json")

    # a spot check with real bytes: the refcheck table's call for the paper's named variants has to
    # agree with what the store actually reads at that position. Cheap, and it catches an off-by-one
    # or a stale table that the aggregate counts would sail straight past.
    spot_check(cfg, seqs, store, check)


def spot_check(cfg: Config, seqs: dict, store, check) -> None:
    """For each `paper_variants` entry, read the reference base from the store and confirm the
    refcheck table called the same allele."""
    wrong, checked = [], 0
    for locus in cfg["paper_variants"]:
        chrom, pos = str(locus).split(":")
        pos = int(pos)
        table = refcheck_table(cfg, chrom)
        if chrom not in seqs or not table.exists():
            continue
        t = pq.read_table(table, columns=["position", "A1", "A2", "match"])
        t = t.filter(pc.equal(t["position"], pos))
        if not t.num_rows:
            wrong.append(f"{chrom}:{pos} absent from {table.name}")
            continue
        base = store.get_substring(seqs[chrom]["digest"], pos - 1, pos).upper()
        for row in t.to_pylist():
            called, a1, a2 = row["match"], row["A1"], row["A2"]
            want = {"a1": a1, "a2": a2, "both": a1}
            if called in want and want[called][:1] != base:
                wrong.append(f"{chrom}:{pos} {a1}/{a2} called {called} but the store reads {base}")
            elif called == "neither" and base in (a1, a2):
                wrong.append(f"{chrom}:{pos} {a1}/{a2} called neither but the store reads {base}")
            checked += 1
    check(not wrong, f"refcheck calls match reference bytes for {checked} paper variants ({wrong[:3]})")
