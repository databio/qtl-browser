"""Reference values for `npm run store-check`: the Python decoders' output for one collection of a
local qtlstore (format version 2), written as JSON for the TS decoder (src/lib/store-decode.ts) to match.

    uv run python ui-topchef/scripts/store_reference.py <store dir> <out.json> \
        [--anchor topchef_ge_cis] [--gwas dcm_jurgens2024_biobanks] [--sample 40]

Run from the repo root (it imports `pipeline`). The experiments are the anchor's collection siblings,
the ones the app opens with VITE_COLLECTION_OF, plus the GWAS named like VITE_GWAS. A genome-wide
store is too large to dump whole, so every object family gets a seeded sample decoded in full and
exact totals over every row, both of which the TS side recomputes from its own decode:

- variants files (`catalog.decode_file`): per chromosome a sample of pages (columns in full) and
  integer sums over every site; the variant index (`decode_vidx`) whole;
- the rsID index (`decode_rsid`): a sample of blocks in full, sums over every record, and
  `rsid_lookup` for a sample of rs numbers;
- hits (`results.decode_hits`, `packfmt_v1.hits_frame_table`), per experiment and chromosome: the
  frame table, a sample of frames in full, and sums over every record;
- cis blocks (`results.read_block` by seek): a seeded sample per cis experiment, blocks with credible
  sets (found through the hits files' kind 1 records) and the largest block;
- trans frames (`results.read_trans`): a seeded sample per trans experiment and the largest frame;
- GWAS (`gwas.read_window`): the windows of the first sampled blocks, the smallest chromosome whole,
  and sums over every chromosome; every row of the bin summary (`gwas.read_bins`);
- what a first gene page reads: gene lookup rows for a sample of keys (`annotation.lookup_bucket`),
  each chromosome's genes and exon models with a sample of rows, every search index part (and the
  trans experiments' unplaced part) with its ords, ids and offsets, and each experiment's counts.

Floats go out as JSON numbers (shortest round-trip repr, so exact); NaN as null, infinities as strings.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline import annotation as an, catalog as cat, gwas as gw, packfmt_v1 as pf, qtlstore as qs, results as res  # noqa: E402


def f(x):
    x = float(x)
    if math.isnan(x):
        return None
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x


def floats(a):
    return [f(x) for x in np.asarray(a, dtype=np.float64).tolist()]


def ints(a):
    return [int(x) for x in np.asarray(a).tolist()]


def isum(a) -> int:
    return int(np.asarray(a, dtype=np.int64).sum())


def fsum_finite(a) -> tuple[float, int]:
    """(exact sum of the finite entries, how many were not finite)."""
    x = np.asarray(a, dtype=np.float64)
    ok = np.isfinite(x)
    return math.fsum(x[ok].tolist()), int((~ok).sum())


def spread(n: int, k: int) -> list[int]:
    """k indices spread over range(n), first and last included."""
    if n <= k:
        return list(range(n))
    return sorted({round(i * (n - 1) / (k - 1)) for i in range(k)})


def clean_row(row: dict) -> dict:
    return {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in row.items()}


class Blocks:
    """`results.read_block` by seek: the codec from each object's header, read once per file."""

    def __init__(self, store: qs.Store):
        self.store = store
        self.codec: dict[str, object] = {}

    def read(self, doc: dict, row: dict) -> dict:
        name = doc["files"][row["chr"]]
        with open(self.store.immutable / name, "rb") as fh:
            if name not in self.codec:
                self.codec[name] = qs.codec_for(qs.parse_file_header(fh.read(qs.HEADER_LEN))["version"])
                if self.codec[name] is pf:
                    raise SystemExit(f"{name}: a version 1 results object; this reference is for version 2 stores")
            fh.seek(row["blk_off"])
            blk = fh.read(row["blk_len"])
        return self.codec[name].decode_gene_block(blk, expect_blk_len=row["blk_len"], details_version=res.DETAILS_VERSION)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("store")
    ap.add_argument("out")
    ap.add_argument("--anchor", default="topchef_ge_cis", help="an experiment id; its collection's experiments are checked")
    ap.add_argument("--gwas", default="dcm_jurgens2024_biobanks", help="the GWAS experiment id, '' for none")
    ap.add_argument("--sample", type=int, default=40)
    ap.add_argument("--smoke-genes", default="ENSG00000099956,ENSG00000243156",
                    help="gene ids whose every trans frame is included, for bench/smoke_trans_tab.mjs (SMARCB1, MICAL3)")
    ap.add_argument("--smoke-variants", default="rs4819361",
                    help="rs numbers whose trans associations' frames are included, for the same suite")
    a = ap.parse_args()
    st = qs.Store(Path(a.store))
    sdoc = json.loads((st.root / "store.json").read_text())
    if sdoc.get("format_version") != 2:
        raise SystemExit(f"{a.store}: format_version {sdoc.get('format_version')}; this reference is for version 2 stores")
    exps = sdoc["experiments"]
    if a.anchor not in exps:
        raise SystemExit(f"{a.store}: no experiment {a.anchor} in store.json ({', '.join(exps)})")
    cname = exps[a.anchor]["collection"]
    ids = [e for e, x in exps.items() if x["collection"] == cname and x["kind"] != "gwas"]
    docs = {e: st.load("experiments", e) for e in ids}
    gdoc = st.load("experiments", a.gwas) if a.gwas else None
    if gdoc is not None and gdoc["kind"] != "gwas":
        raise SystemExit(f"{a.gwas}: kind {gdoc['kind']}, not gwas")
    cis = [e for e in ids if docs[e]["kind"] == "qtl_cis"]
    trans = [e for e in ids if docs[e]["kind"] == "qtl_trans"]
    cdoc = st.load(qs.CATALOGS, docs[a.anchor]["catalog"])
    names = [c["name"] for c in cdoc["chromosomes"]]
    rng = np.random.default_rng(0)
    out = {"format_version": 2, "anchor": a.anchor, "collection": {"name": cname, "doc": json.loads((st.immutable / cname).read_text())},
           "experiments": {e: {"kind": docs[e]["kind"], "phenotype_type": docs[e].get("phenotype_type")} for e in ids},
           "gwas_experiment": a.gwas or None, "catalog": cdoc["id"], "catalog_dir": qs.CATALOGS}

    # ---- variant index, rsID index ----
    vx = cat.decode_vidx((st.immutable / cdoc["vidx"]).read_bytes(), names)
    out["vidx"] = {"page_size": vx["page_size"], "rsid_block_records": vx["rsid_block_records"], "rsid_n": vx["rsid_n"],
                   "rsid_first": ints(vx["rsid_first"]),
                   "chroms": {n: {"n_cis": int(c["n_cis"]), "n_trans": int(c["n_trans"]), "page_off": ints(c["page_off"]),
                                  "page_first_position": ints(c["page_first_position"])} for n, c in vx["chroms"].items()}}
    rbuf = (st.immutable / cdoc["rsid"]).read_bytes()
    rs = cat.decode_rsid(rbuf)
    B = vx["rsid_block_records"]
    n_blocks = -(-len(rs) // B)
    out["rsid"] = {"file": cdoc["rsid"], "n": int(len(rs)),
                   "totals": {"rs_number": isum(rs["rs_number"]), "vidx": isum(rs["vidx"]), "ordinal": isum(rs["ordinal"])},
                   "blocks": {str(b): {"rs_number": ints(rs["rs_number"][b * B:(b + 1) * B]), "vidx": ints(rs["vidx"][b * B:(b + 1) * B]),
                                       "ordinal": ints(rs["ordinal"][b * B:(b + 1) * B])} for b in spread(n_blocks, 12)}}
    uniq = np.unique(rs["rs_number"])
    tries = sorted(set(uniq[:: max(1, len(uniq) // 400)].tolist()) | {1, int(uniq[-1]) + 1, int(uniq[0]) - 1} | set(int(x) for x in vx["rsid_first"]))
    out["rsid_lookup"] = {str(t): [list(x) for x in cat.rsid_lookup(lambda o, n: rbuf[o:o + n], vx["rsid_first"], vx["rsid_n"], t)]
                          for t in tries if t > 0}

    # ---- variants files: sampled pages whole, sums over every site ----
    out["chroms"] = {}
    for c in cdoc["chromosomes"]:
        d = cat.decode_file((st.immutable / c["file"]).read_bytes())
        flags = np.asarray(d["flags"])
        match = (flags & cat.FLAG_MATCH_MASK) >> cat.FLAG_MATCH_SHIFT
        minor = flags & cat.FLAG_ALT_IS_MINOR
        ref, alt = list(d["ref"]), list(d["alt"])
        n_cis, count, ps = int(c["n_cis"]), int(c["count"]), int(cdoc["page_size"])
        n_pages_cis = -(-n_cis // ps)
        n_pages = n_pages_cis + -(-(count - n_cis) // ps)
        pages = {}
        for k in spread(n_pages, 8):
            lo = k * ps if k < n_pages_cis else n_cis + (k - n_pages_cis) * ps
            hi = min(lo + ps, n_cis if k < n_pages_cis else count)
            s = slice(lo, hi)
            pages[str(k)] = {"first_vidx": lo, "n": hi - lo, "pos": ints(d["pos"][s]), "ref": ref[s], "alt": alt[s],
                             "af_code": ints(d["af_code"][s]), "af": floats(d["af"][s]), "rs_number": ints(d["rs_number"][s]),
                             "ma_samples": ints(d["ma_samples"][s]), "ma_count": ints(d["ma_count"][s]),
                             "match": ints(match[s]), "alt_is_minor": ints(minor[s])}
        out["chroms"][c["name"]] = {
            "file": c["file"], "seq_digest": c["seq_digest"], "header": d["header"], "count": count, "n_cis": n_cis,
            "n_pages": n_pages, "pages": pages,
            "totals": {"pos": isum(d["pos"]), "af_code": isum(d["af_code"]), "rs_number": isum(d["rs_number"]),
                       "ma_samples": isum(d["ma_samples"]), "ma_count": isum(d["ma_count"]), "match": isum(match),
                       "alt_is_minor": isum(minor), "allele_len": sum(map(len, ref)) + sum(map(len, alt)),
                       "allele_first": sum(ord(x[0]) for x in ref) + sum(ord(x[0]) for x in alt)},
        }

    # ---- hits: per experiment and chromosome ----
    out["hits"] = {}
    cs_ords: dict[str, set[int]] = {e: set() for e in cis}
    for e in ids:
        out["hits"][e] = {}
        for c in cdoc["chromosomes"]:
            name = docs[e]["hits"][c["name"]]
            hbuf = (st.immutable / name).read_bytes()
            hh = qs.parse_file_header(hbuf)
            offs = pf.hits_frame_table(hbuf, hh["n_cis"], hh["page_size"])
            h = res.decode_hits(hbuf)
            if docs[e]["kind"] == "qtl_cis":
                cs_ords[e].update(int(x) for x in np.unique(h["ord"][h["kind"] == 1]).tolist())
            nonempty = [g for g in range(len(offs) - 1) if offs[g + 1] > offs[g]]
            frames = {}
            for g in [nonempty[i] for i in spread(len(nonempty), 6)]:
                fr = pf.decode_hits_frame(hbuf[offs[g]:offs[g + 1]], g * hh["page_size"], hh["page_size"], f"{name} frame {g}")
                frames[str(g)] = {"vidx": ints(fr["vidx"]), "ord": ints(fr["ord"]), "value": floats(fr["value"]), "beta": floats(fr["beta"]),
                                  "kind": ints(fr["kind"]), "cs_id": ints(fr["cs_id"]), "flags": ints(fr["flags"])}
            value, value_bad = fsum_finite(h["value"])
            beta, beta_bad = fsum_finite(np.abs(np.asarray(h["beta"], dtype=np.float64)))
            out["hits"][e][c["name"]] = {
                "file": name, "header": hh, "frame_off": ints(offs), "n": int(len(h)),
                "kinds": [int((h["kind"] == k).sum()) for k in (0, 1, 2)], "frames": frames,
                "totals": {"vidx": isum(h["vidx"]), "ord": isum(h["ord"]), "cs_id": isum(h["cs_id"]), "flags": isum(h["flags"]),
                           "value": value, "value_not_finite": value_bad, "abs_beta": beta, "beta_not_finite": beta_bad}}

    # ---- cis blocks ----
    blocks = Blocks(st)
    out["blocks"] = []
    index_rows = {e: res.load_index(st, docs[e]).to_pylist() for e in ids}
    for e in cis:
        g = [x for x in index_rows[e] if x["blk_off"] is not None]
        pick = set(rng.choice(len(g), size=min(a.sample, len(g)), replace=False).tolist())
        rows = [row for i, row in enumerate(g) if i in pick]
        chosen = {x["ord"] for x in rows}
        with_cs = [x for x in g if x["ord"] in cs_ords[e] and x["ord"] not in chosen]
        rows += [with_cs[i] for i in spread(len(with_cs), a.sample // 2)] + [max(g, key=lambda x: x["blk_len"])]
        seen = set()
        for row in rows:
            if row["ord"] in seen:
                continue
            seen.add(row["ord"])
            b = blocks.read(docs[e], row)
            out["blocks"].append({
                "experiment": e, "row": clean_row(row), "file": docs[e]["files"][row["chr"]],
                "n_rows": b["n_rows"], "var_start": b["var_start"], "pos_first": b["pos_first"], "pos_last": b["pos_last"], "n_cs": b["n_cs"],
                "nlp_max": f(b["nlp_max"]), "beta_max": f(b["beta_max"]), "lse_min": f(b["lse_min"]), "lse_max": f(b["lse_max"]),
                "details": b["details"],
                "nlp_code": ints(b["nlp_code"]), "beta_code": ints(b["beta_code"]), "se_code": ints(b["se_code"]),
                "nlp": floats(b["nlp"]), "pval": floats(b["pval_nominal"]), "beta": floats(b["beta"]), "se": floats(b["se"]),
                "cs_row": ints(b["cs_row"]), "cs_pip": floats(b["cs_pip"]), "cs_id": ints(b["cs_id"]),
            })

    # ---- trans frames: a sample, the largest, and the smoke suite's genes and variants ----
    smoke_genes = set(filter(None, a.smoke_genes.split(",")))
    smoke_vidx: dict[str, set[int]] = {}
    for rsid in filter(None, a.smoke_variants.split(",")):
        for ordinal, vidx in cat.rsid_lookup(lambda o, n: rbuf[o:o + n], vx["rsid_first"], vx["rsid_n"], int(rsid.lstrip("rs"))):
            smoke_vidx.setdefault(names[ordinal - 1], set()).add(int(vidx))
    out["trans"] = []
    for e in trans:
        g = [x for x in index_rows[e] if x["trans_off"] is not None]
        pick = set(rng.choice(len(g), size=min(a.sample, len(g)), replace=False).tolist())
        rows = [row for i, row in enumerate(g) if i in pick] + [max(g, key=lambda x: x["trans_len"])]
        rows += [x for x in g if x["gene_id"] in smoke_genes]
        by_ord = {x["ord"]: x for x in g}
        for chrom, vs in smoke_vidx.items():
            h = res.decode_hits((st.immutable / docs[e]["hits"][chrom]).read_bytes())
            rows += [by_ord[int(o)] for o in np.unique(h["ord"][(h["kind"] == 2) & np.isin(h["vidx"], list(vs))]).tolist() if int(o) in by_ord]
        seen = set()
        for row in rows:
            if row["ord"] in seen:
                continue
            seen.add(row["ord"])
            fr = res.read_trans(st, docs[e], row)
            out["trans"].append({"experiment": e, "file": docs[e]["file"], "dof": docs[e]["dof"], "phenotype_id": row["phenotype_id"],
                                 "phenotype_type": row["phenotype_type"], "gene_id": row["gene_id"],
                                 "trans_off": int(row["trans_off"]), "trans_len": int(row["trans_len"]), "n_trans": int(row["n_trans"]),
                                 "n_cs": int(fr["n_cs"]), "nlp_max": f(fr["nlp_max"]), "beta_max": f(fr["beta_max"]),
                                 "lse_min": f(fr["lse_min"]), "lse_max": f(fr["lse_max"]), "chr": fr["chr"],
                                 "pos": ints(fr["pos"]), "rs_number": ints(fr["rs_number"]), "af_code": ints(fr["af_code"]), "af": floats(fr["af"]),
                                 "nlp_code": ints(fr["nlp_code"]), "nlp": floats(fr["nlp"]), "p": floats(fr["p"]),
                                 "beta_code": ints(fr["beta_code"]), "beta": floats(fr["beta"]), "se_code": ints(fr["se_code"]), "se": floats(fr["se"]),
                                 "cs_row": ints(fr["cs_row"]), "cs_pip": floats(fr["cs_pip"]), "cs_id": ints(fr["cs_id"]),
                                 "ref": fr["ref"], "alt": fr["alt"]})

    # ---- GWAS: a few block windows whole, the smallest chromosome whole, sums over every chromosome ----
    out["gwas"] = None
    out["gwas_bins"] = None
    if gdoc is not None:
        def window(chrom, lo, hi):
            w = gw.read_window(st, gdoc, chrom, int(lo), int(hi))
            return {"chr": chrom, "lo": int(lo), "hi": int(hi), "pos": ints(w["pos"]), "ref": list(w["ref"]), "alt": list(w["alt"]),
                    "beta": floats(w["beta"]), "se": floats(w["se"]), "af": floats(w["af"]), "p": floats(w["p"]), "n": ints(w["n"]),
                    "rs_number": ints(w["rs_number"]), "cs_row": ints(w["cs_row"]), "cs_pip": floats(w["cs_pip"]), "cs_id": ints(w["cs_id"])}
        wins = []
        for b in out["blocks"]:
            r = b["row"]
            if r.get("w_lo") is not None and (r["chr"], r["w_lo"], r["w_hi"]) not in wins:
                wins.append((r["chr"], r["w_lo"], r["w_hi"]))
        wins = wins[:12]
        smallest = min(gdoc["rows_by_chrom"], key=lambda k: gdoc["rows_by_chrom"][k])
        wins.append((smallest, 1, 2**31 - 1))
        totals = {}
        for chrom in gdoc["files"]:
            w = gw.read_window(st, gdoc, chrom, 1, 2**31 - 1)
            nlp = -np.log10(np.asarray(w["p"], dtype=np.float64))
            af, af_bad = fsum_finite(w["af"])
            nlp_sum, nlp_bad = fsum_finite(nlp)
            totals[chrom] = {"rows": int(len(w["pos"])), "pos": isum(w["pos"]), "n": isum(w["n"]), "rs_number": isum(w["rs_number"]),
                             "abs_beta": math.fsum(np.abs(np.asarray(w["beta"], dtype=np.float64)).tolist()),
                             "se": math.fsum(np.asarray(w["se"], dtype=np.float64).tolist()),
                             "af": af, "af_not_finite": af_bad, "nlp": nlp_sum, "nlp_not_finite": nlp_bad, "n_cs": int(len(w["cs_row"]))}
        out["gwas"] = {"experiment": a.gwas, "windows": [window(*w) for w in wins], "chrom_totals": totals}
        bins = gw.read_bins(st, gdoc)
        out["gwas_bins"] = None if bins is None else {**gdoc["bins"], "rows": bins.to_pylist()}

    # ---- what a first gene page reads ----
    adoc = st.load("annotations", docs[a.anchor]["annotation"])
    genes = an.decode((st.immutable / adoc["genes"]).read_bytes())
    buckets = an.decode_lookup((st.immutable / adoc["lookup"]).read_bytes())[1]
    gl = genes.to_pylist()
    pick = [gl[i] for i in sorted(rng.choice(len(gl), size=min(300, len(gl)), replace=False).tolist())]
    keys = sorted({g["gene_id"] for g in pick} | {g["name"] for g in pick if g["name"]} | {g["name"].lower() for g in pick[:20] if g["name"]}
                  | {"NOT-A-GENE", "ensg00000000000"})
    out["lookup"] = {"file": adoc["lookup"], "identity": adoc["identity_digest"], "n_buckets": len(buckets),
                     "keys": {k: [r for r in buckets[an.lookup_bucket(an.lookup_key(k), len(buckets))].to_pylist()
                                  if r["key"] == an.lookup_key(k)] for k in keys}}
    # the coloc lists for the collection's QTL experiments against the GWAS (SPEC section 20): each
    # pointer, its rows, and every listed gene's lookup rows, which is how the Home track places them
    out["colocs"] = []
    for p in sorted((st.root / "colocs").glob("*.json")):
        cdoc = json.loads(p.read_text())
        if cdoc["qtl"] in ids and cdoc["gwas"] == a.gwas:
            rows = json.loads((st.immutable / cdoc["object"]).read_bytes())
            out["colocs"].append({"doc": cdoc, "rows": rows, "lookup": {
                g: [r for r in buckets[an.lookup_bucket(an.lookup_key(g), len(buckets))].to_pylist() if r["key"] == an.lookup_key(g)]
                for g in sorted({r["gene_id"] for r in rows})}})
    out["annotation_chroms"] = {}
    for c, parts in adoc["chroms"].items():
        g = an.decode((st.immutable / parts["genes"]).read_bytes())
        m = an.decode((st.immutable / parts["exon_models"]).read_bytes()).to_pylist()
        out["annotation_chroms"][c] = {**parts, "n": g.num_rows, "gene_id": g.column("gene_id").to_pylist(),
                                       "tss": g.column("tss").to_pylist(), "models": m[:: max(1, len(m) // 20)]}
    out["index_parts"] = {}
    out["counts"] = {}
    for e in ids:
        doc = docs[e]
        parts = list(doc["search_index_parts"].items())
        if doc.get("search_index_unplaced"):
            parts.append(("unplaced", doc["search_index_unplaced"]))
        out["index_parts"][e] = {}
        for c, p in parts:
            t = res.decode_arrow((st.immutable / p["file"]).read_bytes(), f"{e} index part {c}").to_pylist()
            part = {**p, "ord": [r["ord"] for r in t], "phenotype_id": [r["phenotype_id"] for r in t], "chr": [r["chr"] for r in t]}
            if doc["kind"] == "qtl_cis":
                part["blk_off"] = [r["blk_off"] for r in t]
            else:
                part["trans_off"] = [r["trans_off"] for r in t]
                part["trans_len"] = [r["trans_len"] for r in t]
                part["n_trans"] = [r["n_trans"] for r in t]
                part["gene_id"] = [r["gene_id"] for r in t]
            out["index_parts"][e][c] = part
        out["counts"][e] = {"n_phenotypes": doc["n_phenotypes"], "counts": doc.get("counts"), "precision": doc["precision"],
                            "n_rows": doc.get("n_rows"), "dof": doc.get("dof"), "significance": doc.get("significance")}

    Path(a.out).write_text(json.dumps(out, allow_nan=False, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    n_sites = sum(c["count"] for c in out["chroms"].values())
    n_hits = sum(x["n"] for e in out["hits"].values() for x in e.values())
    print(f"{a.out}: {len(ids)} experiments of {cname} + GWAS {a.gwas or 'none'}; {len(out['blocks'])} blocks, "
          f"{n_sites} sites in {len(out['chroms'])} chromosomes, {out['rsid']['n']} rsID records, {n_hits} hits, "
          f"{len(out['trans'])} trans frames, {len(out['gwas']['windows']) if out['gwas'] else 0} GWAS windows, "
          f"{'no' if out['gwas_bins'] is None else len(out['gwas_bins']['rows'])} GWAS bins, {len(out['lookup']['keys'])} lookup keys, "
          f"{len(out['colocs'])} coloc lists, "
          f"{len(out['annotation_chroms'])} annotation chromosomes, "
          f"{sum(len(p) for p in out['index_parts'].values())} index parts")
    return 0


if __name__ == "__main__":
    sys.exit(main())
