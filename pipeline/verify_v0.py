"""The v1 store against the v0 packs, the migration check for TOPCHeF (store.sbatch runs it).

A bridge module: it reads v0 packs with `packfmt_v0` and `packtool`, and the v1 store with
`catalog`, `results`, `gwas` and `packfmt_v1`. Checks (b) variant catalog cis sites, (c) cis block
codes, (d) trans rows, (e) GWAS rows, (f) the GWAS bin summary against v0's `gwas_dcm_bins.json`. It is the one place allowed to import both codecs
besides `bench_store`, and it goes away with the v0 packs.

    python -m pipeline.verify_v0 --store <store> --id topchef [--genes 40]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from . import catalog as cat
from . import packfmt_v0 as pf0
from . import packfmt_v1 as pf
from . import qtlstore as qs
from .results import load_index, read_block


def verify_v0(store: qs.Store, exp_id: str, chroms: list[str], n_sample: int = 40) -> int:
    """(b) variant catalog cis sites == v0 variants file (A2 = ref, A1 = alt, af code unchanged);
    (c) for sampled genes and introns, u16 nlp and SE codes (sign bit included) and scales equal the
    v0 blocks', matched variant by variant. Returns the number of failures."""
    from . import packtool as pt
    from .common import Config, pack_file
    cfg = Config()
    doc = store.load("experiments", exp_id)
    cdoc = store.load("variant_catalogs", doc["catalog"])
    idx = load_index(store, doc).to_pylist()
    si_v0 = pack_file(cfg, "search_index", "arrow.zst")
    bad = 0
    for c in chroms:
        v0 = pf0.decode_variants_file(pack_file(cfg, f"variants/{c}", "qbv").read_bytes())
        n_cis = v0["n_cis"]
        mine = cat.load_chrom(store, cdoc, c)
        m_cis = cdoc_n_cis(cdoc, c)
        a = sorted(zip(v0["position"][:n_cis].tolist(), v0["A2"][:n_cis], v0["A1"][:n_cis],
                       np.where(np.isnan(v0["af"][:n_cis]), pf0.AF_NULL,
                                np.rint(v0["af"][:n_cis] * pf0.AF_MAXQ)).astype(int).tolist()))
        b = sorted(zip(mine["pos"][:m_cis].tolist(), mine["ref"][:m_cis], mine["alt"][:m_cis],
                       mine["af_code"][:m_cis].astype(int).tolist()))
        ok = a == b
        bad += not ok
        print(f"{'OK  ' if ok else 'FAIL'} (b) {c}: {len(b):,} variant catalog cis sites vs {len(a):,} v0 cis variants "
              f"(pos, ref=A2, alt=A1, af code){'' if ok else f'; {len(set(a) ^ set(b)):,} differ'}", flush=True)
        v0key = list(zip(v0["position"].tolist(), v0["A2"], v0["A1"]))
        mkey = list(zip(mine["pos"].tolist(), mine["ref"], mine["alt"]))
        eq = pack_file(cfg, f"eqtl/{c}", "qbe")
        sq = pack_file(cfg, f"sqtl/{c}", "qbs")
        rows = [r for r in idx if r["chr"] == c and r["has_nominal"]]
        ge = sorted((r for r in rows if r["phenotype_type"] == "ge"), key=lambda r: -(r["n_var"] or 0))[:n_sample]
        tot = {"blocks": 0, "rows": 0, "nlp": 0, "se": 0, "sign": 0, "scales": 0, "missing": 0, "introns": 0}
        by_id = {(r["phenotype_type"], r["phenotype_id"]): r for r in rows}
        for r in ge:
            try:
                off, ln = pt.find_gene_block(eq, r["phenotype_id"], si_v0)
                old = pt.read_block(eq, off, ln, 435)
            except Exception as e:                    # noqa: BLE001
                tot["missing"] += 1
                print(f"  no v0 block for {r['phenotype_id']}: {e}")
                continue
            _cmp(read_block(store, doc, r), old, mkey, v0key, tot, r["phenotype_id"])
            for sp in (old["details"] or {}).get("splice", []):
                mr = by_id.get(("leafcutter", sp["phenotype_id"]))
                if mr is None:
                    continue
                sold = pt.read_block(sq, sp["blk_off"], sp["blk_len"], 480, kind=pf0.KIND_SQTL)
                _cmp(read_block(store, doc, mr), sold, mkey, v0key, tot, sp["phenotype_id"])
                tot["introns"] += 1
        fails = tot["nlp"] + tot["se"] + tot["scales"] + tot["missing"]
        bad += fails
        print(f"{'OK  ' if not fails else 'FAIL'} (c) {c}: {tot['blocks']} blocks ({tot['introns']} introns), "
              f"{tot['rows']:,} rows: nlp code diffs {tot['nlp']}, SE code diffs {tot['se']} "
              f"(sign bit {tot['sign']}), scale diffs {tot['scales']}, missing {tot['missing']}", flush=True)
    return bad


def cdoc_n_cis(cdoc: dict, chrom: str) -> int:
    return next(c["n_cis"] for c in cdoc["chromosomes"] if c["name"] == chrom)


def _cmp(new: dict, old: dict, mkey: list, v0key: list, tot: dict, label: str) -> None:
    tot["blocks"] += 1
    # an all-null row (NLP_NULL and SE_NULL) is dropped on both sides: in v1 it is a vidx the phenotype
    # did not test, and in v0 it is a source row with null slope, SE and p (a monomorphic variant,
    # CONTRACT.md "n_variants"); both formats store those identically, as null codes
    o = {v0key[old["var_start"] + i]: (int(a), int(b)) for i, (a, b) in enumerate(zip(old["nlp_code"], old["se_code"]))
         if not (a == pf0.NLP_NULL and b == pf0.SE_NULL)}
    n = {mkey[new["var_start"] + i]: (int(a), int(b)) for i, (a, b) in enumerate(zip(new["nlp_code"], new["se_code"]))
         if not (a == pf0.NLP_NULL and b == pf0.SE_NULL)}
    tot["rows"] += len(o)
    if set(o) != set(n):
        tot["nlp"] += 1
        print(f"  {label}: variant sets differ ({len(o)} v0, {len(n)} v1)")
        return
    dn = sum(o[k][0] != n[k][0] for k in o)
    ds = sum(o[k][1] != n[k][1] for k in o)
    tot["sign"] += sum((o[k][1] ^ n[k][1]) & pf0.SE_SIGN != 0 for k in o)
    tot["nlp"] += dn > 0
    tot["se"] += ds > 0
    if (new["nlp_max"], new["lse_min"], new["lse_max"]) != (old["nlp_max"], old["lse_min"], old["lse_max"]):
        tot["scales"] += 1
    if dn or ds:
        print(f"  {label}: {dn} nlp and {ds} SE codes differ")


def verify_trans_v0(store: qs.Store, exp_id: str, chroms: list[str], n_sample: int = 40, seed: int = 0) -> int:
    """(d) trans: for sampled genes with a v0 trans frame, the gene's v1 `ge` frame and each intron's v1
    `leafcutter` frame hold v0's rows (same variant chromosome and position), minus exactly the rows whose
    position has no site in the variant catalog (the allele-less trans eQTL variants, left out on
    purpose); -log10 p and beta agree within the two formats' rounding. Returns failures."""
    from .common import Config, pack_file
    from .results import read_trans
    cfg = Config()
    doc = store.load("experiments", exp_id)
    if not any(r.get("trans") for r in doc["results"]):
        print("SKIP (d) trans: the experiment has no trans objects")
        return 0
    cdoc = store.load(qs.CATALOGS, doc["catalog"])
    v1 = {(r["phenotype_type"], r["phenotype_id"]): r for r in load_index(store, doc).to_pylist()}
    from .annotation import decode as decode_arrow
    si0 = decode_arrow(pack_file(cfg, "search_index", "arrow.zst").read_bytes(), "v0 search index").to_pylist()
    names0 = pf0.TRANS_VARIANT_CHROMS
    sites = {}                                   # chromosome -> set of positions with a site
    for c in cdoc["chromosomes"]:
        sites[c["name"]] = set(cat.load_chrom(store, cdoc, c["name"])["pos"].tolist())
    rng = np.random.default_rng(seed)
    bad = 0
    for c in chroms:
        genes = [r for r in si0 if r["chr"] == c and r.get("trans_off") is not None]
        tot = {"genes": 0, "phenotypes": 0, "rows_v0": 0, "rows_v1": 0, "variant_outside_build": 0, "left_out_allele_less": 0,
               "missing_unexplained": 0, "extra_in_v1": 0, "nlp_out_of_bound": 0, "beta_out_of_bound": 0,
               "beta_sign_differs": 0}
        path0 = pack_file(cfg, f"trans/{c}", "qbt")
        with open(path0, "rb") as f0:
            for i in rng.permutation(len(genes))[:n_sample]:
                g = genes[int(i)]
                f0.seek(g["trans_off"])
                fr = pf0.decode_trans_frame(f0.read(g["trans_len"]), 435, 480)
                ids = pf0.trans_phenotype_ids(fr, c, g["gene_id"], g["gene_version"])
                types = ["ge"] * fr["n_e"] + ["leafcutter"] * fr["n_s"]
                e0 = {}
                for t, pid, vc, pos, nlp, beta in zip(types, ids, fr["variant_chr"].tolist(), fr["position"].tolist(),
                                                      fr["nlp"].tolist(), fr["beta"].tolist()):
                    e0.setdefault((t, pid), {})[(names0[vc - 1], int(pos))] = (nlp, beta, fr["nlp_max"], fr["beta_max"])
                tot["genes"] += 1
                for key, rows0 in e0.items():
                    tot["phenotypes"] += 1
                    tot["rows_v0"] += len(rows0)
                    row1 = v1.get(key)
                    fr1 = read_trans(store, doc, row1) if row1 is not None else None
                    rows1 = {} if fr1 is None else {(ch, int(p)): (n, b, fr1["nlp_max"], fr1["beta_max"]) for ch, p, n, b in
                                                    zip(fr1["chr"], fr1["pos"].tolist(), fr1["nlp"].tolist(), fr1["beta"].tolist())}
                    tot["rows_v1"] += len(rows1)
                    tot["extra_in_v1"] += len(set(rows1) - set(rows0))
                    for k0, (n0, b0, nm0, bm0) in rows0.items():
                        if k0 not in rows1:
                            if k0[0] not in sites:
                                tot["variant_outside_build"] += 1
                            elif k0[1] in sites[k0[0]]:
                                tot["missing_unexplained"] += 1
                            else:
                                tot["left_out_allele_less"] += 1
                            continue
                        n1, b1, nm1, bm1 = rows1[k0]
                        if abs(n1 - n0) > (nm0 + nm1) / 131066 + 1e-9:
                            tot["nlp_out_of_bound"] += 1
                        if abs(abs(b1) - abs(b0)) > (bm0 + bm1) / 65534 + 1e-7:
                            tot["beta_out_of_bound"] += 1
                        if (b0 > 0) != (b1 > 0) and abs(b0) > (bm0 + bm1) / 65534:
                            tot["beta_sign_differs"] += 1
        fails = tot["missing_unexplained"] + tot["extra_in_v1"] + tot["nlp_out_of_bound"] + tot["beta_out_of_bound"]
        bad += fails
        print(f"{'OK  ' if not fails else 'FAIL'} (d) {c} trans: {tot}", flush=True)
    return bad


def verify_gwas_v0(store: qs.Store, exp_id: str, chroms: list[str]) -> int:
    """(e) GWAS: every v1 row of a chromosome is a v0 row, as given (ref = NEA, alt = EA) or swapped (ref = EA,
    beta and af mirrored), with identical codes; the v0 rows not in v1 are the ones orientation dropped
    (the reference reads neither allele) and the second copy of an identical pair. Returns failures."""
    from collections import Counter
    from .common import Config, pack_file
    from .gwas import decode_index
    cfg = Config()
    doc = store.load("experiments", exp_id)
    g1 = doc.get("gwas")
    if not g1:
        print("SKIP (e) GWAS: the experiment has no GWAS object")
        return 0
    i0 = pf0.decode_gwas_index(pack_file(cfg, "gwas_index", "bin").read_bytes())
    i1 = decode_index((store.immutable / g1["index"]).read_bytes())
    bad, grand = 0, Counter()
    for c in chroms:
        if c not in i0["chroms"] and c not in g1["files"]:
            continue

        def rows(path, idx, first, dec):
            fp, eo = idx["chroms"][c]
            buf = path.read_bytes()
            out = []
            for k in range(len(eo)):
                b = dec(buf[(first if k == 0 else int(eo[k - 1])):int(eo[k])], idx["n_values"])
                cd = b["codes"]
                al = (b["ea"], b["nea"]) if "ea" in b else (b["ref"], b["alt"])
                out += list(zip(b["position"].tolist(), al[0], al[1], cd["beta"].tolist(), cd["se"].tolist(),
                                (cd["eaf"] if "eaf" in cd else cd["af"]).tolist(), cd["p_mant"].tolist(),
                                cd["p_exp"].tolist(), b["n"].tolist(), b["rs_number"].tolist()))
            return out
        r0 = rows(pack_file(cfg, f"gwas/{c}", "qbg"), i0, pf0.FILE_HEADER_LEN, pf0.decode_gwas_block)
        r1 = rows(store.immutable / g1["files"][c], i1, qs.HEADER_LEN, pf.decode_gwas_block)
        want1 = Counter(r1)
        as_is = Counter((p, nea, ea, b, se, af, m, e, n, rs) for p, ea, nea, b, se, af, m, e, n, rs in r0)
        swap = Counter((p, ea, nea, -b, se, 10000 - af, m, e, n, rs) for p, ea, nea, b, se, af, m, e, n, rs in r0)
        t = Counter()
        for k, v in want1.items():
            a = min(v, as_is.get(k, 0))
            sw = min(v - a, swap.get(k, 0))
            t["as_is"] += a
            t["swapped"] += sw
            t["v1_unmatched"] += v - a - sw
        t["v0_rows"], t["v1_rows"] = len(r0), len(r1)
        t["v0_not_in_v1"] = len(r0) - t["as_is"] - t["swapped"]
        grand += t
        fails = t["v1_unmatched"]
        bad += fails
        print(f"{'OK  ' if not fails else 'FAIL'} (e) {c} GWAS: {dict(t)}", flush=True)
    want_dropped = (g1.get("orientation") or {}).get("dropped")
    print(f"(e) GWAS total: {dict(grand)}; orientation dropped {want_dropped}", flush=True)
    return bad


def verify_gwas_bins_v0(store: qs.Store, exp_id: str, chroms: list[str]) -> int:
    """(f) GWAS bins: the v1 bin summary equals v0's `gwas_dcm_bins.json` on every chromosome the store has
    GWAS files for: the same set of (chr, bin_start) and every value exactly equal (JSON numbers against
    the decoded float64 and uint32; strings and nulls as they are). Returns failures."""
    import json
    from .common import Config
    from .gwas import read_bins
    doc = store.load("experiments", exp_id)
    g1 = doc.get("gwas") or {}
    if not g1.get("bins"):
        print("SKIP (f) GWAS bins: the experiment has no bin summary")
        return 0
    v0 = json.loads((Config().derived / "gwas_dcm_bins.json").read_text())
    cols = v0["columns"]
    keep = set(chroms) & set(g1["files"])
    old = {(cols["chr"][i], cols["bin_start"][i]): {k: v[i] for k, v in cols.items()}
           for i in range(v0["n"]) if cols["chr"][i] in keep}
    new = {(r["chr"], r["bin_start"]): r for r in read_bins(store, doc).to_pylist()}
    bad = [k for k in sorted(set(old) | set(new)) if old.get(k) != new.get(k)]
    fields = {k: sum(1 for key in old if key in new and old[key][k] != new[key][k]) for k in cols}
    print(f"{'OK  ' if not bad else 'FAIL'} (f) GWAS bins: {len(new)} v1 bins, {len(old)} v0 bins on {len(keep)} chromosomes, "
          f"{len(bad)} differ (by field: { {k: v for k, v in fields.items() if v} or 'none'})"
          + (f"; first {bad[:3]}: v0 {[old.get(k) for k in bad[:1]]} v1 {[new.get(k) for k in bad[:1]]}" if bad else ""),
          flush=True)
    return len(bad)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--store", required=True, type=Path)
    ap.add_argument("--id", required=True)
    ap.add_argument("--genes", type=int, default=40)
    args = ap.parse_args(argv)
    from .common import CHROMS
    st = qs.Store(args.store)
    bad = verify_v0(st, args.id, CHROMS, args.genes)
    bad += verify_trans_v0(st, args.id, CHROMS, args.genes)
    bad += verify_gwas_v0(st, args.id, CHROMS)
    bad += verify_gwas_bins_v0(st, args.id, CHROMS)
    print("all comparisons agree" if not bad else f"FAILURES: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
