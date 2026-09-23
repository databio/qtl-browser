"""Reference values for `npm run store-check`: the Python decoders' output for one experiment of a
local qtlstore, written as JSON for the TS decoder (src/lib/store-decode.ts) to match.

    uv run python ui/scripts/store_reference.py <store dir> <out.json> [--experiment topchef] [--sample 60]

Run from the repo root (it imports `pipeline`; the codec is `pipeline/packfmt_v1.py`). Per chromosome:
every site of the variants file (`catalog.decode_file`), the variant index (`decode_vidx`), every rsID
record (`decode_rsid`) plus `rsid_lookup` for a sample of rs numbers, and every hits record frame by
frame (`results.decode_hits`) with its frame table (`packfmt_v1.hits_frame_table`). Per phenotype type:
a seeded sample of blocks plus blocks with credible sets, each through `results.read_block`; every
trans frame through `results.read_trans`; GWAS windows through `gwas.read_window`; and every row of
the GWAS bin summary through `gwas.read_bins`. Floats
go out as JSON numbers (shortest round-trip repr, so exact); NaN as null and infinities as strings.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline import catalog as cat, gwas as gw, packfmt_v1 as pf, qtlstore as qs, results as res  # noqa: E402


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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("store")
    ap.add_argument("out")
    ap.add_argument("--experiment", default="topchef")
    ap.add_argument("--sample", type=int, default=60)
    a = ap.parse_args()
    st = qs.Store(Path(a.store))
    doc = st.load("experiments", a.experiment)
    # the variant catalogs' pointer directory is `variant_catalogs/` (older stores: `catalogs/`)
    cat_dir = qs.CATALOGS
    cdoc = json.loads((st.root / cat_dir / f"{doc['catalog']}.json").read_text())
    names = [c["name"] for c in cdoc["chromosomes"]]
    out = {"experiment": a.experiment, "catalog": cdoc["id"], "catalog_dir": cat_dir, "chroms": {}, "blocks": []}

    vx = cat.decode_vidx((st.immutable / cdoc["vidx"]).read_bytes(), names)
    out["vidx"] = {"page_size": vx["page_size"], "rsid_block_records": vx["rsid_block_records"], "rsid_n": vx["rsid_n"],
                   "rsid_first": ints(vx["rsid_first"]),
                   "chroms": {n: {"n_cis": int(c["n_cis"]), "n_trans": int(c["n_trans"]), "page_off": ints(c["page_off"]),
                                  "page_first_position": ints(c["page_first_position"])} for n, c in vx["chroms"].items()}}
    rs = cat.decode_rsid((st.immutable / cdoc["rsid"]).read_bytes())
    out["rsid"] = {"file": cdoc["rsid"], "rs_number": ints(rs["rs_number"]), "vidx": ints(rs["vidx"]), "ordinal": ints(rs["ordinal"])}
    rbuf = (st.immutable / cdoc["rsid"]).read_bytes()
    uniq = np.unique(rs["rs_number"])
    tries = sorted(set(uniq[::499].tolist()) | {1, int(uniq[-1]) + 1, int(uniq[0]) - 1} | set(int(x) for x in vx["rsid_first"]))
    out["rsid_lookup"] = {str(t): [list(x) for x in cat.rsid_lookup(lambda o, n: rbuf[o:o + n], vx["rsid_first"], vx["rsid_n"], t)]
                          for t in tries if t > 0}

    for c in cdoc["chromosomes"]:
        d = cat.decode_file((st.immutable / c["file"]).read_bytes())
        flags = np.asarray(d["flags"])
        hbuf = (st.immutable / doc["hits"][c["name"]]).read_bytes()
        h = res.decode_hits(hbuf)
        hh = qs.parse_file_header(hbuf)
        table = pf.hits_frame_table(hbuf, hh["n_cis"], hh["page_size"])
        out["chroms"][c["name"]] = {
            "file": c["file"], "seq_digest": c["seq_digest"], "header": d["header"],
            "pos": ints(d["pos"]), "ref": list(d["ref"]), "alt": list(d["alt"]), "af_code": ints(d["af_code"]),
            "af": floats(d["af"]), "rs_number": ints(d["rs_number"]), "ma_samples": ints(d["ma_samples"]),
            "ma_count": ints(d["ma_count"]), "match": ints((flags & cat.FLAG_MATCH_MASK) >> cat.FLAG_MATCH_SHIFT),
            "alt_is_minor": ints(flags & cat.FLAG_ALT_IS_MINOR),
            "hits": {"file": doc["hits"][c["name"]], "header": hh, "frame_off": ints(table),
                     "vidx": ints(h["vidx"]), "ord": ints(h["ord"]), "value": floats(h["value"]), "beta": floats(h["beta"]),
                     "kind": ints(h["kind"]), "cs_id": ints(h["cs_id"]), "flags": ints(h["flags"])},
        }

    # rows as plain Python values (to_pylist keeps None for null integers; pandas would make floats)
    index_rows = res.load_index(st, doc).to_pylist()
    rng = np.random.default_rng(0)
    for ptype in dict.fromkeys(x["phenotype_type"] for x in index_rows):
        # trans-only phenotypes have no block
        g = [x for x in index_rows if x["phenotype_type"] == ptype and x["blk_off"] is not None]
        r = next(x for x in doc["results"] if x["phenotype_type"] == ptype)
        pick = set(rng.choice(len(g), size=min(a.sample, len(g)), replace=False).tolist())
        rows = [row for i, row in enumerate(g) if i in pick]
        # blocks with credible sets and the largest block, so those paths are always covered
        with_cs = []
        for row in g:
            if len(with_cs) >= a.sample // 2:
                break
            if row["ord"] in {x["ord"] for x in rows}:
                continue
            b = res.read_block(st, doc, row)
            if b["n_cs"]:
                with_cs.append(row)
        rows += with_cs + [max(g, key=lambda x: x["blk_len"])]
        seen = set()
        for row in rows:
            if row["ord"] in seen:
                continue
            seen.add(row["ord"])
            b = res.read_block(st, doc, row)
            clean = {k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in row.items()}
            out["blocks"].append({
                "row": clean, "file": r["files"][row["chr"]], "dof": r["dof"],
                "n_rows": b["n_rows"], "var_start": b["var_start"], "pos_first": b["pos_first"], "pos_last": b["pos_last"],
                "nlp_max": f(b["nlp_max"]), "lse_min": f(b["lse_min"]), "lse_max": f(b["lse_max"]), "details": b["details"],
                "nlp_code": ints(b["nlp_code"]), "se_code": ints(b["se_code"]), "nlp": floats(b["nlp"]),
                "pval": floats(b["pval_nominal"]), "se": floats(b["slope_se"]), "slope": floats(b["slope"]),
                "cs_row": ints(b["cs_row"]), "cs_pip": floats(b["cs_pip"]), "cs_id": ints(b["cs_id"]),
            })
    # trans: every frame of every phenotype type
    out["trans"] = []
    for row in (x for x in index_rows if x["trans_off"] is not None):
        fr = res.read_trans(st, doc, row)
        r = next(x for x in doc["results"] if x["phenotype_type"] == row["phenotype_type"])
        out["trans"].append({"file": r["trans"]["file"], "dof": r["dof"], "phenotype_id": row["phenotype_id"],
                             "phenotype_type": row["phenotype_type"], "gene_id": row["gene_id"],
                             "trans_off": int(row["trans_off"]), "trans_len": int(row["trans_len"]), "n_trans": int(row["n_trans"]),
                             "nlp_max": f(fr["nlp_max"]), "beta_max": f(fr["beta_max"]), "chr": fr["chr"],
                             "pos": ints(fr["pos"]), "rs_number": ints(fr["rs_number"]), "af_code": ints(fr["af_code"]),
                             "af": floats(fr["af"]), "nlp_code": ints(fr["nlp_code"]), "nlp": floats(fr["nlp"]), "p": floats(fr["p"]),
                             "beta": floats(fr["beta"]), "se": floats(fr["se"]), "ref": fr["ref"], "alt": fr["alt"]})
    # GWAS: the window of every sampled block's run, and one whole chromosome
    out["gwas"] = []
    if doc.get("gwas"):
        wins = {(b["row"]["chr"], b["row"]["w_lo"], b["row"]["w_hi"]) for b in out["blocks"] if b["row"].get("w_lo") is not None}
        g0 = next(iter(doc["gwas"]["files"]))
        wins.add((g0, 1, 2**31 - 1))
        for chrom, lo, hi in sorted(wins):
            w = gw.read_window(st, doc, chrom, int(lo), int(hi))
            out["gwas"].append({"chr": chrom, "lo": int(lo), "hi": int(hi), "pos": ints(w["pos"]), "ref": list(w["ref"]),
                                "alt": list(w["alt"]), "beta": floats(w["beta"]), "se": floats(w["se"]), "af": floats(w["af"]),
                                "p": floats(w["p"]), "n": ints(w["n"]), "rs_number": ints(w["rs_number"])})
    # the GWAS bin summary: every row
    bins = gw.read_bins(st, doc)
    out["gwas_bins"] = None if bins is None else {**doc["gwas"]["bins"], "rows": bins.to_pylist()}
    Path(a.out).write_text(json.dumps(out, allow_nan=False, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(f"{a.out}: {len(out['blocks'])} blocks, {sum(len(c['pos']) for c in out['chroms'].values())} sites, "
          f"{len(out['rsid']['rs_number'])} rsID records, {sum(len(c['hits']['vidx']) for c in out['chroms'].values())} hits, "
          f"{len(out['trans'])} trans frames, {len(out['gwas'])} GWAS windows, "
          f"{'no' if bins is None else bins.num_rows} GWAS bins")
    return 0


if __name__ == "__main__":
    sys.exit(main())
