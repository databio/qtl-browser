"""Round-trip tests for pipeline/packtool.py: rows -> pack -> rows, through the API and the CLI.

    uv run python -m pipeline.test_packtool               # synthetic cases, then a chr21 smoke test when the packs exist
    uv run python -m pipeline.test_packtool --synthetic   # synthetic cases only

Plain asserts, no test framework, in the style of test_packfmt.py. Every case writes into a fresh
temporary directory. The generators are small copies of test_packfmt's, so a change there cannot
break these tests by accident.
"""
from __future__ import annotations

import contextlib
import io
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy.special import stdtr, stdtrit

from . import packfmt as pf
from . import packtool as pt
from .common import Config, addressed_files

CASES = []
DOF, DOF_S = 435, 480


def case(fn):
    CASES.append(fn)
    return fn


def raises(rule: str, fn, *args, **kwargs) -> None:
    try:
        fn(*args, **kwargs)
    except (pt.PackToolError, ValueError) as e:
        assert rule in str(e), f"expected {rule!r} in error, got {e!r}"
        return
    raise AssertionError(f"no error raised (expected {rule!r})")


# ---- generators -------------------------------------------------------------------------------------
NONSNP = [("AT", "A"), ("A", "AT"), ("N", "A"), ("*", "G"), ("A", "A"), ("C", "N"), ("TT", "T"), ("G", "GACGTACGTA")]


def synth_variants(n: int, seed: int = 0) -> pa.Table:
    """n unique cis variants sorted by (position, A1, A2) with nulls in every value column, as a table."""
    rng = np.random.default_rng(seed)
    snps = list(pf.SNP_CODES)
    keys, pos = set(), 10_000
    while len(keys) < n:
        pos += int(rng.integers(0, 40))
        pair = snps[rng.integers(len(snps))] if rng.random() < 0.9 else NONSNP[rng.integers(len(NONSNP))]
        keys.add((pos, *pair))
    keys = sorted(keys)[:n]
    rs = rng.integers(1, 1_600_000_000, n)
    af = rng.random(n)
    ms = rng.integers(0, 517, n)
    mc = rng.integers(0, 517, n)
    null = lambda p: rng.random(n) < p
    return pa.table({
        "position": pa.array([k[0] for k in keys], pa.int64()), "A1": pa.array([k[1] for k in keys]), "A2": pa.array([k[2] for k in keys]),
        "rs_number": pa.array(rs, mask=null(0.05)), "af": pa.array(af, mask=null(0.03)),
        "ma_samples": pa.array(ms, mask=null(0.02)), "ma_count": pa.array(mc, mask=null(0.02)),
        "match": pa.array([["none", "exact", "position"][i] for i in rng.integers(0, 3, n)]),
    })


def synth_trans_only(n: int, start: int = 200_000) -> pa.Table:
    """n trans-only variants at distinct positions; every fifth has no alleles (flags bit 2), the first has no rsID."""
    return pa.table({"position": pa.array([start + 7 * i for i in range(n)], pa.int64()),
                     "A1": pa.array([None if i % 5 == 0 else "A" for i in range(n)], pa.string()),
                     "A2": pa.array([None if i % 5 == 0 else ("G" if i % 3 else "GT") for i in range(n)], pa.string()),
                     "rs_number": pa.array([i * 11 for i in range(n)], pa.int64()), "af": pa.array(np.linspace(0.01, 0.99, n)),
                     "match": pa.array(["position" if i % 5 == 0 else "exact" for i in range(n)])})


def p_of(slope, se, dof: int) -> np.ndarray:
    return 2.0 * stdtr(dof, -np.abs(np.asarray(slope, dtype=np.float64)) / np.asarray(se, dtype=np.float64))


def synth_results(variants: pa.Table, genes: list[tuple[str, int, int]], seed: int = 1, dof: int = DOF) -> pa.Table:
    """Nominal rows for `genes` = [(phenotype_id, first vidx, n rows)], each gene's rows shuffled
    (the writer sorts them by vidx; genes keep their order, which sets the block order), with a
    null row per gene and credible-set rows on the three smallest p."""
    rng = np.random.default_rng(seed)
    v = variants.to_pydict()
    cols = {k: [] for k in ("phenotype_id", "position", "A1", "A2", "tss_distance", "pval_nominal", "slope", "slope_se", "pip", "cs_id")}
    parts = []
    for pid, first, n in genes:
        se = rng.uniform(0.03, 0.3, n)
        sl = rng.normal(0, 0.25, n)
        p = p_of(sl, se, dof)
        p[n // 2] = sl[n // 2] = se[n // 2] = np.nan
        anchor = v["position"][first] - 1000
        pip = [None] * n
        cs = [None] * n
        for j, r in enumerate(np.argsort(np.nan_to_num(p, nan=2.0))[:3]):
            pip[r], cs[r] = float(rng.uniform(0.1, 1)), 1 + (j == 2)
        for i in range(n):
            cols["phenotype_id"].append(pid); cols["position"].append(v["position"][first + i])
            cols["A1"].append(v["A1"][first + i]); cols["A2"].append(v["A2"][first + i])
            cols["tss_distance"].append(v["position"][first + i] - anchor)
            cols["pval_nominal"].append(None if math.isnan(p[i]) else float(p[i]))
            cols["slope"].append(None if math.isnan(sl[i]) else float(sl[i])); cols["slope_se"].append(None if math.isnan(se[i]) else float(se[i]))
            cols["pip"].append(pip[i]); cols["cs_id"].append(cs[i])
        t = pa.table(cols)
        parts.append(t.take(pa.array(rng.permutation(t.num_rows))))
        cols = {k: [] for k in cols}
    return pa.concat_tables(parts)


def synth_gwas(n: int, seed: int = 3) -> pa.Table:
    """GWAS rows at the source's precision on two chromosomes, shuffled."""
    rng = np.random.default_rng(seed)
    inc = rng.integers(0, 40, n)
    inc[2048:2050] = 0
    pos = np.cumsum(inc) + 10_000
    mant, exp = rng.integers(1000, 10000, n), rng.integers(-43, -3, n)
    p = np.array([float(f"{m}e{e}") for m, e in zip(mant, exp)])
    p[:3] = [1.0, 6.396e-40, 0.05]
    pairs = [("A", "G"), ("C", "T"), ("AT", "A"), ("G", "GACGT" * 40), ("T", "C")]
    ea, nea = zip(*[pairs[i] for i in rng.integers(0, len(pairs), n)])
    rs = np.where(rng.random(n) < 0.125, 0, rng.integers(1, 2**31, n))
    t = pa.table({"chr": ["chr21"] * (n * 3 // 5) + ["chr22"] * (n - n * 3 // 5), "position": pos,
                  "ea": list(ea), "nea": list(nea), "rs_number": pa.array(rs, mask=rs == 0),
                  "beta": np.array([float(f"{x:.4f}") for x in rng.normal(0, 0.2, n)]),
                  "se": np.array([float(f"{x:.4f}") for x in rng.uniform(0.0201, 3.3318, n)]),
                  "eaf": np.array([float(f"{x:.4f}") for x in rng.uniform(0, 1, n)]), "p": p,
                  "n": rng.choice([42637, 422920, 937963], n)})
    return t.take(pa.array(rng.permutation(n)))


TRANS_GENES = {"ENSG_A": ("ALPHA", 3), "ENSG_B": ("BETA", 1), "ENSG_C": ("GAMMA", 12)}   # symbol, gene_version


def synth_trans(seed: int = 13) -> pa.Table:
    """Trans rows for three genes on chrT: A has eQTL rows on two variant chromosomes and two introns,
    B has eQTL rows only, C one intron only. Each gene's rows are shuffled; genes keep their order."""
    rng = np.random.default_rng(seed)
    parts = []

    def run(gene, qt, pid, chroms, n):
        rows = {k: [] for k in ("gene_id", "qtl_type", "phenotype_id", "variant_chr", "position", "rs_number", "af", "pval", "beta")}
        for c in chroms:
            pos = np.sort(rng.choice(np.arange(1_000_000, 1_000_000 + 50 * n * 4), n, replace=False))
            for p in pos.tolist():
                rows["gene_id"].append(gene); rows["qtl_type"].append(qt); rows["phenotype_id"].append(pid); rows["variant_chr"].append(c)
                rows["position"].append(int(p)); rows["rs_number"].append(None if rng.random() < 0.2 else int(rng.integers(1, 2**31)))
                rows["af"].append(float(np.float32(rng.uniform(0.01, 0.99)))); rows["pval"].append(float(10.0 ** -rng.uniform(5, 40)))
                rows["beta"].append(float(np.float32(rng.normal(0, 0.4))))
        return pa.table(rows)
    for gene, runs in (("ENSG_A", [("e", "ENSG_A", ["chr1", "chrX"], 6), ("s", "chrT:100:200:clu_5_+:ENSG_A.3", ["chr2"], 4),
                                   ("s", "chrT:100:300:clu_5_+:ENSG_A.3", ["chr2", "chr3"], 3)]),
                       ("ENSG_B", [("e", "ENSG_B", ["chr7"], 5)]),
                       ("ENSG_C", [("s", "chrT:5:9:clu_77_-:ENSG_C.12", ["chrX"], 4)])):
        t = pa.concat_tables([run(gene, *r) for r in runs])
        parts.append(t.take(pa.array(rng.permutation(t.num_rows))))
    return pa.concat_tables(parts)


def tmpdir() -> Path:
    return Path(tempfile.mkdtemp(prefix="packtool_test_"))


def cli(*args) -> tuple[int, str, str]:
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        rc = pt.main([str(a) for a in args])
    return rc, out.getvalue(), err.getvalue()


def tables_equal(a: pa.Table, b: pa.Table, cols: list[str], float_tol: float = 0.0) -> None:
    for c in cols:
        x, y = a[c].to_pylist(), b[c].to_pylist()
        assert len(x) == len(y), (c, len(x), len(y))
        for i, (p, q) in enumerate(zip(x, y)):
            if isinstance(p, float) and isinstance(q, float):
                assert abs(p - q) <= float_tol, (c, i, p, q)
            else:
                assert p == q, (c, i, p, q)


def write_search_index(path: Path, chrom: str, pointers: pa.Table) -> None:
    """A search_index.parquet with the columns the tool reads, from a trans pointer table."""
    gids = pointers["gene_id"].to_pylist()
    pq.write_table(pa.table({"gene_id": gids, "symbol": [TRANS_GENES[g][0] for g in gids], "chr": [chrom] * len(gids),
                             "tss": list(range(len(gids))), "gene_version": [TRANS_GENES[g][1] for g in gids],
                             "trans_off": pointers["trans_off"], "trans_len": pointers["trans_len"],
                             "blk_off": pa.nulls(len(gids), pa.int64()), "blk_len": pa.nulls(len(gids), pa.int64())}), path)


# ---- cases --------------------------------------------------------------------------------------------
@case
def variants_round_trip():
    d = tmpdir()
    src = synth_variants(1300, seed=2)
    shuffled = src.take(pa.array(np.random.default_rng(0).permutation(src.num_rows)))
    qbv = d / "chrT.qbv"
    pages = pt.write_variants_pack(shuffled, qbv, "chrT")
    h = pt.read_header(qbv)
    assert h["kind"] == 1 and h["chrom"] == "chrT" and h["count"] == 1300 and h["n_cis"] == 1300 and h["page_size"] == 512 and h["extension_matches"]
    assert pages["n"].to_pylist() == [512, 512, 276] and pages["first_vidx"].to_pylist() == [0, 512, 1024] and pages["section"].to_pylist() == ["cis"] * 3
    assert pages.equals(pt.page_table(qbv))
    whole = pt.variant_rows(qbv)
    assert whole["vidx"].to_pylist() == list(range(1300)) and not any(whole["no_alleles"].to_pylist())
    tables_equal(whole, src, ["position", "A1", "A2", "rs_number", "ma_samples", "ma_count", "match"])
    tables_equal(whole, src, ["af"], float_tol=0.5 / pf.AF_MAXQ + 1e-12)
    run = pt.variant_rows(qbv, vidx=1020, n=10)
    assert run["vidx"].to_pylist() == list(range(1020, 1030)) and run["position"].to_pylist() == src["position"].to_pylist()[1020:1030]
    cov = pt.variant_rows(qbv, vidx=1020, n=10, whole_pages=True)
    assert cov.num_rows == 512 + 276 and cov["vidx"][0].as_py() == 512, "the covering pages are pages 1 and 2"
    assert pt.covering_pages(pages, 0, 1) == (32, pages["len"][0].as_py()) and pt.covering_pages(pages, 1299, 1)[0] == pages["off"][2].as_py()
    off, ln = pt.covering_pages(pages, 1024, 276)
    rng = pt.variant_rows(qbv, off=off, length=ln)
    assert rng.num_rows == 276 and rng["vidx"][0].as_py() == 1024
    assert pt.variant_rows(qbv, section="trans").num_rows == 0 and pt.variant_rows(qbv, section="cis").num_rows == 1300
    raises("outside", pt.covering_pages, pages, 1300, 1)
    raises("give both", pt.variant_rows, qbv, vidx=3)
    # every table format carries the rows losslessly: the pack rebuilt from it is byte-identical
    first = qbv.read_bytes()
    for ext in (".tsv", ".csv", ".parquet", ".arrow", ".json"):
        t = d / f"variants{ext}"
        pt.write_table(whole.drop_columns(["vidx", "allele_code", "no_alleles"]), t)
        again = d / f"again{ext}.qbv"
        pt.write_variants_pack(pt.read_table(t), again, "chrT")
        assert again.read_bytes() == first, ext
    # a raw-codec file, other page size
    raw = d / "raw.qbv"
    pt.write_variants_pack(src, raw, "chrT", page_size=1000, codec="raw", level=1)
    assert pt.page_table(raw)["codec"].to_pylist() == ["raw", "raw"] and pt.variant_rows(raw).equals(whole)
    assert pt.check_file(raw)["heap_records"] == int((whole["allele_code"].to_numpy() == 0).sum())
    raises("missing the column", pt.write_variants_pack, src.drop_columns(["A2"]), d / "x.qbv", "chrT")
    raises("strict (A1, A2) order", pt.write_variants_pack, pa.concat_tables([src.slice(0, 2), src.slice(0, 1)]), d / "x.qbv", "chrT")
    raises("expected 2 (eqtl), 3 (sqtl)", pt.list_blocks, qbv)

    # the trans-only section (SPEC section 4): appended after the cis pages, which do not change
    tr = synth_trans_only(700)
    both = d / "both.qbv"
    pages2 = pt.write_variants_pack(shuffled, both, "chrT", trans_only=tr.take(pa.array(np.random.default_rng(1).permutation(700))))
    h2 = pt.read_header(both)
    assert (h2["count"], h2["n_cis"]) == (2000, 1300)
    assert pages2["section"].to_pylist() == ["cis"] * 3 + ["trans"] * 2 and pages2["first_vidx"].to_pylist() == [0, 512, 1024, 1300, 1812]
    cis_end = pages2["off"][3].as_py()
    assert both.read_bytes()[32:cis_end] == first[32:], "adding the section changes no cis byte"
    assert pt.variant_rows(both, section="cis").equals(whole)
    t2 = pt.variant_rows(both, section="trans")
    assert t2.num_rows == 700 and t2["vidx"][0].as_py() == 1300 and t2["position"].to_pylist() == tr["position"].to_pylist()
    assert t2["no_alleles"].to_pylist() == [i % 5 == 0 for i in range(700)] and t2["A1"].null_count == 140 and t2["A2"].null_count == 140
    assert t2["ma_samples"].null_count == 700 and t2["ma_count"].null_count == 700, "the trans files report no sample counts"
    assert t2["rs_number"][0].as_py() is None and t2["rs_number"].to_pylist()[1:] == [i * 11 for i in range(1, 700)]
    assert t2["match"].to_pylist() == tr["match"].to_pylist() and t2["allele_code"].to_pylist()[:4] == [0, 2, 2, 0]
    tables_equal(t2, tr, ["af"], float_tol=0.5 / pf.AF_MAXQ + 1e-12)
    assert pt.variant_rows(both).num_rows == 2000
    straddle = pt.variant_rows(both, vidx=1298, n=4)
    assert straddle["vidx"].to_pylist() == [1298, 1299, 1300, 1301] and straddle["no_alleles"].to_pylist() == [False, False, True, False]
    c = pt.check_file(both)
    assert (c["records"], c["n_cis"], c["trans_only"], c["no_alleles"], c["pages"]) == (2000, 1300, 700, 140, 5)
    # the pipeline's shape: one table with in_cis, and the same through a TSV
    tr_full = tr.append_column("ma_samples", pa.nulls(700, pa.int64())).append_column("ma_count", pa.nulls(700, pa.int64()))
    combined = pa.concat_tables([src.append_column("in_cis", pa.array([True] * 1300)),
                                 tr_full.select(src.column_names).append_column("in_cis", pa.array([False] * 700))], promote_options="default")
    pt.write_variants_pack(combined, d / "combined.qbv", "chrT")
    assert (d / "combined.qbv").read_bytes() == both.read_bytes()
    pt.write_table(combined, d / "combined.tsv")
    pt.write_variants_pack(pt.read_table(d / "combined.tsv"), d / "combined2.qbv", "chrT")
    assert (d / "combined2.qbv").read_bytes() == both.read_bytes()
    raises("strictly increasing", pt.write_variants_pack, src, d / "x.qbv", "chrT", trans_only=pa.concat_tables([tr.slice(0, 2), tr.slice(0, 1)]))
    bad = src.set_column(src.schema.get_field_index("A1"), "A1", pa.array([None] + src["A1"].to_pylist()[1:], pa.string()))
    raises("null allele", pt.write_variants_pack, bad, d / "x.qbv", "chrT")


def build_results(d: Path):
    src = synth_variants(1300, seed=2)
    qbv = d / "chrT.qbv"
    pt.write_variants_pack(src, qbv, "chrT", trans_only=synth_trans_only(50))     # a trans-only section never gets in the way
    introns = [("chrT:11000:12000:clu_1_+:ENSG_A.1", 10, 80), ("chrT:11500:13000:clu_2_+:ENSG_A.1", 40, 700)]
    srows = synth_results(src, introns, seed=5, dof=DOF_S)         # sQTL rows: the source relation uses the sQTL dof
    sptr = pt.write_results_pack(srows, qbv, d / "chrT.qbs", pf.KIND_SQTL)
    genes = [("ENSG_A", 0, 100), ("ENSG_B", 50, 600)]
    erows = synth_results(src, genes, seed=7)
    splice = [{"phenotype_id": pid, "intron_start": int(pid.split(":")[1]), "intron_end": int(pid.split(":")[2]),
               "blk_off": int(sptr["blk_off"][i].as_py()), "blk_len": int(sptr["blk_len"][i].as_py())} for i, (pid, _, _) in enumerate(introns)]
    details = {"ENSG_A": {"v": 0, "gene": {"gene_id": "ENSG_A", "symbol": "ALPHA", "tss": 11000}, "exons": [[11000, 11200], [11800, 12000]], "splice": splice},
               "ENSG_C": {"gene": {"gene_id": "ENSG_C", "symbol": "GAMMA"}, "exons": [], "splice": []}}
    eptr = pt.write_results_pack(erows, qbv, d / "chrT.qbe", pf.KIND_EQTL, details=details)
    return src, qbv, srows, sptr, erows, eptr, details


def check_rows_against_source(t: pa.Table, rows: pa.Table, pid: str, blk: dict, dof: int) -> None:
    s = rows.filter(pc.equal(rows["phenotype_id"], pid)).sort_by([("position", "ascending"), ("A1", "ascending"), ("A2", "ascending")])
    assert t.num_rows == s.num_rows
    tables_equal(t, s, ["position", "A1", "A2", "tss_distance"])
    p_src = np.array([np.nan if x is None else x for x in s["pval_nominal"].to_pylist()])
    p_out = np.array([np.nan if x is None else x for x in t["pval_nominal"].to_pylist()])
    assert np.array_equal(np.isnan(p_src), np.isnan(p_out))
    ok = ~np.isnan(p_src)
    assert np.all(np.abs(-np.log10(p_out[ok]) + np.log10(p_src[ok])) <= blk["nlp_max"] / 131066 + 1e-9), "p within half a step"
    sl_src = np.array([np.nan if x is None else x for x in s["slope"].to_pylist()])
    se_src = np.array([np.nan if x is None else x for x in s["slope_se"].to_pylist()])
    sl_out = np.array([np.nan if x is None else x for x in t["slope"].to_pylist()])
    se_out = np.array([np.nan if x is None else x for x in t["slope_se"].to_pylist()])
    se_b, sl_b = pf.error_bounds(p_src, sl_src, se_src, blk["nlp_max"], blk["lse_min"], blk["lse_max"], dof)
    dom = ~np.isnan(sl_b)
    f32 = 6e-8
    assert np.all(np.abs(se_out[dom] - se_src[dom]) <= pf.BOUND_FACTOR * se_b[dom] + f32 * se_src[dom]), "SE within SPEC section 9"
    assert np.all(np.abs(sl_out[dom] - sl_src[dom]) <= pf.BOUND_FACTOR * sl_b[dom] + f32 * np.abs(sl_src[dom]) + 1e-12), "slope within SPEC section 9"
    assert np.array_equal(np.isnan(sl_out), ~dom) and np.array_equal(np.isnan(se_out), np.isnan(se_src))
    pip_src = s["pip"].to_pylist()
    pip_out = t["pip"].to_pylist()
    assert [x is None for x in pip_src] == [x is None for x in pip_out]
    assert all(abs(a - b) < 1e-6 for a, b in zip(pip_src, pip_out) if a is not None)
    assert t["cs_id"].to_pylist() == s["cs_id"].to_pylist()


@case
def results_round_trip():
    d = tmpdir()
    src, qbv, srows, sptr, erows, eptr, details = build_results(d)
    qbe, qbs = d / "chrT.qbe", d / "chrT.qbs"
    assert pt.read_header(qbe)["count"] == 3 and pt.read_header(qbs)["count"] == 2
    lb = pt.list_blocks(qbe, gene_ids=True)
    assert lb["gene_id"].to_pylist() == ["ENSG_A", "ENSG_B", "ENSG_C"] and lb["symbol"].to_pylist() == ["ALPHA", None, "GAMMA"]
    assert lb["n_rows"].to_pylist() == [100, 600, 0] and lb["var_start"].to_pylist() == [0, 50, None] and lb["n_cs"].to_pylist() == [3, 3, 0]
    assert lb["blk_off"].to_pylist() == eptr["blk_off"].to_pylist() and lb["blk_len"].to_pylist() == eptr["blk_len"].to_pylist()
    assert eptr["phenotype_id"].to_pylist() == ["ENSG_A", "ENSG_B", "ENSG_C"] and eptr["n_var"].to_pylist() == [100, 600, 0]
    assert eptr["var_start"].to_pylist() == [0, 50, None] and eptr["anchor"].to_pylist()[2] is None
    assert pt.list_blocks(qbe, limit=1).num_rows == 1 and "gene_id" not in pt.list_blocks(qbs).column_names
    for gid in ("ENSG_A", "ALPHA"):
        assert pt.find_gene_block(qbe, gid) == (int(eptr["blk_off"][0].as_py()), int(eptr["blk_len"][0].as_py()))
    raises("no block whose details", pt.find_gene_block, qbe, "ENSG_Z")
    for i, gid in enumerate(("ENSG_A", "ENSG_B")):
        off, ln = pt.find_gene_block(qbe, gid)
        b = pt.read_block(qbe, off, ln, DOF)
        assert b["details"]["gene"]["gene_id"] == gid and b["details"]["v"] == 0
        t = pt.block_rows(b, qbv)
        assert t["vidx"].to_pylist() == list(range(b["var_start"], b["var_start"] + b["n_rows"]))
        check_rows_against_source(t, erows, gid, b, DOF)
        codes = pt.block_codes(b)
        assert codes.num_rows == b["n_rows"] and codes["pval_nominal"].null_count == 1 and codes["nlp_code"].to_pylist().count(pf.NLP_MAXQ) >= 1
        cs = pt.block_credible_sets(b)
        assert cs.num_rows == 3 and cs["vidx"].to_pylist() == (cs["row"].to_numpy() + b["var_start"]).tolist()
    off, ln = pt.find_gene_block(qbe, "ENSG_C")
    c = pt.read_block(qbe, off, ln, DOF)
    assert c["n_rows"] == 0 and c["details"]["gene"]["symbol"] == "GAMMA" and pt.block_rows(c, qbv).num_rows == 0
    assert pt.block_rows(c, qbv).column_names[0] == "vidx"
    # introns: found through the gene's details, decoded as kind 3
    for i, pid in enumerate(sptr["phenotype_id"].to_pylist()):
        off, ln = pt.find_intron_block(qbe, pid)
        assert (off, ln) == (int(sptr["blk_off"][i].as_py()), int(sptr["blk_len"][i].as_py()))
        b = pt.read_block(qbs, off, ln, DOF_S)
        assert b["details"] is None and b["n_rows"] == [80, 700][i]
        check_rows_against_source(pt.block_rows(b, qbv), srows, pid, b, DOF_S)
    raises("does not end in an ENSG", pt.find_intron_block, qbe, "chrT:1:2:clu_9_+:XYZ")
    raises("none named", pt.find_intron_block, qbe, "chrT:1:2:clu_9_+:ENSG_A.1")
    raises("needs a details frame", pt.read_block, qbs, off, ln, DOF_S, kind=pf.KIND_EQTL)
    # whole-file checks
    ce = pt.check_file(qbe, variants=qbv, dof=DOF)
    assert (ce["blocks"], ce["rows"], ce["credible_set_records"], ce["empty_blocks"]) == (3, 700, 6, 1)
    cs_ = pt.check_file(qbs, variants=qbv)
    assert (cs_["blocks"], cs_["rows"], cs_["credible_set_records"]) == (2, 780, 6)
    # writer rules
    bad = erows.set_column(erows.schema.get_field_index("position"), "position", pa.array([1] * erows.num_rows, pa.int64()))
    raises("is not in the cis section", pt.write_results_pack, bad, qbv, d / "x.qbe", pf.KIND_EQTL)
    gap = erows.filter(pc.not_equal(erows["position"], src["position"][3]))
    raises("unbroken run", pt.write_results_pack, gap, qbv, d / "x.qbe", pf.KIND_EQTL)
    tss = erows["tss_distance"].to_numpy().copy()
    tss[0] += 1
    raises("not one value", pt.write_results_pack, erows.set_column(erows.schema.get_field_index("tss_distance"), "tss_distance", pa.array(tss)), qbv, d / "x.qbe", pf.KIND_EQTL)
    raises("has no details", pt.write_results_pack, srows, qbv, d / "x.qbs", pf.KIND_SQTL, details={"a": {}})
    nocs = erows.set_column(erows.schema.get_field_index("cs_id"), "cs_id", pa.nulls(erows.num_rows, pa.int64()))
    raises("no cs_id", pt.write_results_pack, nocs, qbv, d / "x.qbe", pf.KIND_EQTL)
    raises("not a results kind", pt.write_results_pack, erows, qbv, d / "x.qbe", pf.KIND_GWAS)
    raises("missing the column", pt.write_results_pack, erows.drop_columns(["slope_se"]), qbv, d / "x.qbe", pf.KIND_EQTL)


def check_trans_rows(t: pa.Table, src: pa.Table, f: dict, gene_id: str) -> None:
    """Decoded frame rows against the source rows of one gene, matched by (qtl_type, phenotype_id, variant_chr, position)."""
    s = src.filter(pc.equal(src["gene_id"], gene_id))
    assert t.num_rows == s.num_rows == f["n_e"] + f["n_s"]
    key = lambda tb, i: (tb["qtl_type"][i].as_py(), tb["phenotype_id"][i].as_py(), tb["variant_chr"][i].as_py(), tb["position"][i].as_py())
    src_rows = {key(s, i): i for i in range(s.num_rows)}
    assert t["qtl_type"].to_pylist() == ["e"] * f["n_e"] + ["s"] * f["n_s"], "eQTL rows first"
    for i in range(t.num_rows):
        j = src_rows.pop(key(t, i))
        r, o = s.slice(j, 1).to_pylist()[0], t.slice(i, 1).to_pylist()[0]
        assert o["gene_id"] == gene_id and o["rs_number"] == r["rs_number"]
        assert abs(o["af"] - r["af"]) <= 0.5 / pf.AF_MAXQ + 1e-12
        assert abs(-math.log10(o["pval"]) + math.log10(r["pval"])) <= f["nlp_max"] / 131066 + 1e-9, "p within half a step"
        assert abs(o["beta"] - r["beta"]) <= f["beta_max"] / 65534 + 1.2e-7 * abs(r["beta"]), "beta within half a step"
        dof = DOF if r["qtl_type"] == "e" else DOF_S
        t_src = -stdtrit(dof, r["pval"] / 2)
        assert abs(o["beta_se"] - abs(r["beta"]) / t_src) <= 0.01 * abs(r["beta"]) / t_src, "beta_se within 1%"
        assert abs(o["r2"] - t_src * t_src / (t_src * t_src + dof)) <= 0.002
        if r["qtl_type"] == "s":
            st, en, clu, strand = pt.parse_phenotype_id(r["phenotype_id"])
            assert (o["intron_start"], o["intron_end"], o["cluster"], o["strand"]) == (st, en, clu, strand)
        else:
            assert o["intron_start"] is None and o["strand"] is None
    assert not src_rows


@case
def trans_round_trip():
    d = tmpdir()
    src = synth_trans()
    qbt = d / "chrT.qbt"
    ptr = pt.write_trans_pack(src, qbt, "chrT", level=3)
    h = pt.read_header(qbt)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"], h["n_cis"]) == (6, "chrT", 3, 0, None) and h["extension_matches"]
    assert ptr["gene_id"].to_pylist() == ["ENSG_A", "ENSG_B", "ENSG_C"] and ptr["trans_off"][0].as_py() == 32
    assert ptr["n_e"].to_pylist() == [12, 5, 0] and ptr["n_s"].to_pylist() == [10, 0, 4]
    assert h["bytes"] == 32 + sum(ptr["trans_len"].to_pylist())
    si = d / "search_index.parquet"
    write_search_index(si, "chrT", ptr)
    walked = pt.trans_frames(qbt)
    assert walked["trans_off"].to_pylist() == ptr["trans_off"].to_pylist() and walked["trans_len"].to_pylist() == ptr["trans_len"].to_pylist()
    listed = pt.trans_frames(qbt, si)
    assert listed["gene_id"].to_pylist() == ["ENSG_A", "ENSG_B", "ENSG_C"] and listed["symbol"].to_pylist() == ["ALPHA", "BETA", "GAMMA"]
    assert pt.trans_frames(qbt, limit=2).num_rows == 2 and pt.trans_frames(qbt, si, limit=1)["gene_version"].to_pylist() == [3]
    for gene, (symbol, version) in TRANS_GENES.items():
        for name in (gene, symbol):
            info = pt.find_trans_frame(qbt, name, si)
            assert info["gene_id"] == gene and info["gene_version"] == version
        f = pt.read_trans_frame(qbt, info["trans_off"], info["trans_len"], DOF, DOF_S)
        t = pt.trans_table(f, "chrT", gene, version)
        check_trans_rows(t, src, f, gene)
        assert t["phenotype_id"].to_pylist()[:f["n_e"]] == [gene] * f["n_e"]
        assert all(p.startswith("chrT:") and p.endswith(f"{gene}.{version}") for p in t["phenotype_id"].to_pylist()[f["n_e"]:])
        anon = pt.trans_table(f, "chrT")
        assert anon["phenotype_id"].null_count == t.num_rows and anon["gene_id"].null_count == t.num_rows
    fa = pt.read_trans_frame(qbt, ptr["trans_off"][0].as_py(), ptr["trans_len"][0].as_py(), DOF, DOF_S)
    assert fa["k"] == 2 and fa["intron_start"].tolist() == [100, 100] and fa["intron_end"].tolist() == [200, 300] and fa["strand"] == ["+", "+"]
    raises("no trans rows", pt.find_trans_frame, qbt, "ENSG_Z", pq.write_table(pa.table({"gene_id": ["ENSG_Z"], "symbol": ["Z"], "chr": ["chrT"],
           "trans_off": pa.nulls(1, pa.int64()), "trans_len": pa.nulls(1, pa.int64()), "gene_version": [1]}), d / "si2.parquet") or d / "si2.parquet")
    raises("no gene", pt.find_trans_frame, qbt, "NOPE", si)
    c = pt.check_file(qbt)
    assert (c["frames"], c["rows"], c["eqtl_rows"], c["sqtl_rows"], c["max_introns"], c["found_through"]) == (3, 31, 17, 14, 2, "zstd frame boundaries")
    c2 = pt.check_file(qbt, search_index=si)
    assert c2["rows"] == 31 and c2["found_through"] == "search_index"
    # writer rules
    raises("must be in (0, 1]", pt.write_trans_pack, src.set_column(src.schema.get_field_index("pval"), "pval", pa.array([1.5] + src["pval"].to_pylist()[1:])), d / "x.qbt", "chrT")
    raises("chr1..chr22 or chrX", pt.write_trans_pack, src.set_column(src.schema.get_field_index("variant_chr"), "variant_chr", pa.array(["chrM"] + src["variant_chr"].to_pylist()[1:])), d / "x.qbt", "chrT")
    raises("not chr:start:end", pt.write_trans_pack, src.set_column(src.schema.get_field_index("phenotype_id"), "phenotype_id", pa.array(["bad:id"] * src.num_rows)), d / "x.qbt", "chrT")
    raises("need phenotype_id", pt.write_trans_pack, src.drop_columns(["phenotype_id"]), d / "x.qbt", "chrT")
    raises("missing the column", pt.write_trans_pack, src.drop_columns(["af"]), d / "x.qbt", "chrT")
    raises("qtl_type must be", pt.write_trans_pack, src.set_column(src.schema.get_field_index("qtl_type"), "qtl_type", pa.array(["x"] * src.num_rows)), d / "x.qbt", "chrT")
    # the same rows through a TSV, and with the four intron columns instead of phenotype_id, are byte-identical
    pt.write_table(src, d / "trans.tsv")
    pt.write_trans_pack(pt.read_table(d / "trans.tsv"), d / "again.qbt", "chrT", level=3)
    assert (d / "again.qbt").read_bytes() == qbt.read_bytes()
    parsed = [pt.parse_phenotype_id(p) if q == "s" else (None, None, None, None) for p, q in zip(src["phenotype_id"].to_pylist(), src["qtl_type"].to_pylist())]
    four = src.drop_columns(["phenotype_id"])
    for i, name in enumerate(("intron_start", "intron_end", "cluster", "strand")):
        four = four.append_column(name, pa.array([x[i] for x in parsed]))
    pt.write_trans_pack(four, d / "four.qbt", "chrT", level=3)
    assert (d / "four.qbt").read_bytes() == qbt.read_bytes()


@case
def gwas_round_trip():
    d = tmpdir()
    src = synth_gwas(5000)
    blocks = pt.write_gwas_packs(src, d / "gwas", block_rows=2048, level=3)
    assert sorted(p.name for p in (d / "gwas").iterdir()) == ["chr21.qbg", "chr22.qbg", "gwas_index.bin"]
    assert blocks["chr"].to_pylist() == ["chr21", "chr21", "chr22"] and blocks["block"].to_pylist() == [0, 1, 0]
    assert blocks.equals(pt.gwas_index_table(d / "gwas" / "gwas_index.bin"))
    idx = pt.gwas_index(d / "gwas" / "gwas_index.bin")
    assert idx["block_rows"] == 2048 and idx["n_values"] == [42637, 422920, 937963]
    h = pt.read_header(d / "gwas" / "chr21.qbg")
    assert h["kind"] == 4 and h["count"] == 3000 and h["page_size"] == 2048
    rs0 = pa.array(np.nan_to_num(src["rs_number"].cast(pa.float64()).to_numpy(zero_copy_only=False), nan=0).astype(np.int64))
    ordered = src.set_column(src.schema.get_field_index("rs_number"), "rs_number", rs0)     # the writer's order: null rs_number sorts as 0
    ordered = ordered.sort_by([("chr", "ascending"), ("position", "ascending"), ("ea", "ascending"), ("nea", "ascending"),
                               ("rs_number", "ascending"), ("p", "ascending")])
    s21 = ordered.filter(pc.equal(ordered["chr"], "chr21"))
    pos = s21["position"].to_numpy()
    dup = int(pos[2047]) if pos[2047] == pos[2048] else None
    assert dup is not None, "the generator repeats a position across the first block boundary"
    for lo, hi in [(dup, dup + 50), (dup - 1, dup), (int(pos[0]), int(pos[0])), (int(pos[-1]) - 5, int(pos[-1]) + 100), (1, 5), (10**9, 10**9 + 1)]:
        got = pt.gwas_rows(d / "gwas" / "chr21.qbg", d / "gwas" / "gwas_index.bin", lo=lo, hi=hi)
        got = got.set_column(got.schema.get_field_index("rs_number"), "rs_number", pc.fill_null(got["rs_number"], 0))
        want = s21.filter(pa.array((pos >= lo) & (pos <= hi)))
        assert got.num_rows == want.num_rows, (lo, hi, got.num_rows, want.num_rows)
        tables_equal(got, want, ["position", "ea", "nea", "rs_number", "beta", "se", "eaf", "n"])
        gp, wp = np.array(got["p"].to_pylist()), np.array(want["p"].to_pylist())
        assert np.all(np.abs(gp - wp) <= 1e-12 * wp)
    b1 = pt.gwas_rows(d / "gwas" / "chr21.qbg", d / "gwas" / "gwas_index.bin", block=1)
    assert b1.num_rows == 3000 - 2048 and b1["block"].to_pylist() == [1] * b1.num_rows
    raises("outside 0..1", pt.gwas_rows, d / "gwas" / "chr21.qbg", d / "gwas" / "gwas_index.bin", block=2)
    raises("give --lo", pt.gwas_rows, d / "gwas" / "chr21.qbg", d / "gwas" / "gwas_index.bin")
    c = pt.check_file(d / "gwas" / "chr22.qbg", index=d / "gwas" / "gwas_index.bin")
    assert (c["blocks"], c["rows"]) == (1, 2000)
    raises("give --index", pt.check_file, d / "gwas" / "chr22.qbg")
    ci = pt.check_file(d / "gwas" / "gwas_index.bin")
    assert ci["chromosomes"] == ["chr21", "chr22"] and ci["blocks"] == {"chr21": 2, "chr22": 1}
    # lossless rules come from packfmt and name the row
    bad = src.set_column(src.schema.get_field_index("p"), "p", pa.array([0.12345] + src["p"].to_pylist()[1:]))
    raises("significant digits", pt.write_gwas_packs, bad, d / "bad")
    raises("missing the column", pt.write_gwas_packs, src.drop_columns(["eaf"]), d / "bad")
    # the pack rebuilt from a TSV of the source is byte-identical
    pt.write_table(src, d / "gwas.tsv")
    pt.write_gwas_packs(pt.read_table(d / "gwas.tsv"), d / "gwas2", block_rows=2048, level=3)
    for name in ("chr21.qbg", "chr22.qbg", "gwas_index.bin"):
        assert (d / "gwas2" / name).read_bytes() == (d / "gwas" / name).read_bytes(), name


@case
def cli_round_trip():
    d = tmpdir()
    src, qbv, srows, sptr, erows, eptr, details = build_results(d)
    qbe, qbs = d / "chrT.qbe", d / "chrT.qbs"
    rc, out, _ = cli("header", qbe)
    assert rc == 0 and json.loads(out)["count"] == 3
    rc, out, _ = cli("blocks", qbe, "--gene-ids", "-o", d / "blocks.parquet")
    assert rc == 0 and pt.read_table(d / "blocks.parquet").equals(pt.list_blocks(qbe, gene_ids=True))
    rc, out, err = cli("block", qbe, "--gene", "ENSG_B")
    assert rc == 1 and "degrees of freedom" in err, "no manifest anywhere near the temp dir"
    (d / "manifest.json").write_text(json.dumps({"packs": {"dof": {"eqtl": DOF, "sqtl": DOF_S}}}))
    rc, out, _ = cli("block", qbe, "--gene", "ENSG_B", "--variants", qbv, "-o", d / "b.tsv")
    assert rc == 0
    off, ln = pt.find_gene_block(qbe, "ENSG_B")
    api = pt.block_rows(pt.read_block(qbe, off, ln, DOF), qbv)
    got = pt.read_table(d / "b.tsv")
    assert got.num_rows == 600 and got.column_names == api.column_names
    tables_equal(got, api, ["vidx", "position", "A1", "A2", "rs_number", "tss_distance", "ma_samples", "ma_count", "cs_id"])
    tables_equal(got, api, ["pval_nominal", "slope", "slope_se", "pip"], float_tol=1e-6)   # float32 columns through their shortest decimal form
    rc, out, _ = cli("block", qbe, "--index", "2", "--details", "--dof", DOF)
    assert rc == 0 and json.loads(out)["gene"]["symbol"] == "GAMMA"
    rc, out, _ = cli("block", qbe, "--off", off, "--len", ln, "--header", "--manifest", d / "manifest.json")
    assert rc == 0 and json.loads(out)["n_rows"] == 600 and json.loads(out)["dof"] == DOF
    rc, out, _ = cli("block", qbe, "--gene", "ENSG_A", "--cs", "--dof", DOF)
    assert rc == 0 and out.startswith("row\tvidx\tpip\tcs_id\n") and out.count("\n") == 4
    rc, out, _ = cli("block", qbe, "--gene", "ENSG_A", "--dof", DOF)
    assert rc == 0 and out.startswith("row\tvidx\tnlp_code\tse_code\t") and out.count("\n") == 101, "codes without --variants"
    pid = sptr["phenotype_id"][1].as_py()
    rc, out, _ = cli("block", qbs, "--intron", pid, "--eqtl", qbe, "--variants", qbv, "--dof", DOF_S, "-o", d / "i.parquet")
    assert rc == 0 and pt.read_table(d / "i.parquet").num_rows == 700
    rc, _, err = cli("block", qbs, "--intron", pid, "--dof", DOF_S)
    assert rc == 1 and "--eqtl" in err
    rc, out, _ = cli("variants", qbv, "--vidx", 1290, "--n", 10, "-o", d / "v.json")
    assert rc == 0 and [r["vidx"] for r in json.loads((d / "v.json").read_text())] == list(range(1290, 1300))
    rc, out, _ = cli("variants", qbv, "--pages")
    assert rc == 0 and out.splitlines()[0] == "page\toff\tlen\tfirst_vidx\tn\tcodec\tsection" and len(out.splitlines()) == 5
    rc, out, _ = cli("variants", qbv, "--section", "trans")
    assert rc == 0 and len(out.splitlines()) == 51 and out.splitlines()[1].startswith("1300\t")
    rc, out, _ = cli("check", qbs, "--variants", qbv)
    assert rc == 0 and json.loads(out)["rows"] == 780
    # writers from tables on disk
    cis = pt.variant_rows(qbv, section="cis").drop_columns(["vidx", "allele_code", "no_alleles"])
    pt.write_table(cis, d / "variants.tsv")
    pt.write_table(synth_trans_only(50), d / "trans_only.tsv")
    rc, out, _ = cli("pack-variants", d / "variants.tsv", "-o", d / "cli.qbv", "--chrom", "chrT", "--trans-only", d / "trans_only.tsv", "--pages", d / "pages.tsv")
    assert rc == 0 and (d / "cli.qbv").read_bytes() == qbv.read_bytes() and json.loads(out)["pages"] == 4 and json.loads(out)["n_cis"] == 1300
    assert pt.read_table(d / "pages.tsv").num_rows == 4
    pt.write_table(erows, d / "eqtl_rows.parquet")
    (d / "details.json").write_text(json.dumps(details))
    rc, out, _ = cli("pack-results", d / "eqtl_rows.parquet", "--variants", qbv, "-o", d / "cli.qbe", "--details", d / "details.json", "--pointers", d / "ptr.tsv")
    assert rc == 0 and (d / "cli.qbe").read_bytes() == qbe.read_bytes() and json.loads(out)["rows"] == 700
    tables_equal(pt.read_table(d / "ptr.tsv"), eptr, ["phenotype_id", "blk_off", "blk_len", "var_start", "n_var", "anchor", "n_cs"])
    pt.write_table(srows, d / "sqtl_rows.tsv")
    rc, out, _ = cli("pack-results", d / "sqtl_rows.tsv", "--variants", qbv, "-o", d / "cli.qbs")
    assert rc == 0 and (d / "cli.qbs").read_bytes() == qbs.read_bytes()
    rc, _, err = cli("pack-results", d / "sqtl_rows.tsv", "--variants", qbv, "-o", d / "cli.pack")
    assert rc == 1 and "cannot tell the kind" in err
    # trans pack
    trans = synth_trans()
    pt.write_table(trans, d / "trans.tsv")
    rc, out, _ = cli("pack-trans", d / "trans.tsv", "-o", d / "chrT.qbt", "--chrom", "chrT", "--pointers", d / "tptr.tsv", "--level", 3)
    assert rc == 0 and json.loads(out)["frames"] == 3 and json.loads(out)["sqtl_rows"] == 14
    tptr = pt.read_table(d / "tptr.tsv")
    write_search_index(d / "si.parquet", "chrT", tptr)
    rc, out, _ = cli("frames", d / "chrT.qbt", "--search-index", d / "si.parquet")
    assert rc == 0 and out.splitlines()[0] == "frame\tgene_id\tsymbol\tgene_version\ttrans_off\ttrans_len" and len(out.splitlines()) == 4
    rc, out, _ = cli("frames", d / "chrT.qbt", "--limit", 1)
    assert rc == 0 and len(out.splitlines()) == 2
    rc, out, _ = cli("trans", d / "chrT.qbt", "--gene", "ALPHA", "--search-index", d / "si.parquet", "-o", d / "alpha.parquet")
    assert rc == 0
    got = pt.read_table(d / "alpha.parquet")
    f = pt.read_trans_frame(d / "chrT.qbt", tptr["trans_off"][0].as_py(), tptr["trans_len"][0].as_py(), DOF, DOF_S)
    check_trans_rows(got, trans, f, "ENSG_A")
    rc, out, _ = cli("trans", d / "chrT.qbt", "--off", tptr["trans_off"][2].as_py(), "--len", tptr["trans_len"][2].as_py(), "--gene-id", "ENSG_C", "--gene-version", 12,
                     "--dof-eqtl", DOF, "--dof-sqtl", DOF_S, "-o", d / "gamma.tsv")
    assert rc == 0 and pt.read_table(d / "gamma.tsv")["phenotype_id"].to_pylist() == ["chrT:5:9:clu_77_-:ENSG_C.12"] * 4
    rc, out, _ = cli("trans", d / "chrT.qbt", "--gene", "ENSG_A", "--search-index", d / "si.parquet", "--header")
    assert rc == 0 and json.loads(out)["k"] == 2 and json.loads(out)["introns"][1]["intron_end"] == 300 and json.loads(out)["dof"] == {"eqtl": DOF, "sqtl": DOF_S}
    rc, _, err = cli("trans", d / "chrT.qbt", "--gene", "ENSG_A")
    assert rc == 1 and "--search-index" in err
    rc, _, err = cli("trans", d / "chrT.qbt", "--off", 32, "--len", 10, "--dof-eqtl", DOF)
    assert rc == 1 and "both --dof-eqtl and --dof-sqtl" in err
    rc, out, _ = cli("check", d / "chrT.qbt", "--search-index", d / "si.parquet")
    assert rc == 0 and json.loads(out)["rows"] == 31
    # GWAS
    g = synth_gwas(3000)
    pt.write_table(g, d / "gwas.parquet")
    rc, out, _ = cli("pack-gwas", d / "gwas.parquet", "-o", d / "g", "--level", 3, "--blocks", d / "gblocks.tsv")
    assert rc == 0 and json.loads(out)["files"] == ["chr21.qbg", "chr22.qbg", "gwas_index.bin"]
    rc, out, _ = cli("gwas", d / "g" / "chr22.qbg", "--index", d / "g" / "gwas_index.bin", "--block", 0, "-o", d / "g22.arrow")
    assert rc == 0 and pt.read_table(d / "g22.arrow").num_rows == 1200
    rc, out, _ = cli("gwas-index", d / "g" / "gwas_index.bin")
    assert rc == 0 and out.splitlines()[0] == "chr\tblock\tfirst_position\tbyte_start\tbyte_end"
    rc, _, err = cli("header", d / "nope.qbv")
    assert rc == 1 and "no such file" in err


# ---- real data ------------------------------------------------------------------------------------
# ---- hits pack, rsID index, variant index (kinds 7, 8, 9) ----------------------------------------------
HITS_GENES = {0: ("ENSG_A", "ALPHA", 3), 1: ("ENSG_B", "BETA", 1), 2: ("ENSG_C", "GAMMA", 12)}   # ord -> gene_id, symbol, version
HITS_INTRONS = {0: (11000, 12000, 5, "+"), 2: (5, 9, 77, "-")}                                   # ord -> intron of its sQTL rows
RESERVED_VIDX = (1, 100, 300, 2100, 2599)   # the hand-written rows own these; 1 must stay empty, so the spread avoids them


def synth_hits(n_variants: int = 2600, seed: int = 21) -> pa.Table:
    """Hit rows over one chromosome: every kind, variants with no rows, two credible sets of one
    phenotype, and rows in the last (short) frame. Shuffled, as the writer sorts."""
    rng = np.random.default_rng(seed)
    rows = {k: [] for k in ("vidx", "kind", "gene", "pval", "beta", "slope", "slope_se", "pip", "cs_id",
                            "significant", "phenotype_id")}

    def add(vidx, kind, ord_, **kw):
        gid, _, ver = HITS_GENES[ord_]
        pid = gid
        if kind % 2 == 1:
            s, e, c, st = HITS_INTRONS[ord_]
            pid = f"chrT:{s}:{e}:clu_{c}_{st}:{gid}.{ver}"
        rows["vidx"].append(int(vidx)); rows["kind"].append(int(kind)); rows["gene"].append(int(ord_))
        rows["phenotype_id"].append(pid)
        for c in ("pval", "beta", "slope", "slope_se", "pip", "cs_id"):
            rows[c].append(kw.get(c))
        rows["significant"].append(bool(kw.get("significant", False)))

    # one variant with all six kinds, in the first frame
    add(100, 0, 0, pval=1e-30, beta=0.9); add(100, 0, 1, pval=1e-12, beta=-0.42)
    add(100, 1, 0, pval=1e-8, beta=0.11)
    add(100, 2, 0, pval=1e-5, slope=0.55, slope_se=0.08, significant=True)
    add(100, 3, 2, pval=2e-4, slope=-1.2, slope_se=0.3)
    add(100, 4, 0, pip=0.9, cs_id=1); add(100, 5, 2, pip=0.25, cs_id=2)
    # a variant in two credible sets of one phenotype keeps both rows
    add(300, 4, 1, pip=0.7, cs_id=1); add(300, 4, 1, pip=0.3, cs_id=2)
    # rows in the last, short frame (2600 variants = two full frames of 1024 and a 552-variant tail)
    add(2599, 0, 1, pval=1.0, beta=0.01)
    add(2100, 3, 0, pval=3e-9, slope=0.4, slope_se=0.05, significant=True)
    # a spread of trans rows, so frames 0 and 1 both carry rows and many variants carry none
    pool = np.setdiff1d(np.arange(n_variants), np.array(RESERVED_VIDX))
    for v in rng.choice(pool, 120, replace=False).tolist():
        ord_ = int(rng.integers(0, 3))
        kind = 1 if (ord_ in HITS_INTRONS and rng.random() < 0.5) else 0
        add(v, kind, ord_, pval=float(10.0 ** -rng.uniform(5, 40)), beta=float(np.float32(rng.normal(0, 0.4))))
    t = pa.table({"vidx": pa.array(rows["vidx"], pa.int64()), "kind": pa.array(rows["kind"], pa.int64()),
                  "gene": pa.array(rows["gene"], pa.int64()),
                  "pval": pa.array([np.nan if x is None else x for x in rows["pval"]], pa.float64()),
                  "beta": pa.array([np.nan if x is None else x for x in rows["beta"]], pa.float64()),
                  "slope": pa.array([np.nan if x is None else x for x in rows["slope"]], pa.float64()),
                  "slope_se": pa.array([np.nan if x is None else x for x in rows["slope_se"]], pa.float64()),
                  "pip": pa.array([np.nan if x is None else x for x in rows["pip"]], pa.float64()),
                  "cs_id": pa.array(rows["cs_id"], pa.int64()),
                  "significant": pa.array(rows["significant"], pa.bool_()),
                  "phenotype_id": pa.array(rows["phenotype_id"], pa.string())})
    return t.take(pa.array(rng.permutation(t.num_rows)))


def write_hits_search_index(path: Path) -> Path:
    """A search_index.parquet whose (chr, tss, gene_id) order is the ord of HITS_GENES."""
    pq.write_table(pa.table({"gene_id": [HITS_GENES[i][0] for i in range(3)],
                             "symbol": [HITS_GENES[i][1] for i in range(3)], "chr": ["chrT"] * 3,
                             "tss": [0, 1, 2], "gene_version": [HITS_GENES[i][2] for i in range(3)]}), path)
    return path


@case
def hits_round_trip():
    d = tmpdir()
    n_variants = 2600
    src = synth_hits(n_variants)
    qbh = d / "chrT.qbh"
    ptr = pt.write_hits_pack(src, qbh, "chrT", n_variants=n_variants, level=3)
    h = pt.read_header(qbh)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"], h["n_cis"]) == (7, "chrT", 3, 1024, None) and h["extension_matches"]
    assert ptr["frame"].to_pylist() == [0, 1, 2] and ptr["first_vidx"].to_pylist() == [0, 1024, 2048]
    assert ptr["n_variants"].to_pylist() == [1024, 1024, 552], "the last frame is short"
    assert sum(ptr["rows"].to_pylist()) == src.num_rows and ptr["hits_off"][0].as_py() == 32
    assert h["bytes"] == 32 + sum(ptr["hits_len"].to_pylist())
    # frames walked by their zstd boundaries agree with the writer's pointers
    walked = pt.hits_frames(qbh)
    assert walked["hits_off"].to_pylist() == ptr["hits_off"].to_pylist()
    assert walked["hits_len"].to_pylist() == ptr["hits_len"].to_pylist()
    assert pt.hits_frames(qbh, limit=2).num_rows == 2
    si = write_hits_search_index(d / "search_index.parquet")
    # every row comes back, per variant, with its gene and phenotype rebuilt
    got = []
    for v in sorted(set(src["vidx"].to_pylist())):
        t = pt.hits_rows(qbh, vidx=v, search_index=si)
        assert t["vidx"].to_pylist() == [v] * t.num_rows
        assert t["kind"].to_pylist() == sorted(t["kind"].to_pylist()), "a variant's rows go by kind"
        got.append(t)
    out = pa.concat_tables(got)
    assert out.num_rows == src.num_rows
    key = [("vidx", "ascending"), ("kind", "ascending"), ("gene", "ascending"), ("cs_id", "ascending")]
    s = src.sort_by(key)
    o = out.sort_by(key)
    tables_equal(o, s, ["vidx", "kind", "gene", "phenotype_id", "significant"])
    assert o["cs_id"].to_pylist() == [None if x is None else x for x in s["cs_id"].to_pylist()]
    assert o["gene_id"].to_pylist() == [HITS_GENES[g][0] for g in s["gene"].to_pylist()]
    assert o["symbol"].to_pylist() == [HITS_GENES[g][1] for g in s["gene"].to_pylist()]
    assert o["qtl_type"].to_pylist() == ["s" if k % 2 else "e" for k in s["kind"].to_pylist()]
    # values within SPEC section 13's limits, per frame scale
    frames = {i: pt.read_hits_frame(qbh, off, ln, i * 1024)
              for i, (off, ln) in enumerate(zip(ptr["hits_off"].to_pylist(), ptr["hits_len"].to_pylist()))}
    for i in range(s.num_rows):
        v, k = s["vidx"][i].as_py(), s["kind"][i].as_py()
        f = frames[v // 1024]
        if k <= 1:
            assert abs(-math.log10(o["pval"][i].as_py()) + math.log10(s["pval"][i].as_py())) <= f["trans_nlp_max"] / 131066 + 1e-9
            assert abs(o["beta"][i].as_py() - s["beta"][i].as_py()) <= f["trans_beta_max"] / 65534 + 1e-12
            assert o["slope"][i].as_py() is None and o["pip"][i].as_py() is None
        elif k <= 3:
            assert abs(-math.log10(o["pval"][i].as_py()) + math.log10(s["pval"][i].as_py())) <= f["perm_nlp_max"] / 131066 + 1e-9
            assert abs(o["slope_se"][i].as_py() - s["slope_se"][i].as_py()) <= f["se_max"] / 131070 + 1e-12
            assert abs(o["slope"][i].as_py() - s["slope"][i].as_py()) <= f["slope_max"] / 65534 + 1e-12
            assert o["beta"][i].as_py() is None
        else:
            assert abs(o["pip"][i].as_py() - s["pip"][i].as_py()) <= 0.5 / 65535 + 1e-12
            assert o["pval"][i].as_py() is None and o["beta"][i].as_py() is None
    # intron fields are null on the even kinds and rebuilt on the odd ones
    for i in range(o.num_rows):
        k = o["kind"][i].as_py()
        want = HITS_INTRONS[o["gene"][i].as_py()] if k % 2 else None
        if want is None:
            assert o["intron_start"][i].as_py() is None and o["strand"][i].as_py() is None
        else:
            assert (o["intron_start"][i].as_py(), o["intron_end"][i].as_py(), o["cluster"][i].as_py(), o["strand"][i].as_py()) == want
    # a variant with no rows, and a whole frame
    empty = pt.hits_rows(qbh, vidx=1, search_index=si)
    assert empty.num_rows == 0 and "vidx" in empty.column_names
    whole = pt.hits_rows(qbh, frame=0, search_index=si)
    assert whole.num_rows == ptr["rows"][0].as_py()
    # without a search_index the ord stays, with no gene id
    anon = pt.hits_rows(qbh, vidx=100)
    assert anon["gene_id"].null_count == anon.num_rows and anon["phenotype_id"].null_count == anon.num_rows
    assert anon["gene"].to_pylist() == [0, 1, 0, 0, 2, 0, 2]
    # check_file
    c = pt.check_file(qbh)
    assert (c["frames"], c["rows"], c["variants"], c["frame_variants"]) == (3, src.num_rows, n_variants, 1024)
    assert c["found_through"] == "zstd frame boundaries"
    assert sum(c["rows_by_kind"].values()) == src.num_rows and c["rows_by_kind"]["credible set eQTL"] == 3
    # the four intron columns instead of phenotype_id give byte-identical output
    parsed = [pt.parse_phenotype_id(p) if k % 2 else (None, None, None, None)
              for p, k in zip(src["phenotype_id"].to_pylist(), src["kind"].to_pylist())]
    four = src.drop_columns(["phenotype_id"])
    for i, name in enumerate(("intron_start", "intron_end", "cluster", "strand")):
        four = four.append_column(name, pa.array([x[i] for x in parsed]))
    pt.write_hits_pack(four, d / "four.qbh", "chrT", n_variants=n_variants, level=3)
    assert (d / "four.qbh").read_bytes() == qbh.read_bytes()
    # the same rows through a TSV are byte-identical too
    pt.write_table(src, d / "hits.tsv")
    pt.write_hits_pack(pt.read_table(d / "hits.tsv"), d / "again.qbh", "chrT", n_variants=n_variants, level=3)
    assert (d / "again.qbh").read_bytes() == qbh.read_bytes()
    # writer rules
    raises("needs at least one variant", pt.write_hits_pack, src.slice(0, 0), d / "x.qbh", "chrT")
    raises("outside 0..", pt.write_hits_pack, src, d / "x.qbh", "chrT", n_variants=200)
    raises("missing the column", pt.write_hits_pack, src.drop_columns(["gene"]), d / "x.qbh", "chrT", n_variants=n_variants)
    raises("needs a phenotype_id", pt.write_hits_pack, src.set_column(
        src.schema.get_field_index("phenotype_id"), "phenotype_id", pa.nulls(src.num_rows, pa.string())),
        d / "x.qbh", "chrT", n_variants=n_variants)
    raises("need phenotype_id", pt.write_hits_pack, src.drop_columns(["phenotype_id"]), d / "x.qbh", "chrT", n_variants=n_variants)
    raises("p must be in", pt.write_hits_pack, src.set_column(
        src.schema.get_field_index("pval"), "pval", pa.array([2.0] * src.num_rows, pa.float64())),
        d / "x.qbh", "chrT", n_variants=n_variants)
    raises("exactly one of --vidx and --frame", pt.hits_rows, qbh)
    raises("exactly one of --vidx and --frame", pt.hits_rows, qbh, vidx=1, frame=0)
    raises("is negative", pt.hits_rows, qbh, vidx=-1)
    raises("outside the file", pt.hits_rows, qbh, frame=9)
    # the CLI
    rc, out, _ = cli("header", qbh)
    assert rc == 0 and json.loads(out)["kind_name"] == "hits"
    rc, out, _ = cli("hits", qbh, "--vidx", 100, "--search-index", si, "-o", d / "v.parquet")
    assert rc == 0 and pt.read_table(d / "v.parquet").num_rows == 7
    rc, out, _ = cli("hits", qbh, "--frames", "-o", "-")
    assert rc == 0 and out.splitlines()[0] == "frame\tfirst_vidx\thits_off\thits_len" and len(out.splitlines()) == 4
    rc, out, _ = cli("check", qbh)
    assert rc == 0 and json.loads(out)["rows"] == src.num_rows
    rc, out, _ = cli("pack-hits", d / "hits.tsv", "-o", d / "cli.qbh", "--chrom", "chrT",
                     "--n-variants", n_variants, "--level", 3, "--pointers", d / "ptr.tsv")
    assert rc == 0 and (d / "cli.qbh").read_bytes() == qbh.read_bytes()
    return d, qbh, ptr, n_variants


@case
def rsid_index_round_trip():
    d = tmpdir()
    rng = np.random.default_rng(31)
    n = 3 * 64 + 5                                   # four blocks, the last one short
    rs = np.sort(rng.choice(np.arange(1, 2_153_660_727), n, replace=False))
    chrom = [pf.VARIANT_CHROMS[i % 23] for i in range(n)]
    vidx = rng.integers(0, 700_000, n)
    src = pa.table({"rs_number": pa.array(rs, pa.int64()), "chr": pa.array(chrom, pa.string()),
                    "vidx": pa.array(vidx, pa.int64())})
    qbr = d / "rsid_index.qbr"
    blocks = pt.write_rsid_index(src.take(pa.array(rng.permutation(n))), qbr, block_records=64)   # rows are sorted by the writer
    h = pt.read_header(qbr)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (8, "all", n, 64) and h["extension_matches"]
    assert h["bytes"] == 32 + 8 * n
    assert blocks["block"].to_pylist() == [0, 1, 2, 3] and blocks["first_rs_number"].to_pylist() == rs[::64].tolist()
    assert blocks["byte_start"].to_pylist() == [32, 32 + 512, 32 + 1024, 32 + 1536]
    assert blocks["byte_end"].to_pylist()[-1] == h["bytes"]
    # every record comes back
    all_rows = pt.rsid_records(qbr)
    assert all_rows.num_rows == n and all_rows["rs_number"].to_pylist() == rs.tolist()
    assert all_rows["chr"].to_pylist() == chrom and all_rows["vidx"].to_pylist() == vidx.tolist()
    assert all_rows["block"].to_pylist() == [i // 64 for i in range(n)]
    last = pt.rsid_records(qbr, block=3)
    assert last.num_rows == 5 and last["record"].to_pylist() == list(range(192, 197))
    assert pt.rsid_records(qbr, block=0, limit=2).num_rows == 2
    # the two-level lookup, with and without the startup file's samples
    for i in (0, 1, 63, 64, n - 1):
        got = pt.rsid_lookup(qbr, int(rs[i]))
        assert got == {"rs_number": int(rs[i]), "chr": chrom[i], "vidx": int(vidx[i]), "block": i // 64,
                       "byte_start": 32 + (i // 64) * 512, "byte_end": 32 + min((i // 64 + 1) * 512, 8 * n)}
    assert pt.rsid_lookup(qbr, int(rs[0]) - 1) is None, "below every block"
    miss = int(rs[0]) + 1
    if miss in set(rs.tolist()):
        miss = int(rs[-1]) + 5
    assert pt.rsid_lookup(qbr, miss) is None, "a miss inside a block that exists"
    assert pt.rsid_lookup(qbr, int(rs[-1]) + 10_000) is None, "above every record"
    c = pt.check_file(qbr)
    assert (c["records"], c["blocks"], c["block_records"], c["last_rs_number"]) == (n, 4, 64, int(rs[-1]))
    # writer rules
    raises("chr must be chr1..chr22 or chrX", pt.write_rsid_index,
           src.set_column(src.schema.get_field_index("chr"), "chr", pa.array(["chrM"] * n)), d / "x.qbr")
    raises("must strictly increase", pt.write_rsid_index,
           src.set_column(src.schema.get_field_index("rs_number"), "rs_number", pa.array([7] * n, pa.int64())), d / "x.qbr")
    raises("vidx", pt.write_rsid_index,
           src.set_column(src.schema.get_field_index("vidx"), "vidx", pa.array([1 << 27] * n, pa.int64())), d / "x.qbr")
    raises("missing the column", pt.write_rsid_index, src.drop_columns(["vidx"]), d / "x.qbr")
    # chr_ordinal instead of chr
    ordinals = src.drop_columns(["chr"]).append_column(
        "chr_ordinal", pa.array([pf.VARIANT_CHROMS.index(c) + 1 for c in chrom], pa.int64()))
    pt.write_rsid_index(ordinals, d / "ord.qbr", block_records=64)
    assert (d / "ord.qbr").read_bytes() == qbr.read_bytes()
    # the CLI
    rc, out, _ = cli("rsid", qbr, "--rs", int(rs[70]))
    assert rc == 0 and json.loads(out)["chr"] == chrom[70] and json.loads(out)["block"] == 1
    rc, out, _ = cli("rsid", qbr, "--rs", int(rs[0]) - 1)
    assert rc == 1, "a miss is a non-zero exit"
    rc, out, _ = cli("rsid", qbr, "--block", 3, "-o", d / "b.tsv")
    assert rc == 0 and pt.read_table(d / "b.tsv").num_rows == 5
    pt.write_table(src, d / "rsid.tsv")
    rc, out, _ = cli("pack-rsid", d / "rsid.tsv", "-o", d / "cli.qbr", "--block-records", 64, "--blocks", d / "bl.tsv")
    assert rc == 0 and (d / "cli.qbr").read_bytes() == qbr.read_bytes()
    return d, qbr, src


@case
def variant_index_round_trip():
    """A whole little build: 23 variants files, 23 hits packs, an rsID index, and the startup file."""
    d = tmpdir()
    (d / "variants").mkdir(); (d / "hits").mkdir()
    P, F, B = 8, 4, 64
    rng = np.random.default_rng(41)
    want, rsid_rows = {}, {"rs_number": [], "chr": [], "vidx": []}
    for i, name in enumerate(pf.VARIANT_CHROMS):
        n_cis, n_tr = 30 + i, 5 + (i % 3)
        cis = synth_variants(n_cis, seed=100 + i)
        tr = synth_trans_only(n_tr, start=500_000 + 1000 * i)
        qbv = d / "variants" / f"{name}.qbv"
        pt.write_variants_pack(cis, qbv, name, page_size=P, trans_only=tr, level=3)
        total = n_cis + n_tr
        # hit rows on a few variants of each chromosome, including the last one
        rows = {"vidx": [], "kind": [], "gene": [], "pval": [], "beta": [], "phenotype_id": []}
        for v in sorted(set(rng.choice(np.arange(total), 6, replace=False).tolist() + [total - 1])):
            rows["vidx"].append(v); rows["kind"].append(0); rows["gene"].append(0)
            rows["pval"].append(float(10.0 ** -rng.uniform(5, 30))); rows["beta"].append(float(rng.normal(0, 0.3)))
            rows["phenotype_id"].append("ENSG_A")
        pt.write_hits_pack(pa.table(rows), d / "hits" / f"{name}.qbh", name, n_variants=total, frame_variants=F, level=3)
        want[name] = {"n_cis": n_cis, "n_trans_only": n_tr, "total": total}
        # rsIDs come from what the file itself decodes to, so a record always points at its own variant
        recs = pa.concat_tables([pt.variant_rows(qbv, section="cis"), pt.variant_rows(qbv, section="trans")])
        for vv, rr in zip(recs["vidx"].to_pylist(), recs["rs_number"].to_pylist()):
            if rr is not None:
                rsid_rows["rs_number"].append(int(rr)); rsid_rows["chr"].append(name); rsid_rows["vidx"].append(int(vv))
    rt = pa.table({k: pa.array(v, pa.int64() if k != "chr" else pa.string()) for k, v in rsid_rows.items()})
    rt = rt.filter(pc.is_valid(rt["rs_number"]))
    # rs_number must not repeat across the whole index
    seen, keep = set(), []
    for j, r in enumerate(rt["rs_number"].to_pylist()):
        if r not in seen:
            seen.add(r); keep.append(j)
    rt = rt.take(pa.array(keep))
    qbr = d / "rsid_index.qbr"
    pt.write_rsid_index(rt, qbr, block_records=B)
    variants = {c: d / "variants" / f"{c}.qbv" for c in pf.VARIANT_CHROMS}
    hits = {c: d / "hits" / f"{c}.qbh" for c in pf.VARIANT_CHROMS}
    qbx = d / "variant_index.qbx"
    tbl = pt.write_variant_index(qbx, variants, hits, qbr, level=3)
    h = pt.read_header(qbx)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (9, "all", 23, P) and h["extension_matches"]
    assert tbl.num_rows == 23 and tbl["chr"].to_pylist() == list(pf.VARIANT_CHROMS)
    idx = pt.read_variant_index(qbx)
    assert (idx["page_size"], idx["frame_variants"], idx["rsid_block_records"]) == (P, F, B)
    assert idx["rsid_n_records"] == rt.num_rows and idx["rsid_n_blocks"] == -(-rt.num_rows // B)
    for name in pf.VARIANT_CHROMS:
        c, w = idx["chroms"][name], want[name]
        assert (c["n_cis"], c["n_trans_only"]) == (w["n_cis"], w["n_trans_only"])
        # the offsets are the files' own
        pages = pt.page_table(variants[name])
        assert c["page_off"].tolist() == pages["off"].to_pylist() + [pt.read_header(variants[name])["bytes"]]
        assert c["n_pages_cis"] + c["n_pages_trans"] == pages.num_rows
        fr = pt.hits_frames(hits[name])
        assert c["hits_off"].tolist() == fr["hits_off"].to_pylist() + [pt.read_header(hits[name])["bytes"]]
        # and the startup file finds the same frames as walking the pack does
        byidx = pt.hits_frames(hits[name], variant_index=qbx)
        assert byidx["hits_off"].to_pylist() == fr["hits_off"].to_pylist()
        assert byidx["hits_len"].to_pylist() == fr["hits_len"].to_pylist()
        # every variant's page and frame resolve, and the page really holds it
        for vidx in (0, w["n_cis"] - 1, w["n_cis"], w["total"] - 1):
            page, off, ln = pf.variant_index_page(idx, name, vidx)
            rows = pt.variant_rows(variants[name], vidx=vidx, n=1)
            assert rows["vidx"][0].as_py() == vidx
            assert (off, ln) == (pages["off"][page].as_py(), pages["len"][page].as_py())
            frame, foff, fln = pf.variant_index_hits(idx, name, vidx)
            assert (foff, fln) == (fr["hits_off"][frame].as_py(), fr["hits_len"][frame].as_py())
            # the position lookup lands on the page that holds this variant
            pos = rows["position"][0].as_py()
            section = "cis" if vidx < w["n_cis"] else "trans"
            assert pf.variant_index_position(idx, name, pos, section)[0] == page
    # the rsID samples are the blocks' first records, and a lookup through the startup file agrees
    recs = pt.rsid_records(qbr)
    assert idx["rsid_first"].tolist() == recs["rs_number"].to_pylist()[::B]
    for j in (0, B, recs.num_rows - 1):
        rs = recs["rs_number"][j].as_py()
        assert pt.rsid_lookup(qbr, rs, variant_index=qbx) == pt.rsid_lookup(qbr, rs)
        got = pt.rsid_lookup(qbr, rs, variant_index=qbx)
        # the record points at a real variant, and the page it names decodes to that rsID
        v = pt.variant_rows(variants[got["chr"]], vidx=got["vidx"], n=1)
        assert v["rs_number"][0].as_py() == rs
    c = pt.check_file(qbx)
    assert c["chromosomes"] == list(pf.VARIANT_CHROMS) and c["rsid_records"] == rt.num_rows
    assert c["variants"] == sum(w["total"] for w in want.values())
    assert c["pages"] == sum(pt.page_table(variants[n]).num_rows for n in pf.VARIANT_CHROMS)
    assert c["frames"] == sum(pt.read_header(hits[n])["count"] for n in pf.VARIANT_CHROMS)
    # a hits pack checked through the startup file agrees with one checked by walking
    assert pt.check_file(hits["chr1"], variant_index=qbx)["rows"] == pt.check_file(hits["chr1"])["rows"]
    assert pt.check_file(hits["chr1"], variant_index=qbx)["found_through"] == "variant_index"
    # writer rules
    raises("no variants file", pt.write_variant_index, d / "x.qbx", {k: v for k, v in list(variants.items())[:5]}, hits, qbr)
    raises("no hits pack", pt.write_variant_index, d / "x.qbx", variants, {k: v for k, v in list(hits.items())[:5]}, qbr)
    raises("header kind", pt.write_variant_index, d / "x.qbx", variants, hits, variants["chr1"])
    # the CLI builds the same bytes
    rc, out, _ = cli("pack-variant-index", "-o", d / "cli.qbx", "--variants", d / "variants",
                     "--hits", d / "hits", "--rsid-index", qbr, "--level", 3)
    assert rc == 0 and (d / "cli.qbx").read_bytes() == qbx.read_bytes()
    rc, out, _ = cli("variant-index", qbx, "-o", "-")
    assert rc == 0 and out.splitlines()[0].startswith("chr\tn_cis") and len(out.splitlines()) == 24
    rc, out, _ = cli("variant-index", qbx, "--chrom", "chr1")
    assert rc == 0 and len(json.loads(out)["hits_off"]) == idx["chroms"]["chr1"]["n_frames"] + 1


def real_chr21() -> None:
    # every file is found by its logical key under `immutable/`, the way a reader finds it through
    # manifest.json: the names carry content hashes, so nothing here may be a fixed path (SPEC section 3)
    cfg = Config()
    root = cfg.derived
    pub = addressed_files(cfg)
    keys = {"qbv": "variants/chr21", "qbe": "eqtl/chr21", "qbs": "sqtl/chr21", "qbt": "trans/chr21",
            "qbtM": "trans/chrM", "qbg": "gwas/chr21", "idx": "gwas_index", "si": "search_index"}
    files = {k: pub.get(key, root / "missing") for k, key in keys.items()}
    files["manifest"] = root / "manifest.json"
    if not all(p.is_file() for p in files.values()):
        print("  real chr21 smoke test skipped (packs not built)")
        return
    if subprocess.run(["pgrep", "-f", "locus-bench"], capture_output=True).returncode == 0:
        print("  real chr21 smoke test skipped (a locus-bench benchmark is running; re-run later)")
        return
    for k in ("qbv", "qbe", "qbs", "qbt", "qbg"):
        h = pt.read_header(files[k])
        assert h["chrom"] == "chr21" and h["extension_matches"], k
    assert pt.read_header(files["idx"])["chrom"] == "all" and pt.read_header(files["qbtM"])["chrom"] == "chrM"
    dof = pt.dof_from_manifest(files["manifest"], "eqtl")
    lb = pt.list_blocks(files["qbe"], limit=5, gene_ids=True)
    gid = lb["gene_id"][1].as_py()
    off, ln = pt.find_gene_block(files["qbe"], gid, files["si"])
    assert (off, ln) == (lb["blk_off"][1].as_py(), lb["blk_len"][1].as_py())
    assert pt.find_gene_block(files["qbe"], lb["symbol"][1].as_py()) == (off, ln), "scan by symbol agrees with search_index"
    b = pt.read_block(files["qbe"], off, ln, dof)
    t = pt.block_rows(b, files["qbv"])
    assert t.num_rows == b["n_rows"] and t.schema.remove(0).equals(pf.READER_SCHEMA)
    splice = b["details"]["splice"]
    if splice:
        so, sl = pt.find_intron_block(files["qbe"], splice[0]["phenotype_id"], files["si"])
        ib = pt.read_block(files["qbs"], so, sl, pt.dof_from_manifest(files["manifest"], "sqtl"))
        assert pt.block_rows(ib, files["qbv"]).num_rows == ib["n_rows"] >= 1
    # the variants file's two sections
    hv = pt.read_header(files["qbv"])
    c = pt.check_file(files["qbv"])
    assert c["records"] == c["header_count"] == hv["count"] and c["n_cis"] == hv["n_cis"] < hv["count"] and c["trans_only"] == hv["count"] - hv["n_cis"]
    pages = pt.page_table(files["qbv"])
    assert pages.filter(pc.equal(pages["section"], "trans"))["first_vidx"][0].as_py() == hv["n_cis"], "the first trans-only variant starts a page"
    tr = pt.variant_rows(files["qbv"], vidx=hv["n_cis"], n=3)
    assert tr["vidx"][0].as_py() == hv["n_cis"] and tr["ma_samples"].null_count == 3
    # the trans pack, found through search_index and by walking
    info = pt.find_trans_frame(files["qbt"], gid, files["si"])
    f = pt.read_trans_frame(files["qbt"], info["trans_off"], info["trans_len"], dof, pt.dof_from_manifest(files["manifest"], "sqtl"))
    tt = pt.trans_table(f, "chr21", info["gene_id"], info["gene_version"])
    assert tt.num_rows == f["n_e"] + f["n_s"] and tt["phenotype_id"].null_count == 0
    fr = pt.trans_frames(files["qbt"], files["si"])
    walked = pt.trans_frames(files["qbt"], limit=5)
    assert fr["trans_off"].to_pylist()[:5] == walked["trans_off"].to_pylist() and fr["trans_len"].to_pylist()[:5] == walked["trans_len"].to_pylist()
    ct = pt.check_file(files["qbt"], search_index=files["si"])
    cm = pt.check_file(files["qbtM"])
    assert ct["frames"] == fr.num_rows and cm["frames"] == pt.read_header(files["qbtM"])["count"]
    idx = pt.gwas_index_table(files["idx"])
    first = idx.filter(pc.equal(idx["chr"], "chr21"))["first_position"][0].as_py()
    g = pt.gwas_rows(files["qbg"], files["idx"], lo=first, hi=first + 100_000)
    assert g.num_rows >= 1 and g["position"][0].as_py() == first
    print(f"  real chr21: {gid} block {off}+{ln}, {t.num_rows} rows; {len(splice)} introns; trans frame {f['n_e']} eQTL + {f['n_s']} sQTL rows; "
          f"trans pack {ct['frames']} frames, {ct['rows']:,} rows (chrM {cm['frames']} frames); GWAS window {g.num_rows} rows; "
          f"variants {c['n_cis']:,} cis + {c['trans_only']:,} trans-only in {c['pages']} pages")


def main() -> int:
    for fn in CASES:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"{len(CASES)} synthetic cases passed")
    if "--synthetic" not in sys.argv:
        real_chr21()
    return 0


if __name__ == "__main__":
    sys.exit(main())
