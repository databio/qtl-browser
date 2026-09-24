"""Tests for the results builder (pipeline/results.py).

    uv run python -m pipeline.test_results      (or pytest)

Small synthetic contract tables over the test_catalog variant catalog: a gene with a gap in its tested
variants and a credible set, a gene with no nominal rows, and a leafcutter cluster of three introns,
only one of which has nominal rows.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import packfmt_v1 as pf
from . import qtlstore as qs
from . import results as rs
from .test_annotation import build_fixture as build_annotation
from .test_catalog import build as build_catalog

NOMINAL = [  # phenotype_type, phenotype_id, gene_id, chr, pos, ref, alt, beta, se, pvalue
    ("ge", "G1", "ENSG1", "chr1", 1, "A", "C", 0.5, 0.1, 1e-6),
    ("ge", "G1", "ENSG1", "chr1", 3, "G", "A", -0.2, 0.2, 0.3),
    ("leafcutter", "I1", "ENSG1", "chr2", 7, "TTG", "T", -0.4, 0.05, 1e-9),
    ("leafcutter", "I1", "ENSG1", "chr2", 9, "T", "A", 0.1, 0.1, 0.5),
    ("leafcutter", "I1", "ENSG1", "chr2", 9, "T", "C", 0.3, 0.2, 0.2),
]


def write_tables(d: Path, with_new_columns: bool = True) -> Path:
    t = d / "tables"
    cols = list(zip(*NOMINAL))
    for c in ("chr1", "chr2"):
        keep = [i for i, r in enumerate(NOMINAL) if r[3] == c]
        part = pa.table({k: pa.array([cols[j][i] for i in keep], typ) for j, (k, typ) in enumerate(
            [("phenotype_type", pa.string()), ("phenotype_id", pa.string()), ("gene_id", pa.string()),
             ("chr", pa.string()), ("pos", pa.int32()), ("ref", pa.string()), ("alt", pa.string()),
             ("beta", pa.float32()), ("se", pa.float32()), ("pvalue", pa.float64())])})
        (t / "nominal" / f"chr={c}").mkdir(parents=True)
        pq.write_table(part, t / "nominal" / f"chr={c}" / "data.parquet")
    ph = {"phenotype_type": ["ge", "ge", "leafcutter", "leafcutter", "leafcutter"],
          "phenotype_id": ["G1", "G2", "I1", "I2", "I3"],
          "gene_id": ["ENSG1", "ENSG2", "ENSG1", None, None],
          "extra": ["{}", "{}", '{"cluster_id": "clu_1"}', '{"cluster_id": "clu_1"}', '{"cluster_id": "clu_1"}']}
    perm = {"phenotype_type": ["ge", "ge", "leafcutter"], "phenotype_id": ["G1", "G2", "I1"],
            "gene_id": ["ENSG1", "ENSG2", "ENSG1"], "n_variants": pa.array([3, 0, 3], pa.int32()),
            "p_perm": [0.001, 0.9, 0.01], "p_beta": [0.002, 0.8, 0.02],
            "lead_chr": ["chr1", "chr1", "chr2"], "lead_pos": pa.array([1, 2, 7], pa.int32()),
            "lead_ref": ["A", "C", "TTG"], "lead_alt": ["C", "CA", "T"]}
    cs = {"phenotype_type": ["ge", "ge", "leafcutter"], "phenotype_id": ["G1", "G1", "I1"],
          "cs_id": pa.array([1, 1, 1], pa.int16()), "chr": ["chr1", "chr1", "chr2"],
          "pos": pa.array([1, 2, 7], pa.int32()), "ref": ["A", "C", "TTG"], "alt": ["C", "CA", "T"],
          "pip": pa.array([0.9, 0.1, 1.0], pa.float32()), "z": pa.array([5.0, 1.0, 6.0], pa.float32()),
          "cs_size": pa.array([2, 2, 1], pa.int32()), "cs_min_r2": pa.array([0.8, 0.8, 1.0], pa.float32())}
    if with_new_columns:
        ph["phenotype_object_id"] = ["G1", "G2", "clu_1", "clu_1", "clu_1"]
        ph["has_nominal"] = [True, False, True, False, False]
        perm["phenotype_object_id"] = ["G1", "G2", "clu_1"]
        cs["phenotype_object_id"] = ["G1", "G1", "clu_1"]
    for name, cols_ in (("phenotypes", ph), ("permuted", perm), ("credible_sets", cs)):
        pq.write_table(pa.table(cols_), t / f"{name}.parquet")
    (t / "ingestion.json").write_text(json.dumps({
        "phenotype_types": ["ge", "leafcutter"], "dof": {"ge": 435, "leafcutter": 480},
        "significance": {"column": "p_perm", "op": "<", "threshold": 0.05},
        "allele_orientation_source": "test"}))
    return t


def build_all(d: Path, with_new_columns: bool = True):
    st, rg, cdoc = build_catalog(d)
    build_annotation(st)
    doc = rs.build(st, "exp", write_tables(d, with_new_columns), "cat1", "ann", ["chr1", "chr2"])
    st.write_store("t", [])
    return st, rg, doc


def rows_by_id(st, doc):
    return {r["phenotype_id"]: r for r in rs.load_index(st, doc).to_pylist()}


def test_pointer_and_validate():
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d))
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and st.notes == []
        assert doc["catalog"] == "cat1" and doc["annotation"] == "ann"
        assert [r["phenotype_type"] for r in doc["results"]] == ["ge", "leafcutter"]
        assert [r["dof"] for r in doc["results"]] == [435, 480]
        assert set(doc["results"][0]["files"]) == {"chr1", "chr2"} and set(doc["hits"]) == {"chr1", "chr2"}
        assert doc["results"][0]["precision"]["neglog10p_max_error"] > 0
        h = qs.parse_file_header((st.immutable / doc["results"][1]["files"]["chr2"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"]) == (rs.KIND_RESULTS, "chr2", 3)


def test_slope_precision_is_measured():
    """`slope_max_error_over_se` is the largest |slope rebuilt from the stored codes - source beta| / se over
    every row, measured, not a bound; null when the results set has no dof."""
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d))
        r = rows_by_id(st, doc)
        src = {(pid, pos, alt): (b, se) for _, pid, _, _, pos, _, alt, b, se, _ in NOMINAL}
        for res, ids in zip(doc["results"], (["G1"], ["I1"])):
            want, n = 0.0, 0
            for pid in ids:
                blk = rs.read_block(st, doc, r[pid])
                sites = catalog.load_chrom(st, st.load(qs.POINTER_DIRS[0], doc["catalog"]), r[pid]["chr"])
                for i, sl in enumerate(blk["slope"]):
                    v = blk["var_start"] + i
                    key = (pid, int(sites["pos"][v]), sites["alt"][v])
                    if key in src and np.isfinite(sl):
                        b, se = (float(np.float32(x)) for x in src[key])
                        want, n = max(want, abs(sl - b) / se), n + 1
            pr = res["precision"]
            assert pr["slope_rows_compared"] == n and n > 0
            assert abs(pr["slope_max_error_over_se"] - want) < 1e-12, (pr, want)
        buf = (st.immutable / doc["results"][0]["files"]["chr1"]).read_bytes()
        blk = buf[r["G1"]["blk_off"]:r["G1"]["blk_off"] + r["G1"]["blk_len"]]
        err, n = rs.slope_error(blk, np.array([0.5, np.nan, -0.2]), np.array([0.1, np.nan, 0.2]), None)
        assert (err, n) == (0.0, 0)                  # no dof, no slope, nothing measured


TRANS = [  # phenotype_type, phenotype_id, gene_id, chr, pos, ref, alt, beta, se, pvalue
    ("ge", "G1", "ENSG1", "chr2", 2, "T", "G", 0.8, 0.1, 1e-8),         # a trans-only site
    ("ge", "G1", "ENSG1", "chr1", 3, "G", "A", -0.3, 0.05, 1e-6),
    ("ge", "G9", "ENSG9", "chr2", 7, "TTG", "T", 0.5, 0.1, 2e-7),       # G9 has trans rows and nothing else
    ("leafcutter", "I2", None, "chr1", 1, "A", "C", -1.2, 0.2, 3e-9),
]
GWAS = [  # chr, pos, ref, alt, beta, se, af, pvalue, n, rs_number
    ("chr1", 2, "C", "CA", 0.0123, 0.0045, 0.3333, 2.5e-3, 1000, 2),
    ("chr1", 1, "A", "C", -0.5, 0.1, 0.25, 1.234e-9, 1000, -1),
    ("chr1", 3, "G", "A", 0.2, 0.05, 0.9, 0.5, 900, 3),
    ("chr2", 9, "T", "A", 0.0001, 0.0002, 0.0, 1.0, 1000, 5),
]
# the adapter's bin table over the source rows (unrounded; the builder rounds as v0 did), 5 bp bins
GWAS_BINS = [  # chr, bin_start, bin_end, min_p, lead_position, lead_rsid, lead_beta, lead_ea, n_gws, n_variants
    ("chr2", 5, 10, 1.0, 9, "rs5", 0.0001, "A", 0, 1),
    ("chr1", 0, 5, 1.23456e-9, 1, None, -0.50049, "C", 1, 4),     # 4: a source row the reference did not read
]
BINS_TYPES = [("chr", pa.string()), ("bin_start", pa.int64()), ("bin_end", pa.int64()), ("min_p", pa.float64()),
              ("lead_position", pa.int64()), ("lead_rsid", pa.string()), ("lead_beta", pa.float64()),
              ("lead_ea", pa.string()), ("n_gws", pa.int32()), ("n_variants", pa.int32())]


def write_bins(t: Path, rows: list) -> None:
    cols = list(zip(*rows))
    pq.write_table(pa.table({k: pa.array(list(cols[j]), typ) for j, (k, typ) in enumerate(BINS_TYPES)}),
                   t / "gwas_bins.parquet")


def add_trans_and_gwas(t: Path) -> None:
    """Contract `trans` and `gwas` tables for the synthetic experiment, plus the trans-only phenotype G9."""
    cols = list(zip(*TRANS))
    pq.write_table(pa.table({k: pa.array(list(cols[j]), typ) for j, (k, typ) in enumerate(
        [("phenotype_type", pa.string()), ("phenotype_id", pa.string()), ("gene_id", pa.string()), ("chr", pa.string()),
         ("pos", pa.int32()), ("ref", pa.string()), ("alt", pa.string()), ("beta", pa.float32()), ("se", pa.float32()),
         ("pvalue", pa.float64())])}), t / "trans.parquet")
    ph = pq.read_table(t / "phenotypes.parquet").to_pylist()
    ph.append({"phenotype_type": "ge", "phenotype_id": "G9", "gene_id": "ENSG9", "extra": "{}",
               "phenotype_object_id": "G9", "has_nominal": False})
    pq.write_table(pa.Table.from_pylist(ph), t / "phenotypes.parquet")
    cols = list(zip(*GWAS))
    pq.write_table(pa.table({k: pa.array(list(cols[j]), typ) for j, (k, typ) in enumerate(
        [("chr", pa.string()), ("pos", pa.int32()), ("ref", pa.string()), ("alt", pa.string()), ("beta", pa.float64()),
         ("se", pa.float64()), ("af", pa.float64()), ("pvalue", pa.float64()), ("n", pa.int64()),
         ("rs_number", pa.int64())])}), t / "gwas.parquet")
    (t / "gwas.json").write_text(json.dumps({"id": "test_gwas", "title": "a test GWAS",
                                             "orientation": {"as_is": 3, "swapped": 1, "dropped": 0}}))
    write_bins(t, GWAS_BINS)


def build_with_trans(d: Path):
    st, rg, cdoc = build_catalog(d)
    build_annotation(st)
    t = write_tables(d)
    add_trans_and_gwas(t)
    doc = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
    st.write_store("t", [])
    return st, rg, doc


def test_trans_objects_index_and_hits():
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_with_trans(Path(d))
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and st.notes == []
        r = rows_by_id(st, doc)
        # G9: trans only, so an index row with no chromosome and no block
        assert r["G9"]["chr"] is None and r["G9"]["blk_off"] is None and r["G9"]["n_trans"] == 1
        assert r["G2"]["trans_off"] is None and r["I1"]["trans_off"] is None
        ge = doc["results"][0]["trans"]
        assert (ge["n_rows"], ge["n_phenotypes"]) == (3, 2) and doc["results"][1]["trans"]["n_rows"] == 1
        h = qs.parse_file_header((st.immutable / ge["file"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"]) == (rs.KIND_TRANS, "all", 2)
        assert doc["trans"]["rows"] == 4 and doc["trans"]["rows_skipped_variant_outside_build"] == 0
        g1 = rs.read_trans(st, doc, r["G1"])
        # sorted by (chromosome ordinal, pos): chr1:3 first, then the chr2 trans-only site; no vidx stored
        assert g1["chr"] == ["chr1", "chr2"] and g1["pos"].tolist() == [3, 2] and "vidx" not in g1
        assert g1["ref"] == ["G", "T"] and g1["alt"] == ["A", "G"] and g1["rs_number"].tolist() == [3, 1]
        np.testing.assert_allclose(g1["beta"], [-0.3, 0.8], atol=0.8 / 65534 + 1e-7)
        np.testing.assert_allclose(-np.log10(g1["p"]), [6, 8], atol=8 / 131066 + 1e-9)
        assert np.all(g1["se"] > 0)
        i2 = rs.read_trans(st, doc, r["I2"])
        assert i2["alt"] == ["C"] and abs(i2["beta"][0] + 1.2) < 1e-6 and r["I2"]["chr"] == "chr2"
        g9 = rs.read_trans(st, doc, r["G9"])
        assert g9["ref"] == ["TTG"] and g9["alt"] == ["T"]
        # variant-keyed: kind 2 records in the hits file of the variant's chromosome
        h2 = rs.decode_hits((st.immutable / doc["hits"]["chr2"]).read_bytes())
        tr = h2[h2["kind"] == rs.HIT_TRANS]
        assert [(int(x["vidx"]), int(x["ord"])) for x in tr] == [(0, r["G9"]["ord"]), (3, r["G1"]["ord"])]
        assert abs(float(tr["value"][1]) - 8.0) < 1e-5 and abs(float(tr["beta"][1]) - 0.8) < 1e-6
        assert np.all(np.isnan(h2[h2["kind"] != rs.HIT_TRANS]["beta"]))


def test_index_parts_and_counts():
    """The browser's copy of the search index: one part per chromosome built, a trans-only part, the
    ord runs each part holds, and the counts the Home page prints. `add_split` on a pointer without
    them writes the same keys `build` did, and validate fails a part that differs from the index."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_with_trans(Path(d))
        r = rows_by_id(st, doc)
        parts = doc["search_index_parts"]
        assert list(parts) == ["chr1", "chr2"]
        got = {c: [x["phenotype_id"] for x in rs.decode_arrow((st.immutable / e["file"]).read_bytes()).to_pylist()]
               for c, e in parts.items()}
        assert got == {"chr1": ["G1", "G2"], "chr2": ["I1", "I2", "I3"]}, got
        assert parts["chr1"]["ords"] == [[r["G1"]["ord"], r["G2"]["ord"]]] and parts["chr2"]["rows"] == 3
        to = doc["search_index_trans_only"]
        assert to["rows"] == 1 and to["ords"] == [[r["G9"]["ord"]] * 2]
        # G2 has no rows; the clu_1 group (p_perm 0.01) makes all three introns significant; I2, I3 name no gene
        assert doc["counts"] == {"ge": {"phenotypes": 2, "with_rows": 1, "significant": 1, "significant_genes": 1},
                                 "leafcutter": {"phenotypes": 3, "with_rows": 1, "significant": 3, "significant_genes": 1}}, doc["counts"]
        assert rs.ord_ranges([0, 1, 2, 5, 7, 8]) == [[0, 2], [5, 5], [7, 8]] and rs.ord_ranges([]) == []
        old = {k: v for k, v in doc.items() if k not in rs.SPLIT_KEYS}
        st.write_pointer("experiments", "exp", old)
        assert any("no search_index_parts" in f for f in st.validate(refget=rg))
        assert rs.add_split(st, "exp") == doc
        assert st.validate(refget=rg) == []
        swapped = dict(doc, search_index_parts={"chr1": parts["chr2"], "chr2": parts["chr1"]})
        st.write_pointer("experiments", "exp", swapped)
        fails = st.validate(refget=rg)
        assert any("part chr1: rows differ" in f for f in fails), fails


def _layout_rows():
    """Search-index rows for the frame-order rule: gene A (an eQTL and two introns), gene B (introns only,
    interleaved with A's by position, so by `ord`), gene C (an eQTL and one intron) and an intron with no
    gene, between them."""
    spec = [("ge", "gA", "A"), ("ge", "gC", "C"), ("leafcutter", "iB1", "B"), ("leafcutter", "iA1", "A"),
            ("leafcutter", "iN", None), ("leafcutter", "iB2", "B"), ("leafcutter", "iA2", "A"), ("leafcutter", "iC1", "C")]
    return [{"ord": i, "phenotype_type": t, "phenotype_id": p, "gene_id": g, "trans_off": None, "trans_len": None}
            for i, (t, p, g) in enumerate(spec)]


def test_trans_frame_order_keeps_genes_contiguous():
    """Frames are grouped by gene (genes in search-index order: the gene's eQTL `ord`, else its smallest),
    `ord` order within a gene, gene-less phenotypes last; `check_trans_layout` passes that layout and
    catches the old `ord` order, where gene A's introns straddle gene B's."""
    rows = _layout_rows()
    lc = [r["ord"] for r in rows if r["phenotype_type"] == "leafcutter"]
    order = rs.trans_frame_order(rows, lc)
    assert [rows[k]["phenotype_id"] for k in order] == ["iA1", "iA2", "iC1", "iB1", "iB2", "iN"]
    assert rs.trans_frame_order(rows, [1, 0]) == [0, 1]

    def lay(order_by_type):
        rr = [dict(r) for r in rows]
        sizes, doc = {}, {"results": []}
        for t, order in order_by_type.items():
            off = qs.HEADER_LEN
            for k in order:
                rr[k]["trans_off"], rr[k]["trans_len"] = off, 10 + k
                off += 10 + k
            sizes[f"{t}.qbt"] = off
            doc["results"].append({"phenotype_type": t, "trans": {"file": f"{t}.qbt", "n_phenotypes": len(order)}})
        return rr, doc, sizes
    rr, doc, sizes = lay({"ge": rs.trans_frame_order(rows, [0, 1]), "leafcutter": order})
    assert rs.check_trans_layout(rr, doc, sizes) == []
    # a gene page reads one range per object: A's two introns, 10 + 3 and 10 + 6 bytes
    assert rs.gene_trans_ranges(rr, "A") == {"ge": (64, 10), "leafcutter": (64, 29)}
    rr, doc, sizes = lay({"ge": [0, 1], "leafcutter": lc})
    fails = rs.check_trans_layout(rr, doc, sizes)
    assert len(fails) == 1 and "2 genes whose frames are not one contiguous range" in fails[0], fails
    rr[3]["trans_off"] += 1                                  # a gap between frames
    assert any("previous frame ends" in f for f in rs.check_trans_layout(rr, doc, sizes))


def test_trans_rejects_orphans():
    with tempfile.TemporaryDirectory() as d:
        st, rg, cdoc = build_catalog(Path(d))
        build_annotation(st)
        t = write_tables(Path(d))
        add_trans_and_gwas(t)
        tr = pq.read_table(t / "trans.parquet").to_pylist()
        pq.write_table(pa.Table.from_pylist(tr + [{**tr[0], "pos": 5}]), t / "trans.parquet")
        try:
            rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
            raise AssertionError("an orphan trans row was accepted")
        except ValueError as e:
            assert "not in the variant catalog" in str(e), e
        pq.write_table(pa.Table.from_pylist(tr + [{**tr[0], "phenotype_id": "NOPE"}]), t / "trans.parquet")
        try:
            rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
            raise AssertionError("a trans row naming an unknown phenotype was accepted")
        except ValueError as e:
            assert "not in `phenotypes`" in str(e), e


def test_gwas_object_and_window():
    from . import gwas
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_with_trans(Path(d))
        g = doc["gwas"]
        assert g["id"] == "test_gwas" and g["n_rows"] == 4 and set(g["files"]) == {"chr1", "chr2"}
        assert g["n_values"] == [900, 1000] and g["orientation"]["swapped"] == 1
        h = qs.parse_file_header((st.immutable / g["files"]["chr1"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (gwas.KIND_GWAS, "chr1", 3, gwas.BLOCK_ROWS)
        w = gwas.read_window(st, doc, "chr1", 1, 2)
        assert w["pos"].tolist() == [1, 2] and w["ref"] == ["A", "C"] and w["alt"] == ["C", "CA"]
        assert w["beta"].tolist() == [-0.5, 0.0123] and w["af"].tolist() == [0.25, 0.3333]
        assert w["p"].tolist() == [1.234e-9, 2.5e-3] and w["rs_number"].tolist() == [0, 2] and w["n"].tolist() == [1000, 1000]
        assert gwas.read_window(st, doc, "chr2", 10, 20)["pos"].size == 0
        assert g["rows_sharing_a_site"] == 0
        t = Path(d) / "tables"
        gw = pq.read_table(t / "gwas.parquet").to_pylist()
        # a site reported twice (the source's two allele orders of one indel): both rows kept, counted
        pq.write_table(pa.Table.from_pylist(gw + [{**gw[1], "beta": 0.3, "n": 900, "pvalue": 0.04}]), t / "gwas.parquet")
        g2 = gwas.build(st, t, st.load(qs.CATALOGS, "cat1"), ["chr1", "chr2"])
        assert g2["rows_sharing_a_site"] == 2 and g2["n_rows"] == 5
        # lossless means lossless: more than 4 decimals is refused, not rounded
        gw[0]["beta"] = 0.12345
        pq.write_table(pa.Table.from_pylist(gw), t / "gwas.parquet")
        try:
            gwas.build(st, t, st.load(qs.CATALOGS, "cat1"), ["chr1", "chr2"])
            raise AssertionError("a 5-decimal beta was accepted")
        except ValueError as e:
            assert "more than 4 decimals" in str(e), e


def test_gwas_bins():
    """The bin summary: v0's values and rounding, catalog table order, named from `gwas.bins`, checked by
    validate; bad bins refused at build and caught by validate."""
    from . import catalog, gwas
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_with_trans(Path(d))
        b = doc["gwas"]["bins"]
        assert (b["bin_bp"], b["n_bins"]) == (5, 2) and b["file"].endswith(".arrow.zst")
        assert b["file"] in qs.object_names(doc)
        t = gwas.read_bins(st, doc)
        assert t.schema == gwas.BINS_SCHEMA
        # chr1 before chr2 (the variant catalog's order); p to 3 significant digits, beta to 3 decimals
        assert t.to_pylist() == [
            {"chr": "chr1", "bin_start": 0, "bin_end": 5, "min_p": 1.23e-9, "lead_position": 1, "lead_rsid": None,
             "lead_beta": -0.5, "lead_ea": "C", "n_gws": 1, "n_variants": 4},
            {"chr": "chr2", "bin_start": 5, "bin_end": 10, "min_p": 1.0, "lead_position": 9, "lead_rsid": "rs5",
             "lead_beta": 0.0, "lead_ea": "A", "n_gws": 0, "n_variants": 1}]
        assert st.validate(refget=rg, sites=catalog.read_sites) == []
        t_dir = Path(d) / "tables"
        cdoc = st.load(qs.CATALOGS, "cat1")
        # a subset build keeps the bins of its chromosomes only
        assert gwas.build(st, t_dir, cdoc, ["chr1"])["bins"]["n_bins"] == 1
        # no bin table: no bin summary
        (t_dir / "gwas_bins.parquet").unlink()
        assert gwas.build(st, t_dir, cdoc, ["chr1", "chr2"])["bins"] is None
        # the builder refuses a bin that is not aligned, a lead outside its bin, more hits than rows
        for bad, rule in (({"bin_start": 1, "bin_end": 6}, "not [k * 5"), ({"lead_position": 12}, "outside the bin"),
                          ({"n_gws": 2}, "n_gws 2, n_variants 1")):
            row = dict(zip([k for k, _ in BINS_TYPES], GWAS_BINS[0]))
            write_bins(t_dir, [GWAS_BINS[1], tuple({**row, **bad}.values())])
            try:
                gwas.build(st, t_dir, cdoc, ["chr1", "chr2"])
                raise AssertionError(f"a bad bin was accepted ({bad})")
            except ValueError as e:
                assert rule in str(e), e
        # validate: a pointer naming a bin summary whose count, bins or chromosomes are wrong
        good = doc["gwas"]["bins"]
        for change, rule in (({"n_bins": 3}, "the pointer says 3"), ({"bin_bp": 10}, "not [k * 10")):
            st.write_pointer("experiments", "exp", {**doc, "gwas": {**doc["gwas"], "bins": {**good, **change}}})
            fails = st.validate(refget=rg, sites=catalog.read_sites)
            assert any("gwas bins" in f and rule in f for f in fails), fails
        write_bins(t_dir, GWAS_BINS)
        only1 = gwas.build(st, t_dir, cdoc, ["chr1", "chr2"])
        st.write_pointer("experiments", "exp", {**doc, "gwas": {**only1, "files": {"chr1": only1["files"]["chr1"]}}})
        assert any("have bins and no GWAS file" in f for f in st.validate(refget=rg, sites=catalog.read_sites))
        st.write_pointer("experiments", "exp", doc)
        assert st.validate(refget=rg, sites=catalog.read_sites) == []


def test_hits_paging():
    """Frames of F variants, an offset table after the header, zero-length frames for empty ranges; one
    frame answers one variant."""
    recs = np.concatenate([pf.hit_records([0, 5, 5, 2049], [3, 1, 0, 7], [0.01, 0.5, 0.9, 0.2], rs.HIT_CS, cs_id=1),
                           pf.hit_records([5], [9], [12.5], rs.HIT_TRANS, beta=[-0.25])])
    buf = rs.encode_hits("chr1", "A" * 32, recs, 2100, frame_variants=1024)
    h = qs.parse_file_header(buf)
    assert (h["count"], h["page_size"], h["n_cis"]) == (5, 1024, 2100)
    offs = pf.hits_frame_table(buf, 2100, 1024)
    assert len(offs) == 4 and offs[0] == 64 + 16 and offs[-1] == len(buf)
    assert offs[2] - offs[1] == 0                             # vidx 1024..2047: no records, no bytes
    f5 = rs.hits_frame(buf, 5)
    assert [(int(x["vidx"]), int(x["kind"]), int(x["ord"])) for x in f5] == [(0, 1, 3), (5, 1, 0), (5, 1, 1), (5, 2, 9)]
    assert rs.hits_frame(buf, 1500).size == 0 and rs.hits_frame(buf, 2049)["ord"].tolist() == [7]
    assert len(rs.decode_hits(buf)) == 5
    empty = rs.encode_hits("chr1", "A" * 32, pf.hit_records([], [], [], rs.HIT_LEAD), 0)
    assert len(empty) == 64 + 4 and len(rs.decode_hits(empty)) == 0


def test_remove_experiment_and_gc():
    """Add an experiment, remove it, gc: the first experiment's objects stay byte for byte, the removed one's
    own objects go, the store validates at every step, and a second gc deletes nothing."""
    import shutil
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, doc = build_with_trans(d)
        before = {n: (st.immutable / n).read_bytes() for n in qs.object_names(doc)}
        other = d / "other"
        shutil.copytree(write_tables(d / "o", with_new_columns=False), other)
        doc2 = rs.build(st, "exp2", other, "cat1", "ann", ["chr1", "chr2"])
        st.write_store("t", [])
        assert st.validate(refget=rg, sites=catalog.read_sites) == []
        only2 = set(qs.object_names(doc2)) - set(qs.object_names(doc)) - set(qs.object_names(st.load(qs.CATALOGS, "cat1")))
        assert only2, "exp2 should own some objects"
        assert qs.gc(st, dry_run=True)["would_delete"] == 0            # everything is referenced
        qs.remove_experiment(st, "exp2")
        assert json.loads((st.root / "store.json").read_text())["experiments"] == ["exp"]
        assert st.validate(refget=rg, sites=catalog.read_sites) == []
        r = qs.gc(st)
        assert r["deleted"] == len(only2) and not any((st.immutable / n).exists() for n in only2)
        assert {n: (st.immutable / n).read_bytes() for n in qs.object_names(doc)} == before
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and qs.gc(st)["deleted"] == 0
        try:
            qs.remove_experiment(st, "exp2")
            raise AssertionError("removing a missing experiment succeeded")
        except ValueError as e:
            assert "does not exist" in str(e)


def test_blocks_codes_and_gaps():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d))
        r = rows_by_id(st, doc)
        assert [r[k]["ord"] for k in ("G1", "G2", "I1", "I2", "I3")] == [0, 1, 2, 3, 4]
        g1 = rs.read_block(st, doc, r["G1"])
        # G1 tests vidx 0 and 2 (chr1 order: 1 A>C, 2 C>CA, 3 G>A); vidx 1 is only in its credible set
        assert (g1["var_start"], g1["n_rows"], g1["anchor"]) == (0, 3, 0)
        nlp_q, _ = pf.quantize_nlp([1e-6, 0.3])
        se_q, _, _ = pf.quantize_se([0.1, 0.2], [0.5, -0.2])
        assert g1["nlp_code"].tolist() == [nlp_q[0], pf.NLP_NULL, nlp_q[1]]
        assert g1["se_code"].tolist() == [se_q[0], pf.SE_NULL, se_q[1]] and g1["se_code"][2] & pf.SE_SIGN
        assert g1["cs_row"].tolist() == [0, 1] and g1["cs_id"].tolist() == [1, 1]
        det = g1["details"]
        assert det["v"] == 1 and det["has_nominal"] and det["n_nominal"] == 2 and det["n_credible_sets"] == 1
        assert det["group"]["p_perm"] == 0.001 and det["group"]["significant"]
        assert not {"symbol", "tss", "biotype", "start", "end", "strand"} & set(det)
        assert (r["G1"]["w_lo"], r["G1"]["w_hi"], r["G1"]["n_var"]) == (1, 3, 3)


def test_no_nominal_and_groups():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d))
        r = rows_by_id(st, doc)
        g2 = rs.read_block(st, doc, r["G2"])
        assert g2["n_rows"] == 0 and g2["details"]["has_nominal"] is False
        assert not r["G2"]["has_nominal"] and r["G2"]["var_start"] is None and not r["G2"]["significant"]
        for k in ("I2", "I3"):                            # non-lead introns: located by the cluster's lead
            assert r[k]["chr"] == "chr2" and not r[k]["has_nominal"] and not r[k]["is_group_lead"]
            b = rs.read_block(st, doc, r[k])
            assert b["n_rows"] == 0 and b["details"]["group"]["lead_phenotype_id"] == "I1"
            assert b["details"]["phenotype_object_id"] == "clu_1" and b["details"]["extra"] == {"cluster_id": "clu_1"}
        i1 = rs.read_block(st, doc, r["I1"])
        assert r["I1"]["is_group_lead"] and r["I1"]["significant"] and i1["n_rows"] == 3


def test_hits():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d))
        h1 = rs.decode_hits((st.immutable / doc["hits"]["chr1"]).read_bytes())
        got = [(int(x["vidx"]), int(x["kind"]), int(x["ord"]), int(x["cs_id"]), int(x["flags"])) for x in h1]
        assert got == [(0, 0, 0, 0, 1), (0, 1, 0, 1, 0), (1, 0, 1, 0, 0), (1, 1, 0, 1, 0)]
        assert np.isclose(h1["value"][1], 0.9)
        h2 = rs.decode_hits((st.immutable / doc["hits"]["chr2"]).read_bytes())
        assert [(int(x["vidx"]), int(x["kind"]), int(x["ord"])) for x in h2] == [(0, 0, 2), (0, 1, 2)]


def test_significance_uses_the_rules_column():
    """A rule on `p_beta` tests `p_beta`, not `p_perm`. With p_beta < 0.015, G1 (p_perm 0.001, p_beta
    0.002) passes and I1 (p_perm 0.01, p_beta 0.02) does not, though I1's p_perm is under 0.015."""
    rule = {"column": "p_beta", "op": "<", "threshold": 0.015}
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        build_annotation(st)
        t = write_tables(d)
        ing = json.loads((t / "ingestion.json").read_text())
        doc = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], ingestion={**ing, "significance": rule})
        assert doc["significance"] == rule
        r = rows_by_id(st, doc)
        assert r["G1"]["significant"] and not r["I1"]["significant"] and not r["I2"]["significant"]
        assert r["I1"]["p_perm"] == 0.01                       # the index still reports p_perm
        assert not rs.read_block(st, doc, r["I1"])["details"]["group"]["significant"]
        assert rs.read_block(st, doc, r["G1"])["details"]["group"]["significant"]
        h2 = rs.decode_hits((st.immutable / doc["hits"]["chr2"]).read_bytes())
        assert [int(x["flags"]) for x in h2 if x["kind"] == 0] == [0]
        for bad in ({"column": "pvalue", "op": "<", "threshold": 0.05}, {"column": "p_beta", "op": ">", "threshold": 1}):
            try:
                rs.build(st, "exp2", t, "cat1", "ann", ["chr1", "chr2"], ingestion={**ing, "significance": bad})
            except ValueError as e:
                assert "significance rule" in str(e), e
            else:
                raise AssertionError(f"no error for {bad}")


def test_defaults_for_old_tables():
    """Tables written before the contract added phenotype_object_id and has_nominal still build,
    with phenotype_object_id = phenotype_id and has_nominal from the nominal rows."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build_all(Path(d), with_new_columns=False)
        r = rows_by_id(st, doc)
        assert r["I1"]["phenotype_object_id"] == "I1" and r["I1"]["has_nominal"] and r["G1"]["has_nominal"]
        # without the group column I2 and I3 are their own groups, with no permuted row: unplaceable,
        # and counted rather than dropped
        assert "I2" not in r and doc["unplaced"] == {"count": 2, "examples": ["I2", "I3"]}


def test_orphan_nominal_row_fails():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        t = write_tables(d)
        p = t / "nominal" / "chr=chr1" / "data.parquet"
        tab = pq.read_table(p, partitioning=None)   # the file alone: no chr= partition column
        pq.write_table(tab.set_column(4, "pos", pa.array([1, 30], pa.int32())), p)
        try:
            rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
        except ValueError as e:
            assert "not in the variant catalog" in str(e), e
        else:
            raise AssertionError("orphan nominal row was accepted")


def test_fitted_dof_and_has_nominal_check():
    assert rs._dof(435) == (435, None)
    d, fit = rs._dof({"dof": 367, "usable": True, "dof_for_manifest": 367, "residual_log10p": 1e-6})
    assert d == 367 and fit["residual_log10p"] == 1e-6
    assert rs._dof({"dof": 12, "usable": False, "reason": "no fit"})[0] is None
    with tempfile.TemporaryDirectory() as d:
        st, rg, cdoc = build_catalog(Path(d))
        build_annotation(st)
        t = write_tables(Path(d))
        ph = pq.read_table(t / "phenotypes.parquet").to_pandas()
        ph["has_nominal"] = True
        pq.write_table(pa.Table.from_pandas(ph, preserve_index=False), t / "phenotypes.parquet")
        try:
            rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
        except ValueError as e:
            assert "has_nominal" in str(e)
        else:
            raise AssertionError("a wrong has_nominal flag was accepted")


def test_read_block_with_null_dof():
    # an unusable dof fit stores `dof: null`; read_block must not invent slopes from a stand-in dof
    with tempfile.TemporaryDirectory() as d:
        st, rg, cdoc = build_catalog(Path(d))
        build_annotation(st)
        t = write_tables(Path(d))
        ing = json.loads((t / "ingestion.json").read_text())
        ing["dof"]["ge"] = {"dof": 12, "usable": False, "reason": "no fit"}
        (t / "ingestion.json").write_text(json.dumps(ing))
        doc = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
        assert [r["dof"] for r in doc["results"]] == [None, 480]
        r = rows_by_id(st, doc)
        g1 = rs.read_block(st, doc, r["G1"])
        assert g1["n_rows"] == 3 and np.isnan(g1["slope"]).all()
        assert np.isfinite(g1["nlp"][[0, 2]]).all() and g1["se_code"][2] & pf.SE_SIGN   # p, SE, sign survive
        i1 = rs.read_block(st, doc, r["I1"])                                             # dof 480: slopes rebuilt
        assert np.isfinite(i1["slope"]).any()


def test_duplicate_credible_set_rows():
    with tempfile.TemporaryDirectory() as d:
        st, rg, cdoc = build_catalog(Path(d))
        build_annotation(st)
        t = write_tables(Path(d))
        cs = pq.read_table(t / "credible_sets.parquet")
        # the adapter removes a source's repeats; the builder accepts none, identical or not
        bad = cs.slice(0, 1).to_pandas()
        bad["pip"] = 0.5
        for extra in (cs.slice(0, 1), pa.Table.from_pandas(bad, schema=cs.schema, preserve_index=False)):
            pq.write_table(pa.concat_tables([cs, extra]), t / "credible_sets.parquet")
            try:
                rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
            except ValueError as e:
                assert "repeat a (phenotype, cs_id, site)" in str(e)
            else:
                raise AssertionError("repeated credible-set rows were accepted")


def main() -> int:
    bad = 0
    tests = [v for k, v in globals().items() if k.startswith("test_")]
    for fn in tests:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001
            import traceback
            traceback.print_exc()
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(tests) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
