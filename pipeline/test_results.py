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


def build_all(d: Path, with_new_columns: bool = True, version: int = qs.FORMAT_VERSION):
    """A whole build at one format version. The default is the store's default; a test of v1's
    own semantics (the rebuilt slope, the lossless GWAS) names `version=1`."""
    st, rg, cdoc = build_catalog(d, version=version)
    build_annotation(st, version=version)
    docs = rs.build(st, "exp", write_tables(d, with_new_columns), "cat1", "ann", ["chr1", "chr2"],
                    codec=qs.codec_for(version))
    st.write_store("t", [], version)
    return st, rg, docs


def read_row(st, docs, row):
    """`read_block` against whichever cis scan owns this row. A search-index row belongs to an
    experiment now, so a test holding a merged view has to say which."""
    return rs.read_block(st, one(docs, "qtl_cis", row["phenotype_type"]), row)


def read_row_trans(st, docs, row):
    return rs.read_trans(st, one(docs, "qtl_trans", row["phenotype_type"]), row)


def rows_by_id(st, doc):
    """`phenotype_id` -> its search-index row.

    Accepts one experiment, or a whole build — in which case every cis scan's rows are merged, which
    is what a test asking about several phenotype types means. Note `ord` is per experiment now, so
    a merged view has repeated ordinals by design."""
    if "kind" not in doc:                     # a {id: pointer} build
        out = {}
        for d in doc.values():
            if d["kind"] == "qtl_cis":
                out |= {r["phenotype_id"]: r for r in rs.load_index(st, d).to_pylist()}
        return out
    return {r["phenotype_id"]: r for r in rs.load_index(st, doc).to_pylist()}


def test_pointer_and_validate():
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d))
        doc = one(docs)
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and st.notes == []
        assert doc["catalog"] == "cat1" and doc["annotation"] == "ann"
        # one experiment per scan: two phenotype types, cis and trans each, and the GWAS
        # one experiment per scan; these tables have no trans or GWAS, so two cis scans
        assert sorted(docs) == ["exp_ge_cis", "exp_leafcutter_cis"]
        assert {d["kind"] for d in docs.values()} == {"qtl_cis"}
        # they all share one collection digest, which is what makes them one collection
        assert len({d["collection"] for d in docs.values()}) == 1
        assert [one(docs, "qtl_cis", p)["dof"] for p in ("ge", "leafcutter")] == [435, 480]
        assert set(doc["files"]) == {"chr1", "chr2"} and set(doc["hits"]) == {"chr1", "chr2"}
        assert doc["precision"]["theoretical_limit"]["neglog10p"] > 0
        lc = one(docs, "qtl_cis", "leafcutter")
        h = qs.parse_file_header((st.immutable / lc["files"]["chr2"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"]) == (rs.KIND_RESULTS, "chr2", 3)


def test_slope_precision_is_measured():
    """`measured_worst.beta_over_se` is the largest |slope rebuilt from the stored codes - source
    beta| / se over every row, measured, not a limit; null when the results set has no dof.

    Under v1 this is the only quantity that can be measured at all: `nlp` and `se` are what v1
    stores, so there is nothing to compare them against, and it has no `beta_max` from which to
    derive a ceiling for the effect size. `theoretical_limit.beta_over_se` is therefore 0.0 here --
    which is why v1's budget breach went unnoticed for so long."""
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d), version=1)
        doc = one(docs)
        r = rows_by_id(st, docs)
        src = {(pid, pos, alt): (b, se) for _, pid, _, _, pos, _, alt, b, se, _ in NOMINAL}
        for res, ids in ((one(docs, "qtl_cis", "ge"), ["G1"]), (one(docs, "qtl_cis", "leafcutter"), ["I1"])):
            want, n = 0.0, 0
            for pid in ids:
                blk = read_row(st, docs, r[pid])
                sites = catalog.load_chrom(st, st.load(qs.POINTER_DIRS[0], doc["catalog"]), r[pid]["chr"])
                for i, sl in enumerate(blk["slope"]):
                    v = blk["var_start"] + i
                    key = (pid, int(sites["pos"][v]), sites["alt"][v])
                    if key in src and np.isfinite(sl):
                        b, se = (float(np.float32(x)) for x in src[key])
                        want, n = max(want, abs(sl - b) / se), n + 1
            pr = res["precision"]
            assert pr["measured_worst"]["rows_compared"] == n and n > 0
            assert abs(pr["measured_worst"]["beta_over_se"] - want) < 1e-12, (pr, want)
            assert pr["theoretical_limit"]["beta_over_se"] == 0.0        # v1 stores no beta_max
            assert pr["theoretical_limit"]["neglog10p"] > 0 and pr["theoretical_limit"]["se_rel"] > 0
        buf = (st.immutable / doc["files"]["chr1"]).read_bytes()
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


def one(docs: dict, kind: str = "qtl_cis", ptype: str = "ge") -> dict:
    """One experiment out of a build's several. `build` returns `{id: pointer}` now -- a cohort with
    two phenotype types, cis and trans, plus a GWAS is five experiments -- so a test that means "the
    eQTL cis scan" has to say so."""
    if kind == "gwas":
        return next(d for d in docs.values() if d["kind"] == "gwas")
    return next(d for d in docs.values() if d["kind"] == kind and d.get("phenotype_type") == ptype)


def build_with_trans(d: Path, version: int = qs.FORMAT_VERSION):
    st, rg, cdoc = build_catalog(d, version=version)
    build_annotation(st, version=version)
    t = write_tables(d)
    add_trans_and_gwas(t)
    docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], codec=qs.codec_for(version))
    st.write_store("t", [], version)
    return st, rg, docs


def test_trans_objects_index_and_hits():
    from . import catalog
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_with_trans(Path(d))
        doc = one(docs)
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and st.notes == []
        r = rows_by_id(st, docs)
        ge, lc = one(docs, "qtl_trans", "ge"), one(docs, "qtl_trans", "leafcutter")
        tr_rows = {d["phenotype_type"]: {x["phenotype_id"]: x for x in rs.load_index(st, d).to_pylist()}
                   for d in (ge, lc)}
        # Membership, not null columns, now says which scan a phenotype is in. G9 has trans rows and
        # no cis result, so it is in the trans index and absent from the cis one -- which is what
        # `search_index_trans_only` used to hold, and why it no longer needs to exist.
        assert "G9" not in r and tr_rows["ge"]["G9"]["n_trans"] == 1
        assert "G2" not in tr_rows["ge"] and "I1" not in tr_rows["leafcutter"]
        assert set(tr_rows["ge"]) == {"G1", "G9"} and set(tr_rows["leafcutter"]) == {"I2"}
        assert (ge["n_rows"], ge["n_phenotypes"]) == (3, 2) and lc["n_rows"] == 1
        h = qs.parse_file_header((st.immutable / ge["file"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"]) == (rs.KIND_TRANS, "all", 2)
        # each trans experiment reports its own type's rows, not the cohort's, and carries no
        # permutation rule, since that is a cis quantity
        assert ge["placement"] == {"rows": 3, "rows_in_source": 3, "rows_skipped_variant_outside_build": 0,
                                   "rows_skipped_phenotype_outside_build": 0}, ge["placement"]
        assert lc["placement"]["rows"] == lc["placement"]["rows_in_source"] == 1
        assert "significance" not in ge and "significance" not in lc and ge["excluded"] is None
        r = {**r, **tr_rows["ge"], **tr_rows["leafcutter"]}      # both scopes, for the reads below
        g1 = read_row_trans(st, docs, r["G1"])
        # sorted by (chromosome ordinal, pos): chr1:3 first, then the chr2 trans-only site; no vidx stored
        assert g1["chr"] == ["chr1", "chr2"] and g1["pos"].tolist() == [3, 2] and "vidx" not in g1
        assert g1["ref"] == ["G", "T"] and g1["alt"] == ["A", "G"] and g1["rs_number"].tolist() == [3, 1]
        np.testing.assert_allclose(g1["beta"], [-0.3, 0.8], atol=0.8 / 65534 + 1e-7)
        np.testing.assert_allclose(-np.log10(g1["p"]), [6, 8], atol=8 / 131066 + 1e-9)
        assert np.all(g1["se"] > 0)
        i2 = read_row_trans(st, docs, r["I2"])
        assert i2["alt"] == ["C"] and abs(i2["beta"][0] + 1.2) < 1e-6 and r["I2"]["chr"] == "chr2"
        g9 = read_row_trans(st, docs, r["G9"])
        assert g9["ref"] == ["TTG"] and g9["alt"] == ["T"]
        # Variant-keyed, and now split by scan: the trans experiment's hits hold kind-2 records only,
        # the cis experiment's hold leads and credible sets only. `ord` is in each experiment's own
        # index, which is why the expected ordinals are read from the trans index rather than written.
        h2 = rs.decode_hits((st.immutable / ge["hits"]["chr2"]).read_bytes())
        assert set(h2["kind"].tolist()) == {rs.HIT_TRANS}
        assert [(int(x["vidx"]), int(x["ord"])) for x in h2] == [(0, tr_rows["ge"]["G9"]["ord"]),
                                                                 (3, tr_rows["ge"]["G1"]["ord"])]
        assert abs(float(h2["value"][1]) - 8.0) < 1e-5 and abs(float(h2["beta"][1]) - 0.8) < 1e-6
        cis2 = rs.decode_hits((st.immutable / one(docs, "qtl_cis", "leafcutter")["hits"]["chr2"]).read_bytes())
        assert rs.HIT_TRANS not in set(cis2["kind"].tolist())
        assert np.all(np.isnan(cis2["beta"]))        # only a trans record carries a beta
        # Each experiment owns its ordinal space, so every index numbers its own rows 0..n-1. That
        # is the invariant the renumbering exists to hold, and it is why an `ord` read from one
        # experiment's index means nothing in another's.
        for dd in docs.values():
            if dd["kind"] != "gwas":
                ords = [x["ord"] for x in rs.load_index(st, dd).to_pylist()]
                assert ords == list(range(len(ords))), (dd["id"], ords)


def test_index_parts_and_counts():
    """The browser's copy of the search index: one part per chromosome built, a trans-only part, the
    ord runs each part holds, and the counts the Home page prints. `add_split` on a pointer without
    them writes the same keys `build` did, and validate fails a part that differs from the index."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_with_trans(Path(d))
        doc = one(docs)
        r = rows_by_id(st, docs)
        # Parts are per experiment now, so each one holds only its own scan's phenotypes. A
        # chromosome with none of them still gets a part, empty, as it always did.
        lc = one(docs, "qtl_cis", "leafcutter")
        got = {eid: {c: [x["phenotype_id"] for x in
                         rs.decode_arrow((st.immutable / e["file"]).read_bytes()).to_pylist()]
                     for c, e in d["search_index_parts"].items()}
               for eid, d in (("ge", doc), ("leafcutter", lc))}
        assert got == {"ge": {"chr1": ["G1", "G2"], "chr2": []},
                       "leafcutter": {"chr1": [], "chr2": ["I1", "I2", "I3"]}}, got
        assert doc["search_index_parts"]["chr1"]["ords"] == [[r["G1"]["ord"], r["G2"]["ord"]]]
        assert lc["search_index_parts"]["chr2"]["rows"] == 3
        # The trans index is split the same way, by the *phenotype's* chromosome, since that is what
        # the gene page asks under. G1 is on chr1 by its cis result; G9 has no cis result and its
        # gene ENSG9 is not in the annotation, so it is the one row nothing can place, and it goes
        # in `search_index_unplaced` rather than being dropped.
        assert "search_index_trans_only" not in doc
        tr = one(docs, "qtl_trans", "ge")
        assert [x["phenotype_id"] for x in rs.load_index(st, tr).to_pylist()] == ["G1", "G9"]
        tparts = {c: [x["phenotype_id"] for x in
                      rs.decode_arrow((st.immutable / e["file"]).read_bytes()).to_pylist()]
                  for c, e in tr["search_index_parts"].items()}
        assert tparts == {"chr1": ["G1"], "chr2": []}, tparts
        unp = rs.decode_arrow((st.immutable / tr["search_index_unplaced"]["file"]).read_bytes())
        assert [x["phenotype_id"] for x in unp.to_pylist()] == ["G9"] and tr["search_index_unplaced"]["rows"] == 1
        assert "search_index_unplaced" not in doc and "search_index_unplaced" not in lc
        # validate checks the trans parts too: dropping the unplaced part is a failure, not a shrug
        st.write_pointer("experiments", tr["id"], {k: v for k, v in tr.items() if k != "search_index_unplaced"})
        assert any("no search_index_unplaced" in f for f in st.validate(refget=rg)), st.validate(refget=rg)
        st.write_pointer("experiments", tr["id"], tr)
        assert st.validate(refget=rg) == []
        # counts are this scan's own: G2 has no rows, and the clu_1 group (p_perm 0.01) makes all
        # three introns significant while I2 and I3 name no gene
        assert doc["counts"] == {"phenotypes": 2, "with_rows": 1, "significant": 1, "significant_genes": 1}, doc["counts"]
        assert lc["counts"] == {"phenotypes": 3, "with_rows": 1, "significant": 3, "significant_genes": 1}, lc["counts"]
        assert rs.ord_ranges([0, 1, 2, 5, 7, 8]) == [[0, 2], [5, 5], [7, 8]] and rs.ord_ranges([]) == []
        old = {k: v for k, v in doc.items() if k not in rs.SPLIT_KEYS}
        st.write_pointer("experiments", doc["id"], old)
        assert any("no search_index_parts" in f for f in st.validate(refget=rg))
        assert rs.add_split(st, doc["id"]) == doc
        assert st.validate(refget=rg) == []
        p = doc["search_index_parts"]
        swapped = dict(doc, search_index_parts={"chr1": p["chr2"], "chr2": p["chr1"]})
        st.write_pointer("experiments", doc["id"], swapped)
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
        """One trans experiment per phenotype type, as `build` writes them: the layout check now
        takes a single scan and that scan's own rows, every one of which has a frame."""
        rr = [dict(r) for r in rows]
        per = {}
        for t, order in order_by_type.items():
            off = qs.HEADER_LEN
            for k in order:
                rr[k]["trans_off"], rr[k]["trans_len"] = off, 10 + k
                off += 10 + k
            per[t] = ({"kind": "qtl_trans", "phenotype_type": t, "file": f"{t}.qbt"}, {f"{t}.qbt": off})
        return rr, per

    def check(rr, per):
        return [f for t, (doc, sizes) in per.items()
                for f in rs.check_trans_layout([r for r in rr if r["phenotype_type"] == t], doc, sizes)]
    rr, per = lay({"ge": rs.trans_frame_order(rows, [0, 1]), "leafcutter": order})
    assert check(rr, per) == []
    # a gene page reads one range per object: A's two introns, 10 + 3 and 10 + 6 bytes
    assert rs.gene_trans_ranges(rr, "A") == {"ge": (64, 10), "leafcutter": (64, 29)}
    rr, per = lay({"ge": [0, 1], "leafcutter": lc})
    fails = check(rr, per)
    assert len(fails) == 1 and "2 genes whose frames are not one contiguous range" in fails[0], fails
    rr[3]["trans_off"] += 1                                  # a gap between frames
    assert any("previous frame ends" in f for f in check(rr, per))


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
    """The v1 GWAS object: lossless at the source's four decimals, and refusing more."""
    from . import gwas
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_with_trans(Path(d), version=1)
        doc = one(docs)
        g = one(docs, "gwas")
        assert g["id"] == "test_gwas" and g["n_rows"] == 4 and set(g["files"]) == {"chr1", "chr2"}
        assert g["n_values"] == [900, 1000] and g["orientation"]["swapped"] == 1
        h = qs.parse_file_header((st.immutable / g["files"]["chr1"]).read_bytes())
        assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (gwas.KIND_GWAS, "chr1", 3, gwas.BLOCK_ROWS)
        w = gwas.read_window(st, g, "chr1", 1, 2)
        assert w["pos"].tolist() == [1, 2] and w["ref"] == ["A", "C"] and w["alt"] == ["C", "CA"]
        assert w["beta"].tolist() == [-0.5, 0.0123] and w["af"].tolist() == [0.25, 0.3333]
        assert w["p"].tolist() == [1.234e-9, 2.5e-3] and w["rs_number"].tolist() == [0, 2] and w["n"].tolist() == [1000, 1000]
        assert gwas.read_window(st, g, "chr2", 10, 20)["pos"].size == 0
        assert g["rows_sharing_a_site"] == 0
        t = Path(d) / "tables"
        gw = pq.read_table(t / "gwas.parquet").to_pylist()
        # a site reported twice (the source's two allele orders of one indel): both rows kept, counted
        pq.write_table(pa.Table.from_pylist(gw + [{**gw[1], "beta": 0.3, "n": 900, "pvalue": 0.04}]), t / "gwas.parquet")
        g2 = gwas.build(st, t, st.load(qs.CATALOGS, "cat1"), ["chr1", "chr2"], codec=pf)
        assert g2["rows_sharing_a_site"] == 2 and g2["n_rows"] == 5 and g2["precision"] is None
        # lossless means lossless: more than 4 decimals is refused, not rounded
        gw[0]["beta"] = 0.12345
        pq.write_table(pa.Table.from_pylist(gw), t / "gwas.parquet")
        try:
            gwas.build(st, t, st.load(qs.CATALOGS, "cat1"), ["chr1", "chr2"], codec=pf)
            raise AssertionError("a 5-decimal beta was accepted")
        except ValueError as e:
            assert "more than 4 decimals" in str(e), e


def test_gwas_v2_quantizes_and_declares_its_version():
    """The v2 GWAS path, which no store held until the codec was threaded into `gwas.build`.

    Three things have to hold at once. The header says 2, which is the only discriminator a `.qbg`
    has. The rows are smaller -- 18 bytes against v1's 21 -- because 9 bytes of lossless p/beta/se
    became 6 codes. And the 5-decimal beta v1 *refuses* is now accepted and comes back within the
    block's limit, which is the whole trade in plan D2: quantize and record the error, rather than
    refuse the source.
    """
    from . import gwas
    from . import packfmt_v2 as p2
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        t = write_tables(d)
        add_trans_and_gwas(t)
        gw = pq.read_table(t / "gwas.parquet").to_pylist()
        five_dp = next(r for r in gw if r["chr"] == "chr1" and r["pos"] == 1)
        five_dp["beta"] = 0.12345                   # v1 refuses this; v2 quantizes it
        pq.write_table(pa.Table.from_pylist(gw), t / "gwas.parquet")
        raised = ""
        try:
            gwas.build(st, t, cdoc, ["chr1", "chr2"], codec=pf)
        except ValueError as e:
            raised = str(e)
        assert "more than 4 decimals" in raised, f"v1 accepted a 5-decimal beta: {raised!r}"

        g = gwas.build(st, t, cdoc, ["chr1", "chr2"])          # the default codec is v2
        # v2 measures what quantizing cost, in the same block shape as a QTL scan's
        pr = g["precision"]
        assert pr["measured_worst"]["rows_compared"] == g["n_rows"] > 0
        for k in ("neglog10p", "beta_over_se", "se_rel"):
            assert pr["measured_worst"][k] <= pr["theoretical_limit"][k] * (1 + 1e-9), (k, pr)
        h = qs.parse_file_header((st.immutable / g["files"]["chr1"]).read_bytes())
        assert h["version"] == 2 and h["kind"] == gwas.KIND_GWAS
        assert qs.parse_file_header((st.immutable / g["index"]).read_bytes())["version"] == 2
        assert p2.GWAS_ROW_BYTES == 18 and pf.GWAS_ROW_BYTES == 21

        # reading it with the v1 codec would decode 18-byte rows at a 21-byte stride, so a caller
        # asserting v1 is refused; with no codec named, the index header's version picks it
        doc2 = g      # `read_window` takes the GWAS experiment itself now
        try:
            gwas.read_window(st, doc2, "chr1", 1, 2, codec=pf)
            raise AssertionError("a v2 GWAS file was read with the v1 codec")
        except ValueError as e:
            assert "format version 2" in str(e), e
        assert gwas.read_window(st, doc2, "chr1", 1, 2)["pos"].tolist() == [1, 2]

        w = gwas.read_window(st, doc2, "chr1", 1, 2, codec=p2)
        assert w["pos"].tolist() == [1, 2] and w["ref"] == ["A", "C"] and w["alt"] == ["C", "CA"]
        assert w["n"].tolist() == [1000, 1000] and w["rs_number"].tolist() == [0, 2]
        src = {r["pos"]: r for r in gw if r["chr"] == "chr1"}
        lim = p2.theoretical_limit(*_gwas_scales(st, g, "chr1"))
        for i, pos in enumerate(w["pos"].tolist()):
            assert abs(w["beta"][i] - src[pos]["beta"]) <= lim["beta_over_se"] * w["se"][i] + 1e-12
            assert abs(w["se"][i] / src[pos]["se"] - 1) <= lim["se_rel"] + 1e-12
        assert abs(w["beta"][0] - 0.12345) < 1e-4, w["beta"][0]


def _gwas_scales(st, g, chrom):
    """The four f64 scales out of a v2 GWAS block's own header, which is what a reader has."""
    from . import packfmt_v2 as p2
    buf = (st.immutable / g["files"][chrom]).read_bytes()
    payload = p2.zstd_unframe(buf[qs.HEADER_LEN:], None, "gwas block")
    _, _, _, nlp_max, beta_max, lse_min, lse_max = p2._GWAS_HEADER.unpack_from(payload, 0)
    return nlp_max, beta_max, lse_min, lse_max


def test_v2_trans_is_one_frame_per_phenotype_with_se_stored():
    """A v2 trans object has v1's framing -- one zstd frame per phenotype, located by the search
    index, a gene's frames contiguous -- and differs only inside the frame: the `QTT3` magic, the
    two SE scales, and a stored `se` column that no longer depends on a dof.

    The framing was merged with the GWAS block's once (`plans/2026-10-07-experiment-model.md` D20)
    and reversed: cutting a phenotype's rows at every chromosome turned 100,000 frames into 1.7
    million blocks, because a trans phenotype's hits sit on a median of 16 chromosomes. So a trans
    document carries no index, and `read_trans` picks the codec from the object's own header."""
    from . import packfmt_v2 as p2
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d, version=2)
        build_annotation(st, version=2)
        t = write_tables(d)
        add_trans_and_gwas(t)
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], codec=p2)
        st.write_store("t", [], 2)
        ge = one(docs, "qtl_trans", "ge")
        assert "index" not in ge and "n_blocks" not in ge
        body = (st.immutable / ge["file"]).read_bytes()
        assert qs.parse_file_header(body)["version"] == 2
        rows = {x["phenotype_id"]: x for x in rs.load_index(st, ge).to_pylist()}
        # one frame per phenotype: its span is exactly one v2 frame, which opens with the v2 magic
        r = rows["G1"]
        span = body[r["trans_off"]:r["trans_off"] + r["trans_len"]]
        assert p2.zstd_unframe(span, None, "frame")[:4] == p2.MAGIC_TRANS
        assert ge["n_phenotypes"] == qs.parse_file_header(body)["count"]
        # G1's rows span two chromosomes inside the one frame, each row naming its own
        g1 = rs.read_trans(st, ge, r)
        assert g1["chr"] == ["chr1", "chr2"] and g1["pos"].tolist() == [3, 2]
        assert g1["ref"] == ["G", "T"] and g1["alt"] == ["A", "G"]
        np.testing.assert_allclose(g1["beta"], [-0.3, 0.8], atol=0.8 / 32766 + 1e-7)
        # se is *stored*, not rebuilt from p and dof, so it survives a null dof
        np.testing.assert_allclose(g1["se"], [0.05, 0.1], rtol=2e-4)
        ge["dof"] = None
        np.testing.assert_allclose(rs.read_trans(st, ge, r)["se"], [0.05, 0.1], rtol=2e-4)
        # the cis reader picks the codec from the object's header the same way: a v2 block decodes
        # to a stored beta, with no slope rebuilt and no dof consulted
        cis = one(docs, "qtl_cis", "ge")
        g1 = next(x for x in rs.load_index(st, cis).to_pylist() if x["phenotype_id"] == "G1")
        blk = rs.read_block(st, dict(cis, dof=None), g1)
        assert "slope" not in blk and np.isfinite(blk["beta"]).any() and blk["details"]["phenotype_id"] == "G1"
        assert st.validate(refget=rg) == [], st.validate(refget=rg)


def test_collection_metadata_comes_from_config_when_given():
    """`collections` (config.yaml's map) names the publication; it wins over what the adapter
    recorded and a missing entry falls back to it. Keyed by the cohort's stem and by a GWAS's id."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        build_annotation(st)
        t = write_tables(d)
        add_trans_and_gwas(t)
        meta = {"exp": {"title": "A cohort paper", "citation": {"doi": "10.1/abc", "year": 2026},
                        "source": {"record": "zenodo_1"}},
                "test_gwas": {"title": "A GWAS paper", "citation": {"doi": "10.1/xyz", "year": 2024}}}
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], collections=meta)
        qc = json.loads((st.immutable / one(docs, "qtl_cis")["collection"]).read_text())
        gc_ = json.loads((st.immutable / one(docs, "gwas")["collection"]).read_text())
        assert qc["title"] == "A cohort paper" and qc["citation"]["doi"] == "10.1/abc" and qc["source"] == {"record": "zenodo_1"}
        assert gc_["title"] == "A GWAS paper" and gc_["citation"]["year"] == 2024
        assert gc_["source"] == {}               # no source given here or by the fixture's adapter
        assert len({d["collection"] for d in docs.values() if d["kind"] != "gwas"}) == 1
        # no entry at all: the adapter's title and source, as before
        docs2 = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"])
        assert json.loads((st.immutable / one(docs2, "gwas")["collection"]).read_text())["title"] == "a test GWAS"


def test_gwas_is_its_own_collection():
    """The defect this restructure exists to remove, in the form it can come back in.

    The DCM GWAS is Jurgens et al. 2024; TOPCHeF colocalized against it and did not produce it. It
    used to be a *field* on the TOPCHeF experiment, so `remove_experiment("topchef")` plus `gc`
    would have deleted another publication's data. Making it its own experiment fixes that only if
    it also gets its own collection: handing it the cohort's digest would move the misattribution
    from a field into a digest, and make one publication's citation describe another's."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_with_trans(Path(d))
        gw, qtl = one(docs, "gwas"), one(docs, "qtl_cis")
        assert gw["collection"] != qtl["collection"]
        gc_ = json.loads((st.immutable / gw["collection"]).read_text())
        qc = json.loads((st.immutable / qtl["collection"]).read_text())
        assert gc_["title"] == "a test GWAS" and qc["title"] != gc_["title"]
        # the GWAS collection records the GWAS source, not the cohort's archives
        assert "zenodo" not in json.dumps(gc_)
        # and every QTL scan still shares one collection with the others
        assert len({d["collection"] for d in docs.values() if d["kind"] != "gwas"}) == 1


def test_gwas_bins():
    """The bin summary: v0's values and rounding, catalog table order, named from `gwas.bins`, checked by
    validate; bad bins refused at build and caught by validate."""
    from . import catalog, gwas
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_with_trans(Path(d))
        doc = one(docs)
        b = one(docs, "gwas")["bins"]
        assert (b["bin_bp"], b["n_bins"]) == (5, 2) and b["file"].endswith(".arrow.zst")
        assert b["file"] in qs.object_names(one(docs, "gwas"))
        gw = one(docs, "gwas")
        t = gwas.read_bins(st, gw)
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
        # validate: the GWAS experiment naming a bin summary whose count, bins or chromosomes are wrong
        good = gw["bins"]
        for change, rule in (({"n_bins": 3}, "the pointer says 3"), ({"bin_bp": 10}, "not [k * 10")):
            st.write_pointer("experiments", gw["id"], {**gw, "bins": {**good, **change}})
            fails = st.validate(refget=rg, sites=catalog.read_sites)
            assert any("gwas bins" in f and rule in f for f in fails), fails
        write_bins(t_dir, GWAS_BINS)
        only1 = gwas.build(st, t_dir, cdoc, ["chr1", "chr2"])
        st.write_pointer("experiments", gw["id"],
                         {**gw, **only1, "files": {"chr1": only1["files"]["chr1"]}, "id": gw["id"]})
        assert any("have bins and no GWAS file" in f for f in st.validate(refget=rg, sites=catalog.read_sites))
        st.write_pointer("experiments", gw["id"], gw)
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
        st, rg, docs = build_with_trans(d)
        kept = {n for dd in docs.values() for n in qs.object_names(dd)}
        before = {n: (st.immutable / n).read_bytes() for n in kept}
        other = d / "other"
        shutil.copytree(write_tables(d / "o", with_new_columns=False), other)
        docs2 = rs.build(st, "exp2", other, "cat1", "ann", ["chr1", "chr2"])
        st.write_store("t", [])
        assert st.validate(refget=rg, sites=catalog.read_sites) == []
        own2 = {n for dd in docs2.values() for n in qs.object_names(dd)}
        only2 = own2 - kept - set(qs.object_names(st.load(qs.CATALOGS, "cat1")))
        assert only2, "exp2 should own some objects"
        # its collection object among them: a collection is a field, so the only thing keeping that
        # document alive is the experiments naming its digest
        coll2 = one(docs2)["collection"]
        assert coll2 in only2 and coll2 != one(docs)["collection"]
        assert qs.gc(st, dry_run=True)["would_delete"] == 0            # everything is referenced
        for eid in sorted(docs2):                  # one scan is one experiment, so removing a
            qs.remove_experiment(st, eid)          # cohort is removing each of its scans
        assert sorted(json.loads((st.root / "store.json").read_text())["experiments"]) == sorted(docs)
        assert st.validate(refget=rg, sites=catalog.read_sites) == []
        r = qs.gc(st)
        assert r["deleted"] == len(only2) and not any((st.immutable / n).exists() for n in only2)
        assert {n: (st.immutable / n).read_bytes() for n in kept} == before
        assert st.validate(refget=rg, sites=catalog.read_sites) == [] and qs.gc(st)["deleted"] == 0
        try:
            qs.remove_experiment(st, "exp2")
            raise AssertionError("removing a missing experiment succeeded")
        except ValueError as e:
            assert "does not exist" in str(e)


def test_blocks_codes_and_gaps():
    """The v1 block's codes: the interleaved (nlp, se) pair with the slope's sign in the SE field."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d), version=1)
        doc = one(docs)
        r = rows_by_id(st, docs)
        # each scan numbers its own rows from 0, so the merged view repeats ordinals: ge has G1, G2
        # and leafcutter has I1, I2, I3, where one shared index numbered all five 0..4
        assert [r[k]["ord"] for k in ("G1", "G2", "I1", "I2", "I3")] == [0, 1, 0, 1, 2]
        g1 = read_row(st, docs, r["G1"])
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


def test_builds_with_no_credible_sets():
    """An experiment with no fine-mapping writes `credible_sets` with zero rows (CONTRACT.md), and a
    chromosome of a normal experiment can have none either. Both used to crash: an empty pandas
    slice loses its string dtypes, so `_vidx`'s merge raised "str and float64" on `ref` rather than
    returning nothing. ARIC's plasma pQTL has no fine-mapping at all, so every chromosome hit it."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        build_annotation(st)
        t = write_tables(d)
        cs = pq.read_table(t / "credible_sets.parquet")
        pq.write_table(cs.slice(0, 0), t / "credible_sets.parquet")      # keep the table, drop the rows
        ing = json.loads((t / "ingestion.json").read_text())
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], ingestion=ing)
        doc = one(docs)
        r = rows_by_id(st, docs)
        assert doc["counts"], doc
        assert read_row(st, docs, r["G1"])["details"]["n_credible_sets"] == 0
        # the hits files still exist and still carry the lead records, just no credible-set kind
        h = rs.decode_hits((st.immutable / doc["hits"]["chr1"]).read_bytes())
        assert len(h) and not any(int(x["kind"]) == rs.HIT_CS for x in h)


def test_no_nominal_and_groups():
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d))
        doc = one(docs)
        r = rows_by_id(st, docs)
        g2 = read_row(st, docs, r["G2"])
        assert g2["n_rows"] == 0 and g2["details"]["has_nominal"] is False
        assert not r["G2"]["has_nominal"] and r["G2"]["var_start"] is None and not r["G2"]["significant"]
        for k in ("I2", "I3"):                            # non-lead introns: located by the cluster's lead
            assert r[k]["chr"] == "chr2" and not r[k]["has_nominal"] and not r[k]["is_group_lead"]
            b = read_row(st, docs, r[k])
            assert b["n_rows"] == 0 and b["details"]["group"]["lead_phenotype_id"] == "I1"
            assert b["details"]["phenotype_object_id"] == "clu_1" and b["details"]["extra"] == {"cluster_id": "clu_1"}
        i1 = read_row(st, docs, r["I1"])
        assert r["I1"]["is_group_lead"] and r["I1"]["significant"] and i1["n_rows"] == 3


def test_hits():
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d))
        doc = one(docs)
        lc = one(docs, "qtl_cis", "leafcutter")
        hits = lambda d, c: rs.decode_hits((st.immutable / d["hits"][c]).read_bytes())    # noqa: E731
        h1 = hits(doc, "chr1")
        got = [(int(x["vidx"]), int(x["kind"]), int(x["ord"]), int(x["cs_id"]), int(x["flags"])) for x in h1]
        assert got == [(0, 0, 0, 0, 1), (0, 1, 0, 1, 0), (1, 0, 1, 0, 0), (1, 1, 0, 1, 0)]
        assert np.isclose(h1["value"][1], 0.9)
        # Hits are per scan: the eQTL leads and credible sets are all on chr1, the sQTL ones all on
        # chr2, so each experiment's other chromosome is empty where one object held both.
        assert len(hits(doc, "chr2")) == 0 and len(hits(lc, "chr1")) == 0
        assert [(int(x["vidx"]), int(x["kind"]), int(x["ord"])) for x in hits(lc, "chr2")] \
            == [(0, 0, 0), (0, 1, 0)]


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
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], ingestion={**ing, "significance": rule})
        doc = one(docs)
        assert doc["significance"] == rule
        r = rows_by_id(st, docs)
        assert r["G1"]["significant"] and not r["I1"]["significant"] and not r["I2"]["significant"]
        # the index reports the column the rule tested, not p_perm (I1's p_perm is 0.01)
        assert r["I1"]["sig_value"] == 0.02 and r["G1"]["sig_value"] == 0.002
        # the block's details keep the source's own numbers under their own names
        gi = read_row(st, docs, r["I1"])["details"]["group"]
        assert gi["p_perm"] == 0.01 and gi["p_beta"] == 0.02 and not gi["significant"]
        assert read_row(st, docs, r["G1"])["details"]["group"]["significant"]
        # I1 is a leafcutter phenotype, so its lead hit is in that scan's object, not the eQTL one
        lc = one(docs, "qtl_cis", "leafcutter")
        h2 = rs.decode_hits((st.immutable / lc["hits"]["chr2"]).read_bytes())
        assert [int(x["flags"]) for x in h2 if x["kind"] == 0] == [0]
        # the hit's value is the tested quantity too
        assert [round(float(x["value"]), 3) for x in h2 if x["kind"] == 0] == [0.02]
        for bad in ({"column": "pvalue", "op": "<", "threshold": 0.05}, {"column": "p_beta", "op": "~", "threshold": 1}):
            try:
                rs.build(st, "exp2", t, "cat1", "ann", ["chr1", "chr2"], ingestion={**ing, "significance": bad})
            except ValueError as e:
                assert "significance rule" in str(e), e
            else:
                raise AssertionError(f"no error for {bad}")


def test_significance_not_assessed():
    """`significance: null` -- a source that ran no permutation pass. Every phenotype's `significant`
    is null, not false: "never assessed" and "tested and not significant" are different claims. The
    type counts say null rather than 0, the lead hits still exist with the flag clear, and a rule on
    a column `permuted` does not have is an error rather than a silently false one."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st, rg, cdoc = build_catalog(d)
        build_annotation(st)
        t = write_tables(d)
        ing = json.loads((t / "ingestion.json").read_text())
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], ingestion={**ing, "significance": None})
        doc = one(docs)
        assert doc["significance"] is None
        r = rows_by_id(st, docs)
        assert all(v["significant"] is None for v in r.values()), {k: v["significant"] for k, v in r.items()}
        assert all(v["sig_value"] is None for v in r.values())
        for dd in docs.values():
            c = dd["counts"]
            assert c["significant"] is None and c["significant_genes"] is None and c["phenotypes"] > 0
        assert read_row(st, docs, r["G1"])["details"]["group"]["significant"] is None
        h2 = rs.decode_hits((st.immutable / one(docs, "qtl_cis", "leafcutter")["hits"]["chr2"]).read_bytes())
        assert [int(x["flags"]) for x in h2 if x["kind"] == 0] == [0]
        try:
            rs.build(st, "exp2", t, "cat1", "ann", ["chr1", "chr2"],
                     ingestion={**ing, "significance": {"column": "pval_perm", "op": "<", "threshold": 0.05}})
        except ValueError as e:
            assert "has no column 'pval_perm'" in str(e), e
        else:
            raise AssertionError("a rule naming a column `permuted` lacks must fail")


def test_migrate_index_renames_p_perm():
    """`migrate-index` renames a stored index's `p_perm` to `sig_value` and redoes the parts and
    counts, leaving every other object alone. Idempotent, and it refuses an experiment whose rule
    tested something other than `p_perm`, where the stored column is not the tested value."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d))
        doc = one(docs)
        before = {k: v for k, v in doc.items() if k in ("results", "hits", "catalog", "annotation")}
        # an index as a store built before the rename holds it
        idx = rs.load_index(st, doc)
        old = idx.rename_columns(["p_perm" if n == "sig_value" else n for n in idx.schema.names])
        stale = {**doc, "search_index": st.put(rs.encode_arrow(old), rs.EXT_INDEX)}
        st.write_pointer("experiments", doc["id"], stale)
        assert "p_perm" in rs.load_index(st, st.load("experiments", doc["id"])).schema.names

        out = rs.migrate_index(st, doc["id"])
        names = rs.load_index(st, out).schema.names
        assert "sig_value" in names and "p_perm" not in names
        assert {k: v for k, v in out.items() if k in before} == before   # nothing else was touched
        assert out["counts"] == doc["counts"] and len(out["search_index_parts"]) == len(doc["search_index_parts"])
        for c, p in out["search_index_parts"].items():
            assert "sig_value" in rs.decode_arrow((st.immutable / p["file"]).read_bytes()).schema.names
        assert rs.migrate_index(st, doc["id"])["search_index"] == out["search_index"]   # idempotent

        st.write_pointer("experiments", "exp_beta", {**stale, "id": "exp_beta",
                                                     "significance": {"column": "p_beta", "op": "<", "threshold": 0.05}})
        try:
            rs.migrate_index(st, "exp_beta")
        except ValueError as e:
            assert "cannot be renamed" in str(e), e
        else:
            raise AssertionError("a rule on another column must not be silently renamed")


def test_defaults_for_old_tables():
    """Tables written before the contract added phenotype_object_id and has_nominal still build,
    with phenotype_object_id = phenotype_id and has_nominal from the nominal rows."""
    with tempfile.TemporaryDirectory() as d:
        st, rg, docs = build_all(Path(d), with_new_columns=False)
        doc = one(docs)
        r = rows_by_id(st, docs)
        assert r["I1"]["phenotype_object_id"] == "I1" and r["I1"]["has_nominal"] and r["G1"]["has_nominal"]
        # without the group column I2 and I3 are their own groups, with no permuted row: unplaceable,
        # and counted rather than dropped -- on the leafcutter scan they belong to, not the eQTL one
        assert "I2" not in r and doc["unplaced"] == {"count": 0, "examples": []}
        assert one(docs, "qtl_cis", "leafcutter")["unplaced"] == {"count": 2, "examples": ["I2", "I3"]}


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
    # an unusable dof fit stores `dof: null`; under v1 read_block must not invent slopes from a
    # stand-in dof (v2 stores beta and never rebuilds a slope, so this is a v1 case)
    with tempfile.TemporaryDirectory() as d:
        st, rg, cdoc = build_catalog(Path(d), version=1)
        build_annotation(st, version=1)
        t = write_tables(Path(d))
        ing = json.loads((t / "ingestion.json").read_text())
        ing["dof"]["ge"] = {"dof": 12, "usable": False, "reason": "no fit"}
        (t / "ingestion.json").write_text(json.dumps(ing))
        docs = rs.build(st, "exp", t, "cat1", "ann", ["chr1", "chr2"], codec=pf)
        doc = one(docs)
        assert [one(docs, "qtl_cis", p)["dof"] for p in ("ge", "leafcutter")] == [None, 480]
        r = rows_by_id(st, docs)
        g1 = read_row(st, docs, r["G1"])
        assert g1["n_rows"] == 3 and np.isnan(g1["slope"]).all()
        assert np.isfinite(g1["nlp"][[0, 2]]).all() and g1["se_code"][2] & pf.SE_SIGN   # p, SE, sign survive
        i1 = read_row(st, docs, r["I1"])                                             # dof 480: slopes rebuilt
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
