"""Tests for the catalog overlap index (pipeline/overlap.py, SPEC.md section 19).

    uv run python -m pipeline.test_overlap      (or pytest)

Two catalogs over the tiny real refgetstore from test_qtlstore, built so that one of them has a
trans-only site on a chromosome where the other does not. That is the case `vidx_in` exists for:
a rank over the presence mask is *not* a vidx once a chromosome carries trans-only sites, and an
implementation that returns the rank passes every test built without one.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pyarrow as pa

from . import catalog as cat
from . import overlap as ov
from . import qtlstore as qs
from .test_qtlstore import SEQ, LEN, refgetstore     # noqa: F401 -- LEN is used via the fixture


def table(rows) -> pa.Table:
    """(chr, pos, ref, alt, in_cis) -> a contract `sites` table with the required columns."""
    c = list(zip(*rows))
    return pa.table({"chr": pa.array(c[0]), "pos": pa.array(c[1], pa.int32()), "ref": pa.array(c[2]),
                     "alt": pa.array(c[3]), "af": pa.array([0.25] * len(rows), pa.float32()),
                     "ma_samples": pa.array([2] * len(rows), pa.int32()),
                     "in_cis": pa.array(c[4])})


# A: chr1 {1 A>C, 3 G>A, 5 T>G}, chr2 {2 C>T}
# B: chr1 {1 A>C, 2 C>CA, 5 T>G}, chr2 {2 C>T, 4 G>A}
# A's chr1 site at 3 is trans-only, so A's chr1 vidx order is (1, 5) then (3): the exact case where
# union order and vidx order disagree.
A_ROWS = [("chr1", 1, "A", "C", True), ("chr1", 5, "T", "G", True), ("chr1", 3, "G", "A", False),
          ("chr2", 2, "C", "T", True)]
B_ROWS = [("chr1", 1, "A", "C", True), ("chr1", 2, "C", "CA", True), ("chr1", 5, "T", "G", True),
          ("chr2", 2, "C", "T", True), ("chr2", 4, "G", "A", True)]


def build(d: Path, a_rows=None, b_rows=None, page_size: int = 2, chunk: int = 4):
    rg = refgetstore(d / "refget")
    coll = rg.list_collections(page_size=10)["results"][0].digest
    st = qs.Store(d / "store")
    cat.build(st, "cat_a", table(a_rows or A_ROWS), rg, coll, ["chr1", "chr2"], page_size=page_size)
    cat.build(st, "cat_b", table(b_rows or B_ROWS), rg, coll, ["chr1", "chr2"], page_size=page_size)
    doc = ov.build(st, "ovl", ["cat_a", "cat_b"], page_size=page_size, chunk=chunk)
    return st, doc, ov.decode((st.immutable / doc["object"]).read_bytes())


def test_union_and_counts():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        assert doc["union"] == 6, doc["union"]           # 1,2,3,5 on chr1; 2,4 on chr2
        assert doc["mask_width"] == 1 and idx["k"] == 2
        assert [c["id"] for c in doc["catalogs"]] == ["cat_a", "cat_b"]
        p = doc["pairs"][0]
        assert (p["shared"], p["a_only"], p["b_only"]) == (3, 1, 2), p
        assert abs(p["jaccard"] - 3 / 6) < 1e-9
        assert doc["by_chrom"]["chr1"]["union"] == 4 and doc["by_chrom"]["chr2"]["union"] == 2


def test_union_sites_are_in_canonical_order():
    """Chromosomes by seq_digest in ASCII order, then pos, ref, alt (SPEC section 5) -- *not* by
    chromosome name and not in the catalog's own table order. In this fixture chr2's digest sorts
    before chr1's, so chr2's sites come first; a reader that assumes name order reads the wrong
    sites."""
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        got = [(p, r, a) for _, p, r, a in ov.union_sites(idx)]
        assert got == [(2, "C", "T"), (4, "G", "A"),                              # chr2
                       (1, "A", "C"), (2, "C", "CA"), (3, "G", "A"), (5, "T", "G")], got
        assert [u for u, *_ in ov.union_sites(idx)] == list(range(6))
        digests = [c["seq_digest"] for c in idx["chroms"]]
        assert digests == sorted(digests)
        names = {c["seq_digest"]: c["name"] for c in st.load(qs.CATALOGS, "cat_a")["chromosomes"]}
        assert [names[x] for x in digests] == ["chr2", "chr1"], "the fixture no longer exercises digest order"


def test_membership_bits_follow_the_pointer_order():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        by_key = {(p, r, a): u for u, p, r, a in ov.union_sites(idx)}
        both = ov.membership(idx, by_key[(1, "A", "C")])
        only_a = ov.membership(idx, by_key[(3, "G", "A")])
        only_b = ov.membership(idx, by_key[(2, "C", "CA")])
        assert both == 0b11 and only_a == 0b01 and only_b == 0b10


def test_vidx_is_not_a_rank():
    """The trap. cat_a's chr1 vidx order is 1 A>C (0), 5 T>G (1), then the trans-only 3 G>A (2),
    because a catalog numbers a chromosome's cis sites before its trans-only ones. Union order
    interleaves them, so at 5 T>G cat_a's rank is 3 and its vidx is 1: an implementation that returns
    the rank is wrong here and right on every catalog with no trans-only section."""
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        by_key = {(p, r, a): u for u, p, r, a in ov.union_sites(idx)}
        u_three, u_five = by_key[(3, "G", "A")], by_key[(5, "T", "G")]
        assert ov.rank(idx, 0, u_five) == 3 and ov.vidx_in(idx, 0, u_five) == 1, "rank was returned as a vidx"
        assert ov.rank(idx, 0, u_three) == 2 and ov.vidx_in(idx, 0, u_three) == 2
        # cat_b has no trans-only site, so for it rank and vidx do agree within a chromosome: the
        # reason this bug hides
        assert ov.rank(idx, 1, u_five) == 4 and ov.vidx_in(idx, 1, u_five) == 2
        # and every site's vidx round-trips against the catalog's own pages
        for i, cid in enumerate(["cat_a", "cat_b"]):
            cdoc = st.load(qs.CATALOGS, cid)
            for c in cdoc["chromosomes"]:
                d0 = cat.decode_file((st.immutable / c["file"]).read_bytes())
                for v, (p, r, a) in enumerate(zip(d0["pos"].tolist(), d0["ref"], d0["alt"])):
                    assert ov.vidx_in(idx, i, by_key[(p, r, a)]) == v, (cid, c["name"], p, r, a, v)


def test_vidx_none_when_absent():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        by_key = {(p, r, a): u for u, p, r, a in ov.union_sites(idx)}
        assert ov.vidx_in(idx, 1, by_key[(3, "G", "A")]) is None      # cat_b has no site there
        assert ov.vidx_in(idx, 0, by_key[(2, "C", "CA")]) is None


def test_rank_across_a_chunk_boundary():
    """chunk=4 puts the union's last two sites in chunk 1, so rank has to add a prefix to a scan."""
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d), chunk=4)
        assert idx["n_chunks"] == 2, idx["n_chunks"]
        assert ov.rank(idx, 0, idx["n_union"]) == 4     # cat_a has 4 sites
        assert ov.rank(idx, 1, idx["n_union"]) == 5     # cat_b has 5
        # chunk 0 is union 0..3: cat_a holds 2 of them, cat_b all 4
        assert int(idx["prefix"][0, 1]) == 2 and int(idx["prefix"][1, 1]) == 4


def test_find_locates_a_site_by_key():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        sd1 = st.load(qs.CATALOGS, "cat_a")["chromosomes"][0]["seq_digest"]
        u = ov.find(idx, sd1, 5, "T", "G")
        assert u is not None and ov.membership(idx, u) == 0b11
        assert ov.find(idx, sd1, 4, "A", "T") is None
        assert ov.find(idx, "z" * 32, 1, "A", "C") is None


def test_mask_width_bounds():
    assert (ov.mask_width(1), ov.mask_width(8), ov.mask_width(9), ov.mask_width(64)) == (1, 1, 2, 8)
    for k in (0, 65):
        try:
            ov.mask_width(k)
        except ValueError:
            continue
        raise AssertionError(f"mask_width({k}) should refuse")


def test_validate_passes_and_catches_a_stale_catalog():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        st.write_store("t", [])
        assert st.validate() == [], st.validate()
        # rebuilding a catalog with different sites changes its identity; the index is now stale
        rg = refgetstore(Path(d) / "refget")
        coll = rg.list_collections(page_size=10)["results"][0].digest
        cat.build(st, "cat_a", table(A_ROWS + [("chr2", 8, "A", "T", True)]), rg, coll, ["chr1", "chr2"], page_size=2)
        fails = st.validate()
        assert any("rebuild the index" in f for f in fails), fails


def test_conflicts_and_suspects_are_reported():
    with tempfile.TemporaryDirectory() as d:
        # same position, disjoint allele pairs -> a conflict; two deletions of one base a few bp
        # apart, in different catalogs -> a normalisation suspect
        a = [("chr1", 1, "A", "C", True), ("chr1", 10, "GT", "G", True)]
        b = [("chr1", 1, "A", "G", True), ("chr1", 12, "CT", "C", True)]
        st, doc, idx = build(Path(d), a, b)
        assert doc["conflicts"]["shared_positions_no_shared_allele_pair"] == 1, doc["conflicts"]
        assert doc["conflicts"]["examples"][0]["pos"] == 1
        assert doc["normalisation_suspects"]["count"] == 1, doc["normalisation_suspects"]
        ex = doc["normalisation_suspects"]["examples"][0]
        assert {ex["a"]["pos"], ex["b"]["pos"]} == {10, 12}


def test_refuses_catalogs_on_different_collections():
    with tempfile.TemporaryDirectory() as d:
        st, doc, idx = build(Path(d))
        bad = st.load(qs.CATALOGS, "cat_b")
        bad["collection_digest"] = "Z" * 32
        (st.root / qs.CATALOGS / "cat_b.json").write_text(__import__("json").dumps(bad))
        try:
            ov.build(st, "ovl2", ["cat_a", "cat_b"])
        except ValueError as e:
            assert "different references" in str(e), e
            return
        raise AssertionError("an index across two sequence collections should be refused")


def main() -> int:
    bad = 0
    tests = [v for k, v in globals().items() if k.startswith("test_")]
    for fn in tests:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(tests) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
