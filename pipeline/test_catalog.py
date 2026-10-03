"""Tests for the variant catalog builder (pipeline/catalog.py).

    uv run python -m pipeline.test_catalog      (or pytest)

A tiny real refgetstore (from test_qtlstore) and a hand-written `sites` table.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa

from . import catalog as cat
from . import packfmt_v1 as pf
from . import qtlstore as qs
from .test_qtlstore import SEQ, LEN, refgetstore


def sites_table() -> pa.Table:
    # deliberately unsorted; chr2 has a multi-allelic position, an indel and a trans-only site
    rows = [
        ("chr2", 9, "T", "A", 0.7, 3, 4, True, "rs5", 5, "exact"),
        ("chr1", 3, "G", "A", 0.1, 10, 12, True, "rs3", 3, "exact"),
        ("chr1", 1, "A", "C", 0.4, 5, 5, True, None, -1, "none"),
        ("chr2", 9, "T", "C", float("nan"), -1, -1, True, "rs5", 5, "position"),
        ("chr2", 7, "TTG", "T", 0.25, 2, 2, True, "rs9", 9, "exact"),
        ("chr2", 2, "T", "G", float("nan"), -1, -1, False, "rs1", 1, "position"),
        ("chr1", 2, "C", "CA", 0.5, 7, 9, True, "rs2", 2, "exact"),
    ]
    cols = list(zip(*rows))
    return pa.table({"chr": pa.array(cols[0]), "pos": pa.array(cols[1], pa.int32()), "ref": pa.array(cols[2]),
                     "alt": pa.array(cols[3]), "af": pa.array(cols[4], pa.float32()),
                     "ma_samples": pa.array(cols[5], pa.int32()), "ma_count": pa.array(cols[6], pa.int32()),
                     "in_cis": pa.array(cols[7]), "rsid": pa.array(cols[8]), "rs_number": pa.array(cols[9], pa.int64()),
                     "match": pa.array(cols[10])})


def build(d: Path, page_size: int = 2, table: pa.Table | None = None, chroms=("chr1", "chr2")):
    rg = refgetstore(d / "refget")
    coll = rg.list_collections(page_size=10)["results"][0].digest
    st = qs.Store(d / "store")
    doc = cat.build(st, "cat1", table if table is not None else sites_table(), rg, coll, list(chroms),
                    page_size=page_size)
    return st, rg, doc


def raises(rule: str, fn, *args, **kwargs) -> None:
    try:
        fn(*args, **kwargs)
    except ValueError as e:
        assert rule in str(e), f"expected {rule!r} in error, got {e!r}"
        return
    raise AssertionError(f"no ValueError raised (expected {rule!r})")


def test_layout_and_order():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build(Path(d))
        assert [c["name"] for c in doc["chromosomes"]] == ["chr1", "chr2"]
        assert doc["attributes"] == ["af", "ma_samples", "ma_count", "rs_number", "match"]   # no rsid text
        c1, c2 = doc["chromosomes"]
        assert (c1["seq_digest"], c1["length"], c1["count"], c1["n_cis"]) == (SEQ["chr1"], LEN["chr1"], 3, 3)
        assert (c2["count"], c2["n_cis"]) == (4, 3)
        assert c2["file_digest"] == c2["file"].split(".")[0]
        d2 = cat.load_chrom(st, doc, "chr2")
        # cis by (pos, ref, alt), then the trans-only section restarting at position 2
        assert list(zip(d2["pos"].tolist(), d2["ref"], d2["alt"])) == [(7, "TTG", "T"), (9, "T", "A"), (9, "T", "C"),
                                                                        (2, "T", "G")]
        h = qs.parse_file_header((st.immutable / c2["file"]).read_bytes())
        assert (h["kind"], h["count"], h["n_cis"], h["page_size"], h["seq_digest"]) == (1, 4, 3, 2, SEQ["chr2"])
        # af codes are v0's, ALT-relative; flags bit 0 alt_is_minor, bits 1-2 the match code
        assert d2["af_code"].tolist() == [round(0.25 * pf.AF_MAXQ), round(np.float32(0.7) * pf.AF_MAXQ), pf.AF_NULL,
                                          pf.AF_NULL]
        assert (d2["flags"] & 1).tolist() == [1, 0, 0, 0]
        assert (d2["flags"] >> 1).tolist() == [1, 1, 2, 2]
        assert d2["ma_samples"].tolist() == [2, 3, pf.COUNT_NULL, pf.COUNT_NULL]
        d1 = cat.load_chrom(st, doc, "chr1")
        assert d1["rs_number"].tolist() == [0, 2, 3]            # -1 is stored as 0 = none
        assert (d1["flags"] & 1).tolist() == [1, 0, 1]          # af 0.5 is not minor


def test_identity_and_validate():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build(Path(d))
        d1, d2 = (cat.load_chrom(st, doc, c) for c in ("chr1", "chr2"))
        want = qs.catalog_identity([(SEQ["chr1"], d1["pos"], d1["ref"], d1["alt"]),
                                    (SEQ["chr2"], d2["pos"], d2["ref"], d2["alt"])])
        assert doc["identity_digest"] == want
        st.write_store("t", [])
        assert st.validate(refget=rg, sites=cat.read_sites) == [] and st.notes == []
        # input order and page size do not change the identity
        t = sites_table()
        _, _, doc2 = build(Path(d) / "b", page_size=512, table=t.take(list(reversed(range(t.num_rows)))))
        assert doc2["identity_digest"] == doc["identity_digest"]
        assert doc2["chromosomes"][0]["file"] != doc["chromosomes"][0]["file"]


def test_identity_is_the_set_of_sites():
    """The same sites give the same identity whatever the cis/trans-only split and the chromosome table
    order, although both change the files (SPEC section 5, canonical order)."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        _, _, doc = build(d / "a")
        t = sites_table()
        flipped = t.set_column(t.column_names.index("in_cis"), "in_cis",
                               pa.array([not x for x in t["in_cis"].to_pylist()]))
        _, _, doc2 = build(d / "b", table=flipped, chroms=("chr2", "chr1"))
        assert [c["name"] for c in doc2["chromosomes"]] == ["chr2", "chr1"]
        assert doc2["chromosomes"][0]["n_cis"] == 1 and doc2["chromosomes"][1]["n_cis"] == 0
        assert doc2["identity_digest"] == doc["identity_digest"]
        assert {c["file"] for c in doc2["chromosomes"]}.isdisjoint(c["file"] for c in doc["chromosomes"])
        # one site fewer is a different identity
        _, _, doc3 = build(d / "c", table=t.slice(1))
        assert doc3["identity_digest"] != doc["identity_digest"]


def test_match_code_3_is_reserved():
    page = lambda m: cat.encode_page(0, [5], [0], [0.1], [1], [1], ["A"], ["G"], [m], codec="raw")  # noqa: E731
    raises("match code", page, 3)
    raw = bytearray(page(2))
    flags_at = pf.PAGE_HEADER_LEN + 4 + 15                      # payload: heap_len, then 15n bytes, then flags
    assert raw[flags_at] == 2 << 1 | 1
    head = qs.file_header(1, "chr1", 1, 512, 1, SEQ["chr1"])
    assert cat.decode_file(head + bytes(raw))["flags"].tolist() == [5]
    raw[flags_at] = 3 << 1
    raises("match code 3", cat.decode_file, head + bytes(raw))


def test_rsid_must_be_rs_number():
    """The variant catalog stores `rs_number`, never the rsID text, so a sites table's `rsid` has to be exactly
    "rs<rs_number>" (null where there is none) or the text would be lost."""
    with tempfile.TemporaryDirectory() as d:
        t = sites_table()
        i = t.column_names.index("rsid")
        wrong = t.set_column(i, "rsid", pa.array(["rs5", "rs3", None, "rs5", "rs9", "rs1", "rs22"]))
        raises("rsid is not", build, Path(d) / "a", table=wrong)
        unnamed = t.set_column(i, "rsid", pa.array(["rs5", "rs3", "rs0", "rs5", "rs9", "rs1", "rs2"]))
        raises("rsid is not", build, Path(d) / "b", table=unnamed)
        raises("without `rs_number`", build, Path(d) / "c", table=t.drop_columns(["rs_number"]))
        _, _, doc = build(Path(d) / "d", table=t.drop_columns(["rsid"]))
        assert "rs_number" in doc["attributes"]


def test_indexes_use_catalog_ordinals():
    with tempfile.TemporaryDirectory() as d:
        st, rg, doc = build(Path(d))
        names = [c["name"] for c in doc["chromosomes"]]
        vx = cat.decode_vidx((st.immutable / doc["vidx"]).read_bytes(), names)
        c2 = vx["chroms"]["chr2"]
        assert (c2["n_cis"], c2["n_trans"]) == (3, 1)
        assert c2["page_first_position"].tolist() == [7, 9, 2]      # pages of 2, trans-only starts a page
        size = (st.immutable / doc["chromosomes"][1]["file"]).stat().st_size
        assert c2["page_off"][0] == qs.HEADER_LEN and c2["page_off"][-1] == size
        r = cat.decode_rsid((st.immutable / doc["rsid"]).read_bytes())
        got = list(zip(r["rs_number"].tolist(), r["ordinal"].tolist(), r["vidx"].tolist()))
        # rs5 names two alleles at one position: both records kept
        assert got == [(1, 2, 3), (2, 1, 1), (3, 1, 2), (5, 2, 1), (5, 2, 2), (9, 2, 0)]
        assert vx["rsid_n"] == 6 and vx["rsid_first"].tolist() == [1]


def test_rejects_bad_sites():
    with tempfile.TemporaryDirectory() as d:
        t = sites_table()
        dup = pa.concat_tables([t, t.slice(0, 1)])
        for bad, rule in ((dup, "not unique"), (t.drop_columns(["in_cis"]), "required columns"),
                          (t.set_column(1, "pos", pa.array([10**6] * t.num_rows, pa.int32())), "outside")):
            try:
                build(Path(d) / rule.replace(" ", "_"), table=bad)
            except ValueError as e:
                assert rule in str(e), e
            else:
                raise AssertionError(f"no error for {rule}")


def test_rsid_run_across_a_block_boundary():
    """A repeated rs_number whose run crosses a block boundary: starting at the last block whose first
    number is <= rs would miss the run's first records; the lookup starts at the last block whose first
    number is below it."""
    B = cat.RSID_BLOCK_RECORDS
    rs = np.arange(1, B + 11, dtype=np.int64) * 10
    rs[B - 2:B + 1] = rs[B - 2]                  # a run of three: the last two of block 0, the first of block 1
    R = int(rs[B - 2])
    vidx = np.arange(len(rs))
    buf, first = cat.encode_rsid(rs, vidx, np.ones(len(rs), dtype=np.int64), "A" * 32)
    assert first[1] == R and first[0] < R
    read = lambda off, ln: buf[off:off + ln]                             # noqa: E731
    got = cat.rsid_lookup(read, first, len(rs), R)
    want = sorted((1, int(v)) for v in vidx[rs == R])
    assert sorted(got) == want and len(got) == 3, got
    # the naive rule (last block whose first number <= rs) finds only the record in block 1
    b = int(np.searchsorted(first, R, side="right")) - 1
    rec = np.frombuffer(buf[qs.HEADER_LEN + 12 * b * B:], cat.RSID_DTYPE)
    assert int((rec["rs_number"] == R).sum()) == 1
    assert cat.rsid_lookup(read, first, len(rs), 10) == [(1, 0)]
    assert cat.rsid_lookup(read, first, len(rs), int(rs[-1])) == [(1, len(rs) - 1)]
    assert cat.rsid_lookup(read, first, len(rs), 5) == [] and cat.rsid_lookup(read, first, len(rs), 10 ** 9) == []


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
