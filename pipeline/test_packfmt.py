"""Round-trip tests for the qtlb reference codec (pipeline/packfmt.py, SPEC.md v0).

    uv run python -m pipeline.test_packfmt

Plain asserts, no test framework. Synthetic cases first, then a real-data round trip for FLNC
(chr7) and SYNPO2L (chr10) from the derived bucket mirror (read only).
"""
from __future__ import annotations

import math
import struct
import sys
import time

import numpy as np
import pyarrow as pa

from . import packfmt as pf
from .common import variants_sql

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


def raises(rule: str, fn, *args, **kwargs) -> None:
    try:
        fn(*args, **kwargs)
    except ValueError as e:
        assert rule in str(e), f"expected {rule!r} in error, got {e!r}"
        return
    raise AssertionError(f"no ValueError raised (expected {rule!r})")


NONSNP = [("AT", "A"), ("A", "AT"), ("N", "A"), ("*", "G"), ("A", "A"), ("C", "N"), ("T", "*"),
          ("ACGT" * 15, "A"), ("C", "G" * 52), ("TT", "T"), ("G", "GACGTACGTA")]


def synth_variants(n: int, seed: int = 0) -> dict:
    """n unique variants sorted by (position, A1, A2), with nulls sprinkled in every value column."""
    rng = np.random.default_rng(seed)
    snps = list(pf.SNP_CODES)
    keys = set()
    pos = 10_000
    while len(keys) < n:
        pos += int(rng.integers(0, 40))           # zero steps give several variants at one position
        pair = snps[rng.integers(len(snps))] if rng.random() < 0.9 else NONSNP[rng.integers(len(NONSNP))]
        keys.add((pos, *pair))
    keys = sorted(keys)[:n]
    rs = rng.integers(1, 1_600_000_000, n).astype(object)
    rs[rng.random(n) < 0.05] = None
    af = rng.random(n)
    af[rng.random(n) < 0.03] = np.nan
    ms = rng.integers(0, 517, n).astype(object)
    ms[rng.random(n) < 0.02] = None
    mc = rng.integers(0, 517, n).astype(object)
    mc[rng.random(n) < 0.02] = None
    match = [["none", "exact", "position"][i] for i in rng.integers(0, 3, n)]
    return {"position": np.array([k[0] for k in keys]), "rs_number": rs, "af": af, "ma_samples": ms,
            "ma_count": mc, "A1": [k[1] for k in keys], "A2": [k[2] for k in keys], "match": match}


def check_variants_equal(dec: dict, src: dict, sl: slice) -> None:
    assert np.array_equal(dec["position"].astype(np.int64), np.asarray(src["position"][sl], dtype=np.int64))
    assert dec["A1"] == list(src["A1"][sl]) and dec["A2"] == list(src["A2"][sl])
    want_rs = [0 if v is None else int(v) for v in src["rs_number"][sl]]
    assert dec["rs_number"].filled(0).astype(np.int64).tolist() == want_rs
    af = np.asarray(src["af"][sl], dtype=np.float64)
    assert np.array_equal(np.isnan(dec["af"]), np.isnan(af))
    ok = ~np.isnan(af)
    assert np.all(np.abs(dec["af"][ok] - af[ok]) <= 0.5 / pf.AF_MAXQ + 1e-15)
    for k in ("ma_samples", "ma_count"):
        want = [None if v is None else int(v) for v in src[k][sl]]
        got = [None if m else int(v) for v, m in zip(dec[k].data, np.ma.getmaskarray(dec[k]))]
        assert got == want, k
    assert dec["match"] == list(src["match"][sl])


# ---- file header ------------------------------------------------------------------------------
@case
def file_header_round_trip():
    for kind, ps in ((pf.KIND_VARIANTS, 1024), (pf.KIND_EQTL, 0), (pf.KIND_SQTL, 0), (pf.KIND_GWAS, 2048), (pf.KIND_GWAS_INDEX, 2048),
                     (pf.KIND_TRANS, 0)):
        n_cis = 120_000 if kind == pf.KIND_VARIANTS else None
        b = pf.file_header(kind, "chrX", 123_456, ps, n_cis)
        assert len(b) == 32 and b[:4] == b"QTLB" and b[16:32] == (123_456).to_bytes(4, "little") + ps.to_bytes(4, "little") + (n_cis or 0).to_bytes(4, "little") + bytes(4)
        h = pf.parse_file_header(b)
        assert h == {"kind": kind, "version": 0, "header_len": 32, "chrom": "chrX", "count": 123_456, "page_size": ps, "n_cis": n_cis}
    assert pf.parse_file_header(pf.file_header(pf.KIND_TRANS, "chrM", 16, 0))["chrom"] == "chrM"
    b = bytearray(pf.file_header(1, "chr7", 1, 1024, 1))
    raises("magic", pf.parse_file_header, b"XTLB" + bytes(b[4:]))
    raises("version", pf.parse_file_header, bytes(b[:5]) + b"\x01" + bytes(b[6:]))
    raises("reserved", pf.parse_file_header, bytes(b[:31]) + b"\x01")
    raises("exceeds count", pf.parse_file_header, bytes(b[:24]) + (2).to_bytes(4, "little") + bytes(4))
    raises("reserved", pf.parse_file_header, pf.file_header(2, "chr1", 1, 0)[:24] + (1).to_bytes(4, "little") + bytes(4))
    raises("needs n_cis", pf.file_header, 1, "chr1", 5, 512)
    raises("needs n_cis", pf.file_header, 1, "chr1", 5, 512, 6)
    raises("only for a variants file", pf.file_header, 6, "chr1", 5, 0, 3)
    raises("1-8 ASCII", pf.file_header, 1, "chr_toolong", 1, 1024, 1)
    raises("page size must be 0", pf.file_header, 2, "chr1", 1, 1024)
    raises("page size must be 0", pf.file_header, 3, "chr1", 1, 512)
    raises("page size must be 0", pf.file_header, 6, "chr1", 1, 512)
    raises("page size 0 not in", pf.file_header, pf.KIND_GWAS, "chr1", 1, 0)
    raises("unknown kind", pf.file_header, 10, "chr1", 1, 0)   # 1-9 are taken; the next free number
    raises("unknown kind", pf.parse_file_header, bytearray(pf.file_header(1, "chr1", 1, 512, 0))[:4] + bytes([10]) + pf.file_header(1, "chr1", 1, 512, 0)[5:])


# ---- variant pages ------------------------------------------------------------------------------
@case
def snp_codes_all_twelve():
    pairs = list(pf.SNP_CODES)
    n = len(pairs)
    for codec in ("raw", "zstd"):
        page = pf.encode_variant_page(7, np.arange(100, 100 + n), [1] * n, [0.5] * n, [1] * n, [1] * n,
                                      [a for a, _ in pairs], [b for _, b in pairs], ["exact"] * n, codec, 19)
        assert len(page) % 4 == 0
        dec = pf.decode_variant_pages(page)
        assert list(zip(dec["A1"], dec["A2"])) == pairs
        assert dec["allele_code"].tolist() == list(range(1, 13))
        assert dec["vidx"].tolist() == list(range(7, 7 + n))
        if codec == "raw":
            assert page[12:16] == bytes(4), "heap_len must be 0 when every record is a SNP"


@case
def non_snp_alleles_in_heap():
    pairs = NONSNP + [("A" * k, "C") for k in range(1, 61)] + [("G", "T" * k) for k in (1, 2, 60)]
    pairs = [p for p in pairs if p not in pf.SNP_CODES]
    n = len(pairs)
    for codec in ("raw", "zstd"):
        page = pf.encode_variant_page(0, np.arange(1, n + 1), [0] * n, [None] * n, [None] * n, [None] * n,
                                      [a for a, _ in pairs], [b for _, b in pairs], ["none"] * n, codec, 3)
        dec = pf.decode_variant_pages(page)
        assert list(zip(dec["A1"], dec["A2"])) == pairs
        assert np.all(dec["allele_code"] == 0)
        if codec == "raw":
            heap_len = int.from_bytes(page[12:16], "little")
            heap = page[12 + 4 + 16 * n:12 + 4 + 16 * n + heap_len]
            assert heap == "".join(f"{a}\t{b}\n" for a, b in pairs).encode()


@case
def value_edges():
    n = 6
    rs = [0, None, 4_000_000_000, 1_562_829_120, 1, 2**32 - 1]
    af = [0.0, 1.0, None, float("nan"), 0.5, 7.6e-6]
    ms = [0, 516, None, 65534, 3, 0]
    mc = [None, 0, 516, 1, 65534, 2]
    match = ["none", "exact", "position", "none", "exact", "position"]
    src = {"position": np.array([5, 5, 6, 7, 7, 8]), "rs_number": np.array(rs, dtype=object), "af": af,
           "ma_samples": np.array(ms, dtype=object), "ma_count": np.array(mc, dtype=object),
           "A1": ["A", "C", "G", "A", "AT", "T"], "A2": ["C", "A", "T", "G", "A", "C"], "match": match}
    for codec in ("raw", "zstd"):
        page = pf.encode_variant_page(0, src["position"], src["rs_number"], np.array([np.nan if v is None else v for v in af]),
                                      src["ma_samples"], src["ma_count"], src["A1"], src["A2"], match, codec, 19)
        dec = pf.decode_variant_pages(page)
        assert dec["rs_number"].mask.tolist() == [True, True, False, False, False, False]
        assert dec["rs_number"].data[2] == 4_000_000_000 and dec["rs_number"].data[5] == 2**32 - 1
        assert dec["af"][0] == 0.0 and dec["af"][1] == 1.0 and math.isnan(dec["af"][2]) and math.isnan(dec["af"][3])
        src_af = {**src, "af": np.array([np.nan if v is None else v for v in af])}
        check_variants_equal(dec, src_af, slice(None))


@case
def pages_ranges_both_codecs():
    n, P = 3000, 1024
    src = synth_variants(n, seed=1)
    for codec in ("raw", "zstd"):
        buf, offs = pf.encode_variants_file("chr7", src["position"], src["rs_number"], src["af"], src["ma_samples"],
                                            src["ma_count"], src["A1"], src["A2"], src["match"], P, codec, 19)
        assert len(offs) == 4 and offs[-1] == len(buf) and np.all(offs % 4 == 0)
        whole = pf.decode_variants_file(buf)
        assert [p["n"] for p in whole["pages"]] == [1024, 1024, 952], "last page shorter than the page size"
        check_variants_equal(whole, src, slice(None))
        # a gene whose range crosses the page 0 / page 1 boundary
        vs, nv = 1000, 100
        off, ln = pf.variant_range(offs, P, vs, nv)
        assert (off, ln) == (offs[0], offs[2] - offs[0])
        pos = src["position"]
        dec = pf.decode_variant_pages(buf[off:off + ln], var_start=vs, n_var=nv,
                                      pos_first=int(pos[vs]), pos_last=int(pos[vs + nv - 1]))
        assert [p["first_vidx"] for p in dec["pages"]] == [0, 1024]
        check_variants_equal(dec, src, slice(0, 2048))
        # a gene inside the short last page
        off2, ln2 = pf.variant_range(offs, P, 2100, 50)
        dec2 = pf.decode_variant_pages(buf[off2:off2 + ln2], var_start=2100, n_var=50)
        assert len(dec2["pages"]) == 1 and dec2["pages"][0]["n"] == 952
        # expectation failures
        raises("does not hold var_start", pf.decode_variant_pages, buf[off:off + ln], var_start=1030, n_var=10)
        raises("does not hold var_start + n_var - 1", pf.decode_variant_pages, buf[off:off + ln], var_start=1000, n_var=1500)
        raises("!= pos_first", pf.decode_variant_pages, buf[off:off + ln], var_start=vs, n_var=nv, pos_first=int(pos[vs]) + 1)
        # structural failures
        p0 = buf[offs[0]:offs[1]]
        p2 = buf[offs[2]:offs[3]]
        raises("does not follow the previous page", pf.decode_variant_pages, p0 + p2)
        raises("runs past the end", pf.decode_variant_pages, buf[off:off + ln - 4])
        bad = bytearray(p0)
        bad[11] = 1
        raises("reserved header byte", pf.decode_variant_pages, bytes(bad))
        bad = bytearray(p0)
        bad[8:10] = (1023).to_bytes(2, "little")
        raises("payload", pf.decode_variant_pages, bytes(bad))
    # padding and mixed codecs
    raw = pf.encode_variant_page(0, [1, 2, 3], [1, 2, 3], [0.1] * 3, [1] * 3, [1] * 3, ["AT"] * 3, ["A", "C", "G"], ["none"] * 3, "raw", 1)
    stored_len = int.from_bytes(raw[:4], "little")
    assert (12 + stored_len) % 4 != 0, "test page must carry padding"
    bad = bytearray(raw)
    bad[-1] = 7
    raises("padding is not zero", pf.decode_variant_pages, bytes(bad))
    z = pf.encode_variant_page(3, [4, 5], [1, 2], [0.1] * 2, [1] * 2, [1] * 2, ["A", "A"], ["C", "G"], ["none"] * 2, "zstd", 1)
    raises("differs from the first page's codec", pf.decode_variant_pages, raw + z)
    bad = bytearray(raw)
    bad[12 + 4 + 16 * 3 - 1] = 8           # flags byte of the last record: bit 3 set
    raises("flags bits 3-7", pf.decode_variant_pages, bytes(bad))
    bad[12 + 4 + 16 * 3 - 1] = 4           # bit 2 on a record that still has a heap record
    raises("heap holds 3 records for 2", pf.decode_variant_pages, bytes(bad))


@case
def encoder_rejects_bad_input():
    ok = dict(position=[1, 2], rs_number=[1, 2], af=[0.1, 0.2], ma_samples=[1, 2], ma_count=[1, 2],
              A1=["A", "A"], A2=["C", "G"], match=["none", "exact"])

    def enc(**kw):
        a = {**ok, **kw}
        return pf.encode_variant_page(0, a["position"], a["rs_number"], a["af"], a["ma_samples"], a["ma_count"],
                                      a["A1"], a["A2"], a["match"], "raw", 1)
    enc()
    raises("decreasing", enc, position=[2, 1])
    raises("rs_number", enc, rs_number=[1, 2**32])
    raises("ma_samples", enc, ma_samples=[65535, 1])
    raises("ma_count", enc, ma_count=[-1, 1])
    raises("af", enc, af=[1.5, 0.1])
    raises("not ASCII without tab", enc, A1=["A\tT", "A"])
    raises("not ASCII without tab", enc, A2=["C", "é"])
    raises("not a string", enc, A1=[None, "A"])
    raises("unknown value", enc, match=["fuzzy", "none"])
    m = 65536
    raises("must be 1..65535", pf.encode_variant_page, 0, np.arange(1, m + 1), np.zeros(m), np.zeros(m), np.zeros(m),
           np.zeros(m), ["A"] * m, ["C"] * m, ["none"] * m, "raw", 1)
    raises("strict (A1, A2) order", pf.encode_variants_file, "chr1", [5, 5], [0, 0], [0.1, 0.1], [1, 1], [1, 1],
           ["A", "A"], ["C", "C"], ["none", "none"], 1024, "raw", 1)
    raises("strict (A1, A2) order", pf.encode_variants_file, "chr1", [5, 5], [0, 0], [0.1, 0.1], [1, 1], [1, 1],
           ["C", "A"], ["A", "C"], ["none", "none"], 1024, "raw", 1)
    raises("null allele", pf.encode_variants_file, "chr1", [5, 6], [0, 0], [0.1, 0.1], [1, 1], [1, 1],
           [None, "A"], [None, "C"], ["none", "none"], 1024, "raw", 1)


def trans_only_rows(n: int, seed: int = 7, start: int = 50) -> dict:
    """Trans-only variants: one per position, a quarter without alleles, some non-SNPs, no counts."""
    rng = np.random.default_rng(seed)
    pos = start + np.cumsum(rng.integers(1, 30, n))
    snps = list(pf.SNP_CODES)
    A1, A2 = [], []
    for _ in range(n):
        r = rng.random()
        pair = (None, None) if r < 0.25 else NONSNP[rng.integers(len(NONSNP))] if r < 0.35 else snps[rng.integers(len(snps))]
        A1.append(pair[0])
        A2.append(pair[1])
    rs = rng.integers(1, 2_000_000_000, n).astype(object)
    rs[rng.random(n) < 0.05] = None
    return {"position": pos, "rs_number": rs, "af": rng.uniform(0.01, 0.99, n).astype(np.float32), "A1": A1, "A2": A2,
            "match": [["none", "exact", "position"][i] for i in rng.integers(0, 3, n)]}


@case
def variants_file_trans_only_section():
    P, n, m = 512, 1300, 700                   # the last cis page is short (276), then 512 + 188 trans-only records
    src = synth_variants(n, seed=4)
    tr = trans_only_rows(m)
    args = (src["position"], src["rs_number"], src["af"], src["ma_samples"], src["ma_count"], src["A1"], src["A2"], src["match"], P)
    for codec in ("raw", "zstd"):
        cis_buf, cis_offs = pf.encode_variants_file("chr3", *args, codec, 19)
        buf, offs = pf.encode_variants_file("chr3", *args, codec, 19, trans_only=tr)
        h = pf.parse_file_header(buf)
        assert (h["count"], h["n_cis"]) == (n + m, n) and pf.parse_file_header(cis_buf)["n_cis"] == n
        assert len(offs) == 6 and np.array_equal(offs[:4], cis_offs) and offs[-1] == len(buf)
        assert buf[32:cis_offs[-1]] == cis_buf[32:], "the trans-only section leaves every cis byte in place"
        whole = pf.decode_variants_file(buf)
        assert whole["n_cis"] == n and [p["n"] for p in whole["pages"]] == [512, 512, 276, 512, 188]
        assert [p["first_vidx"] for p in whole["pages"]] == [0, 512, 1024, 1300, 1812]
        check_variants_equal({k: whole[k][:n] for k in ("position", "A1", "A2", "rs_number", "af", "ma_samples", "ma_count", "match")},
                             src, slice(0, n))
        t = slice(n, n + m)
        assert np.array_equal(whole["position"][t].astype(np.int64), tr["position"])
        assert whole["A1"][t] == tr["A1"] and whole["A2"][t] == tr["A2"] and whole["match"][t] == tr["match"]
        assert whole["no_alleles"][t].tolist() == [a is None for a in tr["A1"]] and not whole["no_alleles"][:n].any()
        assert np.all(whole["allele_code"][t][whole["no_alleles"][t]] == 0)
        assert whole["rs_number"].filled(0)[t].astype(np.int64).tolist() == [0 if v is None else int(v) for v in tr["rs_number"]]
        assert np.all(np.abs(whole["af"][t] - tr["af"].astype(np.float64)) <= 0.5 / pf.AF_MAXQ + 1e-12)
        assert whole["ma_samples"].mask[t].all() and whole["ma_count"].mask[t].all(), "trans-only counts are null"
        # a trans-only range decodes alone; a range across both sections needs n_cis
        tr_only = pf.decode_variant_pages(buf[offs[3]:offs[5]])
        assert tr_only["A1"] == tr["A1"] and tr_only["vidx"][0] == n
        both = buf[offs[2]:offs[4]]
        raises("positions decrease across pages", pf.decode_variant_pages, both)
        assert len(pf.decode_variant_pages(both, n_cis=n)["position"]) == 276 + 512
        raises("reaches n_cis", pf.decode_variant_pages, buf[offs[2]:offs[3]], var_start=1250, n_var=51, n_cis=n)
        pf.decode_variant_pages(buf[offs[2]:offs[3]], var_start=1200, n_var=100, n_cis=n)
        raises("straddles n_cis", pf.decode_variant_pages, buf[offs[0]:offs[1]], n_cis=100)
        raises("on a cis record", pf.decode_variant_pages, buf[offs[3]:offs[4]], n_cis=n + 1000)
        # an empty trans-only section is the cis-only file
        empty = {k: v[:0] for k, v in tr.items()}
        eb, eo = pf.encode_variants_file("chr3", *args, codec, 19, trans_only=empty)
        assert eb == cis_buf and np.array_equal(eo, cis_offs)
    # bit 2 rules, on a raw page
    page = pf.encode_variant_page(0, [7, 9], [1, 0], [0.5, 0.25], [None, None], [None, None], [None, "AT"], [None, "A"],
                                  ["position", "none"], "raw", 1)
    dec = pf.decode_variant_pages(page)
    assert dec["A1"] == [None, "AT"] and dec["A2"] == [None, "A"] and dec["no_alleles"].tolist() == [True, False]
    assert page[12 + 4 + 16 * 2 - 2] == 0x04 | 2 and page[12 + 4 + 14 * 2] == 0, "flags bit 2 plus match code, allele code 0"
    assert page[12 + 4 + 16 * 2:12 + 4 + 16 * 2 + 5] == b"AT\tA\n", "no heap record for the record without alleles"
    bad = bytearray(page)
    bad[12 + 4 + 14 * 2] = 2                  # allele code on the bit-2 record
    raises("nonzero allele code", pf.decode_variant_pages, bytes(bad))
    raises("not a string", pf.encode_variant_page, 0, [7], [1], [0.5], [None], [None], [None], ["A"], ["none"], "raw", 1)
    bad_tr = {**tr, "position": tr["position"].copy()}
    bad_tr["position"][5] = bad_tr["position"][4]
    raises("strictly increasing", pf.encode_variants_file, "chr3", *args, "raw", 1, trans_only=bad_tr)


# ---- trans frames (SPEC section 12) --------------------------------------------------------------
DOF_E, DOF_S = 435, 480
CHR_CODE = {c: i + 1 for i, c in enumerate(pf.TRANS_VARIANT_CHROMS)}


def trans_source(spec: list, seed: int = 9) -> list[dict]:
    """One gene's source rows: spec is (qtl_type, intron (start, end, cluster, strand) or None, variant_chr, positions).
    p from -log10 p 5 to 60, beta = sign * se * t (float32), r2 = t^2 / (t^2 + dof), as the Zenodo trans files relate them."""
    from scipy.special import stdtrit
    rng = np.random.default_rng(seed)
    rows = []
    for qt, intron, chrom, positions in spec:
        dof = DOF_E if qt == "e" else DOF_S
        for x in positions:
            p = 10.0 ** -rng.uniform(5.0, 60.0)
            t = -stdtrit(dof, p / 2)
            se = rng.uniform(0.02, 0.3)
            rows.append({"qtl_type": qt, "intron": intron, "variant_chr": chrom, "position": int(x),
                         "rs_number": None if rng.random() < 0.1 else int(rng.integers(1, 2**32 - 1)),
                         "af": float(np.float32(rng.uniform(0.01, 0.99))), "pval": p,
                         "beta": float(np.float32((1 if rng.random() < 0.5 else -1) * se * t)),
                         "beta_se": float(np.float32(se)), "r2": float(np.float32(t * t / (t * t + dof)))})
    rng.shuffle(rows)                          # the encoder sorts
    return rows


def encode_rows(rows: list[dict], level: int = 19) -> bytes:
    col = lambda k: [r[k] for r in rows]       # noqa: E731
    intr = [r["intron"] or (None, None, None, None) for r in rows]
    return pf.encode_trans_frame(col("qtl_type"), col("variant_chr"), col("position"), col("rs_number"), col("af"), col("pval"),
                                 col("beta"), [x[0] for x in intr], [x[1] for x in intr], [x[2] for x in intr], [x[3] for x in intr], level)


def frame_order(rows: list[dict]) -> list[dict]:
    return sorted(rows, key=lambda r: (r["qtl_type"], r["intron"] or (0, 0, 0, ""), CHR_CODE[r["variant_chr"]], r["position"]))


def check_trans(dec: dict, rows: list[dict], gene_chr: str = "chr7", gene_id: str = "ENSG00000128591", version: int = 16) -> None:
    want = frame_order(rows)
    assert dec["n_e"] + dec["n_s"] == len(want) and dec["n_e"] == sum(r["qtl_type"] == "e" for r in want)
    assert dec["qtl_type"] == [r["qtl_type"] for r in want]
    assert [pf.TRANS_VARIANT_CHROMS[c - 1] for c in dec["variant_chr"].tolist()] == [r["variant_chr"] for r in want]
    assert dec["position"].tolist() == [r["position"] for r in want]
    assert dec["rs_number"].tolist() == [r["rs_number"] or 0 for r in want]
    assert pf.trans_phenotype_ids(dec, gene_chr, gene_id, version) == [
        gene_id if r["qtl_type"] == "e" else f"{gene_chr}:{r['intron'][0]}:{r['intron'][1]}:clu_{r['intron'][2]}_{r['intron'][3]}:{gene_id}.{version}"
        for r in want]
    src = {k: np.array([r[k] for r in want], dtype=np.float64) for k in ("af", "pval", "beta", "beta_se", "r2")}
    assert np.all(np.abs(dec["af"] - src["af"]) <= 0.5 / pf.AF_MAXQ + 1e-12)
    assert np.all(np.abs(dec["nlp"] + np.log10(src["pval"])) <= dec["nlp_max"] / 131066 + 1e-9)
    assert np.all(np.abs(dec["beta"] - src["beta"]) <= dec["beta_max"] / 65534 + 1e-12)
    assert np.all(np.abs(dec["beta_se"] - src["beta_se"]) / src["beta_se"] <= 0.01)
    assert np.all(np.abs(dec["r2"] - src["r2"]) <= 0.002)
    assert int(dec["nlp_code"].max()) == pf.NLP_MAXQ and int(np.abs(dec["beta_code"].astype(np.int32)).max()) == pf.BETA_MAXQ


def trans_payload(frame: bytes) -> bytearray:
    return bytearray(pf.zstd_unframe(frame, None, "test"))


def trans_offsets(n: int, n_s: int, k: int) -> dict:
    r = 32 + 13 * k + pf.pad4(k)
    c = r + 14 * n + pf.pad4(14 * n)
    return {"strand_pad": 32 + 13 * k, "position": r, "rs": r + 4 * n, "af": r + 8 * n, "nlp": r + 10 * n, "beta": r + 12 * n,
            "chr": c, "intron": c + n, "end": c + n + n_s + pf.pad4(n + n_s)}


@case
def trans_frame_both_types_and_chr_change():
    i1, i2 = (1000, 2000, 5, "+"), (900, 5000, 7, "-")
    spec = [("e", None, "chr1", [100, 200]), ("e", None, "chr2", [50, 60]),
            ("s", i1, "chr1", [300]), ("s", i1, "chrX", [10, 20]), ("s", i2, "chr22", [5])]
    rows = trans_source(spec)
    frame = encode_rows(rows)
    fp = pf.zstandard.get_frame_parameters(frame)
    assert fp.has_checksum and fp.content_size > 0
    dec = pf.decode_trans_frame(frame, DOF_E, DOF_S)
    assert (dec["n_e"], dec["n_s"], dec["k"]) == (4, 4, 2)
    assert dec["intron_start"].tolist() == [900, 1000] and dec["strand"] == ["-", "+"] and dec["cluster"].tolist() == [7, 5]
    assert dec["intron"].tolist() == [0, 1, 1, 1]
    check_trans(dec, rows)
    raw = trans_payload(frame)
    o = trans_offsets(8, 4, 2)
    assert len(raw) == o["end"] == 184 and raw[:4] == b"QTT0"
    # absolute at each run start and each chromosome change, else the difference
    assert np.frombuffer(bytes(raw), "<u4", 8, o["position"]).tolist() == [100, 100, 50, 10, 5, 300, 10, 10]
    assert list(raw[o["chr"]:o["chr"] + 8]) == [1, 1, 2, 2, 22, 1, 23, 23]
    assert pf.encode_trans_frame(*[[r[k] for r in rows] for k in ("qtl_type", "variant_chr", "position", "rs_number", "af", "pval", "beta")],
                                 *[[(r["intron"] or (None,) * 4)[j] for r in rows] for j in range(4)], 19) == frame, "deterministic"

    def mutated(off: int, value: bytes) -> bytes:
        x = bytearray(raw)
        x[off:off + len(value)] = value
        return pf.zstd_frame(bytes(x), 1)
    raises("padding", pf.decode_trans_frame, mutated(o["strand_pad"], b"\x01"), DOF_E, DOF_S)
    raises("intron index", pf.decode_trans_frame, mutated(o["intron"] + 3, b"\x02"), DOF_E, DOF_S)
    raises("nlp code", pf.decode_trans_frame, mutated(o["nlp"] + 2, (65534).to_bytes(2, "little")), DOF_E, DOF_S)
    raises("af code", pf.decode_trans_frame, mutated(o["af"], (65535).to_bytes(2, "little")), DOF_E, DOF_S)
    raises("beta code", pf.decode_trans_frame, mutated(o["beta"], (-32768).to_bytes(2, "little", signed=True)), DOF_E, DOF_S)
    raises("variant_chr", pf.decode_trans_frame, mutated(o["chr"], b"\x18"), DOF_E, DOF_S)
    raises("variant_chr decreases", pf.decode_trans_frame, mutated(o["chr"] + 3, b"\x01"), DOF_E, DOF_S)
    raises("position", pf.decode_trans_frame, mutated(o["position"] + 4, bytes(4)), DOF_E, DOF_S)
    raises("never decrease", pf.decode_trans_frame, mutated(o["intron"], b"\x01\x00"), DOF_E, DOF_S)
    raises("magic", pf.decode_trans_frame, mutated(0, b"QTT1"), DOF_E, DOF_S)
    raises("reserved", pf.decode_trans_frame, mutated(14, b"\x01"), DOF_E, DOF_S)
    raises("k", pf.decode_trans_frame, mutated(12, bytes(2)), DOF_E, DOF_S)
    raises("intron table", pf.decode_trans_frame, mutated(32, (1000).to_bytes(4, "little")), DOF_E, DOF_S)
    raises("length", pf.decode_trans_frame, pf.zstd_frame(bytes(raw[:-4]), 1), DOF_E, DOF_S)
    raises("checksum", pf.decode_trans_frame, pf.zstandard.ZstdCompressor().compress(bytes(raw)), DOF_E, DOF_S)


@case
def trans_frame_one_type():
    rows = trans_source([("e", None, "chr5", [10, 20, 30])], seed=2)
    frame = encode_rows(rows)
    dec = pf.decode_trans_frame(frame, DOF_E, DOF_S)
    assert (dec["n_e"], dec["n_s"], dec["k"]) == (3, 0, 0) and dec["intron"].size == 0
    assert len(trans_payload(frame)) == trans_offsets(3, 0, 0)["end"] == 80
    check_trans(dec, rows)
    rows = trans_source([("s", (5, 50, 3, "-"), "chr9", [7, 8]), ("s", (5, 60, 3, "-"), "chr9", [7])], seed=3)
    dec = pf.decode_trans_frame(encode_rows(rows), DOF_E, DOF_S)
    assert (dec["n_e"], dec["n_s"], dec["k"]) == (0, 3, 2) and dec["intron"].tolist() == [0, 0, 1]
    assert dec["position"].tolist() == [7, 8, 7], "a new run starts from an absolute position"
    check_trans(dec, rows, gene_chr="chrM", gene_id="ENSG00000198888", version=2)


@case
def trans_frame_encoder_rejects():
    rows = trans_source([("e", None, "chr1", [10, 20]), ("s", (1, 9, 1, "+"), "chr2", [5])], seed=4)
    encode_rows(rows)
    raises("p must be in (0, 1]", encode_rows, [{**rows[0], "pval": 0.0}] + rows[1:])
    raises("p must be in (0, 1]", encode_rows, [{**rows[0], "pval": float("nan")}] + rows[1:])
    raises("af", encode_rows, [{**rows[0], "af": None}] + rows[1:])
    raises("strictly increasing", encode_rows, rows + [{**frame_order(rows)[0]}])
    raises("intron", encode_rows, [{**r, "intron": None} if r["qtl_type"] == "s" else r for r in rows])
    raises("qtl_type", encode_rows, [{**rows[0], "qtl_type": "t"}] + rows[1:])
    raises("variant_chr", encode_rows, [{**rows[0], "variant_chr": "chrY"}] + rows[1:])
    raises("at least one row", encode_rows, [])
    many = [{**rows[0], "qtl_type": "s", "intron": (i, i + 5, i, "+")} for i in range(256)]
    raises("255", encode_rows, many)
    encode_rows(many[:255])


# ---- gene blocks ---------------------------------------------------------------------------------
DETAILS = {"v": 0, "gene": {"gene_id": "ENSG0", "symbol": "TEST", "slope": float("nan"), "tested": True, "note": "é"},
           "exons": [[1, 2]], "splice": []}
DETAILS_JSON = {"v": 0, "gene": {"gene_id": "ENSG0", "symbol": "TEST", "slope": None, "tested": True, "note": "é"},
                "exons": [[1, 2]], "splice": []}


def dummy_variants(first: int, n: int) -> dict:
    return pf.decode_variant_pages(pf.encode_variant_page(first, np.arange(1000, 1000 + n), np.arange(1, n + 1), np.full(n, 0.25),
                                                          np.full(n, 10), np.full(n, 12), ["A"] * n, ["G"] * n, ["exact"] * n, "raw", 1))


def p_of(slope, se, dof: int) -> np.ndarray:
    """Two-sided Student-t p for the source relation |slope| / se = t (NaN in, NaN out)."""
    from scipy.special import stdtr
    return 2.0 * stdtr(dof, -np.abs(np.asarray(slope, dtype=np.float64)) / np.asarray(se, dtype=np.float64))


def within_bounds(d: dict, p, slope, se, dof: int) -> tuple[float, float]:
    """Assert decoded slope_se and slope are within BOUND_FACTOR x the SPEC per-row bound; return the max ratios."""
    se_b, sl_b = pf.error_bounds(p, slope, se, d["nlp_max"], d["lse_min"], d["lse_max"], dof)
    ok = ~np.isnan(se_b)
    se_r = np.abs(d["slope_se"][ok] - np.asarray(se, dtype=np.float64)[ok]) / se_b[ok]
    sl_r = np.abs(d["slope"][ok] - np.asarray(slope, dtype=np.float64)[ok]) / sl_b[ok]
    assert np.all(np.isfinite(d["slope"][ok])) and np.all(np.isfinite(d["slope_se"][ok])), "non-null rows decode finite"
    assert np.all(se_r <= pf.BOUND_FACTOR) and np.all(sl_r <= pf.BOUND_FACTOR), (se_r.max(), sl_r.max())
    return (float(se_r.max()) if se_r.size else 0.0, float(sl_r.max()) if sl_r.size else 0.0)


@case
def block_no_rows():
    b = pf.encode_gene_block(DETAILS, None, None, [], [], [], None, None, None, 19)
    assert len(b) % 4 == 0 and b[:4] == b"QGB0"
    d = pf.decode_gene_block(b, 435, expect_blk_len=len(b), expect_n_var=None, expect_var_start=None)
    assert d["n_rows"] == 0 and d["var_start"] is None and d["details"] == DETAILS_JSON
    assert (d["nlp_max"], d["lse_min"], d["lse_max"]) == (0.0, 0.0, 0.0)
    assert int.from_bytes(b[12:16], "little") == 0xFFFFFFFF
    t = pf.gene_page_rows(d, dummy_variants(0, 4))
    assert t.num_rows == 0 and t.schema.equals(pf.READER_SCHEMA)
    raises("n_rows 0 != search_index n_var 5", pf.decode_gene_block, b, 435, expect_n_var=5)
    raises("rows", pf.encode_gene_block, DETAILS, 3, None, [], [], [], None, None, None, 1)


@case
def block_nulls_zero_one():
    # 0 plain; 1 null row; 2 p = 0 (slope cannot be derived); 3 p = 1, slope 0; 4 steepest; 5 negative;
    # 6 slope null but se present (no sign, so SE is null); 7 p null but slope and se present
    p = np.array([0.5, np.nan, 0.0, 1.0, 1e-20, 0.03, 0.2, np.nan])
    slope = np.array([0.1, np.nan, 0.3, 0.0, 0.4, -0.05, np.nan, 0.2])
    se = np.array([0.15, np.nan, 0.02, 0.07, 0.04, 0.023, 0.05, 0.05])
    b = pf.encode_gene_block(DETAILS, 100, 5000, p, slope, se, [4], [0.9], [1], 19, pos_first=1002, pos_last=1009)
    d = pf.decode_gene_block(b, 435, expect_n_var=8, expect_var_start=100)
    assert d["nlp_code"].tolist()[1:5] == [pf.NLP_NULL, pf.NLP_ZERO, 0, pf.NLP_MAXQ]
    assert d["se_code"][1] == pf.SE_NULL and d["se_code"][6] == pf.SE_NULL
    assert d["se_code"][2] == 0 and d["se_code"][0] == pf.SE_MAXQ, "smallest se holds code 0, largest 32766"
    assert d["se_code"][5] & pf.SE_SIGN and not d["se_code"][0] & pf.SE_SIGN and d["negative"].tolist()[5]
    assert d["lse_min"] == math.log(0.02) and d["lse_max"] == math.log(0.15)
    assert math.isnan(d["pval_nominal"][1]) and d["pval_nominal"][2] == 0.0 and d["pval_nominal"][3] == 1.0
    assert np.all(np.isnan(d["slope"][[1, 2, 6, 7]])), "null row, p = 0, null se, null p give no slope"
    assert np.all(np.isnan(d["slope_se"][[1, 6]])) and abs(d["slope_se"][2] - 0.02) < 1e-15 and abs(d["slope_se"][7] - 0.05) < 1e-5
    assert abs(d["slope"][3]) < 1e-15, "p = 1 derives a slope of about 0, never a blow-up"
    assert np.sign(d["slope"][5]) == -1 and np.sign(d["slope"][0]) == 1
    t = pf.gene_page_rows(d, dummy_variants(98, 12))
    assert t.schema.equals(pf.READER_SCHEMA)
    assert t["pval_nominal"].to_pylist()[1:4] == [None, 0.0, 1.0] and t["pval_nominal"].to_pylist()[7] is None
    sl, se_out = t["slope"].to_pylist(), t["slope_se"].to_pylist()
    assert [sl[i] is None for i in range(8)] == [False, True, True, False, False, False, True, True]
    assert [se_out[i] is None for i in range(8)] == [False, True, False, False, False, False, True, False]
    assert t["tss_distance"].to_pylist() == [x - 5000 for x in range(1002, 1010)]
    assert t["pip"].to_pylist()[4] == np.float32(0.9) and t["cs_id"].to_pylist() == [None] * 4 + [1, None, None, None]
    raises("var_start 100 != search_index var_start 101", pf.decode_gene_block, b, 435, expect_var_start=101)
    raises("pos_first", pf.gene_page_rows, {**d, "pos_first": 999}, dummy_variants(98, 12))
    raises("do not cover", pf.gene_page_rows, d, dummy_variants(102, 12))
    raises("se must be finite and > 0", pf.encode_gene_block, DETAILS, 0, 0, [0.1], [0.1], [0.0], None, None, None, 1,
           pos_first=1, pos_last=1)


@case
def p_one_small_slope_and_null_row():
    """The rows that broke the earlier slope layout: p = 1 or p reading back as 1, slope about 0."""
    dof = 435
    slope = np.array([0.0, -1e-9, 3e-12, -0.002, np.nan, 0.15, -1.2, 0.0, 4e-5])
    se = np.array([0.05, 0.08, 0.2, 0.11, np.nan, 0.04, 0.3, 0.06, 0.09])
    p = p_of(slope, se, dof)
    assert p[0] == 1.0 and p[7] == 1.0 and np.isnan(p[4]) and 1 - 1e-7 < p[1] < 1
    b = pf.encode_gene_block(DETAILS, 0, 0, p, slope, se, None, None, None, 1, pos_first=1000, pos_last=1008)
    d = pf.decode_gene_block(b, dof)
    code0 = d["nlp_code"] == 0
    assert code0[[0, 1, 2, 7]].all(), "p = 1 and p within a half step of 1 read back as exactly 1"
    assert np.all(np.abs(d["slope"][code0]) <= 1e-15), d["slope"][code0]
    assert d["nlp_code"][4] == pf.NLP_NULL and d["se_code"][4] == pf.SE_NULL
    assert d["negative"].tolist() == [False, True, False, True, False, False, True, False, False]
    within_bounds(d, p, slope, se, dof)
    t = pf.gene_page_rows(d, dummy_variants(0, 9))
    assert [t[c].to_pylist()[4] for c in ("pval_nominal", "slope", "slope_se")] == [None, None, None], "null row"
    se_out = np.array(t["slope_se"].to_pylist()[:4] + t["slope_se"].to_pylist()[5:], dtype=np.float64)
    assert np.all(np.isfinite(se_out)) and se_out.max() < 0.31, "no SE blow-up"


@case
def block_all_p_one_single_row_scales():
    b = pf.encode_gene_block(DETAILS, 0, -10, [1.0, 1.0, np.nan], [0.0, -0.0, 0.0], [0.05, 0.05, np.nan], [], [], [], 1,
                             pos_first=1, pos_last=3)
    d = pf.decode_gene_block(b, 435)
    assert d["nlp_max"] == 0.0 and d["nlp_code"].tolist() == [0, 0, pf.NLP_NULL] and d["n_cs"] == 0
    assert d["lse_min"] == d["lse_max"] == math.log(0.05) and d["se_code"].tolist() == [0, 0, pf.SE_NULL]
    assert d["pval_nominal"][:2].tolist() == [1.0, 1.0] and np.all(np.abs(d["slope_se"][:2] - 0.05) < 1e-16)
    assert np.all(np.abs(d["slope"][:2]) < 1e-15), "-0.0 is stored as positive; slope about 0"
    sl1, se1 = -0.7, 0.13
    p1 = p_of([sl1], [se1], 480)
    b1 = pf.encode_gene_block(DETAILS, 9, 0, p1, [sl1], [se1], [0], [1.0], [3], 1, pos_first=42, pos_last=42)
    d1 = pf.decode_gene_block(b1, 480)
    assert d1["nlp_code"].tolist() == [pf.NLP_MAXQ] and d1["se_code"].tolist() == [pf.SE_SIGN]
    assert abs(d1["slope_se"][0] - se1) < 1e-15 and abs(d1["slope"][0] - sl1) < 1e-9
    b2 = pf.encode_gene_block(DETAILS, 0, 0, [0.1, 0.2], [np.nan, np.nan], [np.nan, np.nan], None, None, None, 1, pos_first=1, pos_last=2)
    d2 = pf.decode_gene_block(b2, 435)
    assert (d2["lse_min"], d2["lse_max"]) == (0.0, 0.0) and d2["se_code"].tolist() == [pf.SE_NULL] * 2
    se = np.array([0.01, 0.1, 1.0, 0.1])
    b3 = pf.encode_gene_block(DETAILS, 0, 0, [0.1] * 4, [-0.2, 0.2, 0.3, 0.0], se, None, None, None, 1, pos_first=1, pos_last=4)
    d3 = pf.decode_gene_block(b3, 435)
    assert d3["se_code"].tolist() == [pf.SE_SIGN | 0, 16383, pf.SE_MAXQ, 16383]
    assert np.all(np.abs(np.log(d3["slope_se"]) - np.log(se)) <= (d3["lse_max"] - d3["lse_min"]) / 65532 * (1 + 1e-9))


@case
def block_credible_sets():
    # memberships from SuSiE: row 3 is in two sets; the block keeps both and row output picks best
    raw = [(3, 0.2, 1), (3, 0.7, 2), (8, 0.9, 1), (0, 0.05, 2)]
    raw.sort(key=lambda x: (x[0], x[2]))
    n = 10
    p = np.full(n, 0.01)
    se = np.full(n, 0.04)
    b = pf.encode_gene_block(DETAILS, 0, 0, p, np.full(n, 0.1), se, [x[0] for x in raw], [x[1] for x in raw],
                             [x[2] for x in raw], 1, pos_first=1000, pos_last=1009)
    d = pf.decode_gene_block(b, 435)
    assert d["cs_row"].tolist() == [0, 3, 3, 8] and d["cs_id"].tolist() == [2, 1, 2, 1]
    assert d["cs_pip"].tolist() == np.array([0.05, 0.2, 0.7, 0.9], dtype=np.float32).tolist()
    t = pf.gene_page_rows(d, dummy_variants(0, n))
    assert t["cs_id"].to_pylist() == [2, None, None, 2, None, None, None, None, 1, None]
    raises("strictly ascending", pf.encode_gene_block, DETAILS, 0, 0, p, np.full(n, 0.1), se, [3, 3], [0.2, 0.7], [2, 1], 1,
           pos_first=1000, pos_last=1009)
    raises("cs_row", pf.encode_gene_block, DETAILS, 0, 0, p, np.full(n, 0.1), se, [10], [0.2], [1], 1, pos_first=1000, pos_last=1009)


@case
def sqtl_blocks_kind3():
    # an intron block: no details frame, one null row, row 4 in two credible sets (as LINC01954's 12 variants)
    n = 12
    rng = np.random.default_rng(11)
    se = rng.uniform(0.03, 0.2, n)
    slope = rng.normal(0.0, 0.3, n)
    p = p_of(slope, se, 480)
    p[5] = slope[5] = se[5] = np.nan
    cs = [(0, 0.1, 1), (4, 0.3, 1), (4, 0.6, 2), (11, 0.9, 3)]
    b = pf.encode_gene_block(None, 40, 900, p, slope, se, [x[0] for x in cs], [x[1] for x in cs], [x[2] for x in cs], 19,
                             pos_first=1000, pos_last=1000 + n - 1)
    body = 64 + 4 * n + 12 * len(cs)
    assert len(b) == body + pf.pad4(body) and b[56:64] == bytes(8), "64 + 4 n_rows + 12 n_cs, details_zlen = details_len = 0"
    d = pf.decode_gene_block(b, 480, kind=pf.KIND_SQTL, expect_blk_len=len(b), expect_n_var=n, expect_var_start=40)
    assert d["details"] is None and d["details_zlen"] == 0 and d["details_len"] == 0
    assert d["cs_row"].tolist() == [0, 4, 4, 11] and d["cs_id"].tolist() == [1, 1, 2, 3]
    within_bounds(d, p, slope, se, 480)
    t = pf.gene_page_rows(d, dummy_variants(40, n))
    assert t["cs_id"].to_pylist()[4] == 2 and abs(t["pip"].to_pylist()[4] - 0.6) < 1e-6, "row output keeps the higher-PIP set"
    assert t["tss_distance"].to_pylist() == list(range(100, 100 + n))
    assert t["pval_nominal"].null_count == 1 and t["slope"].null_count == 1 and t["slope_se"].null_count == 1
    # a one-row intron
    b1 = pf.encode_gene_block(None, 0, 1000, [0.03], [0.2], [0.09], [0], [1.0], [0], 1, pos_first=1000, pos_last=1000)
    d1 = pf.decode_gene_block(b1, 480, kind=pf.KIND_SQTL, expect_blk_len=80, expect_n_var=1, expect_var_start=0)
    assert d1["n_rows"] == 1 and d1["n_cs"] == 1 and len(b1) == 80
    # the details rule per kind
    eq = pf.encode_gene_block(DETAILS, 40, 900, p, slope, se, None, None, None, 1, pos_first=1000, pos_last=1000 + n - 1)
    raises("needs a details frame", pf.decode_gene_block, b, 480)
    raises("no details frame", pf.decode_gene_block, eq, 480, kind=pf.KIND_SQTL)
    x = bytearray(b)
    x[60:64] = (5).to_bytes(4, "little")
    raises("no details frame", pf.decode_gene_block, bytes(x), 480, kind=pf.KIND_SQTL)
    raises("at least one row", pf.encode_gene_block, None, None, None, [], [], [], None, None, None, 1)
    raises("not a results kind", pf.decode_gene_block, b, 480, kind=pf.KIND_VARIANTS)
    # a kind 3 file: blocks back to back after the header
    buf, offs = pf.encode_eqtl_file("chr2", [b, b1], pf.KIND_SQTL)
    h, walked = pf.walk_eqtl_file(buf, pf.KIND_SQTL)
    assert h["kind"] == pf.KIND_SQTL and h["count"] == 2 and h["page_size"] == 0
    assert walked == [(int(offs[0]), len(b)), (int(offs[1]), len(b1))] and offs[-1] == len(buf)
    raises("header kind", pf.walk_eqtl_file, buf, pf.KIND_EQTL)


# ---- GWAS pack (SPEC section 11) ------------------------------------------------------------------
def gwas_rows(n: int, seed: int = 3) -> dict:
    """Synthetic GWAS rows at the source's precision: 4-decimal beta/se/eaf, 4-significant-digit p, repeated positions."""
    rng = np.random.default_rng(seed)
    inc = rng.integers(0, 40, n)
    inc[2048:2050] = 0                    # a position repeated across the first block boundary
    pos = np.cumsum(inc) + 10_000
    mant = rng.integers(1000, 10000, n)
    exp = rng.integers(-43, -3, n)
    p = np.array([float(f"{m}e{e}") for m, e in zip(mant, exp)])
    p[:3] = [1.0, 6.396e-40, 0.05]
    pairs = [("A", "G"), ("C", "T"), ("AT", "A"), ("G", "GACGT" * 40), ("T", "C")]
    ea, nea = zip(*[pairs[i] for i in rng.integers(0, len(pairs), n)])
    return {"position": pos, "beta": np.array([float(f"{x:.4f}") for x in rng.normal(0, 0.2, n)]),
            "se": np.array([float(f"{x:.4f}") for x in rng.uniform(0.0201, 3.3318, n)]),
            "eaf": np.array([float(f"{x:.4f}") for x in rng.uniform(0, 1, n)]), "p": p,
            "rs_number": np.where(rng.random(n) < 0.125, 0, rng.integers(1, 2**31, n)),
            "n": rng.choice([42637, 422920, 937963], n), "ea": list(ea), "nea": list(nea)}


@case
def gwas_p_codes_edges():
    p = np.array([1.0, 0.001, 0.5, 6.396e-40, 1e-5, 9.999e-5, 1.234e-19, 1.234e-23, 0.1234, 5e-8])
    m, e = pf.gwas_p_codes(p)
    assert m.tolist()[:4] == [1000, 1000, 5000, 6396] and e.tolist()[:4] == [-3, -6, -4, -43]
    back = pf.gwas_p_value(m, e)
    assert np.all(np.abs(back - p) <= 1e-12 * p)
    assert back[6] == p[6] and back[8] == p[8], "p_exp down to -22 reads back as the source double"
    raises("more than 4 significant digits", pf.gwas_p_codes, [0.12345])
    raises("not in (0, 1]", pf.gwas_p_codes, [0.0])
    raises("not in (0, 1]", pf.gwas_p_codes, [1.5])
    raises("not in (0, 1]", pf.gwas_p_codes, [np.nan])
    raises("more than 4 decimals", pf.gwas_scaled, [0.12345], "beta", -10, 10**6)
    raises("outside 0..65535", pf.gwas_scaled, [6.5536], "se", 0, 65535)
    raises("is null", pf.gwas_scaled, [np.nan], "eaf", 0, 10000)


@case
def gwas_block_and_index_round_trip():
    n, rows_per = 5000, 2048
    r = gwas_rows(n)
    nv = [42637, 422920, 937963]
    codes = pf.gwas_codes(r["position"], r["beta"], r["se"], r["eaf"], r["p"], r["rs_number"], r["n"], nv)
    frames, fp, eo, off = [], [], [], pf.FILE_HEADER_LEN
    for s in range(0, n, rows_per):
        e = min(s + rows_per, n)
        f = pf.encode_gwas_block({k: v[s:e] for k, v in codes.items()}, r["ea"][s:e], r["nea"][s:e], 19)
        frames.append(f)
        off += len(f)
        fp.append(int(r["position"][s]))
        eo.append(off)
    buf = pf.file_header(pf.KIND_GWAS, "chr22", n, rows_per) + b"".join(frames)
    idx = pf.decode_gwas_index(pf.encode_gwas_index(nv, [("chr22", fp, eo)], rows_per, 19))
    assert idx["block_rows"] == rows_per and idx["n_values"] == nv and idx["chroms"]["chr22"][1].tolist() == eo
    got = {k: [] for k in ("position", "ea", "nea", "rs_number", "beta", "se", "eaf", "p", "n")}
    for k in range(len(frames)):
        start = pf.FILE_HEADER_LEN if k == 0 else eo[k - 1]
        d = pf.decode_gwas_block(buf[start:eo[k]], nv, expect_rows=min(rows_per, n - k * rows_per))
        body = zstandard_payload(buf[start:eo[k]])
        assert len(body) == 8 + 21 * d["rows"] + struct.unpack_from("<I", body, 4)[0]
        for c in got:
            got[c].extend(list(d[c]))
    for c in ("position", "ea", "nea", "rs_number", "n", "beta", "se", "eaf"):
        assert list(got[c]) == list(r[c]) if c in ("ea", "nea") else np.array_equal(np.array(got[c]), np.asarray(r[c])), c
    assert np.all(np.abs(np.array(got["p"]) - r["p"]) <= 1e-12 * r["p"])
    # window rule: every row with lo <= position <= hi is inside the range, including a position repeated across a boundary
    fpa, eoa = np.array(fp), np.array(eo)
    dup = next(k for k in range(1, len(fp)) if r["position"][k * rows_per - 1] == fp[k])
    rng = np.random.default_rng(5)
    wins = [(fp[dup], fp[dup] + 50), (fp[dup] - 1, fp[dup]), (1, 5), (int(r["position"][-1]), int(r["position"][-1]) + 9),
            (int(r["position"][0]), int(r["position"][0])), (10**9, 10**9 + 1)]
    wins += [(int(a), int(a) + int(w)) for a, w in zip(rng.integers(9000, int(r["position"][-1]) + 100, 40), rng.integers(0, 60_000, 40))]
    for lo, hi in wins:
        want = np.flatnonzero((r["position"] >= lo) & (r["position"] <= hi))
        w = pf.gwas_window(fpa, eoa, lo, hi)
        if w is None:
            assert want.size == 0 and hi < fp[0]
            continue
        ks, ke, b0, b1 = w
        pos = np.concatenate([pf.decode_gwas_block(buf[(pf.FILE_HEADER_LEN if k == 0 else eo[k - 1]):eo[k]], nv)["position"] for k in range(ks, ke + 1)])
        first_row = ks * rows_per
        sel = np.flatnonzero((pos >= lo) & (pos <= hi)) + first_row
        assert np.array_equal(sel, want), (lo, hi)
        assert b0 == (pf.FILE_HEADER_LEN if ks == 0 else eo[ks - 1]) and b1 == eo[ke]
    raises("lo", pf.gwas_window, fpa, eoa, 10, 9)
    # rules
    raises("not in the n table", pf.gwas_codes, r["position"], r["beta"], r["se"], r["eaf"], r["p"], r["rs_number"], np.full(n, 7), nv)
    raises("decreasing", pf.gwas_codes, r["position"][::-1], r["beta"], r["se"], r["eaf"], r["p"], r["rs_number"], r["n"], nv)
    raises("not an ASCII string", pf.encode_gwas_block, {k: v[:2] for k, v in codes.items()}, ["A\tC", "A"], ["G", "T"], 1)
    raises("at least one row", pf.encode_gwas_block, {k: v[:0] for k, v in codes.items()}, [], [], 1)
    x = pf.zstd_frame(struct.pack("<III", 1, 0, 5) + bytes(17), 1)     # position 5, p_mant 0
    raises("mantissa", pf.decode_gwas_block, x, nv)
    raises("checksum", pf.decode_gwas_index, pf.file_header(pf.KIND_GWAS_INDEX, "all", 0, 2048) + pf.zstandard.ZstdCompressor().compress(b"x" * 8))


def zstandard_payload(frame: bytes) -> bytes:
    return pf.zstd_unframe(frame, None, "test")


@case
def block_validation():
    n = 5
    b = pf.encode_gene_block(DETAILS, 7, 0, np.full(n, 0.2), np.linspace(-1, 1, n), np.linspace(0.05, 0.2, n), [1], [0.5], [1], 19,
                             pos_first=10, pos_last=20)
    pf.decode_gene_block(b, 435, expect_blk_len=len(b), expect_n_var=n, expect_var_start=7)

    def mutate(off, val):
        x = bytearray(b)
        x[off:off + len(val)] = val
        return bytes(x)
    raises("magic", pf.decode_gene_block, mutate(0, b"X"), 435)
    raises("length field", pf.decode_gene_block, b + bytes(4), 435)
    raises("search_index blk_len", pf.decode_gene_block, b, 435, expect_blk_len=len(b) + 4)
    raises("n_rows 5 != search_index n_var 6", pf.decode_gene_block, b, 435, expect_n_var=6)
    raises("n_rows 5 != search_index n_var None", pf.decode_gene_block, b, 435, expect_n_var=None)
    raises("lse_min", pf.decode_gene_block, mutate(40, struct.pack("<d", 5.0)), 435)
    raises("lse_min", pf.decode_gene_block, mutate(48, struct.pack("<d", float("inf"))), 435)
    raises("0x7FFF", pf.decode_gene_block, mutate(64 + 2, b"\xff\x7f"), 435)
    ones = bytearray(b)
    for i in range(n):
        ones[64 + 4 * i + 2:64 + 4 * i + 4] = (1 | (int.from_bytes(b[64 + 4 * i + 2:64 + 4 * i + 4], "little") & pf.SE_SIGN)).to_bytes(2, "little")
    raises("no row holds log(SE) code 0", pf.decode_gene_block, bytes(ones), 435)
    cs_off = 64 + 4 * n
    raises("credible-set padding", pf.decode_gene_block, mutate(cs_off + 9, b"\x01"), 435)
    raises("credible-set row", pf.decode_gene_block, mutate(cs_off, b"\x09"), 435)
    dz_off = cs_off + 12
    dzlen = int.from_bytes(b[56:60], "little")
    raises("block details", pf.decode_gene_block, mutate(dz_off + dzlen - 1, bytes([b[dz_off + dzlen - 1] ^ 0xFF])), 435)
    body = dz_off + dzlen
    if body < len(b):
        raises("padding", pf.decode_gene_block, mutate(len(b) - 1, b"\x01"), 435)
    lowered = bytearray(b)
    for i in range(n):
        lowered[64 + 4 * i + 1] = 0        # high byte of every nlp code: 65533 -> 253
    raises("no row holds code 65533", pf.decode_gene_block, bytes(lowered), 435)


@case
def quantization_bounds():
    rng = np.random.default_rng(3)
    for trial in range(50):
        m = int(rng.integers(1, 5000))
        p = 10.0 ** -rng.uniform(0, rng.uniform(0.01, 299), m)
        p[rng.random(m) < 0.01] = 1.0
        q, nlp_max = pf.quantize_nlp(p)
        x = -np.log10(p)
        step = nlp_max / pf.NLP_MAXQ
        err = np.abs(pf.dequantize_nlp(q, nlp_max) - x)
        assert np.all(err <= step / 2 * (1 + 1e-9) + 1e-300), (trial, err.max(), step)
        se = np.exp(rng.normal(-2.5, rng.uniform(0.01, 1.5), m))
        sl = rng.normal(0, rng.uniform(1e-3, 3), m)
        sl[rng.random(m) < 0.01] = 0.0
        qs, lo, hi = pf.quantize_se(se, sl)
        se_hat, neg = pf.dequantize_se(qs, lo, hi)
        h = (hi - lo) / (2 * pf.SE_MAXQ)
        assert np.all(np.abs(np.log(se_hat) - np.log(se)) <= h * (1 + 1e-9) + 1e-12), trial
        assert np.array_equal(neg, sl < 0) and int((qs & pf.SE_QMASK).max()) <= pf.SE_MAXQ
    assert pf.p_from_nlp(np.array([np.inf]))[0] == 0.0
    raises("exceeds the format limit", pf.quantize_nlp, np.array([1e-301]))
    t = pf.t_from_p(np.array([1.0, 0.0, np.nan, 0.05]), 435)
    assert 0 <= t[0] < 1e-15 and t[1] == np.inf and np.isnan(t[2]) and abs(t[3] - 1.9654) < 1e-3


@case
def reader_schema_matches():
    b = pf.encode_gene_block(DETAILS, 2, 900, [0.1, 0.2], [0.3, -0.3], [0.18, 0.23], [0], [0.5], [1], 1, pos_first=1002, pos_last=1003)
    t = pf.gene_page_rows(pf.decode_gene_block(b, 435), dummy_variants(0, 8))
    assert t.schema.equals(pf.READER_SCHEMA)
    assert [(f.name, str(f.type)) for f in t.schema] == [
        ("position", "int32"), ("A1", "string"), ("A2", "string"), ("rs_number", "int64"), ("tss_distance", "int32"),
        ("af", "float"), ("ma_samples", "int16"), ("ma_count", "int16"), ("pval_nominal", "double"),
        ("slope", "float"), ("slope_se", "float"), ("pip", "float"), ("cs_id", "int8")]


@case
def details_json_rules():
    j = pf.details_bytes({"v": 0, "a": float("nan"), "b": np.float64(0.1), "c": np.int64(3), "d": np.bool_(True), "e": "é"})
    assert j == '{"v":0,"a":null,"b":0.1,"c":3,"d":true,"e":"é"}'.encode("utf-8")
    raises("Out of range float", pf.details_bytes, {"x": float("inf")})


# ---- real data ------------------------------------------------------------------------------------
T_BANDS = [(0.0, 0.01), (0.01, 0.1), (0.1, 1.0), (1.0, math.inf)]


# ---- hits pack, rsID index, variant index (SPEC sections 13 to 15) ----------------------------
def reframe(payload) -> bytes:
    """A mutated hits/variant-index payload back into one legal zstd frame."""
    return pf.zstd_frame(bytes(payload), 3)


def hits_frame(rows: list[dict], first_vidx: int, n_variants: int, level: int = 3) -> bytes:
    """encode_hits_frame from a list of row dicts; every optional column is passed explicitly."""
    g = lambda k, d=None: [r.get(k, d) for r in rows]
    fl = lambda k: np.array([np.nan if r.get(k) is None else float(r[k]) for r in rows], dtype=np.float64)
    return pf.encode_hits_frame(
        first_vidx, n_variants, np.array(g("vidx"), dtype=np.int64), np.array(g("kind"), dtype=np.int64),
        np.array(g("gene"), dtype=np.int64), pval=fl("pval"), beta=fl("beta"), slope=fl("slope"),
        slope_se=fl("slope_se"), pip=fl("pip"), cs_id=np.array(g("cs_id"), dtype=object),
        intron_start=np.array(g("intron_start"), dtype=object), intron_end=np.array(g("intron_end"), dtype=object),
        cluster=np.array(g("cluster"), dtype=object), strand=g("strand"),
        significant=[bool(r.get("significant", False)) for r in rows], level=level)


def intron(s: int, e: int, c: int, strand: str) -> dict:
    return {"intron_start": s, "intron_end": e, "cluster": c, "strand": strand}


SHORT_ROW = [{"vidx": 3072, "kind": 2, "gene": 1, "pval": 1e-9, "slope": 0.2, "slope_se": 0.05}]

HITS_ROWS = [
    # slot 0 (vidx 100): all six kinds, so every value column and both flag bits are exercised
    {"vidx": 100, "kind": 0, "gene": 7, "pval": 1e-30, "beta": 0.90},
    {"vidx": 100, "kind": 0, "gene": 3, "pval": 1e-12, "beta": -0.42},
    {"vidx": 100, "kind": 1, "gene": 9, "pval": 1e-8, "beta": 0.11, **intron(11000, 12000, 5, "+")},
    {"vidx": 100, "kind": 2, "gene": 7, "pval": 1e-5, "slope": 0.55, "slope_se": 0.08, "significant": True},
    {"vidx": 100, "kind": 3, "gene": 9, "pval": 2e-4, "slope": -1.20, "slope_se": 0.30, **intron(11000, 12000, 5, "-")},
    {"vidx": 100, "kind": 4, "gene": 7, "pip": 0.90, "cs_id": 1},
    {"vidx": 100, "kind": 5, "gene": 9, "pip": 0.25, "cs_id": 2, **intron(11000, 12000, 5, "+")},
    # slot 1 (vidx 101) has no rows at all
    # slot 2: one variant in two credible sets of one phenotype, kept as two rows, PIP descending
    {"vidx": 102, "kind": 4, "gene": 3, "pip": 0.70, "cs_id": 1},
    {"vidx": 102, "kind": 4, "gene": 3, "pip": 0.30, "cs_id": 2},
    # slot 3: the largest ord a u16 holds, and a p of exactly 1
    {"vidx": 103, "kind": 0, "gene": 65535, "pval": 1.0, "beta": 0.01},
]


@case
def hits_frame_all_kinds():
    src = HITS_ROWS
    b = hits_frame(src, 100, 4)
    d = pf.decode_hits_frame(b, 100, "t")
    assert d["n_variants"] == 4 and d["n_rows"] == len(src)
    assert d["count"].tolist() == [7, 0, 2, 1] and d["start"].tolist() == [0, 7, 7, 9]
    assert d["vidx"].tolist() == [100] * 7 + [102, 102, 103]
    # rows of one variant go by kind, then p ascending (0-3) or PIP descending (4-5), then ord
    assert d["kind"].tolist() == [0, 0, 1, 2, 3, 4, 5, 4, 4, 0]
    assert d["gene"].tolist() == [7, 3, 9, 7, 9, 7, 9, 3, 3, 65535]
    assert pf.hits_slice(d, 100) == (0, 7) and pf.hits_slice(d, 101) == (7, 7) and pf.hits_slice(d, 103) == (9, 10)
    raises("outside the frame", pf.hits_slice, d, 104)
    # every value comes back inside SPEC section 13's limits
    for i, r in ((0, src[0]), (1, src[1]), (2, src[2]), (9, src[9])):
        assert abs(-math.log10(d["pval"][i]) + math.log10(r["pval"])) <= d["trans_nlp_max"] / 131066 + 1e-9, i
        assert abs(d["beta"][i] - r["beta"]) <= d["trans_beta_max"] / 65534 + 1e-12, i
    for i, r in ((3, src[3]), (4, src[4])):
        assert abs(-math.log10(d["pval"][i]) + math.log10(r["pval"])) <= d["perm_nlp_max"] / 131066 + 1e-9, i
        assert abs(d["slope_se"][i] - r["slope_se"]) <= d["se_max"] / 131070 + 1e-12, i
        assert abs(d["slope"][i] - r["slope"]) <= d["slope_max"] / 65534 + 1e-12, i
    for i, r in ((5, src[5]), (6, src[6]), (7, src[7]), (8, src[8])):
        assert abs(d["pip"][i] - r["pip"]) <= 0.5 / 65535 + 1e-12, i
        assert d["cs_id"][i] == r["cs_id"], i
    # kinds carry only their own values; the others are NaN or -1
    trans, lead, cs = np.array([0, 1, 2, 9]), np.array([3, 4]), np.array([5, 6, 7, 8])
    assert np.all(np.isnan(d["beta"][np.r_[lead, cs]])) and np.all(np.isnan(d["slope"][np.r_[trans, cs]]))
    assert np.all(np.isnan(d["pip"][np.r_[trans, lead]])) and np.all(d["cs_id"][np.r_[trans, lead]] == -1)
    assert np.all(np.isnan(d["pval"][cs])) and np.all(np.isnan(d["slope_se"][np.r_[trans, cs]]))
    # flags: bit 0 is the intron strand, bit 1 the significant lead
    assert d["strand"] == [None, None, "+", None, "-", None, "+", None, None, None]
    assert d["significant"].tolist() == [False, False, False, True, False, False, False, False, False, False]
    assert d["flags"].tolist() == [0, 0, 0, 2, 1, 0, 0, 0, 0, 0], "bit 0 only on a '-' intron, bit 1 only on a significant lead"
    assert d["intron_start"].tolist() == [0, 0, 11000, 0, 11000, 0, 11000, 0, 0, 0]
    assert d["intron_end"].tolist() == [0, 0, 12000, 0, 12000, 0, 12000, 0, 0, 0]
    assert d["cluster"].tolist() == [0, 0, 5, 0, 5, 0, 5, 0, 0, 0]
    assert np.all(d["v2"][trans] == 0) and np.all(d["v3"][cs] == 0), "v2 is 0 on kinds 0-1, v3 on kinds 4-5"
    # every typed column starts at a multiple of 4 and the payload is exactly the layout's length
    off = pf.hits_offsets(4, len(src))
    assert d["payload_len"] == off["end"] and all(off[k] % 4 == 0 for k in ("gene", "v1", "v2", "v3", "intron_start", "intron_end", "cluster"))
    # order is independent of the input order: the encoder sorts
    assert hits_frame(list(reversed(src)), 100, 4) == b


@case
def hits_frame_empty_and_short_last():
    # a frame with no rows at all is valid, and so is a short last frame
    e = pf.decode_hits_frame(pf.encode_hits_frame(2048, 1024, None, None, None, level=3), 2048, "t")
    assert e["n_rows"] == 0 and e["n_variants"] == 1024 and e["count"].sum() == 0
    assert e["trans_nlp_max"] == e["trans_beta_max"] == e["perm_nlp_max"] == e["se_max"] == e["slope_max"] == 0.0
    assert pf.hits_slice(e, 2048) == (0, 0) and e["payload_len"] == pf.hits_offsets(1024, 0)["end"]
    short = hits_frame(SHORT_ROW, 3072, 5)
    s = pf.decode_hits_frame(short, 3072, "t")
    assert s["n_variants"] == 5 and s["n_rows"] == 1 and s["count"].tolist() == [1, 0, 0, 0, 0]
    assert s["trans_nlp_max"] == 0.0 and s["trans_beta_max"] == 0.0 and s["perm_nlp_max"] > 0
    # a frame decoded without first_vidx carries no absolute vidx, and hits_slice then takes the slot
    nv = pf.decode_hits_frame(short, None, "t")
    assert "vidx" not in nv and pf.hits_slice(nv, 0) == (0, 1)


@case
def hits_frame_encoder_rejects():
    ok = [{"vidx": 10, "kind": 2, "gene": 1, "pval": 1e-9, "slope": 0.2, "slope_se": 0.05}]
    bad = lambda **kw: [{**ok[0], **kw}]
    raises("must be 1..65535", hits_frame, ok, 10, 0)
    raises("must be 1..65535", hits_frame, ok, 10, 65536)
    raises("vidx", hits_frame, bad(vidx=99), 10, 5)
    raises("p must be in (0, 1]", hits_frame, bad(pval=0.0), 10, 5)
    raises("p must be in (0, 1]", hits_frame, bad(pval=1.5), 10, 5)
    raises("p must be in (0, 1]", hits_frame, bad(pval=None), 10, 5)
    raises("need a finite slope", hits_frame, bad(slope_se=0.0), 10, 5)
    raises("need a finite slope", hits_frame, bad(slope=None), 10, 5)
    raises("pip must be in [0, 1]", hits_frame, [{"vidx": 10, "kind": 4, "gene": 1, "pip": 1.5, "cs_id": 0}], 10, 5)
    raises("cs_id must be", hits_frame, [{"vidx": 10, "kind": 4, "gene": 1, "pip": 0.5}], 10, 5)
    raises("needs a strand", hits_frame, [{"vidx": 10, "kind": 1, "gene": 1, "pval": 0.1, "beta": 0.1,
                                          "intron_start": 1, "intron_end": 2, "cluster": 3}], 10, 5)
    raises("a strand on a kind 0, 2 or 4 row", hits_frame, bad(strand="+"), 10, 5)
    raises("an intron field on a kind 0, 2 or 4 row", hits_frame, bad(intron_start=5), 10, 5)
    raises("significant flag is only for lead rows", hits_frame, [{"vidx": 10, "kind": 0, "gene": 1, "pval": 0.1,
                                                                  "beta": 0.1, "significant": True}], 10, 5)
    raises("kind", hits_frame, bad(kind=6), 10, 5)
    raises("gene", hits_frame, bad(gene=65536), 10, 5)


@case
def hits_frame_reader_rejects():
    b = hits_frame(HITS_ROWS, 100, 4)
    raw = bytearray(pf.zstd_unframe(b, None, "t"))
    off = pf.hits_offsets(4, len(HITS_ROWS))

    def mutate(i, val, n=1):
        x = bytearray(raw)
        x[i:i + n] = int(val).to_bytes(n, "little")
        return reframe(x)

    raises("magic", pf.decode_hits_frame, mutate(0, 0x30485651 ^ 0xFF, 4), 100, "t")
    raises("reserved header field", pf.decode_hits_frame, mutate(6, 1, 2), 100, "t")
    raises("n_variants", pf.decode_hits_frame, mutate(4, 0, 2), 100, "t")
    raises("must be finite and at least 0", pf.decode_hits_frame, reframe(raw[:8] + struct.pack("<d", float("nan")) + raw[16:]), 100, "t")
    raises("must be finite and at least 0", pf.decode_hits_frame, reframe(raw[:8] + struct.pack("<d", -1.0) + raw[16:]), 100, "t")
    raises("decompressed length", pf.decode_hits_frame, reframe(raw + b"\0\0\0\0"), 100, "t")
    raises("decompressed length", pf.decode_hits_frame, mutate(pf.HITS_HEADER_LEN, 9, 2), 100, "t")   # a count that no longer sums
    raises("shorter than the 48-byte header", pf.decode_hits_frame, reframe(raw[:40]), 100, "t")
    raises("a kind above 5", pf.decode_hits_frame, mutate(off["kind"], 6), 100, "t")
    raises("reserved flag bits", pf.decode_hits_frame, mutate(off["flags"], 0x04), 100, "t")
    raises("flags bit 0 on a row without an intron", pf.decode_hits_frame, mutate(off["flags"], 1), 100, "t")
    raises("flags bit 0 on a row without an intron", pf.decode_hits_frame, mutate(off["flags"] + 5, 2), 100, "t")
    raises("an nlp code above", pf.decode_hits_frame, mutate(off["v1"], 65534, 2), 100, "t")
    raises("v2 must be 0 on kinds 0-1", pf.decode_hits_frame, mutate(off["v2"], 1, 2), 100, "t")
    raises("v2 must be 0 on kinds 0-1", pf.decode_hits_frame, mutate(off["v3"] + 2 * 5, 1, 2), 100, "t")  # v3 on a kind 4 row
    raises("a signed code of -32768", pf.decode_hits_frame, mutate(off["v3"], 0x8000, 2), 100, "t")
    raises("an intron field on a kind 0, 2 or 4 row", pf.decode_hits_frame, mutate(off["intron_start"], 1, 4), 100, "t")
    raises("an intron field on a kind 0, 2 or 4 row", pf.decode_hits_frame, mutate(off["cluster"], 1, 4), 100, "t")
    # padding after the counts and after the flags must be zero. 4 variants and 10 rows are both
    # already aligned, so this frame has no padding at all; the short frame leaves 2 bytes of each.
    assert pf.pad4(2 * 4) == 0 and pf.pad4(2 * len(HITS_ROWS)) == 0
    s_raw = bytearray(pf.zstd_unframe(hits_frame(SHORT_ROW, 3072, 5), None, "t"))
    s_off = pf.hits_offsets(5, 1)
    assert pf.pad4(2 * 5) == 2 and pf.pad4(2 * 1) == 2
    for i in (pf.HITS_HEADER_LEN + 2 * 5, s_off["flags"] + 1):
        x = bytearray(s_raw)
        x[i] = 1
        raises("padding is not zero", pf.decode_hits_frame, reframe(x), 3072, "t")
    # the scale rule: no row holds the largest code any more
    raises("the scale rule", pf.decode_hits_frame, mutate(off["v1"], 1, 2), 100, "t")
    raises("the scale rule", pf.decode_hits_frame, reframe(raw[:24] + struct.pack("<d", 0.0) + raw[32:]), 100, "t")
    # order inside one variant: kind must not decrease and v1 must not increase inside one kind
    raises("must go by kind", pf.decode_hits_frame, mutate(off["kind"], 1), 100, "t")
    raises("must go by kind", pf.decode_hits_frame, mutate(off["v1"] + 2 * 8, 65535, 2), 100, "t")


@case
def rsid_index_block_math():
    B = 8
    n = 3 * B + 5                      # a short last block
    rs = (np.arange(n, dtype=np.int64) + 1) * 3
    co = np.array([(i % 23) + 1 for i in range(n)], dtype=np.int64)
    vi = np.arange(n, dtype=np.int64) * 7
    buf, first = pf.encode_rsid_index(rs, co, vi, block_records=B)
    h = pf.parse_file_header(buf)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (pf.KIND_RSID, "all", n, B)
    assert len(buf) == 32 + 8 * n and first.tolist() == rs[::B].tolist() and first.size == 4
    for b, want in ((0, B), (1, B), (2, B), (3, 5)):
        o, ln = pf.rsid_block_range(b, n, B)
        assert (o, ln) == (32 + b * B * 8, want * 8)
        d = pf.decode_rsid_block(buf[o:o + ln], "t")
        assert d["rs_number"].tolist() == rs[b * B:b * B + want].tolist()
        assert d["vidx"].tolist() == vi[b * B:b * B + want].tolist()
        assert d["chrom"] == [pf.VARIANT_CHROMS[c - 1] for c in co[b * B:b * B + want].tolist()]
    raises("outside a file", pf.rsid_block_range, 4, n, B)
    raises("outside a file", pf.rsid_block_range, -1, n, B)
    # rsid_block_of picks the last block whose first record is at or below the wanted number
    assert pf.rsid_block_of(first, int(rs[0]) - 1) is None
    for b in range(4):
        assert pf.rsid_block_of(first, int(first[b])) == b
        assert pf.rsid_block_of(first, int(first[b]) + 1) == b
    assert pf.rsid_block_of(first, int(rs[-1]) + 10_000) == 3
    # a full lookup through the block math, hit and miss
    for i in (0, 1, B - 1, B, n - 1):
        b = pf.rsid_block_of(first, int(rs[i]))
        o, ln = pf.rsid_block_range(b, n, B)
        assert pf.rsid_find(pf.decode_rsid_block(buf[o:o + ln], "t"), int(rs[i])) == (pf.VARIANT_CHROMS[co[i] - 1], int(vi[i]))
    b = pf.rsid_block_of(first, int(rs[0]) + 1)
    o, ln = pf.rsid_block_range(b, n, B)
    assert pf.rsid_find(pf.decode_rsid_block(buf[o:o + ln], "t"), int(rs[0]) + 1) is None, "a miss inside the right block"
    # rejects
    raises("must strictly increase", pf.encode_rsid_index, np.array([5, 5]), np.array([1, 1]), np.array([0, 1]))
    raises("must strictly increase", pf.encode_rsid_index, np.array([5, 4]), np.array([1, 1]), np.array([0, 1]))
    raises("chr_ordinal", pf.encode_rsid_index, np.array([1, 2]), np.array([0, 1]), np.array([0, 1]))
    raises("chr_ordinal", pf.encode_rsid_index, np.array([1, 2]), np.array([1, 24]), np.array([0, 1]))
    raises("vidx", pf.encode_rsid_index, np.array([1, 2]), np.array([1, 1]), np.array([0, 1 << 27]))
    raises("rs_number", pf.encode_rsid_index, np.array([0, 2]), np.array([1, 1]), np.array([0, 1]))
    raises("records per block", pf.encode_rsid_index, rs, co, vi, 0)
    raises("whole number of 8-byte records", pf.decode_rsid_block, buf[32:32 + 12], "t")
    raises("whole number of 8-byte records", pf.decode_rsid_block, b"", "t")
    raises("chromosome ordinal outside", pf.decode_rsid_block, struct.pack("<II", 5, 0), "t")
    raises("must strictly increase", pf.decode_rsid_block, struct.pack("<IIII", 9, 1 << 27, 9, 1 << 27), "t")


def synth_variant_index(page_size: int = 512, frame_variants: int = 1024, block_records: int = 4096, seed: int = 9):
    rng = np.random.default_rng(seed)
    chroms = []
    for i in range(len(pf.VARIANT_CHROMS)):
        n_cis, n_tr = 1500 + 100 * i, 600 + 7 * i      # two trans-only pages, so n_cis + P is a real index
        pc, ptr = -(-n_cis // page_size), -(-n_tr // page_size)
        nf = -(-(n_cis + n_tr) // frame_variants)
        chroms.append({
            "n_cis": n_cis, "n_trans_only": n_tr,
            "page_off": 32 + np.cumsum(np.r_[0, rng.integers(900, 1200, pc + ptr)]),
            "page_first_position": np.sort(rng.choice(np.arange(1, 10 ** 8), pc + ptr, replace=False)),
            "hits_off": 32 + np.cumsum(np.r_[0, rng.integers(500, 800, nf)])})
    rsid_n = 3 * block_records + 17
    rsid_first = (np.arange(4, dtype=np.int64) + 1) * 1000
    return chroms, rsid_first, rsid_n


@case
def variant_index_round_trip():
    P, F, B = 512, 1024, 4096
    chroms, rsid_first, rsid_n = synth_variant_index(P, F, B)
    buf = pf.encode_variant_index(chroms, P, F, rsid_first, rsid_n, B, level=3)
    h = pf.parse_file_header(buf)
    assert (h["kind"], h["chrom"], h["count"], h["page_size"]) == (pf.KIND_VARIANT_INDEX, "all", 23, P)
    d = pf.decode_variant_index(buf, "t")
    assert (d["page_size"], d["frame_variants"], d["rsid_block_records"]) == (P, F, B)
    assert d["rsid_n_records"] == rsid_n and d["rsid_n_blocks"] == 4 and d["rsid_first"].tolist() == rsid_first.tolist()
    assert list(d["chroms"]) == list(pf.VARIANT_CHROMS)
    for name, src in zip(pf.VARIANT_CHROMS, chroms):
        c = d["chroms"][name]
        assert (c["n_cis"], c["n_trans_only"]) == (src["n_cis"], src["n_trans_only"])
        assert c["n_pages_cis"] == -(-src["n_cis"] // P) and c["n_pages_trans"] == -(-src["n_trans_only"] // P)
        assert c["n_frames"] == -(-(src["n_cis"] + src["n_trans_only"]) // F)
        assert c["page_off"].tolist() == list(src["page_off"]) and c["hits_off"].tolist() == list(src["hits_off"])
        assert c["page_first_position"].tolist() == list(src["page_first_position"])
    # the three reader algorithms of SPEC section 15
    name = "chr1"
    c, src = d["chroms"][name], chroms[0]
    n_cis = src["n_cis"]
    for vidx, want in ((0, 0), (P - 1, 0), (P, 1), (n_cis - 1, c["n_pages_cis"] - 1), (n_cis, c["n_pages_cis"]),
                       (n_cis + P, c["n_pages_cis"] + 1)):
        page, off, ln = pf.variant_index_page(d, name, vidx)
        assert page == want and off == c["page_off"][page] and ln == c["page_off"][page + 1] - c["page_off"][page]
    raises("outside", pf.variant_index_page, d, name, n_cis + src["n_trans_only"])
    raises("outside", pf.variant_index_page, d, name, -1)
    for vidx in (0, F - 1, F, n_cis + src["n_trans_only"] - 1):
        fr, off, ln = pf.variant_index_hits(d, name, vidx)
        assert fr == vidx // F and off == c["hits_off"][fr] and ln == c["hits_off"][fr + 1] - c["hits_off"][fr]
    raises("outside", pf.variant_index_hits, d, name, n_cis + src["n_trans_only"])
    first = c["page_first_position"]
    for k in (0, 1, c["n_pages_cis"] - 1):
        assert pf.variant_index_position(d, name, int(first[k]), "cis")[0] == k
        assert pf.variant_index_position(d, name, int(first[k]) + 1, "cis")[0] == k
    assert pf.variant_index_position(d, name, int(first[c["n_pages_cis"]]), "trans")[0] == c["n_pages_cis"]
    raises("below the", pf.variant_index_position, d, name, int(first[0]) - 1, "cis")
    # encoder rejects
    raises("chromosomes, expected", pf.encode_variant_index, chroms[:5], P, F, rsid_first, rsid_n, B, 3)
    raises("block samples", pf.encode_variant_index, chroms, P, F, rsid_first[:3], rsid_n, B, 3)
    raises("block samples", pf.encode_variant_index, chroms, P, F, np.array([5, 4, 3, 2]), rsid_n, B, 3)
    bad = [dict(c) for c in chroms]
    bad[0] = {**bad[0], "hits_off": bad[0]["hits_off"][:-1]}
    raises("do not follow", pf.encode_variant_index, bad, P, F, rsid_first, rsid_n, B, 3)
    bad2 = [dict(c) for c in chroms]
    bad2[0] = {**bad2[0], "page_off": np.r_[bad2[0]["page_off"][:2][::-1], bad2[0]["page_off"][2:]]}
    raises("must increase from byte 32", pf.encode_variant_index, bad2, P, F, rsid_first, rsid_n, B, 3)
    bad3 = [dict(c) for c in chroms]
    bad3[0] = {**bad3[0], "page_first_position": np.r_[0, bad3[0]["page_first_position"][1:]]}
    raises("page_first_position", pf.encode_variant_index, bad3, P, F, rsid_first, rsid_n, B, 3)
    raises("page size and frame variants", pf.encode_variant_index, chroms, 0, F, rsid_first, rsid_n, B, 3)
    # reader rejects, by mutating the payload inside the single zstd frame
    p = bytearray(pf.zstd_unframe(buf[32:], None, "t"))

    def vx(payload):
        return buf[:32] + reframe(payload)

    raises("magic", pf.decode_variant_index, vx(bytearray(b"XXXX") + p[4:]), "t")
    raises("a reserved header field", pf.decode_variant_index, vx(p[:6] + b"\x01\x00" + p[8:]), "t")
    raises("a reserved header field", pf.decode_variant_index, vx(p[:28] + b"\x01\x00\x00\x00" + p[32:]), "t")
    raises("chromosomes", pf.decode_variant_index, vx(p[:4] + struct.pack("<H", 22) + p[6:]), "t")
    raises("page size", pf.decode_variant_index, vx(p[:8] + struct.pack("<I", 256) + p[12:]), "t")
    raises("rsID blocks", pf.decode_variant_index, vx(p[:24] + struct.pack("<I", 9) + p[28:]), "t")
    raises("bytes after the rsID block samples", pf.decode_variant_index, vx(p + b"\0\0\0\0"), "t")
    raises("payload ends early", pf.decode_variant_index, vx(p[:-4]), "t")
    raises("header kind", pf.decode_variant_index, pf.file_header(pf.KIND_RSID, "all", 23, P) + buf[32:], "t")
    raises("payload is shorter", pf.decode_variant_index, vx(p[:16]), "t")


def real_gene(con, cfg, symbol: str, dof: int) -> dict:
    d = cfg.derived
    gid, chrom, gbin = con.execute(f"SELECT gene_id, chr, bin FROM '{cfg.tables / 'genes.parquet'}' WHERE symbol = ?", [symbol]).fetchone()
    vp = variants_sql(cfg, chrom)
    con.execute("DROP TABLE IF EXISTS v")
    con.execute(f"""
        CREATE TABLE v AS
        WITH vals AS (
            SELECT position, A1, A2, min(af) AS af, max(af) <> min(af) AS af_differs,
                   min(ma_samples) AS ma_samples, min(ma_count) AS ma_count
            FROM read_parquet('{cfg.tables}/cis_eqtl_nominal/chr={chrom}/*/*.parquet', hive_partitioning=false) GROUP BY 1, 2, 3)
        SELECT (row_number() OVER (ORDER BY vv.position, vv.A1, vv.A2) - 1) AS vidx,
               vv.position, vv.A1, vv.A2, vv.rs_number, vv.match, vals.af, vals.ma_samples, vals.ma_count
        FROM (SELECT position, A1, A2, rs_number, match FROM {vp} WHERE in_cis) vv
        LEFT JOIN vals USING (position, A1, A2)
        ORDER BY vidx
    """)
    vt = con.execute("SELECT * FROM v ORDER BY vidx").to_arrow_table()
    nom = con.execute(f"""
        SELECT v.vidx, n.position, n.A1, n.A2, n.rs_number, n.tss_distance, n.af, n.ma_samples, n.ma_count,
               n.pval_nominal, n.slope, n.slope_se, n.pip, n.cs_id
        FROM read_parquet('{cfg.tables}/cis_eqtl_nominal/chr={chrom}/bin={gbin}/data.parquet', hive_partitioning=false) n
        LEFT JOIN v USING (position, A1, A2)
        WHERE n.gene_id = ? ORDER BY v.vidx
    """, [gid]).to_arrow_table()
    gene = con.execute(f"SELECT * FROM '{cfg.tables / 'genes.parquet'}' WHERE gene_id = ?", [gid]).to_arrow_table().to_pylist()[0]
    exons = [(r[0], r[1]) for r in con.execute(f"SELECT start, \"end\" FROM '{cfg.tables / 'exons.parquet'}' WHERE gene_id = ? ORDER BY start, \"end\"", [gid]).fetchall()]
    splice = con.execute(f"""SELECT s.* EXCLUDE (gene_id, symbol, chr, tss), p.blk_off, p.blk_len FROM '{cfg.tables / 'splice_phenotypes.parquet'}' s
        LEFT JOIN '{d}/_tmp/pack_pointers/sqtl_{chrom}.parquet' p USING (phenotype_id)
        WHERE s.gene_id = ? ORDER BY s.cluster_id, s.intron_start, s.intron_end""", [gid]).to_arrow_table().to_pylist()

    vidx = nom["vidx"].to_numpy()
    n = nom.num_rows
    assert vidx.min() == vidx[0] and vidx.max() - vidx.min() + 1 == n and len(np.unique(vidx)) == n, "contiguous run"
    var_start = int(vidx[0])
    npos = nom["position"].to_numpy().astype(np.int64)
    anchors = npos - nom["tss_distance"].to_numpy()
    assert anchors.min() == anchors.max(), "constant anchor"

    P = cfg["packs"]["variant_page_size"]
    out = {"symbol": symbol, "chrom": chrom, "rows": n, "variants": vt.num_rows}
    for codec in ("raw", "zstd"):
        t0 = time.perf_counter()
        vbuf, offs = pf.encode_variants_file(chrom, vt["position"], vt["rs_number"], vt["af"], vt["ma_samples"], vt["ma_count"],
                                             vt["A1"], vt["A2"], vt["match"], P, codec, cfg["packs"]["zstd_level"])
        out[f"encode_variants_{codec}_s"] = time.perf_counter() - t0
        out[f"variants_{codec}_MB"] = len(vbuf) / 1e6
        t0 = time.perf_counter()
        whole = pf.decode_variants_file(vbuf)
        out[f"decode_variants_{codec}_s"] = time.perf_counter() - t0
        assert np.array_equal(whole["position"], vt["position"].to_numpy())
    var_off, var_len = pf.variant_range(offs, P, var_start, n)     # zstd file from the last loop

    cs = nom.filter(pa.compute.is_valid(nom["pip"]))
    rows_idx = (cs["vidx"].to_numpy() - var_start).astype(np.int64)
    details = pf.details_from_tables(gene, exons, splice)
    t0 = time.perf_counter()
    blk = pf.encode_gene_block(details, var_start, int(anchors[0]), nom["pval_nominal"], nom["slope"], nom["slope_se"],
                               rows_idx, cs["pip"], cs["cs_id"], cfg["packs"]["zstd_level"],
                               pos_first=int(npos[0]), pos_last=int(npos[-1]))
    out["encode_block_ms"] = (time.perf_counter() - t0) * 1e3
    t0 = time.perf_counter()
    b = pf.decode_gene_block(blk, dof, expect_blk_len=len(blk), expect_n_var=n, expect_var_start=var_start)
    pages = pf.decode_variant_pages(vbuf[var_off:var_off + var_len], var_start=var_start, n_var=n,
                                    pos_first=int(npos[0]), pos_last=int(npos[-1]))
    t = pf.gene_page_rows(b, pages)
    out["decode_block_pages_rows_ms"] = (time.perf_counter() - t0) * 1e3
    out["block_bytes"], out["variant_range_bytes"] = len(blk), var_len
    assert b["details"] == __import__("json").loads(pf.details_bytes(details))
    assert t.schema.equals(pf.READER_SCHEMA) and t.num_rows == n

    for col in ("position", "A1", "A2", "rs_number", "tss_distance", "ma_samples", "ma_count", "pip", "cs_id"):
        assert t[col].to_pylist() == nom[col].to_pylist(), col
    af_src = nom["af"].to_numpy(zero_copy_only=False).astype(np.float64)
    af_out = t["af"].to_numpy(zero_copy_only=False).astype(np.float64)
    out["af_max_err"] = float(np.nanmax(np.abs(af_out - af_src)))
    assert out["af_max_err"] <= 7.7e-6

    p_src = nom["pval_nominal"].to_numpy(zero_copy_only=False)
    s_src = nom["slope"].to_numpy(zero_copy_only=False).astype(np.float64)
    se_src = nom["slope_se"].to_numpy(zero_copy_only=False).astype(np.float64)
    fin = ~np.isnan(p_src) & (p_src > 0)
    nstep = b["nlp_max"] / pf.NLP_MAXQ
    nlp_err = np.abs(-np.log10(t["pval_nominal"].to_numpy(zero_copy_only=False)[fin]) - (-np.log10(p_src[fin])))
    out["nlp_max_err"], out["nlp_half_step"] = float(nlp_err.max()), nstep / 2
    assert out["nlp_max_err"] <= nstep / 2 + 1e-9
    ok = ~np.isnan(s_src) & ~np.isnan(se_src)
    assert np.array_equal(np.isnan(b["slope_se"]), ~ok), "SE is null exactly where the source slope or SE is null"
    se_b, sl_b = pf.error_bounds(p_src, s_src, se_src, b["nlp_max"], b["lse_min"], b["lse_max"], dof)
    dom = ~np.isnan(sl_b)
    se_err = np.abs(b["slope_se"] - se_src)
    sl_err = np.abs(b["slope"] - s_src)
    assert np.all(np.isfinite(b["slope"][dom])) and np.all(np.isfinite(b["slope_se"][dom]))
    assert np.all(se_err[dom] <= pf.BOUND_FACTOR * se_b[dom]) and np.all(sl_err[dom] <= pf.BOUND_FACTOR * sl_b[dom])
    s_out = t["slope"].to_numpy(zero_copy_only=False).astype(np.float64)
    se_out = t["slope_se"].to_numpy(zero_copy_only=False).astype(np.float64)
    f32 = 6e-8                                             # the float32 output column adds half an ulp
    assert np.all(np.abs(s_out[dom] - s_src[dom]) <= pf.BOUND_FACTOR * sl_b[dom] + f32 * np.abs(s_src[dom]) + 1e-12)
    assert np.all(np.abs(se_out[dom] - se_src[dom]) <= pf.BOUND_FACTOR * se_b[dom] + f32 * se_src[dom])
    tt = np.abs(s_src / se_src)
    out["bands"] = {}
    for lo, hi in T_BANDS:
        m = dom & (tt >= lo) & (tt < hi)
        if m.any():
            out["bands"][f"[{lo}, {hi})"] = (int(m.sum()), float(se_err[m].max()), float(sl_err[m].max()),
                                             float((se_err[m] / se_b[m]).max()), float((sl_err[m] / sl_b[m]).max()))
    out["slope_max_err"], out["se_max_err"] = float(sl_err[dom].max()), float(se_err[dom].max())
    out["slope_err_over_se"] = float((sl_err[dom] / se_src[dom]).max())
    out["rows_p_reads_1"] = int((b["nlp_code"] == 0).sum())
    out["null_rows"] = int((~dom).sum())
    return out


def real_data() -> None:
    import duckdb
    from .common import Config
    cfg = Config()
    con = duckdb.connect()
    assert con.execute("SELECT current_setting('default_collation')").fetchone()[0] == "", "byte-wise VARCHAR order"
    dof = cfg["packs"]["dof"]["eqtl"]
    for symbol in ("FLNC", "SYNPO2L"):
        r = real_gene(con, cfg, symbol, dof)
        print(f"  real {symbol} ({r['chrom']}): {r['rows']:,} rows, {r['variants']:,} cis variants on the chromosome")
        print(f"    variants file raw {r['variants_raw_MB']:.2f} MB (encode {r['encode_variants_raw_s']:.2f} s, decode {r['decode_variants_raw_s']:.2f} s); "
              f"zstd {r['variants_zstd_MB']:.2f} MB (encode {r['encode_variants_zstd_s']:.2f} s, decode {r['decode_variants_zstd_s']:.2f} s)")
        print(f"    block {r['block_bytes']:,} B, variant range {r['variant_range_bytes']:,} B (zstd); encode block {r['encode_block_ms']:.1f} ms; "
              f"decode block + pages + rows {r['decode_block_pages_rows_ms']:.1f} ms")
        print(f"    af max err {r['af_max_err']:.2e}; -log10 p max err {r['nlp_max_err']:.2e} (half step {r['nlp_half_step']:.2e}); "
              f"slope_se max err {r['se_max_err']:.2e}; derived slope max err {r['slope_max_err']:.2e} ({r['slope_err_over_se']:.2e} of SE)")
        for band, (cnt, se_e, sl_e, se_r, sl_r) in r["bands"].items():
            print(f"    |t| in {band}: {cnt:,} rows, max abs err SE {se_e:.2e}, slope {sl_e:.2e}; max err / bound SE {se_r:.3f}, slope {sl_r:.3f}")
        print(f"    rows whose p reads back as 1: {r['rows_p_reads_1']}; null rows: {r['null_rows']}")


def test_synthetic_cases():
    """pytest entry (`uv run pytest pipeline/test_packfmt.py`): every synthetic case. The real-data round
    trip runs from `uv run python -m pipeline.test_packfmt`."""
    for fn in CASES:
        fn()


def main() -> int:
    for fn in CASES:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"{len(CASES)} synthetic cases passed")
    if "--synthetic" not in sys.argv:
        real_data()
        print("real-data round trip passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
