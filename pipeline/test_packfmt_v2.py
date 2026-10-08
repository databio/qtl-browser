"""Tests for the qtlstore v2 codec: cis blocks, trans frames, GWAS blocks.

    uv run python -m pipeline.test_packfmt_v2

Synthetic throughout, so it needs no data. Two kinds of case:

- **round trip**: encode, decode, and check every stored value comes back within the limit the
  header itself declares -- which is also a check that `theoretical_limit` is computable from the
  header alone, since that is all the test gives it.
- **the v1 failures**: a `p` of 0 and a `-log10 p` over 300 both broke v1, the first by losing the
  effect size and the second by failing the build. They are the reason this format version exists,
  so each gets a case naming what it used to do.
"""
from __future__ import annotations

import math
import sys

import numpy as np
from scipy.special import stdtr

from . import packfmt_v1 as v1
from . import packfmt_v2 as v2

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


DET = {"v": 1, "phenotype_type": "ge", "phenotype_id": "G1"}


def rows(n=4000, dof=435, seed=0):
    """A phenotype's nominal rows where p is exactly the two-sided Student-t tail of beta/se, the
    relation a real source satisfies."""
    rng = np.random.default_rng(seed)
    se = np.exp(rng.normal(-2.0, 0.3, n))
    t = rng.standard_t(dof, n)
    beta = t * se
    p = 2.0 * stdtr(float(dof), -np.abs(t))
    return p, beta, se


def limits_hold(d: dict, p, beta, se) -> dict:
    """`measured_worst` against `theoretical_limit`, the latter rebuilt from the decoded header
    alone. Returns both, and raises if the measurement exceeds the limit -- which the arithmetic
    forbids, so it would mean the encoder or the scales are wrong."""
    tl = v2.theoretical_limit(d["nlp_max"], d["beta_max"], d["lse_min"], d["lse_max"])
    mw = v2.measured_worst(p, beta, se, d["nlp"], d["beta"], d["se"])
    over = [k for k in ("neglog10p", "beta_over_se", "se_rel") if mw[k] > tl[k] * (1 + 1e-9)]
    if over:
        raise AssertionError(f"measured_worst exceeds theoretical_limit for {over}: {mw} vs {tl}")
    return {"limit": tl, "worst": mw}


# ---- the codes ---------------------------------------------------------------------------------
@case
def beta_codes_round_trip_and_reserve_null():
    b = np.array([0.0, 1.5, -1.5, 0.3, np.nan, np.inf])
    q, bmax = v2.quantize_beta(b)
    assert bmax == 1.5, bmax
    assert q[4] == v2.BETA_NULL and q[5] == v2.BETA_NULL
    assert abs(int(q[1])) == v2.BETA_MAXQ and int(q[2]) == -v2.BETA_MAXQ
    out = v2.dequantize_beta(q, bmax)
    assert np.isnan(out[4]) and np.isnan(out[5])
    for i in (0, 1, 2, 3):
        assert abs(out[i] - b[i]) <= bmax / (2 * v2.BETA_MAXQ) * (1 + 1e-12), (i, out[i], b[i])
    assert out[0] == 0.0 and out[1] > 0 and out[2] < 0        # sign is intrinsic to the i16


@case
def se_codes_use_all_sixteen_bits():
    """v1 spent bit 15 of the SE field on the slope's sign, leaving 15 bits of magnitude. In v2 beta
    carries the sign, so se gets 16 bits and its relative error halves. Checked against v1 on the
    same input rather than asserted."""
    se = np.exp(np.linspace(-4, -1, 500))
    q2, lo2, hi2 = v2.quantize_se(se)
    d2 = v2.dequantize_se(q2, lo2, hi2)
    worst2 = float(np.max(np.abs(d2 - se) / se))
    q1, lo1, hi1 = v1.quantize_se(se, np.ones_like(se))
    d1, _ = v1.dequantize_se(q1, lo1, hi1)
    worst1 = float(np.max(np.abs(d1 - se) / se))
    assert worst2 < worst1, (worst2, worst1)
    assert 1.7 < worst1 / worst2 < 2.3, worst1 / worst2
    assert q2.max() <= v2.SE_MAXQ and lo2 == math.log(se.min())


@case
def se_null_and_a_constant_block():
    q, lo, hi = v2.quantize_se(np.array([np.nan, -1.0, 0.0]))
    assert (q == v2.SE_NULL).all() and (lo, hi) == (0.0, 0.0)
    q, lo, hi = v2.quantize_se(np.array([0.5, 0.5, 0.5]))
    assert lo == hi and (q == 0).all()
    assert np.allclose(v2.dequantize_se(q, lo, hi), 0.5)


@case
def nlp_has_no_ceiling_where_v1_had_300():
    """v1's `quantize_nlp` raises above NLP_LIMIT = 300, which fails the build for 4 of ARIC's 118
    chr22 blocks (worst 317.86). That limit guarded v1's inverse-t slope rebuild, not its storage.
    v2 inverts nothing, so it encodes, and the resolution degrades gracefully."""
    p = np.array([1e-400, 1e-320, 1e-10, 0.5])    # 1e-400 underflows a double to exactly 0.0
    assert p[0] == 0.0
    try:
        v1.quantize_nlp(np.array([1e-320, 0.5]))
    except ValueError as e:
        assert "exceeds the format limit" in str(e), e
    else:
        raise AssertionError("v1 is expected to refuse -log10 p over 300")
    q, nlp_max = v2.quantize_nlp(p)
    assert q[0] == v2.NLP_ZERO                     # the underflowed one keeps its reserved code
    assert nlp_max > 319, nlp_max
    nlp = v2.dequantize_nlp(q, nlp_max)
    assert math.isinf(nlp[0])
    assert abs(nlp[1] - -math.log10(1e-320)) <= nlp_max / (2 * v2.NLP_MAXQ)


# ---- cis block ---------------------------------------------------------------------------------
@case
def cis_block_round_trips_within_its_own_limits():
    p, beta, se = rows()
    blk = v2.encode_gene_block(DET, 100, p, beta, se, [0, 0, 7], [0.2, 0.3, 0.9], [1, 2, 1], 19,
                               pos_first=1000, pos_last=99000)
    assert len(blk) == v2.BLOCK_HEADER_LEN + 6 * len(p) + 12 * 3 + (len(blk) - v2.BLOCK_HEADER_LEN
                                                                    - 6 * len(p) - 36)
    d = v2.decode_gene_block(blk, expect_blk_len=len(blk), expect_n_var=len(p), expect_var_start=100)
    assert d["details"] == DET and d["n_rows"] == len(p) and d["var_start"] == 100 and d["flags"] == 0
    assert d["cs_row"].tolist() == [0, 0, 7] and d["cs_id"].tolist() == [1, 2, 1]
    limits_hold(d, p, beta, se)


@case
def cis_block_is_six_bytes_a_row():
    p, beta, se = rows(n=10_000)
    blk = v2.encode_gene_block(DET, 0, p, beta, se, None, None, None, 19, pos_first=1, pos_last=10_000)
    overhead = len(blk) - 6 * len(p)
    assert overhead < 400, overhead             # header + details frame + padding only
    d = v2.decode_gene_block(blk)
    assert d["nlp_code"].dtype == np.uint16 and d["beta_code"].dtype == np.int16
    assert d["se_code"].dtype == np.uint16


@case
def cis_block_keeps_beta_when_p_underflowed():
    """The v1 failure this format exists for. With `p = 0` v1 stores NLP_ZERO and the reader gets
    `slope = NaN`, because the slope is rebuilt through `t(p, dof)` and no t exists. v2 stores beta,
    so the row keeps its effect size and only the p stays honestly unknown."""
    p, beta, se = rows()
    p[5] = 0.0
    blk1 = v1.encode_gene_block(DET, 0, 0, p, beta, se, None, None, None, 19, pos_first=1, pos_last=9)
    d1 = v1.decode_gene_block(blk1, 435)
    assert d1["pval_nominal"][5] == 0.0 and math.isnan(d1["slope"][5])

    blk2 = v2.encode_gene_block(DET, 0, p, beta, se, None, None, None, 19, pos_first=1, pos_last=9)
    d2 = v2.decode_gene_block(blk2)
    assert d2["nlp_code"][5] == v2.NLP_ZERO              # p is still unknown, truthfully
    assert d2["pval_nominal"][5] == 0.0 and math.isinf(d2["nlp"][5])
    assert np.isfinite(d2["beta"][5]) and np.isfinite(d2["se"][5])
    assert abs(d2["beta"][5] - beta[5]) / se[5] < 1e-3, (d2["beta"][5], beta[5])


@case
def cis_block_with_no_rows():
    blk = v2.encode_gene_block(DET, None, [], [], [], None, None, None, 19)
    d = v2.decode_gene_block(blk, expect_n_var=None, expect_var_start=None)
    assert d["n_rows"] == 0 and d["var_start"] is None and d["details"] == DET
    assert d["nlp_max"] == 0.0 and d["beta_max"] == 0.0


@case
def cis_block_rejects_a_damaged_header():
    p, beta, se = rows(n=50)
    blk = bytearray(v2.encode_gene_block(DET, 0, p, beta, se, None, None, None, 19,
                                         pos_first=1, pos_last=50))
    for name, mutate in (("magic", lambda b: b.__setitem__(slice(0, 4), b"XXXX")),
                         ("blk_len", lambda b: b.__setitem__(slice(4, 8), (len(b) + 4).to_bytes(4, "little")))):
        bad = bytearray(blk)
        mutate(bad)
        try:
            v2.decode_gene_block(bytes(bad))
        except ValueError:
            pass
        else:
            raise AssertionError(f"a damaged {name} must raise")
    try:
        v2.decode_gene_block(bytes(blk), expect_n_var=len(p) + 1)
    except ValueError as e:
        assert "n_rows" in str(e), e
    else:
        raise AssertionError("an n_var mismatch with the search index must raise")


@case
def cis_limit_is_computable_from_the_header_alone():
    """`theoretical_limit` takes only the four scales, so a reader can recompute it from a stored
    block. `se_min` comes from `exp(lse_min)`, which must be the real smallest se."""
    p, beta, se = rows()
    blk = v2.encode_gene_block(DET, 0, p, beta, se, None, None, None, 19, pos_first=1, pos_last=9)
    d = v2.decode_gene_block(blk)
    assert abs(math.exp(d["lse_min"]) - float(se.min())) / float(se.min()) < 1e-12
    tl = v2.theoretical_limit(d["nlp_max"], d["beta_max"], d["lse_min"], d["lse_max"])
    assert abs(tl["beta_over_se"] - d["beta_max"] / (2 * v2.BETA_MAXQ) / math.exp(d["lse_min"])) < 1e-15
    assert tl["af"] == 0.5 / v1.AF_MAXQ


# ---- trans frame -------------------------------------------------------------------------------
def trans_rows(n=600, seed=1):
    rng = np.random.default_rng(seed)
    ordinal = np.sort(rng.integers(1, 4, n))
    pos = np.zeros(n, dtype=np.int64)
    for o in np.unique(ordinal):
        m = ordinal == o
        pos[m] = np.sort(rng.integers(1000, 10_000_000, m.sum()))
    se = np.exp(rng.normal(-2.0, 0.3, n))
    beta = rng.standard_t(400, n) * se
    p = 10.0 ** -rng.uniform(1, 30, n)
    af = rng.uniform(0.01, 0.5, n)
    af_code = np.rint(af * v1.AF_MAXQ).astype("<u2")
    rs = rng.integers(0, 10**8, n)
    ref = ["A"] * n
    alt = ["C"] * n
    ref[0], alt[0] = "AT", "A"          # an indel, so it lands in the heap
    return ordinal, pos, rs, af_code, ref, alt, p, beta, se


@case
def trans_frame_round_trips_and_stores_se():
    """v1's trans frame stored `(nlp, beta)` and had the reader rebuild `se = |beta| / t(nlp, dof)` --
    the mirror of the cis side's missing beta. v2 stores se, so neither side derives anything."""
    ordinal, pos, rs, af_code, ref, alt, p, beta, se = trans_rows()
    fr = v2.encode_trans_frame(ordinal, pos, rs, af_code, ref, alt, p, beta, se, 19)
    d = v2.decode_trans_frame(fr)
    assert d["n"] == len(p)
    assert d["ordinal"].tolist() == ordinal.tolist()
    assert d["pos"].tolist() == pos.tolist()                 # delta coding, resets per chromosome
    assert d["rs_number"].tolist() == rs.tolist()
    assert d["ref"][0] == "AT" and d["alt"][0] == "A"        # out of the heap
    assert d["ref"][1] == "A" and d["alt"][1] == "C"         # out of the SNP code table
    assert np.allclose(d["af"], af_code / v1.AF_MAXQ)
    limits_hold(d, p, beta, se)
    assert "se_code" in d and np.all(np.isfinite(d["se"]))


@case
def trans_frame_rejects_unsorted_rows_and_null_p():
    ordinal, pos, rs, af_code, ref, alt, p, beta, se = trans_rows(n=20)
    bad = pos.copy()
    bad[5], bad[6] = bad[6], bad[5]
    if ordinal[5] == ordinal[6]:
        try:
            v2.encode_trans_frame(ordinal, bad, rs, af_code, ref, alt, p, beta, se, 19)
        except ValueError as e:
            assert "sorted" in str(e), e
        else:
            raise AssertionError("unsorted positions must raise")
    p2 = p.copy(); p2[3] = np.nan
    try:
        v2.encode_trans_frame(ordinal, pos, rs, af_code, ref, alt, p2, beta, se, 19)
    except ValueError as e:
        assert "null p" in str(e), e
    else:
        raise AssertionError("a null p must raise in a trans frame")


# ---- GWAS block --------------------------------------------------------------------------------
def gwas_rows(n=2048, seed=2):
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.choice(np.arange(1, 50_000_000), n, replace=False))
    se = np.exp(rng.normal(-2.0, 0.8, n))        # a GWAS block spans more se than a cis block
    beta = rng.standard_normal(n) * se
    p = 10.0 ** -rng.uniform(0.1, 40, n)
    af = rng.uniform(0.001, 0.5, n)
    n_values = np.array([42637, 422920, 937963], dtype=np.int64)
    n_code = rng.integers(0, len(n_values), n).astype(np.uint8)
    rs = rng.integers(0, 10**9, n)
    ref = ["A"] * n
    alt = ["G"] * n
    ref[3], alt[3] = "A", "ATTT"
    return pos, p, beta, se, af, n_code, rs, ref, alt, n_values


@case
def gwas_block_round_trips():
    pos, p, beta, se, af, n_code, rs, ref, alt, n_values = gwas_rows()
    fr = v2.encode_gwas_block(pos, p, beta, se, af, n_code, rs, ref, alt, 19)
    d = v2.decode_gwas_block(fr, n_values, expect_rows=len(pos))
    assert d["position"].tolist() == pos.tolist()
    assert d["rs_number"].tolist() == rs.tolist()
    assert d["rows"] == len(pos)
    # `n` is the sample size and `rows` the row count, the same way round as v1's GWAS decode
    assert d["n"].tolist() == n_values[n_code.astype(np.int64)].tolist()
    assert d["ref"][3] == "A" and d["alt"][3] == "ATTT"
    assert np.max(np.abs(d["af"] - af)) <= 0.5 / v1.AF_MAXQ * (1 + 1e-9)
    limits_hold(d, p, beta, se)


@case
def gwas_block_is_smaller_than_v1s():
    """v1's GWAS row is 21 bytes and lossless; v2's is 18 and quantized. The payload shrinks because
    9 bytes of p/beta/se become 6, which more than pays for the 32 bytes of scales once per block."""
    B = 2048
    v2_payload = v2.GWAS_HEADER_LEN + v2.GWAS_ROW_BYTES * B
    v1_payload = v1.GWAS_HEADER_LEN + v1.GWAS_ROW_BYTES * B
    assert v1_payload == 43_016 and v2_payload == 36_904, (v1_payload, v2_payload)
    assert v2_payload < v1_payload
    assert v2.GWAS_HEADER_LEN / B < 0.02                  # the scales cost under 0.02 B/row


@case
def gwas_block_rejects_bad_input():
    pos, p, beta, se, af, n_code, rs, ref, alt, n_values = gwas_rows(n=100)
    for label, kw in (("decreasing position", {"position": pos[::-1]}),
                      ("af over 1", {"af": np.where(np.arange(100) == 4, 1.5, af)}),
                      ("non-finite beta", {"beta": np.where(np.arange(100) == 2, np.inf, beta)})):
        args = {"position": pos, "p": p, "beta": beta, "se": se, "af": af, "n_code": n_code,
                "rs_number": rs, "ref": ref, "alt": alt, **kw}
        try:
            v2.encode_gwas_block(args["position"], args["p"], args["beta"], args["se"], args["af"],
                                 args["n_code"], args["rs_number"], args["ref"], args["alt"], 19)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{label} must raise")
    fr = v2.encode_gwas_block(pos, p, beta, se, af, n_code, rs, ref, alt, 19)
    try:
        v2.decode_gwas_block(fr, n_values[:1])            # an n_code past the table
    except ValueError as e:
        assert "n_code" in str(e), e
    else:
        raise AssertionError("an n_code outside the index table must raise")


# ---- across all three --------------------------------------------------------------------------
@case
def every_kind_uses_the_same_triple_and_scales():
    """The point of v2: one encoding of (p, beta, se) in all three object kinds, each against four
    f64 scales of its own scope. Asserted on the headers, so a future divergence fails here."""
    assert (v2.ROW_BYTES, v2.TRANS_ROW_BYTES, v2.GWAS_ROW_BYTES) == (6, 18, 18)
    # identity bytes: cis addresses rows by a contiguous vidx run, so it carries none
    assert v2.TRANS_ROW_BYTES - v2.ROW_BYTES == 12 and v2.GWAS_ROW_BYTES - v2.ROW_BYTES == 12
    p, beta, se = rows(n=200)
    cis = v2.decode_gene_block(v2.encode_gene_block(DET, 0, p, beta, se, None, None, None, 19,
                                                    pos_first=1, pos_last=200))
    o, ps, rs, ac, rf, al, tp, tb, ts = trans_rows(n=200)
    tr = v2.decode_trans_frame(v2.encode_trans_frame(o, ps, rs, ac, rf, al, tp, tb, ts, 19))
    gp, gpv, gb, gs, gaf, gn, grs, gr, ga, nv = gwas_rows(n=200)
    gw = v2.decode_gwas_block(v2.encode_gwas_block(gp, gpv, gb, gs, gaf, gn, grs, gr, ga, 19), nv)
    for name, d in (("cis", cis), ("trans", tr), ("gwas", gw)):
        for k in ("nlp_max", "beta_max", "lse_min", "lse_max"):
            assert k in d, (name, k)
        for k in ("nlp_code", "beta_code", "se_code", "nlp", "beta", "se"):
            assert k in d, (name, k)
        tl = v2.theoretical_limit(d["nlp_max"], d["beta_max"], d["lse_min"], d["lse_max"])
        assert set(tl) == {"neglog10p", "beta_over_se", "se_rel", "af"}, (name, tl)


def main() -> int:
    bad = 0
    for fn in CASES:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001 - report every case
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(CASES) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
