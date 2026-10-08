"""qtlstore v2 codec: p, beta and se all stored, nothing derived.

`plans/2026-10-07-qtlstore-v2.md` is the design and the evidence. This module is the cis results
block; trans frames and GWAS blocks follow once this is verified against a real build.

Everything that does not change is imported from `packfmt_v1` rather than copied, so the diff
between the formats is exactly the delta: zstd framing, the details JSON rules, the nlp quantizer and
its reserved codes, the credible-set record, and the column coercion helpers all stand.

What changes, and why (every number measured, see the plan):

- **beta is stored**, as an `i16` against a per-block `beta_max`, instead of being rebuilt as
  `sign * se * t(nlp, dof)`. That rebuild is the source of every defect this format version exists
  to remove: a `p` of 0 lost the effect size, `-log10 p` over 300 failed the build outright, a null
  `dof` cost an experiment every beta, and the error grew as `t^3` so it breached its own 4.47e-03
  budget -- 4.466e-03 on TOPCHeF's leafcutter set, over 499,596,989 rows, in the store live today.
- **se takes the freed sign bit.** The sign now rides in beta's `i16`, so se uses all 16 bits
  instead of 15 and its error halves.
- **the row array is columnar**, `nlp[n]` then `beta[n]` then `se[n]`, not interleaved. A reader that
  wants only `-log10 p` -- the locus plot -- can sub-range the first `2n` bytes rather than all `6n`.
- **the block header grows 64 -> 72** for `beta_max`, and v1's dead `anchor` i32 (always 0) becomes
  `flags`, reserved 0, which is where a future per-block `i32` beta would be signalled (plan G1).

The rulers are not interchangeable and the reasons are measured over 170 ARIC blocks:

- `se` spans a median factor of 5.2 within one block, is strictly positive, and gets a **log** ruler:
  constant relative error, 1.3e-05 at 16 bits.
- `|beta|` spans a median factor of 28,264 and crosses zero, so no log ruler exists; a **linear**
  ruler bounds the error at `(beta_max / se_min) / 65,532` in units of a row's own SE. Note the
  denominator: `beta_max` and `se_min` are generally different rows, so the ratio exceeds any row's
  own `|z|` -- 141 against a `z_max` of 9.8 on a DCM GWAS block. The 4.47e-03 budget is reached at a
  ratio of 293.
"""
from __future__ import annotations

import json
import math
import struct

import numpy as np

from . import packfmt_v1 as pf
# unchanged from v1: framing, details rules, the nlp reserved codes and dequantizer, the
# credible-set record, column coercion
from .packfmt_v1 import (  # noqa: F401 - re-exported as the v2 surface
    CS_DTYPE, CS_RECORD_LEN, NLP_MAXQ, NLP_NULL, NLP_ZERO, U32_MAX, UNCHECKED,
    details_bytes, dequantize_nlp, p_from_nlp, zstd_frame, zstd_unframe, pad4,
)
# Also unchanged, and re-exported so a builder handed `codec=packfmt_v2` never reaches past it to
# `pf` for part of an object. Hit records hold `value`/`beta` as `f4` and the GWAS index holds
# offsets, so neither carries a quantized statistic and neither has a v2 layout to differ from.
from .packfmt_v1 import (  # noqa: F401
    AF_MAXQ, AF_NULL, HIT_CS, HIT_DTYPE, HIT_LEAD, HIT_TRANS, HITS_FRAME_VARIANTS, MATCH_CODES,
    SNP_CODES, decode_gwas_index_payload, decode_hits_frame, encode_gwas_index_payload,
    encode_hits_body, gwas_shared_codes, gwas_window, hit_records, hits_frame_table, t_from_p,
)

FORMAT_VERSION = 2

# ---- codes ------------------------------------------------------------------------------------
BETA_NULL = -32768                  # i16; the rest of the range is magnitude with sign
BETA_MAXQ = 32766
SE_NULL = 0xFFFF                    # u16; no sign bit any more, so 16 bits of magnitude
SE_MAXQ = 65534
# the 4.47e-03 slope budget of SPEC section 9 is reached at this beta_max/se_min ratio
BETA_BUDGET_RATIO = 293

MAGIC_BLOCK = b"QGB2"
BLOCK_HEADER_LEN = 72
_BLOCK_HEADER = struct.Struct("<4sIIIIIIIddddII")
ROW_BYTES = 6
assert _BLOCK_HEADER.size == BLOCK_HEADER_LEN
NO_VAR_START = 0xFFFFFFFF


def quantize_nlp(p) -> tuple[np.ndarray, float]:
    """`u16` codes and `nlp_max`, with **no ceiling on nlp_max**.

    Identical to v1's quantizer except that v1 raises above `NLP_LIMIT = 300`. That limit was never
    a storage limit -- `nlp_max` is an f64 and the code is relative to it -- it existed because v1's
    reader inverted `p` to rebuild the slope and `10**-nlp` underflows a double past ~308. v2
    inverts nothing. 4 of ARIC's 118 chr22 blocks exceed 300 (worst 317.86) and fail the v1 build
    outright; here they encode, and the resolution degrades gracefully: at `nlp_max = 10,000` the
    step is 0.15 in -log10 p, which is `dz ~ 0.0016` at `z ~ 214`.

    Reserved codes are v1's: 65535 for a null or NaN p (never tested), 65534 for p == 0 (underflowed
    in the source, so no finite -log10 p exists). Both still mean what they meant -- but a `p` of 0
    no longer costs the row its effect size, because beta is stored.
    """
    p = np.asarray(p, dtype=np.float64)
    null = np.isnan(p)
    zero = p == 0
    ok = ~null & ~zero
    if np.any((p[ok] < 0) | (p[ok] > 1)):
        raise ValueError("p outside [0, 1]")
    nlp = np.zeros(p.shape, dtype=np.float64)
    nlp[ok] = np.maximum(-np.log10(p[ok]), 0.0)          # -log10(1) is -0.0
    nlp_max = float(nlp[ok].max()) if np.any(ok) else 0.0
    q = np.zeros(p.shape, dtype=np.uint16)
    if nlp_max > 0:
        q[ok] = np.rint(nlp[ok] / nlp_max * NLP_MAXQ).astype(np.uint16)
    q[zero] = NLP_ZERO
    q[null] = NLP_NULL
    return q, nlp_max


def quantize_beta(beta) -> tuple[np.ndarray, float]:
    """`i16` codes and `beta_max` for one scope: `q = rint(beta / beta_max * 32766)`.

    Linear, because beta crosses zero and `log(0)` does not exist. Null or non-finite -> -32768.
    `beta_max` is the largest `|beta|` over the non-null rows, 0.0 when every row is null."""
    b = np.asarray(beta, dtype=np.float64)
    null = ~np.isfinite(b)
    beta_max = float(np.abs(b[~null]).max()) if np.any(~null) else 0.0
    q = np.zeros(b.shape, dtype=np.int16)
    if beta_max > 0:
        scaled = np.rint(b[~null] / beta_max * BETA_MAXQ)
        if np.any(np.abs(scaled) > BETA_MAXQ):           # rint can land one past on the extreme
            scaled = np.clip(scaled, -BETA_MAXQ, BETA_MAXQ)
        q[~null] = scaled.astype(np.int16)
    q[null] = BETA_NULL
    return q, beta_max


def dequantize_beta(q, beta_max: float) -> np.ndarray:
    """float64 beta = `q * beta_max / 32766`; NaN for the null code."""
    q = np.asarray(q, dtype=np.int16)
    out = q.astype(np.float64) * (float(beta_max) / BETA_MAXQ)
    out[q == BETA_NULL] = np.nan
    return out


def quantize_se(se) -> tuple[np.ndarray, float, float]:
    """`u16` log codes and `(lse_min, lse_max)`: `q = rint((ln se - lse_min) / span * 65534)`.

    A log ruler gives constant *relative* error, which suits a quantity whose within-block spread is
    a small factor and which is strictly positive. No sign bit: beta carries the sign in v2, so all
    16 bits are magnitude. Null, non-finite or non-positive se -> 0xFFFF."""
    s = np.asarray(se, dtype=np.float64)
    null = ~np.isfinite(s) | (s <= 0)
    q = np.full(s.shape, SE_NULL, dtype=np.uint16)
    if not np.any(~null):
        return q, 0.0, 0.0
    lse = np.log(s[~null])
    lse_min, lse_max = float(lse.min()), float(lse.max())
    span = lse_max - lse_min
    q[~null] = 0 if span == 0 else np.rint((lse - lse_min) / span * SE_MAXQ).astype(np.uint16)
    return q, lse_min, lse_max


def dequantize_se(q, lse_min: float, lse_max: float) -> np.ndarray:
    """float64 se = `exp(lse_min + q * span / 65534)`; NaN for 0xFFFF."""
    q = np.asarray(q, dtype=np.uint16)
    step = (float(lse_max) - float(lse_min)) / SE_MAXQ
    out = np.exp(float(lse_min) + q.astype(np.float64) * step)
    out[q == SE_NULL] = np.nan
    return out


# ---- precision ---------------------------------------------------------------------------------
def theoretical_limit(nlp_max: float, beta_max: float, lse_min: float, lse_max: float) -> dict:
    """The worst error the codes admit: half a quantization step in each quantity's own unit.

    Computed from the four scales and nothing else -- no row is decoded and no source value is
    consulted -- so **any reader can recompute it from a stored block header**, which is what makes
    it a property of the format rather than of a build. It holds for every row in the scope,
    including rows a later build with the same scales would write.

    Rounding to the nearest code leaves the true value at most half a step away, so each entry is
    `step / 2`:

    - `neglog10p`: the step is `nlp_max / 65533`. Absolute, in -log10 p.
    - `beta_over_se`: the step is `beta_max / 32766` and it is **absolute**, so the row it hurts most
      is the one with the smallest standard error -- generally **not** the row that set `beta_max`,
      which is why this is a ratio across different rows and runs well above any row's own `|z|`.
      `se_min` comes from `exp(lse_min)`, since `lse_min` is the log of the smallest stored se.
      Expressed per row's own SE so it is comparable to SPEC section 9's 4.47e-03 budget, to the
      trans side and to the GWAS side.
    - `se_rel`: the ruler is linear in `ln se`, so a half-step in log space is a constant
      *relative* error, `expm1(span / (2 * 65534))`.
    - `af`: the variant catalog's own `u16` code, `0.5 / 65534`.
    """
    se_min = math.exp(float(lse_min))
    return {
        "neglog10p": float(nlp_max) / (2 * NLP_MAXQ),
        "beta_over_se": (float(beta_max) / (2 * BETA_MAXQ) / se_min) if se_min > 0 else 0.0,
        "se_rel": math.expm1((float(lse_max) - float(lse_min)) / (2 * SE_MAXQ)),
        "af": 0.5 / pf.AF_MAXQ,
    }


def measured_worst(p, beta, se, nlp_d, beta_d, se_d) -> dict:
    """The worst error that actually occurred: every row decoded and compared against the source
    values the adapter supplied, with `rows_compared` as the n of the comparison.

    Named `measured_worst` and not `measured` because every entry is a maximum, and the pair should
    say so on both sides: `theoretical_limit` announces that it is an extreme, so this must too, or
    a reader takes it for a typical error.

    A limit is arithmetic; this is fact. Keeping both shows when a limit is loose -- on real ARIC
    blocks these run at 91% to 100% of the limit, so the limit is a real prediction rather than a
    vacuous ceiling. And `measured_worst` exceeding `theoretical_limit` is a **detectable defect**:
    the arithmetic says it cannot happen, so the encoder or the scales are wrong."""
    p = np.asarray(p, dtype=np.float64)
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    out = {"neglog10p": 0.0, "beta_over_se": 0.0, "se_rel": 0.0, "rows_compared": 0}
    fin = np.isfinite(p) & (p > 0) & np.isfinite(nlp_d)
    if np.any(fin):
        out["neglog10p"] = float(np.max(np.abs(nlp_d[fin] - -np.log10(p[fin]))))
    ok = np.isfinite(beta) & np.isfinite(beta_d) & np.isfinite(se) & (se > 0)
    if np.any(ok):
        out["beta_over_se"] = float(np.max(np.abs(beta_d[ok] - beta[ok]) / se[ok]))
        out["rows_compared"] = int(ok.sum())
    sok = np.isfinite(se) & (se > 0) & np.isfinite(se_d)
    if np.any(sok):
        out["se_rel"] = float(np.max(np.abs(se_d[sok] - se[sok]) / se[sok]))
    return out


# ---- cis block ---------------------------------------------------------------------------------
def encode_gene_block(details: dict, var_start: int | None, p, beta, se,
                      cs_row, cs_pip, cs_id, level: int, *, flags: int = 0,
                      pos_first: int | None = None, pos_last: int | None = None) -> bytes:
    """One phenotype's block: `72 + 6n + 12k + details_zlen`, padded to 4.

    `p`, `beta` and `se` are the phenotype's nominal rows in vidx order, all three stored. The row
    array is columnar: `nlp[n]`, `beta[n]`, `se[n]`. A phenotype with no rows passes empty arrays and
    None for `var_start`, `pos_first` and `pos_last`.
    """
    if not isinstance(details, dict):
        raise ValueError("block: a v2 block needs a details object")
    p = pf._float_column(p, "p")
    n = len(p)
    b = pf._float_column(beta, "beta", n)
    s = pf._float_column(se, "se", n)
    if n == 0:
        if any(v is not None for v in (var_start, pos_first, pos_last)):
            raise ValueError("block: a block with no rows takes var_start, pos_first, pos_last = None")
        h_var_start, h_first, h_last = NO_VAR_START, 0, 0
    else:
        if any(v is None for v in (var_start, pos_first, pos_last)):
            raise ValueError("block: var_start, pos_first, pos_last are required when there are rows")
        if var_start < 0 or var_start + n - 1 >= NO_VAR_START:
            raise ValueError(f"block: rows {var_start}..{var_start + n - 1} do not fit below 0xFFFFFFFF")
        if not 1 <= pos_first <= pos_last <= U32_MAX:
            raise ValueError(f"block: pos_first {pos_first}, pos_last {pos_last} out of order or range")
        h_var_start, h_first, h_last = var_start, pos_first, pos_last
    if not 0 <= flags <= U32_MAX:
        raise ValueError(f"block: flags {flags} does not fit u32")

    nlp_q, nlp_max = quantize_nlp(p)
    beta_q, beta_max = quantize_beta(b)
    se_q, lse_min, lse_max = quantize_se(s)

    k = 0 if cs_row is None else len(cs_row)
    if k and n == 0:
        raise ValueError("block: credible-set records on a block with no rows")
    cs = np.zeros(k, dtype=CS_DTYPE)
    if k:
        rows = pf._int_column(cs_row, "cs_row", 0, n - 1, -1, k)
        ids = pf._int_column(cs_id, "cs_id", 0, 127, -1, k)
        pip = pf._float_column(cs_pip, "cs_pip", k)
        if np.any(rows < 0) or np.any(ids < 0):
            raise ValueError("block: null credible-set row or cs_id")
        if any((int(rows[i]), int(ids[i])) >= (int(rows[i + 1]), int(ids[i + 1])) for i in range(k - 1)):
            raise ValueError("block: credible-set (row, cs_id) pairs must be strictly ascending")
        if not np.all((pip >= 0) & (pip <= 1)):
            raise ValueError("block: pip outside [0, 1] or null")
        cs["row"], cs["pip"], cs["cs_id"] = rows, pip.astype(np.float32), ids

    dj = details_bytes(details)
    dz = zstd_frame(dj, level)
    body = BLOCK_HEADER_LEN + ROW_BYTES * n + CS_RECORD_LEN * k + len(dz)
    blk_len = body + pad4(body)
    if blk_len > U32_MAX:
        raise ValueError("block: longer than 4 GiB")
    header = _BLOCK_HEADER.pack(MAGIC_BLOCK, blk_len, n, h_var_start, flags, h_first, h_last, k,
                                nlp_max, lse_min, lse_max, beta_max, len(dz), len(dj))
    return b"".join((header, nlp_q.astype("<u2").tobytes(), beta_q.astype("<i2").tobytes(),
                     se_q.astype("<u2").tobytes(), cs.tobytes(), dz, bytes(blk_len - body)))


def decode_gene_block(buf: bytes, *, expect_blk_len: int | None = None,
                      expect_n_var=UNCHECKED, expect_var_start=UNCHECKED,
                      details_version: int = 1) -> dict:
    """Header fields, details dict, and per-row arrays. **Nothing is derived**: `beta` and `se` are
    decoded, not rebuilt, so there is no `dof` argument and no row is lost to a `p` of 0."""
    buf = bytes(buf)
    total = len(buf)
    if total < BLOCK_HEADER_LEN:
        raise ValueError(f"block: {total} bytes is shorter than the {BLOCK_HEADER_LEN}-byte header")
    (magic, blk_len, n, var_start, flags, pos_first, pos_last, n_cs,
     nlp_max, lse_min, lse_max, beta_max, dzlen, dlen) = _BLOCK_HEADER.unpack_from(buf, 0)
    if magic != MAGIC_BLOCK:
        raise ValueError(f"block: magic {magic!r} is not {MAGIC_BLOCK!r}")
    if blk_len != total:
        raise ValueError(f"block: length field {blk_len} != range length {total}")
    if expect_blk_len is not None and blk_len != expect_blk_len:
        raise ValueError(f"block: length field {blk_len} != search_index blk_len {expect_blk_len}")
    body = BLOCK_HEADER_LEN + ROW_BYTES * n + CS_RECORD_LEN * n_cs + dzlen
    if body + pad4(body) != blk_len:
        raise ValueError(f"block: 72 + 6*n_rows + 12*n_cs + details_zlen padded to 4 = "
                         f"{body + pad4(body)} != block length {blk_len}")
    if any(buf[body:]):
        raise ValueError("block: padding is not zero")
    if dzlen == 0:
        raise ValueError("block: a v2 block needs a details frame")
    for name, v in (("nlp_max", nlp_max), ("lse_min", lse_min), ("lse_max", lse_max), ("beta_max", beta_max)):
        if not math.isfinite(v):
            raise ValueError(f"block: {name} is not finite")
    if nlp_max < 0 or beta_max < 0 or lse_min > lse_max:
        raise ValueError("block: bad scales")
    if n == 0 and var_start != NO_VAR_START:
        raise ValueError("block: no rows but var_start is set")
    if expect_n_var is not UNCHECKED:
        want = 0 if expect_n_var is None else int(expect_n_var)
        if n != want:
            raise ValueError(f"block: n_rows {n} != search_index n_var {want}")
    if expect_var_start is not UNCHECKED:
        want = NO_VAR_START if expect_var_start is None else int(expect_var_start)
        if var_start != want:
            raise ValueError(f"block: var_start {var_start} != search_index var_start {want}")

    o = BLOCK_HEADER_LEN
    nlp_code = np.frombuffer(buf, "<u2", n, o)
    beta_code = np.frombuffer(buf, "<i2", n, o + 2 * n)
    se_code = np.frombuffer(buf, "<u2", n, o + 4 * n)
    o += ROW_BYTES * n
    cs = np.frombuffer(buf, CS_DTYPE, n_cs, o)
    o += CS_RECORD_LEN * n_cs
    dj = zstd_unframe(buf[o:o + dzlen], dlen, "block details")
    if dj.startswith(b"\xef\xbb\xbf"):
        raise ValueError("block details: JSON starts with a BOM")
    try:
        details = json.loads(dj.decode("utf-8"), parse_constant=pf._reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        raise ValueError(f"block details: {e}") from None
    if not isinstance(details, dict) or details.get("v") != details_version:
        raise ValueError(f"block details: not an object with v = {details_version}")
    nlp = dequantize_nlp(nlp_code, nlp_max)
    return {
        "blk_len": blk_len, "n_rows": n, "var_start": None if n == 0 else var_start, "flags": flags,
        "pos_first": pos_first, "pos_last": pos_last, "n_cs": n_cs,
        "nlp_max": nlp_max, "lse_min": lse_min, "lse_max": lse_max, "beta_max": beta_max,
        "details": details,
        "nlp_code": nlp_code, "beta_code": beta_code, "se_code": se_code,
        "nlp": nlp, "pval_nominal": p_from_nlp(nlp),
        "beta": dequantize_beta(beta_code, beta_max),
        "se": dequantize_se(se_code, lse_min, lse_max),
        "cs_row": cs["row"], "cs_pip": cs["pip"], "cs_id": cs["cs_id"],
    }


# ---- trans frame -------------------------------------------------------------------------------
MAGIC_TRANS = b"QTT3"      # v1 is QTT2; the version digit is the only difference, as QGB0 -> QGB2
TRANS_HEADER_LEN = 48
_TRANS_HEADER = struct.Struct("<4sIddddII")
TRANS_ROW_BYTES = 18
assert _TRANS_HEADER.size == TRANS_HEADER_LEN
# columns in payload order; 4+4+2+2+2+2+1+1 = 18
TRANS_COLUMNS = (("delta", "<u4", 4), ("rs_number", "<u4", 4), ("af_code", "<u2", 2),
                 ("nlp_code", "<u2", 2), ("beta_code", "<i2", 2), ("se_code", "<u2", 2),
                 ("ordinal", "u1", 1), ("allele", "u1", 1))
assert sum(w for _, _, w in TRANS_COLUMNS) == TRANS_ROW_BYTES


def encode_trans_frame(ordinal, pos, rs_number, af_code, ref, alt, p, beta, se, level: int) -> bytes:
    """One phenotype's trans rows as one zstd frame: `48 + 18n + heap_len`.

    The v1 frame stored `(nlp, beta)` and had the reader rebuild `se = |beta| / t(nlp, dof)` -- the
    mirror image of the cis side's missing beta, and with the same consequences: a null `dof` cost
    trans its standard errors and its r2, and nothing measured the error of the one value it
    reconstructed. v2 stores `se` too, so the trans and cis rows now hold the same triple and neither
    derives anything.

    A trans row carries its own variant identity because it can be anywhere in the genome: 12 of its
    18 bytes. Position is delta coded, absolute on the first row of each chromosome and from the row
    before otherwise, which is why rows must be sorted by `(ordinal, pos, ref, alt)`.
    """
    ordinal = np.asarray(ordinal, dtype=np.int64)
    n = len(ordinal)
    if n < 1:
        raise ValueError("trans frame: no rows")
    pos = np.asarray(pos, dtype=np.int64)
    if ordinal.min() < 1 or ordinal.max() > 255:
        raise ValueError("trans frame: chromosome ordinal outside 1..255")
    new_chrom = np.r_[True, ordinal[1:] != ordinal[:-1]]
    step = np.r_[0, np.diff(pos)]
    if np.any(np.diff(ordinal) < 0) or np.any(~new_chrom & (step < 0)):
        raise ValueError("trans frame: rows not sorted by (ordinal, pos)")
    if pos.min() < 1 or pos.max() > U32_MAX:
        raise ValueError("trans frame: position out of range")
    delta = np.where(new_chrom, pos, step).astype("<u4")

    b = pf._float_column(beta, "beta", n)
    if not np.all(np.isfinite(b)):
        raise ValueError("trans frame: beta must be finite")
    s = pf._float_column(se, "se", n)
    nlp_q, nlp_max = quantize_nlp(pf._float_column(p, "p", n))
    if np.any(nlp_q == NLP_NULL):
        raise ValueError("trans frame: null p")
    beta_q, beta_max = quantize_beta(b)
    se_q, lse_min, lse_max = quantize_se(s)

    codes = np.fromiter((pf.SNP_CODES.get((r, t), 0) for r, t in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(codes == 0):
        r, t = ref[i], alt[i]
        if not (r and t and r.isascii() and t.isascii()) or "\t" in r + t or "\n" in r + t:
            raise ValueError(f"trans frame: allele pair {r!r}/{t!r} is not two non-empty ASCII strings")
        heap.append(f"{r}\t{t}\n")
    heap_b = "".join(heap).encode("ascii")
    cols = {"delta": delta, "rs_number": np.where(np.asarray(rs_number, dtype=np.int64) < 0, 0,
                                                  np.asarray(rs_number, dtype=np.int64)),
            "af_code": np.asarray(af_code), "nlp_code": nlp_q, "beta_code": beta_q, "se_code": se_q,
            "ordinal": ordinal, "allele": codes}
    payload = b"".join(
        [_TRANS_HEADER.pack(MAGIC_TRANS, n, nlp_max, beta_max, lse_min, lse_max, len(heap_b), 0)]
        + [np.asarray(cols[c]).astype(t).tobytes() for c, t, _ in TRANS_COLUMNS] + [heap_b])
    assert len(payload) == TRANS_HEADER_LEN + TRANS_ROW_BYTES * n + len(heap_b)
    return zstd_frame(payload, level)


def decode_trans_frame(frame: bytes, what: str = "trans frame") -> dict:
    """One trans frame's rows. `se` is decoded now, not rebuilt, so there is no `dof` argument."""
    p = zstd_unframe(frame, None, what)
    if len(p) < TRANS_HEADER_LEN:
        raise ValueError(f"{what}: shorter than its header")
    magic, n, nlp_max, beta_max, lse_min, lse_max, heap_len, reserved = _TRANS_HEADER.unpack_from(p, 0)
    if magic != MAGIC_TRANS or reserved or n < 1:
        raise ValueError(f"{what}: bad magic, reserved field or row count")
    if len(p) != TRANS_HEADER_LEN + TRANS_ROW_BYTES * n + heap_len:
        raise ValueError(f"{what}: payload length {len(p)} != {TRANS_HEADER_LEN} + "
                         f"{TRANS_ROW_BYTES}n + heap_len")
    for name, v in (("nlp_max", nlp_max), ("beta_max", beta_max), ("lse_min", lse_min), ("lse_max", lse_max)):
        if not math.isfinite(v):
            raise ValueError(f"{what}: {name} is not finite")
    if nlp_max < 0 or beta_max < 0 or lse_min > lse_max:
        raise ValueError(f"{what}: bad scales")
    o, col = TRANS_HEADER_LEN, {}
    for name, t, w in TRANS_COLUMNS:
        col[name] = np.frombuffer(p, t, n, o)
        o += w * n
    ordinal = col["ordinal"].astype(np.int64)
    if ordinal.min() < 1 or np.any(np.diff(ordinal) < 0):
        raise ValueError(f"{what}: chromosome ordinals not >= 1 and non-decreasing")
    if np.any(col["allele"] > 12) or np.any(col["beta_code"] == BETA_NULL) or np.any(col["nlp_code"] == NLP_NULL):
        raise ValueError(f"{what}: reserved allele, beta or nlp code")
    new_chrom = np.r_[True, ordinal[1:] != ordinal[:-1]]
    d = col["delta"].astype(np.int64)
    run = np.cumsum(new_chrom) - 1
    starts = np.flatnonzero(new_chrom)
    csum = np.cumsum(np.where(new_chrom, 0, d))
    pos = d[starts][run] + csum - csum[starts][run]
    pairs = [pf.SNP_ALLELES.get(c, (None, None)) for c in col["allele"].tolist()]
    zero = np.flatnonzero(col["allele"] == 0)
    recs = p[o:].decode("ascii").split("\n")[:-1] if heap_len else []
    if len(recs) != zero.size:
        raise ValueError(f"{what}: heap holds {len(recs)} records for {zero.size} code-0 rows")
    for i, rec in zip(zero.tolist(), recs):
        pairs[i] = tuple(rec.split("\t"))
    nlp = dequantize_nlp(col["nlp_code"], nlp_max)
    af = np.where(col["af_code"] == pf.AF_NULL, np.nan, col["af_code"] / pf.AF_MAXQ)
    return {"n": n, "nlp_max": nlp_max, "beta_max": beta_max, "lse_min": lse_min, "lse_max": lse_max,
            "ordinal": ordinal, "pos": pos, "rs_number": col["rs_number"].astype(np.int64),
            "af_code": col["af_code"], "af": af, "nlp_code": col["nlp_code"], "nlp": nlp,
            "p": p_from_nlp(nlp), "beta_code": col["beta_code"],
            "beta": dequantize_beta(col["beta_code"], beta_max),
            "se_code": col["se_code"], "se": dequantize_se(col["se_code"], lse_min, lse_max),
            "ref": [x[0] for x in pairs], "alt": [x[1] for x in pairs]}


# ---- GWAS block --------------------------------------------------------------------------------
GWAS_HEADER_LEN = 40
_GWAS_HEADER = struct.Struct("<IIdddd")
GWAS_ROW_BYTES = 18
assert _GWAS_HEADER.size == GWAS_HEADER_LEN
GWAS_COLUMNS = (("position_delta", "<u4", 4), ("rs_number", "<u4", 4), ("af_code", "<u2", 2),
                ("nlp_code", "<u2", 2), ("beta_code", "<i2", 2), ("se_code", "<u2", 2),
                ("n_code", "u1", 1), ("allele", "u1", 1))
assert sum(w for _, _, w in GWAS_COLUMNS) == GWAS_ROW_BYTES
GWAS_MAX_N_VALUES = pf.GWAS_MAX_N_VALUES          # u8 n_code; gap G4 -- a per-variant-N meta-analysis breaks it


def gwas_codes(position, beta, se, af, p, rs_number, n, n_values) -> dict[str, np.ndarray]:
    """`gwas.build`'s precode step, same signature as v1's. v1 returns finished codes here, because
    its scales are fixed constants; v2's four scales are per block, so p/beta/se/af pass through
    **unquantized** and `encode_gwas_columns` does the quantizing once it knows the block's rows.

    So the only work left is the three columns neither codec quantizes -- position, rs_number and
    `n_code` -- which is exactly `pf.gwas_shared_codes`. Validating them here rather than per block
    means a chromosome with a bad N or a position out of order still fails before any bytes are
    written, as it did in v1."""
    shared = pf.gwas_shared_codes(position, rs_number, n, n_values)
    k = shared["position"].size
    for name, v in (("beta", beta), ("se", se), ("af", af), ("p", p)):
        if len(v) != k:
            raise ValueError(f"gwas {name}: {len(v)} rows, position has {k}")
    return {**shared, "beta": np.asarray(beta, dtype=np.float64), "se": np.asarray(se, dtype=np.float64),
            "af": np.asarray(af, dtype=np.float64), "p": np.asarray(p, dtype=np.float64)}


def encode_gwas_columns(codes: dict, ref: list, alt: list, level: int) -> bytes:
    """`gwas_codes` output sliced to one block, encoded. The counterpart of v1's name of the same
    shape, so `gwas.build` runs one loop over either codec."""
    return encode_gwas_block(codes["position"], codes["p"], codes["beta"], codes["se"], codes["af"],
                             codes["n_code"], codes["rs_number"], ref, alt, level)


def encode_gwas_block(position, p, beta, se, af, n_code, rs_number, ref, alt, level: int) -> bytes:
    """One GWAS block of B consecutive rows as one zstd frame: `40 + 18n + heap_len`.

    v1 stored this losslessly -- `i32` beta at `rint(b * 1e4)`, `u16` se likewise, and p as a `u16`
    mantissa with an `i8` exponent -- and `gwas_codes` refused any source carrying more precision
    than that holds. v2 quantizes it like everything else (plan D2), which **loses** that exactness
    and the build-time guard with it; the `precision` block is the replacement. The README's
    principle already said the store is a view and the source releases are the archive, and GWAS was
    its exception.

    The uncompressed row gets smaller -- 9 bytes of lossless p/beta/se become 6 codes, so at
    B = 2048 the payload goes 43,016 -> 36,904 bytes, -14.2%, and the four f64 scales cost
    0.0156 B/row -- but **the stored object does not**. Measured on the DCM GWAS, chr21/22, the same
    349,950 rows: v1 13.290 B/row against v2 13.585, +2.2%. Quantized codes spread their entropy
    across all 16 bits, where v1's `rint(x * 1e4)` leaves structure for zstd to find, and that costs
    more than the three bytes saved. The reason to do it anyway is uniformity, not size.

    `af` is the raw frequency; it is coded `rint(af * 65534)` here rather than v1's `rint(af * 1e4)`,
    which aligns the GWAS side with the variant catalog and improves its af limit from 5.0e-05 to
    7.63e-06 (plan D8).
    """
    pos = np.asarray(position, dtype=np.int64)
    n = pos.size
    if n < 1:
        raise ValueError("gwas block: needs at least one row")
    if np.any(np.diff(pos) < 0):
        raise ValueError("gwas block: decreasing position")
    if pos.min() < 1 or pos.max() > U32_MAX:
        raise ValueError("gwas block: position out of range")
    for name, v in (("p", p), ("beta", beta), ("se", se), ("af", af), ("n_code", n_code),
                    ("rs_number", rs_number), ("ref", ref), ("alt", alt)):
        if len(v) != n:
            raise ValueError(f"gwas block: {name} has {len(v)} rows, position has {n}")
    deltas = np.empty(n, dtype="<u4")
    deltas[0] = pos[0]
    deltas[1:] = np.diff(pos)

    b = pf._float_column(beta, "beta", n)
    if not np.all(np.isfinite(b)):
        raise ValueError("gwas block: beta must be finite")
    nlp_q, nlp_max = quantize_nlp(pf._float_column(p, "p", n))
    if np.any(nlp_q == NLP_NULL):
        raise ValueError("gwas block: null p")
    beta_q, beta_max = quantize_beta(b)
    se_q, lse_min, lse_max = quantize_se(pf._float_column(se, "se", n))
    a = pf._float_column(af, "af", n)
    if np.any(np.isfinite(a) & ((a < 0) | (a > 1))):
        raise ValueError("gwas block: af outside [0, 1]")
    af_q = np.where(np.isfinite(a), np.rint(np.nan_to_num(a) * pf.AF_MAXQ), pf.AF_NULL).astype("<u2")

    allele = np.fromiter((pf.SNP_CODES.get((x, y), 0) for x, y in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(allele == 0).tolist():
        for s, nm in ((ref[i], "ref"), (alt[i], "alt")):
            if not isinstance(s, str) or not s or not s.isascii() or "\t" in s or "\n" in s:
                raise ValueError(f"gwas {nm}: row {i} allele {s!r} is not an ASCII string without tab or newline")
        heap.append(f"{ref[i]}\t{alt[i]}\n")
    heap_b = "".join(heap).encode("ascii")
    rs = np.asarray(rs_number, dtype=np.int64)
    cols = {"position_delta": deltas, "rs_number": np.where(rs < 0, 0, rs), "af_code": af_q,
            "nlp_code": nlp_q, "beta_code": beta_q, "se_code": se_q,
            "n_code": np.asarray(n_code), "allele": allele}
    payload = b"".join([_GWAS_HEADER.pack(n, len(heap_b), nlp_max, beta_max, lse_min, lse_max)]
                       + [np.asarray(cols[c]).astype(t).tobytes() for c, t, _ in GWAS_COLUMNS] + [heap_b])
    assert len(payload) == GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + len(heap_b)
    return zstd_frame(payload, level)


def decode_gwas_block(frame: bytes, n_values, *, expect_rows: int | None = None,
                      what: str = "gwas block") -> dict:
    """One GWAS block's rows. `n_values` is the index's table of distinct n, for `n_code`."""
    buf = zstd_unframe(frame, None, what)
    if len(buf) < GWAS_HEADER_LEN:
        raise ValueError(f"{what}: shorter than its header")
    n, heap_len, nlp_max, beta_max, lse_min, lse_max = _GWAS_HEADER.unpack_from(buf, 0)
    if n < 1:
        raise ValueError(f"{what}: row count {n}")
    if expect_rows is not None and n != expect_rows:
        raise ValueError(f"{what}: {n} rows, index says {expect_rows}")
    if len(buf) != GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + heap_len:
        raise ValueError(f"{what}: payload length {len(buf)} != {GWAS_HEADER_LEN} + "
                         f"{GWAS_ROW_BYTES}n + heap_len")
    for name, v in (("nlp_max", nlp_max), ("beta_max", beta_max), ("lse_min", lse_min), ("lse_max", lse_max)):
        if not math.isfinite(v):
            raise ValueError(f"{what}: {name} is not finite")
    if nlp_max < 0 or beta_max < 0 or lse_min > lse_max:
        raise ValueError(f"{what}: bad scales")
    o, col = GWAS_HEADER_LEN, {}
    for name, t, w in GWAS_COLUMNS:
        col[name] = np.frombuffer(buf, t, n, o)
        o += w * n
    nv = np.asarray(n_values, dtype=np.int64)
    if np.any(col["n_code"] >= nv.size):
        raise ValueError(f"{what}: n_code outside the index's table of {nv.size} values")
    if np.any(col["allele"] > 12) or np.any(col["nlp_code"] == NLP_NULL) or np.any(col["beta_code"] == BETA_NULL):
        raise ValueError(f"{what}: reserved allele, nlp or beta code")
    pos = np.cumsum(col["position_delta"].astype(np.int64))
    pairs = [pf.SNP_ALLELES.get(c, (None, None)) for c in col["allele"].tolist()]
    zero = np.flatnonzero(col["allele"] == 0)
    recs = buf[o:].decode("ascii").split("\n")[:-1] if heap_len else []
    if len(recs) != zero.size:
        raise ValueError(f"{what}: heap holds {len(recs)} records for {zero.size} code-0 rows")
    for i, rec in zip(zero.tolist(), recs):
        pairs[i] = tuple(rec.split("\t"))
    nlp = dequantize_nlp(col["nlp_code"], nlp_max)
    # `rows` is the row count and `n` the per-row sample size, as in v1's GWAS decode. The trans and
    # cis decoders use `n` for the row count, but a GWAS row is the one place both quantities exist,
    # and a caller that got the wrong one here would index a scalar or mask with a count.
    return {"rows": n, "nlp_max": nlp_max, "beta_max": beta_max, "lse_min": lse_min, "lse_max": lse_max,
            "position": pos, "rs_number": col["rs_number"].astype(np.int64),
            "af_code": col["af_code"],
            "af": np.where(col["af_code"] == pf.AF_NULL, np.nan, col["af_code"] / pf.AF_MAXQ),
            "nlp_code": col["nlp_code"], "nlp": nlp, "p": p_from_nlp(nlp),
            "beta_code": col["beta_code"], "beta": dequantize_beta(col["beta_code"], beta_max),
            "se_code": col["se_code"], "se": dequantize_se(col["se_code"], lse_min, lse_max),
            "n_code": col["n_code"], "n": nv[col["n_code"].astype(np.int64)],
            "ref": [x[0] for x in pairs], "alt": [x[1] for x in pairs]}
