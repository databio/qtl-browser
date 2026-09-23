"""Codec for qtlb **v1** (SPEC.md at the repo root): the byte-level pieces every v1 object shares.

The v1 store modules (`qtlstore`, `catalog`, `annotation`, `results`) import this module and never
`packfmt_v0`; the v0 build (`steps_pack*`, `packcheck`, `packtool`, `upload`) imports only
`packfmt_v0`. The two share code by copy for now: where v1 reuses a v0 layout unchanged (variant
pages, result blocks, quantizers, zstd framing, SNP codes) the functions here are the v0 ones,
byte for byte, and a v1 store rebuilt after the split is identical object for object.

What lives here: the quantizers and the slope derivation (SPEC.md section 13), zstd framing (section
1), the variant page header and allele codes (section 5), the results block (section 8), and the
trans frame, hits frame and GWAS block codecs (sections 8-10). The 64-byte v1 file header and object
names are `qtlstore`'s.

Float64 throughout. Only the Python writer rounds (`np.rint`, round half to even); readers only
multiply.
"""
from __future__ import annotations

import json
import math
import struct

import numpy as np
import pyarrow as pa
import zstandard
from scipy.special import gammaln, stdtrit

# ---- quantization codes ---------------------------------------------------------------------
NLP_NULL, NLP_ZERO, NLP_MAXQ = 65535, 65534, 65533      # u16 codes for -log10 p
SE_NULL, SE_SIGN, SE_MAXQ = 0xFFFF, 0x8000, 32766        # u16 SE field: bit 15 slope sign, bits 0-14 log(SE) code
SE_QMASK = 0x7FFF
REL_F32 = 1.2e-7            # relative mismatch of |slope| = se * t in the float32 source (two half-ulps)
BOUND_FACTOR = 1.25         # SPEC.md section 13: a decoded value may miss its source by BOUND_FACTOR x its per-row bound
NLP_LIMIT = 300.0           # keeps inverse-t reference behavior inside the tested double range


def quantize_nlp(p: np.ndarray) -> tuple[np.ndarray, float]:
    """u16 codes and nlp_max for one phenotype. Null or NaN p -> 65535, p == 0 -> 65534,
    else q = rint(-log10(p) / nlp_max * 65533). nlp_max is the largest finite -log10 p
    (0.0 when every p is 1 or null)."""
    p = np.asarray(p, dtype=np.float64)
    null = np.isnan(p)
    zero = p == 0
    ok = ~null & ~zero
    if np.any((p[ok] < 0) | (p[ok] > 1)):
        raise ValueError("p outside [0, 1]")
    nlp = np.zeros(p.shape, dtype=np.float64)
    nlp[ok] = -np.log10(p[ok])
    nlp[ok] = np.maximum(nlp[ok], 0.0)          # -log10(1) is -0.0
    nlp_max = float(nlp[ok].max()) if np.any(ok) else 0.0
    if nlp_max > NLP_LIMIT:
        raise ValueError(f"-log10 p {nlp_max} exceeds the format limit {NLP_LIMIT:g}")
    q = np.zeros(p.shape, dtype=np.uint16)
    if nlp_max > 0:
        q[ok] = np.rint(nlp[ok] / nlp_max * NLP_MAXQ).astype(np.uint16)
    q[zero] = NLP_ZERO
    q[null] = NLP_NULL
    return q, nlp_max


def dequantize_nlp(q: np.ndarray, nlp_max: float) -> np.ndarray:
    """float64 -log10 p = q * nlp_max / 65533; NaN for 65535, +inf for 65534."""
    q = np.asarray(q, dtype=np.uint16)
    out = q.astype(np.float64) * (float(nlp_max) / NLP_MAXQ)
    out[q == NLP_ZERO] = np.inf
    out[q == NLP_NULL] = np.nan
    return out


def p_from_nlp(nlp: np.ndarray) -> np.ndarray:
    """p = 10^-nlp; +inf gives 0, NaN stays NaN."""
    return np.power(10.0, -np.asarray(nlp, dtype=np.float64))


def quantize_se(se: np.ndarray, slope: np.ndarray) -> tuple[np.ndarray, float, float]:
    """u16 SE-field codes and (lse_min, lse_max) for one phenotype.

    A row is null (0xFFFF) when se or slope is null or NaN: without the slope there is no sign.
    Otherwise bit 15 is set when slope < 0 (so -0.0 and 0.0 are positive) and bits 0-14 hold
    q = rint((ln se - lse_min) / (lse_max - lse_min) * 32766), or 0 when lse_min == lse_max.
    lse_min and lse_max are the smallest and largest ln se over non-null rows (both 0.0 when every
    row is null). Raises on a non-null se that is not finite and positive, or an infinite slope."""
    se = np.asarray(se, dtype=np.float64)
    s = np.asarray(slope, dtype=np.float64)
    if se.shape != s.shape:
        raise ValueError(f"se and slope shapes differ: {se.shape} vs {s.shape}")
    if np.any(np.isinf(s)):
        raise ValueError("infinite slope")
    ok = ~np.isnan(se) & ~np.isnan(s)
    if np.any(~np.isfinite(se[ok]) | (se[ok] <= 0)):
        raise ValueError("se must be finite and > 0 where present")
    q = np.full(se.shape, SE_NULL, dtype=np.uint16)
    if not np.any(ok):
        return q, 0.0, 0.0
    lse = np.log(se[ok])
    lse_min, lse_max = float(lse.min()), float(lse.max())
    code = np.zeros(lse.shape, dtype=np.uint16)
    if lse_max > lse_min:
        code = np.rint((lse - lse_min) / (lse_max - lse_min) * SE_MAXQ).astype(np.uint16)
    q[ok] = code | np.where(s[ok] < 0, SE_SIGN, 0).astype(np.uint16)
    return q, lse_min, lse_max


def dequantize_se(q: np.ndarray, lse_min: float, lse_max: float) -> tuple[np.ndarray, np.ndarray]:
    """(se, negative): se = exp(lse_min + code * ((lse_max - lse_min) / 32766)) in float64 with
    code = bits 0-14, NaN for 0xFFFF; negative is bit 15 (False for null rows)."""
    q = np.asarray(q, dtype=np.uint16)
    null = q == SE_NULL
    step = (float(lse_max) - float(lse_min)) / SE_MAXQ
    se = np.exp(float(lse_min) + (q & SE_QMASK).astype(np.float64) * step)
    se[null] = np.nan
    return se, ((q & SE_SIGN) != 0) & ~null


def t_from_p(p: np.ndarray, dof: int) -> np.ndarray:
    """|t| = max(0, -scipy.special.stdtrit(dof, p / 2)); +inf for p = 0, NaN for NaN.
    At p = 1 scipy returns a signed zero or about 6.6e-17 depending on the build; the clamp keeps
    t >= 0 there for every reader."""
    p = np.asarray(p, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        t = np.maximum(-stdtrit(dof, p / 2.0), 0.0)
    return np.where(p == 0, np.inf, t)


def slope_from_se(se: np.ndarray, negative: np.ndarray, p: np.ndarray, dof: int) -> np.ndarray:
    """slope = sign * se * t(p, dof); NaN where se or p is NaN, or p = 0 (t is infinite)."""
    se = np.asarray(se, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    out = np.full(np.broadcast(se, p).shape, np.nan)
    ok = ~np.isnan(se) & ~np.isnan(p) & (p > 0)
    if np.any(ok):
        out[ok] = np.where(np.asarray(negative)[ok], -1.0, 1.0) * se[ok] * t_from_p(p[ok], dof)
    return out


def _log_t_pdf(t: np.ndarray, dof: int) -> np.ndarray:
    nu = float(dof)
    return (gammaln((nu + 1) / 2) - gammaln(nu / 2) - 0.5 * math.log(nu * math.pi)
            - (nu + 1) / 2 * np.log1p(t * t / nu))


def error_bounds(p, slope, se, nlp_max: float, lse_min: float, lse_max: float, dof: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-row worst-case absolute errors (SPEC.md section 13) of the decoded slope_se and the derived
    slope against the source values, from the phenotype's two rounding steps:

        se_bound    = se * (expm1(h_lse) + 1e-12)                    h_lse = (lse_max - lse_min) / 65532
        slope_bound = se_bound * (t + dt) + se * dt + 1.2e-7 * |slope|    dt = h_nlp / D, h_nlp = nlp_max / 131066

    t = |slope| / se from the source; D = d(-log10 p)/dt = 2 f(t; dof) / (p ln 10), f the Student-t
    density (so dt is how far t moves when -log10 p moves half a step). The last term is the
    float32 source: |slope| and se * t agree only to two half-ulps. NaN where p, slope, or se is
    null, se <= 0, or p = 0."""
    p = np.asarray(p, dtype=np.float64)
    s = np.asarray(slope, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    se_b = np.full(p.shape, np.nan)
    sl_b = np.full(p.shape, np.nan)
    ok = ~np.isnan(p) & (p > 0) & ~np.isnan(s) & ~np.isnan(se) & (se > 0)
    if not np.any(ok):
        return se_b, sl_b
    h_lse = (float(lse_max) - float(lse_min)) / (2 * SE_MAXQ)
    h_nlp = float(nlp_max) / (2 * NLP_MAXQ)
    t = np.abs(s[ok]) / se[ok]
    d = 2.0 * np.exp(_log_t_pdf(t, dof) - np.log(p[ok])) / math.log(10.0)
    dt = h_nlp / d
    se_b[ok] = se[ok] * (math.expm1(h_lse) + 1e-12)
    sl_b[ok] = se_b[ok] * (t + dt) + se[ok] * dt + REL_F32 * np.abs(s[ok])
    return se_b, sl_b


# ---- format constants (SPEC.md sections 1, 5, 8) ---------------------------------------------
DETAILS_VERSION = 1                # `v` of every v1 block's details JSON
PAGE_HEADER_LEN, BLOCK_HEADER_LEN = 12, 64
MAGIC_BLOCK = b"QGB0"
CODECS = {"raw": 0, "zstd": 1}
CODEC_NAMES = {v: k for k, v in CODECS.items()}
CS_RECORD_LEN = 12
NO_VAR_START = 0xFFFFFFFF          # block var_start when n_rows is 0
U32_MAX = 0xFFFFFFFF
AF_NULL, AF_MAXQ = 65535, 65534    # af code = rint(af * 65534)
COUNT_NULL = 65535                 # ma_samples, ma_count
MAX_PAGE_RECORDS = 65535           # page header n is a u16

SNP_CODES = {("A", "C"): 1, ("A", "G"): 2, ("A", "T"): 3, ("C", "A"): 4, ("C", "G"): 5, ("C", "T"): 6,
             ("G", "A"): 7, ("G", "C"): 8, ("G", "T"): 9, ("T", "A"): 10, ("T", "C"): 11, ("T", "G"): 12}
SNP_ALLELES = {code: pair for pair, code in SNP_CODES.items()}
MATCH_CODES = {"none": 0, "exact": 1, "position": 2}
MATCH_NAMES = {v: k for k, v in MATCH_CODES.items()}

_PAGE_HEADER = struct.Struct("<IIHBB")
_BLOCK_HEADER = struct.Struct("<4sIIIiIIIdddII")
PAIR_DTYPE = np.dtype([("nlp", "<u2"), ("se", "<u2")])
CS_DTYPE = np.dtype([("row", "<u4"), ("pip", "<f4"), ("cs_id", "u1"), ("pad", "u1", (3,))])
assert _PAGE_HEADER.size == PAGE_HEADER_LEN
assert _BLOCK_HEADER.size == BLOCK_HEADER_LEN and PAIR_DTYPE.itemsize == 4 and CS_DTYPE.itemsize == CS_RECORD_LEN

UNCHECKED = object()   # sentinel for decoder expectations that were not supplied


def pad4(n: int) -> int:
    """Zero bytes needed to bring a length of n to a multiple of 4."""
    return -n % 4


# ---- zstd framing (SPEC.md section 1) ---------------------------------------------------------
_ZC: dict[int, zstandard.ZstdCompressor] = {}
_ZD = zstandard.ZstdDecompressor()


def zstd_frame(data: bytes, level: int) -> bytes:
    """Exactly one zstd frame: Frame_Content_Size present, content checksum on, no dictionary."""
    c = _ZC.get(level)
    if c is None:
        c = _ZC[level] = zstandard.ZstdCompressor(level=level, write_checksum=True, write_content_size=True,
                                                  write_dict_id=False)
    return c.compress(data)


def zstd_unframe(frame: bytes, expect_len: int | None, what: str) -> bytes:
    """Decompress one frame, enforcing the SPEC framing rules; ValueError names the rule."""
    try:
        fp = zstandard.get_frame_parameters(frame)
    except zstandard.ZstdError as e:
        raise ValueError(f"{what}: not a zstd frame ({e})") from None
    if fp.content_size in (zstandard.CONTENTSIZE_UNKNOWN, zstandard.CONTENTSIZE_ERROR):
        raise ValueError(f"{what}: zstd frame has no content size")
    if not fp.has_checksum:
        raise ValueError(f"{what}: zstd frame has no content checksum")
    if fp.dict_id != 0:
        raise ValueError(f"{what}: zstd frame names a dictionary")
    if expect_len is not None and fp.content_size != expect_len:
        raise ValueError(f"{what}: zstd content size {fp.content_size} != expected {expect_len}")
    try:
        out = _ZD.decompress(frame, max_output_size=fp.content_size, allow_extra_data=False)
    except zstandard.ZstdError as e:
        raise ValueError(f"{what}: {e}") from None
    if len(out) != fp.content_size:
        raise ValueError(f"{what}: decompressed {len(out)} bytes, frame says {fp.content_size}")
    return out


# ---- input coercion ---------------------------------------------------------------------------
def _to_numpy(x):
    if isinstance(x, (pa.Array, pa.ChunkedArray)):
        return x.to_numpy(zero_copy_only=False)
    return x


def _int_column(x, name: str, lo: int, hi: int, null_code: int, n: int) -> np.ndarray:
    """int64 array; nulls (None, NaN, masked) become `null_code`; non-null values outside [lo, hi] raise."""
    x = _to_numpy(x)
    if np.ma.isMaskedArray(x):
        mask, a = np.ma.getmaskarray(x), np.ma.getdata(x)
    else:
        a = np.asarray(x)
        mask = np.zeros(a.shape, dtype=bool)
    if a.shape != (n,):
        raise ValueError(f"{name}: expected {n} values, got shape {a.shape}")
    out = np.full(n, null_code, dtype=np.int64)
    if a.dtype.kind == "O":
        for i in range(n):
            v = a[i]
            if mask[i] or v is None or (isinstance(v, float) and math.isnan(v)):
                continue
            if isinstance(v, (bool, np.bool_, str, bytes)) or int(v) != v:
                raise ValueError(f"{name}: non-integer value {v!r}")
            if not lo <= int(v) <= hi:
                raise ValueError(f"{name}: value {v} outside {lo}..{hi}")
            out[i] = int(v)
        return out
    if a.dtype.kind == "f":
        null = mask | np.isnan(a)
        f = a[~null]
        if not np.all(np.isfinite(f)) or np.any(f != np.rint(f)):
            raise ValueError(f"{name}: non-integer value")
    elif a.dtype.kind in "iu":
        null = mask
        f = a[~null]
    else:
        raise ValueError(f"{name}: unsupported dtype {a.dtype}")
    if f.size and (f.min() < lo or f.max() > hi):
        raise ValueError(f"{name}: value outside {lo}..{hi} (min {f.min()}, max {f.max()})")
    out[~null] = f.astype(np.int64)
    return out


def _float_column(x, name: str, n: int | None = None) -> np.ndarray:
    """float64 array; None and masked entries become NaN."""
    x = _to_numpy(x)
    if np.ma.isMaskedArray(x):
        a = np.ma.filled(x.astype(np.float64), np.nan)
    else:
        a = np.asarray(x)
        if a.dtype.kind == "O":
            a = np.array([np.nan if v is None else float(v) for v in a], dtype=np.float64)
        else:
            a = a.astype(np.float64)
    if a.ndim != 1 or (n is not None and a.shape != (n,)):
        raise ValueError(f"{name}: expected {n} values, got shape {a.shape}")
    return a


def _str_list(x, name: str, n: int) -> list:
    if isinstance(x, (pa.Array, pa.ChunkedArray)):
        x = x.to_pylist()
    x = list(x)
    if len(x) != n:
        raise ValueError(f"{name}: expected {n} values, got {len(x)}")
    return x


# ---- results block (SPEC.md section 8) --------------------------------------------------------
def _json_clean(x):
    if isinstance(x, dict):
        return {str(k): _json_clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_json_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_json_clean(v) for v in x.tolist()]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (float, np.floating)):
        f = float(x)
        return None if math.isnan(f) else f
    return x


def details_bytes(details: dict) -> bytes:
    """Gene details JSON: UTF-8 without BOM, no whitespace, NaN -> null, floats in Python's
    shortest round-trip form. Infinite floats raise."""
    return json.dumps(_json_clean(details), allow_nan=False, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def encode_gene_block(details: dict, var_start: int | None, anchor: int | None, p, slope, slope_se,
                      cs_row, cs_pip, cs_id, level: int, *, pos_first: int | None = None,
                      pos_last: int | None = None) -> bytes:
    """Block header + raw pairs + credible-set records + zstd details frame (SPEC.md section 8,
    "Results file"), padded to a multiple of 4. Every v1 block has a details frame.

    `p`, `slope`, and `slope_se` are the phenotype's nominal rows in vidx order (None/NaN = null).
    The pair stores -log10 p and slope_se with the slope's sign; the slope itself is not stored. A
    phenotype with no rows passes empty p/slope/slope_se and None for var_start, anchor, pos_first,
    pos_last. Credible-set memberships are sorted by (`cs_row`, `cs_id`) with unique pairs;
    `cs_row` is the 0-based row within the block."""
    if not isinstance(details, dict):
        raise ValueError("block: a v1 block needs a details object")
    p = _float_column(p, "p")
    n = len(p)
    s = _float_column(slope, "slope", n)
    se = _float_column(slope_se, "slope_se", n)
    if n == 0:
        if any(v is not None for v in (var_start, anchor, pos_first, pos_last)):
            raise ValueError("block: a block with no rows takes var_start, anchor, pos_first, pos_last = None")
        h_var_start, h_anchor, h_first, h_last = NO_VAR_START, 0, 0, 0
    else:
        if any(v is None for v in (var_start, anchor, pos_first, pos_last)):
            raise ValueError("block: var_start, anchor, pos_first, pos_last are required when there are rows")
        if var_start < 0 or var_start + n - 1 >= NO_VAR_START:
            raise ValueError(f"block: rows {var_start}..{var_start + n - 1} do not fit below 0xFFFFFFFF")
        if not -2**31 <= anchor < 2**31:
            raise ValueError(f"block: anchor {anchor} does not fit i32")
        if not 1 <= pos_first <= pos_last <= U32_MAX:
            raise ValueError(f"block: pos_first {pos_first}, pos_last {pos_last} out of order or range")
        h_var_start, h_anchor, h_first, h_last = var_start, anchor, pos_first, pos_last
    nlp_q, nlp_max = quantize_nlp(p)
    se_q, lse_min, lse_max = quantize_se(se, s)
    pairs = np.empty(n, dtype=PAIR_DTYPE)
    pairs["nlp"], pairs["se"] = nlp_q, se_q

    k = 0 if cs_row is None else len(cs_row)
    if k and n == 0:
        raise ValueError("block: credible-set records on a block with no rows")
    cs = np.zeros(k, dtype=CS_DTYPE)
    if k:
        rows = _int_column(cs_row, "cs_row", 0, n - 1, -1, k)
        ids = _int_column(cs_id, "cs_id", 0, 127, -1, k)
        pip = _float_column(cs_pip, "cs_pip", k)
        if np.any(rows < 0) or np.any(ids < 0):
            raise ValueError("block: null credible-set row or cs_id")
        if any((int(rows[i]), int(ids[i])) >= (int(rows[i + 1]), int(ids[i + 1])) for i in range(k - 1)):
            raise ValueError("block: credible-set (row, cs_id) pairs must be strictly ascending")
        if not np.all((pip >= 0) & (pip <= 1)):
            raise ValueError("block: pip outside [0, 1] or null")
        cs["row"], cs["pip"], cs["cs_id"] = rows, pip.astype(np.float32), ids

    dj = details_bytes(details)
    dz = zstd_frame(dj, level)
    body = BLOCK_HEADER_LEN + 4 * n + CS_RECORD_LEN * k + len(dz)
    blk_len = body + pad4(body)
    if blk_len > U32_MAX:
        raise ValueError("block: longer than 4 GiB")
    header = _BLOCK_HEADER.pack(MAGIC_BLOCK, blk_len, n, h_var_start, h_anchor, h_first, h_last, k,
                                nlp_max, lse_min, lse_max, len(dz), len(dj))
    return b"".join((header, pairs.tobytes(), cs.tobytes(), dz, bytes(blk_len - body)))


def _reject_constant(name):
    raise ValueError(f"block details: JSON constant {name} is not allowed")


def decode_gene_block(buf: bytes, dof: int, *, expect_blk_len: int | None = None,
                      expect_n_var=UNCHECKED, expect_var_start=UNCHECKED, details_version: int = DETAILS_VERSION) -> dict:
    """Header fields, details dict, and per-row arrays: nlp_code, se_code, nlp, pval_nominal,
    slope_se (NaN for the null code), negative (the sign bit), and slope derived with
    `slope_from_se` and `dof` (NaN where p or slope_se is null, or p = 0). Credible sets stay sparse:
    cs_row, cs_pip, cs_id. Every v1 block carries a details frame whose `v` is `details_version`.

    Expectations from the search index: `expect_blk_len`; `expect_n_var` (None means the phenotype
    has no rows, so n_rows must be 0); `expect_var_start` (None likewise). Raises ValueError on the
    first rule that fails."""
    buf = bytes(buf)
    total = len(buf)
    if total < BLOCK_HEADER_LEN:
        raise ValueError(f"block: {total} bytes is shorter than the 64-byte header")
    (magic, blk_len, n, var_start, anchor, pos_first, pos_last, n_cs, nlp_max, lse_min, lse_max,
     dzlen, dlen) = _BLOCK_HEADER.unpack_from(buf, 0)
    if magic != MAGIC_BLOCK:
        raise ValueError(f"block: magic {magic!r} is not {MAGIC_BLOCK!r}")
    if blk_len != total:
        raise ValueError(f"block: length field {blk_len} != range length {total}")
    if expect_blk_len is not None and blk_len != expect_blk_len:
        raise ValueError(f"block: length field {blk_len} != search_index blk_len {expect_blk_len}")
    body = BLOCK_HEADER_LEN + 4 * n + CS_RECORD_LEN * n_cs + dzlen
    if body + pad4(body) != blk_len:
        raise ValueError(f"block: 64 + 4*n_rows + 12*n_cs + details_zlen padded to 4 = {body + pad4(body)} != block length {blk_len}")
    if any(buf[body:]):
        raise ValueError("block: padding is not zero")
    if dzlen == 0:
        raise ValueError("block: a v1 block needs a details frame")
    if not (math.isfinite(nlp_max) and nlp_max >= 0 and math.isfinite(lse_min) and math.isfinite(lse_max)
            and lse_min <= lse_max):
        raise ValueError(f"block: scales nlp_max {nlp_max} (finite, >= 0), lse_min {lse_min} <= lse_max {lse_max} (finite) do not hold")
    if n == 0:
        if (var_start, anchor, pos_first, pos_last, n_cs, nlp_max, lse_min, lse_max) != (NO_VAR_START, 0, 0, 0, 0, 0.0, 0.0, 0.0):
            raise ValueError("block: n_rows 0 needs var_start 0xFFFFFFFF, anchor, pos_first, pos_last, n_cs 0 and zero scales")
    else:
        if var_start + n - 1 >= NO_VAR_START:
            raise ValueError(f"block: var_start {var_start} + n_rows reaches 0xFFFFFFFF")
        if not 1 <= pos_first <= pos_last:
            raise ValueError(f"block: pos_first {pos_first}, pos_last {pos_last} out of order")
    if expect_n_var is not UNCHECKED and n != (expect_n_var or 0):
        raise ValueError(f"block: n_rows {n} != search_index n_var {expect_n_var}")
    if expect_var_start is not UNCHECKED and (None if n == 0 else var_start) != expect_var_start:
        raise ValueError(f"block: var_start {var_start} != search_index var_start {expect_var_start}")

    pairs = np.frombuffer(buf, PAIR_DTYPE, n, BLOCK_HEADER_LEN)
    nlp_q, se_q = pairs["nlp"].copy(), pairs["se"].copy()
    fin = nlp_q[nlp_q <= NLP_MAXQ]
    if nlp_max == 0 and np.any(fin != 0):
        raise ValueError("block: nlp_max is 0 but a finite -log10 p code is nonzero")
    if nlp_max > 0 and (fin.size == 0 or int(fin.max()) != NLP_MAXQ):
        raise ValueError("block: nlp_max > 0 but no row holds code 65533")
    lq = (se_q[se_q != SE_NULL] & SE_QMASK).astype(np.int32)
    if np.any(lq > SE_MAXQ):
        raise ValueError("block: SE code 0x7FFF is not allowed (log(SE) codes stop at 32766; 0xFFFF is null)")
    if lq.size == 0 and (lse_min, lse_max) != (0.0, 0.0):
        raise ValueError("block: no row has an SE but lse_min, lse_max are not both 0")
    if lq.size and lse_min == lse_max and np.any(lq != 0):
        raise ValueError("block: lse_min == lse_max but a log(SE) code is nonzero")
    if lse_min < lse_max and (int(lq.min()) != 0 or int(lq.max()) != SE_MAXQ):
        raise ValueError("block: lse_min < lse_max but no row holds log(SE) code 0 or no row holds 32766")

    cs = np.frombuffer(buf, CS_DTYPE, n_cs, BLOCK_HEADER_LEN + 4 * n)
    rows = cs["row"].astype(np.int64)
    if n_cs:
        ids = cs["cs_id"].astype(np.int64)
        if any((int(rows[i]), int(ids[i])) >= (int(rows[i + 1]), int(ids[i + 1])) for i in range(n_cs - 1)):
            raise ValueError("block: credible-set (row, cs_id) pairs are not strictly ascending")
        if rows[-1] >= n:
            raise ValueError(f"block: credible-set row {rows[-1]} >= n_rows {n}")
        if np.any(cs["pad"] != 0):
            raise ValueError("block: credible-set padding bytes are not zero")
        if not np.all((cs["pip"] >= 0) & (cs["pip"] <= 1)):
            raise ValueError("block: pip outside [0, 1]")
        if np.any(cs["cs_id"] > 127):
            raise ValueError("block: cs_id above 127 does not fit the reader's TINYINT")

    dj = zstd_unframe(buf[BLOCK_HEADER_LEN + 4 * n + CS_RECORD_LEN * n_cs:body], dlen, "block details")
    if dj.startswith(b"\xef\xbb\xbf"):
        raise ValueError("block details: JSON starts with a BOM")
    try:
        details = json.loads(dj.decode("utf-8"), parse_constant=_reject_constant)
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        raise ValueError(f"block details: {e}") from None
    if not isinstance(details, dict) or details.get("v") != details_version:
        raise ValueError(f"block details: not an object with v = {details_version}")

    nlp = dequantize_nlp(nlp_q, nlp_max)
    pval = p_from_nlp(nlp)
    se, negative = dequantize_se(se_q, lse_min, lse_max)
    return {
        "blk_len": blk_len, "n_rows": n, "var_start": None if n == 0 else var_start,
        "anchor": None if n == 0 else anchor, "pos_first": pos_first, "pos_last": pos_last,
        "n_cs": n_cs, "nlp_max": nlp_max, "lse_min": lse_min, "lse_max": lse_max, "details_zlen": dzlen,
        "details_len": dlen, "dof": dof, "details": details,
        "nlp_code": nlp_q, "se_code": se_q, "nlp": nlp, "pval_nominal": pval, "slope_se": se,
        "negative": negative, "slope": slope_from_se(se, negative, pval, dof),
        "cs_row": rows, "cs_pip": cs["pip"].copy(), "cs_id": cs["cs_id"].astype(np.int8),
    }


BLOCK_FIELDS = ("blk_len", "n_rows", "var_start", "anchor", "pos_first", "pos_last", "n_cs", "nlp_max", "lse_min", "lse_max",
                "details_zlen", "details_len")


def parse_block_header(buf: bytes, off: int = 0) -> dict:
    """The 64-byte block header at `off` as a dict of BLOCK_FIELDS (SPEC.md section 8), checking only
    the magic; `decode_gene_block` applies every rule. For listings that need the header fields of
    many blocks without decoding their rows."""
    if len(buf) - off < BLOCK_HEADER_LEN:
        raise ValueError(f"block header at byte {off}: truncated")
    vals = _BLOCK_HEADER.unpack_from(buf, off)
    if vals[0] != MAGIC_BLOCK:
        raise ValueError(f"block header at byte {off}: magic {vals[0]!r} is not {MAGIC_BLOCK!r}")
    return dict(zip(BLOCK_FIELDS, vals[1:]))


# ---- hits frames (SPEC.md section 8, "Hits file") ---------------------------------------------
# One 20-byte record per (variant, fact): u32 vidx, u32 ord (search-index row), f32 value, f32 beta,
# u8 kind, u8 cs_id, u8 flags, u8 0. Kinds: 0 lead of a group (value p_perm), 1 credible-set member
# (value pip), 2 trans association (value -log10 p, beta the ALT effect). `beta` is NaN on kinds 0-1.
HIT_LEAD, HIT_CS, HIT_TRANS = 0, 1, 2
HIT_MAX_KIND = 2
HIT_FLAG_SIG = 0x01                   # kind 0: the group passes the experiment's significance rule
HIT_DTYPE = np.dtype([("vidx", "<u4"), ("ord", "<u4"), ("value", "<f4"), ("beta", "<f4"), ("kind", "u1"),
                      ("cs_id", "u1"), ("flags", "u1"), ("pad", "u1")])
HITS_FRAME_VARIANTS = 1024            # variants per hits frame (header page size)
assert HIT_DTYPE.itemsize == 20


def hit_records(vidx, ord_, value, kind, cs_id=0, flags=0, beta=np.nan) -> np.ndarray:
    """Build hit records from columns (scalars broadcast)."""
    vidx = np.asarray(vidx, dtype=np.int64)
    a = np.zeros(len(vidx), dtype=HIT_DTYPE)
    a["vidx"], a["ord"], a["value"], a["beta"] = vidx, ord_, value, beta
    a["kind"], a["cs_id"], a["flags"] = kind, cs_id, flags
    return a


def sort_hits(a: np.ndarray) -> np.ndarray:
    """Records in file order: (vidx, kind, ord, cs_id)."""
    return a[np.lexsort((a["cs_id"], a["ord"], a["kind"], a["vidx"]))]


def encode_hits_body(recs: np.ndarray, n_variants: int, frame_variants: int, level: int) -> tuple[bytes, int]:
    """The bytes after the 64-byte header: u32 frame_off[n_frames + 1], then the frames. Frame g holds the
    records with g*F <= vidx < (g+1)*F as one zstd frame; a frame with no records is zero bytes long, so
    `frame_off[g] == frame_off[g+1]` means "no hits" without a request. Offsets are absolute file offsets.
    Returns (body, n_frames)."""
    from_header = 64
    if frame_variants < 1:
        raise ValueError("hits: frame_variants must be >= 1")
    recs = sort_hits(recs)
    if len(recs) and int(recs["vidx"][-1]) >= n_variants:
        raise ValueError(f"hits: vidx {int(recs['vidx'][-1])} beyond the chromosome's {n_variants} variants")
    if np.any(recs["kind"] > HIT_MAX_KIND) or np.any(recs["pad"]):
        raise ValueError("hits: kind above 2 or nonzero pad")
    n_frames = -(-n_variants // frame_variants)
    cut = np.searchsorted(recs["vidx"], np.arange(n_frames + 1, dtype=np.int64) * frame_variants)
    off = from_header + 4 * (n_frames + 1)
    offs, frames = [off], []
    for g in range(n_frames):
        part = recs[cut[g]:cut[g + 1]]
        f = zstd_frame(part.tobytes(), level) if len(part) else b""
        frames.append(f)
        off += len(f)
        offs.append(off)
    if off > U32_MAX:
        raise ValueError("hits: file over 4 GiB")
    return np.asarray(offs, dtype="<u4").tobytes() + b"".join(frames), n_frames


def hits_frame_table(buf: bytes, n_variants: int, frame_variants: int) -> np.ndarray:
    """frame_off[n_frames + 1] from the bytes after a hits file's header (at least 64 + 4 (n_frames + 1))."""
    n_frames = -(-n_variants // frame_variants)
    offs = np.frombuffer(buf, "<u4", n_frames + 1, 64).astype(np.int64)
    if offs[0] != 64 + 4 * (n_frames + 1) or np.any(np.diff(offs) < 0):
        raise ValueError("hits: frame offset table does not start after itself or decreases")
    return offs


def decode_hits_frame(frame: bytes, first_vidx: int, frame_variants: int, what: str = "hits frame") -> np.ndarray:
    """One frame's records (empty for a zero-length frame), checked: sorted, inside the frame's vidx range,
    kinds 0-2, pad zero, flags only on kind 0, beta NaN except on kind 2."""
    if not frame:
        return np.zeros(0, dtype=HIT_DTYPE)
    raw = zstd_unframe(frame, None, what)
    if len(raw) % HIT_DTYPE.itemsize or not raw:
        raise ValueError(f"{what}: {len(raw)} bytes is not a positive multiple of 20")
    a = np.frombuffer(raw, HIT_DTYPE)
    v = a["vidx"].astype(np.int64)
    if v.min() < first_vidx or v.max() >= first_vidx + frame_variants:
        raise ValueError(f"{what}: vidx outside {first_vidx}..{first_vidx + frame_variants - 1}")
    if np.any(a["kind"] > HIT_MAX_KIND) or np.any(a["pad"]) or np.any((a["flags"] != 0) & (a["kind"] != HIT_LEAD)) \
            or np.any(a["flags"] & (0xFF ^ HIT_FLAG_SIG)):
        raise ValueError(f"{what}: bad kind, flags or pad")
    if np.any(np.isfinite(a["beta"]) & (a["kind"] != HIT_TRANS)):
        raise ValueError(f"{what}: beta set on a lead or credible-set record")
    order = np.lexsort((a["cs_id"], a["ord"], a["kind"], a["vidx"]))
    if np.any(order != np.arange(len(a))):
        raise ValueError(f"{what}: records not sorted by (vidx, kind, ord, cs_id)")
    return a


# ---- trans frames (SPEC.md section 9) -----------------------------------------------------------
MAGIC_TRANS = b"QTT2"
_TRANS_HEADER = struct.Struct("<4sIddII")      # magic, n, nlp_max, beta_max, heap_len, reserved
TRANS_HEADER_LEN = 32
TRANS_ROW_BYTES = 16                           # pos 4, rs 4, af 2, nlp 2, beta 2, ordinal 1, allele 1
BETA_MAXQ = 32767
assert _TRANS_HEADER.size == TRANS_HEADER_LEN


def encode_trans_frame(ordinal, pos, rs_number, af_code, ref, alt, p, beta, level: int) -> bytes:
    """One phenotype's trans rows as one zstd frame. Rows must be sorted by (ordinal, pos, ref, alt).
    `ordinal` is the variant chromosome's 1-based position in the variant catalog's table; `af_code`
    the variant catalog's u16 af code; `p` in (0, 1] (0 allowed: code 65534);
    `beta` the ALT effect, finite. Position is stored as a delta: absolute on the first row and
    wherever the ordinal changes, else the difference from the row before (>= 0)."""
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
    b = _float_column(beta, "beta", n)
    if not np.all(np.isfinite(b)):
        raise ValueError("trans frame: beta must be finite")
    nlp_q, nlp_max = quantize_nlp(_float_column(p, "p", n))
    if np.any(nlp_q == NLP_NULL):
        raise ValueError("trans frame: null p")
    beta_max = float(np.abs(b).max())
    bq = np.zeros(n, dtype="<i2") if beta_max == 0 else np.rint(b / beta_max * BETA_MAXQ).astype("<i2")
    rs = np.asarray(rs_number, dtype=np.int64)
    rs = np.where(rs < 0, 0, rs).astype("<u4")
    codes = np.fromiter((SNP_CODES.get((r, t), 0) for r, t in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(codes == 0):
        r, t = ref[i], alt[i]
        if not (r and t and r.isascii() and t.isascii()) or "\t" in r + t or "\n" in r + t:
            raise ValueError(f"trans frame: allele pair {r!r}/{t!r} is not two non-empty ASCII strings")
        heap.append(f"{r}\t{t}\n")
    heap_b = "".join(heap).encode("ascii")
    payload = b"".join((_TRANS_HEADER.pack(MAGIC_TRANS, n, nlp_max, beta_max, len(heap_b), 0),
                        delta.tobytes(), rs.tobytes(),
                        np.asarray(af_code, dtype="<u2").tobytes(), nlp_q.astype("<u2").tobytes(), bq.tobytes(),
                        ordinal.astype("u1").tobytes(), codes.tobytes(), heap_b))
    return zstd_frame(payload, level)


def decode_trans_frame(frame: bytes, what: str = "trans frame") -> dict:
    """Rows of one trans frame: ordinal, pos, rs_number (0 = none), af (NaN null), nlp, p, beta, ref,
    alt, plus the codes and scales. No vidx: a row's (ordinal, pos, ref, alt) names its site."""
    p = zstd_unframe(frame, None, what)
    if len(p) < TRANS_HEADER_LEN:
        raise ValueError(f"{what}: shorter than its header")
    magic, n, nlp_max, beta_max, heap_len, reserved = _TRANS_HEADER.unpack_from(p, 0)
    if magic != MAGIC_TRANS or reserved or n < 1:
        raise ValueError(f"{what}: bad magic, reserved field or row count")
    if len(p) != TRANS_HEADER_LEN + TRANS_ROW_BYTES * n + heap_len:
        raise ValueError(f"{what}: payload length {len(p)} != 32 + 16n + heap_len")
    if not (math.isfinite(nlp_max) and nlp_max >= 0 and math.isfinite(beta_max) and beta_max >= 0):
        raise ValueError(f"{what}: bad scales")
    o = TRANS_HEADER_LEN
    col = {}
    for name, t, w in (("delta", "<u4", 4), ("rs_number", "<u4", 4), ("af_code", "<u2", 2),
                       ("nlp_code", "<u2", 2), ("beta_code", "<i2", 2), ("ordinal", "u1", 1), ("allele", "u1", 1)):
        col[name] = np.frombuffer(p, t, n, o)
        o += w * n
    ordinal = col["ordinal"].astype(np.int64)
    if ordinal.min() < 1 or np.any(np.diff(ordinal) < 0):
        raise ValueError(f"{what}: chromosome ordinals not >= 1 and non-decreasing")
    if np.any(col["allele"] > 12) or np.any(col["beta_code"] == -32768) or np.any(col["nlp_code"] == NLP_NULL):
        raise ValueError(f"{what}: reserved allele, beta or nlp code")
    new_chrom = np.r_[True, ordinal[1:] != ordinal[:-1]]
    pos = np.empty(n, dtype=np.int64)
    d = col["delta"].astype(np.int64)
    run = np.cumsum(new_chrom) - 1
    starts = np.flatnonzero(new_chrom)
    csum = np.cumsum(np.where(new_chrom, 0, d))
    pos[:] = d[starts][run] + csum - csum[starts][run]
    pairs = [SNP_ALLELES.get(c, (None, None)) for c in col["allele"].tolist()]
    zero = np.flatnonzero(col["allele"] == 0)
    heap = p[o:]
    recs = heap.decode("ascii").split("\n")[:-1] if heap_len else []
    if len(recs) != zero.size:
        raise ValueError(f"{what}: heap holds {len(recs)} records for {zero.size} code-0 rows")
    for i, rec in zip(zero.tolist(), recs):
        pairs[i] = tuple(rec.split("\t"))
    nlp = dequantize_nlp(col["nlp_code"], nlp_max)
    af = np.where(col["af_code"] == AF_NULL, np.nan, col["af_code"] / AF_MAXQ)
    return {"n": n, "nlp_max": nlp_max, "beta_max": beta_max, "ordinal": ordinal,
            "pos": pos, "rs_number": col["rs_number"].astype(np.int64),
            "af_code": col["af_code"], "af": af, "nlp_code": col["nlp_code"], "nlp": nlp, "p": p_from_nlp(nlp),
            "beta_code": col["beta_code"], "beta": col["beta_code"] * (beta_max / BETA_MAXQ),
            "ref": [x[0] for x in pairs], "alt": [x[1] for x in pairs]}


# ---- GWAS blocks and index (SPEC.md section 10) ------------------------------------------------
GWAS_SCALE = 10_000                  # beta, se, eaf are stored as rint(x * 1e4)
GWAS_DECIMAL_TOL = 1e-6              # |x * 1e4 - rint(x * 1e4)| below this counts as at most 4 decimals
GWAS_P_TOL = 1e-12                   # |p_mant * 10^p_exp - p| <= this * p: at most 4 significant digits
GWAS_P_MANT = (1000, 9999)
GWAS_P_EXP_MIN = -128                # i8
GWAS_MAX_N_VALUES = 255              # n codes are u8
GWAS_COLUMNS = (("position_delta", "<u4"), ("beta", "<i4"), ("rs_number", "<u4"), ("se", "<u2"), ("af", "<u2"),
                ("p_mant", "<u2"), ("p_exp", "i1"), ("n_code", "u1"), ("allele", "u1"))
GWAS_ROW_BYTES = sum(np.dtype(t).itemsize for _, t in GWAS_COLUMNS)
GWAS_HEADER_LEN = 8                  # payload: u32 n, u32 heap_len
assert GWAS_ROW_BYTES == 21
_POW10 = np.array([float(f"1e{k}") for k in range(0, 1 - GWAS_P_EXP_MIN)])   # the double nearest 10^k, k = 0..128


def _first(mask: np.ndarray) -> int:
    return int(np.flatnonzero(mask)[0])


def gwas_scaled(x, name: str, lo: int, hi: int) -> np.ndarray:
    """rint(x * 1e4) as int64; ValueError names the first row that is null, has more than 4 decimals, or codes outside lo..hi."""
    x = np.asarray(x, dtype=np.float64)
    s = x * GWAS_SCALE
    bad = ~np.isfinite(s)
    if bad.any():
        raise ValueError(f"gwas {name}: row {_first(bad)} is null or not finite")
    q = np.rint(s)
    bad = np.abs(s - q) >= GWAS_DECIMAL_TOL
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas {name}: row {i} value {x[i]!r} has more than 4 decimals")
    bad = (q < lo) | (q > hi)
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas {name}: row {i} code {int(q[i])} outside {lo}..{hi}")
    return q.astype(np.int64)


def gwas_p_codes(p) -> tuple[np.ndarray, np.ndarray]:
    """p = p_mant * 10^p_exp with p_mant in 1000..9999 (SPEC.md section 10). ValueError names the first row whose p is not in
    (0, 1], or does not fit 4 significant digits within GWAS_P_TOL."""
    p = np.asarray(p, dtype=np.float64)
    bad = ~((p > 0) & (p <= 1))
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas p: row {i} value {p[i]!r} not in (0, 1]")
    e = np.floor(np.log10(p)).astype(np.int64) - 3
    bad = e - 1 < GWAS_P_EXP_MIN
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas p: row {i} value {p[i]!r} is below the i8 exponent range")
    m = np.rint(p * _POW10[-e])
    up = m > GWAS_P_MANT[1]                 # mant rounds to 10000: the next exponent
    e[up] += 1
    m[up] = np.rint(p[up] * _POW10[-e[up]])
    down = m < GWAS_P_MANT[0]               # log10 landed one decade high
    e[down] -= 1
    m[down] = np.rint(p[down] * _POW10[-e[down]])
    bad = (m < GWAS_P_MANT[0]) | (m > GWAS_P_MANT[1])
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas p: row {i} value {p[i]!r} gives mantissa {m[i]}")
    back = m / _POW10[-e]
    bad = np.abs(back - p) > GWAS_P_TOL * p
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas p: row {i} value {p[i]!r} has more than 4 significant digits")
    return m.astype("<u2"), e.astype("i1")


def gwas_p_value(p_mant, p_exp) -> np.ndarray:
    """Reader form: p_mant / 10^-p_exp, dividing by the double nearest the power of ten."""
    return np.asarray(p_mant, dtype=np.float64) / _POW10[-np.asarray(p_exp, dtype=np.int64)]


def gwas_codes(position, beta, se, af, p, rs_number, n, n_values) -> dict[str, np.ndarray]:
    """Whole columns (a chromosome, or one block) to SPEC.md section 10 codes, enforcing the lossless rules. ValueError names
    the column and the first bad row. `n_values` is the index's sorted table of distinct n. Alleles are coded per block."""
    k = len(position)
    pos = np.asarray(position, dtype=np.int64)
    if pos.size and (pos.min() < 1 or pos.max() > U32_MAX):
        raise ValueError("gwas position: outside 1..2^32-1")
    if np.any(np.diff(pos) < 0):
        raise ValueError(f"gwas position: decreasing at row {_first(np.diff(pos) < 0) + 1}")
    rs = np.asarray(rs_number, dtype=np.int64)
    if rs.size != k or (rs.size and (rs.min() < 0 or rs.max() > U32_MAX)):
        raise ValueError("gwas rs_number: outside 0..2^32-1")
    nv = np.asarray(n_values, dtype=np.int64)
    if not 1 <= nv.size <= GWAS_MAX_N_VALUES or np.any(np.diff(nv) <= 0) or nv[0] < 0 or nv[-1] > U32_MAX:
        raise ValueError(f"gwas n values: need 1 to {GWAS_MAX_N_VALUES} ascending distinct u32 values, got {nv.size}")
    nn = np.asarray(n, dtype=np.int64)
    code = np.searchsorted(nv, nn)
    bad = (code >= nv.size) | (nv[np.minimum(code, nv.size - 1)] != nn)
    if bad.any():
        i = _first(bad)
        raise ValueError(f"gwas n: row {i} value {nn[i]} is not in the n table")
    mant, exp = gwas_p_codes(p)
    return {"position": pos, "beta": gwas_scaled(beta, "beta", -(2**31 - 1), 2**31 - 1), "rs_number": rs,
            "se": gwas_scaled(se, "se", 0, 65535), "af": gwas_scaled(af, "af", 0, GWAS_SCALE),
            "p_mant": mant, "p_exp": exp, "n_code": code.astype(np.uint8)}


def encode_gwas_block(codes: dict, ref: list, alt: list, level: int) -> bytes:
    """One GWAS block: `codes` from gwas_codes sliced to the block's rows, as one zstd frame of
    u32 n, u32 heap_len, the GWAS_COLUMNS back to back, then the allele heap (SPEC.md section 5 rules)."""
    pos = np.asarray(codes["position"], dtype=np.int64)
    n = pos.size
    if n < 1:
        raise ValueError("gwas block: needs at least one row")
    if len(ref) != n or len(alt) != n or any(len(np.asarray(codes[c])) != n for c, _ in GWAS_COLUMNS[1:-1]):
        raise ValueError("gwas block: columns differ in length")
    if np.any(np.diff(pos) < 0):
        raise ValueError("gwas block: decreasing position")
    deltas = np.empty(n, dtype="<u4")
    deltas[0] = pos[0]
    deltas[1:] = np.diff(pos)
    allele = np.fromiter((SNP_CODES.get((a, b), 0) for a, b in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(allele == 0).tolist():
        for s, nm in ((ref[i], "ref"), (alt[i], "alt")):
            if not isinstance(s, str) or not s.isascii() or "\t" in s or "\n" in s:
                raise ValueError(f"gwas {nm}: row {i} allele {s!r} is not an ASCII string without tab or newline")
        heap.append(f"{ref[i]}\t{alt[i]}\n")
    heap_b = "".join(heap).encode("ascii")
    cols = {**codes, "position_delta": deltas, "allele": allele}
    payload = b"".join([struct.pack("<II", n, len(heap_b))] + [np.asarray(cols[c]).astype(t).tobytes() for c, t in GWAS_COLUMNS] + [heap_b])
    assert len(payload) == GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + len(heap_b)
    return zstd_frame(payload, level)


def decode_gwas_block(frame: bytes, n_values, *, expect_rows: int | None = None, what: str = "gwas block") -> dict:
    """One block's frame to row arrays: position, ref, alt, rs_number (0 = none), beta, se, af (ALT), p (float64), n; plus the codes."""
    payload = zstd_unframe(frame, None, what)
    if len(payload) < GWAS_HEADER_LEN:
        raise ValueError(f"{what}: payload of {len(payload)} bytes is shorter than its 8-byte header")
    n, heap_len = struct.unpack_from("<II", payload, 0)
    if n < 1:
        raise ValueError(f"{what}: no rows")
    if expect_rows is not None and n != expect_rows:
        raise ValueError(f"{what}: {n} rows, expected {expect_rows}")
    if len(payload) != GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + heap_len:
        raise ValueError(f"{what}: payload length {len(payload)} != 8 + 21n + heap_len = {GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + heap_len}")
    c, off = {}, GWAS_HEADER_LEN
    for name, t in GWAS_COLUMNS:
        c[name] = np.frombuffer(payload, dtype=t, count=n, offset=off)
        off += np.dtype(t).itemsize * n
    pos = np.cumsum(c["position_delta"], dtype=np.int64)
    if pos[0] < 1 or pos[-1] > U32_MAX:
        raise ValueError(f"{what}: position outside 1..2^32-1")
    nv = np.asarray(n_values, dtype=np.int64)
    if np.any(c["allele"] > 12):
        raise ValueError(f"{what}: reserved allele code")
    if np.any(c["af"] > GWAS_SCALE):
        raise ValueError(f"{what}: af code above 10000")
    if np.any((c["p_mant"] < GWAS_P_MANT[0]) | (c["p_mant"] > GWAS_P_MANT[1])):
        raise ValueError(f"{what}: p mantissa outside 1000..9999")
    if np.any((c["p_exp"] > -3) | ((c["p_exp"] == -3) & (c["p_mant"] != 1000))):
        raise ValueError(f"{what}: p above 1")
    if np.any(c["n_code"] >= nv.size):
        raise ValueError(f"{what}: n code beyond the {nv.size}-value n table")
    pairs = [SNP_ALLELES.get(x) for x in c["allele"].tolist()]
    zero = np.flatnonzero(c["allele"] == 0)
    heap = payload[off:]
    if zero.size == 0:
        if heap_len:
            raise ValueError(f"{what}: heap has {heap_len} bytes but no code-0 rows")
    else:
        if max(heap, default=0) >= 0x80 or not heap.endswith(b"\n"):
            raise ValueError(f"{what}: heap is not ASCII records ending in newline")
        recs = heap[:-1].decode("ascii").split("\n")
        if len(recs) != zero.size:
            raise ValueError(f"{what}: heap holds {len(recs)} records for {zero.size} code-0 rows")
        for i, rec in zip(zero.tolist(), recs):
            f = rec.split("\t")
            if len(f) != 2:
                raise ValueError(f"{what}: heap record for row {i} does not hold exactly one tab")
            if (f[0], f[1]) in SNP_CODES:
                raise ValueError(f"{what}: SNP {f[0]}/{f[1]} stored in the heap instead of a code")
            pairs[i] = (f[0], f[1])
    return {"rows": n, "position": pos, "ref": [x[0] for x in pairs], "alt": [x[1] for x in pairs],
            "rs_number": c["rs_number"].astype(np.int64), "beta": c["beta"] / GWAS_SCALE, "se": c["se"] / GWAS_SCALE,
            "af": c["af"] / GWAS_SCALE, "p": gwas_p_value(c["p_mant"], c["p_exp"]), "n": nv[c["n_code"]], "codes": c}


def encode_gwas_index_payload(n_values, chroms: list[tuple[str, np.ndarray, np.ndarray]], first_block: int) -> bytes:
    """The GWAS index payload (before its zstd frame): u32 n_values_count, u32 n_values[], u32 n_chroms, then
    per chromosome an 8-byte name, u32 n_blocks, u32 first_position[n_blocks], u32 end_offset[n_blocks].
    `chroms` is (name, first_position, end_offset) in file order; `first_block` is where block 0 starts
    in every GWAS file (64, the v1 header length)."""
    nv = np.asarray(n_values, dtype=np.int64)
    if not 1 <= nv.size <= GWAS_MAX_N_VALUES or np.any(np.diff(nv) <= 0) or nv[0] < 0 or nv[-1] > U32_MAX:
        raise ValueError("gwas index: n values must be 1 to 255 ascending distinct u32 values")
    parts = [struct.pack("<I", nv.size), nv.astype("<u4").tobytes(), struct.pack("<I", len(chroms))]
    for name, fp, eo in chroms:
        fp, eo = np.asarray(fp, dtype=np.int64), np.asarray(eo, dtype=np.int64)
        nm = name.encode("ascii")
        if not 1 <= len(nm) <= 8 or b"\0" in nm:
            raise ValueError(f"gwas index: chromosome {name!r} must be 1-8 ASCII bytes")
        if fp.size < 1 or fp.size != eo.size:
            raise ValueError(f"gwas index {name}: first_position and end_offset need the same nonzero length")
        if fp[0] < 1 or fp[-1] > U32_MAX or np.any(np.diff(fp) < 0):
            raise ValueError(f"gwas index {name}: first positions must be non-decreasing u32 values from 1")
        if eo[0] <= first_block or eo[-1] > U32_MAX or np.any(np.diff(eo) <= 0):
            raise ValueError(f"gwas index {name}: end offsets must increase from past byte {first_block} and fit u32")
        parts += [nm.ljust(8, b"\0"), struct.pack("<I", fp.size), fp.astype("<u4").tobytes(), eo.astype("<u4").tobytes()]
    return b"".join(parts)


def decode_gwas_index_payload(p: bytes, first_block: int) -> dict:
    """{n_values, chroms: {name: (first_position, end_offset)}} from an index payload."""
    off = 0

    def u32s(k: int) -> np.ndarray:
        nonlocal off
        if off + 4 * k > len(p):
            raise ValueError("gwas index: payload ends early")
        a = np.frombuffer(p, "<u4", k, off).astype(np.int64)
        off += 4 * k
        return a

    nv = u32s(int(u32s(1)[0]))
    if not 1 <= nv.size <= GWAS_MAX_N_VALUES or np.any(np.diff(nv) <= 0):
        raise ValueError("gwas index: n values must be 1 to 255 ascending distinct values")
    n_chroms = int(u32s(1)[0])
    chroms = {}
    for _ in range(n_chroms):
        if off + 8 > len(p):
            raise ValueError("gwas index: payload ends early")
        raw = p[off:off + 8]
        off += 8
        name = raw.rstrip(b"\0")
        if not name or b"\0" in name or max(name) >= 0x80 or name.decode() in chroms:
            raise ValueError(f"gwas index: bad or repeated chromosome name {raw!r}")
        nb = int(u32s(1)[0])
        fp, eo = u32s(nb), u32s(nb)
        if nb < 1 or fp[0] < 1 or np.any(np.diff(fp) < 0) or eo[0] <= first_block or np.any(np.diff(eo) <= 0):
            raise ValueError(f"gwas index {name.decode()}: blocks must be non-empty, positions non-decreasing, offsets increasing")
        chroms[name.decode()] = (fp, eo)
    if off != len(p):
        raise ValueError(f"gwas index: {len(p) - off} bytes after the last chromosome")
    return {"n_values": nv.tolist(), "chroms": chroms}


def gwas_window(first_position, end_offset, lo: int, hi: int, first_block: int) -> tuple[int, int, int, int] | None:
    """The window rule (SPEC.md section 10): (start block, end block, first byte, end byte exclusive) of the
    blocks holding every row with lo <= position <= hi, or None when no block starts at or before hi.
    Readers still filter decoded rows."""
    if lo > hi:
        raise ValueError(f"gwas window: lo {lo} > hi {hi}")
    end = int(np.searchsorted(first_position, hi, side="right")) - 1
    if end < 0:
        return None
    start = max(int(np.searchsorted(first_position, lo, side="left")) - 1, 0)
    return start, end, (first_block if start == 0 else int(end_offset[start - 1])), int(end_offset[end])
