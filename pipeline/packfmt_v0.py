"""Reference codec for the qtlb binary pack format (SPEC.md, version 0).

Not a pipeline step, and it writes no files. `packcheck` uses the quantizers, the slope derivation,
and the per-row error bounds, so its fidelity numbers come from the exact code the spec defines;
`steps_pack`, `steps_pack_trans`, and `steps_gwas` import the encoders, `validate` uses the error
bounds, and `packtool` (the command-line inspector and converter) uses the decoders.

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
BOUND_FACTOR = 1.25         # SPEC section 9: a decoded value may miss its source by BOUND_FACTOR x its per-row bound
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
    """Per-row worst-case absolute errors (SPEC section 9) of the decoded slope_se and the derived
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


# ---- format constants (SPEC.md sections 2-5 and 11) ------------------------------------------
MAGIC_FILE = b"QTLB"
KIND_VARIANTS, KIND_EQTL, KIND_SQTL, KIND_GWAS, KIND_GWAS_INDEX, KIND_TRANS = 1, 2, 3, 4, 5, 6
KIND_HITS, KIND_RSID, KIND_VARIANT_INDEX = 7, 8, 9
RESULT_KINDS = (KIND_EQTL, KIND_SQTL)   # block files: header page size 0
PAGED_KINDS = (KIND_VARIANTS, KIND_GWAS, KIND_GWAS_INDEX, KIND_HITS, KIND_RSID, KIND_VARIANT_INDEX)
# header page size 1..65535: variants per page (1), GWAS rows per block (4, 5), variants per hits frame (7),
# records per rsID block (8), variants per page again (9, the value kind 1 uses)
ZERO_PAGE_KINDS = (*RESULT_KINDS, KIND_TRANS)               # header page size 0
VERSION = 0
FILE_HEADER_LEN, PAGE_HEADER_LEN, BLOCK_HEADER_LEN = 32, 12, 64
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
FLAG_MATCH, FLAG_NO_ALLELES, FLAG_RESERVED = 0x03, 0x04, 0xF8   # variant record flags: bits 0-1, bit 2, bits 3-7

_FILE_HEADER = struct.Struct("<4sBBH8sIIII")
_PAGE_HEADER = struct.Struct("<IIHBB")
_BLOCK_HEADER = struct.Struct("<4sIIIiIIIdddII")
PAIR_DTYPE = np.dtype([("nlp", "<u2"), ("se", "<u2")])
CS_DTYPE = np.dtype([("row", "<u4"), ("pip", "<f4"), ("cs_id", "u1"), ("pad", "u1", (3,))])
assert _FILE_HEADER.size == FILE_HEADER_LEN and _PAGE_HEADER.size == PAGE_HEADER_LEN
assert _BLOCK_HEADER.size == BLOCK_HEADER_LEN and PAIR_DTYPE.itemsize == 4 and CS_DTYPE.itemsize == CS_RECORD_LEN

# gene details JSON (SPEC section 5): key order is fixed here, not taken from the source row
GENE_FIELDS = ("gene_id", "gene_id_version", "symbol", "start", "end", "strand", "tss", "biotype", "tested",
               "num_var", "lead_position", "lead_A1", "lead_A2", "lead_rsid", "lead_af", "lead_tss_distance",
               "slope", "slope_se", "pval_nominal", "pval_perm", "pval_beta", "qval", "is_egene",
               "n_credible_sets", "n_trans_pairs")
SPLICE_FIELDS = ("phenotype_id", "intron_start", "intron_end", "cluster_id", "strand", "num_var", "lead_position",
                 "lead_A1", "lead_A2", "lead_rsid", "lead_af", "lead_tss_distance", "slope", "slope_se",
                 "pval_nominal", "pval_perm", "pval_beta", "qval", "is_sqtl", "n_credible_sets", "blk_off", "blk_len")

# SPEC section 8: the per-gene table a reader produces, row order vidx
READER_SCHEMA = pa.schema([
    ("position", pa.int32()), ("A1", pa.string()), ("A2", pa.string()), ("rs_number", pa.int64()),
    ("tss_distance", pa.int32()), ("af", pa.float32()), ("ma_samples", pa.int16()), ("ma_count", pa.int16()),
    ("pval_nominal", pa.float64()), ("slope", pa.float32()), ("slope_se", pa.float32()),
    ("pip", pa.float32()), ("cs_id", pa.int8()),
])

UNCHECKED = object()   # sentinel for decoder expectations that were not supplied


def pad4(n: int) -> int:
    """Zero bytes needed to bring a length of n to a multiple of 4."""
    return -n % 4


# ---- zstd framing ---------------------------------------------------------------------------
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


# ---- file header (SPEC section 3) -----------------------------------------------------------
def file_header(kind: int, chrom: str, count: int, page_size: int, n_cis: int | None = None) -> bytes:
    """32 bytes: magic QTLB, kind, version 0, header length 32, chromosome (ASCII, zero-padded to 8),
    count (variants, blocks, GWAS rows, index chromosomes, or trans frames), page size (kinds 1, 4, 5) or 0
    (kinds 2, 3, 6), then u32 `n_cis` (kind 1, required: the cis variants, at most count; zero otherwise) and
    4 reserved zero bytes."""
    if kind not in (*PAGED_KINDS, *ZERO_PAGE_KINDS):
        raise ValueError(f"file header: unknown kind {kind}")
    try:
        name = chrom.encode("ascii")
    except UnicodeEncodeError:
        raise ValueError(f"file header: chromosome {chrom!r} is not ASCII") from None
    if not 1 <= len(name) <= 8 or b"\0" in name:
        raise ValueError(f"file header: chromosome {chrom!r} must be 1-8 ASCII bytes without NUL")
    if not 0 <= count <= U32_MAX:
        raise ValueError(f"file header: count {count} does not fit u32")
    if kind in PAGED_KINDS and not 1 <= page_size <= MAX_PAGE_RECORDS:
        raise ValueError(f"file header: page size {page_size} not in 1..65535")
    if kind in ZERO_PAGE_KINDS and page_size != 0:
        raise ValueError("file header: page size must be 0 for a results or trans file (kinds 2, 3, 6)")
    if kind == KIND_VARIANTS:
        if n_cis is None or not 0 <= n_cis <= count:
            raise ValueError(f"file header: a variants file needs n_cis in 0..count ({count}), got {n_cis}")
    elif n_cis:
        raise ValueError(f"file header: n_cis is only for a variants file (kind 1), not kind {kind}")
    return _FILE_HEADER.pack(MAGIC_FILE, kind, VERSION, FILE_HEADER_LEN, name.ljust(8, b"\0"), count, page_size, n_cis or 0, 0)


def parse_file_header(b: bytes) -> dict:
    """Header fields; `n_cis` is the kind 1 value and None for the other kinds (whose field must be zero)."""
    if len(b) < FILE_HEADER_LEN:
        raise ValueError(f"file header: {len(b)} bytes, need 32")
    magic, kind, version, hlen, name, count, page_size, n_cis, reserved = _FILE_HEADER.unpack_from(b, 0)
    if magic != MAGIC_FILE:
        raise ValueError(f"file header: magic {magic!r} is not {MAGIC_FILE!r}")
    if kind not in (*PAGED_KINDS, *ZERO_PAGE_KINDS):
        raise ValueError(f"file header: unknown kind {kind}")
    if version != VERSION:
        raise ValueError(f"file header: version {version}, this reader supports {VERSION}")
    if hlen != FILE_HEADER_LEN:
        raise ValueError(f"file header: header length {hlen} != 32")
    if reserved != 0 or (kind != KIND_VARIANTS and n_cis != 0):
        raise ValueError("file header: reserved bytes are not zero")
    if kind == KIND_VARIANTS and n_cis > count:
        raise ValueError(f"file header: n_cis {n_cis} exceeds count {count}")
    chrom = name.rstrip(b"\0")
    if not chrom or b"\0" in chrom or max(chrom) >= 0x80:
        raise ValueError(f"file header: bad chromosome field {name!r}")
    if kind in PAGED_KINDS and not 1 <= page_size <= MAX_PAGE_RECORDS:
        raise ValueError(f"file header: kind {kind} with page size {page_size}, not 1..65535")
    if kind in ZERO_PAGE_KINDS and page_size != 0:
        raise ValueError("file header: results or trans file with nonzero page size")
    return {"kind": kind, "version": version, "header_len": hlen, "chrom": chrom.decode("ascii"),
            "count": count, "page_size": page_size, "n_cis": n_cis if kind == KIND_VARIANTS else None}


# ---- variants file (SPEC section 4) ---------------------------------------------------------
def encode_variant_page(first_vidx: int, position, rs_number, af, ma_samples, ma_count,
                        A1: list[str], A2: list[str], match: list[str], codec: str, level: int) -> bytes:
    """Page header + payload (SPEC section 'Variants file'), padded to a multiple of 4.

    Nulls: rs_number None/NaN/0 -> 0; af None/NaN -> 65535; counts None/NaN/masked -> 65535. A record
    whose A1 and A2 are both None has flags bit 2 (alleles not reported), allele code 0, and no heap
    record; only trans-only records may (encode_variants_file enforces that).
    Raises on decreasing positions, rs_number > 2^32-1, counts > 65534, af outside [0, 1],
    non-ASCII alleles or alleles holding a tab or newline, one null allele, more than 65535 records,
    unknown match."""
    n = len(position)
    if not 1 <= n <= MAX_PAGE_RECORDS:
        raise ValueError(f"page: {n} records, must be 1..65535")
    if codec not in CODECS:
        raise ValueError(f"page: unknown codec {codec!r}")
    if first_vidx < 0 or first_vidx + n - 1 > U32_MAX:
        raise ValueError(f"page: vidx {first_vidx}..{first_vidx + n - 1} does not fit u32")
    pos = _int_column(position, "position", 1, U32_MAX, -1, n)
    if np.any(pos < 0):
        raise ValueError("position: null")
    d = np.diff(pos)
    if np.any(d < 0):
        raise ValueError("position: decreasing within a page")
    deltas = np.empty(n, dtype="<u4")
    deltas[0] = pos[0]
    deltas[1:] = d
    rs = _int_column(rs_number, "rs_number", 0, U32_MAX, 0, n).astype("<u4")
    a = _float_column(af, "af", n)
    ok = ~np.isnan(a)
    if np.any((a[ok] < 0) | (a[ok] > 1)):
        raise ValueError("af: value outside [0, 1]")
    afq = np.full(n, AF_NULL, dtype="<u2")
    afq[ok] = np.rint(a[ok] * AF_MAXQ).astype("<u2")
    ms = _int_column(ma_samples, "ma_samples", 0, COUNT_NULL - 1, COUNT_NULL, n).astype("<u2")
    mc = _int_column(ma_count, "ma_count", 0, COUNT_NULL - 1, COUNT_NULL, n).astype("<u2")
    a1, a2, m = _str_list(A1, "A1", n), _str_list(A2, "A2", n), _str_list(match, "match", n)
    no_alleles = np.fromiter((x is None and y is None for x, y in zip(a1, a2)), dtype=bool, count=n)
    codes = np.fromiter((SNP_CODES.get((x, y), 0) for x, y in zip(a1, a2)), dtype=np.uint8, count=n)
    heap_parts = []
    for i in np.flatnonzero((codes == 0) & ~no_alleles):
        for s, nm in ((a1[i], "A1"), (a2[i], "A2")):
            if not isinstance(s, str):
                raise ValueError(f"{nm}: record {i} allele is {s!r}, not a string")
            if not s.isascii() or "\t" in s or "\n" in s:
                raise ValueError(f"{nm}: record {i} allele {s!r} is not ASCII without tab or newline")
        heap_parts.append(f"{a1[i]}\t{a2[i]}\n")
    heap = "".join(heap_parts).encode("ascii")
    try:
        flags = np.fromiter((MATCH_CODES[v] for v in m), dtype=np.uint8, count=n)
    except KeyError as e:
        raise ValueError(f"match: unknown value {e.args[0]!r}") from None
    flags |= np.where(no_alleles, FLAG_NO_ALLELES, 0).astype(np.uint8)
    payload = b"".join((struct.pack("<I", len(heap)), deltas.tobytes(), rs.tobytes(), afq.tobytes(),
                        ms.tobytes(), mc.tobytes(), codes.tobytes(), flags.tobytes(), heap))
    stored = payload if codec == "raw" else zstd_frame(payload, level)
    page = _PAGE_HEADER.pack(len(stored), first_vidx, n, CODECS[codec], 0) + stored
    return page + bytes(pad4(len(page)))


def decode_variant_pages(buf: bytes, *, var_start: int | None = None, n_var: int | None = None,
                         pos_first: int | None = None, pos_last: int | None = None, n_cis: int | None = None) -> dict:
    """Walk consecutive pages from byte 0 of `buf`; returns vidx, position, rs_number (masked where
    0, i.e. None), af (NaN for null), ma_samples and ma_count (masked where null), A1, A2 (None where
    flags bit 2 says the alleles are not reported), no_alleles, match, allele_code, pages (offset,
    length, first_vidx, n, codec), codec.

    Without `n_cis` positions must never decrease over the whole range, so a range that crosses
    from the cis section into the trans-only section fails. With `n_cis` (the file header value):
    no page straddles it, positions never decrease within each section and may restart at vidx
    n_cis, and bit 2 appears only at vidx >= n_cis.

    With var_start and n_var (SPEC section 7, step 4): the first page must hold var_start, the last
    page must hold var_start + n_var - 1 (so the range is exactly the covering pages), and with
    pos_first / pos_last the positions at those two rows must match the block header. With n_cis too,
    the run must end before n_cis."""
    buf = bytes(buf)
    total = len(buf)
    if total == 0:
        raise ValueError("variant pages: empty range")
    off, next_vidx, codec0 = 0, None, None
    pages: list[dict] = []
    cols: dict[str, list] = {k: [] for k in ("position", "rs", "af", "ms", "mc", "allele", "flags")}
    alleles: list[tuple[str, str]] = []
    while off < total:
        where = f"page at byte {off}"
        if total - off < PAGE_HEADER_LEN:
            raise ValueError(f"{where}: header truncated")
        stored_len, first_vidx, n, codec, reserved = _PAGE_HEADER.unpack_from(buf, off)
        if reserved != 0:
            raise ValueError(f"{where}: reserved header byte is not zero")
        if codec not in CODEC_NAMES:
            raise ValueError(f"{where}: unknown codec {codec}")
        if codec0 is None:
            codec0 = codec
        elif codec != codec0:
            raise ValueError(f"{where}: codec {codec} differs from the first page's codec {codec0}")
        if n == 0:
            raise ValueError(f"{where}: zero records")
        if next_vidx is not None and first_vidx != next_vidx:
            raise ValueError(f"{where}: first_vidx {first_vidx} does not follow the previous page (expected {next_vidx})")
        if first_vidx + n - 1 > U32_MAX:
            raise ValueError(f"{where}: vidx overflows u32")
        end = off + PAGE_HEADER_LEN + stored_len
        stop = end + pad4(end - off)
        if stop > total:
            raise ValueError(f"{where}: stored length {stored_len} runs past the end of the range")
        if any(buf[end:stop]):
            raise ValueError(f"{where}: padding is not zero")
        stored = buf[off + PAGE_HEADER_LEN:end]
        payload = stored if codec == CODECS["raw"] else zstd_unframe(stored, None, where)
        fixed = 4 + 16 * n
        if len(payload) < fixed:
            raise ValueError(f"{where}: payload {len(payload)} bytes is shorter than 4 + 16n = {fixed}")
        (heap_len,) = struct.unpack_from("<I", payload, 0)
        if len(payload) != fixed + heap_len:
            raise ValueError(f"{where}: payload length {len(payload)} != 4 + 16n + heap_len = {fixed + heap_len}")
        pos = np.cumsum(np.frombuffer(payload, "<u4", n, 4), dtype=np.int64)
        if pos[0] < 1 or pos[-1] > U32_MAX:
            raise ValueError(f"{where}: position outside 1..2^32-1")
        allele = np.frombuffer(payload, "u1", n, 4 + 14 * n)
        flags = np.frombuffer(payload, "u1", n, 4 + 15 * n)
        if np.any(allele > 12):
            raise ValueError(f"{where}: reserved allele code {int(allele.max())}")
        if np.any(flags & FLAG_RESERVED):
            raise ValueError(f"{where}: flags bits 3-7 are not zero")
        if np.any((flags & FLAG_MATCH) == 3):
            raise ValueError(f"{where}: reserved rsID match code 3")
        no_al = (flags & FLAG_NO_ALLELES) != 0
        if np.any(no_al & (allele != 0)):
            raise ValueError(f"{where}: flags bit 2 (alleles not reported) on a record with a nonzero allele code")
        if n_cis is not None:
            if first_vidx < n_cis < first_vidx + n:
                raise ValueError(f"{where}: page {first_vidx}..{first_vidx + n - 1} straddles n_cis {n_cis}")
            if first_vidx < n_cis and np.any(no_al):
                raise ValueError(f"{where}: flags bit 2 (alleles not reported) on a cis record (vidx below n_cis {n_cis})")
        zero = np.flatnonzero((allele == 0) & ~no_al)
        heap = payload[fixed:]
        pairs = [SNP_ALLELES.get(c, (None, None)) for c in allele.tolist()]
        if zero.size == 0:
            if heap_len:
                raise ValueError(f"{where}: heap has {heap_len} bytes but no code-0 records")
        else:
            if max(heap, default=0) >= 0x80 or not heap.endswith(b"\n"):
                raise ValueError(f"{where}: heap is not ASCII records ending in newline")
            recs = heap[:-1].decode("ascii").split("\n")
            if len(recs) != zero.size:
                raise ValueError(f"{where}: heap holds {len(recs)} records for {zero.size} code-0 records")
            for i, rec in zip(zero.tolist(), recs):
                f = rec.split("\t")
                if len(f) != 2:
                    raise ValueError(f"{where}: heap record for row {i} does not hold exactly one tab")
                if (f[0], f[1]) in SNP_CODES:
                    raise ValueError(f"{where}: SNP {f[0]}/{f[1]} stored in the heap instead of a code")
                pairs[i] = (f[0], f[1])
        cols["position"].append(pos)
        cols["rs"].append(np.frombuffer(payload, "<u4", n, 4 + 4 * n))
        cols["af"].append(np.frombuffer(payload, "<u2", n, 4 + 8 * n))
        cols["ms"].append(np.frombuffer(payload, "<u2", n, 4 + 10 * n))
        cols["mc"].append(np.frombuffer(payload, "<u2", n, 4 + 12 * n))
        cols["allele"].append(allele)
        cols["flags"].append(flags)
        alleles.extend(pairs)
        pages.append({"offset": off, "length": stop - off, "first_vidx": first_vidx, "n": n,
                      "codec": CODEC_NAMES[codec]})
        next_vidx = first_vidx + n
        off = stop
    c = {k: np.concatenate(v) for k, v in cols.items()}
    first = pages[0]["first_vidx"]
    vidx = np.arange(first, next_vidx, dtype=np.int64)
    down = np.diff(c["position"]) < 0
    if n_cis is not None:
        down &= vidx[1:] != n_cis                 # positions restart at the first trans-only variant
    if np.any(down):
        raise ValueError("variant pages: positions decrease across pages" + ("" if n_cis is None else " within a section"))
    if var_start is not None or n_var is not None:
        if var_start is None or n_var is None or n_var < 1:
            raise ValueError("variant pages: var_start and n_var >= 1 must be given together")
        last = var_start + n_var - 1
        if n_cis is not None and last >= n_cis:
            raise ValueError(f"variant pages: run {var_start}..{last} reaches n_cis {n_cis}")
        if not first <= var_start < first + pages[0]["n"]:
            raise ValueError(f"variant pages: first page ({first}..{first + pages[0]['n'] - 1}) does not hold var_start {var_start}")
        if not pages[-1]["first_vidx"] <= last < next_vidx:
            raise ValueError(f"variant pages: last page ({pages[-1]['first_vidx']}..{next_vidx - 1}) does not hold var_start + n_var - 1 = {last}")
        if pos_first is not None and c["position"][var_start - first] != pos_first:
            raise ValueError(f"variant pages: position at var_start {c['position'][var_start - first]} != pos_first {pos_first}")
        if pos_last is not None and c["position"][last - first] != pos_last:
            raise ValueError(f"variant pages: position at the last row {c['position'][last - first]} != pos_last {pos_last}")
    return {
        "vidx": vidx.astype(np.uint32),
        "position": c["position"].astype(np.uint32),
        "rs_number": np.ma.MaskedArray(c["rs"], mask=c["rs"] == 0),
        "af": np.where(c["af"] == AF_NULL, np.nan, c["af"] / AF_MAXQ),
        "ma_samples": np.ma.MaskedArray(c["ms"], mask=c["ms"] == COUNT_NULL),
        "ma_count": np.ma.MaskedArray(c["mc"], mask=c["mc"] == COUNT_NULL),
        "A1": [p[0] for p in alleles],
        "A2": [p[1] for p in alleles],
        "no_alleles": (c["flags"] & FLAG_NO_ALLELES) != 0,
        "match": [MATCH_NAMES[f & FLAG_MATCH] for f in c["flags"].tolist()],
        "allele_code": c["allele"],
        "pages": pages,
        "codec": CODEC_NAMES[codec0],
    }


def encode_variants_file(chrom: str, position, rs_number, af, ma_samples, ma_count, A1, A2, match,
                         page_size: int, codec: str, level: int, *, trans_only: dict | None = None) -> tuple[bytes, np.ndarray]:
    """File header + pages of `page_size` records in vidx order. The positional arguments are the
    chromosome's cis variants, sorted by position, A1, A2 byte-wise (checked; no null allele).

    `trans_only` is the trans-only section (SPEC section 4): a dict of `position`, `rs_number`, `af`,
    `A1`, `A2`, `match` for the chromosome's trans-only variants, one per position, positions strictly
    increasing. A1 and A2 are both None where the alleles are not reported (flags bit 2); ma_samples
    and ma_count are written as null. Its records take vidx n_cis.. and start a new page, so the cis
    pages are the same bytes as without it. The header's n_cis is the cis count.

    Returns the file bytes and page start offsets over both sections (n_pages + 1 entries; the last
    is the file length). `variant_range` on these offsets is valid for cis runs only."""
    n = len(position)
    pos = _int_column(position, "position", 1, U32_MAX, -1, n)
    if np.any(pos < 0):
        raise ValueError("position: null")
    a1, a2 = _str_list(A1, "A1", n), _str_list(A2, "A2", n)
    if any(x is None or y is None for x, y in zip(a1, a2)):
        raise ValueError("variants file: a cis variant has a null allele")
    same = np.flatnonzero(np.diff(pos) == 0)
    if np.any(np.diff(pos) < 0):
        raise ValueError("variants file: positions are not ascending")
    for i in same.tolist():
        if not (a1[i], a2[i]) < (a1[i + 1], a2[i + 1]):   # str order is byte order for ASCII
            raise ValueError(f"variants file: rows {i}, {i + 1} at position {pos[i]} are not in strict (A1, A2) order")
    rs = _int_column(rs_number, "rs_number", 0, U32_MAX, 0, n)
    afv = _float_column(af, "af", n)
    counts = []
    for x, nm in ((ma_samples, "ma_samples"), (ma_count, "ma_count")):
        v = _int_column(x, nm, 0, COUNT_NULL - 1, COUNT_NULL, n)
        counts.append(np.ma.MaskedArray(v, mask=v == COUNT_NULL))
    m = _str_list(match, "match", n)
    t = trans_only or {k: [] for k in ("position", "rs_number", "af", "A1", "A2", "match")}
    nt = len(t["position"])
    tpos = _int_column(t["position"], "trans-only position", 1, U32_MAX, -1, nt)
    if np.any(tpos < 0) or np.any(np.diff(tpos) <= 0):
        raise ValueError("variants file: trans-only positions must be non-null and strictly increasing")
    trs = _int_column(t["rs_number"], "trans-only rs_number", 0, U32_MAX, 0, nt)
    taf = _float_column(t["af"], "trans-only af", nt)
    tnull = np.ma.MaskedArray(np.full(nt, COUNT_NULL), mask=np.ones(nt, dtype=bool))
    ta1, ta2, tm = _str_list(t["A1"], "trans-only A1", nt), _str_list(t["A2"], "trans-only A2", nt), _str_list(t["match"], "trans-only match", nt)
    parts = [file_header(KIND_VARIANTS, chrom, n + nt, page_size, n)]
    offsets = [FILE_HEADER_LEN]
    for first, k, cols in ((0, n, (pos, rs, afv, counts[0], counts[1], a1, a2, m)),
                           (n, nt, (tpos, trs, taf, tnull, tnull, ta1, ta2, tm))):
        for s in range(0, k, page_size):
            e = min(k, s + page_size)
            pg = encode_variant_page(first + s, *(c[s:e] for c in cols), codec, level)
            parts.append(pg)
            offsets.append(offsets[-1] + len(pg))
    if offsets[-1] > U32_MAX:
        raise ValueError(f"variants file: {offsets[-1]} bytes exceeds 4 GiB")
    return b"".join(parts), np.array(offsets, dtype=np.int64)


def variant_range(page_offsets: np.ndarray, page_size: int, var_start: int, n_var: int) -> tuple[int, int]:
    """`search_index` var_off, var_len: from the first page holding var_start to the end (padding
    included) of the page holding var_start + n_var - 1."""
    if n_var < 1:
        raise ValueError("variant_range: n_var must be >= 1")
    first, last = var_start // page_size, (var_start + n_var - 1) // page_size
    if last + 1 >= len(page_offsets):
        raise ValueError("variant_range: range runs past the last page")
    return int(page_offsets[first]), int(page_offsets[last + 1] - page_offsets[first])


def decode_variants_file(buf: bytes) -> dict:
    """Whole-file decode for tooling and validate: header, then every page under the header's n_cis
    (the `decode_variant_pages` section rules), checking that pages start at vidx 0 and hold `count`
    records, and that within each section every page but the last holds exactly `page_size` records
    (so the first trans-only variant starts a page)."""
    h = parse_file_header(buf)
    if h["kind"] != KIND_VARIANTS:
        raise ValueError(f"variants file: header kind {h['kind']}")
    if h["count"] == 0:
        if len(buf) != FILE_HEADER_LEN:
            raise ValueError("variants file: bytes after the header of an empty file")
        return {"header": h, "n_cis": 0, "pages": []}
    P, n_cis, count = h["page_size"], h["n_cis"], h["count"]
    out = decode_variant_pages(buf[FILE_HEADER_LEN:], n_cis=n_cis)
    sizes = [p["n"] for p in out["pages"]]
    want = [min(P, n_cis - s) for s in range(0, n_cis, P)] + [min(P, count - s) for s in range(n_cis, count, P)]
    if out["pages"][0]["first_vidx"] != 0 or sizes != want:
        raise ValueError(f"variants file: {len(sizes)} pages of {sum(sizes)} records from vidx {out['pages'][0]['first_vidx']} do not follow "
                         f"header count {count}, n_cis {n_cis}, page size {P} ({len(want)} pages)")
    out["header"] = h
    out["n_cis"] = n_cis
    return out


def walk_variants_file(buf: bytes) -> tuple[dict, list[dict]]:
    """Header and every page's {offset, length, first_vidx, n, codec} from the page headers alone,
    without decompressing anything: for tooling that needs page offsets (a run's covering pages, the
    section boundary at n_cis). Checks the header kind, that first_vidx chains from 0 across both
    sections, and that the records total the header count; `decode_variants_file` is the full check."""
    h = parse_file_header(buf)
    if h["kind"] != KIND_VARIANTS:
        raise ValueError(f"variants file: header kind {h['kind']}")
    out, off, next_vidx = [], FILE_HEADER_LEN, 0
    while off < len(buf):
        if len(buf) - off < PAGE_HEADER_LEN:
            raise ValueError(f"variants file: page header at byte {off} is truncated")
        stored_len, first_vidx, n, codec, reserved = _PAGE_HEADER.unpack_from(buf, off)
        end = off + PAGE_HEADER_LEN + stored_len
        stop = end + pad4(end - off)
        if stop > len(buf) or n == 0 or reserved != 0 or codec not in CODEC_NAMES or first_vidx != next_vidx:
            raise ValueError(f"variants file: bad page header at byte {off} (length {stored_len}, first_vidx {first_vidx}, "
                             f"n {n}, codec {codec}, reserved {reserved})")
        out.append({"offset": off, "length": stop - off, "first_vidx": first_vidx, "n": n, "codec": CODEC_NAMES[codec]})
        next_vidx, off = first_vidx + n, stop
    if next_vidx != h["count"]:
        raise ValueError(f"variants file: pages hold {next_vidx} records, header count {h['count']}")
    return h, out


# ---- trans pack (SPEC section 12) -----------------------------------------------------------
MAGIC_TRANS = b"QTT0"
TRANS_HEADER_LEN = 32
TRANS_VARIANT_CHROMS = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX",)    # variant_chr code = index + 1
STRANDS = ("+", "-")                                                          # strand code = index
BETA_MAXQ = 32767
MAX_TRANS_INTRONS = 255
_TRANS_HEADER = struct.Struct("<4sIIHHdd")
assert _TRANS_HEADER.size == TRANS_HEADER_LEN


def _trans_run_starts(n_e: int, intron: np.ndarray) -> np.ndarray:
    """True where a row starts a run: row 0, the first sQTL row, and each sQTL row whose intron differs from the row before."""
    n = n_e + intron.size
    start = np.zeros(n, dtype=bool)
    if n:
        start[0] = True
    if intron.size:
        start[n_e] = True
        start[n_e + 1:] |= intron[1:] != intron[:-1]
    return start


def encode_trans_frame(qtl_type, variant_chr, position, rs_number, af, pval, beta,
                       intron_start, intron_end, cluster, strand, level: int) -> bytes:
    """One gene's trans rows as one zstd frame (SPEC section 12).

    Rows may come in any order: the frame holds the eQTL rows, then the sQTL rows grouped by intron in
    intron-table order (start, end, cluster, strand), and each run by variant chromosome, then position.
    `qtl_type` is 'e' or 's'; `variant_chr` is 'chr1'..'chr22' or 'chrX'; the intron fields are read
    for sQTL rows only (strand '+' or '-'). nlp_max is the gene's largest -log10 p and beta_max its
    largest |beta|. Raises on no rows, an unknown qtl_type or chromosome, p not in (0, 1], an af that
    is null or outside [0, 1], a non-finite beta, a null intron field on an sQTL row, more than 255
    introns, or two rows with the same run, chromosome, and position."""
    n = len(position)
    if n < 1:
        raise ValueError("trans frame: needs at least one row")
    qt = _str_list(qtl_type, "qtl_type", n)
    if any(q not in ("e", "s") for q in qt):
        raise ValueError("trans frame: qtl_type must be 'e' or 's'")
    is_s = np.array([q == "s" for q in qt], dtype=bool)
    code_of = {c: i + 1 for i, c in enumerate(TRANS_VARIANT_CHROMS)}
    try:
        chr_code = np.array([code_of[v] for v in _str_list(variant_chr, "variant_chr", n)], dtype=np.int64)
    except (KeyError, TypeError):
        raise ValueError("trans frame: variant_chr must be chr1..chr22 or chrX") from None
    pos = _int_column(position, "position", 1, U32_MAX, -1, n)
    if np.any(pos < 0):
        raise ValueError("position: null")
    rs = _int_column(rs_number, "rs_number", 0, U32_MAX, 0, n)
    a = _float_column(af, "af", n)
    if not np.all((a >= 0) & (a <= 1)):
        raise ValueError("trans frame: af must be non-null and in [0, 1]")
    p = _float_column(pval, "pval", n)
    if not np.all((p > 0) & (p <= 1)):
        raise ValueError("trans frame: p must be in (0, 1]")
    b = _float_column(beta, "beta", n)
    if not np.all(np.isfinite(b)):
        raise ValueError("trans frame: beta must be finite")
    istart = _int_column(intron_start, "intron_start", 0, U32_MAX, -1, n)
    iend = _int_column(intron_end, "intron_end", 0, U32_MAX, -1, n)
    iclu = _int_column(cluster, "cluster", 0, U32_MAX, -1, n)
    istr = np.array([STRANDS.index(x) if x in STRANDS else -1 for x in _str_list(strand, "strand", n)], dtype=np.int64)
    if np.any(is_s & ((istart < 0) | (iend < 0) | (iclu < 0) | (istr < 0))):
        raise ValueError("trans frame: an sQTL row has a null intron field or a strand other than '+' or '-'")
    si = np.flatnonzero(is_s)
    n_s, n_e = int(si.size), n - int(si.size)
    table, inv = np.unique(np.stack([istart[si], iend[si], iclu[si], istr[si]], axis=1), axis=0, return_inverse=True)
    k = len(table)
    if k > MAX_TRANS_INTRONS:
        raise ValueError(f"trans frame: {k} introns, at most {MAX_TRANS_INTRONS}")
    iidx = np.zeros(n, dtype=np.int64)
    iidx[si] = np.asarray(inv).reshape(-1)
    order = np.lexsort((pos, chr_code, iidx, is_s))
    pos, rs, a, p, b, chr_code, iidx = (x[order] for x in (pos, rs, a, p, b, chr_code, iidx))
    seg = _trans_run_starts(n_e, iidx[n_e:])
    seg[1:] |= chr_code[1:] != chr_code[:-1]
    d = np.diff(pos)
    if np.any(~seg[1:] & (d <= 0)):
        raise ValueError("trans frame: positions must be strictly increasing within a run and chromosome")
    pcol = np.where(seg, pos, np.concatenate(([0], d)))
    nlp_q, nlp_max = quantize_nlp(p)
    beta_max = float(np.abs(b).max())
    bq = np.zeros(n, dtype="<i2") if beta_max == 0 else np.rint(b / beta_max * BETA_MAXQ).astype("<i2")
    payload = b"".join((
        _TRANS_HEADER.pack(MAGIC_TRANS, n_e, n_s, k, 0, nlp_max, beta_max),
        table[:, 0].astype("<u4").tobytes(), table[:, 1].astype("<u4").tobytes(), table[:, 2].astype("<u4").tobytes(),
        table[:, 3].astype("u1").tobytes(), bytes(pad4(k)),
        pcol.astype("<u4").tobytes(), rs.astype("<u4").tobytes(), np.rint(a * AF_MAXQ).astype("<u2").tobytes(),
        nlp_q.astype("<u2").tobytes(), bq.tobytes(), bytes(pad4(14 * n)),
        chr_code.astype("u1").tobytes(), iidx[n_e:].astype("u1").tobytes(), bytes(pad4(n + n_s))))
    return zstd_frame(payload, level)


def decode_trans_frame(buf: bytes, dof_e: int, dof_s: int, what: str = "trans frame") -> dict:
    """One trans frame under SPEC section 12's reader rules (ValueError names the rule). Returns n_e, n_s,
    k, nlp_max, beta_max, the intron table (intron_start, intron_end, cluster, strand as '+'/'-'), and per
    row: qtl_type ('e'/'s'), variant_chr (codes 1..23), position (absolute), rs_number (0 = none), af_code,
    af, nlp_code, nlp, pval, beta_code, beta, t, beta_se = |beta| / t, r2 = t^2 / (t^2 + dof); and intron
    (the n_s sQTL rows' table index). t = max(0, -stdtrit(dof, p / 2)) with dof_e for eQTL rows, dof_s for sQTL."""
    raw = zstd_unframe(bytes(buf), None, what)
    if len(raw) < TRANS_HEADER_LEN:
        raise ValueError(f"{what}: {len(raw)} bytes is shorter than the 32-byte header")
    magic, n_e, n_s, k, reserved, nlp_max, beta_max = _TRANS_HEADER.unpack_from(raw, 0)
    if magic != MAGIC_TRANS:
        raise ValueError(f"{what}: magic {magic!r} is not {MAGIC_TRANS!r}")
    if reserved:
        raise ValueError(f"{what}: reserved header field is not zero")
    n = n_e + n_s
    if n == 0:
        raise ValueError(f"{what}: no rows")
    if (n_s == 0 and k != 0) or (n_s > 0 and not 1 <= k <= MAX_TRANS_INTRONS):
        raise ValueError(f"{what}: k = {k} introns with {n_s} sQTL rows (k must be 0 without sQTL rows, else 1..255)")
    if not (math.isfinite(nlp_max) and nlp_max >= 0 and math.isfinite(beta_max) and beta_max >= 0):
        raise ValueError(f"{what}: scales nlp_max {nlp_max} and beta_max {beta_max} must be finite and >= 0")
    r = TRANS_HEADER_LEN + 13 * k + pad4(k)
    c = r + 14 * n + pad4(14 * n)
    end = c + n + n_s + pad4(n + n_s)
    if len(raw) != end:
        raise ValueError(f"{what}: decompressed length {len(raw)} != {end} for n_e {n_e}, n_s {n_s}, k {k}")
    if any(raw[TRANS_HEADER_LEN + 13 * k:r]) or any(raw[r + 14 * n:c]) or any(raw[c + n + n_s:end]):
        raise ValueError(f"{what}: padding is not zero")
    istart = np.frombuffer(raw, "<u4", k, 32)
    iend = np.frombuffer(raw, "<u4", k, 32 + 4 * k)
    clu = np.frombuffer(raw, "<u4", k, 32 + 8 * k)
    strand = np.frombuffer(raw, "u1", k, 32 + 12 * k)
    if np.any(strand > 1):
        raise ValueError(f"{what}: strand code above 1")
    keys = list(zip(istart.tolist(), iend.tolist(), clu.tolist(), strand.tolist()))
    if any(keys[i] >= keys[i + 1] for i in range(k - 1)):
        raise ValueError(f"{what}: intron table is not strictly ascending by (start, end, cluster, strand)")
    pcol = np.frombuffer(raw, "<u4", n, r).astype(np.int64)
    rs = np.frombuffer(raw, "<u4", n, r + 4 * n)
    afq = np.frombuffer(raw, "<u2", n, r + 8 * n)
    nq = np.frombuffer(raw, "<u2", n, r + 10 * n)
    bq = np.frombuffer(raw, "<i2", n, r + 12 * n)
    chr_code = np.frombuffer(raw, "u1", n, c)
    intron = np.frombuffer(raw, "u1", n_s, c + n)
    if np.any(nq > NLP_MAXQ):
        raise ValueError(f"{what}: nlp code above 65533")
    if np.any(afq > AF_MAXQ):
        raise ValueError(f"{what}: af code above 65534")
    if np.any(bq == -32768):
        raise ValueError(f"{what}: beta code -32768")
    if np.any((chr_code < 1) | (chr_code > 23)):
        raise ValueError(f"{what}: variant_chr code outside 1..23")
    if np.any(intron >= k):
        raise ValueError(f"{what}: intron index at or above k = {k}")
    if n_s and (np.any(np.diff(intron.astype(np.int64)) < 0) or np.unique(intron).size != k):
        raise ValueError(f"{what}: intron indices must never decrease and must use every intron of the table")
    seg = _trans_run_starts(n_e, intron)
    if np.any(~seg[1:] & (np.diff(chr_code.astype(np.int64)) < 0)):
        raise ValueError(f"{what}: variant_chr decreases inside a run")
    seg[1:] |= chr_code[1:] != chr_code[:-1]
    if np.any(pcol == 0):
        raise ValueError(f"{what}: a position entry is 0 (positions start at 1 and increase within a run and chromosome)")
    csum = np.cumsum(pcol)
    starts = np.flatnonzero(seg)
    pos = csum - np.repeat(csum[starts] - pcol[starts], np.diff(np.r_[starts, n]))
    if pos.max() > U32_MAX:
        raise ValueError(f"{what}: position above 2^32-1")
    nlp = nq.astype(np.float64) * (nlp_max / NLP_MAXQ)
    pval = np.power(10.0, -nlp)
    beta = bq.astype(np.float64) * (beta_max / BETA_MAXQ)
    dof = np.where(np.arange(n) < n_e, float(dof_e), float(dof_s))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.maximum(-stdtrit(dof, pval / 2.0), 0.0)
        beta_se = np.abs(beta) / t
    return {"n_e": n_e, "n_s": n_s, "k": k, "nlp_max": nlp_max, "beta_max": beta_max,
            "intron_start": istart.copy(), "intron_end": iend.copy(), "cluster": clu.copy(), "strand": [STRANDS[s] for s in strand.tolist()],
            "qtl_type": ["e"] * n_e + ["s"] * n_s, "variant_chr": chr_code.copy(), "position": pos.astype(np.uint32),
            "rs_number": rs.copy(), "af_code": afq.copy(), "af": afq / AF_MAXQ, "nlp_code": nq.copy(), "nlp": nlp, "pval": pval,
            "beta_code": bq.copy(), "beta": beta, "t": t, "beta_se": beta_se, "r2": t * t / (t * t + dof), "intron": intron.copy(),
            "payload_len": len(raw)}


def trans_phenotype_ids(frame: dict, gene_chr: str, gene_id: str, gene_version: int) -> list[str]:
    """Each decoded row's phenotype_id: the gene_id for eQTL rows, `chr:start:end:clu_<cluster>_<strand>:gene_id.version` for sQTL rows."""
    introns = [f"{gene_chr}:{s}:{e}:clu_{c}_{st}:{gene_id}.{gene_version}" for s, e, c, st in
               zip(frame["intron_start"].tolist(), frame["intron_end"].tolist(), frame["cluster"].tolist(), frame["strand"])]
    return [gene_id] * frame["n_e"] + [introns[i] for i in frame["intron"].tolist()]


# ---- eQTL results pack (SPEC section 5) -----------------------------------------------------
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


def details_from_tables(gene: dict, exons: list[tuple[int, int]], splice_rows: list[dict]) -> dict:
    """Build details directly from genes, collapsed exons, and splice phenotypes.

    Unknown or missing fields raise so schema changes cannot be dropped silently.
    """
    extra = set(gene) - set(GENE_FIELDS) - {"chr", "bin"}
    missing = [k for k in GENE_FIELDS if k not in gene]
    if extra or missing:
        raise ValueError(f"gene row: unexpected columns {sorted(extra)}, missing {missing}")
    splice = []
    for s in splice_rows:
        extra = set(s) - set(SPLICE_FIELDS)
        missing = [k for k in SPLICE_FIELDS if k not in s]
        if extra or missing:
            raise ValueError(f"splice row: unexpected fields {sorted(extra)}, missing {missing}")
        splice.append({k: s[k] for k in SPLICE_FIELDS})
    return {"v": VERSION,
            "gene": {k: gene[k] for k in GENE_FIELDS},
            "exons": [[int(start), int(end)] for start, end in exons],
            "splice": splice}


def encode_gene_block(details: dict | None, var_start: int | None, anchor: int | None, p, slope, slope_se,
                      cs_row, cs_pip, cs_id, level: int, *, pos_first: int | None = None,
                      pos_last: int | None = None) -> bytes:
    """Block header + raw pairs + credible-set records + zstd details frame (SPEC section 'eQTL
    results pack'), padded to a multiple of 4. `details=None` writes a kind 3 (sQTL) block: no
    details frame, details_zlen = details_len = 0, and at least one row.

    `p`, `slope`, and `slope_se` are the gene's nominal rows in vidx order (None/NaN = null). The
    pair stores -log10 p and slope_se with the slope's sign; the slope itself is not stored. A gene
    with no eQTL rows passes empty p/slope/slope_se and None for var_start, anchor, pos_first,
    pos_last. Credible-set memberships are sorted by (`cs_row`, `cs_id`) with unique pairs;
    `cs_row` is the 0-based row within the gene."""
    p = _float_column(p, "p")
    n = len(p)
    s = _float_column(slope, "slope", n)
    se = _float_column(slope_se, "slope_se", n)
    if details is None and n == 0:
        raise ValueError("block: an sQTL block (no details) needs at least one row")
    if n == 0:
        if any(v is not None for v in (var_start, anchor, pos_first, pos_last)):
            raise ValueError("block: a gene with no rows takes var_start, anchor, pos_first, pos_last = None")
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
        raise ValueError("block: credible-set records on a gene with no rows")
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

    dj = b"" if details is None else details_bytes(details)
    dz = b"" if details is None else zstd_frame(dj, level)
    body = BLOCK_HEADER_LEN + 4 * n + CS_RECORD_LEN * k + len(dz)
    blk_len = body + pad4(body)
    if blk_len > U32_MAX:
        raise ValueError("block: longer than 4 GiB")
    header = _BLOCK_HEADER.pack(MAGIC_BLOCK, blk_len, n, h_var_start, h_anchor, h_first, h_last, k,
                                nlp_max, lse_min, lse_max, len(dz), len(dj))
    return b"".join((header, pairs.tobytes(), cs.tobytes(), dz, bytes(blk_len - body)))


def _reject_constant(name):
    raise ValueError(f"block details: JSON constant {name} is not allowed")


def decode_gene_block(buf: bytes, dof: int, *, kind: int = KIND_EQTL, expect_blk_len: int | None = None,
                      expect_n_var=UNCHECKED, expect_var_start=UNCHECKED, details_version: int = VERSION) -> dict:
    """Header fields, details dict, and per-row arrays: nlp_code, se_code, nlp, pval_nominal,
    slope_se (NaN for the null code), negative (the sign bit), and slope derived with
    `slope_from_se` and `dof` (NaN where p or slope_se is null, or p = 0). tss_distance needs
    positions, so the caller gets `anchor`. Credible sets stay sparse: cs_row, cs_pip, cs_id.

    `kind` is the file's kind: a kind 2 (eQTL) block must carry a details frame, and a kind 3
    (sQTL) block must have none (details_zlen = details_len = 0, details returned as None) and at
    least one row.

    Expectations from `search_index` (SPEC section 7, step 3): `expect_blk_len`; `expect_n_var`
    (None means the gene is not eQTL-tested, so n_rows must be 0); `expect_var_start` (None when
    not eQTL-tested). `details_version` is the `v` the details JSON must carry: 0 here, 1 for the
    qtlb v1 results objects (`pipeline/results.py`), which reuse this block layout. Raises ValueError
    on the first rule that fails."""
    if kind not in RESULT_KINDS:
        raise ValueError(f"block: kind {kind} is not a results kind")
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
    if kind == KIND_EQTL and dzlen == 0:
        raise ValueError("block: a kind 2 (eQTL) block needs a details frame")
    if kind == KIND_SQTL and (dzlen or dlen or n == 0):
        raise ValueError(f"block: a kind 3 (sQTL) block has no details frame and at least one row (details_zlen {dzlen}, details_len {dlen}, n_rows {n})")
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

    details = None
    if kind == KIND_EQTL:
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


def encode_eqtl_file(chrom: str, blocks: list[bytes], kind: int) -> tuple[bytes, np.ndarray]:
    """File header of `kind` (2 eQTL, 3 sQTL) + blocks back to back. Returns the bytes and block
    start offsets (len(blocks) + 1 entries; the last is the file length): eQTL blk_off in
    search_index, sQTL blk_off in the gene details."""
    if kind not in RESULT_KINDS:
        raise ValueError(f"results file: kind {kind} is not a results kind")
    offsets = [FILE_HEADER_LEN]
    for b in blocks:
        if len(b) % 4:
            raise ValueError("results file: block length is not a multiple of 4")
        offsets.append(offsets[-1] + len(b))
    if offsets[-1] > U32_MAX:
        raise ValueError(f"results file: {offsets[-1]} bytes exceeds 4 GiB")
    return file_header(kind, chrom, len(blocks), 0) + b"".join(blocks), np.array(offsets, dtype=np.int64)


def walk_eqtl_file(buf: bytes, kind: int) -> tuple[dict, list[tuple[int, int]]]:
    """Header and (blk_off, blk_len) of every block, following block length fields; checks the
    header kind, magic, the block count, and that the last block ends at the end of the file."""
    h = parse_file_header(buf)
    if h["kind"] != kind or kind not in RESULT_KINDS:
        raise ValueError(f"results file: header kind {h['kind']}, expected {kind}")
    out, off = [], FILE_HEADER_LEN
    while off < len(buf):
        if len(buf) - off < 8 or buf[off:off + 4] != MAGIC_BLOCK:
            raise ValueError(f"results file: no block magic at byte {off}")
        (ln,) = struct.unpack_from("<I", buf, off + 4)
        if ln < BLOCK_HEADER_LEN or ln % 4 or off + ln > len(buf):
            raise ValueError(f"results file: bad block length {ln} at byte {off}")
        out.append((off, ln))
        off += ln
    if len(out) != h["count"]:
        raise ValueError(f"results file: {len(out)} blocks, header count {h['count']}")
    return h, out


BLOCK_FIELDS = ("blk_len", "n_rows", "var_start", "anchor", "pos_first", "pos_last", "n_cs", "nlp_max", "lse_min", "lse_max",
                "details_zlen", "details_len")


def parse_block_header(buf: bytes, off: int = 0) -> dict:
    """The 64-byte block header at `off` as a dict of BLOCK_FIELDS (SPEC section 5), checking only
    the magic; `decode_gene_block` applies every rule. For listings that need the header fields of
    many blocks without decoding their rows."""
    if len(buf) - off < BLOCK_HEADER_LEN:
        raise ValueError(f"block header at byte {off}: truncated")
    vals = _BLOCK_HEADER.unpack_from(buf, off)
    if vals[0] != MAGIC_BLOCK:
        raise ValueError(f"block header at byte {off}: magic {vals[0]!r} is not {MAGIC_BLOCK!r}")
    return dict(zip(BLOCK_FIELDS, vals[1:]))


# ---- GWAS pack and index (SPEC section 11) ---------------------------------------------------
GWAS_SCALE = 10_000                  # beta, se, eaf are stored as rint(x * 1e4)
GWAS_DECIMAL_TOL = 1e-6              # |x * 1e4 - rint(x * 1e4)| below this counts as at most 4 decimals
GWAS_P_TOL = 1e-12                   # |p_mant * 10^p_exp - p| <= this * p: at most 4 significant digits
GWAS_P_MANT = (1000, 9999)
GWAS_P_EXP_MIN = -128                # i8
GWAS_MAX_N_VALUES = 255              # n codes are u8
GWAS_INDEX_CHROM = "all"
GWAS_COLUMNS = (("position_delta", "<u4"), ("beta", "<i4"), ("rs_number", "<u4"), ("se", "<u2"), ("eaf", "<u2"),
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
    """p = p_mant * 10^p_exp with p_mant in 1000..9999 (SPEC section 11). ValueError names the first row whose p is not in
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


def gwas_codes(position, beta, se, eaf, p, rs_number, n, n_values) -> dict[str, np.ndarray]:
    """Whole columns (a chromosome, or one block) to SPEC section 11 codes, enforcing the lossless rules. ValueError names
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
            "se": gwas_scaled(se, "se", 0, 65535), "eaf": gwas_scaled(eaf, "eaf", 0, GWAS_SCALE),
            "p_mant": mant, "p_exp": exp, "n_code": code.astype(np.uint8)}


def encode_gwas_block(codes: dict, ea: list, nea: list, level: int) -> bytes:
    """One GWAS block: `codes` from gwas_codes sliced to the block's rows, as one zstd frame of
    u32 n, u32 heap_len, the GWAS_COLUMNS back to back, then the allele heap (SPEC section 4 rules)."""
    pos = np.asarray(codes["position"], dtype=np.int64)
    n = pos.size
    if n < 1:
        raise ValueError("gwas block: needs at least one row")
    if len(ea) != n or len(nea) != n or any(len(np.asarray(codes[c])) != n for c, _ in GWAS_COLUMNS[1:-1]):
        raise ValueError("gwas block: columns differ in length")
    if np.any(np.diff(pos) < 0):
        raise ValueError("gwas block: decreasing position")
    deltas = np.empty(n, dtype="<u4")
    deltas[0] = pos[0]
    deltas[1:] = np.diff(pos)
    allele = np.fromiter((SNP_CODES.get((a, b), 0) for a, b in zip(ea, nea)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(allele == 0).tolist():
        for s, nm in ((ea[i], "ea"), (nea[i], "nea")):
            if not isinstance(s, str) or not s.isascii() or "\t" in s or "\n" in s:
                raise ValueError(f"gwas {nm}: row {i} allele {s!r} is not an ASCII string without tab or newline")
        heap.append(f"{ea[i]}\t{nea[i]}\n")
    heap_b = "".join(heap).encode("ascii")
    cols = {**codes, "position_delta": deltas, "allele": allele}
    payload = b"".join([struct.pack("<II", n, len(heap_b))] + [np.asarray(cols[c]).astype(t).tobytes() for c, t in GWAS_COLUMNS] + [heap_b])
    assert len(payload) == GWAS_HEADER_LEN + GWAS_ROW_BYTES * n + len(heap_b)
    return zstd_frame(payload, level)


def decode_gwas_block(frame: bytes, n_values, *, expect_rows: int | None = None, what: str = "gwas block") -> dict:
    """One block's frame to row arrays: position, ea, nea, rs_number (0 = none), beta, se, eaf, p (float64), n; plus the codes."""
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
    if np.any(c["eaf"] > GWAS_SCALE):
        raise ValueError(f"{what}: eaf code above 10000")
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
    return {"rows": n, "position": pos, "ea": [x[0] for x in pairs], "nea": [x[1] for x in pairs],
            "rs_number": c["rs_number"].astype(np.int64), "beta": c["beta"] / GWAS_SCALE, "se": c["se"] / GWAS_SCALE,
            "eaf": c["eaf"] / GWAS_SCALE, "p": gwas_p_value(c["p_mant"], c["p_exp"]), "n": nv[c["n_code"]], "codes": c}


def encode_gwas_index(n_values, chroms: list[tuple[str, np.ndarray, np.ndarray]], block_rows: int, level: int) -> bytes:
    """gwas_index.bin: a kind 5 header (chromosome `all`, count = chromosomes, page size = rows per block), then one zstd
    frame: u32 n_values_count, u32 n_values[], u32 n_chroms, then per chromosome an 8-byte name, u32 n_blocks,
    u32 first_position[n_blocks], u32 end_offset[n_blocks]. `chroms` is (name, first_position, end_offset) in file order."""
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
        if eo[0] <= FILE_HEADER_LEN or eo[-1] > U32_MAX or np.any(np.diff(eo) <= 0):
            raise ValueError(f"gwas index {name}: end offsets must increase from past byte 32 and fit u32")
        parts += [nm.ljust(8, b"\0"), struct.pack("<I", fp.size), fp.astype("<u4").tobytes(), eo.astype("<u4").tobytes()]
    return file_header(KIND_GWAS_INDEX, GWAS_INDEX_CHROM, len(chroms), block_rows) + zstd_frame(b"".join(parts), level)


def decode_gwas_index(buf: bytes) -> dict:
    """{block_rows, n_values, chroms: {name: (first_position, end_offset)}} from gwas_index.bin."""
    h = parse_file_header(buf)
    if h["kind"] != KIND_GWAS_INDEX or h["chrom"] != GWAS_INDEX_CHROM:
        raise ValueError(f"gwas index: header kind {h['kind']} chromosome {h['chrom']!r}, expected 5 and 'all'")
    p = zstd_unframe(buf[FILE_HEADER_LEN:], None, "gwas index")
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
    if n_chroms != h["count"]:
        raise ValueError(f"gwas index: {n_chroms} chromosomes, header count {h['count']}")
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
        if nb < 1 or fp[0] < 1 or np.any(np.diff(fp) < 0) or eo[0] <= FILE_HEADER_LEN or np.any(np.diff(eo) <= 0):
            raise ValueError(f"gwas index {name.decode()}: blocks must be non-empty, positions non-decreasing, offsets increasing")
        chroms[name.decode()] = (fp, eo)
    if off != len(p):
        raise ValueError(f"gwas index: {len(p) - off} bytes after the last chromosome")
    return {"block_rows": h["page_size"], "n_values": nv.tolist(), "chroms": chroms}


def gwas_window(first_position, end_offset, lo: int, hi: int) -> tuple[int, int, int, int] | None:
    """SPEC section 11 window rule: (start block, end block, first byte, end byte exclusive) of the blocks that hold every
    row with lo <= position <= hi, or None when no block starts at or before hi. Readers still filter decoded rows."""
    if lo > hi:
        raise ValueError(f"gwas window: lo {lo} > hi {hi}")
    end = int(np.searchsorted(first_position, hi, side="right")) - 1
    if end < 0:
        return None
    start = max(int(np.searchsorted(first_position, lo, side="left")) - 1, 0)
    return start, end, (FILE_HEADER_LEN if start == 0 else int(end_offset[start - 1])), int(end_offset[end])


# ---- reader output (SPEC section 8) ---------------------------------------------------------
def gene_page_rows(block: dict, variants: dict) -> pa.Table:
    """Join a decoded block to decoded pages by vidx: the Arrow table a reader must produce
    (SPEC section 'Reader output', `READER_SCHEMA`). Nulls: rs_number for 0, af for code 65535,
    counts for 65535, pval_nominal for the null code, slope_se for the SE null code, slope where p
    or slope_se is null or p = 0 (t is infinite there), pip and cs_id where the row has no
    credible-set record."""
    n = block["n_rows"]
    if n == 0:
        return READER_SCHEMA.empty_table()
    first = int(variants["vidx"][0])
    i0 = block["var_start"] - first
    if i0 < 0 or i0 + n > len(variants["vidx"]):
        raise ValueError(f"reader: variant pages {first}..{first + len(variants['vidx']) - 1} do not cover rows {block['var_start']}..{block['var_start'] + n - 1}")
    sl = slice(i0, i0 + n)
    pos = variants["position"][sl].astype(np.int64)
    if pos[0] != block["pos_first"] or pos[-1] != block["pos_last"]:
        raise ValueError(f"reader: positions {pos[0]}..{pos[-1]} != block pos_first {block['pos_first']}, pos_last {block['pos_last']}")
    tss = pos - block["anchor"]
    if tss.min() < -2**31 or tss.max() >= 2**31:
        raise ValueError("reader: tss_distance does not fit i32")

    def masked(arr, typ, cast):
        data, mask = np.ma.getdata(arr), np.ma.getmaskarray(arr)
        if np.any(data[~mask] > np.iinfo(cast).max):
            raise ValueError(f"reader: value above {np.iinfo(cast).max} for {typ}")
        return pa.array(np.where(mask, 0, data).astype(cast), mask=mask, type=typ)

    af = variants["af"][sl]
    p = block["pval_nominal"]
    p_null = block["nlp_code"] == NLP_NULL
    se_null = block["se_code"] == SE_NULL
    s_null = p_null | (block["nlp_code"] == NLP_ZERO) | se_null
    pip = np.zeros(n, dtype=np.float32)
    cs_id = np.zeros(n, dtype=np.int8)
    no_cs = np.ones(n, dtype=bool)
    for row, value, cid in zip(block["cs_row"], block["cs_pip"], block["cs_id"]):
        if no_cs[row] or value > pip[row] or (value == pip[row] and cid < cs_id[row]):
            pip[row], cs_id[row], no_cs[row] = value, cid, False
    arrays = [
        pa.array(pos.astype(np.int32), type=pa.int32()),
        pa.array(variants["A1"][sl], type=pa.string()),
        pa.array(variants["A2"][sl], type=pa.string()),
        masked(variants["rs_number"][sl], pa.int64(), np.int64),
        pa.array(tss.astype(np.int32), type=pa.int32()),
        pa.array(np.where(np.isnan(af), 0, af).astype(np.float32), mask=np.isnan(af), type=pa.float32()),
        masked(variants["ma_samples"][sl], pa.int16(), np.int16),
        masked(variants["ma_count"][sl], pa.int16(), np.int16),
        pa.array(np.where(p_null, 0, p), mask=p_null, type=pa.float64()),
        pa.array(np.where(s_null, 0, block["slope"]).astype(np.float32), mask=s_null, type=pa.float32()),
        pa.array(np.where(se_null, 0, block["slope_se"]).astype(np.float32), mask=se_null, type=pa.float32()),
        pa.array(pip, mask=no_cs, type=pa.float32()),
        pa.array(cs_id, mask=no_cs, type=pa.int8()),
    ]
    return pa.Table.from_arrays(arrays, schema=READER_SCHEMA)


# ---- hits pack (SPEC section 13) -------------------------------------------------------------
MAGIC_HITS = b"QVH0"
HITS_HEADER_LEN = 48
HITS_MAX_KIND = 5
HITS_KIND_NAMES = ("trans eQTL", "trans sQTL", "lead eQTL", "lead sQTL", "credible set eQTL", "credible set sQTL")
HITS_FRAME_VARIANTS = 1024           # variants per frame; the file header carries the build's value
MAX_HITS_FRAME_VARIANTS = 65535      # the frame header's n_variants is a u16
PIP_MAXQ = 65535                     # kinds 4-5: v1 = rint(pip * 65535)
SE_MAXQ_LIN = 65535                  # kinds 2-3: v2 = rint(slope_se / se_max * 65535)
SIGNED_MAXQ = 32767                  # kinds 0-1 beta and kinds 2-3 slope: v3 = rint(x / max * 32767)
HIT_FLAG_MINUS, HIT_FLAG_SIG, HIT_FLAG_RESERVED = 0x01, 0x02, 0xFC
_HITS_HEADER = struct.Struct("<4sHHddddd")
assert _HITS_HEADER.size == HITS_HEADER_LEN


def hits_offsets(n_variants: int, n_rows: int) -> dict:
    """Byte offset of every column in a decoded hits frame, and `end`, its payload length (SPEC section 13)."""
    a = HITS_HEADER_LEN + 2 * n_variants
    a += pad4(2 * n_variants)
    b = a + 2 * n_rows
    b += pad4(2 * n_rows)
    return {"count": HITS_HEADER_LEN, "kind": a, "flags": a + n_rows, "gene": b, "v1": b + 2 * n_rows,
            "v2": b + 4 * n_rows, "v3": b + 6 * n_rows, "intron_start": b + 8 * n_rows,
            "intron_end": b + 12 * n_rows, "cluster": b + 16 * n_rows, "end": b + 20 * n_rows}


def _hits_scale(values: np.ndarray, sel: np.ndarray, maxq: int, signed: bool) -> tuple[np.ndarray, float]:
    """Codes for one scaled column and its scale: `scale` is the largest |value| over `sel` (0.0 when
    `sel` is empty), and the code is `rint(value / scale * maxq)` (0 when the scale is 0)."""
    out = np.zeros(values.size, dtype=np.int64)
    if not sel.any():
        return out, 0.0
    v = values[sel]
    if not np.all(np.isfinite(v)):
        raise ValueError("hits frame: a scaled value is null or not finite")
    scale = float(np.abs(v).max())
    if scale > 0:
        out[sel] = np.rint(v / scale * maxq).astype(np.int64)
    if not signed and np.any(out < 0):
        raise ValueError("hits frame: a negative value in an unsigned column")
    return out, scale


def encode_hits_frame(first_vidx, n_variants, vidx, kind, gene_ord, pval=None, beta=None, slope=None, slope_se=None,
                      pip=None, cs_id=None, intron_start=None, intron_end=None, cluster=None, strand=None,
                      significant=None, level: int = 19) -> bytes:
    """One frame of a variant-keyed hits pack as one zstd frame (SPEC section 13).

    The frame covers variant indices `[first_vidx, first_vidx + n_variants)`; every row's `vidx` must
    lie inside it, and a frame with no rows is valid. Rows may come in any order: the encoder sorts
    them by `vidx`, `kind`, then p ascending (kinds 0-3) or PIP descending (kinds 4-5), then `gene`
    (the `search_index` ord), then the intron fields. Per kind, the columns read are: 0-1 (trans)
    `pval` and `beta`; 2-3 (lead) `pval` (the permutation p), `slope_se` and `slope`, with
    `significant` as flags bit 1; 4-5 (credible set) `pip` and `cs_id`. The odd kinds are the sQTL
    ones and need `intron_start`, `intron_end`, `cluster`, and `strand` ('+' or '-', which becomes
    flags bit 0); the even kinds must leave those null or zero. Raises on a value outside its range,
    a null where the kind needs one, an intron field or a `significant` flag on a kind that has none,
    or more than 65535 rows for one variant."""
    n_variants, first_vidx = int(n_variants), int(first_vidx)
    if not 1 <= n_variants <= MAX_HITS_FRAME_VARIANTS:
        raise ValueError(f"hits frame: {n_variants} variants, must be 1..{MAX_HITS_FRAME_VARIANTS}")
    if first_vidx < 0 or first_vidx + n_variants - 1 > U32_MAX:
        raise ValueError(f"hits frame: first_vidx {first_vidx} with {n_variants} variants does not fit u32")
    n = 0 if vidx is None else len(vidx)
    zero = np.zeros(n, dtype=np.int64)
    if n == 0:
        counts = np.zeros(n_variants, dtype="<u2")
        payload = b"".join((_HITS_HEADER.pack(MAGIC_HITS, n_variants, 0, 0.0, 0.0, 0.0, 0.0, 0.0),
                            counts.tobytes(), bytes(pad4(2 * n_variants))))
        return zstd_frame(payload, level)
    vi = _int_column(vidx, "vidx", first_vidx, first_vidx + n_variants - 1, -1, n)
    kd = _int_column(kind, "kind", 0, HITS_MAX_KIND, -1, n)
    go = _int_column(gene_ord, "gene", 0, 65535, -1, n)
    if np.any(vi < 0) or np.any(kd < 0) or np.any(go < 0):
        raise ValueError("hits frame: vidx, kind and gene must be non-null on every row")
    is_trans, is_lead, is_cs = kd <= 1, (kd == 2) | (kd == 3), kd >= 4
    is_s = (kd % 2) == 1
    p = _float_column(pval, "pval", n) if pval is not None else np.full(n, np.nan)
    need_p = is_trans | is_lead
    if np.any(~np.isfinite(p[need_p])) or np.any(p[need_p] <= 0) or np.any(p[need_p] > 1):
        raise ValueError("hits frame: p must be in (0, 1] on kinds 0-3")
    b = _float_column(beta, "beta", n) if beta is not None else np.zeros(n)
    sl = _float_column(slope, "slope", n) if slope is not None else np.zeros(n)
    se = _float_column(slope_se, "slope_se", n) if slope_se is not None else np.zeros(n)
    pp = _float_column(pip, "pip", n) if pip is not None else np.zeros(n)
    if np.any(~np.isfinite(sl[is_lead])) or np.any(~np.isfinite(se[is_lead])) or np.any(se[is_lead] <= 0):
        raise ValueError("hits frame: kinds 2-3 need a finite slope and a finite slope_se above 0")
    if np.any(~np.isfinite(pp[is_cs])) or np.any(pp[is_cs] < 0) or np.any(pp[is_cs] > 1):
        raise ValueError("hits frame: pip must be in [0, 1] on kinds 4-5")
    cs = _int_column(cs_id, "cs_id", 0, 65535, -1, n) if cs_id is not None else np.full(n, -1, dtype=np.int64)
    if np.any(cs[is_cs] < 0):
        raise ValueError("hits frame: cs_id must be 0..65535 on kinds 4-5")
    istart = _int_column(intron_start, "intron_start", 0, U32_MAX, 0, n) if intron_start is not None else zero.copy()
    iend = _int_column(intron_end, "intron_end", 0, U32_MAX, 0, n) if intron_end is not None else zero.copy()
    clu = _int_column(cluster, "cluster", 0, U32_MAX, 0, n) if cluster is not None else zero.copy()
    st = _str_list(strand, "strand", n) if strand is not None else [None] * n
    sig = np.array(_str_list(significant, "significant", n), dtype=bool) if significant is not None else np.zeros(n, dtype=bool)
    minus = np.array([x == "-" for x in st], dtype=bool)
    if np.any(np.array([x not in ("+", "-") for x in st], dtype=bool) & is_s):
        raise ValueError("hits frame: an sQTL row (kinds 1, 3, 5) needs a strand of '+' or '-'")
    if np.any(np.array([x is not None for x in st], dtype=bool) & ~is_s):
        raise ValueError("hits frame: a strand on a kind 0, 2 or 4 row, which has no intron")
    no_intron = ~is_s
    if np.any((istart[no_intron] != 0) | (iend[no_intron] != 0) | (clu[no_intron] != 0)):
        raise ValueError("hits frame: an intron field on a kind 0, 2 or 4 row, which has no intron")
    if np.any(sig & ~is_lead):
        raise ValueError("hits frame: the significant flag is only for lead rows (kinds 2-3)")
    key = np.zeros(n)
    key[need_p] = p[need_p]
    key[is_cs] = -pp[is_cs]
    order = np.lexsort((clu, iend, istart, go, key, kd, vi))
    vi, kd, go, p, b, sl, se, pp, cs, istart, iend, clu, minus, sig = (
        x[order] for x in (vi, kd, go, p, b, sl, se, pp, cs, istart, iend, clu, minus, sig))
    is_trans, is_lead, is_cs, is_s = kd <= 1, (kd == 2) | (kd == 3), kd >= 4, (kd % 2) == 1
    counts = np.bincount(vi - first_vidx, minlength=n_variants)
    if counts.max(initial=0) > 65535:
        raise ValueError(f"hits frame: variant {first_vidx + int(counts.argmax())} has {int(counts.max())} rows, at most 65535")
    v1 = np.zeros(n, dtype=np.int64)
    v2 = np.zeros(n, dtype=np.int64)
    trans_nlp_max = perm_nlp_max = 0.0
    for sel, out in ((is_trans, "trans"), (is_lead, "lead")):
        if sel.any():
            q, m = quantize_nlp(p[sel])
            if int(q.max()) > NLP_MAXQ:
                raise ValueError("hits frame: a p of 0 or null reached the encoder")
            v1[sel] = q.astype(np.int64)
            if out == "trans":
                trans_nlp_max = m
            else:
                perm_nlp_max = m
    v1[is_cs] = np.rint(pp[is_cs] * PIP_MAXQ).astype(np.int64)
    v2[is_cs] = cs[is_cs]
    se_codes, se_max = _hits_scale(se, is_lead, SE_MAXQ_LIN, False)
    v2[is_lead] = se_codes[is_lead]
    v3 = np.zeros(n, dtype=np.int64)
    beta_codes, trans_beta_max = _hits_scale(b, is_trans, SIGNED_MAXQ, True)
    slope_codes, slope_max = _hits_scale(sl, is_lead, SIGNED_MAXQ, True)
    v3[is_trans] = beta_codes[is_trans]
    v3[is_lead] = slope_codes[is_lead]
    flags = (minus * HIT_FLAG_MINUS) | (sig * HIT_FLAG_SIG)
    payload = b"".join((
        _HITS_HEADER.pack(MAGIC_HITS, n_variants, 0, trans_nlp_max, trans_beta_max, perm_nlp_max, se_max, slope_max),
        counts.astype("<u2").tobytes(), bytes(pad4(2 * n_variants)),
        kd.astype("u1").tobytes(), flags.astype("u1").tobytes(), bytes(pad4(2 * n)),
        go.astype("<u2").tobytes(), v1.astype("<u2").tobytes(), v2.astype("<u2").tobytes(), v3.astype("<i2").tobytes(),
        istart.astype("<u4").tobytes(), iend.astype("<u4").tobytes(), clu.astype("<u4").tobytes()))
    return zstd_frame(payload, level)


def decode_hits_frame(buf: bytes, first_vidx: int | None = None, what: str = "hits frame") -> dict:
    """One hits frame under SPEC section 13's reader rules (ValueError names the rule). Returns
    `n_variants`, the five scales, `count` and `start` per variant slot, and per row: `kind`, `flags`,
    `gene` (the `search_index` ord), the raw codes `v1`, `v2`, `v3`, `intron_start`, `intron_end`,
    `cluster`, `strand` ('+'/'-' on kinds 1, 3, 5, else None), `significant`, and the derived
    `pval` (kinds 0-3), `beta` (0-1), `slope` and `slope_se` (2-3), `pip` and `cs_id` (4-5), each NaN
    or None where its kind does not carry it. With `first_vidx` the rows also carry absolute `vidx`."""
    raw = zstd_unframe(bytes(buf), None, what)
    if len(raw) < HITS_HEADER_LEN:
        raise ValueError(f"{what}: {len(raw)} bytes is shorter than the 48-byte header")
    magic, n_variants, reserved, trans_nlp_max, trans_beta_max, perm_nlp_max, se_max, slope_max = _HITS_HEADER.unpack_from(raw, 0)
    if magic != MAGIC_HITS:
        raise ValueError(f"{what}: magic {magic!r} is not {MAGIC_HITS!r}")
    if reserved:
        raise ValueError(f"{what}: reserved header field is not zero")
    if not 1 <= n_variants <= MAX_HITS_FRAME_VARIANTS:
        raise ValueError(f"{what}: n_variants {n_variants} is not 1..{MAX_HITS_FRAME_VARIANTS}")
    scales = {"trans_nlp_max": trans_nlp_max, "trans_beta_max": trans_beta_max, "perm_nlp_max": perm_nlp_max,
              "se_max": se_max, "slope_max": slope_max}
    for name, v in scales.items():
        if not (math.isfinite(v) and v >= 0):
            raise ValueError(f"{what}: scale {name} = {v} must be finite and at least 0")
    o0 = hits_offsets(n_variants, 0)
    if len(raw) < o0["kind"]:
        raise ValueError(f"{what}: {len(raw)} bytes cannot hold {n_variants} counts")
    count = np.frombuffer(raw, "<u2", n_variants, HITS_HEADER_LEN).astype(np.int64)
    n = int(count.sum())
    off = hits_offsets(n_variants, n)
    if len(raw) != off["end"]:
        raise ValueError(f"{what}: decompressed length {len(raw)} != {off['end']} for {n_variants} variants and {n} rows")
    if any(raw[HITS_HEADER_LEN + 2 * n_variants:off["kind"]]) or any(raw[off["flags"] + n:off["gene"]]):
        raise ValueError(f"{what}: padding is not zero")
    kd = np.frombuffer(raw, "u1", n, off["kind"]).astype(np.int64)
    flags = np.frombuffer(raw, "u1", n, off["flags"]).astype(np.int64)
    gene = np.frombuffer(raw, "<u2", n, off["gene"]).astype(np.int64)
    v1 = np.frombuffer(raw, "<u2", n, off["v1"]).astype(np.int64)
    v2 = np.frombuffer(raw, "<u2", n, off["v2"]).astype(np.int64)
    v3 = np.frombuffer(raw, "<i2", n, off["v3"]).astype(np.int64)
    istart = np.frombuffer(raw, "<u4", n, off["intron_start"]).astype(np.int64)
    iend = np.frombuffer(raw, "<u4", n, off["intron_end"]).astype(np.int64)
    clu = np.frombuffer(raw, "<u4", n, off["cluster"]).astype(np.int64)
    if np.any(kd > HITS_MAX_KIND):
        raise ValueError(f"{what}: a kind above {HITS_MAX_KIND}")
    is_trans, is_lead, is_cs, is_s = kd <= 1, (kd == 2) | (kd == 3), kd >= 4, (kd % 2) == 1
    if np.any(flags & HIT_FLAG_RESERVED):
        raise ValueError(f"{what}: reserved flag bits 2-7 are set")
    if np.any(((flags & HIT_FLAG_MINUS) != 0) & ~is_s) or np.any(((flags & HIT_FLAG_SIG) != 0) & ~is_lead):
        raise ValueError(f"{what}: flags bit 0 on a row without an intron, or bit 1 on a row that is not a lead")
    if np.any(v1[is_trans | is_lead] > NLP_MAXQ):
        raise ValueError(f"{what}: an nlp code above {NLP_MAXQ} on kinds 0-3")
    if np.any(v2[is_trans] != 0) or np.any(v3[is_cs] != 0):
        raise ValueError(f"{what}: v2 must be 0 on kinds 0-1 and v3 must be 0 on kinds 4-5")
    if np.any(v3 == -32768):
        raise ValueError(f"{what}: a signed code of -32768")
    if np.any((istart[~is_s] != 0) | (iend[~is_s] != 0) | (clu[~is_s] != 0)):
        raise ValueError(f"{what}: an intron field on a kind 0, 2 or 4 row")
    for scale, sel, codes, maxq in ((trans_nlp_max, is_trans, v1, NLP_MAXQ), (perm_nlp_max, is_lead, v1, NLP_MAXQ),
                                    (se_max, is_lead, v2, SE_MAXQ_LIN), (trans_beta_max, is_trans, np.abs(v3), SIGNED_MAXQ),
                                    (slope_max, is_lead, np.abs(v3), SIGNED_MAXQ)):
        top = int(codes[sel].max(initial=0))
        if (scale > 0 and top != maxq) or (scale == 0 and top != 0):
            raise ValueError(f"{what}: the scale rule (scale {scale} with largest code {top}, expected {maxq if scale > 0 else 0})")
    start = np.r_[0, np.cumsum(count)[:-1]]
    slot = np.repeat(np.arange(n_variants), count)
    same = np.r_[False, (slot[1:] == slot[:-1])] if n else np.zeros(0, dtype=bool)
    if n and (np.any(same & (np.diff(kd, prepend=kd[0]) < 0))
              or np.any(same & (np.diff(kd, prepend=kd[0]) == 0) & (np.diff(v1, prepend=v1[0]) > 0))):
        raise ValueError(f"{what}: rows of one variant must go by kind, then p ascending or PIP descending (v1 never increases)")
    nlp = np.where(is_trans, v1 * (trans_nlp_max / NLP_MAXQ), v1 * (perm_nlp_max / NLP_MAXQ))
    pval = np.where(is_trans | is_lead, np.power(10.0, -nlp), np.nan)
    out = {"n_variants": n_variants, "n_rows": n, **scales, "count": count, "start": start,
           "kind": kd, "flags": flags, "gene": gene, "v1": v1, "v2": v2, "v3": v3,
           "intron_start": istart, "intron_end": iend, "cluster": clu,
           "strand": [("-" if f & HIT_FLAG_MINUS else "+") if s else None for f, s in zip(flags.tolist(), is_s.tolist())],
           "significant": (flags & HIT_FLAG_SIG) != 0,
           "pval": pval, "nlp": np.where(is_trans | is_lead, nlp, np.nan),
           "beta": np.where(is_trans, v3 * (trans_beta_max / SIGNED_MAXQ), np.nan),
           "slope": np.where(is_lead, v3 * (slope_max / SIGNED_MAXQ), np.nan),
           "slope_se": np.where(is_lead, v2 * (se_max / SE_MAXQ_LIN), np.nan),
           "pip": np.where(is_cs, v1 / PIP_MAXQ, np.nan),
           "cs_id": np.where(is_cs, v2, -1), "payload_len": len(raw)}
    if first_vidx is not None:
        out["first_vidx"] = int(first_vidx)
        out["vidx"] = slot + int(first_vidx)
    return out


def hits_slice(frame: dict, vidx: int) -> tuple[int, int]:
    """(start, stop) of one variant's rows in a decoded frame. `vidx` is absolute when the frame was
    decoded with `first_vidx`, otherwise the slot inside the frame."""
    i = int(vidx) - int(frame.get("first_vidx", 0))
    if not 0 <= i < frame["n_variants"]:
        raise ValueError(f"hits frame: variant {vidx} is outside the frame")
    return int(frame["start"][i]), int(frame["start"][i] + frame["count"][i])


# ---- rsID index (SPEC section 14) -------------------------------------------------------------
RSID_CHROM = "all"
RSID_BLOCK_RECORDS = 4096            # records per block; the file header carries the build's value
RSID_RECORD_LEN = 8
RSID_VIDX_BITS = 27
RSID_MAX_VIDX = (1 << RSID_VIDX_BITS) - 1
VARIANT_CHROMS = TRANS_VARIANT_CHROMS        # chr1..chr22, chrX: ordinal = index + 1 (sections 13-15)
RSID_DTYPE = np.dtype([("rs_number", "<u4"), ("ref", "<u4")])
assert RSID_DTYPE.itemsize == RSID_RECORD_LEN


def encode_rsid_index(rs_number, chr_ordinal, vidx, block_records: int = RSID_BLOCK_RECORDS) -> tuple[bytes, np.ndarray]:
    """The whole rsID index: a kind 8 header (chromosome `all`, count = records, page size = records
    per block) and the uncompressed 8-byte records `(rs_number, (chr_ordinal << 27) | vidx)`, which
    must already be sorted by `rs_number`. Returns the bytes and each block's first `rs_number`.
    Raises when `rs_number` does not strictly increase, a chromosome ordinal is outside 1..23, or a
    `vidx` does not fit 27 bits."""
    n = len(rs_number)
    rs = _int_column(rs_number, "rs_number", 1, U32_MAX, -1, n)
    co = _int_column(chr_ordinal, "chr_ordinal", 1, len(VARIANT_CHROMS), -1, n)
    vi = _int_column(vidx, "vidx", 0, RSID_MAX_VIDX, -1, n)
    if np.any(rs < 0) or np.any(co < 0) or np.any(vi < 0):
        raise ValueError("rsid index: rs_number, chr_ordinal and vidx must be non-null")
    if n > 1 and np.any(np.diff(rs) <= 0):
        raise ValueError("rsid index: rs_number must strictly increase")
    if not 1 <= block_records <= MAX_PAGE_RECORDS:
        raise ValueError(f"rsid index: {block_records} records per block, must be 1..{MAX_PAGE_RECORDS}")
    rec = np.empty(n, dtype=RSID_DTYPE)
    rec["rs_number"] = rs
    rec["ref"] = (co << RSID_VIDX_BITS) | vi
    return file_header(KIND_RSID, RSID_CHROM, n, block_records) + rec.tobytes(), rs[::block_records].copy()


def decode_rsid_block(buf: bytes, what: str = "rsid block") -> dict:
    """One block, or any whole number of records, from an rsID index: `rs_number`, `chr_ordinal`,
    `vidx`, and `chrom`. Checks the record length, the ordinals, and that `rs_number` increases."""
    if len(buf) == 0 or len(buf) % RSID_RECORD_LEN:
        raise ValueError(f"{what}: {len(buf)} bytes is not a whole number of 8-byte records")
    rec = np.frombuffer(buf, dtype=RSID_DTYPE)
    rs = rec["rs_number"].astype(np.int64)
    ref = rec["ref"].astype(np.int64)
    co = ref >> RSID_VIDX_BITS
    if np.any((co < 1) | (co > len(VARIANT_CHROMS))):
        raise ValueError(f"{what}: a chromosome ordinal outside 1..{len(VARIANT_CHROMS)}")
    if rs.size > 1 and np.any(np.diff(rs) <= 0):
        raise ValueError(f"{what}: rs_number must strictly increase")
    return {"rs_number": rs, "chr_ordinal": co, "vidx": ref & RSID_MAX_VIDX,
            "chrom": [VARIANT_CHROMS[c - 1] for c in co.tolist()]}


def rsid_block_of(rsid_first: np.ndarray, rs_number: int) -> int | None:
    """The block that can hold `rs_number`: the last one whose first record is at or below it, or
    None when the value is below every block's first record."""
    i = int(np.searchsorted(np.asarray(rsid_first), int(rs_number), side="right")) - 1
    return None if i < 0 else i


def rsid_block_range(block: int, n_records: int, block_records: int = RSID_BLOCK_RECORDS) -> tuple[int, int]:
    """(byte offset, length) of one block of an rsID index file."""
    start = block * block_records
    if block < 0 or start >= n_records:
        raise ValueError(f"rsid index: block {block} is outside a file of {n_records} records")
    return FILE_HEADER_LEN + start * RSID_RECORD_LEN, min(block_records, n_records - start) * RSID_RECORD_LEN


def rsid_find(block: dict, rs_number: int) -> tuple[str, int] | None:
    """(chromosome, vidx) of `rs_number` in a decoded block, or None when the block does not hold it."""
    i = int(np.searchsorted(block["rs_number"], int(rs_number)))
    if i >= block["rs_number"].size or int(block["rs_number"][i]) != int(rs_number):
        return None
    return block["chrom"][i], int(block["vidx"][i])


# ---- variant index (SPEC section 15) ----------------------------------------------------------
MAGIC_VARIANT_INDEX = b"QVX0"
VARIANT_INDEX_HEADER_LEN = 32
_VX_HEADER = struct.Struct("<4sHHIIIIII")
assert _VX_HEADER.size == VARIANT_INDEX_HEADER_LEN


def encode_variant_index(chroms: list[dict], page_size: int, frame_variants: int, rsid_first,
                         rsid_n_records: int, rsid_block_records: int = RSID_BLOCK_RECORDS, level: int = 19) -> bytes:
    """The startup file: a kind 9 header (chromosome `all`, count = chromosomes, page size = variants
    per page) and one zstd frame holding, for chr1..chr22 then chrX, each variants file's page offsets
    and first positions and each hits file's frame offsets, then every rsID block's first `rs_number`
    (SPEC section 15).

    `chroms` is one dict per chromosome, in that order, with `n_cis`, `n_trans_only`, `page_off`
    (`n_pages + 1` byte offsets, the last the variants file's size), `page_first_position` (`n_pages`),
    and `hits_off` (`n_frames + 1`, the last the hits file's size). Raises when a page or frame count
    does not follow from the variant counts, an offset run does not increase from byte 32, or a first
    position is 0."""
    if len(chroms) != len(VARIANT_CHROMS):
        raise ValueError(f"variant index: {len(chroms)} chromosomes, expected {len(VARIANT_CHROMS)}")
    if not 1 <= page_size <= MAX_PAGE_RECORDS or not 1 <= frame_variants <= MAX_HITS_FRAME_VARIANTS:
        raise ValueError("variant index: page size and frame variants must be 1..65535")
    rf = _int_column(rsid_first, "rsid_first", 1, U32_MAX, -1, len(rsid_first))
    n_blocks = -(-int(rsid_n_records) // int(rsid_block_records))
    if rf.size != n_blocks or (rf.size > 1 and np.any(np.diff(rf) <= 0)):
        raise ValueError(f"variant index: {rf.size} rsID block samples for {n_blocks} blocks, or they do not increase")
    parts = [_VX_HEADER.pack(MAGIC_VARIANT_INDEX, len(chroms), 0, page_size, frame_variants,
                             rsid_block_records, int(rsid_n_records), n_blocks, 0)]
    for name, c in zip(VARIANT_CHROMS, chroms):
        n_cis, n_tr = int(c["n_cis"]), int(c["n_trans_only"])
        pc, pt = -(-n_cis // page_size), -(-n_tr // page_size)
        nf = -(-(n_cis + n_tr) // frame_variants)
        po = _int_column(c["page_off"], f"{name} page_off", FILE_HEADER_LEN, U32_MAX, -1, len(c["page_off"]))
        pf = _int_column(c["page_first_position"], f"{name} page_first_position", 1, U32_MAX, -1, len(c["page_first_position"]))
        ho = _int_column(c["hits_off"], f"{name} hits_off", FILE_HEADER_LEN, U32_MAX, -1, len(c["hits_off"]))
        if po.size != pc + pt + 1 or pf.size != pc + pt or ho.size != nf + 1:
            raise ValueError(f"variant index {name}: {po.size - 1} pages and {ho.size - 1} frames do not follow "
                             f"{n_cis} cis + {n_tr} trans-only variants ({pc + pt} pages, {nf} frames)")
        if np.any(np.diff(po) <= 0) or np.any(np.diff(ho) <= 0) or po[0] != FILE_HEADER_LEN or ho[0] != FILE_HEADER_LEN:
            raise ValueError(f"variant index {name}: page and frame offsets must increase from byte {FILE_HEADER_LEN}")
        parts += [struct.pack("<IIIII", n_cis, n_tr, pc, pt, nf), po.astype("<u4").tobytes(),
                  pf.astype("<u4").tobytes(), ho.astype("<u4").tobytes()]
    parts.append(rf.astype("<u4").tobytes())
    return file_header(KIND_VARIANT_INDEX, RSID_CHROM, len(chroms), page_size) + zstd_frame(b"".join(parts), level)


def decode_variant_index(buf: bytes, what: str = "variant index") -> dict:
    """The startup file under SPEC section 15's reader rules: `page_size`, `frame_variants`,
    `rsid_block_records`, `rsid_n_records`, `rsid_first`, and `chroms`, a dict of chromosome name to
    `n_cis`, `n_trans_only`, `n_pages_cis`, `n_pages_trans`, `n_frames`, `page_off`,
    `page_first_position`, and `hits_off`."""
    h = parse_file_header(buf)
    if h["kind"] != KIND_VARIANT_INDEX or h["chrom"] != RSID_CHROM:
        raise ValueError(f"{what}: header kind {h['kind']} chromosome {h['chrom']!r}, expected {KIND_VARIANT_INDEX} and {RSID_CHROM!r}")
    p = zstd_unframe(buf[FILE_HEADER_LEN:], None, what)
    if len(p) < VARIANT_INDEX_HEADER_LEN:
        raise ValueError(f"{what}: payload is shorter than the 32-byte header")
    magic, n_chrom, reserved, page_size, frame_variants, rsid_block_records, rsid_n, rsid_blocks, reserved2 = _VX_HEADER.unpack_from(p, 0)
    if magic != MAGIC_VARIANT_INDEX:
        raise ValueError(f"{what}: magic {magic!r} is not {MAGIC_VARIANT_INDEX!r}")
    if reserved or reserved2:
        raise ValueError(f"{what}: a reserved header field is not zero")
    if n_chrom != h["count"] or n_chrom != len(VARIANT_CHROMS):
        raise ValueError(f"{what}: {n_chrom} chromosomes, header count {h['count']}, expected {len(VARIANT_CHROMS)}")
    if page_size != h["page_size"] or not 1 <= page_size <= MAX_PAGE_RECORDS:
        raise ValueError(f"{what}: page size {page_size} differs from the file header's {h['page_size']} or is out of range")
    if not 1 <= frame_variants <= MAX_HITS_FRAME_VARIANTS or not 1 <= rsid_block_records <= MAX_PAGE_RECORDS:
        raise ValueError(f"{what}: frame variants {frame_variants} or rsID block records {rsid_block_records} out of range")
    if rsid_blocks != -(-rsid_n // rsid_block_records):
        raise ValueError(f"{what}: {rsid_blocks} rsID blocks for {rsid_n} records of {rsid_block_records}")
    off = VARIANT_INDEX_HEADER_LEN
    chroms = {}

    def u32s(k: int) -> np.ndarray:
        nonlocal off
        if off + 4 * k > len(p):
            raise ValueError(f"{what}: payload ends early")
        a = np.frombuffer(p, "<u4", k, off).astype(np.int64)
        off += 4 * k
        return a

    for name in VARIANT_CHROMS:
        n_cis, n_tr, pc, pt, nf = (int(x) for x in u32s(5))
        if pc != -(-n_cis // page_size) or pt != -(-n_tr // page_size) or nf != -(-(n_cis + n_tr) // frame_variants):
            raise ValueError(f"{what} {name}: {pc} + {pt} pages and {nf} frames do not follow {n_cis} cis + {n_tr} trans-only variants")
        po, pf, ho = u32s(pc + pt + 1), u32s(pc + pt), u32s(nf + 1)
        if po.size and (po[0] != FILE_HEADER_LEN or np.any(np.diff(po) <= 0)):
            raise ValueError(f"{what} {name}: page offsets must increase from byte {FILE_HEADER_LEN}")
        if ho.size and (ho[0] != FILE_HEADER_LEN or np.any(np.diff(ho) <= 0)):
            raise ValueError(f"{what} {name}: hits frame offsets must increase from byte {FILE_HEADER_LEN}")
        if np.any(pf < 1) or np.any(np.diff(pf[:pc]) < 0) or np.any(np.diff(pf[pc:]) < 0):
            raise ValueError(f"{what} {name}: page first positions must be at least 1 and never decrease within a section")
        chroms[name] = {"n_cis": n_cis, "n_trans_only": n_tr, "n_pages_cis": pc, "n_pages_trans": pt, "n_frames": nf,
                        "page_off": po, "page_first_position": pf, "hits_off": ho}
    rf = u32s(rsid_blocks)
    if rf.size > 1 and np.any(np.diff(rf) <= 0):
        raise ValueError(f"{what}: rsID block samples must strictly increase")
    if off != len(p):
        raise ValueError(f"{what}: {len(p) - off} bytes after the rsID block samples")
    return {"page_size": page_size, "frame_variants": frame_variants, "rsid_block_records": rsid_block_records,
            "rsid_n_records": rsid_n, "rsid_n_blocks": rsid_blocks, "rsid_first": rf, "chroms": chroms}


def variant_index_page(index: dict, chrom: str, vidx: int) -> tuple[int, int, int]:
    """(page, byte offset, length) of the variants-file page holding `vidx` (SPEC section 15)."""
    c = index["chroms"][chrom]
    P, n_cis = index["page_size"], c["n_cis"]
    if not 0 <= vidx < n_cis + c["n_trans_only"]:
        raise ValueError(f"variant index {chrom}: vidx {vidx} is outside {n_cis + c['n_trans_only']} variants")
    page = vidx // P if vidx < n_cis else c["n_pages_cis"] + (vidx - n_cis) // P
    return page, int(c["page_off"][page]), int(c["page_off"][page + 1] - c["page_off"][page])


def variant_index_position(index: dict, chrom: str, position: int, section: str = "cis") -> tuple[int, int, int]:
    """(page, byte offset, length) of the page that can hold `position` in one section, or a
    ValueError when the section has no page at or below it. The caller still checks the decoded
    records: a position the section does not hold lands in the page before it."""
    c = index["chroms"][chrom]
    lo, hi = (0, c["n_pages_cis"]) if section == "cis" else (c["n_pages_cis"], c["n_pages_cis"] + c["n_pages_trans"])
    first = c["page_first_position"][lo:hi]
    i = int(np.searchsorted(first, int(position), side="right")) - 1
    if i < 0:
        raise ValueError(f"variant index {chrom}: position {position} is below the {section} section's first page")
    page = lo + i
    return page, int(c["page_off"][page]), int(c["page_off"][page + 1] - c["page_off"][page])


def variant_index_hits(index: dict, chrom: str, vidx: int) -> tuple[int, int, int]:
    """(frame, byte offset, length) of the hits frame holding `vidx` (SPEC section 15)."""
    c = index["chroms"][chrom]
    if not 0 <= vidx < c["n_cis"] + c["n_trans_only"]:
        raise ValueError(f"variant index {chrom}: vidx {vidx} is outside {c['n_cis'] + c['n_trans_only']} variants")
    frame = vidx // index["frame_variants"]
    return frame, int(c["hits_off"][frame]), int(c["hits_off"][frame + 1] - c["hits_off"][frame])


def variants_page_firsts(buf: bytes) -> tuple[dict, list[dict], np.ndarray]:
    """Header, every page (`walk_variants_file`), and each page's first position: the position column
    starts 4 bytes into the payload and its first entry is absolute (section 4). What the
    `variant_index` build reads; it decompresses a zstd page but decodes nothing else."""
    h, pages = walk_variants_file(buf)
    first = np.zeros(len(pages), dtype=np.int64)
    for i, pg in enumerate(pages):
        stored = _PAGE_HEADER.unpack_from(buf, pg["offset"])[0]
        body = buf[pg["offset"] + PAGE_HEADER_LEN:pg["offset"] + PAGE_HEADER_LEN + stored]
        if pg["codec"] == "zstd":
            body = zstd_unframe(body, None, f"variants file, page at byte {pg['offset']}")
        if len(body) < 8:
            raise ValueError(f"variants file: page at byte {pg['offset']} has no position column")
        first[i] = int(np.frombuffer(body, "<u4", 1, 4)[0])
    return h, pages, first
