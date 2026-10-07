"""Measure a unified p/beta/se encoding against the current one, on real TOPCHeF numbers.

    uv run python -m pipeline.schema_probe cis  --kind eqtl --chrom chr21
    uv run python -m pipeline.schema_probe gwas --chrom chr21
    uv run python -m pipeline.schema_probe cis  --kind sqtl --chrom chr22 --json out.json

`plans/2026-10-07-unified-schema.md` is the proposal. This module only measures it: it writes no
store objects, touches neither UI, and never speaks to the bucket. It reads the local source tables
(`data/derived/cis_{e,s}qtl_nominal`, `data/derived/gwas_dcm`) -- the same numbers a v1 build
ingests -- encodes them both ways, decodes both, and compares each against the source.

What "current" means per object kind, from `packfmt_v1`:

- cis: `nlp` u16 against the block's `nlp_max`, `se` as a 15-bit log code with the slope's sign in
  bit 15, and **beta not stored** -- rebuilt as `sign * se * t(nlp, dof)`.
- GWAS: lossless. `beta` `i32` = `rint(b * 1e4)`, `se`/`af` `u16` likewise, `p` as a `u16` mantissa
  in 1000..9999 with an `i8` exponent. 21 bytes a row, and `gwas_codes` refuses a source carrying
  more precision than that holds.

What "unified" means, the same for both: `u16` nlp / `i16` beta / `u16` log-coded se, 6 bytes, four
f64 scales per scope. The sign rides in beta's `i16`, which frees bit 15 of the SE field, so se gets
16 bits rather than 15.

Three things are reported, because they are three different questions:

1. **Bytes a row**, and what that projects to over the real row counts.
2. **Accuracy against the source**, per field, in the format's own currency: beta error in units of
   the row's own SE (what SPEC section 9's 4.47e-03 budget is denominated in), se as a relative
   error, nlp as an absolute error in -log10 p.
3. **What each scheme cannot represent at all.** The current scheme loses beta wherever the rebuild
   cannot run -- `p = 0`, a null dof, or `-log10 p` over `NLP_LIMIT`, which fails the build outright.
   The unified scheme loses beta only where the source had none. This is the comparison the whole
   proposal rests on, so it is counted rather than argued.

It also measures the risk the proposal introduces (plan decision D4): with all three stored
independently, `beta/se` and `t(nlp, dof)` can disagree by more than quantization noise, and no
stored value is authoritative. `consistency` reports the worst observed disagreement.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from . import packfmt_v1 as pf
from .common import Config, log

# unified codes. nlp reuses the existing quantizer untouched; se drops the sign and takes bit 15.
BETA_MAXQ = 32766          # i16, -32768 reserved for null
SE16_MAXQ = 65534          # u16, 65535 reserved for null
UNIFIED_ROW_BYTES = 6
CIS_ROW_BYTES_NOW = 4
GWAS_ROW_BYTES_NOW = pf.GWAS_ROW_BYTES          # 21
GWAS_IDENTITY_BYTES = 12                        # pos delta, rs, af, n_code, allele
GWAS_ROW_BYTES_UNIFIED = GWAS_IDENTITY_BYTES + UNIFIED_ROW_BYTES

# real row counts, for projecting a store total from a sampled chromosome
ROWS = {"eqtl": 123_458_053, "sqtl": 499_596_989, "gwas": 12_504_079}
OBJECT_MB = {"eqtl_sqtl": 3284.2, "gwas": 167.4}


# ---- unified codec -----------------------------------------------------------------------------
def quantize_beta(beta: np.ndarray) -> tuple[np.ndarray, float]:
    """`i16` codes and `beta_max` for one scope. Null or non-finite -> -32768; a linear ruler,
    because beta crosses zero and a log one is undefined there."""
    b = np.asarray(beta, dtype=np.float64)
    null = ~np.isfinite(b)
    beta_max = float(np.abs(b[~null]).max()) if np.any(~null) else 0.0
    q = np.zeros(b.shape, dtype=np.int16)
    if beta_max > 0:
        q[~null] = np.rint(b[~null] / beta_max * BETA_MAXQ).astype(np.int16)
    q[null] = -32768
    return q, beta_max


def dequantize_beta(q: np.ndarray, beta_max: float) -> np.ndarray:
    q = np.asarray(q, dtype=np.int16)
    out = q.astype(np.float64) * (float(beta_max) / BETA_MAXQ)
    out[q == -32768] = np.nan
    return out


def quantize_se16(se: np.ndarray) -> tuple[np.ndarray, float, float]:
    """`u16` log codes and (lse_min, lse_max). No sign bit: beta carries the sign now, so all 16
    bits are magnitude and the error is half what `packfmt_v1.quantize_se` gives."""
    s = np.asarray(se, dtype=np.float64)
    null = ~np.isfinite(s) | (s <= 0)
    if np.any(~null):
        lse = np.log(s[~null])
        lse_min, lse_max = float(lse.min()), float(lse.max())
    else:
        lse_min = lse_max = 0.0
    q = np.full(s.shape, 65535, dtype=np.uint16)
    if np.any(~null):
        span = lse_max - lse_min
        if span == 0:
            q[~null] = 0
        else:
            q[~null] = np.rint((np.log(s[~null]) - lse_min) / span * SE16_MAXQ).astype(np.uint16)
    return q, lse_min, lse_max


def dequantize_se16(q: np.ndarray, lse_min: float, lse_max: float) -> np.ndarray:
    q = np.asarray(q, dtype=np.uint16)
    step = (float(lse_max) - float(lse_min)) / SE16_MAXQ
    out = np.exp(float(lse_min) + q.astype(np.float64) * step)
    out[q == 65535] = np.nan
    return out


# ---- comparison --------------------------------------------------------------------------------
def _pct(x: np.ndarray, qs=(50, 95, 99, 100)) -> dict:
    x = x[np.isfinite(x)]
    if not x.size:
        return {}
    return {f"p{q}": float(np.percentile(x, q)) for q in qs}


def compare_scope(p: np.ndarray, beta: np.ndarray, se: np.ndarray, dof: int | None) -> dict:
    """One scope's rows under both schemes, each against the source. Errors only over rows the
    scheme can actually represent; the rows it cannot are counted separately."""
    p = np.asarray(p, dtype=np.float64)
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    n = len(p)
    out: dict = {"rows": n}

    # ---- current: nlp + se(+sign), beta rebuilt through the inverse t
    over_limit = False
    try:
        nlp_q, nlp_max = pf.quantize_nlp(p)
    except ValueError as e:                     # NLP_LIMIT: the build would fail here
        over_limit = True
        out["current"] = {"build_fails": str(e)}
        nlp_q, nlp_max = None, None
    if not over_limit:
        se_q, lse_min, lse_max = pf.quantize_se(se, beta)
        nlp_d = pf.dequantize_nlp(nlp_q, nlp_max)
        p_d = pf.p_from_nlp(nlp_d)
        se_d, neg = pf.dequantize_se(se_q, lse_min, lse_max)
        beta_d = pf.slope_from_se(se_d, neg, p_d, dof) if dof else np.full(n, np.nan)
        ok = np.isfinite(beta_d) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
        out["current"] = {
            "bytes_per_row": CIS_ROW_BYTES_NOW,
            "nlp_max": nlp_max, "beta_stored": False,
            "beta_err_in_se": _pct(np.abs(beta_d[ok] - beta[ok]) / se[ok]),
            "se_rel_err": _pct(np.abs(se_d - se) / se),
            "nlp_abs_err": _pct(np.abs(nlp_d - (-np.log10(np.where(p > 0, p, np.nan))))),
            "beta_unrepresentable": int(n - ok.sum()),
            "why": {"p_zero": int((p == 0).sum()), "no_dof": 0 if dof else n,
                    "p_null": int(np.isnan(p).sum())},
        }

    # ---- unified: nlp + beta + se, nothing rebuilt
    try:
        u_nlp_q, u_nlp_max = pf.quantize_nlp(p)
        u_over = False
    except ValueError:
        # the unified scheme has no inverse t, so the 300 limit does not apply to it; quantize the
        # same way without the ceiling to show what it would hold
        nlp = np.where(p > 0, -np.log10(np.where(p > 0, p, 1.0)), np.nan)
        u_nlp_max = float(np.nanmax(nlp))
        u_nlp_q = np.zeros(n, dtype=np.uint16)
        fin = np.isfinite(nlp)
        u_nlp_q[fin] = np.rint(nlp[fin] / u_nlp_max * pf.NLP_MAXQ).astype(np.uint16)
        u_nlp_q[p == 0] = pf.NLP_ZERO
        u_nlp_q[np.isnan(p)] = pf.NLP_NULL
        u_over = True
    b_q, beta_max = quantize_beta(beta)
    s_q, u_lse_min, u_lse_max = quantize_se16(se)
    b_d = dequantize_beta(b_q, beta_max)
    s_d = dequantize_se16(s_q, u_lse_min, u_lse_max)
    nlp_u = u_nlp_q.astype(np.float64) * (u_nlp_max / pf.NLP_MAXQ)
    nlp_u[u_nlp_q == pf.NLP_ZERO] = np.inf
    nlp_u[u_nlp_q == pf.NLP_NULL] = np.nan
    ok_b = np.isfinite(b_d) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
    out["unified"] = {
        "bytes_per_row": UNIFIED_ROW_BYTES,
        "nlp_max": u_nlp_max, "beta_max": beta_max, "beta_stored": True,
        "exceeds_current_nlp_limit": bool(u_over),
        "beta_err_in_se": _pct(np.abs(b_d[ok_b] - beta[ok_b]) / se[ok_b]),
        "se_rel_err": _pct(np.abs(s_d - se) / se),
        "nlp_abs_err": _pct(np.abs(nlp_u - (-np.log10(np.where(p > 0, p, np.nan))))),
        "beta_unrepresentable": int(n - ok_b.sum()),
        "why": {"source_beta_null": int((~np.isfinite(beta)).sum())},
    }

    # ---- the risk the proposal adds (plan D4): do the three stored values agree?
    if dof:
        ok_c = np.isfinite(b_d) & np.isfinite(s_d) & np.isfinite(nlp_u) & (s_d > 0)
        if np.any(ok_c):
            t_stored = np.abs(b_d[ok_c]) / s_d[ok_c]
            t_from_p = pf.t_from_p(pf.p_from_nlp(nlp_u[ok_c]), dof)
            out["consistency"] = {
                "rows": int(ok_c.sum()),
                "abs_t_disagreement": _pct(np.abs(t_stored - t_from_p)),
                "note": "|beta|/se against t(nlp, dof); both are stored now, neither is authoritative",
            }
    return out


# ---- sources -----------------------------------------------------------------------------------
def cis_scopes(cfg: Config, kind: str, chrom: str, limit: int | None, tables: Path | None = None):
    """(phenotype id, p, beta, se) per phenotype for one chromosome. A phenotype is the cis scope,
    so it is also the quantization scope.

    Two input layouts, because the two things worth measuring live in different places. TOPCHeF's v1
    contract tables are on Rivanna, so its numbers come from the local v0-era tables
    (`gene_id`/`phenotype_id`, `pval_nominal`, `slope`, `slope_se`) -- that shows what unification
    *costs* on a study the current scheme handles perfectly. `--tables` reads a contract tree
    instead (`phenotype_id`, `pvalue`, `beta`, `se`), which is how ARIC gets measured -- and ARIC is
    the study that shows what unification *buys*, since it is the one with `p = 0` rows and blocks
    over `NLP_LIMIT`.
    """
    if tables is not None:
        files = sorted(glob.glob(str(Path(tables) / "nominal" / f"chr={chrom}" / "*.parquet")))
        if not files:
            raise SystemExit(f"no contract nominal files for {chrom} under {tables}")
        return _scopes(files, "phenotype_id", ("pvalue", "beta", "se"), limit)
    key = {"eqtl": "cis_eqtl_nominal", "sqtl": "cis_sqtl_nominal"}[kind]
    files = sorted(glob.glob(str(cfg.derived / key / f"chr={chrom}" / "bin=*" / "*.parquet")))
    if not files:
        raise SystemExit(f"no {key} files for {chrom} under {cfg.derived}")
    idc = "gene_id" if kind == "eqtl" else "phenotype_id"
    return _scopes(files, idc, ("pval_nominal", "slope", "slope_se"), limit)


def _scopes(files: list[str], idc: str, cols: tuple[str, str, str], limit: int | None):
    seen = 0
    pc, bc, sc = cols
    for f in files:
        t = pq.read_table(f, columns=[idc, pc, bc, sc])
        ids = np.array(t[idc].to_pylist())
        p = t[pc].to_numpy(zero_copy_only=False).astype(np.float64)
        b = t[bc].to_numpy(zero_copy_only=False).astype(np.float64)
        s = t[sc].to_numpy(zero_copy_only=False).astype(np.float64)
        o = np.argsort(ids, kind="stable")
        ids, p, b, s = ids[o], p[o], b[o], s[o]
        cuts = np.r_[0, np.flatnonzero(ids[1:] != ids[:-1]) + 1, len(ids)]
        for i in range(len(cuts) - 1):
            a, z = cuts[i], cuts[i + 1]
            yield ids[a], p[a:z], b[a:z], s[a:z]
            seen += 1
            if limit and seen >= limit:
                return


def gwas_scopes(cfg: Config, chrom: str, rows_per_block: int, limit: int | None):
    """(label, p, beta, se) per block of one chromosome's GWAS rows. The GWAS block is the scope.

    `beta` and `se` are rounded back to 4 decimals on the way in. The local table stores them as
    **float32**, so Jurgens' published `-0.0827` reads back as `-0.08269999921321869`, which
    `gwas_codes` rightly refuses as carrying more than 4 decimals -- its whole job is to catch a
    source it cannot hold losslessly. Rounding restores the number the source actually printed,
    which is both what the real adapter encodes and the only honest reference to score against: an
    error measured against float32 noise is not an error in the codec.

    This is not the same as the cis path, where float32 **is** the source's precision -- tensorQTL
    writes float32, which is why `packfmt_v1.REL_F32` exists.
    """
    f = cfg.derived / "gwas_dcm" / f"chr={chrom}" / "data.parquet"
    if not f.exists():
        raise SystemExit(f"no GWAS table at {f}")
    t = pq.read_table(f, columns=["position", "beta", "se", "p"])
    p = t["p"].to_numpy(zero_copy_only=False).astype(np.float64)
    b = np.round(t["beta"].to_numpy(zero_copy_only=False).astype(np.float64), 4)
    s = np.round(t["se"].to_numpy(zero_copy_only=False).astype(np.float64), 4)
    for k, a in enumerate(range(0, len(p), rows_per_block)):
        z = min(a + rows_per_block, len(p))
        yield f"{chrom}:block{k}", p[a:z], b[a:z], s[a:z]
        if limit and k + 1 >= limit:
            return


def gwas_current(p: np.ndarray, beta: np.ndarray, se: np.ndarray) -> dict:
    """The GWAS side's current lossless codes, and whether the source fits them. `gwas_codes`
    *raises* on a source carrying more precision than `rint(x * 1e4)` holds, so a failure here is
    the build-time guard that unifying would remove."""
    out: dict = {"bytes_per_row": GWAS_ROW_BYTES_NOW, "lossless": True}
    try:
        b_q = pf.gwas_scaled(beta, "beta", -(2**31 - 1), 2**31 - 1)
        s_q = pf.gwas_scaled(se, "se", 0, 65535)
        mant, exp = pf.gwas_p_codes(p)
    except ValueError as e:
        out |= {"lossless": False, "build_fails": str(e)[:160]}
        return out
    b_d = b_q.astype(np.float64) / pf.GWAS_SCALE
    s_d = s_q.astype(np.float64) / pf.GWAS_SCALE
    p_d = mant.astype(np.float64) * np.power(10.0, exp.astype(np.float64))
    ok = np.isfinite(s_d) & (s_d > 0)
    out |= {
        "beta_err_in_se": _pct(np.abs(b_d[ok] - beta[ok]) / se[ok]),
        "se_rel_err": _pct(np.abs(s_d - se) / se),
        "nlp_abs_err": _pct(np.abs(-np.log10(p_d) - (-np.log10(np.where(p > 0, p, np.nan))))),
        "beta_unrepresentable": 0,
    }
    return out


# ---- roll-up -----------------------------------------------------------------------------------
def roll(scopes: list[dict]) -> dict:
    """Worst case and median-of-medians across scopes, plus the counts summed."""
    def gather(path):
        vals = []
        for s in scopes:
            d = s
            for k in path:
                d = (d or {}).get(k) if isinstance(d, dict) else None
            if isinstance(d, (int, float)):
                vals.append(float(d))
        return vals
    out: dict = {"scopes": len(scopes), "rows": int(sum(s.get("rows", 0) for s in scopes))}
    for scheme in ("current", "unified"):
        o: dict = {}
        for field in ("beta_err_in_se", "se_rel_err", "nlp_abs_err"):
            med = gather([scheme, field, "p50"])
            worst = gather([scheme, field, "p100"])
            if med:
                o[field] = {"median_of_scopes": float(np.median(med)), "worst_scope": float(max(worst or med))}
        o["beta_unrepresentable"] = int(sum(int((s.get(scheme) or {}).get("beta_unrepresentable", 0)) for s in scopes))
        o["scopes_whose_build_fails"] = sum(1 for s in scopes if "build_fails" in (s.get(scheme) or {}))
        o["bytes_per_row"] = next((int((s.get(scheme) or {}).get("bytes_per_row", 0)) for s in scopes
                                   if (s.get(scheme) or {}).get("bytes_per_row")), None)
        out[scheme] = o
    dis = gather(["consistency", "abs_t_disagreement", "p100"])
    if dis:
        out["consistency"] = {"worst_abs_t_disagreement": float(max(dis)),
                              "median_scope_worst": float(np.median(dis))}
    return out


def project(kind: str, cur_bytes: int, uni_bytes: int) -> dict:
    """Store-level projection from the per-row widths and the known row counts."""
    if kind == "gwas":
        rows, now_mb = ROWS["gwas"], OBJECT_MB["gwas"]
        raw_now = rows * GWAS_ROW_BYTES_NOW
        raw_uni = rows * GWAS_ROW_BYTES_UNIFIED
        return {"rows": rows, "objects_now_mb": now_mb,
                "raw_now_mb": raw_now / 1e6, "raw_unified_mb": raw_uni / 1e6,
                "objects_unified_mb_est": now_mb * raw_uni / raw_now,
                "note": "GWAS blocks are zstd framed, so the estimate scales the measured object size"}
    rows = ROWS["eqtl"] + ROWS["sqtl"]
    return {"rows": rows, "objects_now_mb": OBJECT_MB["eqtl_sqtl"],
            "pairs_now_mb": rows * cur_bytes / 1e6, "pairs_unified_mb": rows * uni_bytes / 1e6,
            "added_mb": rows * (uni_bytes - cur_bytes) / 1e6,
            "note": "cis blocks store pairs uncompressed, so the delta is exact"}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("what", choices=("cis", "gwas"))
    ap.add_argument("--kind", choices=("eqtl", "sqtl"), default="eqtl")
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--limit", type=int, help="stop after this many scopes")
    ap.add_argument("--rows-per-block", type=int, default=4096, help="GWAS block size")
    ap.add_argument("--tables", type=Path, help="a contract tables tree (ARIC), instead of the v0 tables")
    ap.add_argument("--dof", type=int, help="override the dof used to rebuild beta under the current scheme")
    ap.add_argument("--json", type=Path)
    args = ap.parse_args(argv)
    cfg = Config()

    scopes: list[dict] = []
    if args.what == "cis":
        dof = args.dof or cfg["packs"]["dof"][args.kind]
        src = args.tables or f"v0 {args.kind}"
        log(f"schema_probe cis {src} {args.chrom}: dof {dof}")
        for _id, p, b, s in cis_scopes(cfg, args.kind, args.chrom, args.limit, args.tables):
            scopes.append(compare_scope(p, b, s, dof))
    else:
        log(f"schema_probe gwas {args.chrom}: blocks of {args.rows_per_block} rows")
        for label, p, b, s in gwas_scopes(cfg, args.chrom, args.rows_per_block, args.limit):
            d = compare_scope(p, b, s, None)
            d["current"] = gwas_current(p, b, s)        # the real GWAS codec, not the cis one
            scopes.append(d)

    out = roll(scopes)
    cur = out["current"]["bytes_per_row"] or (GWAS_ROW_BYTES_NOW if args.what == "gwas" else CIS_ROW_BYTES_NOW)
    uni = GWAS_ROW_BYTES_UNIFIED if args.what == "gwas" else UNIFIED_ROW_BYTES
    out["bytes_per_row"] = {"current": cur, "unified": uni}
    out["projection"] = project("gwas" if args.what == "gwas" else args.kind, cur, uni)
    out["what"] = f"{args.what}:{args.kind if args.what == 'cis' else 'dcm'}:{args.chrom}"
    print(json.dumps(out, indent=1, default=str))
    if args.json:
        args.json.write_text(json.dumps({"summary": out, "scopes": scopes}, indent=1, default=str) + "\n")
        log(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
