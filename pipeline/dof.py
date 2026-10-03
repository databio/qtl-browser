"""Fit the Student-t degrees of freedom that ties `beta`, `se` and `pvalue` together.

    uv run python -m pipeline.dof 'data/derived/_tables/cis_eqtl_nominal/chr=*/bin=*/data.parquet' \
        --n-samples 516 --beta slope --se slope_se --p pval_nominal --json out.json

SPEC section 5 stores `-log10 p` and `slope_se` and derives `slope = sign * slope_se * t`, with
`t = -stdtrit(dof, p / 2)` (decision D13). That trade is only sound when `dof` is right: one wrong
dof silently bends every slope in the experiment. TOPCHeF publishes its dof (435 eQTL, 480 sQTL).
The eQTL Catalogue does not, so it has to come out of the data.

The fit is a grid search over integer dof minimising the median `|log10 p ratio|` between the
published p and `2 * stdtr(dof, -|beta/se|)`. It works because a correct dof reproduces the
published p to the float32 precision of the source itself (packcheck check 3: 0 of 122,561,970
eQTL and 0 of 495,708,407 sQTL rows above 1e-5 relative error, max 5.96e-08), so the objective at
the right dof sits at the source's noise floor while the nearest wrong dof sits orders of magnitude
above it. `fit()` reports that separation as `margin` rather than asking the caller to trust it.

Two fits are run on every table. The plan's fit is a uniform random sample of rows. The second uses
the same number of rows taken from the large-|t| tail, because the sensitivity of `log10 p` to dof
grows roughly with `t^3`: a uniform sample is dominated by rows near p = 1, where every candidate
dof gives nearly the same p. On a source that prints p to four significant figures the uniform fit
cannot separate dof +/- 1 from its own printing noise while the tail fit separates them by an order
of magnitude, so running both is how the fit knows which case it is in.

The result therefore answers two questions, and they are not the same question. `usable` is the
plan's rule and decides what the experiment JSON stores: does this dof rebuild the source's p
closely enough. `identified` says whether the data pins the integer at all. A fit can be usable and
not identified, and that is worth reporting rather than hiding behind one boolean.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import stdtr

# Above this residual the experiment stores `dof: null` and the reader uses the stored -log10(p)
# directly, which is what it does today anyway (plan 2, "Degrees of freedom"; CONTRACT.md).
#
# The limit gates the right thing, which is not "is this integer the source's dof" but "does this
# integer rebuild the source's p", because the reader turns p back into t and t into the slope.
# 0.01 in log10 p is a rebuilt slope wrong by about 1.7e-3 of itself around t = 5, a few times the
# codec's own quantization of the slope (4.47e-3 of slope_se, SPEC section 9). So it is a real
# ceiling, not a formality, and a fit near it is worth refusing.
RESIDUAL_LIMIT = 0.01
# The runner-up (always dof +/- 1) must fit at least this much worse before the integer itself is
# called identified. This does NOT gate storing the dof, because the cost of being one out is small
# and confined: reading TOPCHeF eQTL's 435 as 434 moves the rebuilt slope by 2.6e-6 of slope_se at
# |t| = 1, 1.7e-4 at |t| = 5, and 2.3e-2 at |t| = 30, so it stays inside the codec's own 4.47e-3
# budget up to about |t| = 15 and exceeds it only on the steepest rows. Being twenty out does not:
# the same row at |t| = 30 moves by 0.49 of slope_se. `identified` says which case the fit is in.
IDENTIFIED_MARGIN = 2.0
DEFAULT_SAMPLE = 100_000

# Grid half-width below n_samples. The plan proposed n-60 to n-1. TOPCHeF has n = 516 (2n = 1032
# from ma_count/maf on every sampled row) and dof 435 for eQTL, which is n-81: 70 expression PCs,
# 2, and 9 other covariates. So the proposed grid, 456 to 515, excludes the published answer on the
# one dataset where the answer is known, and the search lands on its own lower edge. sQTL's 480 is
# n-36 and does fall inside, which is exactly how a too-narrow grid hides: it is right for the type
# with fewer phenotype PCs and silently wrong for the other.
#
# 200 covers any QTL model that is sane for n in the hundreds and costs 200 medians. `search`
# widens further when the minimum lands near an edge, so this is a starting box, not a claim about
# how many covariates a pipeline may fit.
DEFAULT_WIDTH = 200
EDGE = 3                     # argmin this close to a grid end means the box was too small
MIN_DOF = 2                  # t with dof < 2 has no variance; no QTL model produces one
MAX_DOF = 100_000

# Row filter. p at or below P_FLOOR has underflowed in the source and carries no t. p above P_CEIL
# is |t| below ~0.67, where the candidate dof values are indistinguishable; including those rows
# would drag the median toward zero for every dof and flatten the very minimum being measured.
P_FLOOR = 1e-300
P_CEIL = 0.5
TINY = 1e-323                # keeps log10 finite when a candidate dof underflows the computed p

NEIGHBOURS = (-20, -5, -1, 1, 5, 20)


def usable(beta, se, pvalue) -> np.ndarray:
    """Rows that carry information about dof."""
    beta = np.asarray(beta, dtype=np.float64)
    se = np.asarray(se, dtype=np.float64)
    p = np.asarray(pvalue, dtype=np.float64)
    return np.isfinite(beta) & np.isfinite(se) & (se > 0) & (p > P_FLOOR) & (p < P_CEIL)


def statistics(beta, se, pvalue) -> tuple[np.ndarray, np.ndarray]:
    """(|t|, log10 p) over the usable rows."""
    ok = usable(beta, se, pvalue)
    beta = np.asarray(beta, dtype=np.float64)[ok]
    se = np.asarray(se, dtype=np.float64)[ok]
    p = np.asarray(pvalue, dtype=np.float64)[ok]
    return np.abs(beta / se), np.log10(p)


def objective(t: np.ndarray, lp: np.ndarray, grid: np.ndarray, chunk: int = 64) -> np.ndarray:
    """Median |log10 p ratio| for each dof in `grid`. Chunked over dof so the working array stays
    `chunk * len(t)` doubles however wide the grid is."""
    grid = np.asarray(grid, dtype=np.float64)
    obj = np.empty(len(grid))
    with np.errstate(divide="ignore"):
        for i in range(0, len(grid), chunk):
            d = grid[i:i + chunk, None]
            p = 2.0 * stdtr(d, -t[None, :])
            obj[i:i + chunk] = np.median(np.abs(np.log10(np.maximum(p, TINY)) - lp[None, :]), axis=1)
    return obj


def grid_for(n_samples: int, width: int = DEFAULT_WIDTH) -> np.ndarray:
    """Integer dof from n_samples - width to n_samples - 1."""
    if n_samples is None or n_samples < MIN_DOF + 1:
        raise ValueError(f"n_samples must be at least {MIN_DOF + 1}, got {n_samples}")
    return np.arange(max(MIN_DOF, n_samples - width), n_samples, dtype=np.int64)


def search(t: np.ndarray, lp: np.ndarray, grid: np.ndarray, hi_limit: int = MAX_DOF,
           rounds: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate `grid`, widening it while the minimum sits within EDGE of an end. The objective is
    V-shaped in dof (the p a candidate dof predicts moves monotonically with 1/dof at fixed t), so
    a minimum against an edge means the box was drawn too small, not that the answer is there.
    Returns the final (grid, objective)."""
    grid = np.asarray(grid, dtype=np.int64)
    obj = objective(t, lp, grid)
    for _ in range(rounds):
        k = int(np.argmin(obj))
        lo, hi = int(grid[0]), int(grid[-1])
        width = max(len(grid), 16)
        if k <= EDGE and lo > MIN_DOF:
            add = np.arange(max(MIN_DOF, lo - width), lo, dtype=np.int64)
        elif k >= len(grid) - 1 - EDGE and hi < hi_limit:
            add = np.arange(hi + 1, min(hi_limit, hi + width) + 1, dtype=np.int64)
        else:
            break
        if len(add) == 0:
            break
        grid = np.concatenate([add, grid]) if add[0] < lo else np.concatenate([grid, add])
        obj = np.concatenate([objective(t, lp, add), obj]) if add[0] < lo else \
            np.concatenate([obj, objective(t, lp, add)])
    return grid, obj


def _subsample(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    return np.arange(n) if n <= k else rng.choice(n, size=k, replace=False)


def _report(grid: np.ndarray, obj: np.ndarray, t: np.ndarray, lp: np.ndarray) -> dict:
    """One grid search turned into an answer plus the evidence that it is the answer."""
    k = int(np.argmin(obj))
    dof = int(grid[k])
    best = float(obj[k])
    order = np.argsort(obj)
    second = int(order[1]) if len(order) > 1 else k
    by_dof = {int(d): float(o) for d, o in zip(grid, obj)}
    # How much worse each neighbour is. This is the number that says whether the minimum is a point
    # or a basin: on a dataset whose dof is unknown, a margin near 1 means the fit cannot tell.
    neighbours = {f"{o:+d}": by_dof.get(dof + o) for o in NEIGHBOURS}
    return {
        "dof": dof,
        "residual_log10p": best,
        "n": len(t),
        "grid": [int(grid[0]), int(grid[-1])],
        "at_edge": k <= EDGE or k >= len(grid) - 1 - EDGE,
        "runner_up": {"dof": int(grid[second]), "residual_log10p": float(obj[second])},
        "margin": float(obj[second] / best) if best > 0 else float("inf"),
        "neighbours": neighbours,
        "neighbour_ratios": {kk: (None if v is None else (float(v / best) if best > 0 else float("inf")))
                             for kk, v in neighbours.items()},
        "min_abs_t": float(np.min(t)) if len(t) else None,
        "median_abs_t": float(np.median(t)) if len(t) else None,
    }


def fit(beta, se, pvalue, n_samples: int | None = None, grid=None, sample: int = DEFAULT_SAMPLE,
        seed: int = 0, width: int = DEFAULT_WIDTH) -> dict:
    """Fit the integer dof of one nominal statistics table.

    `grid` overrides the default box drawn from `n_samples`; one of the two is required. The result
    carries `dof` and `residual_log10p` for the experiment JSON, and `usable` for the decision the
    plan asks for: a residual above RESIDUAL_LIMIT means store `dof: null`.
    """
    if grid is None:
        grid = grid_for(n_samples, width)
    grid = np.asarray(grid, dtype=np.int64)
    n_rows = len(np.asarray(pvalue))
    t_all, lp_all = statistics(beta, se, pvalue)
    hi_limit = min(MAX_DOF, n_samples - 1) if n_samples else MAX_DOF
    out: dict = {"dof": None, "residual_log10p": None, "usable": False, "rows": n_rows,
                 "rows_usable": len(t_all), "n_samples": n_samples, "sample_requested": sample}
    if len(t_all) == 0:
        out["reason"] = "no rows with a finite beta, a positive se and 0 < p < 0.5"
        return out

    rng = np.random.default_rng(seed)
    pick = _subsample(len(t_all), sample, rng)
    uniform_grid, uniform_obj = search(t_all[pick], lp_all[pick], grid, hi_limit)
    uni = _report(uniform_grid, uniform_obj, t_all[pick], lp_all[pick])

    # The same count of rows from the large-|t| tail, where dof actually bites.
    if len(t_all) <= sample:
        tail_idx = np.arange(len(t_all))
    else:
        tail_idx = np.argpartition(t_all, len(t_all) - sample)[len(t_all) - sample:]
    tail_grid, tail_obj = search(t_all[tail_idx], lp_all[tail_idx], grid, hi_limit)
    tail = _report(tail_grid, tail_obj, t_all[tail_idx], lp_all[tail_idx])

    # Score the answer on the tail rows too, at the dof the plan's uniform sample chose. The tail
    # is where a wrong dof does its damage, so this is the residual that bounds the rebuilt slope.
    tail_at_uni = float(objective(t_all[tail_idx], lp_all[tail_idx], [uni["dof"]])[0])
    residual = max(uni["residual_log10p"], tail_at_uni)

    notes = []
    if uni["at_edge"] or tail["at_edge"]:
        notes.append("the minimum sits against a grid edge even after widening")
    if uni["dof"] != tail["dof"]:
        notes.append(f"the uniform sample fits dof {uni['dof']}, the large-|t| tail fits {tail['dof']}")
    if tail["margin"] < IDENTIFIED_MARGIN:
        notes.append(f"dof {tail['runner_up']['dof']} fits only {tail['margin']:.3g}x worse on the tail; "
                     "the data does not pin the integer")
    if residual > RESIDUAL_LIMIT:
        notes.append(f"residual {residual:.3g} is above the {RESIDUAL_LIMIT} log10 p limit")

    out.update({
        "dof": uni["dof"],
        "residual_log10p": residual,
        "residual_uniform": uni["residual_log10p"],
        "residual_tail_at_dof": tail_at_uni,
        # Store the dof when it rebuilds p well enough, which is the plan's rule ...
        "usable": residual <= RESIDUAL_LIMIT and not (uni["at_edge"] or tail["at_edge"]),
        # ... and separately say whether the integer is pinned, which is a different question.
        "identified": uni["dof"] == tail["dof"] and tail["margin"] >= IDENTIFIED_MARGIN,
        "margin": tail["margin"],
        "uniform": uni,
        "tail": tail,
    })
    if n_samples:
        out["implied_covariates"] = n_samples - uni["dof"]
    if notes:
        out["reason"] = "; ".join(notes)
    return out


def dof_for_manifest(result: dict) -> int | None:
    """What the experiment JSON stores: the fitted dof, or None when the fit is not trustworthy."""
    return result["dof"] if result.get("usable") else None


def sharpness_table(result: dict, which: str = "uniform") -> str:
    """The objective around the optimum, as the plan's step-4 evidence."""
    r = result[which]
    head = (f"  dof {r['dof']:>6d}  residual {r['residual_log10p']:.3e}   "
            f"(n = {r['n']}, grid {r['grid'][0]}-{r['grid'][1]})")
    rows = [head]
    for off in NEIGHBOURS:
        v = r["neighbours"][f"{off:+d}"]
        ratio = r["neighbour_ratios"][f"{off:+d}"]
        rows.append(f"  dof {r['dof'] + off:>6d}  residual {'-' if v is None else f'{v:.3e}'}"
                    f"   x{'-' if ratio is None else f'{ratio:.4g}'} worse")
    return "\n".join(rows)


# ---- CLI -------------------------------------------------------------------------------------

def sample_parquet(patterns: list[str], cols: tuple[str, str, str], max_rows: int, seed: int,
                   groups_wanted: int = 40):
    """Rows from row groups drawn at random across every matching file.

    Row groups because parquet cannot skip within one; at random across files so the sample is not
    one end of the genome. TOPCHeF's raw nominal files hold about 800k rows per group, so taking
    whole groups would fill a 4M-row budget from two files and call that genome-wide. Each group
    contributes at most `max_rows / groups_wanted` rows instead, which costs more read and buys a
    sample spread over tens of windows on many chromosomes."""
    import glob as _glob

    import pyarrow as pa
    import pyarrow.parquet as pq

    files = sorted({f for pat in patterns for f in _glob.glob(pat)})
    if not files:
        sys.exit(f"dof: no parquet files matched {patterns}")
    groups = []
    for f in files:
        groups += [(f, i) for i in range(pq.ParquetFile(f).metadata.num_row_groups)]
    rng = np.random.default_rng(seed)
    rng.shuffle(groups)
    # Spread over `groups_wanted` groups, or over every group there is when the file is small;
    # dividing by the constant would leave the budget unfilled on a table with few row groups.
    per_group = max(1, max_rows // max(1, min(groups_wanted, len(groups))))
    parts, n = [], 0
    schema = pa.schema([(c, pa.float64()) for c in cols])
    for f, i in groups:
        tb = pq.ParquetFile(f).read_row_group(i, columns=list(cols)).cast(schema)
        if tb.num_rows > per_group:
            tb = tb.take(np.sort(rng.choice(tb.num_rows, size=per_group, replace=False)))
        parts.append(tb)
        n += tb.num_rows
        if n >= max_rows:
            break
    tb = pa.concat_tables(parts)
    arrays = [tb[c].to_numpy(zero_copy_only=False) for c in cols]
    return arrays, {"files": len(files), "files_sampled": len({f for f, _ in groups[:len(parts)]}),
                    "row_groups_read": len(parts), "rows_read": tb.num_rows,
                    "path": str(Path(files[0]).parent)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parquet", nargs="+", help="parquet file or glob holding nominal statistics")
    ap.add_argument("--n-samples", type=int, help="samples in the experiment; sets the default grid")
    ap.add_argument("--dof-range", nargs=2, type=int, metavar=("LO", "HI"), help="explicit dof grid, inclusive")
    ap.add_argument("--width", type=int, default=DEFAULT_WIDTH, help=f"grid width below n_samples (default {DEFAULT_WIDTH})")
    ap.add_argument("--beta", default="beta", help="effect size column (default beta)")
    ap.add_argument("--se", default="se", help="standard error column (default se)")
    ap.add_argument("--p", default="pvalue", help="p-value column (default pvalue)")
    ap.add_argument("--sample", type=int, default=DEFAULT_SAMPLE, help=f"rows per fit (default {DEFAULT_SAMPLE})")
    ap.add_argument("--max-rows", type=int, default=4_000_000, help="rows to read before subsampling")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", help="write the full result here")
    args = ap.parse_args(argv)

    if args.n_samples is None and args.dof_range is None:
        sys.exit("dof: give --n-samples or --dof-range")
    (beta, se, p), where = sample_parquet(args.parquet, (args.beta, args.se, args.p), args.max_rows, args.seed)
    grid = np.arange(args.dof_range[0], args.dof_range[1] + 1) if args.dof_range else None
    r = fit(beta, se, p, n_samples=args.n_samples, grid=grid, sample=args.sample,
            seed=args.seed, width=args.width)
    r["source"] = where
    print(f"{where['path']}: {where['rows_read']} rows read from {where['row_groups_read']} row groups "
          f"of {where['files']} files, {r['rows_usable']} usable")
    if r["dof"] is None:
        print(f"no fit: {r['reason']}")
    else:
        print(f"dof {r['dof']}  residual {r['residual_log10p']:.3e}  usable={r['usable']}  "
              f"identified={r['identified']}  margin x{r['margin']:.4g}"
              + (f"  ({r['reason']})" if r.get("reason") else ""))
        if r.get("implied_covariates") is not None:
            print(f"implied covariates (n_samples - dof): {r['implied_covariates']}")
        for which in ("uniform", "tail"):
            print(f"{which} sample:")
            print(sharpness_table(r, which))
    print(f"stored dof: {dof_for_manifest(r)}")
    if args.json:
        Path(args.json).write_text(json.dumps(r, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
