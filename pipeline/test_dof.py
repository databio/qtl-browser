"""Tests for the degrees-of-freedom fit (pipeline/dof.py).

    uv run python -m pipeline.test_dof

Plain asserts, no test framework, synthetic data only. The real-data check (TOPCHeF eQTL must fit
435 and sQTL 480) is not here: it needs the nominal tables, and lives in the CLI.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.special import stdtr

from . import dof

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


def synthetic(true_dof: int, n: int = 400_000, seed: int = 1):
    """A nominal statistics table whose p is exactly the two-sided Student-t tail of beta/se.
    This is the relation SPEC section 5 assumes; the fit has to find `true_dof` back out of it."""
    rng = np.random.default_rng(seed)
    t = rng.standard_t(true_dof, size=n)
    se = np.exp(rng.normal(-2.0, 0.6, size=n))          # SEs spread over a couple of decades
    beta = t * se
    p = 2.0 * stdtr(float(true_dof), -np.abs(t))
    return beta, se, p


def round_sig(x: np.ndarray, digits: int) -> np.ndarray:
    """x printed to `digits` significant figures, as a source's text output would hold it."""
    e = np.floor(np.log10(np.abs(x)))
    return np.round(x / 10.0 ** e, digits - 1) * 10.0 ** e


@case
def recovers_a_known_dof():
    beta, se, p = synthetic(419)
    r = dof.fit(beta, se, p, n_samples=500, sample=20_000)
    assert r["dof"] == 419, r
    assert r["usable"] and r["identified"], r
    assert dof.dof_for_manifest(r) == 419
    assert r["implied_covariates"] == 81
    assert r["uniform"]["dof"] == r["tail"]["dof"] == 419, r
    assert not r["uniform"]["at_edge"], r["uniform"]


@case
def the_minimum_is_sharp():
    """A wrong dof has to stand out, or the fit means nothing on a dataset with no known answer.
    Measured on float32 data, because exact doubles put the residual at the optimum at 0.0 and any
    ratio against it is infinite; float32 is the noise floor real sources carry."""
    beta, se, p = synthetic(419)
    r = dof.fit(beta.astype(np.float32), se.astype(np.float32), p.astype(np.float32),
                n_samples=500, sample=20_000)
    for which in ("uniform", "tail"):
        g = r[which]
        assert g["margin"] > 10, f"{which}: runner-up only {g['margin']:.2f}x worse: {g}"
        for off in ("-20", "-5", "-1", "+1", "+5", "+20"):
            ratio = g["neighbour_ratios"][off]
            assert ratio is not None and ratio > 10, f"{which} dof{off}: only {ratio}x worse: {g}"
        # The neighbourhood is monotone away from the optimum: a V, not a rough basin.
        n = g["neighbours"]
        assert n["+20"] > n["+5"] > n["+1"] > g["residual_log10p"], g
        assert n["-20"] > n["-5"] > n["-1"] > g["residual_log10p"], g


@case
def float32_source_still_recovers():
    """TOPCHeF's tables carry float32 beta, se and p. That rounding is the noise floor packcheck
    measured at 5.96e-08 relative, and it must not move the answer."""
    beta, se, p = synthetic(419)
    b32 = beta.astype(np.float32).astype(np.float64)
    s32 = se.astype(np.float32).astype(np.float64)
    p32 = p.astype(np.float32).astype(np.float64)
    r = dof.fit(b32, s32, p32, n_samples=500, sample=20_000)
    assert r["dof"] == 419 and r["usable"], r
    assert r["uniform"]["margin"] > 2, r["uniform"]


@case
def a_coarsely_printed_p_says_so():
    """A source that prints p to four significant figures buries the dof +/- 1 separation in its
    own rounding. The fit is still usable - a unit either way rebuilds p far inside the limit - but
    it must not claim the integer is identified."""
    beta, se, p = synthetic(419)
    r = dof.fit(beta, se, round_sig(p, 4), n_samples=500, sample=20_000)
    assert abs(r["dof"] - 419) <= 1, r
    assert r["residual_log10p"] < dof.RESIDUAL_LIMIT and r["usable"], r
    assert not r["identified"], r
    assert "does not pin" in r["reason"], r
    # the large-|t| tail is the half of the fit that still carries the signal
    assert abs(r["tail"]["dof"] - 419) <= abs(r["uniform"]["dof"] - 419), r


@case
def six_significant_figures_still_identifies():
    """Six digits is enough for the tail to separate neighbours by a wide margin."""
    beta, se, p = synthetic(419)
    r = dof.fit(beta, se, round_sig(p, 6), n_samples=500, sample=20_000)
    assert r["dof"] == 419 and r["usable"] and r["identified"], r
    assert r["tail"]["margin"] > 10, r["tail"]["margin"]


@case
def widens_past_the_plans_grid():
    """The plan proposed n-60 to n-1. TOPCHeF is n-81, so a fitter that trusted those bounds would
    have returned the edge. `search` has to walk out of a box drawn too small."""
    beta, se, p = synthetic(419)
    narrow = np.arange(500 - 60, 500)                  # 440..499, the true 419 is outside it
    r = dof.fit(beta, se, p, n_samples=500, grid=narrow, sample=20_000)
    assert r["dof"] == 419, r
    assert r["usable"], r
    assert r["uniform"]["grid"][0] <= 419, r["uniform"]["grid"]


@case
def an_unrelated_p_is_refused():
    """p that does not come from beta/se at all: the residual blows past the limit and the
    experiment JSON stores `dof: null`."""
    beta, se, _ = synthetic(419, n=100_000)
    rng = np.random.default_rng(7)
    p = rng.uniform(1e-8, 0.4, size=len(beta))
    r = dof.fit(beta, se, p, n_samples=500, sample=10_000)
    assert r["residual_log10p"] > dof.RESIDUAL_LIMIT, r
    assert not r["usable"], r
    assert dof.dof_for_manifest(r) is None
    assert "residual" in r["reason"] or "does not pin" in r["reason"], r


@case
def row_filter():
    beta = np.array([1.0, 1.0, np.nan, 1.0, 1.0, 1.0, 1.0])
    se = np.array([0.1, 0.0, 0.1, -0.1, 0.1, 0.1, 0.1])
    p = np.array([0.01, 0.01, 0.01, 0.01, np.nan, 0.9, 0.0])
    assert list(dof.usable(beta, se, p)) == [True, False, False, False, False, False, False]
    t, lp = dof.statistics(beta, se, p)
    assert len(t) == 1 and abs(t[0] - 10.0) < 1e-12 and abs(lp[0] + 2.0) < 1e-12


@case
def no_usable_rows():
    r = dof.fit(np.array([1.0]), np.array([0.1]), np.array([0.9]), n_samples=500)
    assert r["dof"] is None and not r["usable"] and r["rows_usable"] == 0, r
    assert dof.dof_for_manifest(r) is None


@case
def grid_bounds():
    g = dof.grid_for(500)
    assert g[0] == 300 and g[-1] == 499 and len(g) == 200
    assert dof.grid_for(500, width=60)[0] == 440
    assert dof.grid_for(10)[0] == dof.MIN_DOF                 # never below a t with no variance
    raises("n_samples must be at least", dof.grid_for, 2)
    raises("n_samples must be at least", dof.grid_for, None)


@case
def objective_shape():
    """The objective is V-shaped in dof, which is what lets `search` trust an interior minimum."""
    beta, se, p = synthetic(419, n=50_000, seed=3)
    t, lp = dof.statistics(beta, se, p)
    grid = np.arange(380, 461)
    obj = dof.objective(t, lp, grid)
    k = int(np.argmin(obj))
    assert grid[k] == 419, grid[k]
    assert np.all(np.diff(obj[:k]) < 0), "not decreasing up to the minimum"
    assert np.all(np.diff(obj[k:]) > 0), "not increasing after the minimum"
    # chunking must not change the answer
    assert np.allclose(obj, dof.objective(t, lp, grid, chunk=7))


@case
def cli_sampler_spreads_over_row_groups():
    """A 4M-row budget must not come out of two row groups: TOPCHeF's raw files hold ~800k rows
    per group, so whole groups would sample two windows and call the answer genome-wide."""
    import tempfile
    from pathlib import Path

    import pyarrow as pa
    import pyarrow.parquet as pq

    beta, se, p = synthetic(419, n=200_000, seed=5)
    with tempfile.TemporaryDirectory() as d:
        for i in range(4):                       # four "chromosomes", 10 row groups each
            a, b = i * 50_000, (i + 1) * 50_000
            pq.write_table(pa.table({"beta": beta[a:b], "se": se[a:b], "pvalue": p[a:b]}),
                           Path(d) / f"chr{i}.parquet", row_group_size=5_000)
        (arrays, where) = dof.sample_parquet([str(Path(d) / "*.parquet")], ("beta", "se", "pvalue"),
                                             max_rows=40_000, seed=0)
        assert where["files"] == 4 and where["files_sampled"] == 4, where
        assert where["row_groups_read"] >= 20, where     # not one group per file
        assert 40_000 <= where["rows_read"] < 45_000, where
        r = dof.fit(*arrays, n_samples=500, sample=20_000)
        assert r["dof"] == 419, r

        # a budget larger than the table returns the whole table, not a fraction of it
        (_, all_of_it) = dof.sample_parquet([str(Path(d) / "*.parquet")], ("beta", "se", "pvalue"),
                                            max_rows=10_000_000, seed=0)
        assert all_of_it["rows_read"] == 200_000, all_of_it


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
