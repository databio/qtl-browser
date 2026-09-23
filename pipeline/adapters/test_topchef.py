"""Tests for the TOPCHeF adapter (pipeline/adapters/topchef.py).

    uv run python -m pipeline.adapters.test_topchef

Plain asserts, no test framework, synthetic data only. The rules worth pinning here are the ones a
genome-wide build cannot show you: which allele becomes `ref`, whether anything swaps, and what a
swap would do to `beta` and `af` if one ever happened.
"""
from __future__ import annotations

import sys

import numpy as np

from .. import qtlstore as qs
from . import topchef as tc

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


def _rows(rows):
    """(match, A1, A2) triples as the three object arrays `reference_allele` takes."""
    return (np.array([r[0] for r in rows], dtype=object),
            np.array([r[1] for r in rows], dtype=object),
            np.array([r[2] for r in rows], dtype=object))


@case
def snps_take_the_refcheck_call():
    """At a SNP the refcheck class is a reading of the reference, so it is followed either way."""
    match, a1, a2 = _rows([("a2", "G", "A"), ("a1", "G", "A"), ("both", "A", "A"),
                           ("neither", "C", "T"), ("unchecked", "", "")])
    got = tc.reference_allele(match, a1, a2)
    assert list(got) == ["A", "G", "A", None, None], list(got)


@case
def disputed_is_only_the_a1_indels():
    match, a1, a2 = _rows([("a1", "G", "A"),        # SNP: a reading, not a tie-break
                           ("a1", "GT", "G"),       # indel called a1: the tie-break fired
                           ("a2", "GT", "G"),       # indel called a2: unambiguous
                           ("neither", "GT", "C")])
    assert list(tc.disputed(match, a1, a2)) == [False, True, False, False]


@case
def a_disputed_indel_takes_a2_only_when_a2_reads():
    """The tie-break is settled by asking the reference, not by preferring an answer.

    `classify` takes the longer allele when both prefix-match, which puts `ref` on A1 at 97,433 cis
    indels. `af` says A1 is the minor allele at those sites (mean 0.1504 where dbSNP even puts REF
    on A1), so A1 is ALT and A2 is the reference -- but only where the reference really does read
    A2, which is what `a2_reads` carries. Where it does not, the refcheck call stands.
    """
    match, a1, a2 = _rows([("a1", "GT", "G"), ("a1", "GT", "G")])
    reads = np.array([True, False])
    assert list(tc.reference_allele(match, a1, a2, reads)) == ["G", "GT"]
    # with no check run, nothing is overridden: every disputed indel keeps the refcheck call
    assert list(tc.reference_allele(match, a1, a2, None)) == ["GT", "GT"]


@case
def topchef_orientation_swaps_nothing_and_flips_no_sign():
    """The adapter's own call through `orient_to_ref`, on the shapes the release actually has.

    A2 is REF and `af`/`slope` are A1's, so every row must come out as-is: `ref` = A2, `alt` = A1,
    `beta` untouched, `af` untouched. A build whose slopes come out negated has inverted every
    effect size in the store, and nothing in a p-value or an SE would show it.
    """
    rows = [("a2", "G", "A", 0.5, 0.10),            # SNP
            ("a2", "T", "C", -0.2, 0.25),           # SNP, negative slope
            ("a1", "GT", "G", 0.7, 0.08),           # disputed indel, reference also reads A2
            ("a2", "G", "GTT", 0.1, 0.40),          # indel, unambiguous
            ("both", "A", "A", -0.3, 0.80)]         # A1 == A2
    match, a1, a2 = _rows([(r[0], r[1], r[2]) for r in rows])
    beta = np.array([r[3] for r in rows])
    af = np.array([r[4] for r in rows])
    reads = np.array([False, False, True, False, False])          # only the disputed row is asked
    r = qs.orient_to_ref(tc.reference_allele(match, a1, a2, reads), a1, a2, beta, af)
    assert r["counts"] == {"as_is": 5, "swapped": 0, "dropped": 0}, r["counts"]
    assert list(r["ref"]) == list(a2) and list(r["alt"]) == list(a1)
    assert np.array_equal(r["beta"], beta), "a slope was negated: the orientation is inverted"
    assert np.array_equal(r["af"], af), "af was mirrored: the orientation is inverted"


@case
def a_swap_would_negate_beta_and_mirror_af():
    """What the SQL in `_swap` and `_swap_af` has to do, checked against `orient_to_ref` itself.

    TOPCHeF never takes this path. It is here so that a source that does -- the next adapter, or a
    TOPCHeF rebuild where the indel rule changes -- cannot pick up a half-applied swap.
    """
    match, a1, a2 = _rows([("a1", "G", "A")])       # a SNP whose A1 is the reference: a real swap
    beta, af = np.array([0.5]), np.array([0.1])
    r = qs.orient_to_ref(tc.reference_allele(match, a1, a2), a1, a2, beta, af)
    assert r["counts"] == {"as_is": 0, "swapped": 1, "dropped": 0}
    assert list(r["ref"]) == ["G"] and list(r["alt"]) == ["A"]
    assert r["beta"][0] == -0.5 and abs(r["af"][0] - 0.9) < 1e-12
    assert tc._swap("x") == "CASE WHEN o.swapped THEN -x ELSE x END"
    assert tc._swap_af("x") == "CASE WHEN o.swapped THEN 1.0 - x ELSE x END"


@case
def a_site_with_no_reference_allele_is_dropped():
    """`neither` and `unchecked` have no anchored `ref`, so they get no site row."""
    match, a1, a2 = _rows([("neither", "C", "T"), ("unchecked", "", ""), ("a2", "G", "A")])
    r = qs.orient_to_ref(tc.reference_allele(match, a1, a2), a1, a2, np.zeros(3), np.zeros(3))
    assert list(r["keep"]) == [False, False, True]
    assert r["counts"]["dropped"] == 2


@case
def a_variant_without_source_alleles_is_excluded_not_inferred():
    """A trans eQTL variant the release names only as chr:pos gets no reference allele and no site,
    whatever the other arrays hold: alleles are never guessed for it."""
    match, a1, a2 = _rows([("no_source_alleles", "", ""), ("no_source_alleles", "G", "A"), ("a2", "G", "A")])
    ref = tc.reference_allele(match, a1, a2, np.array([True, True, True]))
    assert ref[0] is None and ref[1] is None and ref[2] == "A"
    r = qs.orient_to_ref(ref, a1, a2, np.ones(3), np.zeros(3))
    assert list(r["keep"]) == [False, False, True]


@case
def source_paths_and_the_experiment_facts():
    """The file names and study constants the v0 steps now read from here instead of spelling out."""
    class Cfg:
        def __init__(self):
            self.cfg = {"packs": {"dof": {"eqtl": 435, "sqtl": 480}},
                        "sig_column": "pval_perm", "sig_threshold": 0.05}

        def __getitem__(self, k):
            return self.cfg[k]

        def raw_dir(self, name):
            from pathlib import Path
            return Path("/raw") / name

        def raw_glob(self, name):
            return f"/raw/{name}/*.parquet"

    cfg = Cfg()
    assert str(tc.source_file(cfg, "e", "nominal", "chr7")).endswith(
        "cis_eQTL_nominal/topchef_chr7_MaxPC70.cis_qtl_pairs.chr7.parquet")
    assert str(tc.source_file(cfg, "s", "susie", "chr7")).endswith(
        "cis_sQTL_SuSiE/topchefSplice_chr7_MaxPC25.SuSiE_summary.parquet")
    assert len(tc.cis_sources()) == 6 and "trans_eQTL" not in tc.cis_sources()
    assert tc.dof(cfg) == {"ge": 435, "leafcutter": 480}
    assert tc.significance(cfg) == {"column": "p_perm", "op": "<", "threshold": 0.05}
    assert tc.PHENOTYPE_TYPE == {"e": "ge", "s": "leafcutter"}
    # the contract leafcutter nominal reads the very file the v0 sQTL packs stream
    from ..steps_pack import _raw_sqtl
    assert _raw_sqtl(cfg, "chr7") == tc.source_file(cfg, "s", "nominal", "chr7")


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
