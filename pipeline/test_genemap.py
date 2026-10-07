"""Tests for `pipeline.genemap`: symbols and UniProt accessions -> ENSG. Synthetic, no data needed.

    uv run python -m pipeline.test_genemap
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pyarrow as pa

from . import genemap

# the annotation stands in for GENCODE: one plain gene, one renamed, and two sharing a symbol
GENES = pa.table({
    "gene_id": ["ENSG001", "ENSG002", "ENSG003", "ENSG004", "ENSG005"],
    "name": ["FLNC", "VSIR", "DUP", "DUP", "MZB1"],
})

HGNC_HEADER = "hgnc_id\tsymbol\tstatus\talias_symbol\tprev_symbol\tuniprot_ids\tensembl_gene_id"
HGNC_ROWS = [
    # approved, current symbol, accession
    "HGNC:1\tFLNC\tApproved\t\t\tQ14315\tENSG001",
    # renamed: C10orf54 is its previous symbol
    "HGNC:2\tVSIR\tApproved\tVISTA\tC10orf54\tQ9H7M9\tENSG002",
    # names one of the two DUP genes, so it can break that tie
    "HGNC:3\tDUP\tApproved\t\t\tP00001\tENSG003",
    # PACAP is a previous symbol of MZB1 *and* an alias of a gene the annotation lacks
    "HGNC:4\tMZB1\tApproved\t\tPACAP\tQ8WU39\tENSG005",
    "HGNC:5\tADCYAP1\tApproved\tPACAP\t\tP18509\tENSG_ABSENT",
    # withdrawn records are skipped even though they carry a symbol
    "HGNC:6\tGONE\tWithdrawn\t\tFLNC\tQ99999\tENSG001",
]


def _hgnc(d: Path, rows=HGNC_ROWS) -> genemap.Hgnc:
    p = d / "hgnc.txt"
    p.write_text("\n".join([HGNC_HEADER, *rows]) + "\n")
    return genemap.Hgnc(p)


def resolver(d: Path, overrides=None) -> genemap.Resolver:
    return genemap.Resolver(GENES, _hgnc(d), overrides={} if overrides is None else overrides)


def test_unique_symbol_wins_without_help():
    """A symbol exactly one annotation gene carries needs no outside authority."""
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        assert r.resolve("FLNC") == ("ENSG001", "gencode")
        assert r.resolve("FLNC", "Q14315") == ("ENSG001", "gencode")


def test_accession_outranks_a_name():
    """PACAP is MZB1's previous symbol and ADCYAP1's alias. The accession decides, and it decides
    before either name rule runs -- which is the whole reason the rule exists."""
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        assert r.resolve("PACAP", "Q8WU39") == ("ENSG005", "uniprot")
        # without the accession it is a previous symbol of exactly one gene the annotation has
        assert r.resolve("PACAP") == ("ENSG005", "hgnc_prev")


def test_retired_symbol_via_prev():
    with tempfile.TemporaryDirectory() as d:
        assert resolver(Path(d)).resolve("C10orf54") == ("ENSG002", "hgnc_prev")


def test_alias_resolves_only_within_the_annotation():
    """VISTA is an alias of VSIR, which the annotation has. An HGNC record naming a gene the
    annotation lacks resolves to nothing: a `gene_id` the annotation does not carry would dangle."""
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        assert r.resolve("VISTA") == ("ENSG002", "hgnc_alias")
        assert r.resolve("ADCYAP1", "P18509") == (None, "absent")


def test_shared_symbol_is_ambiguous_until_something_picks():
    """Two annotation genes named DUP. HGNC's own record names one, so the tie breaks; an accession
    HGNC does not tie to either leaves it ambiguous rather than guessing."""
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        assert r.resolve("DUP", "P00001") == ("ENSG003", "uniprot")
        assert r.resolve("DUP") == ("ENSG003", "hgnc_current")
        r2 = genemap.Resolver(GENES, _hgnc(Path(d), [x for x in HGNC_ROWS if "\tDUP\t" not in x]), overrides={})
        assert r2.resolve("DUP") == (None, "ambiguous")
        assert r2.resolve("DUP", "P00001") == (None, "ambiguous")


def test_withdrawn_hgnc_records_are_skipped():
    """HGNC:6 is Withdrawn and lists FLNC as a previous symbol. Honouring it would let a retired
    record speak for a symbol the annotation already resolves."""
    with tempfile.TemporaryDirectory() as d:
        h = _hgnc(Path(d))
        assert h.withdrawn == 1
        assert "Q99999" not in h.ensg_of_uniprot


def test_unknown_and_empty():
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        assert r.resolve("NOPE") == (None, "absent")
        assert r.resolve(None) == (None, "absent")
        assert r.resolve("") == (None, "absent")


def test_override_wins_and_a_stale_one_fails_loudly():
    """An override is a person's call where the annotation itself is ambiguous, so it outranks every
    rule -- and naming a gene the annotation lacks raises at construction rather than at use."""
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d), overrides={("DUP", None): "ENSG004"})
        assert r.resolve("DUP", "P00001") == ("ENSG004", "override")
        try:
            resolver(Path(d), overrides={("DUP", None): "ENSG_NOT_HERE"})
        except ValueError as e:
            assert "does not have" in str(e), e
        else:
            raise AssertionError("a stale override must fail the build")


def test_report_counts_every_rule_and_lists_what_failed():
    with tempfile.TemporaryDirectory() as d:
        r = resolver(Path(d))
        rep = r.report([("FLNC", "Q14315"), ("C10orf54", None), ("PACAP", "Q8WU39"),
                        ("NOPE", None), ("FLNC", "Q14315")])          # the repeat is counted once
        assert rep["keys"] == 4 and rep["resolved"] == 3
        assert rep["counts"] == {**dict.fromkeys(genemap.HOW, 0),
                                 "gencode": 1, "hgnc_prev": 1, "uniprot": 1, "absent": 1}
        assert rep["unresolved"] == {"ambiguous": [], "absent": ["NOPE"]}
        assert rep["map"][("FLNC", "Q14315")] == "ENSG001"


def test_shipped_overrides_are_keyed_as_documented():
    """The shipped table is annotation-version specific, so its shape is pinned here: every key a
    (symbol, accession) pair and every value an unversioned ENSG."""
    for (sym, up), ensg in genemap.OVERRIDES.items():
        assert isinstance(sym, str) and sym
        assert up is None or isinstance(up, str)
        assert ensg.startswith("ENSG") and "." not in ensg


CASES = [test_unique_symbol_wins_without_help, test_accession_outranks_a_name,
         test_retired_symbol_via_prev, test_alias_resolves_only_within_the_annotation,
         test_shared_symbol_is_ambiguous_until_something_picks,
         test_withdrawn_hgnc_records_are_skipped, test_unknown_and_empty,
         test_override_wins_and_a_stale_one_fails_loudly,
         test_report_counts_every_rule_and_lists_what_failed,
         test_shipped_overrides_are_keyed_as_documented]


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
