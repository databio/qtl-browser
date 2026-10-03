"""Tests for the GTF attribute parser and step 2 (pipeline/steps_gtf.py).

    uv run python -m pipeline.test_gtf

Plain asserts, no test framework, synthetic data only. The parser has one job that is easy to get
wrong and silent when it is: GTF writes string attributes quoted and numeric ones bare, and a
quoted-only regex drops every bare key without ever raising. That is exactly how `exon_number`
became 0 on all 1,378,020 rows of the v0 `exons.parquet`.
"""
from __future__ import annotations

import gzip
import re
import sys
import tempfile
from pathlib import Path

import pyarrow.parquet as pq

from . import steps_gtf as sg

CASES = []

# what `steps_gtf` read before the fix: quoted attributes only. Nothing this regex finds today may
# come out differently from the new one -- the fix adds keys, it does not change existing ones.
OLD_ATTR = re.compile(r'(\S+) "([^"]*)"')


def case(fn):
    CASES.append(fn)
    return fn


def gtf_line(chrom, feature, start, end, strand, attrs: str) -> str:
    return "\t".join([chrom, "HAVANA", feature, str(start), str(end), ".", strand, ".", attrs])


# A GENCODE-shaped attribute string: quoted strings, bare numbers, a repeated `tag`, and a quoted
# value holding a space.
REAL = ('gene_id "ENSG00000000001.5"; transcript_id "ENST1.1"; gene_type "protein_coding"; '
        'gene_name "AAA"; exon_number 1; level 2; transcript_support_level "NA"; '
        'havana_gene "OTTHUMG00000000001.2"; tag "basic"; tag "CCDS"; '
        'remark "one two three";')

G1 = 'gene_id "ENSG00000000001.5"; gene_type "protein_coding"; gene_name "AAA";'
G2 = 'gene_id "ENSG00000000002.12"; gene_type "lncRNA"; gene_name "BBB";'
PAR = 'gene_id "ENSG00000000003.7_PAR_Y"; gene_type "protein_coding"; gene_name "CCC";'

LINES = [
    "##description: a fixture, not a release",
    gtf_line("chr1", "gene", 100, 900, "+", G1),
    gtf_line("chr1", "transcript", 100, 900, "+", G1 + ' transcript_id "ENST1.1";'),   # ignored
    gtf_line("chr1", "exon", 100, 200, "+", G1 + ' transcript_id "ENST1.1"; exon_number 1;'),
    gtf_line("chr1", "exon", 800, 900, "+", G1 + ' transcript_id "ENST1.1"; exon_number 2;'),
    gtf_line("chr1", "exon", 100, 250, "+", G1 + ' transcript_id "ENST2.1"; exon_number 1;'),
    gtf_line("chr2", "gene", 500, 1500, "-", G2),
    gtf_line("chr2", "exon", 1400, 1500, "-", G2 + ' transcript_id "ENST3.2"; exon_number 1;'),
    gtf_line("chrY", "gene", 10, 20, "+", PAR),
    gtf_line("chrY", "exon", 10, 20, "+", PAR + ' transcript_id "ENST4.1_PAR_Y"; exon_number 1;'),
]


@case
def quoted_and_bare_attributes():
    a = sg.parse_attrs(REAL)
    assert a["gene_id"] == "ENSG00000000001.5", a
    assert a["gene_name"] == "AAA", a
    assert a["exon_number"] == "1", a                # bare: the key the old regex never produced
    assert a["level"] == "2", a
    assert a["remark"] == "one two three", a         # quoted values keep their spaces
    assert a["tag"] == "CCDS", a                     # a repeated key keeps the last value, as before
    assert int(a.get("exon_number", 0)) == 1, a      # the call site that used to take the default


@case
def nothing_that_parsed_before_changes():
    """Every key the quoted-only regex found still parses to the same value."""
    old = dict(OLD_ATTR.findall(REAL))
    new = sg.parse_attrs(REAL)
    assert old, "fixture has no quoted attributes"
    diff = {k: (v, new.get(k)) for k, v in old.items() if new.get(k) != v}
    assert not diff, diff


@case
def awkward_terminators():
    """A quoted value holding a `;`, a missing final semicolon, and trailing whitespace."""
    a = sg.parse_attrs('gene_name "A;B"; exon_number 3;   ')
    assert a == {"gene_name": "A;B", "exon_number": "3"}, a
    # the old regex parsed an unterminated final quoted attribute; the new one must too
    assert sg.parse_attrs('gene_id "ENSG1.1"; gene_name "AAA"')["gene_name"] == "AAA"
    assert sg.parse_attrs('gene_id "ENSG1.1"; exon_number 7')["exon_number"] == "7"


@case
def agrees_with_the_annotation_parser():
    """`annotation.parse_gtf` reads the same file for the v1 annotation object. If the two parsers
    disagree, one of the two copies of the gene model in the store is wrong."""
    from . import annotation as an
    for s in (REAL, G1, G2, PAR, G1 + ' transcript_id "ENST1.1"; exon_number 12;'):
        assert sg.parse_attrs(s) == an.parse_attrs(s), s


class _Cfg:
    """Just enough Config for `steps_gtf.run`: the GTF path and an output directory."""

    def __init__(self, d: Path):
        self.gtf = d / "x.gtf.gz"
        self.tables = d / "_tables"
        with gzip.open(self.gtf, "wt") as fh:
            fh.write("\n".join(LINES) + "\n")


@case
def run_writes_real_exon_numbers():
    with tempfile.TemporaryDirectory() as d:
        cfg = _Cfg(Path(d))
        sg.run(cfg)
        ex = pq.read_table(cfg.tables / "exons.parquet").to_pylist()
        assert [e["exon_number"] for e in ex] == [1, 1, 2, 1], ex
        assert all(e["gene_id"] != "ENSG00000000003" for e in ex), "_PAR_Y exon kept"
        assert ex[0] == {"gene_id": "ENSG00000000001", "transcript_id": "ENST1.1", "exon_number": 1,
                         "chr": "chr1", "start": 100, "end": 200, "strand": "+"}, ex[0]
        genes = pq.read_table(cfg.tables / "gene_annotation.parquet").to_pylist()
        assert [g["gene_id"] for g in genes] == ["ENSG00000000001", "ENSG00000000002"], genes
        assert [g["symbol"] for g in genes] == ["AAA", "BBB"], genes
        assert [g["biotype"] for g in genes] == ["protein_coding", "lncRNA"], genes
        assert genes[1]["tss"] == 1500, genes[1]          # minus strand: the TSS is the gene's end


def main() -> int:
    for fn in CASES:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(CASES)} cases passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
