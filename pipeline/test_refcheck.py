"""Tests for the refget allele check (pipeline/steps_refget.py).

    uv run python -m pipeline.test_refcheck

Plain asserts, no test framework, synthetic data only.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from . import steps_refget as rg

# 1-based: A at 1, C at 2, ... repeating ACGT; an N at 20
SEQ = bytearray(b"ACGT" * 50)
SEQ[19] = ord("N")
SEQ = np.frombuffer(bytes(SEQ), dtype=np.uint8)

ROWS = [  # position, A1, A2, expected class
    (1, "A", "G", "a1"),          # SNP, A1 is the reference
    (2, "T", "C", "a2"),          # SNP, A2 is the reference
    (3, "GT", "G", "a1"),         # deletion: both match by prefix, the longer allele wins
    (3, "GAA", "G", "a2"),        # insertion: only the shorter allele reads GT...
    (4, "t", "c", "a1"),          # lower-case alleles
    (5, "A", "A", "both"),        # A1 == A2, allowed into the heap
    (6, "A", "T", "neither"),     # reference is C, no flip either
    (7, "C", "T", "neither"),     # reference is G: C is its complement, a strand flip
    (20, "A", "C", "neither"),    # an N in the reference
    (9, None, None, "unchecked"), # trans eQTL-only position, no alleles
    (500, "A", "C", "neither"),   # past the end of the sequence
    (200, "TTT", "T", "a2"),      # a long allele running off the end does not match
]


def classify_cases():
    pos = np.array([r[0] for r in ROWS], dtype=np.int64)
    a1, a2 = [r[1] for r in ROWS], [r[2] for r in ROWS]
    got = [rg.MATCH[c] for c in rg.classify(SEQ, pos, a1, a2)]
    want = [r[3] for r in ROWS]
    assert got == want, "\n".join(f"{r[:3]}: got {g}, want {w}" for r, g, w in zip(ROWS, got, want) if g != w)
    detail = rg.neither_detail(SEQ, pos, a1, a2, rg.classify(SEQ, pos, a1, a2))
    assert detail == {"strand_flip": 1, "on_n": 1}, detail


def threshold():
    ok = {"a1": 998, "a2": 1, "both": 1, "neither": 0}
    assert rg.check_threshold(ok, 0.999) == 1.0
    low = {"a1": 990, "a2": 0, "both": 0, "neither": 10}
    try:
        rg.check_threshold(low, 0.999)
    except ValueError as e:
        assert "min_match_fraction" in str(e), e
    else:
        raise AssertionError("no ValueError below the minimum")


def store_round_trip():
    """A real (tiny) refgetstore gives the same bases `classify` expects."""
    from gtars.refget import RefgetStore
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / "t.fa"
        fa.write_text(">chrT test\n" + bytes(SEQ).decode().lower() + "\n")
        store = RefgetStore.on_disk(str(Path(d) / "store"))
        store.set_quiet(True)
        meta, _ = store.add_sequence_collection_from_fasta(str(fa))
        (rec,) = list(store.get_collection(meta.digest))
        seq = rg.load_sequence(store, rec.metadata.sha512t24u, rec.metadata.length)
        assert np.array_equal(seq, SEQ), "store bases differ from the FASTA"
        assert len(rec.metadata.sha512t24u) == 32 and len(meta.digest) == 32


class _Cfg:
    """Just enough Config for `spot_check`: a `_tables` dir and `paper_variants`."""

    def __init__(self, derived: Path, paper_variants: dict):
        self.derived = derived
        self.tables = derived / "_tables"
        self._d = {"paper_variants": paper_variants}

    def __getitem__(self, k):
        return self._d[k]


def spot_check_cases():
    """`spot_check` compares a stored call against bytes read back from the store, and complains
    when they disagree."""
    import pyarrow as pa
    from gtars.refget import RefgetStore
    from .common import write_parquet

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fa = d / "t.fa"
        fa.write_text(">chrT test\n" + bytes(SEQ).decode().lower() + "\n")
        store = RefgetStore.on_disk(str(d / "store"))
        store.set_quiet(True)
        meta, _ = store.add_sequence_collection_from_fasta(str(fa))
        (rec,) = list(store.get_collection(meta.digest))
        seqs = {"chrT": {"digest": rec.metadata.sha512t24u, "length": rec.metadata.length}}

        # position 1 is A (A1), position 2 is C (A2); the third row lies about position 6
        rows = [(1, "A", "G", "a1"), (2, "T", "C", "a2"), (6, "A", "T", "a1")]
        write_parquet(pa.table({
            "chr": pa.array(["chrT"] * len(rows), pa.string()),
            "position": pa.array([r[0] for r in rows], pa.int32()),
            "A1": pa.array([r[1] for r in rows], pa.string()),
            "A2": pa.array([r[2] for r in rows], pa.string()),
            "in_cis": pa.array([True] * len(rows), pa.bool_()),
            "match": pa.array([r[3] for r in rows], pa.string()).dictionary_encode(),
        }), rg.refcheck_table(_Cfg(d, {}), "chrT"), 1000)

        results = []
        def check(ok, msg):
            results.append((ok, msg))

        rg.spot_check(_Cfg(d, {"chrT:1": "rs1", "chrT:2": "rs2"}), seqs, store, check)
        assert results[-1][0], results[-1][1]

        results.clear()
        rg.spot_check(_Cfg(d, {"chrT:6": "rs6"}), seqs, store, check)
        assert not results[-1][0], "a wrong call was not caught"
        assert "called a1" in results[-1][1], results[-1][1]

        # a position with no row in the table is a failure, not a silent skip
        results.clear()
        rg.spot_check(_Cfg(d, {"chrT:99": "rs99"}), seqs, store, check)
        assert not results[-1][0], "a missing variant was not caught"


CASES = [classify_cases, threshold, store_round_trip, spot_check_cases]


def main() -> int:
    for fn in CASES:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(CASES)} cases passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
