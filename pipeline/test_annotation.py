"""Tests for the annotation builder (pipeline/annotation.py).

    uv run python -m pipeline.test_annotation

Plain asserts, no test framework, same shape as test_qtlstore.py. The fixture is a five-line GTF
that carries every rule the parser has to get right: a minus-strand gene (TSS is `end`), a versioned
id, a `_PAR_Y` duplicate, and two transcripts of one gene.
"""
from __future__ import annotations

import gzip
import random
import struct
import sys
import tempfile
from pathlib import Path

import pyarrow as pa

from . import annotation as an
from . import qtlstore as qs

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


def gtf_line(chrom, feature, start, end, strand, attrs: str) -> str:
    return "\t".join([chrom, "HAVANA", feature, str(start), str(end), ".", strand, ".", attrs])


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


def write_gtf(path: Path, lines=None) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(LINES if lines is None else lines) + "\n")
    return path


def build_fixture(st: qs.Store, annot_id: str = "ann", genes=(("ENSG1", "chr1", 1, 50), ("ENSG2", "chr1", 60, 90))) -> dict:
    """A real annotation (tables, gene models, lookup) in `st`, from a GTF of one-exon genes
    `(gene_id, chr, start, end)` on `+`, named `SYM<n>`: what the other suites' stores point at."""
    lines = []
    for i, (gid, chrom, start, end) in enumerate(genes):
        attrs = f'gene_id "{gid}.1"; gene_type "protein_coding"; gene_name "SYM{i + 1}";'
        lines += [gtf_line(chrom, "gene", start, end, "+", attrs),
                  gtf_line(chrom, "exon", start, end, "+", attrs + f' transcript_id "T{i + 1}.1"; exon_number 1;')]
    st.root.mkdir(parents=True, exist_ok=True)
    return an.build(st, annot_id, write_gtf(st.root / f"{annot_id}.gtf.gz", lines))


@case
def parse_follows_the_v0_rules():
    with tempfile.TemporaryDirectory() as d:
        t = an.parse_gtf(write_gtf(Path(d) / "x.gtf.gz"))
    genes = t["genes"].to_pylist()
    assert [g["gene_id"] for g in genes] == ["ENSG00000000001", "ENSG00000000002"], genes
    assert genes[0] == {"gene_id": "ENSG00000000001", "version": 5, "name": "AAA", "biotype": "protein_coding",
                        "chr": "chr1", "tss": 100, "strand": "+", "start": 100, "end": 900}
    # minus strand: the TSS is the gene's end, the same rule steps_gtf uses
    assert genes[1]["tss"] == 1500 and genes[1]["strand"] == "-" and genes[1]["biotype"] == "lncRNA"
    ex = t["exons"].to_pylist()
    assert len(ex) == 4, ex                                    # the PAR_Y exon is dropped with its gene
    assert [e["transcript_id"] for e in ex] == ["ENST1.1", "ENST2.1", "ENST1.1", "ENST3.2"]
    assert ex[0] == {"gene_id": "ENSG00000000001", "transcript_id": "ENST1.1", "exon_number": 1,
                     "chr": "chr1", "start": 100, "end": 200, "strand": "+"}
    assert t["genes"].schema == an.GENE_SCHEMA and t["exons"].schema == an.EXON_SCHEMA


@case
def parse_rejects_an_unversioned_id():
    with tempfile.TemporaryDirectory() as d:
        bad = [LINES[0], gtf_line("chr1", "gene", 1, 2, "+", 'gene_id "ENSG1"; gene_type "x"; gene_name "y";')]
        raises("integer version", an.parse_gtf, write_gtf(Path(d) / "x.gtf.gz", bad))


@case
def identity_is_content_not_file_order():
    """Same rows in any GTF order give one identity; a moved coordinate gives another. That is what
    makes the digest reproducible across two machines re-downloading the same release."""
    with tempfile.TemporaryDirectory() as d:
        a = an.identity_digest(an.parse_gtf(write_gtf(Path(d) / "a.gtf.gz")))
        body = LINES[1:]
        random.Random(7).shuffle(body)
        shuffled = an.identity_digest(an.parse_gtf(write_gtf(Path(d) / "b.gtf.gz", ["##other header"] + body)))
        assert a == shuffled, "identity depends on GTF line order"
        moved = [ln.replace("\t900\t", "\t901\t") for ln in LINES]
        assert a != an.identity_digest(an.parse_gtf(write_gtf(Path(d) / "c.gtf.gz", moved)))
        renamed = [ln.replace('gene_name "AAA"', 'gene_name "ZZZ"') for ln in LINES]
        assert a != an.identity_digest(an.parse_gtf(write_gtf(Path(d) / "d.gtf.gz", renamed)))


@case
def arrow_zst_round_trip():
    with tempfile.TemporaryDirectory() as d:
        t = an.parse_gtf(write_gtf(Path(d) / "x.gtf.gz"))
    for k in ("genes", "exons"):
        blob = an.encode(t[k])
        assert blob[:4] == b"\x28\xb5\x2f\xfd", "not a zstd frame"      # zstd magic
        assert an.decode(blob).equals(t[k])
    raises("not a zstd frame", an.decode, b"not zstd at all", "annotation")


@case
def build_writes_objects_then_pointer():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        gtf = write_gtf(d / "gencode.v34.annotation.gtf.gz")
        st = qs.Store(d / "store")
        doc = an.build(st, "gencode_v34", gtf, {"name": "gencode", "version": "v34", "md5": "0" * 32})
        assert doc["n_genes"] == 2 and doc["n_exons"] == 4 and doc["n_transcripts"] == 3
        assert doc["identity_digest"] == an.identity_digest(an.parse_gtf(gtf))
        st.write_store("t", [])
        assert st.validate() == []
        got = an.load(st, "gencode_v34")
        assert got["genes"].equals(an.parse_gtf(gtf)["genes"]) and got["exons"].num_rows == 4
        assert got["doc"]["source"]["md5"] == "0" * 32       # provenance rides in the pointer, not the identity


@case
def same_gtf_gives_the_same_digests():
    """Reproducibility, the reason the object exists: two builds of one release must be the same
    annotation, down to the object names and the pointer bytes."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        src = {"name": "gencode", "version": "v34"}
        a = an.build(qs.Store(d / "s1"), "gencode_v34", write_gtf(d / "a.gtf.gz"), src)
        # a second copy of the same release: different path, different file name, shuffled lines
        body = LINES[1:]
        random.Random(3).shuffle(body)
        b = an.build(qs.Store(d / "s2"), "gencode_v34", write_gtf(d / "b.gtf.gz", ["##date: later"] + body), src)
        assert a == b, (a, b)
        assert (d / "s1" / "annotations" / "gencode_v34.json").read_bytes() == \
               (d / "s2" / "annotations" / "gencode_v34.json").read_bytes()


# (key as given, key normalized, bucket of 1024): the same vectors are checked in ui/scripts/store-check.ts,
# so the Python builder and the browser reader cannot hash a key into different buckets
LOOKUP_VECTORS = [("", "", 453), ("a", "A", 716), ("FLNC", "FLNC", 208), ("flnc", "FLNC", 208),
                  ("ENSG00000128591", "ENSG00000128591", 414), ("HLA-DRB1", "HLA-DRB1", 280),
                  ("Y_RNA", "Y_RNA", 828), ("Ä", "Ä", 834)]


@case
def lookup_key_and_bucket_are_pinned():
    for raw, key, bucket in LOOKUP_VECTORS:
        assert an.lookup_key(raw) == key, (raw, an.lookup_key(raw))
        assert an.lookup_bucket(key) == bucket, (key, an.lookup_bucket(key))
    assert an.lookup_bucket("A", 2 ** 32) == 0xC40BF6CC           # FNV-1a 32 of "A", the published value


@case
def chrom_tables_hold_each_chromosome_and_the_union_of_exons():
    with tempfile.TemporaryDirectory() as d:
        t = an.parse_gtf(write_gtf(Path(d) / "x.gtf.gz"))
    m = an.chrom_tables(t)
    assert list(m) == ["chr1", "chr2"]
    assert all(x["genes"].schema == an.GENE_SCHEMA and x["exon_models"].schema == an.EXON_MODEL_SCHEMA for x in m.values())
    # the genes object is the genes table's rows, unchanged
    assert m["chr1"]["genes"].equals(t["genes"].slice(0, 1)) and m["chr2"]["genes"].equals(t["genes"].slice(1, 1))
    # ENST1.1 [100,200] and ENST2.1 [100,250] overlap and merge; [800,900] stays apart
    assert m["chr1"]["exon_models"].to_pylist() == [{"gene_id": "ENSG00000000001", "exon_starts": [100, 800],
                                                     "exon_ends": [250, 900]}]
    assert m["chr2"]["exon_models"].to_pylist()[0]["exon_starts"] == [1400]
    # touching exons (start = previous end + 1) stay two exons; overlapping ones by one base merge
    ex = an.exon_models(pa.Table.from_pylist(
        [{"gene_id": "G", "transcript_id": "T", "exon_number": i, "chr": "chr1", "start": s, "end": e, "strand": "+"}
         for i, (s, e) in enumerate([(10, 20), (21, 30), (30, 40), (5, 8)])], schema=an.EXON_SCHEMA))
    assert ex == {"G": ([5, 10, 21], [8, 20, 40])}, ex


@case
def lookup_finds_ids_and_symbols_by_range():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st = qs.Store(d / "s")
        doc = an.build(st, "gencode_v34", write_gtf(d / "x.gtf.gz"))
        buf = (st.immutable / doc["lookup"]).read_bytes()
    h, buckets = an.decode_lookup(buf)
    assert (h["kind"], h["chrom"], h["count"], h["seq_digest"]) == (10, "all", an.LOOKUP_BUCKETS, doc["identity_digest"])
    assert sum(b.num_rows for b in buckets) == 4                    # two genes, each by id and by symbol
    for key in ("ENSG00000000002", "BBB"):
        # what a reader does: one u32 pair from the offset table, then that byte range alone
        b = an.lookup_bucket(key)
        lo, hi = struct.unpack_from("<2I", buf, qs.HEADER_LEN + 4 * b)
        rows = [r for r in an.decode(buf[lo:hi]).to_pylist() if r["key"] == key]
        assert rows == [{"key": key, "gene_id": "ENSG00000000002", "name": "BBB", "chr": "chr2", "tss": 1500}], rows
    raises("tile", an.decode_lookup, buf[:-1])


@case
def add_split_matches_build_and_validate_catches_a_bad_part():
    """An annotation built before the split objects gets the same pointer from `add_split` as a fresh
    build, and validate fails a gene models object that does not match the tables."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        gtf = write_gtf(d / "x.gtf.gz")
        st = qs.Store(d / "s")
        doc = an.build(st, "gencode_v34", gtf)
        old = {k: v for k, v in doc.items() if k not in an.SPLIT_KEYS}
        st.write_pointer("annotations", "gencode_v34", old)
        st.write_store("t", [])
        assert any("no chroms or lookup" in f for f in st.validate())
        assert an.add_split(st, "gencode_v34") == doc
        assert st.validate() == []
        t = an.parse_gtf(gtf)
        chr1 = an.chrom_tables(t)["chr1"]["exon_models"]
        bad = dict(doc, chroms=dict(doc["chroms"], chr2=dict(doc["chroms"]["chr2"], exon_models=st.put(an.encode(chr1), an.EXT))))
        st.write_pointer("annotations", "gencode_v34", bad)
        assert any("chr2 exon_models: rows differ" in f for f in st.validate()), st.validate()


@case
def gtf_source_reads_sources_yaml():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "sources.yaml").write_text(
            "dest: data/raw\nsources:\n"
            "  - name: gencode\n    version: v34 (GRCh38.p13)\n    dir: gencode_v34\n    files:\n"
            "      - url: https://example/gencode.v34.annotation.gtf.gz\n        size: 43164654\n"
            "        md5: 53901912df1002aae5a173cb505399d7\n")
        gtf = d / "gencode_v34" / "gencode.v34.annotation.gtf.gz"
        s = an.gtf_source(d / "sources.yaml", gtf)
        assert s == {"file": gtf.name, "name": "gencode", "version": "v34 (GRCh38.p13)",
                     "url": "https://example/gencode.v34.annotation.gtf.gz", "size": 43164654,
                     "md5": "53901912df1002aae5a173cb505399d7"}, s
        assert an.gtf_source(d / "nope.yaml", gtf) == {"file": gtf.name}      # degrades, does not fail
        # a second release in the same directory does not inherit the entry's release label
        v39 = d / "gencode_v34" / "gencode.v39.annotation.gtf.gz"
        v39.parent.mkdir()
        v39.write_bytes(b"x")
        assert an.gtf_source(d / "sources.yaml", v39) == {"file": v39.name, "md5": "9dd4e461268c8034f5c8564e155c67a6",
                                                          "size": 1}


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
