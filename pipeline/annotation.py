"""Step 6: the annotation object -- one gene/transcript annotation per GENCODE/Ensembl release.

qtlb v0 fuses annotation into results. A gene's symbol, biotype, TSS, strand and window bounds live
inside every results file's gene-details JSON, inside `search_index`, and again inside
`genes.parquet`, `splice_phenotypes.parquet` and `credible_sets.parquet`. One experiment can carry
five copies of "where does ENSG00000128591 start"; two experiments built on two GENCODE releases
carry ten, and nothing in the store can tell you when two of them disagree.

This module builds the one object they should all point at instead:

    annotations/<annot_id>.json      mutable pointer: object digests + the release's identity
      genes  <digest>.arrow.zst      gene_id, version, name, biotype, chr, tss, strand, start, end
      exons  <digest>.arrow.zst      gene_id, transcript_id, exon_number, chr, start, end, strand

Experiments name an `annotation` id and join to it on `gene_id` (CONTRACT.md: contract tables carry
**no** annotation columns). The browser loads an annotation once per store, not once per experiment.

Encoding
--------
Arrow IPC **stream**, uncompressed, wrapped in exactly one zstd frame -- the same `arrow.zst` shape
`search_index` already uses, and for the same reason: the Arrow JavaScript reader cannot decode
IPC-internal buffer compression, while the browser already carries a zstd decoder for the packs.
(The framing helper is `packfmt_v1.zstd_frame`, the v1 codec.)

Identity
--------
`identity_digest` is sha512t24u over a canonical text of every stored row, gene lines then exon
lines, in the table's own (sorted) order:

    G\\t<gene_id>\\t<version>\\t<name>\\t<biotype>\\t<chr>\\t<strand>\\t<start>\\t<end>\\t<tss>\\n
    X\\t<gene_id>\\t<transcript_id>\\t<exon_number>\\t<chr>\\t<strand>\\t<start>\\t<end>\\n

**Included: every value the object stores**, derived ones too. `tss` is a function of
start/end/strand, so it is redundant -- deliberately. The whole point of the object is that two
experiments cannot disagree about where a gene starts, and a future change to the TSS rule is
exactly such a disagreement; it must surface as a different annotation, not as the same identity
with different bytes.

**Excluded, and why:**

- *Encoding choices* -- Arrow schema metadata, record-batch layout, zstd level, pyarrow version.
  Two machines re-encoding the same release must agree that it is the same release. This is why the
  identity is not just the object digest: content-addressed bytes are reproducible only within one
  environment, the canonical text is reproducible across them.
- *The GTF file's own bytes* -- its name, its path, its mtime, and its `##date` comment header,
  which changes on a re-download of an unchanged release. The source file's md5 is recorded in the
  pointer as provenance instead; provenance answers "where did this come from", identity answers
  "is this the same annotation", and they are not the same question.
- *The pointer's `id`* -- calling this release `gencode_v34` or `gencode34` does not change the
  annotation. An id collision between two different releases is caught by the identity, not hidden
  by it.
- *GTF fields this object does not store* -- CDS/UTR features, `tag`, `level`, `havana_*`,
  transcript support levels. By construction the identity covers exactly what is stored, so it can
  never claim more than the object delivers. Adding a column later changes the identity, correctly.

The pointer document is a pure function of (GTF content, id, source metadata): same GTF in, byte
identical `annotations/<id>.json` out. No build timestamp, on purpose.

What the browser reads
----------------------
The whole-genome tables are the canonical annotation, but a gene page needs one gene: 11 MB of
exons and 1.2 MB of genes before the first plot is what made the v1 gene page slower than v0. So
the pointer also names objects derived from them (`split`, SPEC.md section 6), which the browser
reads instead:

    chroms   {chr: {genes, exon_models}}  per chromosome: its `genes` rows, and each of those genes'
                                          collapsed exon model (`gene_id`, `exon_starts`, `exon_ends`)
    lookup   <digest>.qgl                 gene_id and symbol -> gene_id, name, chr, tss, hashed into
                                          1024 buckets read by range (kind 10)

They are derived, so they stay out of the identity; `check_split` proves they match the tables.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import re
import struct
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from . import qtlstore as qs
from .packfmt_v1 import zstd_frame, zstd_unframe

# GTF attributes are `key "value";` for strings and bare `key value;` for numbers. `steps_gtf` reads
# only the quoted form, which is why every `exon_number` in v0's exons.parquet is 0: GENCODE writes
# `exon_number 1;` unquoted, the regex misses it, and `a.get("exon_number", 0)` supplies the
# default. Both forms are read here. Quoted values are extracted exactly as steps_gtf extracts them,
# so symbol and biotype cannot differ for reasons that have nothing to do with this design.
ATTR = re.compile(r'(\S+)\s+(?:"([^"]*)"|([^";]*?))\s*;')
EXT = "arrow.zst"
ZSTD_LEVEL = 19                        # matches config.yaml packs.zstd_level

GENE_SCHEMA = pa.schema([("gene_id", pa.string()), ("version", pa.int32()), ("name", pa.string()),
                         ("biotype", pa.string()), ("chr", pa.string()), ("tss", pa.int32()),
                         ("strand", pa.string()), ("start", pa.int32()), ("end", pa.int32())])
EXON_SCHEMA = pa.schema([("gene_id", pa.string()), ("transcript_id", pa.string()), ("exon_number", pa.int32()),
                         ("chr", pa.string()), ("start", pa.int32()), ("end", pa.int32()), ("strand", pa.string())])
# one chromosome's collapsed exon models, one row per gene in genes order (SPEC.md section 6)
EXON_MODEL_SCHEMA = pa.schema([("gene_id", pa.string()), ("exon_starts", pa.list_(pa.int32())),
                               ("exon_ends", pa.list_(pa.int32()))])
# one bucket of the gene lookup (kind 10)
LOOKUP_SCHEMA = pa.schema([("key", pa.string()), ("gene_id", pa.string()), ("name", pa.string()),
                           ("chr", pa.string()), ("tss", pa.int32())])
LOOKUP_EXT = "qgl"
LOOKUP_BUCKETS = 1024


# ---- parsing ----------------------------------------------------------------------------------
def parse_attrs(s: str) -> dict[str, str]:
    # findall gives "" for the branch that did not participate, so `or` picks the one that did
    return {k: (q or bare) for k, q, bare in ATTR.findall(s)}


def _split_version(gid: str) -> tuple[str, int]:
    """`ENSG00000128591.16` -> (`ENSG00000128591`, 16). A non-integer suffix is a hard error rather
    than a silent 0: the unversioned id is the join key every experiment uses, and a release whose
    ids are shaped differently needs a decision, not a default."""
    base, _, ver = gid.partition(".")
    if not ver.isdigit():
        raise ValueError(f"gene id {gid!r} has no integer version suffix")
    return base, int(ver)


def parse_gtf(path: Path) -> dict[str, pa.Table]:
    """A GENCODE/Ensembl GTF -> `{"genes": table, "exons": table}`.

    `_PAR_Y` records are dropped, as in `steps_gtf`: they are a second copy of a pseudoautosomal
    gene placed on chrY under the same unversioned `gene_id`, so keeping them would make `gene_id`
    non-unique and break every join in the store.

    Both tables are sorted to a total order so the identity does not depend on GTF line order:
    genes by (chr, start, end, gene_id), exons by (chr, gene_id, start, end, transcript_id,
    exon_number). Exon order also matches what the v0 pack builder wants -- one gene's exons in one
    tight run.
    """
    genes: list[tuple] = []
    exons: list[tuple] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            a = parse_attrs(f[8])
            gid_v = a["gene_id"]
            if gid_v.endswith("_PAR_Y"):
                continue
            gene_id, version = _split_version(gid_v)
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            if f[2] == "gene":
                tss = start if strand == "+" else end
                genes.append((chrom, start, end, gene_id, version, a.get("gene_name") or "",
                              a.get("gene_type") or "", strand, tss))
            else:
                exons.append((chrom, gene_id, start, end, a["transcript_id"], int(a.get("exon_number", 0)), strand))
    genes.sort()
    exons.sort()
    gt = pa.Table.from_arrays(_columns(genes, (3, 4, 5, 6, 0, 8, 7, 1, 2), GENE_SCHEMA), schema=GENE_SCHEMA)
    et = pa.Table.from_arrays(_columns(exons, (1, 4, 5, 0, 2, 3, 6), EXON_SCHEMA), schema=EXON_SCHEMA)
    ids = gt.column("gene_id").to_pylist()
    if len(set(ids)) != len(ids):
        dup = sorted({g for g in ids if ids.count(g) > 1})[:5]
        raise ValueError(f"{path.name}: gene_id is not unique ({len(ids) - len(set(ids))} duplicates, e.g. {dup})")
    return {"genes": gt, "exons": et}


def _columns(rows: list[tuple], order: tuple[int, ...], schema: pa.Schema) -> list[pa.Array]:
    """Transpose sort-ordered tuples into the schema's column order."""
    cols = list(zip(*rows)) if rows else [()] * len(order)
    return [pa.array(cols[i], type=schema.field(j).type) for j, i in enumerate(order)]


# ---- identity ---------------------------------------------------------------------------------
def identity_digest(tables: dict[str, pa.Table]) -> str:
    """sha512t24u over the canonical text described in the module docstring."""
    h = hashlib.sha512()
    g = tables["genes"]
    for gene_id, version, name, biotype, chrom, tss, strand, start, end in zip(
            *(g.column(c).to_pylist() for c in
              ("gene_id", "version", "name", "biotype", "chr", "tss", "strand", "start", "end"))):
        h.update(f"G\t{gene_id}\t{version}\t{name}\t{biotype}\t{chrom}\t{strand}\t{start}\t{end}\t{tss}\n"
                 .encode("utf-8"))
    e = tables["exons"]
    for gene_id, tx, num, chrom, start, end, strand in zip(
            *(e.column(c).to_pylist() for c in
              ("gene_id", "transcript_id", "exon_number", "chr", "start", "end", "strand"))):
        h.update(f"X\t{gene_id}\t{tx}\t{num}\t{chrom}\t{strand}\t{start}\t{end}\n".encode("utf-8"))
    return base64.urlsafe_b64encode(h.digest()[:24]).decode("ascii")


# ---- encoding ---------------------------------------------------------------------------------
def encode(table: pa.Table, level: int = ZSTD_LEVEL) -> bytes:
    sink = pa.BufferOutputStream()
    with pa.ipc.new_stream(sink, table.schema) as w:
        w.write_table(table.combine_chunks())
    return zstd_frame(sink.getvalue().to_pybytes(), level)


def decode(data: bytes, what: str = "annotation") -> pa.Table:
    return pa.ipc.open_stream(pa.py_buffer(zstd_unframe(data, None, what))).read_all()


# ---- what the browser reads: gene models per chromosome, the gene lookup ----------------------
def exon_models(exons: pa.Table) -> dict[str, tuple[list[int], list[int]]]:
    """Each gene's collapsed exon model: the union of its transcripts' exons as sorted intervals, an
    exon merged into the one before it when it starts at or before that one's end (touching exons,
    `start == end + 1`, stay apart). The gene track draws this model."""
    by_gene: dict[str, list[tuple[int, int]]] = {}
    for g, s, e in zip(*(exons.column(c).to_pylist() for c in ("gene_id", "start", "end"))):
        by_gene.setdefault(g, []).append((s, e))
    out: dict[str, tuple[list[int], list[int]]] = {}
    for g, iv in by_gene.items():
        s_out: list[int] = []
        e_out: list[int] = []
        for s, e in sorted(iv):
            if e_out and s <= e_out[-1]:
                e_out[-1] = max(e_out[-1], e)
            else:
                s_out.append(s)
                e_out.append(e)
        out[g] = (s_out, e_out)
    return out


def chrom_tables(tables: dict[str, pa.Table]) -> dict[str, dict[str, pa.Table]]:
    """Per chromosome, in the order chromosomes first appear in `genes`: `genes`, that chromosome's
    `genes` rows in table order, and `exon_models`, one row per one of those genes with its collapsed
    exon model (empty lists for a gene with no exon records). Two objects rather than one so a reader
    of gene rows alone (the search box, the gene list) does not download exon models."""
    g = tables["genes"]
    models = exon_models(tables["exons"])
    chrom = g.column("chr").to_pylist()
    ids = g.column("gene_id").to_pylist()
    out = {}
    for c in dict.fromkeys(chrom):
        rows = [i for i, x in enumerate(chrom) if x == c]
        m = [models.get(ids[i], ([], [])) for i in rows]
        out[c] = {"genes": g.take(pa.array(rows, pa.int64())),
                  "exon_models": pa.Table.from_arrays([pa.array([ids[i] for i in rows], pa.string()),
                                                       pa.array([x[0] for x in m], pa.list_(pa.int32())),
                                                       pa.array([x[1] for x in m], pa.list_(pa.int32()))],
                                                      schema=EXON_MODEL_SCHEMA)}
    return out


def lookup_key(s: str) -> str:
    """A lookup key: the string with ASCII a-z upper-cased and every other character unchanged (so
    JavaScript and Python agree without Unicode case rules)."""
    return s.encode("utf-8").upper().decode("utf-8")


def lookup_bucket(key: str, n_buckets: int = LOOKUP_BUCKETS) -> int:
    """FNV-1a 32 over the key's UTF-8 bytes, modulo the bucket count."""
    h = 0x811C9DC5
    for b in key.encode("utf-8"):
        h = ((h ^ b) * 0x01000193) & 0xFFFFFFFF
    return h % n_buckets


def lookup_buckets(genes: pa.Table, n_buckets: int = LOOKUP_BUCKETS) -> list[pa.Table]:
    """The gene lookup's rows by bucket: every gene under its `gene_id` key and, when it has one that
    differs, its `name` key. Rows within a bucket sorted by (key, gene_id)."""
    rows: list[list[tuple]] = [[] for _ in range(n_buckets)]
    for gid, name, chrom, tss in zip(*(genes.column(c).to_pylist() for c in ("gene_id", "name", "chr", "tss"))):
        keys = [lookup_key(gid)]
        if name and lookup_key(name) != keys[0]:
            keys.append(lookup_key(name))
        for k in keys:
            rows[lookup_bucket(k, n_buckets)].append((k, gid, name, chrom, tss))
    out = []
    for r in rows:
        r.sort(key=lambda x: (x[0], x[1]))
        cols = list(zip(*r)) if r else [()] * 5
        out.append(pa.Table.from_arrays([pa.array(cols[i], type=LOOKUP_SCHEMA.field(i).type) for i in range(5)],
                                        schema=LOOKUP_SCHEMA))
    return out


def encode_lookup(buckets: list[pa.Table], identity: str, level: int = ZSTD_LEVEL) -> bytes:
    """Kind 10: the v1 header (chromosome `all`, count = buckets, seq_digest = the annotation's
    identity digest), u32 byte offsets of each bucket plus the end, then each bucket as one
    `encode`d Arrow object (an empty bucket is zero bytes)."""
    n = len(buckets)
    frames = [encode(b, level) if b.num_rows else b"" for b in buckets]
    off = qs.HEADER_LEN + 4 * (n + 1)
    offs = []
    for f in frames:
        offs.append(off)
        off += len(f)
    offs.append(off)
    if off > qs.U32_MAX:
        raise ValueError("gene lookup: over 4 GiB")
    return (qs.file_header(qs.KIND_GENE_LOOKUP, qs.ALL, n, 0, 0, identity)
            + struct.pack(f"<{n + 1}I", *offs) + b"".join(frames))


def decode_lookup(buf: bytes, what: str = "gene lookup") -> tuple[dict, list[pa.Table]]:
    """(header, buckets) of a kind-10 object."""
    h = qs.parse_file_header(buf)
    if h["kind"] != qs.KIND_GENE_LOOKUP:
        raise ValueError(f"{what}: header kind {h['kind']}, expected {qs.KIND_GENE_LOOKUP}")
    n = h["count"]
    offs = struct.unpack_from(f"<{n + 1}I", buf, qs.HEADER_LEN)
    if offs[0] != qs.HEADER_LEN + 4 * (n + 1) or offs[-1] != len(buf) or any(b < a for a, b in zip(offs, offs[1:])):
        raise ValueError(f"{what}: bucket offsets do not tile the object")
    out = []
    for b in range(n):
        lo, hi = offs[b], offs[b + 1]
        out.append(decode(buf[lo:hi], f"{what} bucket {b}") if hi > lo else LOOKUP_SCHEMA.empty_table())
    return h, out


SPLIT_KEYS = ("chroms", "lookup")


def split(store: qs.Store, tables: dict[str, pa.Table], identity: str, level: int = ZSTD_LEVEL) -> dict:
    """Put the objects the browser reads and return their pointer keys: `chroms` (chromosome ->
    `{genes, exon_models}`) and `lookup`. `build` and `add_split` both come here, so a store migrated in
    place and a store built from scratch hold the same bytes."""
    return {"chroms": {c: {k: store.put(encode(t, level), EXT) for k, t in parts.items()}
                       for c, parts in chrom_tables(tables).items()},
            "lookup": store.put(encode_lookup(lookup_buckets(tables["genes"]), identity, level), LOOKUP_EXT)}


def add_split(store: qs.Store, annot_id: str, level: int = ZSTD_LEVEL) -> dict:
    """Give an annotation built before the split objects existed its `chroms` and `lookup`, from the
    tables it already names, and rewrite its pointer (objects first). Returns the new pointer."""
    a = load(store, annot_id)
    doc = {k: v for k, v in a["doc"].items() if k not in SPLIT_KEYS}
    parts = split(store, {"genes": a["genes"], "exons": a["exons"]}, doc["identity_digest"], level)
    out = {}
    for k, v in doc.items():            # the split keys go right after `exons`, as `build` writes them
        out[k] = v
        if k == "exons":
            out |= parts
    store.write_pointer("annotations", annot_id, out)
    return out


def check_split(store: qs.Store, doc: dict) -> list[str]:
    """The browser's objects agree with the canonical tables: every chromosome of `genes` has a genes
    object holding exactly its rows and an exon models object holding their exon models, and the
    lookup's header and rows are what `lookup_buckets` makes of `genes`. Failures (empty = pass)."""
    where = f"annotations/{doc.get('id')}"
    if not isinstance(doc.get("chroms"), dict) or not isinstance(doc.get("lookup"), str):
        return [f"{where}: no chroms or lookup (the browser reads these; `python -m pipeline.annotation add-split`)"]
    names = [doc["genes"], doc["exons"], doc["lookup"]] + [n for p in doc["chroms"].values() for n in p.values()]
    if not all(isinstance(n, str) and (store.immutable / n).exists() for n in names):
        return []                       # missing objects are reported by validate already
    genes = decode((store.immutable / doc["genes"]).read_bytes(), f"{where} genes")
    exons = decode((store.immutable / doc["exons"]).read_bytes(), f"{where} exons")
    fails = []
    want = chrom_tables({"genes": genes, "exons": exons})
    if set(want) != set(doc["chroms"]):
        fails.append(f"{where}: chroms {sorted(doc['chroms'])} != the genes table's chromosomes {sorted(want)}")
    for c, parts in doc["chroms"].items():
        if set(parts) != {"genes", "exon_models"}:
            fails.append(f"{where} chroms {c}: keys {sorted(parts)}, expected exon_models and genes")
            continue
        for k, n in parts.items():
            try:
                got = decode((store.immutable / n).read_bytes(), f"{where} {c} {k}")
            except Exception as e:  # noqa: BLE001 -- any decode failure is a validation failure
                fails.append(f"{where} {c} {k}: {e}")
                continue
            if c in want and not got.equals(want[c][k]):
                fails.append(f"{where} {c} {k}: rows differ from the genes and exons tables")
    try:
        h, buckets = decode_lookup((store.immutable / doc["lookup"]).read_bytes(), f"{where} lookup")
    except ValueError as e:
        return fails + [f"{where} lookup: {e}"]
    exp = lookup_buckets(genes, h["count"])
    bad = [b for b in range(h["count"]) if not buckets[b].equals(exp[b])]
    if bad:
        fails.append(f"{where} lookup: {len(bad)} buckets differ from the genes table, e.g. bucket {bad[0]}")
    return fails


# ---- build ------------------------------------------------------------------------------------
def build(store: qs.Store, annot_id: str, gtf: Path, source: dict | None = None,
          level: int = ZSTD_LEVEL) -> dict:
    """Parse `gtf`, put both tables and the browser's split objects in the store, write
    `annotations/<annot_id>.json`, return it.

    Objects first and the pointer last, which `Store.write_pointer` enforces anyway.
    """
    tables = parse_gtf(Path(gtf))
    identity = identity_digest(tables)
    doc = {
        "id": annot_id,
        "identity_digest": identity,
        "genes": store.put(encode(tables["genes"], level), EXT),
        "exons": store.put(encode(tables["exons"], level), EXT),
        **split(store, tables, identity, level),
        "n_genes": tables["genes"].num_rows,
        "n_transcripts": len(set(tables["exons"].column("transcript_id").to_pylist())),
        "n_exons": tables["exons"].num_rows,
        "source": source if source is not None else {"file": Path(gtf).name},
    }
    store.write_pointer("annotations", annot_id, doc)
    return doc


def load(store: qs.Store, annot_id: str) -> dict:
    """`{"doc": ..., "genes": table, "exons": table}` for a stored annotation."""
    doc = store.load("annotations", annot_id)
    out = {"doc": doc}
    for k in ("genes", "exons"):
        out[k] = decode((store.immutable / doc[k]).read_bytes(), f"annotations/{annot_id} {k}")
    return out


def gtf_source(sources_yaml: Path, gtf: Path) -> dict:
    """Provenance for the pointer, read from `data/raw/sources.yaml`: which release this is, where
    it came from and its checksum. Identity does not depend on any of it (see the module docstring),
    so a missing sources.yaml degrades to the file name rather than failing the build.

    A source entry describes the GTF only when one of its files *is* that GTF: `data/raw/gencode_v34`
    is a symlink to the whole annotations brick, so v39 sits in a directory named for v34 and must
    not inherit v34's release label. A GTF no entry lists gets its own md5 and size instead."""
    out = {"file": Path(gtf).name}
    if Path(sources_yaml).exists():
        import yaml
        doc = yaml.safe_load(Path(sources_yaml).read_text())
        for src in doc.get("sources", []):
            if src.get("dir") != Path(gtf).parent.name:
                continue
            for f in src.get("files", []):
                url = f.get("url") or (src.get("url_base") or "").format(file=f.get("file", ""))
                if Path(f.get("file") or url).name == Path(gtf).name:
                    out |= {"name": src.get("name"), "version": src.get("version")}
                    out |= {k: v for k, v in (("url", url), ("md5", f.get("md5")), ("size", f.get("size"))) if v}
                    return out
    if not Path(gtf).exists():
        return out
    data = Path(gtf).read_bytes()
    return out | {"md5": hashlib.md5(data).hexdigest(), "size": len(data)}


# ---- migration check --------------------------------------------------------------------------
# Proves the annotation object holds exactly what the v0 tables hold, which is what licenses step 7
# to delete the annotation columns from `search_index`, the results packs and the small tables.
# It dies with those tables.
def verify(gtf: Path, tables_dir: Path) -> int:
    """Compare a fresh parse against `gene_annotation.parquet`, `genes.parquet` and `exons.parquet`.
    Returns the number of mismatches; prints one line per comparison."""
    t = parse_gtf(Path(gtf))
    bad = 0
    mine = {r["gene_id"]: r for r in t["genes"].to_pylist()}

    for name, cols in (("gene_annotation.parquet", ("gene_id_version", "symbol", "chr", "start", "end", "strand",
                                                    "tss", "biotype")),
                       ("genes.parquet", ("symbol", "chr", "start", "end", "strand", "tss", "biotype"))):
        path = Path(tables_dir) / name
        if not path.exists():
            print(f"SKIP {name}: not present")
            continue
        old = pq.read_table(path, columns=["gene_id", *cols]).to_pylist()
        diffs: dict[str, int] = {}
        example: dict[str, tuple] = {}
        for r in old:
            m = mine.get(r["gene_id"])
            if m is None:
                diffs["missing gene"] = diffs.get("missing gene", 0) + 1
                example.setdefault("missing gene", (r["gene_id"], None, None))
                continue
            for c in cols:
                got = f"{m['gene_id']}.{m['version']}" if c == "gene_id_version" else m["name" if c == "symbol" else c]
                if got != r[c]:
                    diffs[c] = diffs.get(c, 0) + 1
                    example.setdefault(c, (r["gene_id"], r[c], got))
        extra = len(mine) - len(old)
        bad += sum(diffs.values()) + abs(extra)
        detail = "; ".join(f"{c}: {n} differ, e.g. {example[c]}" for c, n in sorted(diffs.items()))
        print(f"{'FAIL' if diffs or extra else 'OK  '} {name}: {len(old)} rows vs {len(mine)} genes"
              f"{f', {extra:+d} row difference' if extra else ''}{'; ' + detail if detail else ''}")

    path = Path(tables_dir) / "exons.parquet"
    if not path.exists():
        print("SKIP exons.parquet: not present")
        return bad
    # exon_number is compared apart from the rest: v0 reads only quoted GTF attributes, so every
    # exon_number in exons.parquet is the 0 default rather than the exon's ordinal (see ATTR).
    cols = ["gene_id", "transcript_id", "chr", "start", "end", "strand"]
    old = pq.read_table(path, columns=[*cols, "exon_number"])
    a = sorted(zip(*(old.column(c).to_pylist() for c in cols)))
    b = sorted(zip(*(t["exons"].column(c).to_pylist() for c in cols)))
    if a == b:
        print(f"OK   exons.parquet: {len(a)} rows identical as a multiset (exon_number excluded)")
    else:
        only_old, only_new = set(a) - set(b), set(b) - set(a)
        bad += len(only_old) + len(only_new)
        print(f"FAIL exons.parquet: {len(a)} old vs {len(b)} new; {len(only_old)} only old, {len(only_new)} only new;"
              f" e.g. {sorted(only_old)[:1] or sorted(only_new)[:1]}")
    zeros = old.column("exon_number").to_pylist().count(0)
    mine_zeros = t["exons"].column("exon_number").to_pylist().count(0)
    print(f"NOTE exons.parquet exon_number: {zeros} of {old.num_rows} are 0 in v0, {mine_zeros} here"
          f"{' (v0 never parsed the unquoted attribute)' if zeros == old.num_rows else ''}")
    return bad


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="parse a GTF into a store's annotation object")
    b.add_argument("--gtf", required=True, type=Path)
    b.add_argument("--store", required=True, type=Path, help="store root (immutable/ and annotations/ live here)")
    b.add_argument("--id", required=True, help="annotation id, e.g. gencode_v34")
    b.add_argument("--sources", type=Path, default=Path("data/raw/sources.yaml"))
    s = sub.add_parser("add-split", help="give a stored annotation its per-chromosome objects and lookup (rewrites its pointer)")
    s.add_argument("--store", required=True, type=Path)
    s.add_argument("--id", required=True)
    v = sub.add_parser("verify", help="compare a fresh parse against the v0 _tables (migration check)")
    v.add_argument("--gtf", required=True, type=Path)
    v.add_argument("--tables", required=True, type=Path)
    args = ap.parse_args(argv)
    if args.cmd == "add-split":
        doc = add_split(qs.Store(args.store), args.id)
        print(json.dumps({"lookup": doc["lookup"], "chroms": len(doc["chroms"])}))
        return 0
    if args.cmd == "verify":
        bad = verify(args.gtf, args.tables)
        print(f"{'MISMATCHES: ' + str(bad) if bad else 'all comparisons agree'}")
        return 1 if bad else 0
    doc = build(qs.Store(args.store), args.id, args.gtf, gtf_source(args.sources, args.gtf))
    print(json.dumps(doc, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
