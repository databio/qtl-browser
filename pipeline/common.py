"""Shared helpers: config, paths, DuckDB connections, parquet writing."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml

ROOT = Path(__file__).resolve().parents[1]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# A full build is tens of GB and most of an hour, which is a bad loop to debug in: every breakage
# costs a full re-run to find the next one. These two variables cut a smoke build down to a couple
# of chromosomes so the whole pipeline can be exercised end to end in minutes.
#
#   QTLB_CHROMS=chr21,chr22 QTLB_DERIVED=/scratch/$USER/qtl-browser/derived-smoke \
#       uv run python -m pipeline build
#
# QTLB_DERIVED is what keeps a smoke run from colliding with a real one: the step markers, the
# tables and the packs all hang off `derived`, so pointing it elsewhere gives the smoke build its
# own everything. Never run a subset build into the real derived directory.
#
# A subset build is for finding breakage, not for release. Two things are deliberately wrong in it:
# `validate` compares egene and sQTL counts against `paper_counts`, which only a genome-wide build
# can meet; and chromosome ordinals in the variant and rsID indexes are positions in CHROMS, so a
# subset's packs are not byte-comparable to a full release.
if os.environ.get("QTLB_CHROMS"):
    CHROMS = [c.strip() for c in os.environ["QTLB_CHROMS"].split(",") if c.strip()]
    unknown = [c for c in CHROMS if not re.fullmatch(r"chr(\d{1,2}|X|Y|M)", c)]
    if unknown:
        sys.exit(f"QTLB_CHROMS: {unknown} are not chromosome names")


def log(msg: str) -> None:
    print(f"[{dt.datetime.now():%H:%M:%S}] {msg}", flush=True)


class Config:
    def __init__(self, path: Path | None = None):
        path = path or ROOT / "pipeline" / "config.yaml"
        self.cfg = yaml.safe_load(path.read_text())
        self.raw = ROOT / self.cfg["raw"]
        # QTLB_DERIVED redirects every output (see CHROMS above); absolute or relative to ROOT
        self.derived = Path(os.environ["QTLB_DERIVED"]) if os.environ.get("QTLB_DERIVED") else ROOT / self.cfg["derived"]
        self.zenodo = self.raw / self.cfg["zenodo_dir"]
        self.gtf = self.raw / self.cfg["gencode_gtf"]
        self.dbsnp_vcf = self.raw / self.cfg["dbsnp_vcf"]
        self.assembly_report = self.raw / self.cfg["dbsnp_assembly_report"]
        self.tmp = self.derived / "_tmp"
        # build intermediates: everything the pipeline reads but the browser never does. `_tables/`
        # is never uploaded (pipeline/README.md), so a table here costs local disk and nothing else.
        self.tables = self.derived / "_tables"
        # every file the browser reads at a byte offset, named by its own content hash
        self.immutable = self.derived / self.cfg["r2"]["immutable_prefix"].rstrip("/")
        self.done_dir = self.derived / ".done"

    def __getitem__(self, k):
        return self.cfg[k]

    def raw_dir(self, name: str) -> Path:
        """Extracted per-chromosome parquet dir for one Zenodo archive."""
        return self.zenodo / name

    def raw_glob(self, name: str) -> str:
        return str(self.raw_dir(name) / "*.parquet")

    def mark_done(self, step: str) -> None:
        self.done_dir.mkdir(parents=True, exist_ok=True)
        (self.done_dir / step).write_text(dt.datetime.now().isoformat())

    def is_done(self, step: str) -> bool:
        return (self.done_dir / step).exists()


def connect(cfg: Config, memory_limit: str | None = None, threads: int | None = None,
            temp_dir: Path | None = None) -> duckdb.DuckDBPyConnection:
    """Open an in-memory DuckDB. Concurrent processes MUST pass their own `temp_dir`: DuckDB
    spill files have fixed names, so two processes sharing one temp directory corrupt each other."""
    con = duckdb.connect()
    temp_dir = temp_dir or cfg.tmp
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = '{temp_dir}'")
    con.execute("SET preserve_insertion_order = false")
    con.execute("SET enable_progress_bar = false")
    if memory_limit:
        con.execute(f"SET memory_limit = '{memory_limit}'")
    if threads:
        con.execute(f"SET threads = {threads}")
    return con


def phenotype_batches(pfile: pq.ParquetFile, cols: list[str], batch_size: int = 1_000_000):
    """Tables of whole phenotypes from a file whose rows are grouped by phenotype_id. A table can be
    empty when one phenotype spans a whole batch; its rows arrive with the next table."""
    carry: pa.Table | None = None
    for batch in pfile.iter_batches(batch_size=batch_size, columns=cols):
        tb = pa.Table.from_batches([batch])
        if carry is not None and carry.num_rows:
            tb = pa.concat_tables([carry, tb])
        tb = tb.combine_chunks()
        ids = tb["phenotype_id"].combine_chunks()
        n = tb.num_rows
        ch = pc.indices_nonzero(pc.not_equal(ids.slice(0, n - 1), ids.slice(1))).to_numpy()
        last = int(ch[-1]) + 1 if len(ch) else 0
        yield tb.slice(0, last)
        carry = tb.slice(last)
    if carry is not None and carry.num_rows:
        yield carry.combine_chunks()


def strip_metadata(table: pa.Table) -> pa.Table:
    return table.replace_schema_metadata(None)


# ---- content-addressed publishing (SPEC.md section 3) ------------------------------------------
# A logical key uses letters, digits, "_" and "/" only, so the "." in a published name is
# unambiguous: <stem>.<sha16>.<ext>, stem being the key with "/" replaced by ".".
NAME = re.compile(r"^(?P<stem>[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*)\.(?P<sha>[0-9a-f]{16})\.(?P<ext>[a-z0-9.]+)$")
# Arrow schema metadata key in search_index: logical key -> sha256 of every pack its offsets reach
PACKS_METADATA_KEY = b"qtl_browser.packs"


def digests(path: Path) -> tuple[str, str]:
    """(sha256 hex, md5 hex) in one pass. MD5 is R2's ETag for a single-part upload."""
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as f:
        while chunk := f.read(8 << 20):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest()


def publish_file(cfg, tmp: Path, key: str, ext: str) -> Path:
    """Move a finished file to `immutable/<stem>.<sha16>.<ext>` and delete older builds of the same
    key. Changing one byte changes the name, so a browser can cache it for a year and a stale pack
    can never be paired with a fresh index. Write `tmp` under `cfg.tmp` so the rename stays on one
    filesystem."""
    stem = key.replace("/", ".")
    sha, _ = digests(tmp)
    out = cfg.immutable / f"{stem}.{sha[:16]}.{ext}"
    out.parent.mkdir(parents=True, exist_ok=True)
    os.replace(tmp, out)
    for old in out.parent.glob(f"{stem}.*.{ext}"):
        m = NAME.match(old.name)
        if old != out and m and m["stem"] == stem:
            old.unlink()
    return out


def stage(cfg, key: str, ext: str) -> Path:
    """Where a writer builds a file before `publish_file` names it by its hash."""
    out = cfg.tmp / "publish" / f"{key.replace('/', '.')}.{ext}"
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def addressed_files(cfg) -> dict[str, Path]:
    """Logical key -> its one current file. Fails on a stray file or two builds of one key."""
    out: dict[str, Path] = {}
    for p in sorted(cfg.immutable.glob("*")):
        m = NAME.match(p.name)
        if not m:
            die(f"unexpected file in {cfg.immutable}: {p.name}")
        key = m["stem"].replace(".", "/")
        if key in out:
            die(f"two builds of {key}: {out[key].name}, {p.name}")
        out[key] = p
    return out


def published(cfg, key: str) -> Path:
    """One published file by its logical key, with a message naming the step that writes it."""
    files = addressed_files(cfg)
    if key not in files:
        die(f"{key} is not published in {cfg.immutable.name}/; run the step that writes it")
    return files[key]


UNBUILT_SHA = "0" * 16


def pack_file(cfg, key: str, ext: str) -> Path:
    """The published path of `key`, or the name it would carry if the step that writes it had run.
    That placeholder does not exist, so a caller can say "missing, run step X" (or report a FAIL)
    instead of dying the way `published` does."""
    f = addressed_files(cfg).get(key)
    return f if f is not None else cfg.immutable / f"{key.replace('/', '.')}.{UNBUILT_SHA}.{ext}"


# the search index is not a pack, but it is content-addressed like one (SPEC section 3)
SEARCH_INDEX_EXT = "arrow.zst"
SEARCH_INDEX_NAME = f"search_index.{SEARCH_INDEX_EXT}"       # the logical name, without its hash


def search_index_path(cfg) -> Path:
    """The browser's one non-pack file (SPEC section 6): an Arrow IPC stream in a single zstd
    frame. Arrow rather than parquet so the browser needs no parquet reader at all."""
    return pack_file(cfg, "search_index", SEARCH_INDEX_EXT)


def read_search_index(cfg, columns: list[str] | None = None) -> pa.Table:
    from . import packfmt_v0 as packfmt
    path = search_index_path(cfg)
    raw = packfmt.zstd_unframe(path.read_bytes(), None, path.name)
    t = pa.ipc.open_stream(pa.py_buffer(raw)).read_all()
    return t.select(columns) if columns else t


def search_index_packs(cfg) -> dict[str, str]:
    """`qtl_browser.packs` from the search index's Arrow schema metadata: the logical key of every
    pack its byte offsets reach, to that file's full SHA-256 (SPEC section 3). Empty when the index
    predates the metadata."""
    from . import packfmt_v0 as packfmt
    path = search_index_path(cfg)
    raw = packfmt.zstd_unframe(path.read_bytes(), None, path.name)
    md = pa.ipc.open_stream(pa.py_buffer(raw)).schema.metadata or {}
    blob = md.get(PACKS_METADATA_KEY)
    return json.loads(blob) if blob else {}


def register_search_index(cfg, con, name: str = "search_index") -> str:
    """Make the search index queryable as `name` on a DuckDB connection."""
    con.register(name, read_search_index(cfg))
    return name


def variants_path(cfg) -> Path:
    """`_tables/variants.parquet`: every tested variant, sorted by (chr, position, A1, A2)."""
    return cfg.tables / "variants.parquet"


def variants_sql(cfg, chrom: str | None = None) -> str:
    """The variant table as a SQL source, optionally narrowed to one chromosome. One file with
    (chr, position) row-group statistics replaced the per-chromosome tree the browser used to read."""
    src = f"'{variants_path(cfg)}'"
    return src if chrom is None else f"(SELECT * FROM {src} WHERE chr = '{chrom}')"


def write_parquet(table: pa.Table, path: Path, row_group_size: int, stats_columns: list[str] | None = None,
                  metadata: dict[bytes, bytes] | None = None) -> None:
    """Plain parquet: zstd, dictionary, fixed row-group size."""
    path.parent.mkdir(parents=True, exist_ok=True)
    table = strip_metadata(table)
    if metadata:
        table = table.replace_schema_metadata(metadata)
    pq.write_table(
        table, path, compression="zstd", compression_level=9, use_dictionary=True,
        row_group_size=row_group_size, write_statistics=stats_columns if stats_columns else True,
    )
def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)
