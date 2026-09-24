"""Shared helpers: config, paths, DuckDB connections, parquet writing."""
from __future__ import annotations

import datetime as dt
import os
import re
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

ROOT = Path(__file__).resolve().parents[1]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# Two variables cut a run down to a few chromosomes in its own tree, so a change can be exercised
# end to end in minutes (adapter.sbatch and store.sbatch pass them through):
#
#   QTLB_CHROMS=chr21,chr22 QTLB_DERIVED=/scratch/$USER/qtl-browser/derived-smoke \
#       uv run python -m pipeline.adapters.topchef
#
# QTLB_DERIVED is what keeps a subset run from colliding with a real one: the step markers and the
# tables all hang off `derived`, so pointing it elsewhere gives the run its own everything. Never
# run a subset into the frozen v0 derived directory.
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
        # build intermediates and contract tables: nothing under `_tables/` is ever deployed
        self.tables = self.derived / "_tables"
        # the frozen v0 packs, each named by its own content hash; read only by the v0 comparisons
        self.immutable = self.derived / "immutable"
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


def strip_metadata(table: pa.Table) -> pa.Table:
    return table.replace_schema_metadata(None)


# ---- the frozen v0 build's content-addressed files (v0 SPEC section 3) ---------------------------
# A logical key uses letters, digits, "_" and "/" only, so the "." in a published name is
# unambiguous: <stem>.<sha16>.<ext>, stem being the key with "/" replaced by ".".
NAME = re.compile(r"^(?P<stem>[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*)\.(?P<sha>[0-9a-f]{16})\.(?P<ext>[a-z0-9.]+)$")


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


UNBUILT_SHA = "0" * 16


def pack_file(cfg, key: str, ext: str) -> Path:
    """The published path of `key`, or the name it would carry if the step that writes it had run.
    That placeholder does not exist, so a caller can report a FAIL instead of dying."""
    f = addressed_files(cfg).get(key)
    return f if f is not None else cfg.immutable / f"{key.replace('/', '.')}.{UNBUILT_SHA}.{ext}"


# the v0 search index is not a pack, but it is content-addressed like one (v0 SPEC section 3)
SEARCH_INDEX_EXT = "arrow.zst"


def search_index_path(cfg) -> Path:
    """The v0 search index (v0 SPEC section 6): an Arrow IPC stream in a single zstd frame."""
    return pack_file(cfg, "search_index", SEARCH_INDEX_EXT)


def read_search_index(cfg, columns: list[str] | None = None) -> pa.Table:
    from . import packfmt_v0 as packfmt
    path = search_index_path(cfg)
    raw = packfmt.zstd_unframe(path.read_bytes(), None, path.name)
    t = pa.ipc.open_stream(pa.py_buffer(raw)).read_all()
    return t.select(columns) if columns else t


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
