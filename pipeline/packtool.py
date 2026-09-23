"""Inspect qtlb pack files and convert rows to and from them, outside the pipeline.

    uv run python -m pipeline.packtool header       PACK
    uv run python -m pipeline.packtool blocks       PACK.qbe|.qbs [--gene-ids] [--limit N] [-o TABLE]
    uv run python -m pipeline.packtool block        PACK.qbe|.qbs (--gene ID | --intron ID --eqtl PACK.qbe | --index I | --off O --len L)
                                                    [--variants PACK.qbv] [--dof N | --manifest manifest.json] [--details | --cs | --header] [-o TABLE]
    uv run python -m pipeline.packtool variants     PACK.qbv [--vidx V --n N | --off O --len L | --section cis|trans | --pages] [-o TABLE]
    uv run python -m pipeline.packtool frames       PACK.qbt [--search-index search_index.parquet] [--limit N] [-o TABLE]
    uv run python -m pipeline.packtool trans        PACK.qbt (--gene ID --search-index SI | --off O --len L [--gene-id G --gene-version V])
                                                    [--dof-eqtl N --dof-sqtl N | --manifest M] [--header] [-o TABLE]
    uv run python -m pipeline.packtool gwas         PACK.qbg --index gwas_index.bin (--lo L --hi H | --block K) [-o TABLE]
    uv run python -m pipeline.packtool gwas-index   gwas_index.bin [-o TABLE]
    uv run python -m pipeline.packtool check        PACK [--variants PACK.qbv] [--index gwas_index.bin] [--search-index SI] [--dof N]
    uv run python -m pipeline.packtool pack-variants TABLE -o PACK.qbv --chrom chr21 [--trans-only TABLE] [--page-size 512] [--codec zstd] [--level 19]
    uv run python -m pipeline.packtool pack-results  TABLE --variants PACK.qbv -o PACK.qbe|.qbs [--details D.json] [--pointers TABLE] [--level 19]
    uv run python -m pipeline.packtool pack-trans    TABLE -o PACK.qbt --chrom chr21 [--pointers TABLE] [--level 19]
    uv run python -m pipeline.packtool pack-gwas     TABLE -o DIR [--block-rows 2048] [--level 19]

`TABLE` is a path ending in `.parquet`, `.arrow` / `.feather` (Arrow IPC file), `.tsv`, `.csv`, or
`.json` (a list of objects); `-` means TSV on stdout. Every command takes plain paths and tables:
nothing here reads `config.yaml` or assumes the `data/derived/` layout. All bytes go through
`packfmt`, the reference codec, so this file adds no second encoder or decoder: it reads tables,
finds the right bytes, calls `packfmt`, and writes tables. `PACKS.md` is the plain-language
overview of the files and `SPEC.md` the byte layout.

Degrees of freedom (`dof`) are needed to derive a slope from a result block, and `beta_se` and `r2`
from a trans frame. Give `--dof` (or `--dof-eqtl`/`--dof-sqtl`), or `--manifest manifest.json` (its
`packs.dof`); with neither, the tool looks for a `manifest.json` in the pack's parent directories,
which finds `data/derived/manifest.json` for pipeline builds.

Trans frames carry no gene id (SPEC section 12): a gene's frame is found through `search_index`
(`trans_off`, `trans_len`), so `trans --gene` needs `--search-index`. `frames` and `check` walk the
frames without it, using the zstd frame boundaries.

Adding a file kind (the variant-page files: hits pack, rsID index, variant index): register it in
`KINDS` below, add a read function that calls the new `packfmt` decoder and returns a table, a
write function that calls the new encoder, a `check_file` branch, and a subcommand in `main`.
`test_packtool.py` shows the round-trip test shape to copy.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.feather as pafeather
import pyarrow.parquet as pq
import zstandard

from . import packfmt_v0 as packfmt

# ---- file kinds ----------------------------------------------------------------------------------
# kind number -> (name, extension, one-line meaning). SPEC.md section 3 numbers the kinds. Extensions
# are a convention (SPEC section 3); the header's kind byte is what the tool trusts.
KINDS = {
    packfmt.KIND_VARIANTS: ("variants", ".qbv", "one chromosome's variants in pages: the cis section, then the trans-only section (SPEC section 4)"),
    packfmt.KIND_EQTL: ("eqtl", ".qbe", "one block per gene: eQTL rows, credible sets, gene details (SPEC section 5)"),
    packfmt.KIND_SQTL: ("sqtl", ".qbs", "one block per tested intron: sQTL rows and credible sets (SPEC section 5)"),
    packfmt.KIND_GWAS: ("gwas", ".qbg", "one chromosome's GWAS rows in zstd blocks (SPEC section 11)"),
    packfmt.KIND_GWAS_INDEX: ("gwas_index", ".bin", "block positions and offsets of every GWAS pack (SPEC section 11)"),
    packfmt.KIND_TRANS: ("trans", ".qbt", "one zstd frame per gene with its trans eQTL and sQTL rows (SPEC section 12)"),
    packfmt.KIND_HITS: ("hits", ".qbh", "one zstd frame per 1024 variants: each variant's trans, lead and credible-set rows (SPEC section 13)"),
    packfmt.KIND_RSID: ("rsid_index", ".qbr", "every rsID as a sorted 8-byte record, in uncompressed blocks (SPEC section 14)"),
    packfmt.KIND_VARIANT_INDEX: ("variant_index", ".qbx", "startup file: page, hits frame and rsID block offsets of every chromosome (SPEC section 15)"),
}
KIND_BY_NAME = {v[0]: k for k, v in KINDS.items()}
KIND_BY_EXT = {v[1]: k for k, v in KINDS.items() if v[1] != ".bin"}
DOF_KEY = {packfmt.KIND_EQTL: "eqtl", packfmt.KIND_SQTL: "sqtl"}    # manifest packs.dof key per results kind
PLACEHOLDER_DOF = 1000   # when only a block's or frame's structure is wanted: derived values are then meaningless and never returned
                         # (large enough that t(p) stays about 37 even at p = 1e-300, so nothing overflows float32)

TABLE_EXTS = (".parquet", ".arrow", ".feather", ".ipc", ".tsv", ".csv", ".json")
_ZSTD_CHUNK = 1 << 20


class PackToolError(ValueError):
    """A bad argument, a missing file, or a file that breaks a SPEC rule (message names which)."""


# ---- tables in and out ---------------------------------------------------------------------------
def read_table(path: str | os.PathLike) -> pa.Table:
    """Read a table by extension: `.parquet`; `.arrow`, `.feather`, `.ipc` (Arrow IPC file);
    `.tsv`, `.csv` (header row, types inferred, empty fields are null); `.json` (a list of objects)."""
    p = Path(path)
    if not p.is_file():
        raise PackToolError(f"table {p}: no such file")
    ext = p.suffix.lower()
    if ext == ".parquet":
        return pq.read_table(p)
    if ext in (".arrow", ".feather", ".ipc"):
        return pafeather.read_table(p)
    if ext in (".tsv", ".csv"):
        delim = "\t" if ext == ".tsv" else ","
        return pacsv.read_csv(p, parse_options=pacsv.ParseOptions(delimiter=delim),
                              convert_options=pacsv.ConvertOptions(strings_can_be_null=True))
    if ext == ".json":
        rows = json.loads(p.read_text())
        if not isinstance(rows, list):
            raise PackToolError(f"table {p}: JSON must be a list of objects")
        return pa.Table.from_pylist(rows)
    raise PackToolError(f"table {p}: unknown extension {ext!r}, expected one of {TABLE_EXTS}")


def write_table(table: pa.Table, path: str | os.PathLike | None) -> None:
    """Write a table by extension (see `read_table`); `None` or `-` writes TSV to stdout. Nulls are
    empty fields in TSV and CSV, and `null` in JSON."""
    if path is None or str(path) == "-":
        sys.stdout.write(_tsv_bytes(table).decode("utf-8"))
        return
    p = Path(path)
    ext = p.suffix.lower()
    p.parent.mkdir(parents=True, exist_ok=True)
    if ext == ".parquet":
        pq.write_table(table, p)
    elif ext in (".arrow", ".feather", ".ipc"):
        pafeather.write_feather(table, p, compression="uncompressed")
    elif ext == ".tsv":
        p.write_bytes(_tsv_bytes(table))
    elif ext == ".csv":
        pacsv.write_csv(table, p, write_options=pacsv.WriteOptions(delimiter=","))
    elif ext == ".json":
        p.write_text(json.dumps(_clean(table.to_pylist()), indent=1))
    else:
        raise PackToolError(f"table {p}: unknown extension {ext!r}, expected one of {TABLE_EXTS} or '-'")


def _tsv_bytes(table: pa.Table) -> bytes:
    """TSV without quotes (pack strings are ASCII alleles and ids, never tabs or newlines)."""
    sink = pa.BufferOutputStream()
    pacsv.write_csv(table, sink, write_options=pacsv.WriteOptions(delimiter="\t", quoting_style="none"))
    b = sink.getvalue().to_pybytes()
    head, _, rest = b.partition(b"\n")             # pyarrow quotes the header names even so
    return head.replace(b'"', b"") + b"\n" + rest


def _clean(x):
    """JSON-safe copy: numpy scalars and arrays to Python, NaN and infinities to null, bytes to str."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_clean(v) for v in x.tolist()]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        f = float(x)
        return f if math.isfinite(f) else None
    if isinstance(x, bytes):
        return x.decode("ascii", "replace")
    return x


def _masked_int(a, dtype=pa.int64()) -> pa.Array:
    data, mask = np.ma.getdata(a), np.ma.getmaskarray(a)
    return pa.array(np.where(mask, 0, data).astype(np.int64), mask=mask, type=dtype)


def _nan_float(a) -> pa.Array:
    a = np.asarray(a, dtype=np.float64)
    m = ~np.isfinite(a)
    return pa.array(np.where(m, 0.0, a), mask=m, type=pa.float64())


def _zero_null_int(a) -> pa.Array:
    a = np.asarray(a).astype(np.int64)
    return pa.array(a, mask=a == 0, type=pa.int64())


def _column(table: pa.Table, name: str, required: bool = True, default=None):
    """One column of `table`, or `default` (a scalar, repeated) when the column is missing and not required."""
    if name in table.column_names:
        return table[name].combine_chunks()
    if required:
        raise PackToolError(f"table is missing the column {name!r}; it has {table.column_names}")
    return pa.array([default] * table.num_rows)


def _floats(col) -> np.ndarray:
    """float64 with NaN for nulls, from any numeric or null-only column."""
    if col.null_count == len(col) and col.type == pa.null():
        return np.full(len(col), np.nan)
    return col.cast(pa.float64()).to_numpy(zero_copy_only=False)


def _ints(col, name: str) -> np.ndarray:
    a = _floats(col)
    if np.any(np.isnan(a)):
        raise PackToolError(f"{name}: null values are not allowed")
    if np.any(a != np.rint(a)):
        raise PackToolError(f"{name}: non-integer values")
    return a.astype(np.int64)


def _strings(col) -> list:
    """Python strings with None for nulls (a null-only column gives all None)."""
    if col.type == pa.null():
        return [None] * len(col)
    return col.cast(pa.string()).to_pylist()


# ---- files -----------------------------------------------------------------------------------------
def read_range(path: str | os.PathLike, off: int, length: int) -> bytes:
    """`length` bytes from byte `off` of the file, exactly (a short read is an error): the local
    stand-in for the browser's `Range: bytes=off-(off+length-1)` request."""
    p = Path(path)
    if not p.is_file():
        raise PackToolError(f"{p}: no such file")
    if off < 0 or length < 0:
        raise PackToolError(f"{p}: negative range {off}+{length}")
    with open(p, "rb") as f:
        f.seek(off)
        b = f.read(length)
    if len(b) != length:
        raise PackToolError(f"{p}: range {off}+{length} runs past the end of the file ({p.stat().st_size} bytes)")
    return b


def read_header(path: str | os.PathLike) -> dict:
    """The 32-byte file header as a dict (`packfmt.parse_file_header`: kind, version, chrom, count,
    page_size, n_cis) plus `path`, `bytes`, `kind_name`, and whether the extension matches the kind."""
    p = Path(path)
    h = packfmt.parse_file_header(read_range(p, 0, packfmt.FILE_HEADER_LEN))
    name, ext, _ = KINDS.get(h["kind"], ("?", "", ""))
    return {"path": str(p), "bytes": p.stat().st_size, "kind_name": name,
            "extension_matches": p.suffix == ext, **h}


def _expect_kind(h: dict, kinds: tuple[int, ...], path) -> None:
    if h["kind"] not in kinds:
        want = ", ".join(f"{k} ({KINDS[k][0]})" for k in kinds)
        raise PackToolError(f"{path}: header kind {h['kind']} ({KINDS.get(h['kind'], ('?',))[0]}), expected {want}")


def find_manifest(start: str | os.PathLike) -> Path | None:
    """The nearest `manifest.json` in `start` or its parents (`data/derived/manifest.json` for pipeline builds)."""
    p = Path(start).resolve()
    for d in [p, *p.parents]:
        m = d / "manifest.json"
        if m.is_file():
            return m
    return None


def dof_from_manifest(manifest: str | os.PathLike, key: str) -> int:
    """`packs.dof.eqtl` or `packs.dof.sqtl` from a pipeline `manifest.json`."""
    m = json.loads(Path(manifest).read_text())
    try:
        return int(m["packs"]["dof"][key])
    except KeyError as e:
        raise PackToolError(f"{manifest}: no packs.dof.{key} ({e})") from None


def _manifest_near(manifest, near) -> Path:
    if manifest is None and near is not None:
        manifest = find_manifest(Path(near).parent)
    if manifest is None:
        raise PackToolError("the derived values need the degrees of freedom: give --dof (or --dof-eqtl and --dof-sqtl) or --manifest manifest.json")
    return Path(manifest)


def resolve_dof(kind: int, dof: int | None, manifest: str | os.PathLike | None, near: str | os.PathLike | None) -> int:
    """`dof` if given, else the manifest's value for a results kind (2 eQTL, 3 sQTL), else the
    manifest nearest `near`; PackToolError otherwise."""
    if dof is not None:
        return int(dof)
    return dof_from_manifest(_manifest_near(manifest, near), DOF_KEY[kind])


def resolve_dofs(dof_e: int | None, dof_s: int | None, manifest: str | os.PathLike | None, near: str | os.PathLike | None) -> tuple[int, int]:
    """(eQTL dof, sQTL dof) for a trans frame: both given, else both from the manifest."""
    if dof_e is not None and dof_s is not None:
        return int(dof_e), int(dof_s)
    if dof_e is not None or dof_s is not None:
        raise PackToolError("give both --dof-eqtl and --dof-sqtl, or --manifest")
    m = _manifest_near(manifest, near)
    return dof_from_manifest(m, "eqtl"), dof_from_manifest(m, "sqtl")


# ---- variants file (kind 1) ------------------------------------------------------------------------
def page_table(path: str | os.PathLike) -> pa.Table:
    """Every page of a variants file from its page headers alone, without decompressing
    (`packfmt.walk_variants_file`): `page`, `off` (file byte offset), `len` (bytes including
    padding), `first_vidx`, `n` records, `codec`, `section` (`cis` for pages below `n_cis`, `trans`
    for the trans-only section). The pages are what a browser range request covers: a run of
    variants `[vidx, vidx + n)` is read as the bytes of the pages holding its first and last record
    (SPEC section 6, `var_off`/`var_len`)."""
    h, pages = packfmt.walk_variants_file(Path(path).read_bytes())
    n_cis = h["n_cis"]
    return pa.table({"page": pa.array(range(len(pages)), pa.int64()),
                     "off": pa.array([p["offset"] for p in pages], pa.int64()), "len": pa.array([p["length"] for p in pages], pa.int64()),
                     "first_vidx": pa.array([p["first_vidx"] for p in pages], pa.int64()), "n": pa.array([p["n"] for p in pages], pa.int64()),
                     "codec": pa.array([p["codec"] for p in pages], pa.string()),
                     "section": pa.array(["cis" if p["first_vidx"] < n_cis else "trans" for p in pages], pa.string())})


def covering_pages(pages: pa.Table, vidx: int, n: int) -> tuple[int, int]:
    """(off, len) of the pages that hold `[vidx, vidx + n)`: the bytes a reader fetches for that run."""
    if n < 1:
        raise PackToolError("a run needs n >= 1")
    first = pages["first_vidx"].to_numpy()
    cnt = pages["n"].to_numpy()
    total = int(first[-1] + cnt[-1]) if len(first) else 0
    if vidx < 0 or vidx + n > total:
        raise PackToolError(f"run {vidx}+{n} is outside the file's {total} records")
    a = int(np.searchsorted(first, vidx, side="right")) - 1
    b = int(np.searchsorted(first, vidx + n - 1, side="right")) - 1
    off = int(pages["off"][a].as_py())
    return off, int(pages["off"][b].as_py() + pages["len"][b].as_py()) - off


def _decode_run(path: str | os.PathLike, vidx: int, n: int, *, pos_first: int | None = None, pos_last: int | None = None) -> dict:
    """`packfmt.decode_variant_pages` over the pages covering a phenotype's run `[vidx, vidx + n)`,
    with SPEC section 7 step 4's checks (the run lies in the cis section, below `n_cis`)."""
    h = read_header(path)
    off, ln = covering_pages(page_table(path), vidx, n)
    return packfmt.decode_variant_pages(read_range(path, off, ln), var_start=vidx, n_var=n, pos_first=pos_first, pos_last=pos_last,
                                        n_cis=h["n_cis"])


def variants_table(decoded: dict, start: int = 0, stop: int | None = None) -> pa.Table:
    """Decoded pages (a `packfmt.decode_variant_pages` dict) as a table, rows `start:stop`: `vidx`,
    `position` (1-based GRCh38), `A1`, `A2` (null when the record has flags bit 2, alleles not
    reported; trans-only records only), `no_alleles`, `rs_number` (null when the file holds 0), `af`
    (A1 frequency, null for code 65535), `ma_samples`, `ma_count` (null for 65535; always null in the
    trans-only section), `match` (rsID match: none, exact, position), `allele_code` (0 = alleles
    from the heap or not reported, 1-12 = SNP table)."""
    sl = slice(start, stop)
    return pa.table({
        "vidx": pa.array(decoded["vidx"][sl].astype(np.int64)),
        "position": pa.array(decoded["position"][sl].astype(np.int64)),
        "A1": pa.array(decoded["A1"][sl], pa.string()), "A2": pa.array(decoded["A2"][sl], pa.string()),
        "no_alleles": pa.array(np.asarray(decoded["no_alleles"])[sl].astype(bool)),
        "rs_number": _masked_int(decoded["rs_number"][sl]),
        "af": _nan_float(decoded["af"][sl]),
        "ma_samples": _masked_int(decoded["ma_samples"][sl]), "ma_count": _masked_int(decoded["ma_count"][sl]),
        "match": pa.array(decoded["match"][sl], pa.string()),
        "allele_code": pa.array(decoded["allele_code"][sl].astype(np.int16)),
    })


def _empty_variants() -> pa.Table:
    return variants_table({"vidx": np.array([], np.uint32), "position": np.array([], np.uint32), "A1": [], "A2": [],
                           "no_alleles": np.array([], bool), "rs_number": np.array([], np.uint32), "af": np.array([], np.float64),
                           "ma_samples": np.array([], np.uint16), "ma_count": np.array([], np.uint16), "match": [], "allele_code": np.array([], np.uint8)})


def variant_rows(path: str | os.PathLike, *, vidx: int | None = None, n: int | None = None,
                 off: int | None = None, length: int | None = None, section: str | None = None,
                 whole_pages: bool = False) -> pa.Table:
    """Rows of a variants file (columns: `variants_table`). With `vidx` and `n`: the records
    `[vidx, vidx + n)` from either section (or, with `whole_pages`, every record of the pages covering
    them, which is what a browser decodes). With `off` and `length`: decode that byte range as pages,
    as a browser would from `search_index.var_off`/`var_len`. With `section` `cis` or `trans`: that
    section. With none of them: the whole file, checked as one (`packfmt.decode_variants_file`)."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_VARIANTS,), p)
    n_cis, count = h["n_cis"], h["count"]
    if section is not None:
        if section not in ("cis", "trans"):
            raise PackToolError(f"section must be cis or trans, not {section!r}")
        vidx, n = (0, n_cis) if section == "cis" else (n_cis, count - n_cis)
        if n == 0:
            return _empty_variants()
    if vidx is not None or n is not None:
        if vidx is None or n is None:
            raise PackToolError("give both --vidx and --n")
        off, ln = covering_pages(page_table(p), vidx, n)
        d = packfmt.decode_variant_pages(read_range(p, off, ln), n_cis=n_cis)
        first = int(d["vidx"][0])
        return variants_table(d) if whole_pages else variants_table(d, vidx - first, vidx - first + n)
    if off is not None or length is not None:
        if off is None or length is None:
            raise PackToolError("give both --off and --len")
        return variants_table(packfmt.decode_variant_pages(read_range(p, off, length), n_cis=n_cis))
    d = packfmt.decode_variants_file(p.read_bytes())
    return variants_table(d) if d["pages"] else _empty_variants()


def write_variants_pack(table: pa.Table, out: str | os.PathLike, chrom: str, *, page_size: int = 512,
                        codec: str = "zstd", level: int = 19, trans_only: pa.Table | None = None) -> pa.Table:
    """Rows to a kind 1 variants file. Cis columns: `position` (int, 1-based), `A1`, `A2` (strings;
    SNPs take a 1-byte code, anything else goes to the page's allele heap); optional `rs_number`
    (int, null or 0 = none), `af` (float in [0, 1], null allowed), `ma_samples`, `ma_count` (int
    0..65534, null allowed), `match` (none, exact, position; default none). Rows are sorted here by
    (position, A1, A2) byte-wise, which defines `vidx`; duplicate keys fail.

    `trans_only` is the trans-only section (SPEC section 4): `position` (one row per position),
    optional `A1`, `A2` (null together where the source reports no alleles, written with flags bit
    2), `rs_number`, `af`, `match`; sample counts are always null there. It is sorted by position
    here. When `trans_only` is not given and `table` has a boolean `in_cis` column (the pipeline's
    variant table shape), the false rows become the trans-only section. Writes to a temporary name
    and renames. Returns the page table of the new file."""
    if not isinstance(table, pa.Table):
        raise PackToolError("write_variants_pack: table must be a pyarrow.Table")
    if trans_only is None and "in_cis" in table.column_names:
        cis = _column(table, "in_cis").cast(pa.bool_()).to_numpy(zero_copy_only=False).astype(bool)
        table, trans_only = table.filter(pa.array(cis)), table.filter(pa.array(~cis))
    for c in ("position", "A1", "A2"):
        _column(table, c)
    t = table.sort_by([("position", "ascending"), ("A1", "ascending"), ("A2", "ascending")])
    if "match" not in t.column_names:
        t = t.append_column("match", pa.array(["none"] * t.num_rows, pa.string()))
    for c in ("rs_number", "af", "ma_samples", "ma_count"):
        if c not in t.column_names:
            t = t.append_column(c, pa.nulls(t.num_rows, pa.float64()))
    tr = None
    if trans_only is not None and trans_only.num_rows:
        _column(trans_only, "position")
        x = trans_only.sort_by([("position", "ascending")])
        a1 = _strings(_column(x, "A1", required=False)) if "A1" in x.column_names else [None] * x.num_rows
        a2 = _strings(_column(x, "A2", required=False)) if "A2" in x.column_names else [None] * x.num_rows
        tr = {"position": _ints(_column(x, "position"), "trans-only position"),
              "rs_number": _floats(_column(x, "rs_number", required=False)) if "rs_number" in x.column_names else np.full(x.num_rows, np.nan),
              "af": _floats(_column(x, "af", required=False)) if "af" in x.column_names else np.full(x.num_rows, np.nan),
              "A1": a1, "A2": a2,
              "match": [m or "none" for m in (_strings(_column(x, "match", required=False)) if "match" in x.column_names else [None] * x.num_rows)]}
    buf, offsets = packfmt.encode_variants_file(chrom, _column(t, "position"), _floats(_column(t, "rs_number")), _floats(_column(t, "af")),
                                                _floats(_column(t, "ma_samples")), _floats(_column(t, "ma_count")),
                                                _strings(_column(t, "A1")), _strings(_column(t, "A2")), _strings(_column(t, "match")),
                                                page_size, codec, level, trans_only=tr)
    _write_bytes(out, buf)
    return page_table(out)


def _write_bytes(out: str | os.PathLike, buf: bytes) -> None:
    p = Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_bytes(buf)
    os.replace(tmp, p)


# ---- results packs (kinds 2 and 3) -------------------------------------------------------------------
def list_blocks(path: str | os.PathLike, *, limit: int | None = None, gene_ids: bool = False) -> pa.Table:
    """Every block of an eQTL or sQTL pack, from the block headers (`packfmt.parse_block_header`):
    `block`, `blk_off`, and the SPEC section 5 header fields (`blk_len`, `n_rows`, `var_start`
    (null when `n_rows` is 0), `anchor`, `pos_first`, `pos_last`, `n_cs`, `nlp_max`, `lse_min`,
    `lse_max`, `details_zlen`, `details_len`). With `gene_ids` on a kind 2 pack, each block's
    details are decoded and `gene_id` and `symbol` added (slower: every block is decoded)."""
    p = Path(path)
    buf = p.read_bytes()
    h = packfmt.parse_file_header(buf)
    _expect_kind(h, packfmt.RESULT_KINDS, p)
    _, blocks = packfmt.walk_eqtl_file(buf, h["kind"])
    blocks = blocks[:limit] if limit is not None else blocks
    heads = [packfmt.parse_block_header(buf, off) for off, _ in blocks]
    vs = np.array([x["var_start"] for x in heads], dtype=np.int64)
    out = {"block": pa.array(range(len(blocks)), pa.int64()), "blk_off": pa.array([off for off, _ in blocks], pa.int64())}
    for k in packfmt.BLOCK_FIELDS:
        vals = [x[k] for x in heads]
        if k == "var_start":
            out[k] = pa.array(np.where(vs == packfmt.NO_VAR_START, 0, vs), mask=vs == packfmt.NO_VAR_START, type=pa.int64())
        elif k in ("nlp_max", "lse_min", "lse_max"):
            out[k] = pa.array(vals, pa.float64())
        else:
            out[k] = pa.array(vals, pa.int64())
    if gene_ids and h["kind"] == packfmt.KIND_EQTL:
        gids, syms = [], []
        for off, ln in blocks:
            g = packfmt.decode_gene_block(buf[off:off + ln], PLACEHOLDER_DOF, kind=h["kind"], expect_blk_len=ln)["details"]["gene"]
            gids.append(g.get("gene_id")); syms.append(g.get("symbol"))
        out = {"block": out["block"], "gene_id": pa.array(gids, pa.string()), "symbol": pa.array(syms, pa.string()), **{k: v for k, v in out.items() if k != "block"}}
    return pa.table(out)


def read_search_index_file(path: str | os.PathLike, columns: list[str] | None = None) -> pa.Table:
    """The pipeline's search index, as either the Arrow IPC stream in one zstd frame it writes
    today (`search_index.arrow.zst`) or a plain parquet file. The tool takes a path, so it must
    accept whatever the caller names."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pq.read_table(path, columns=columns)
    raw = packfmt.zstd_unframe(path.read_bytes(), None, str(path))
    t = pa.ipc.open_stream(pa.py_buffer(raw)).read_all()
    return t.select(columns) if columns else t


def _search_index_rows(search_index: str | os.PathLike, gene: str, columns: list[str]) -> list[dict]:
    si = read_search_index_file(search_index, ["gene_id", "symbol", "chr"] + columns)
    rows = [r for r in si.to_pylist() if r["gene_id"] == gene or r["symbol"] == gene]
    if not rows:
        raise PackToolError(f"{search_index}: no gene {gene!r}")
    return rows


def find_gene_block(path: str | os.PathLike, gene_id: str, search_index: str | os.PathLike | None = None) -> tuple[int, int]:
    """(blk_off, blk_len) of a gene's block in a kind 2 pack. With `search_index` (the pipeline's
    parquet), a direct lookup of `blk_off`/`blk_len`; otherwise every block's details are decoded
    until `gene.gene_id` matches. Also accepts a gene symbol."""
    p = Path(path)
    if search_index is not None:
        chrom = read_header(p)["chrom"]
        rows = _search_index_rows(search_index, gene_id, ["blk_off", "blk_len"])
        rows = [r for r in rows if r["chr"] == chrom] or rows
        r = rows[0]
        if r["chr"] != chrom:
            raise PackToolError(f"{gene_id} is on {r['chr']}, but {p} holds {chrom}")
        if r["blk_off"] is None:
            raise PackToolError(f"{gene_id}: search_index says this gene has no block (not tested for eQTL or sQTL)")
        return int(r["blk_off"]), int(r["blk_len"])
    buf = p.read_bytes()
    h = packfmt.parse_file_header(buf)
    _expect_kind(h, (packfmt.KIND_EQTL,), p)
    for off, ln in packfmt.walk_eqtl_file(buf, h["kind"])[1]:
        g = packfmt.decode_gene_block(buf[off:off + ln], PLACEHOLDER_DOF, kind=h["kind"], expect_blk_len=ln)["details"]["gene"]
        if gene_id in (g.get("gene_id"), g.get("symbol")):
            return off, ln
    raise PackToolError(f"{p}: no block whose details name the gene {gene_id!r}")


def find_intron_block(eqtl_path: str | os.PathLike, phenotype_id: str, search_index: str | os.PathLike | None = None) -> tuple[int, int]:
    """(blk_off, blk_len) of an intron's block in the sQTL pack, read from its gene's details in
    the eQTL pack (`splice[].blk_off`/`blk_len`, SPEC section 5). `phenotype_id` is the leafcutter
    string `chr:start:end:clu_N_strand:ENSG.v`; the gene comes from its last field."""
    gene_id = phenotype_id.rsplit(":", 1)[-1].split(".")[0]
    if not gene_id.startswith("ENSG"):
        raise PackToolError(f"phenotype_id {phenotype_id!r} does not end in an ENSG gene id")
    off, ln = find_gene_block(eqtl_path, gene_id, search_index)
    d = packfmt.decode_gene_block(read_range(eqtl_path, off, ln), PLACEHOLDER_DOF, kind=packfmt.KIND_EQTL, expect_blk_len=ln)
    for s in d["details"].get("splice", []):
        if s.get("phenotype_id") == phenotype_id:
            return int(s["blk_off"]), int(s["blk_len"])
    raise PackToolError(f"{gene_id}'s details list {len(d['details'].get('splice', []))} introns, none named {phenotype_id!r}")


def read_block(path: str | os.PathLike, off: int, length: int, dof: int, *, kind: int | None = None) -> dict:
    """One decoded block (`packfmt.decode_gene_block`: header fields, `details`, per-row `pval_nominal`,
    `slope_se`, `slope`, codes, and sparse credible sets) from bytes `[off, off + length)` of an eQTL
    or sQTL pack. `kind` defaults to the file header's. `dof` is the Student-t degrees of freedom
    that turns p and SE into the slope (`manifest.packs.dof`)."""
    if kind is None:
        h = read_header(path)
        _expect_kind(h, packfmt.RESULT_KINDS, path)
        kind = h["kind"]
    return packfmt.decode_gene_block(read_range(path, off, length), dof, kind=kind, expect_blk_len=length)


def block_rows(block: dict, variants_path: str | os.PathLike) -> pa.Table:
    """The SPEC section 8 reader table for a decoded block (`read_block`), joined to its variants
    from the kind 1 file: `vidx`, then `packfmt.READER_SCHEMA` (`position`, `A1`, `A2`, `rs_number`,
    `tss_distance`, `af`, `ma_samples`, `ma_count`, `pval_nominal`, `slope`, `slope_se`, `pip`,
    `cs_id`). Rows are in `vidx` order; a block with no rows gives an empty table."""
    n = block["n_rows"]
    if n == 0:
        t = packfmt.READER_SCHEMA.empty_table()
        return t.add_column(0, "vidx", pa.array([], pa.int64()))
    pages = _decode_run(variants_path, block["var_start"], n, pos_first=block["pos_first"], pos_last=block["pos_last"])
    t = packfmt.gene_page_rows(block, pages)
    return t.add_column(0, "vidx", pa.array(np.arange(block["var_start"], block["var_start"] + n, dtype=np.int64)))


def block_codes(block: dict) -> pa.Table:
    """A decoded block's rows without variants: `row`, `vidx`, `nlp_code`, `se_code` (the two
    stored u16), and the values they decode to: `pval_nominal`, `slope_se`, `slope` (null where the
    codes say null; `slope` also null where p = 0)."""
    n = block["n_rows"]
    vs = block["var_start"] if block["var_start"] is not None else 0
    return pa.table({
        "row": pa.array(np.arange(n, dtype=np.int64)), "vidx": pa.array(np.arange(vs, vs + n, dtype=np.int64)),
        "nlp_code": pa.array(block["nlp_code"].astype(np.int64)), "se_code": pa.array(block["se_code"].astype(np.int64)),
        "pval_nominal": _nan_float(block["pval_nominal"]), "slope_se": _nan_float(block["slope_se"]), "slope": _nan_float(block["slope"]),
    })


def block_credible_sets(block: dict) -> pa.Table:
    """A decoded block's credible-set records, one per membership: `row` (0-based within the
    block), `vidx`, `pip`, `cs_id`."""
    vs = block["var_start"] if block["var_start"] is not None else 0
    return pa.table({"row": pa.array(block["cs_row"].astype(np.int64)), "vidx": pa.array(block["cs_row"].astype(np.int64) + vs),
                     "pip": pa.array(block["cs_pip"].astype(np.float64)), "cs_id": pa.array(block["cs_id"].astype(np.int64))})


RESULT_COLUMNS = ("phenotype_id", "position", "A1", "A2", "tss_distance", "pval_nominal", "slope", "slope_se")


def write_results_pack(rows: pa.Table, variants_path: str | os.PathLike, out: str | os.PathLike, kind: int, *,
                       level: int = 19, details: dict[str, dict] | None = None) -> pa.Table:
    """Nominal rows to a kind 2 (eQTL) or kind 3 (sQTL) pack, one block per phenotype.

    `rows` columns: `phenotype_id` (gene id for eQTL, leafcutter intron id for sQTL), `position`,
    `A1`, `A2` (must exist in the cis section of the kind 1 file at `variants_path`, which supplies
    `vidx`), `tss_distance` (int; `position - tss_distance` must be one value per phenotype, the
    block's `anchor`), `pval_nominal` (float in [0, 1], null allowed), `slope`, `slope_se` (floats,
    null allowed; a row's slope and SE are stored together or not at all); optional `pip` (float in
    [0, 1], null = no credible set) and `cs_id` (int 0..127, required where `pip` is set). A
    phenotype's variants must be one unbroken run of `vidx`, as SPEC requires (check 1). Blocks are
    written in order of each phenotype's first row.

    `details` (kind 2 only) maps phenotype id to its gene details JSON object (SPEC section 5); a
    gene without an entry gets the minimal `{"v": 0, "gene": {"gene_id": id}, "exons": [],
    "splice": []}`, and a gene in `details` without rows gets an empty block (the sQTL-only case).
    Returns a pointer table: `phenotype_id`, `blk_off`, `blk_len`, `var_start`, `n_var`, `anchor`,
    `n_cs` (the `search_index` columns of SPEC section 6, minus the page range and window)."""
    if kind not in packfmt.RESULT_KINDS:
        raise PackToolError(f"kind {kind} is not a results kind (2 eQTL, 3 sQTL)")
    if kind == packfmt.KIND_SQTL and details:
        raise PackToolError("an sQTL pack (kind 3) has no details; drop the details argument")
    for c in RESULT_COLUMNS:
        _column(rows, c)
    vh = read_header(variants_path)
    _expect_kind(vh, (packfmt.KIND_VARIANTS,), variants_path)
    v = variant_rows(variants_path, section="cis")
    key_to_vidx = {k: i for i, k in enumerate(zip(v["position"].to_pylist(), v["A1"].to_pylist(), v["A2"].to_pylist()))}

    pid = _strings(_column(rows, "phenotype_id"))
    pos = _ints(_column(rows, "position"), "position")
    a1, a2 = _strings(_column(rows, "A1")), _strings(_column(rows, "A2"))
    tssd = _ints(_column(rows, "tss_distance"), "tss_distance")
    p, sl, se = (_floats(_column(rows, c)) for c in ("pval_nominal", "slope", "slope_se"))
    pip = _floats(_column(rows, "pip")) if "pip" in rows.column_names else np.full(rows.num_rows, np.nan)
    csid = _floats(_column(rows, "cs_id")) if "cs_id" in rows.column_names else np.full(rows.num_rows, np.nan)
    vidx = np.empty(rows.num_rows, dtype=np.int64)
    for i, k in enumerate(zip(pos.tolist(), a1, a2)):
        j = key_to_vidx.get(k)
        if j is None:
            raise PackToolError(f"row {i} ({pid[i]}): variant {k[0]} {k[1]}/{k[2]} is not in the cis section of {variants_path}")
        vidx[i] = j

    groups: dict[str, list[int]] = {}
    for i, g in enumerate(pid):
        groups.setdefault(g, []).append(i)
    if details:
        for g in details:
            groups.setdefault(g, [])
    blocks, pointers = [], []
    for g, idx in groups.items():
        idx = np.array(idx, dtype=np.int64)
        d = None
        if kind == packfmt.KIND_EQTL:
            d = dict((details or {}).get(g) or {"gene": {"gene_id": g}, "exons": [], "splice": []})
            d.setdefault("v", packfmt.VERSION)
        if idx.size == 0:
            blocks.append(packfmt.encode_gene_block(d, None, None, [], [], [], None, None, None, level))
            pointers.append((g, None, 0, None, 0))
            continue
        idx = idx[np.argsort(vidx[idx], kind="stable")]
        run = vidx[idx]
        if np.any(np.diff(run) != 1):
            raise PackToolError(f"{g}: its variants are not one unbroken run of vidx (SPEC check 1); first break after vidx {run[np.flatnonzero(np.diff(run) != 1)[0]]}")
        anchor = pos[idx] - tssd[idx]
        if anchor.min() != anchor.max():
            raise PackToolError(f"{g}: position - tss_distance is not one value ({anchor.min()}..{anchor.max()})")
        has = ~np.isnan(pip[idx])
        cs_row = np.flatnonzero(has)
        if cs_row.size and np.any(np.isnan(csid[idx][has])):
            raise PackToolError(f"{g}: a row with pip has no cs_id")
        order = np.lexsort((csid[idx][has], cs_row))
        blocks.append(packfmt.encode_gene_block(d, int(run[0]), int(anchor[0]), p[idx], sl[idx], se[idx],
                                                cs_row[order], pip[idx][has][order], csid[idx][has][order].astype(np.int64), level,
                                                pos_first=int(pos[idx][0]), pos_last=int(pos[idx][-1])))
        pointers.append((g, int(run[0]), int(idx.size), int(anchor[0]), int(cs_row.size)))
    buf, offsets = packfmt.encode_eqtl_file(vh["chrom"], blocks, kind)
    _write_bytes(out, buf)
    return pa.table({
        "phenotype_id": pa.array([x[0] for x in pointers], pa.string()),
        "blk_off": pa.array(offsets[:-1].astype(np.int64)), "blk_len": pa.array(np.diff(offsets).astype(np.int64)),
        "var_start": pa.array([x[1] for x in pointers], pa.int64()), "n_var": pa.array([x[2] for x in pointers], pa.int64()),
        "anchor": pa.array([x[3] for x in pointers], pa.int64()), "n_cs": pa.array([x[4] for x in pointers], pa.int64()),
    })


# ---- trans pack (kind 6) -------------------------------------------------------------------------------
def walk_trans_frames(buf: bytes, *, limit: int | None = None) -> list[tuple[int, int]]:
    """(off, len) of the zstd frames after the header of a trans pack, found from the frame
    boundaries alone (each frame is decompressed to find its end; nothing is decoded). The pipeline
    finds frames through `search_index.trans_off`/`trans_len`; this is for tooling without it."""
    h = packfmt.parse_file_header(buf)
    if h["kind"] != packfmt.KIND_TRANS:
        raise PackToolError(f"trans pack: header kind {h['kind']}")
    out, off, dctx = [], packfmt.FILE_HEADER_LEN, zstandard.ZstdDecompressor()
    while off < len(buf) and (limit is None or len(out) < limit):
        dobj, fed = dctx.decompressobj(), 0
        try:
            while not dobj.eof:
                chunk = buf[off + fed:off + fed + _ZSTD_CHUNK]
                if not chunk:
                    raise PackToolError(f"trans pack: frame at byte {off} runs past the end of the file")
                dobj.decompress(chunk)
                fed += len(chunk)
        except zstandard.ZstdError as e:
            raise PackToolError(f"trans pack: frame at byte {off}: {e}") from None
        ln = fed - len(dobj.unused_data)
        out.append((off, ln))
        off += ln
    if limit is None and len(out) != h["count"]:
        raise PackToolError(f"trans pack: {len(out)} frames, header count {h['count']}")
    return out


def trans_frames(path: str | os.PathLike, search_index: str | os.PathLike | None = None, *, limit: int | None = None) -> pa.Table:
    """Every frame of a trans pack: `frame`, `trans_off`, `trans_len`, and with `search_index` also
    `gene_id`, `symbol`, `gene_version` (the frames carry no gene id; the index does). With the
    index, the pointers must be contiguous from byte 32 to the file end (SPEC section 12); without
    it the frames are walked by their zstd boundaries."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_TRANS,), p)
    if search_index is None:
        fr = walk_trans_frames(p.read_bytes(), limit=limit)
        return pa.table({"frame": pa.array(range(len(fr)), pa.int64()), "trans_off": pa.array([o for o, _ in fr], pa.int64()),
                         "trans_len": pa.array([n for _, n in fr], pa.int64())})
    si = read_search_index_file(search_index, ["gene_id", "symbol", "chr", "trans_off", "trans_len", "gene_version"])
    rows = sorted((r for r in si.to_pylist() if r["chr"] == h["chrom"] and r["trans_off"] is not None), key=lambda r: r["trans_off"])
    off = packfmt.FILE_HEADER_LEN
    for r in rows:
        if r["trans_off"] != off:
            raise PackToolError(f"{search_index}: {r['gene_id']}'s trans_off {r['trans_off']} != {off}, the byte after the previous frame")
        off += r["trans_len"]
    if len(rows) != h["count"] or off != h["bytes"]:
        raise PackToolError(f"{search_index}: {len(rows)} frames ending at byte {off} for {p} (header count {h['count']}, {h['bytes']} bytes)")
    rows = rows[:limit] if limit is not None else rows
    return pa.table({"frame": pa.array(range(len(rows)), pa.int64()), "gene_id": pa.array([r["gene_id"] for r in rows], pa.string()),
                     "symbol": pa.array([r["symbol"] for r in rows], pa.string()), "gene_version": pa.array([r["gene_version"] for r in rows], pa.int64()),
                     "trans_off": pa.array([r["trans_off"] for r in rows], pa.int64()), "trans_len": pa.array([r["trans_len"] for r in rows], pa.int64())})


def find_trans_frame(path: str | os.PathLike, gene: str, search_index: str | os.PathLike) -> dict:
    """A gene's frame in a trans pack, by gene id or symbol through `search_index`: {gene_id, symbol,
    chr, gene_version, trans_off, trans_len}. Fails when the gene has no trans rows or sits on
    another chromosome than the pack's."""
    chrom = read_header(path)["chrom"]
    rows = _search_index_rows(search_index, gene, ["trans_off", "trans_len", "gene_version"])
    rows = [r for r in rows if r["chr"] == chrom] or rows
    r = rows[0]
    if r["chr"] != chrom:
        raise PackToolError(f"{gene} is on {r['chr']}, but {path} holds {chrom}")
    if r["trans_off"] is None:
        raise PackToolError(f"{gene}: search_index says this gene has no trans rows")
    return {k: r[k] for k in ("gene_id", "symbol", "chr", "gene_version", "trans_off", "trans_len")}


def read_trans_frame(path: str | os.PathLike, off: int, length: int, dof_e: int, dof_s: int) -> dict:
    """One decoded frame (`packfmt.decode_trans_frame`) from bytes `[off, off + length)` of a trans
    pack: `n_e`, `n_s`, `k`, `nlp_max`, `beta_max`, the intron table, and per-row arrays (positions
    absolute, `beta_se` and `r2` derived with the two dof)."""
    _expect_kind(read_header(path), (packfmt.KIND_TRANS,), path)
    return packfmt.decode_trans_frame(read_range(path, off, length), dof_e, dof_s, what=f"{path} frame at {off}+{length}")


def trans_table(frame: dict, gene_chr: str, gene_id: str | None = None, gene_version: int | None = None) -> pa.Table:
    """A decoded frame as rows: `gene_id`, `qtl_type` (e or s), `phenotype_id` (the gene id for eQTL
    rows, `chr:start:end:clu_N_strand:gene.version` for sQTL rows; null when `gene_id` or
    `gene_version` is unknown), `variant_chr` (chr1..chr22, chrX), `position`, `rs_number` (null
    for 0), `af`, `pval`, `beta`, `beta_se`, `r2`, and the sQTL row's `intron_start`, `intron_end`,
    `cluster`, `strand` (null on eQTL rows). Rows are in frame order (SPEC section 12)."""
    n_e, n_s = frame["n_e"], frame["n_s"]
    n = n_e + n_s
    known = gene_id is not None and gene_version is not None
    pid = packfmt.trans_phenotype_ids(frame, gene_chr, gene_id, gene_version) if known else [None] * n
    idx = frame["intron"].astype(np.int64)
    mask = np.r_[np.ones(n_e, dtype=bool), np.zeros(n_s, dtype=bool)]

    def intron_col(values, typ):
        full = np.r_[np.zeros(n_e, dtype=np.int64), np.asarray(values, dtype=np.int64)[idx]] if n_s else np.zeros(n, dtype=np.int64)
        return pa.array(full, mask=mask, type=typ)
    strands = [None] * n_e + [frame["strand"][i] for i in idx.tolist()]
    return pa.table({
        "gene_id": pa.array([gene_id] * n, pa.string()), "qtl_type": pa.array(frame["qtl_type"], pa.string()),
        "phenotype_id": pa.array(pid, pa.string()),
        "variant_chr": pa.array([packfmt.TRANS_VARIANT_CHROMS[c - 1] for c in frame["variant_chr"].tolist()], pa.string()),
        "position": pa.array(frame["position"].astype(np.int64)), "rs_number": _zero_null_int(frame["rs_number"]),
        "af": pa.array(np.asarray(frame["af"], dtype=np.float64)), "pval": pa.array(np.asarray(frame["pval"], dtype=np.float64)),
        "beta": pa.array(np.asarray(frame["beta"], dtype=np.float64)), "beta_se": _nan_float(frame["beta_se"]), "r2": _nan_float(frame["r2"]),
        "intron_start": intron_col(frame["intron_start"], pa.int64()), "intron_end": intron_col(frame["intron_end"], pa.int64()),
        "cluster": intron_col(frame["cluster"], pa.int64()), "strand": pa.array(strands, pa.string()),
    })


TRANS_COLUMNS = ("gene_id", "qtl_type", "variant_chr", "position", "af", "pval", "beta")


def parse_phenotype_id(phenotype_id: str) -> tuple[int, int, int, str]:
    """(intron_start, intron_end, cluster, strand) from a leafcutter id `chr:start:end:clu_N_strand:gene.v`."""
    parts = phenotype_id.split(":")
    if len(parts) != 5 or not parts[3].startswith("clu_") or parts[3][-2] != "_" or parts[3][-1] not in packfmt.STRANDS:
        raise PackToolError(f"phenotype_id {phenotype_id!r} is not chr:start:end:clu_N_strand:gene.version")
    try:
        return int(parts[1]), int(parts[2]), int(parts[3][4:-2]), parts[3][-1]
    except ValueError:
        raise PackToolError(f"phenotype_id {phenotype_id!r} has a non-integer field") from None


def write_trans_pack(rows: pa.Table, out: str | os.PathLike, chrom: str, *, level: int = 19) -> pa.Table:
    """Trans rows to a kind 6 pack, one zstd frame per gene, in order of each gene's first row.

    Columns: `gene_id`, `qtl_type` (e or s), `variant_chr` (chr1..chr22, chrX), `position`, `af`
    (float in [0, 1]), `pval` (float in (0, 1]), `beta` (finite float); optional `rs_number` (null or
    0 = none). sQTL rows also need their intron: either `phenotype_id` (the leafcutter id, parsed
    here) or the four columns `intron_start`, `intron_end`, `cluster`, `strand`. `packfmt` sorts each
    gene's rows into frame order and enforces SPEC section 12 (at most 255 introns per gene, one row
    per run, chromosome, and position). `chrom` is the genes' chromosome (chrM allowed). Returns a
    pointer table: `gene_id`, `trans_off`, `trans_len`, `n_e`, `n_s`."""
    for c in TRANS_COLUMNS:
        _column(rows, c)
    n = rows.num_rows
    gid = _strings(_column(rows, "gene_id"))
    qt = _strings(_column(rows, "qtl_type"))
    if any(q not in ("e", "s") for q in qt):
        raise PackToolError("qtl_type must be 'e' or 's' on every row")
    vchr = _strings(_column(rows, "variant_chr"))
    pos = _ints(_column(rows, "position"), "position")
    rs = np.nan_to_num(_floats(_column(rows, "rs_number")), nan=0.0).astype(np.int64) if "rs_number" in rows.column_names else np.zeros(n, dtype=np.int64)
    af, pv, beta = (_floats(_column(rows, c)) for c in ("af", "pval", "beta"))
    is_s = np.array([q == "s" for q in qt], dtype=bool)
    if all(c in rows.column_names for c in ("intron_start", "intron_end", "cluster", "strand")):
        istart, iend, iclu = (_floats(_column(rows, c)) for c in ("intron_start", "intron_end", "cluster"))
        strand = _strings(_column(rows, "strand"))
    elif "phenotype_id" in rows.column_names:
        pid = _strings(_column(rows, "phenotype_id"))
        istart, iend, iclu = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan)
        strand = [None] * n
        for i in np.flatnonzero(is_s).tolist():
            if pid[i] is None:
                raise PackToolError(f"row {i}: an sQTL row needs a phenotype_id")
            istart[i], iend[i], iclu[i], strand[i] = parse_phenotype_id(pid[i])
    elif is_s.any():
        raise PackToolError("sQTL rows need phenotype_id, or intron_start, intron_end, cluster, strand")
    else:
        istart, iend, iclu, strand = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), [None] * n
    groups: dict[str, list[int]] = {}
    for i, g in enumerate(gid):
        if g is None:
            raise PackToolError(f"row {i}: gene_id is null")
        groups.setdefault(g, []).append(i)
    frames, ptr = [], {"gene_id": [], "trans_off": [], "trans_len": [], "n_e": [], "n_s": []}
    off = packfmt.FILE_HEADER_LEN
    for g, idx in groups.items():
        sel = np.array(idx, dtype=np.int64)
        try:
            frame = packfmt.encode_trans_frame([qt[i] for i in idx], [vchr[i] for i in idx], pos[sel], rs[sel], af[sel], pv[sel], beta[sel],
                                               istart[sel], iend[sel], iclu[sel], [strand[i] for i in idx], level)
        except ValueError as e:
            raise PackToolError(f"{g}: {e}") from None
        frames.append(frame)
        ptr["gene_id"].append(g); ptr["trans_off"].append(off); ptr["trans_len"].append(len(frame))
        ptr["n_e"].append(int((~is_s[sel]).sum())); ptr["n_s"].append(int(is_s[sel].sum()))
        off += len(frame)
    if off > packfmt.U32_MAX:
        raise PackToolError(f"trans pack: {off} bytes exceeds 4 GiB")
    _write_bytes(out, packfmt.file_header(packfmt.KIND_TRANS, chrom, len(frames), 0) + b"".join(frames))
    return pa.table({k: pa.array(v, pa.string() if k == "gene_id" else pa.int64()) for k, v in ptr.items()})


# ---- GWAS pack and index (kinds 4 and 5) -----------------------------------------------------------
def gwas_index(path: str | os.PathLike) -> dict:
    """`packfmt.decode_gwas_index`: {block_rows, n_values, chroms: {name: (first_position, end_offset)}}."""
    p = Path(path)
    _expect_kind(read_header(p), (packfmt.KIND_GWAS_INDEX,), p)
    return packfmt.decode_gwas_index(p.read_bytes())


def gwas_index_table(path: str | os.PathLike) -> pa.Table:
    """The GWAS index as one row per block: `chr`, `block`, `first_position`, `byte_start`,
    `byte_end` (exclusive) in that chromosome's pack."""
    idx = gwas_index(path)
    cols = {"chr": [], "block": [], "first_position": [], "byte_start": [], "byte_end": []}
    for name, (fp, eo) in idx["chroms"].items():
        starts = np.r_[packfmt.FILE_HEADER_LEN, eo[:-1]]
        cols["chr"] += [name] * len(fp); cols["block"] += list(range(len(fp)))
        cols["first_position"] += fp.tolist(); cols["byte_start"] += starts.tolist(); cols["byte_end"] += eo.tolist()
    return pa.table({k: pa.array(v, pa.string() if k == "chr" else pa.int64()) for k, v in cols.items()})


def _gwas_table(parts: list[dict], block_ids: list[int]) -> pa.Table:
    cat = lambda k: np.concatenate([np.asarray(d[k]) for d in parts]) if parts else np.array([])
    return pa.table({
        "block": pa.array(np.concatenate([np.full(d["rows"], b) for d, b in zip(parts, block_ids)]).astype(np.int64) if parts else np.array([], np.int64)),
        "position": pa.array(cat("position").astype(np.int64)),
        "ea": pa.array([x for d in parts for x in d["ea"]], pa.string()), "nea": pa.array([x for d in parts for x in d["nea"]], pa.string()),
        "rs_number": _zero_null_int(cat("rs_number")),
        "beta": pa.array(cat("beta").astype(np.float64)), "se": pa.array(cat("se").astype(np.float64)), "eaf": pa.array(cat("eaf").astype(np.float64)),
        "p": pa.array(cat("p").astype(np.float64)), "n": pa.array(cat("n").astype(np.int64)),
    })


def gwas_rows(path: str | os.PathLike, index_path: str | os.PathLike, *, lo: int | None = None, hi: int | None = None,
              block: int | None = None) -> pa.Table:
    """GWAS rows of one chromosome's pack: `block`, `position`, `ea` (effect allele), `nea`,
    `rs_number` (null when none), `beta`, `se`, `eaf`, `p`, `n`. With `lo` and `hi`: the rows with
    `lo <= position <= hi`, fetched as SPEC section 11's window rule does (one byte range through
    the index, then filtered). With `block`: every row of that block."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_GWAS,), p)
    idx = gwas_index(index_path)
    if h["chrom"] not in idx["chroms"]:
        raise PackToolError(f"{index_path}: no chromosome {h['chrom']!r} (it has {list(idx['chroms'])})")
    if idx["block_rows"] != h["page_size"]:
        raise PackToolError(f"{index_path}: {idx['block_rows']} rows per block, but {p} says {h['page_size']}")
    fp, eo = idx["chroms"][h["chrom"]]
    if len(eo) and eo[-1] != h["bytes"]:
        raise PackToolError(f"{index_path}: last end offset {eo[-1]} != file size {h['bytes']}")
    if block is not None:
        if not 0 <= block < len(fp):
            raise PackToolError(f"{p}: block {block} is outside 0..{len(fp) - 1}")
        ks, ke = block, block
        b0, b1 = (packfmt.FILE_HEADER_LEN if block == 0 else int(eo[block - 1])), int(eo[block])
    else:
        if lo is None or hi is None:
            raise PackToolError("give --lo and --hi, or --block")
        w = packfmt.gwas_window(fp, eo, lo, hi)
        if w is None:
            return _gwas_table([], [])
        ks, ke, b0, b1 = w
    buf = read_range(p, b0, b1 - b0)
    parts = []
    B, nblk = h["page_size"], len(fp)
    for k in range(ks, ke + 1):
        s = (packfmt.FILE_HEADER_LEN if k == 0 else int(eo[k - 1])) - b0
        expect = B if k < nblk - 1 else h["count"] - B * (nblk - 1)
        d = packfmt.decode_gwas_block(buf[s:int(eo[k]) - b0], idx["n_values"], expect_rows=expect, what=f"{p} block {k}")
        if d["position"][0] != fp[k]:
            raise PackToolError(f"{p} block {k}: first position {d['position'][0]} != index {fp[k]}")
        parts.append(d)
    t = _gwas_table(parts, list(range(ks, ke + 1)))
    if block is None:
        pos = t["position"].to_numpy()
        t = t.filter(pa.array((pos >= lo) & (pos <= hi)))
    return t


GWAS_COLUMNS = ("chr", "position", "ea", "nea", "rs_number", "beta", "se", "eaf", "p", "n")


def chrom_order(name: str) -> tuple[int, str]:
    """Sort key: chr1..chr22, chrX, chrY, chrM, then anything else by name."""
    s = name[3:] if name.startswith("chr") else name
    if s.isdigit():
        return int(s), ""
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(s, 100), name


def write_gwas_packs(table: pa.Table, out_dir: str | os.PathLike, *, block_rows: int = 2048, level: int = 19) -> pa.Table:
    """GWAS rows to kind 4 packs, `<out_dir>/<chr>.qbg`, and the kind 5 index `<out_dir>/gwas_index.bin`.

    Columns: `chr`, `position` (int), `ea`, `nea` (alleles as in the source), `rs_number` (int, null
    or 0 = none), `beta`, `se`, `eaf` (floats printed to at most 4 decimals), `p` (float in (0, 1]
    with at most 4 significant digits), `n` (int; at most 255 distinct values). The lossless rules
    of SPEC section 11 are enforced by `packfmt.gwas_codes`; a row that breaks one fails the write
    with its row number. Rows are sorted here by chromosome, position, ea, nea, rs_number, p.
    Returns one row per block: `chr`, `block`, `first_position`, `byte_start`, `byte_end`."""
    for c in GWAS_COLUMNS:
        _column(table, c)
    t = table.set_column(table.schema.get_field_index("rs_number"), "rs_number",
                         pa.array(np.nan_to_num(_floats(_column(table, "rs_number")), nan=0.0).astype(np.int64)))
    chroms = _strings(t["chr"].combine_chunks())
    ordv = pa.array([chrom_order(c)[0] for c in chroms], pa.int64())
    t = t.append_column("_ord", ordv).sort_by([("_ord", "ascending"), ("chr", "ascending"), ("position", "ascending"),
                                                ("ea", "ascending"), ("nea", "ascending"), ("rs_number", "ascending"), ("p", "ascending")]).drop_columns(["_ord"])
    n_values = sorted(set(_ints(_column(t, "n"), "n").tolist()))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    index_chroms, rows = [], {"chr": [], "block": [], "first_position": [], "byte_start": [], "byte_end": []}
    chr_col = np.array(_strings(t["chr"].combine_chunks()), dtype=object)
    for name in sorted(set(chr_col.tolist()), key=chrom_order):
        sub = t.filter(pa.array(chr_col == name))
        codes = packfmt.gwas_codes(_ints(sub["position"], "position"), _floats(sub["beta"]), _floats(sub["se"]), _floats(sub["eaf"]),
                                   _floats(sub["p"]), _ints(sub["rs_number"], "rs_number"), _ints(sub["n"], "n"), n_values)
        ea, nea = _strings(sub["ea"].combine_chunks()), _strings(sub["nea"].combine_chunks())
        parts, fp, eo, off = [packfmt.file_header(packfmt.KIND_GWAS, name, sub.num_rows, block_rows)], [], [], packfmt.FILE_HEADER_LEN
        for s in range(0, sub.num_rows, block_rows):
            e = min(s + block_rows, sub.num_rows)
            frame = packfmt.encode_gwas_block({k: v[s:e] for k, v in codes.items()}, ea[s:e], nea[s:e], level)
            parts.append(frame)
            rows["chr"].append(name); rows["block"].append(len(fp)); rows["first_position"].append(int(codes["position"][s]))
            rows["byte_start"].append(off); off += len(frame); rows["byte_end"].append(off)
            fp.append(int(codes["position"][s])); eo.append(off)
        _write_bytes(out / f"{name}{KINDS[packfmt.KIND_GWAS][1]}", b"".join(parts))
        index_chroms.append((name, np.array(fp), np.array(eo)))
    _write_bytes(out / "gwas_index.bin", packfmt.encode_gwas_index(n_values, index_chroms, block_rows, level))
    return pa.table({k: pa.array(v, pa.string() if k == "chr" else pa.int64()) for k, v in rows.items()})


# ---- whole-file check -------------------------------------------------------------------------------
def check_file(path: str | os.PathLike, *, variants: str | os.PathLike | None = None, index: str | os.PathLike | None = None,
               search_index: str | os.PathLike | None = None, dof: int | None = None,
               variant_index: str | os.PathLike | None = None) -> dict:
    """Decode a whole pack with `packfmt` and return counts. Any SPEC rule that fails raises with
    the rule. Kind 1: every page of both sections. Kinds 2 and 3: every block (with `variants`, also
    each block's run against the variants file, SPEC section 7 steps 3 to 5; without `dof` a
    placeholder is used, which only affects the slope values, not the checks). Kind 4: every block
    through the index (`index` required). Kind 5: the index itself. Kind 6: every frame, through
    `search_index` pointers when given (checked contiguous) or by walking the zstd frames. Reads the
    whole file into memory; fine for one chromosome. Kind 7: every frame, through `variant_index`
    offsets when given (checked against the file) or by walking the zstd frames. Kind 8: every block.
    Kind 9: the startup file itself."""
    p = Path(path)
    h = read_header(p)
    k = h["kind"]
    out = {"path": str(p), "kind": k, "kind_name": h["kind_name"], "chrom": h["chrom"], "bytes": h["bytes"], "header_count": h["count"]}
    if k == packfmt.KIND_VARIANTS:
        d = packfmt.decode_variants_file(p.read_bytes())
        n_cis = h["n_cis"]
        out.update(pages=len(d["pages"]), records=int(sum(pg["n"] for pg in d["pages"])), n_cis=n_cis, trans_only=h["count"] - n_cis,
                   page_size=h["page_size"], codec=d.get("codec"),
                   heap_records=int(((d["allele_code"] == 0) & ~d["no_alleles"]).sum()) if d["pages"] else 0,
                   no_alleles=int(d["no_alleles"].sum()) if d["pages"] else 0)
        return out
    if k in packfmt.RESULT_KINDS:
        buf = p.read_bytes()
        _, blocks = packfmt.walk_eqtl_file(buf, k)
        rows = cs = dz = empty = 0
        for off, ln in blocks:
            b = packfmt.decode_gene_block(buf[off:off + ln], dof if dof is not None else PLACEHOLDER_DOF, kind=k, expect_blk_len=ln)
            rows += b["n_rows"]; cs += b["n_cs"]; dz += b["details_zlen"]; empty += b["n_rows"] == 0
            if variants is not None:
                block_rows(b, variants)
        out.update(blocks=len(blocks), rows=rows, credible_set_records=cs, empty_blocks=empty, details_bytes=dz,
                   checked_against_variants=variants is not None)
        return out
    if k == packfmt.KIND_GWAS:
        if index is None:
            raise PackToolError("a GWAS pack's blocks are found only through the index: give --index gwas_index.bin")
        idx = gwas_index(index)
        fp, _ = idx["chroms"].get(h["chrom"], (None, None))
        if fp is None:
            raise PackToolError(f"{index}: no chromosome {h['chrom']!r}")
        total = 0
        for b in range(len(fp)):
            total += gwas_rows(p, index, block=b).num_rows
        if total != h["count"]:
            raise PackToolError(f"{p}: blocks hold {total} rows, header count {h['count']}")
        out.update(blocks=len(fp), rows=total, block_rows=h["page_size"])
        return out
    if k == packfmt.KIND_GWAS_INDEX:
        idx = gwas_index(p)
        out.update(chromosomes=list(idx["chroms"]), blocks={c: len(v[0]) for c, v in idx["chroms"].items()},
                   block_rows=idx["block_rows"], n_values=idx["n_values"])
        return out
    if k == packfmt.KIND_TRANS:
        fr = trans_frames(p, search_index)
        buf = p.read_bytes()
        n_e = n_s = 0
        max_k = max_rows = 0
        for off, ln in zip(fr["trans_off"].to_pylist(), fr["trans_len"].to_pylist()):
            f = packfmt.decode_trans_frame(buf[off:off + ln], PLACEHOLDER_DOF, PLACEHOLDER_DOF, what=f"{p} frame at {off}+{ln}")
            n_e += f["n_e"]; n_s += f["n_s"]; max_k = max(max_k, f["k"]); max_rows = max(max_rows, f["n_e"] + f["n_s"])
        if fr.num_rows != h["count"]:
            raise PackToolError(f"{p}: {fr.num_rows} frames, header count {h['count']}")
        out.update(frames=fr.num_rows, rows=n_e + n_s, eqtl_rows=n_e, sqtl_rows=n_s, max_rows_in_a_frame=max_rows, max_introns=max_k,
                   found_through="search_index" if search_index is not None else "zstd frame boundaries")
        return out
    if k == packfmt.KIND_HITS:
        fr = hits_frames(p, variant_index=variant_index)
        buf = p.read_bytes()
        rows = by_kind = None
        n_variants = 0
        for i, (off, ln) in enumerate(zip(fr["hits_off"].to_pylist(), fr["hits_len"].to_pylist())):
            f = packfmt.decode_hits_frame(buf[off:off + ln], i * h["page_size"], what=f"{p} frame at {off}+{ln}")
            counts = np.bincount(f["kind"], minlength=packfmt.HITS_MAX_KIND + 1)
            by_kind = counts if by_kind is None else by_kind + counts
            rows = (rows or 0) + f["n_rows"]
            n_variants += f["n_variants"]
        if fr.num_rows != h["count"]:
            raise PackToolError(f"{p}: {fr.num_rows} frames, header count {h['count']}")
        out.update(frames=fr.num_rows, variants=n_variants, rows=rows or 0, frame_variants=h["page_size"],
                   rows_by_kind={packfmt.HITS_KIND_NAMES[i]: int(c) for i, c in enumerate(by_kind if by_kind is not None else [])},
                   found_through="variant_index" if variant_index is not None else "zstd frame boundaries")
        return out
    if k == packfmt.KIND_RSID:
        B = h["page_size"]
        n_blocks = -(-h["count"] // B)
        prev = None
        for b in range(n_blocks):
            off, ln = packfmt.rsid_block_range(b, h["count"], B)
            d = packfmt.decode_rsid_block(read_range(p, off, ln), what=f"{p} block {b}")
            if prev is not None and int(d["rs_number"][0]) <= prev:
                raise PackToolError(f"{p}: block {b} starts at rs{int(d['rs_number'][0])}, not above the previous block's last")
            prev = int(d["rs_number"][-1])
        if h["bytes"] != packfmt.FILE_HEADER_LEN + h["count"] * packfmt.RSID_RECORD_LEN:
            raise PackToolError(f"{p}: {h['bytes']} bytes for {h['count']} records")
        out.update(records=h["count"], blocks=n_blocks, block_records=B, last_rs_number=prev)
        return out
    if k == packfmt.KIND_VARIANT_INDEX:
        idx = read_variant_index(p)
        out.update(chromosomes=list(idx["chroms"]), page_size=idx["page_size"], frame_variants=idx["frame_variants"],
                   rsid_records=idx["rsid_n_records"], rsid_blocks=idx["rsid_n_blocks"],
                   pages=sum(c["n_pages_cis"] + c["n_pages_trans"] for c in idx["chroms"].values()),
                   frames=sum(c["n_frames"] for c in idx["chroms"].values()),
                   variants=sum(c["n_cis"] + c["n_trans_only"] for c in idx["chroms"].values()))
        return out
    raise PackToolError(f"{p}: kind {k} has no checker yet")


# ---- CLI ----------------------------------------------------------------------------------------------
def _print_json(obj) -> None:
    print(json.dumps(_clean(obj), indent=2))


def _kind_from_out(out: str, kind: str | None) -> int:
    if kind is not None:
        if kind not in KIND_BY_NAME:
            raise PackToolError(f"unknown kind {kind!r}; one of {sorted(KIND_BY_NAME)}")
        return KIND_BY_NAME[kind]
    ext = Path(out).suffix
    if ext not in KIND_BY_EXT:
        raise PackToolError(f"cannot tell the kind from {out!r}: give --kind, or use an extension in {sorted(KIND_BY_EXT)}")
    return KIND_BY_EXT[ext]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m pipeline.packtool", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    out_help = "TABLE to write (.parquet, .arrow, .feather, .tsv, .csv, .json); default - is TSV on stdout"
    dof_help = "degrees of freedom, for the derived slope (default: packs.dof from the nearest manifest.json above the pack)"
    manifest_help = "manifest.json to read packs.dof from (default: the nearest one above the pack)"

    s = sub.add_parser("header", help="print a pack's 32-byte file header as JSON")
    s.add_argument("pack", help="any pack file (.qbv, .qbe, .qbs, .qbt, .qbg) or gwas_index.bin")

    s = sub.add_parser("blocks", help="list the blocks of an eQTL or sQTL pack from their headers")
    s.add_argument("pack", help="an eQTL (.qbe) or sQTL (.qbs) pack"); s.add_argument("--limit", type=int, help="stop after this many blocks")
    s.add_argument("--gene-ids", action="store_true", help="decode details for gene_id and symbol (kind 2)")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("block", help="decode one eQTL or sQTL block to rows")
    s.add_argument("pack", help="an eQTL (.qbe) or sQTL (.qbs) pack")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene", help="gene id or symbol (kind 2; scans details unless --search-index is given)")
    g.add_argument("--intron", help="intron phenotype_id (kind 3; needs --eqtl for the gene's details)")
    g.add_argument("--index", type=int, help="block number, 0-based, in file order")
    g.add_argument("--off", type=int, help="byte offset (with --len)")
    s.add_argument("--len", type=int, help="byte length (with --off)"); s.add_argument("--eqtl", help="the chromosome's eQTL pack, for --intron")
    s.add_argument("--search-index", help="search_index.parquet for a direct --gene lookup")
    s.add_argument("--variants", help="the chromosome's variants pack; with it, rows carry position, alleles, and tss_distance")
    s.add_argument("--dof", type=int, help=dof_help); s.add_argument("--manifest", help=manifest_help)
    s.add_argument("--details", action="store_true", help="print the gene details JSON instead of rows")
    s.add_argument("--cs", action="store_true", help="print the credible-set records instead of rows")
    s.add_argument("--header", action="store_true", help="print the block header fields as JSON instead of rows")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("variants", help="decode variant pages to rows")
    s.add_argument("pack", help="a variants pack (.qbv)")
    s.add_argument("--vidx", type=int, help="first variant index, its 0-based rank in the chromosome (with --n; a phenotype's var_start)")
    s.add_argument("--n", type=int, help="number of variants from --vidx (a phenotype's n_var)")
    s.add_argument("--off", type=int, help="byte offset of whole pages (with --len; search_index var_off, var_len)")
    s.add_argument("--len", type=int, help="byte length (with --off)")
    s.add_argument("--section", choices=["cis", "trans"], help="every record of one section: the cis-tested variants, or the trans-only ones (with no range: every record)")
    s.add_argument("--whole-pages", action="store_true", help="with --vidx/--n: every record of the covering pages")
    s.add_argument("--pages", action="store_true", help="list the pages (offsets, counts, section) instead of rows")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("frames", help="list the frames of a trans pack")
    s.add_argument("pack", help="a trans pack (.qbt)"); s.add_argument("--search-index", help="search_index.parquet: adds gene_id, symbol, gene_version and checks contiguity")
    s.add_argument("--limit", type=int, help="stop after this many frames"); s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("trans", help="decode one gene's trans frame to rows")
    s.add_argument("pack", help="a trans pack (.qbt)")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene", help="gene id or symbol (needs --search-index; frames carry no gene id)")
    g.add_argument("--off", type=int, help="byte offset (with --len)")
    s.add_argument("--len", type=int, help="byte length (with --off)"); s.add_argument("--search-index", help="search_index.parquet (trans_off, trans_len per gene)")
    s.add_argument("--gene-id", help="with --off: the gene id, for phenotype ids"); s.add_argument("--gene-version", type=int, help="with --off: the gene version, for phenotype ids")
    s.add_argument("--dof-eqtl", type=int, help="eQTL degrees of freedom, for beta_se and r2 (default: packs.dof.eqtl from the nearest manifest.json)")
    s.add_argument("--dof-sqtl", type=int, help="sQTL degrees of freedom (default: packs.dof.sqtl from the nearest manifest.json)")
    s.add_argument("--manifest", help=manifest_help)
    s.add_argument("--header", action="store_true", help="print the frame header and intron table as JSON instead of rows")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("gwas", help="decode GWAS rows in a position window, or one block")
    s.add_argument("pack", help="a GWAS pack (.qbg)"); s.add_argument("--index", required=True, help="gwas_index.bin")
    s.add_argument("--lo", type=int, help="first position of the window, inclusive (with --hi)")
    s.add_argument("--hi", type=int, help="last position of the window, inclusive (with --lo)")
    s.add_argument("--block", type=int, help="one block by number (0-based within the chromosome) instead of a window")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("gwas-index", help="list every GWAS block from the index")
    s.add_argument("index", help="gwas_index.bin"); s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("check", help="decode a whole pack and print counts; fails on the first broken rule")
    s.add_argument("pack", help="any pack file, or gwas_index.bin")
    s.add_argument("--variants", help="eQTL/sQTL packs: the chromosome's variants pack, to also check every block's run against it")
    s.add_argument("--index", help="GWAS packs: gwas_index.bin (required; blocks are found only through it)")
    s.add_argument("--search-index", help="trans packs: search_index.parquet, to find frames by its pointers instead of by walking the zstd frames")
    s.add_argument("--dof", type=int, help="eQTL/sQTL packs: degrees of freedom (affects only the derived slopes, not the checks)")
    s.add_argument("--variant-index", help="hits packs: variant_index.qbx, to find frames by its offsets instead of by walking the file")

    s = sub.add_parser("hits", help="decode one variant's rows, or a whole frame, of a hits pack")
    s.add_argument("pack", help="a hits pack (.qbh)")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--vidx", type=int, help="one variant's index in its chromosome's variants file")
    g.add_argument("--frame", type=int, help="a whole frame by number, 0-based")
    g.add_argument("--frames", action="store_true", help="list the frames (offsets and lengths) instead of rows")
    s.add_argument("--variant-index", help="variant_index.qbx: find the frame by its offset instead of walking the file")
    s.add_argument("--search-index", help="search_index.parquet: adds gene_id, symbol and the rebuilt sQTL phenotype_id")
    s.add_argument("--limit", type=int, help="with --frames: stop after this many")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("rsid", help="look one rsID up in the rsID index, or list a block")
    s.add_argument("pack", help="an rsID index (.qbr)")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--rs", type=int, help="an rs number, without the 'rs' (the block math of SPEC section 14)")
    g.add_argument("--block", type=int, help="list one block's records, 0-based")
    s.add_argument("--variant-index", help="variant_index.qbx: take the block samples from the startup file, as the browser does")
    s.add_argument("--limit", type=int, help="with --block: stop after this many records")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("variant-index", help="print the startup file: one row per chromosome, or its offsets")
    s.add_argument("pack", help="a variant index (.qbx)")
    s.add_argument("--chrom", help="print one chromosome's page and frame offsets as JSON instead")
    s.add_argument("-o", "--out", default="-", help=out_help)

    s = sub.add_parser("pack-hits", help="hit rows (vidx, kind, gene, and the values of each kind) to a hits pack")
    s.add_argument("table", help="TABLE of hit rows (see write_hits_pack for the columns)")
    s.add_argument("-o", "--out", required=True, help="output pack (.qbh)")
    s.add_argument("--chrom", required=True, help="the variants' chromosome, written in the header")
    s.add_argument("--n-variants", type=int, help="the chromosome's variant count (default: one past the largest vidx)")
    s.add_argument("--frame-variants", type=int, default=packfmt.HITS_FRAME_VARIANTS, help="variants per frame (default 1024)")
    s.add_argument("--pointers", help="write the frame table here (TABLE path)")
    s.add_argument("--level", type=int, default=19, help="zstd level (default 19)")

    s = sub.add_parser("pack-rsid", help="rsID rows (rs_number, chr, vidx) to an rsID index")
    s.add_argument("table", help="TABLE of rsID rows (see write_rsid_index for the columns)")
    s.add_argument("-o", "--out", required=True, help="output index (.qbr)")
    s.add_argument("--block-records", type=int, default=packfmt.RSID_BLOCK_RECORDS, help="records per block (default 4096)")
    s.add_argument("--blocks", help="write the block table here (TABLE path)")

    s = sub.add_parser("pack-variant-index", help="variants files, hits packs and an rsID index to the startup file")
    s.add_argument("-o", "--out", required=True, help="output index (.qbx)")
    s.add_argument("--variants", required=True, help="directory of <chr>.qbv files")
    s.add_argument("--hits", required=True, help="directory of <chr>.qbh files")
    s.add_argument("--rsid-index", required=True, help="rsid_index.qbr")
    s.add_argument("--level", type=int, default=19, help="zstd level (default 19)")

    s = sub.add_parser("pack-variants", help="rows (position, A1, A2, rs_number, af, ma_samples, ma_count, match) to a variants pack")
    s.add_argument("table", help="TABLE of cis variants (see write_variants_pack for the columns); a boolean in_cis column splits off the trans-only section instead of --trans-only")
    s.add_argument("-o", "--out", required=True, help="output pack (.qbv)"); s.add_argument("--chrom", required=True, help="chromosome name written in the header, e.g. chr21")
    s.add_argument("--trans-only", help="TABLE of trans-only variants (position, A1, A2, rs_number, af, match) for the second section")
    s.add_argument("--page-size", type=int, default=512, help="variants per page (default 512)"); s.add_argument("--codec", choices=sorted(packfmt.CODECS), default="zstd", help="page compression (default zstd)")
    s.add_argument("--level", type=int, default=19, help="zstd level (default 19)")
    s.add_argument("--pages", help="also write the page table here (TABLE path)")

    s = sub.add_parser("pack-results", help="nominal rows to an eQTL (.qbe) or sQTL (.qbs) pack")
    s.add_argument("table", help="TABLE of nominal rows (see write_results_pack for the columns)")
    s.add_argument("--variants", required=True, help="the chromosome's variants pack, which defines each row's vidx"); s.add_argument("-o", "--out", required=True, help="output pack (.qbe or .qbs)")
    s.add_argument("--kind", choices=["eqtl", "sqtl"], help="default: from the output extension")
    s.add_argument("--details", help="JSON object mapping gene id to its details object (kind 2)")
    s.add_argument("--pointers", help="write the pointer table here (TABLE path)")
    s.add_argument("--level", type=int, default=19, help="zstd level for the details frames (default 19)")

    s = sub.add_parser("pack-trans", help="trans rows (gene_id, qtl_type, variant_chr, position, rs_number, af, pval, beta, phenotype_id) to a trans pack")
    s.add_argument("table", help="TABLE of trans rows (see write_trans_pack for the columns)"); s.add_argument("-o", "--out", required=True, help="output pack (.qbt)")
    s.add_argument("--chrom", required=True, help="the genes' chromosome, written in the header")
    s.add_argument("--pointers", help="write the pointer table here (TABLE path)"); s.add_argument("--level", type=int, default=19, help="zstd level (default 19)")

    s = sub.add_parser("pack-gwas", help="GWAS rows (chr, position, ea, nea, rs_number, beta, se, eaf, p, n) to <dir>/<chr>.qbg and <dir>/gwas_index.bin")
    s.add_argument("table", help="TABLE of GWAS rows, any number of chromosomes (see write_gwas_packs for the columns)"); s.add_argument("-o", "--out", required=True, help="output directory")
    s.add_argument("--block-rows", type=int, default=2048, help="rows per zstd block (default 2048)"); s.add_argument("--level", type=int, default=19, help="zstd level (default 19)")
    s.add_argument("--blocks", help="write the block table here (TABLE path)")

    a = ap.parse_args(argv)
    try:
        return _run(a)
    except (PackToolError, ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    except BrokenPipeError:                       # `| head` closed stdout; not an error
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0


def _run(a) -> int:
    if a.cmd == "header":
        _print_json(read_header(a.pack))
    elif a.cmd == "blocks":
        write_table(list_blocks(a.pack, limit=a.limit, gene_ids=a.gene_ids), a.out)
    elif a.cmd == "block":
        h = read_header(a.pack)
        _expect_kind(h, packfmt.RESULT_KINDS, a.pack)
        if a.gene is not None:
            off, ln = find_gene_block(a.pack, a.gene, a.search_index)
        elif a.intron is not None:
            if a.eqtl is None:
                raise PackToolError("--intron needs --eqtl PACK.qbe (the intron's block pointer lives in its gene's details)")
            off, ln = find_intron_block(a.eqtl, a.intron, a.search_index)
        elif a.index is not None:
            t = list_blocks(a.pack)
            if not 0 <= a.index < t.num_rows:
                raise PackToolError(f"{a.pack}: block {a.index} is outside 0..{t.num_rows - 1}")
            off, ln = int(t["blk_off"][a.index].as_py()), int(t["blk_len"][a.index].as_py())
        else:
            if a.len is None:
                raise PackToolError("--off needs --len")
            off, ln = a.off, a.len
        dof = resolve_dof(h["kind"], a.dof, a.manifest, a.pack)
        b = read_block(a.pack, off, ln, dof, kind=h["kind"])
        if a.header:
            _print_json({"blk_off": off, **{k: b[k] for k in ("blk_len", "n_rows", "var_start", "anchor", "pos_first", "pos_last", "n_cs",
                                                              "nlp_max", "lse_min", "lse_max", "details_zlen", "details_len", "dof")}})
        elif a.details:
            if b["details"] is None:
                raise PackToolError("an sQTL block (kind 3) has no details")
            _print_json(b["details"])
        elif a.cs:
            write_table(block_credible_sets(b), a.out)
        elif a.variants:
            write_table(block_rows(b, a.variants), a.out)
        else:
            write_table(block_codes(b), a.out)
    elif a.cmd == "variants":
        if a.pages:
            write_table(page_table(a.pack), a.out)
        else:
            write_table(variant_rows(a.pack, vidx=a.vidx, n=a.n, off=a.off, length=a.len, section=a.section, whole_pages=a.whole_pages), a.out)
    elif a.cmd == "frames":
        write_table(trans_frames(a.pack, a.search_index, limit=a.limit), a.out)
    elif a.cmd == "trans":
        h = read_header(a.pack)
        _expect_kind(h, (packfmt.KIND_TRANS,), a.pack)
        if a.gene is not None:
            if a.search_index is None:
                raise PackToolError("--gene needs --search-index search_index.parquet (trans frames carry no gene id)")
            info = find_trans_frame(a.pack, a.gene, a.search_index)
            off, ln, gene_id, version = info["trans_off"], info["trans_len"], info["gene_id"], info["gene_version"]
        else:
            if a.len is None:
                raise PackToolError("--off needs --len")
            off, ln, gene_id, version = a.off, a.len, a.gene_id, a.gene_version
        dof_e, dof_s = resolve_dofs(a.dof_eqtl, a.dof_sqtl, a.manifest, a.pack)
        f = read_trans_frame(a.pack, off, ln, dof_e, dof_s)
        if a.header:
            _print_json({"trans_off": off, "trans_len": ln, "gene_id": gene_id, "gene_version": version, "dof": {"eqtl": dof_e, "sqtl": dof_s},
                         **{k: f[k] for k in ("n_e", "n_s", "k", "nlp_max", "beta_max", "payload_len")},
                         "introns": [{"intron_start": s, "intron_end": e, "cluster": c, "strand": st} for s, e, c, st in
                                     zip(f["intron_start"].tolist(), f["intron_end"].tolist(), f["cluster"].tolist(), f["strand"])]})
        else:
            write_table(trans_table(f, h["chrom"], gene_id, version), a.out)
    elif a.cmd == "gwas":
        write_table(gwas_rows(a.pack, a.index, lo=a.lo, hi=a.hi, block=a.block), a.out)
    elif a.cmd == "gwas-index":
        write_table(gwas_index_table(a.index), a.out)
    elif a.cmd == "hits":
        if a.frames:
            write_table(hits_frames(a.pack, limit=a.limit, variant_index=a.variant_index), a.out)
        else:
            write_table(hits_rows(a.pack, vidx=a.vidx, frame=a.frame, variant_index=a.variant_index,
                                  search_index=a.search_index), a.out)
    elif a.cmd == "rsid":
        if a.block is not None:
            write_table(rsid_records(a.pack, block=a.block, limit=a.limit), a.out)
        else:
            got = rsid_lookup(a.pack, a.rs, a.variant_index)
            if got is None:
                raise PackToolError(f"rs{a.rs} is not in {a.pack}")
            _print_json(got)
    elif a.cmd == "variant-index":
        if a.chrom:
            idx = read_variant_index(a.pack)
            if a.chrom not in idx["chroms"]:
                raise PackToolError(f"{a.pack}: no chromosome {a.chrom!r}")
            c = idx["chroms"][a.chrom]
            _print_json({"chr": a.chrom, "page_size": idx["page_size"], "frame_variants": idx["frame_variants"],
                         **{k: c[k] for k in ("n_cis", "n_trans_only", "n_pages_cis", "n_pages_trans", "n_frames")},
                         "page_off": c["page_off"].tolist(), "page_first_position": c["page_first_position"].tolist(),
                         "hits_off": c["hits_off"].tolist()})
        else:
            write_table(variant_index_table(a.pack), a.out)
    elif a.cmd == "pack-hits":
        ptr = write_hits_pack(read_table(a.table), a.out, a.chrom, n_variants=a.n_variants,
                              frame_variants=a.frame_variants, level=a.level)
        if a.pointers:
            write_table(ptr, a.pointers)
        _print_json({**read_header(a.out), "frames": ptr.num_rows, "rows": int(sum(ptr["rows"].to_pylist()))})
    elif a.cmd == "pack-rsid":
        blocks = write_rsid_index(read_table(a.table), a.out, block_records=a.block_records)
        if a.blocks:
            write_table(blocks, a.blocks)
        _print_json({**read_header(a.out), "blocks": blocks.num_rows})
    elif a.cmd == "pack-variant-index":
        vd, hd = Path(a.variants), Path(a.hits)
        variants = {c: vd / f"{c}.qbv" for c in packfmt.VARIANT_CHROMS if (vd / f"{c}.qbv").is_file()}
        hits = {c: hd / f"{c}.qbh" for c in packfmt.VARIANT_CHROMS if (hd / f"{c}.qbh").is_file()}
        write_variant_index(a.out, variants, hits, a.rsid_index, level=a.level)
        _print_json({**read_header(a.out), "chromosomes": len(variants)})
    elif a.cmd == "check":
        _print_json(check_file(a.pack, variants=a.variants, index=a.index, search_index=a.search_index, dof=a.dof,
                               variant_index=a.variant_index))
    elif a.cmd == "pack-variants":
        tr = read_table(a.trans_only) if a.trans_only else None
        pages = write_variants_pack(read_table(a.table), a.out, a.chrom, page_size=a.page_size, codec=a.codec, level=a.level, trans_only=tr)
        if a.pages:
            write_table(pages, a.pages)
        _print_json({**read_header(a.out), "pages": pages.num_rows})
    elif a.cmd == "pack-results":
        kind = _kind_from_out(a.out, a.kind)
        details = json.loads(Path(a.details).read_text()) if a.details else None
        ptr = write_results_pack(read_table(a.table), a.variants, a.out, kind, level=a.level, details=details)
        if a.pointers:
            write_table(ptr, a.pointers)
        _print_json({**read_header(a.out), "blocks": ptr.num_rows, "rows": int(sum(ptr["n_var"].to_pylist())),
                     "credible_set_records": int(sum(ptr["n_cs"].to_pylist()))})
    elif a.cmd == "pack-trans":
        ptr = write_trans_pack(read_table(a.table), a.out, a.chrom, level=a.level)
        if a.pointers:
            write_table(ptr, a.pointers)
        _print_json({**read_header(a.out), "frames": ptr.num_rows, "eqtl_rows": int(sum(ptr["n_e"].to_pylist())),
                     "sqtl_rows": int(sum(ptr["n_s"].to_pylist()))})
    elif a.cmd == "pack-gwas":
        blocks = write_gwas_packs(read_table(a.table), a.out, block_rows=a.block_rows, level=a.level)
        if a.blocks:
            write_table(blocks, a.blocks)
        _print_json({"dir": a.out, "files": sorted(p.name for p in Path(a.out).iterdir() if p.suffix in (".qbg", ".bin")),
                     "blocks": blocks.num_rows, "index": read_header(Path(a.out) / "gwas_index.bin")})
    return 0


# ---- hits pack (kind 7) -------------------------------------------------------------------------------
def walk_hits_frames(buf: bytes, *, limit: int | None = None) -> list[tuple[int, int]]:
    """(off, len) of the zstd frames after the header of a hits pack, from the frame boundaries alone.
    The browser finds a frame through `variant_index.hits_off`; this is for tooling without it."""
    h = packfmt.parse_file_header(buf)
    if h["kind"] != packfmt.KIND_HITS:
        raise PackToolError(f"hits pack: header kind {h['kind']}")
    out, off, dctx = [], packfmt.FILE_HEADER_LEN, zstandard.ZstdDecompressor()
    while off < len(buf) and (limit is None or len(out) < limit):
        dobj, fed = dctx.decompressobj(), 0
        try:
            while not dobj.eof:
                chunk = buf[off + fed:off + fed + _ZSTD_CHUNK]
                if not chunk:
                    raise PackToolError(f"hits pack: frame at byte {off} runs past the end of the file")
                dobj.decompress(chunk)
                fed += len(chunk)
        except zstandard.ZstdError as e:
            raise PackToolError(f"hits pack: frame at byte {off}: {e}") from None
        ln = fed - len(dobj.unused_data)
        out.append((off, ln))
        off += ln
    if limit is None and len(out) != h["count"]:
        raise PackToolError(f"hits pack: {len(out)} frames, header count {h['count']}")
    return out


def hits_frames(path: str | os.PathLike, *, limit: int | None = None, variant_index: str | os.PathLike | None = None) -> pa.Table:
    """Every frame of a hits pack: `frame`, `first_vidx` (frame number times the header's frame size),
    `hits_off`, `hits_len`. With `variant_index` the offsets come from the startup file instead and
    must match the file (SPEC section 15); without it the frames are walked by their zstd boundaries."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_HITS,), p)
    if variant_index is None:
        fr = walk_hits_frames(p.read_bytes(), limit=limit)
        off = [o for o, _ in fr]
        ln = [n for _, n in fr]
    else:
        c = read_variant_index(variant_index)["chroms"].get(h["chrom"])
        if c is None:
            raise PackToolError(f"{variant_index}: no chromosome {h['chrom']!r}")
        ho = c["hits_off"]
        if len(ho) != h["count"] + 1 or int(ho[-1]) != h["bytes"]:
            raise PackToolError(f"{variant_index}: {len(ho) - 1} frames ending at byte {int(ho[-1])} for {p} "
                                f"(header count {h['count']}, {h['bytes']} bytes)")
        off, ln = ho[:-1].tolist(), np.diff(ho).tolist()
        if limit is not None:
            off, ln = off[:limit], ln[:limit]
    n = len(off)
    return pa.table({"frame": pa.array(range(n), pa.int64()),
                     "first_vidx": pa.array([i * h["page_size"] for i in range(n)], pa.int64()),
                     "hits_off": pa.array(off, pa.int64()), "hits_len": pa.array(ln, pa.int64())})


def read_hits_frame(path: str | os.PathLike, off: int, length: int, first_vidx: int | None = None) -> dict:
    """One decoded frame (`packfmt.decode_hits_frame`) from bytes `[off, off + length)` of a hits pack."""
    _expect_kind(read_header(path), (packfmt.KIND_HITS,), path)
    return packfmt.decode_hits_frame(read_range(path, off, length), first_vidx, what=f"{path} frame at {off}+{length}")


def _gene_rows(search_index: str | os.PathLike | None) -> dict[int, dict]:
    """ord -> {gene_id, symbol, chr, gene_version} from `search_index`, whose row order defines `ord`."""
    if search_index is None:
        return {}
    t = read_search_index_file(search_index, ["gene_id", "symbol", "chr", "tss", "gene_version"])
    rows = sorted(t.to_pylist(), key=lambda r: (r["chr"], r["tss"], r["gene_id"]))
    return {i: r for i, r in enumerate(rows)}


def hits_table(frame: dict, chrom: str, genes: dict[int, dict] | None = None, *, start: int = 0, stop: int | None = None) -> pa.Table:
    """A decoded hits frame as rows `start:stop`: `vidx` (when the frame was decoded with
    `first_vidx`), `kind`, `kind_name`, `qtl_type`, `gene` (the `search_index` ord) and, with
    `genes`, `gene_id`, `symbol` and the rebuilt sQTL `phenotype_id`; then `pval`, `beta`, `slope`,
    `slope_se`, `pip`, `cs_id`, `significant`, and the intron fields. Each value is null on the kinds
    that do not carry it (SPEC section 13)."""
    n = frame["n_rows"]
    stop = n if stop is None else stop
    sl = slice(start, stop)
    kd = frame["kind"][sl]
    ord_ = frame["gene"][sl]
    g = genes or {}
    gid = [g.get(int(o), {}).get("gene_id") for o in ord_.tolist()]
    sym = [g.get(int(o), {}).get("symbol") for o in ord_.tolist()]
    ver = [g.get(int(o), {}).get("gene_version") for o in ord_.tolist()]
    pid = []
    for i, k in enumerate(kd.tolist()):
        if gid[i] is None:
            pid.append(None)
        elif k % 2 == 0:
            pid.append(gid[i])
        else:
            j = start + i
            pid.append(f"{chrom}:{int(frame['intron_start'][j])}:{int(frame['intron_end'][j])}:"
                       f"clu_{int(frame['cluster'][j])}_{frame['strand'][j]}:{gid[i]}.{ver[i]}")
    is_s = (kd % 2) == 1
    cols = {
        "kind": pa.array(kd, pa.int8()),
        "kind_name": pa.array([packfmt.HITS_KIND_NAMES[k] for k in kd.tolist()], pa.string()),
        "qtl_type": pa.array(["s" if s else "e" for s in is_s.tolist()], pa.string()),
        "gene": pa.array(ord_, pa.int64()), "gene_id": pa.array(gid, pa.string()), "symbol": pa.array(sym, pa.string()),
        "phenotype_id": pa.array(pid, pa.string()),
        "pval": _nan_float(frame["pval"][sl]), "beta": _nan_float(frame["beta"][sl]),
        "slope": _nan_float(frame["slope"][sl]), "slope_se": _nan_float(frame["slope_se"][sl]),
        "pip": _nan_float(frame["pip"][sl]),
        "cs_id": pa.array(np.where(frame["cs_id"][sl] < 0, 0, frame["cs_id"][sl]), mask=frame["cs_id"][sl] < 0, type=pa.int64()),
        "significant": pa.array(frame["significant"][sl]),
        "intron_start": pa.array(frame["intron_start"][sl], mask=~is_s, type=pa.int64()),
        "intron_end": pa.array(frame["intron_end"][sl], mask=~is_s, type=pa.int64()),
        "cluster": pa.array(frame["cluster"][sl], mask=~is_s, type=pa.int64()),
        "strand": pa.array(frame["strand"][sl], pa.string()),
    }
    if "vidx" in frame:
        cols = {"vidx": pa.array(frame["vidx"][sl], pa.int64()), **cols}
    return pa.table(cols)


def hits_rows(path: str | os.PathLike, *, vidx: int | None = None, frame: int | None = None,
              variant_index: str | os.PathLike | None = None, search_index: str | os.PathLike | None = None) -> pa.Table:
    """One variant's rows (`vidx`) or a whole frame (`frame`) of a hits pack as a table. The frame is
    found through `variant_index` when given, otherwise by walking the file's zstd frames."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_HITS,), p)
    if (vidx is None) == (frame is None):
        raise PackToolError("give exactly one of --vidx and --frame")
    fv = h["page_size"]
    if vidx is not None:
        if vidx < 0:
            raise PackToolError(f"vidx {vidx} is negative")
        frame = vidx // fv
    if variant_index is not None:
        c = read_variant_index(variant_index)["chroms"].get(h["chrom"])
        if c is None or frame + 1 >= len(c["hits_off"]):
            raise PackToolError(f"{variant_index}: no frame {frame} for {h['chrom']}")
        off, ln = int(c["hits_off"][frame]), int(c["hits_off"][frame + 1] - c["hits_off"][frame])
    else:
        fr = walk_hits_frames(p.read_bytes(), limit=frame + 1)
        if frame >= len(fr):
            raise PackToolError(f"{p}: frame {frame} is outside the file's {len(fr)} frames")
        off, ln = fr[frame]
    f = read_hits_frame(p, off, ln, frame * fv)
    genes = _gene_rows(search_index)
    if vidx is None:
        return hits_table(f, h["chrom"], genes)
    a, b = packfmt.hits_slice(f, vidx)
    return hits_table(f, h["chrom"], genes, start=a, stop=b)


HITS_COLUMNS = ("vidx", "kind", "gene")


def write_hits_pack(rows: pa.Table, out: str | os.PathLike, chrom: str, *, n_variants: int | None = None,
                    frame_variants: int = packfmt.HITS_FRAME_VARIANTS, level: int = 19) -> pa.Table:
    """Hit rows to a kind 7 pack, one zstd frame per `frame_variants` variant indices.

    Columns: `vidx`, `kind` (0 trans eQTL, 1 trans sQTL, 2 lead eQTL, 3 lead sQTL, 4 credible set
    eQTL, 5 credible set sQTL), `gene` (the row's `search_index` ord); then per kind `pval` and
    `beta` (0-1), `pval`, `slope_se`, `slope` and `significant` (2-3), `pip` and `cs_id` (4-5). The
    odd kinds need their intron, as `phenotype_id` or as `intron_start`, `intron_end`, `cluster`,
    `strand`. `n_variants` is the chromosome's variant count, which fixes the frame count (default:
    one past the largest `vidx`). Rows may come in any order. Returns a pointer table: `frame`,
    `first_vidx`, `n_variants`, `rows`, `hits_off`, `hits_len`."""
    for c in HITS_COLUMNS:
        _column(rows, c)
    n = rows.num_rows
    vidx = _ints(_column(rows, "vidx"), "vidx")
    kind = _ints(_column(rows, "kind"), "kind")
    gene = _ints(_column(rows, "gene"), "gene")
    total = int(n_variants if n_variants is not None else (vidx.max() + 1 if n else 0))
    if total <= 0:
        raise PackToolError("a hits pack needs at least one variant: give n_variants")
    if n and (vidx.min() < 0 or vidx.max() >= total):
        raise PackToolError(f"vidx {int(vidx.min())}..{int(vidx.max())} is outside 0..{total - 1}")
    is_s = (kind % 2) == 1
    cols = {c: _floats(_column(rows, c, required=False, default=None)) for c in ("pval", "beta", "slope", "slope_se", "pip")}
    cs_id = _floats(_column(rows, "cs_id", required=False, default=None))
    sig = _column(rows, "significant", required=False, default=False)
    sig = np.array([bool(x) for x in (_strings(sig) if sig.type == pa.string() else sig.to_pylist())], dtype=bool)
    if all(c in rows.column_names for c in ("intron_start", "intron_end", "cluster", "strand")):
        istart, iend, iclu = (_floats(_column(rows, c)) for c in ("intron_start", "intron_end", "cluster"))
        strand = _strings(_column(rows, "strand"))
    elif "phenotype_id" in rows.column_names:
        pid = _strings(_column(rows, "phenotype_id"))
        istart, iend, iclu = (np.full(n, np.nan) for _ in range(3))
        strand = [None] * n
        for i in np.flatnonzero(is_s).tolist():
            if pid[i] is None:
                raise PackToolError(f"row {i}: an sQTL row (kind {int(kind[i])}) needs a phenotype_id")
            istart[i], iend[i], iclu[i], strand[i] = parse_phenotype_id(pid[i])
    elif is_s.any():
        raise PackToolError("kinds 1, 3 and 5 need phenotype_id, or intron_start, intron_end, cluster, strand")
    else:
        istart, iend, iclu, strand = np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan), [None] * n
    n_frames = -(-total // frame_variants)
    frames, ptr = [], {"frame": [], "first_vidx": [], "n_variants": [], "rows": [], "hits_off": [], "hits_len": []}
    order = np.argsort(vidx, kind="stable")
    bounds = np.searchsorted(vidx[order], np.arange(n_frames + 1) * frame_variants)
    off = packfmt.FILE_HEADER_LEN
    for g in range(n_frames):
        first = g * frame_variants
        nv = min(frame_variants, total - first)
        sel = order[bounds[g]:bounds[g + 1]]
        try:
            if sel.size:
                frame = packfmt.encode_hits_frame(
                    first, nv, vidx[sel], kind[sel], gene[sel], pval=cols["pval"][sel], beta=cols["beta"][sel],
                    slope=cols["slope"][sel], slope_se=cols["slope_se"][sel], pip=cols["pip"][sel], cs_id=cs_id[sel],
                    intron_start=istart[sel], intron_end=iend[sel], cluster=iclu[sel],
                    strand=[strand[i] for i in sel.tolist()], significant=sig[sel], level=level)
            else:
                frame = packfmt.encode_hits_frame(first, nv, None, None, None, level=level)
        except ValueError as e:
            raise PackToolError(f"frame {g} (variants {first}..{first + nv - 1}): {e}") from None
        frames.append(frame)
        ptr["frame"].append(g); ptr["first_vidx"].append(first); ptr["n_variants"].append(nv)
        ptr["rows"].append(int(sel.size)); ptr["hits_off"].append(off); ptr["hits_len"].append(len(frame))
        off += len(frame)
    if off > packfmt.U32_MAX:
        raise PackToolError(f"hits pack: {off} bytes exceeds 4 GiB")
    _write_bytes(out, packfmt.file_header(packfmt.KIND_HITS, chrom, n_frames, frame_variants) + b"".join(frames))
    return pa.table({k: pa.array(v, pa.int64()) for k, v in ptr.items()})


# ---- rsID index (kind 8) ------------------------------------------------------------------------------
def rsid_records(path: str | os.PathLike, *, block: int | None = None, limit: int | None = None) -> pa.Table:
    """The rsID index as rows: `rs_number`, `chr`, `vidx`, and `block`. One block with `block`, the
    whole file otherwise (9.2M rows genome-wide, so prefer a block or `rsid_lookup`)."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_RSID,), p)
    B = h["page_size"]
    if block is None:
        buf = read_range(p, packfmt.FILE_HEADER_LEN, h["count"] * packfmt.RSID_RECORD_LEN)
        first = 0
    else:
        off, ln = packfmt.rsid_block_range(block, h["count"], B)
        buf = read_range(p, off, ln)
        first = block * B
    d = packfmt.decode_rsid_block(buf, what=f"{p} block {block}" if block is not None else str(p))
    n = d["rs_number"].size if limit is None else min(limit, d["rs_number"].size)
    return pa.table({"block": pa.array((np.arange(n) + first) // B, pa.int64()),
                     "record": pa.array(np.arange(n) + first, pa.int64()),
                     "rs_number": pa.array(d["rs_number"][:n], pa.int64()),
                     "chr": pa.array(d["chrom"][:n], pa.string()), "vidx": pa.array(d["vidx"][:n], pa.int64())})


def rsid_lookup(path: str | os.PathLike, rs_number: int, variant_index: str | os.PathLike | None = None) -> dict | None:
    """One rsID through the block math (SPEC section 14): {rs_number, chr, vidx, block, byte_start,
    byte_end}, or None when the index does not hold it. With `variant_index` the block is chosen from
    the startup file's samples, as the browser does; without it the samples are read from the file."""
    p = Path(path)
    h = read_header(p)
    _expect_kind(h, (packfmt.KIND_RSID,), p)
    B = h["page_size"]
    n_blocks = -(-h["count"] // B)
    if variant_index is not None:
        first = read_variant_index(variant_index)["rsid_first"]
    else:
        first = np.array([int(np.frombuffer(read_range(p, packfmt.FILE_HEADER_LEN + b * B * packfmt.RSID_RECORD_LEN, 4), "<u4")[0])
                          for b in range(n_blocks)], dtype=np.int64)
    b = packfmt.rsid_block_of(first, rs_number)
    if b is None:
        return None
    off, ln = packfmt.rsid_block_range(b, h["count"], B)
    got = packfmt.rsid_find(packfmt.decode_rsid_block(read_range(p, off, ln), what=f"{p} block {b}"), rs_number)
    if got is None:
        return None
    return {"rs_number": int(rs_number), "chr": got[0], "vidx": got[1], "block": b, "byte_start": off, "byte_end": off + ln}


def write_rsid_index(rows: pa.Table, out: str | os.PathLike, *, block_records: int = packfmt.RSID_BLOCK_RECORDS) -> pa.Table:
    """rsID rows to a kind 8 index. Columns: `rs_number`, `vidx`, and the chromosome as `chr`
    (chr1..chr22, chrX) or `chr_ordinal` (1..23). Rows are sorted here by `rs_number`, which must not
    repeat. Returns a block table: `block`, `first_rs_number`, `byte_start`, `byte_end`."""
    rs = _ints(_column(rows, "rs_number"), "rs_number")
    vidx = _ints(_column(rows, "vidx"), "vidx")
    if "chr_ordinal" in rows.column_names:
        co = _ints(_column(rows, "chr_ordinal"), "chr_ordinal")
    else:
        names = _strings(_column(rows, "chr"))
        try:
            co = np.array([packfmt.VARIANT_CHROMS.index(c) + 1 for c in names], dtype=np.int64)
        except ValueError:
            raise PackToolError("chr must be chr1..chr22 or chrX on every row") from None
    order = np.argsort(rs, kind="stable")
    buf, first = packfmt.encode_rsid_index(rs[order], co[order], vidx[order], block_records=block_records)
    _write_bytes(out, buf)
    starts = packfmt.FILE_HEADER_LEN + np.arange(first.size) * block_records * packfmt.RSID_RECORD_LEN
    ends = np.minimum(starts + block_records * packfmt.RSID_RECORD_LEN, len(buf))
    return pa.table({"block": pa.array(range(first.size), pa.int64()), "first_rs_number": pa.array(first, pa.int64()),
                     "byte_start": pa.array(starts, pa.int64()), "byte_end": pa.array(ends, pa.int64())})


# ---- variant index (kind 9) ---------------------------------------------------------------------------
def read_variant_index(path: str | os.PathLike) -> dict:
    """The startup file decoded (`packfmt.decode_variant_index`)."""
    return packfmt.decode_variant_index(Path(path).read_bytes(), what=str(path))


def variant_index_table(path: str | os.PathLike) -> pa.Table:
    """One row per chromosome of the startup file: `chr`, `n_cis`, `n_trans_only`, `n_pages_cis`,
    `n_pages_trans`, `n_frames`, `variants_bytes` (the last page offset), `hits_bytes`."""
    idx = read_variant_index(path)
    rows = {k: [] for k in ("chr", "n_cis", "n_trans_only", "n_pages_cis", "n_pages_trans", "n_frames", "variants_bytes", "hits_bytes")}
    for name, c in idx["chroms"].items():
        rows["chr"].append(name)
        for k in ("n_cis", "n_trans_only", "n_pages_cis", "n_pages_trans", "n_frames"):
            rows[k].append(c[k])
        rows["variants_bytes"].append(int(c["page_off"][-1]))
        rows["hits_bytes"].append(int(c["hits_off"][-1]))
    return pa.table({k: pa.array(v, pa.string() if k == "chr" else pa.int64()) for k, v in rows.items()})


def write_variant_index(out: str | os.PathLike, variants: dict, hits: dict, rsid_index: str | os.PathLike,
                        *, level: int = 19) -> pa.Table:
    """The kind 9 startup file from the files it indexes: `variants` and `hits` map each chromosome
    (chr1..chr22, chrX) to its variants file and hits pack, and `rsid_index` is the kind 8 file. Each
    variants file is walked for its page offsets and first positions, each hits pack for its frame
    offsets, and the rsID blocks for their first `rs_number`. Returns `variant_index_table`."""
    chroms = []
    fv = None
    for name in packfmt.VARIANT_CHROMS:
        for d, what in ((variants, "variants file"), (hits, "hits pack")):
            if name not in d:
                raise PackToolError(f"variant index: no {what} for {name}")
        vh, pages, first = packfmt.variants_page_firsts(Path(variants[name]).read_bytes())
        hh = read_header(hits[name])
        _expect_kind(hh, (packfmt.KIND_HITS,), hits[name])
        if fv is None:
            fv = hh["page_size"]
        elif hh["page_size"] != fv:
            raise PackToolError(f"variant index: {name} has {hh['page_size']} variants per frame, {fv} elsewhere")
        fr = walk_hits_frames(Path(hits[name]).read_bytes())
        chroms.append({"n_cis": vh["n_cis"], "n_trans_only": vh["count"] - vh["n_cis"],
                       "page_off": [p["offset"] for p in pages] + [pages[-1]["offset"] + pages[-1]["length"]] if pages else [packfmt.FILE_HEADER_LEN],
                       "page_first_position": first,
                       "hits_off": [o for o, _ in fr] + [fr[-1][0] + fr[-1][1]] if fr else [packfmt.FILE_HEADER_LEN]})
    rh = read_header(rsid_index)
    _expect_kind(rh, (packfmt.KIND_RSID,), rsid_index)
    B = rh["page_size"]
    n_blocks = -(-rh["count"] // B)
    rf = [int(np.frombuffer(read_range(rsid_index, packfmt.FILE_HEADER_LEN + b * B * packfmt.RSID_RECORD_LEN, 4), "<u4")[0])
          for b in range(n_blocks)]
    buf = packfmt.encode_variant_index(chroms, int(_variants_page_size(variants)), int(fv), rf, rh["count"], B, level=level)
    _write_bytes(out, buf)
    return variant_index_table(out)


def _variants_page_size(variants: dict) -> int:
    """The variants files' page size, which must be the same for every chromosome."""
    sizes = {read_header(p)["page_size"] for p in variants.values()}
    if len(sizes) != 1:
        raise PackToolError(f"variant index: the variants files disagree on the page size {sorted(sizes)}")
    return sizes.pop()

if __name__ == "__main__":
    sys.exit(main())
