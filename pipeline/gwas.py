"""The GWAS object of an experiment (qtlb v1, SPEC.md section 10).

Reads the contract `gwas` table (CONTRACT.md) from an experiment's tables directory and writes

      <digest>.qbg                    one per chromosome: v1 header (kind 4) + zstd blocks of rows
      <digest>.qgi                    the GWAS index (kind 5): n table, per chromosome block first
                                      positions and end offsets, so a window is one range request
      <digest>.arrow.zst              the bin summary, when the tables hold `gwas_bins.parquet`: the
                                      strongest p per 5 Mb window, one row per window (the landing track)

and returns the experiment's `gwas` entry. `results.build` calls it when the table exists.

The rows are the source's, lossless at its printed precision (4 decimals for beta, se and af; 4
significant digits for p): the v0 GWAS block layout with the alleles in the store's orientation. `ref`
is the reference base(s), `alt` the other allele, `beta` and `af` describe ALT. The adapter oriented
them (negating beta and mirroring af where the source's effect allele was the reference, both exact at
4 decimals) and dropped, with a count, rows whose alleles the reference does not read.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import packfmt_v1 as pf
from . import qtlstore as qs

KIND_GWAS, KIND_GWAS_INDEX = qs.KIND_GWAS, qs.KIND_GWAS_INDEX
EXT_GWAS, EXT_GWAS_INDEX = "qbg", "qgi"
BLOCK_ROWS = 2048
ZSTD_LEVEL = 19
COLUMNS = ("chr", "pos", "ref", "alt", "beta", "se", "af", "pvalue", "n", "rs_number")
ORDER = ["pos", "ref", "alt", "rs_number", "pvalue", "n", "beta"]
EXT_BINS = "arrow.zst"
# the bin summary's schema (SPEC.md section 10, "Bin summary"): v0's gwas_dcm_bins.json columns, typed
BINS_SCHEMA = pa.schema([("chr", pa.string()), ("bin_start", pa.uint32()), ("bin_end", pa.uint32()),
                         ("min_p", pa.float64()), ("lead_position", pa.uint32()), ("lead_rsid", pa.string()),
                         ("lead_beta", pa.float64()), ("lead_ea", pa.string()), ("n_gws", pa.uint32()),
                         ("n_variants", pa.uint32())])


def build(store: qs.Store, tables: Path, cdoc: dict, chroms: list[str], block_rows: int = BLOCK_ROWS,
          level: int = ZSTD_LEVEL) -> dict | None:
    """The `gwas` entry of an experiment, or None when `tables` has no `gwas.parquet`."""
    path = Path(tables) / "gwas.parquet"
    if not path.exists():
        return None
    meta_p = Path(tables) / "gwas.json"
    meta = json.loads(meta_p.read_text()) if meta_p.exists() else {}
    t = pq.read_table(path, filters=[("chr", "in", list(chroms))]).to_pandas()
    missing = [c for c in COLUMNS if c not in t.columns]
    if missing:
        raise ValueError(f"gwas: required columns missing: {missing}")
    table = {c["name"]: c for c in cdoc["chromosomes"]}
    stray = sorted(set(t["chr"]) - set(table))
    if stray:
        raise ValueError(f"gwas: chromosomes {stray} are not in the variant catalog's table")
    # A source may report one site twice (the DCM GWAS lists 796,531 indels once per allele order, with
    # different N and statistics: two measurements, not one record). Both rows are kept and counted.
    if t.duplicated(COLUMNS).any():
        raise ValueError("gwas: identical rows repeat")
    sharing = int(t.duplicated(["chr", "pos", "ref", "alt"], keep=False).sum())
    for c in ("beta", "se", "af", "pvalue", "n"):
        if t[c].isna().any():
            raise ValueError(f"gwas: null {c}")
    n_values = sorted(int(x) for x in t["n"].unique())
    files, index_chroms, rows_by_chrom = {}, [], {}
    for c in [x for x in chroms if x in set(t["chr"])]:
        d = t[t["chr"] == c].sort_values(ORDER, kind="mergesort")
        pos = d["pos"].to_numpy(dtype=np.int64)
        if pos.min() < 1 or pos.max() > table[c]["length"]:
            raise ValueError(f"gwas {c}: position outside 1..{table[c]['length']}")
        rs = np.where(d["rs_number"].to_numpy(dtype=np.int64) < 0, 0, d["rs_number"].to_numpy(dtype=np.int64))
        try:
            codes = pf.gwas_codes(pos, d["beta"].to_numpy(), d["se"].to_numpy(), d["af"].to_numpy(),
                                  d["pvalue"].to_numpy(), rs, d["n"].to_numpy(dtype=np.int64), n_values)
        except ValueError as e:
            raise ValueError(f"gwas {c}: {e}") from None
        ref, alt = d["ref"].tolist(), d["alt"].tolist()
        blocks, fp, eo, off = [], [], [], qs.HEADER_LEN
        for s in range(0, len(d), block_rows):
            e = min(s + block_rows, len(d))
            fr = pf.encode_gwas_block({k: v[s:e] for k, v in codes.items()}, ref[s:e], alt[s:e], level)
            blocks.append(fr)
            off += len(fr)
            fp.append(int(pos[s]))
            eo.append(off)
        if off > pf.U32_MAX:
            raise ValueError(f"gwas {c}: over 4 GiB")
        body = qs.file_header(KIND_GWAS, c, len(d), block_rows, 0, table[c]["seq_digest"]) + b"".join(blocks)
        files[c] = store.put(body, EXT_GWAS)
        index_chroms.append((c, np.array(fp), np.array(eo)))
        rows_by_chrom[c] = len(d)
    payload = pf.encode_gwas_index_payload(n_values, index_chroms, qs.HEADER_LEN)
    index = qs.file_header(KIND_GWAS_INDEX, qs.ALL, len(index_chroms), block_rows, 0, cdoc["collection_digest"]) + \
        pf.zstd_frame(payload, level)
    return {"id": meta.get("id"), "title": meta.get("title"), "files": files, "index": store.put(index, EXT_GWAS_INDEX),
            "n_rows": int(sum(rows_by_chrom.values())), "rows_by_chrom": rows_by_chrom,
            "rows_sharing_a_site": sharing, "block_rows": block_rows,
            "n_values": n_values, "bins": build_bins(store, tables, cdoc, [c for c in chroms if c in files], level),
            "orientation": meta.get("orientation"), "source": meta.get("source")}


# ---- the bin summary --------------------------------------------------------------------------
def _round_bins(t):
    """v0's rounding (steps_gwas.run, when it wrote gwas_dcm_bins.json): p to 3 significant digits
    through its shortest decimal text, the lead's beta to 3 decimals (Python `round`)."""
    t = t.copy()
    t["min_p"] = [float(f"{p:.3g}") for p in t["min_p"]]
    t["lead_beta"] = [round(float(b), 3) for b in t["lead_beta"]]
    return t


def build_bins(store: qs.Store, tables: Path, cdoc: dict, chroms: list[str], level: int = ZSTD_LEVEL) -> dict | None:
    """The `gwas.bins` entry: the bin summary object from `gwas_bins.parquet` on `chroms`, or None when
    the tables have no such table. Rows in variant catalog table order, then by `bin_start`."""
    path = Path(tables) / "gwas_bins.parquet"
    if not path.exists() or not chroms:
        return None
    # the build's chromosomes only (as for the GWAS rows): a subset build has a subset variant catalog
    t = pq.read_table(path, filters=[("chr", "in", list(chroms))]).to_pandas()
    missing = [c for c in BINS_SCHEMA.names if c not in t.columns]
    if missing:
        raise ValueError(f"gwas bins: required columns missing: {missing}")
    order = {c["name"]: i for i, c in enumerate(cdoc["chromosomes"])}
    stray = sorted(set(t["chr"]) - set(order))
    if stray:
        raise ValueError(f"gwas bins: chromosomes {stray} are not in the variant catalog's table")
    t = t.assign(_o=t["chr"].map(order)).sort_values(["_o", "bin_start"], kind="mergesort").drop(columns="_o")
    t = _round_bins(t[list(BINS_SCHEMA.names)]).reset_index(drop=True)
    table = pa.Table.from_pandas(t, schema=BINS_SCHEMA, preserve_index=False)
    width = sorted(set((t["bin_end"] - t["bin_start"]).tolist()))
    if len(width) > 1:
        raise ValueError(f"gwas bins: bins of several widths {width}")
    bin_bp = int(width[0]) if width else 0
    fails = check_bins(table, cdoc, bin_bp, len(t))
    if fails:
        raise ValueError(f"gwas bins: {fails[0]}" + (f" (and {len(fails) - 1} more)" if len(fails) > 1 else ""))
    sink = pa.BufferOutputStream()
    with pa.ipc.new_stream(sink, BINS_SCHEMA) as w:
        w.write_table(table.combine_chunks())
    return {"file": store.put(pf.zstd_frame(sink.getvalue().to_pybytes(), level), EXT_BINS),
            "bin_bp": bin_bp, "n_bins": len(t)}


def decode_bins(buf: bytes) -> pa.Table:
    """The bin summary table from its object bytes (one zstd frame around an Arrow IPC stream)."""
    return pa.ipc.open_stream(pa.py_buffer(pf.zstd_unframe(buf, None, "gwas bins"))).read_all()


def read_bins(store: qs.Store, doc: dict) -> pa.Table | None:
    """An experiment's GWAS bin summary, or None when it has none."""
    b = (doc.get("gwas") or {}).get("bins")
    return decode_bins((store.immutable / b["file"]).read_bytes()) if b else None


def check_bins(t: pa.Table, cdoc: dict, bin_bp: int, n_bins: int) -> list[str]:
    """What `Store.validate` checks on a bin summary (SPEC.md section 10): the schema; `n_bins` rows;
    chromosomes of the variant catalog, in its table order, bins ascending and unique; every bin
    `[bin_start, bin_start + bin_bp)` aligned to `bin_bp` and inside the chromosome; the lead position
    inside its bin; `0 < min_p <= 1`; `0 <= n_gws <= n_variants`, `n_variants >= 1`; `lead_ea` present."""
    if t.schema != BINS_SCHEMA:
        return [f"schema {t.schema.to_string(show_schema_metadata=False)!r} is not the bin schema"]
    fails = []
    if t.num_rows != n_bins:
        fails.append(f"{t.num_rows} bins, the pointer says {n_bins}")
    if t.num_rows and bin_bp < 1:
        fails.append(f"bin_bp {bin_bp}")
        return fails
    d = t.to_pydict()
    table = {c["name"]: (i, c["length"]) for i, c in enumerate(cdoc["chromosomes"])}
    prev = None
    for i in range(t.num_rows):
        c, s, e = d["chr"][i], d["bin_start"][i], d["bin_end"][i]
        if c not in table:
            fails.append(f"bin {i}: chromosome {c!r} not in the variant catalog table")
            continue
        key = (table[c][0], s)
        where = f"bin {c}:{s}"
        if prev is not None and key <= prev:
            fails.append(f"{where}: out of order or repeated")
        prev = key
        if s % bin_bp or e != s + bin_bp or s >= table[c][1]:
            fails.append(f"{where}: not [k * {bin_bp}, +{bin_bp}) inside the chromosome ({table[c][1]})")
        if not s <= d["lead_position"][i] < e:
            fails.append(f"{where}: lead position {d['lead_position'][i]} outside the bin")
        p = d["min_p"][i]
        if p is None or not 0 < p <= 1:
            fails.append(f"{where}: min_p {p}")
        g, n = d["n_gws"][i], d["n_variants"][i]
        if g is None or n is None or n < 1 or g > n:
            fails.append(f"{where}: n_gws {g}, n_variants {n}")
        if d["lead_beta"][i] is None or not d["lead_ea"][i]:
            fails.append(f"{where}: lead beta or effect allele missing")
    return fails


def decode_index(buf: bytes) -> dict:
    """{block_rows, n_values, chroms: {name: (first_position, end_offset)}} from a `.qgi` object."""
    h = qs.parse_file_header(buf)
    if h["kind"] != KIND_GWAS_INDEX or h["chrom"] != qs.ALL:
        raise ValueError(f"gwas index: header kind {h['kind']} chromosome {h['chrom']!r}")
    out = pf.decode_gwas_index_payload(pf.zstd_unframe(buf[qs.HEADER_LEN:], None, "gwas index"), qs.HEADER_LEN)
    if len(out["chroms"]) != h["count"]:
        raise ValueError(f"gwas index: {len(out['chroms'])} chromosomes, header count {h['count']}")
    return {"block_rows": h["page_size"], **out}


def read_window(store: qs.Store, doc: dict, chrom: str, lo: int, hi: int) -> dict:
    """The GWAS rows with lo <= pos <= hi on `chrom`, through the index (one byte range), as columns."""
    g = doc["gwas"]
    idx = decode_index((store.immutable / g["index"]).read_bytes())
    empty = {"pos": np.zeros(0, np.int64), "ref": [], "alt": [], "beta": np.zeros(0), "se": np.zeros(0),
             "af": np.zeros(0), "p": np.zeros(0), "n": np.zeros(0, np.int64), "rs_number": np.zeros(0, np.int64)}
    if chrom not in idx["chroms"]:
        return empty
    fp, eo = idx["chroms"][chrom]
    w = pf.gwas_window(fp, eo, lo, hi, qs.HEADER_LEN)
    if w is None:
        return empty
    start, end, a, b = w
    with open(store.immutable / g["files"][chrom], "rb") as f:
        f.seek(a)
        buf = f.read(b - a)
    cols = {k: [] for k in empty}
    base = a
    for k in range(start, end + 1):
        s0 = (qs.HEADER_LEN if k == 0 else int(eo[k - 1])) - base
        r = pf.decode_gwas_block(buf[s0:int(eo[k]) - base], idx["n_values"])
        keep = (r["position"] >= lo) & (r["position"] <= hi)
        cols["pos"].append(r["position"][keep])
        for c in ("beta", "se", "af", "p", "n", "rs_number"):
            cols[c].append(np.asarray(r[c])[keep])
        cols["ref"] += [x for x, m in zip(r["ref"], keep) if m]
        cols["alt"] += [x for x, m in zip(r["alt"], keep) if m]
    return {k: (v if k in ("ref", "alt") else np.concatenate(v)) for k, v in cols.items()}
