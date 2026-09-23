"""Steps pack_hits and pack_variant_index (SPEC.md sections 13 to 15), and their checks.

`pack_hits` publishes `hits/<chr>` for chr1..chr22 and chrX: one zstd frame per
`packs.hits_frame_variants` variant indices, covering both sections of the chromosome's variants
file. A frame holds every hit row of its variants -- trans eQTL and sQTL rows, the phenotypes the
variant is the lead of, and its credible-set memberships -- each carrying a `search_index` `ord`
instead of a gene id. The variant page reads one frame, so all three of its short lists cost one
request. Each frame's `hits_off`/`hits_len` goes to `_tmp/pack_pointers/hits_<chr>.parquet`.

`pack_variant_index` publishes `rsid_index` (kind 8), the uncompressed 8-byte records that
turn an rs number into a (chromosome, vidx), and `variant_index` (kind 9), the startup
file holding every variants-file page offset and first position, every hits frame offset, and each
rsID block's first rs number. Both are rebuilt whenever a variants file or a hits pack changes.

Both steps rest on `_tmp/variant_idx.parquet`, the (chr, position) -> vidx map. `vidx` is defined by
the variants file's own order -- the cis section by (position, A1, A2), then the trans-only section
by position -- and the intermediate is checked against the built files record by record before
anything uses it.

`validate` reads the three files with its own readers, written from SPEC rather than from `packfmt`,
and writes reference files for `npm run pack-check` to `data/derived/_tmp/pack_check/variant/`.
"""
from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from . import packfmt_v0 as packfmt
from . import packtool as pt
from .common import (CHROMS, Config, connect, log, publish_file, read_search_index, register_search_index, stage,
                     variants_sql, write_parquet)
from .common import search_index_path
from .steps_pack import (EXT, PackError, _raw_sqtl, _require, _write_arrow, pack_file, pack_paths, pointer_dir,
                         sqtl_path)

# Row counts every source must contribute, measured on the 2026-09-14 derived tables.
EXPECT = {"variants": 9_215_026, "rsid": 9_214_133, "rows": 16_424_482,
          "trans_eqtl": 2_680_117, "trans_sqtl": 13_182_408,
          "lead_eqtl": 19_423, "lead_sqtl": 80_750, "cs_eqtl": 215_462, "cs_sqtl": 246_322}
KIND_COUNT_KEY = ("trans_eqtl", "trans_sqtl", "lead_eqtl", "lead_sqtl", "cs_eqtl", "cs_sqtl")

# The four intron fields of an sQTL phenotype id, chr:start:end:clu_<n>_<strand>:<gene>.<version>.
_INTRON = """split_part({p}, ':', 2)::UINTEGER AS intron_start,
       split_part({p}, ':', 3)::UINTEGER AS intron_end,
       regexp_extract(split_part({p}, ':', 4), '^clu_([0-9]+)_[+-]$', 1)::UINTEGER AS cluster,
       right(split_part({p}, ':', 4), 1) AS strand"""


# ---- sources: the one place the variant-page stage repoints when the tables move ------------------
def variants_table(cfg: Config, chrom: str | None = None) -> str:
    return variants_sql(cfg, chrom)


def genes_table(cfg: Config) -> str:
    return f"'{cfg.tables / 'genes.parquet'}'"


def introns_table(cfg: Config) -> str:
    return f"'{cfg.tables / 'splice_phenotypes.parquet'}'"


def credible_sets_table(cfg: Config) -> str:
    return f"'{cfg.tables / 'credible_sets.parquet'}'"


def trans_source(cfg: Config, chrom: str) -> Path:
    """One variant chromosome's trans rows: the variant-keyed copy the `trans` step writes."""
    return cfg.tables / "trans_by_variant" / f"chr={chrom}" / "data.parquet"


# ---- outputs ---------------------------------------------------------------------------------------
def hits_path(cfg: Config, chrom: str) -> Path:
    return pack_file(cfg, f"hits/{chrom}", EXT["hits"])


def rsid_index_path(cfg: Config) -> Path:
    return pack_file(cfg, "rsid_index", EXT["rsid_index"])


def variant_index_path(cfg: Config) -> Path:
    return pack_file(cfg, "variant_index", EXT["variant_index"])


def variant_idx_path(cfg: Config) -> Path:
    return cfg.tmp / "variant_idx.parquet"


def _frame_variants(cfg: Config) -> int:
    return int(cfg["packs"]["hits_frame_variants"])


def _block_records(cfg: Config) -> int:
    return int(cfg["packs"]["rsid_block_records"])


# ---- build: the (chr, position) -> vidx intermediate ------------------------------------------------
def _variant_idx(cfg: Config, con) -> dict[str, dict]:
    """`_tmp/variant_idx.parquet`: every variant's index in its chromosome's variants file.

    `vidx` follows the file's own order (SPEC section 4): the cis section by (position, A1, A2), then
    the trans-only section by position. Every record of every built file is decoded and compared, so
    a wrong ordering here cannot reach a hits pack or the rsID index."""
    out = variant_idx_path(cfg)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    t0 = time.time()
    con.execute(f"""COPY (
        SELECT chr, position, rs_number, in_cis,
               (row_number() OVER (PARTITION BY chr ORDER BY NOT in_cis, position, A1, A2) - 1)::UINTEGER AS vidx
        FROM {variants_table(cfg)}
    ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000)""")
    os.replace(tmp, out)
    total, distinct = con.execute(f"SELECT count(*), count(DISTINCT (chr, position)) FROM '{out}'").fetchone()
    _require(total == distinct == EXPECT["variants"],
             f"variant_idx: {total:,} rows, {distinct:,} distinct (chr, position), expected {EXPECT['variants']:,}")
    stats = {}
    for chrom in CHROMS:
        n, n_cis = con.execute(f"SELECT count(*), count(*) FILTER (WHERE in_cis) FROM '{out}' WHERE chr = ?", [chrom]).fetchone()
        h = pt.read_header(pack_paths(cfg, chrom)[0])
        _require(n == h["count"] and n_cis == h["n_cis"],
                 f"variant_idx {chrom}: {n:,} variants ({n_cis:,} cis); the variants file header says {h['count']:,} and {h['n_cis']:,}")
        # every record of the built file, against the order this intermediate defines
        src = con.execute(f"SELECT position FROM '{out}' WHERE chr = ? ORDER BY vidx", [chrom]).fetch_arrow_table()
        dec = pa.concat_tables([pt.variant_rows(pack_paths(cfg, chrom)[0], section=s) for s in ("cis", "trans")])
        _require(dec.num_rows == n and dec["vidx"].to_pylist() == list(range(n)),
                 f"variant_idx {chrom}: the variants file decodes to {dec.num_rows:,} records, not {n:,} in vidx order")
        _require(np.array_equal(dec["position"].to_numpy(), src["position"].to_numpy()),
                 f"variant_idx {chrom}: a decoded position differs from the position at that vidx")
        stats[chrom] = {"variants": n, "n_cis": n_cis, "n_trans_only": n - n_cis}
    log(f"variant_idx: {total:,} variants over {len(CHROMS)} chromosomes, every position checked against its variants file, "
        f"{time.time() - t0:.1f} s")
    return stats


# ---- build: the hits pack ---------------------------------------------------------------------------
def _hit_rows(cfg: Config, con, chrom: str, sample=None) -> pa.Table:
    """Every hit row of one variant chromosome, joined to its `vidx` and to `search_index`'s `ord`.

    With `sample`, only those variant indices: validate compares a sample without rebuilding the
    whole chromosome, and the rows it compares against then come from this very query."""
    src = trans_source(cfg, chrom)
    _require(src.exists(), f"pack_hits {chrom}: {src} is missing; run `build --step trans`")
    if sample is None:
        con.execute(f"CREATE OR REPLACE TABLE vi AS SELECT position, vidx FROM '{variant_idx_path(cfg)}' WHERE chr = '{chrom}'")
        stray = con.execute(f"SELECT count(*) FROM read_parquet('{src}', hive_partitioning = false) WHERE variant_chr <> '{chrom}'").fetchone()[0]
        _require(stray == 0, f"pack_hits {chrom}: {stray:,} rows of {src.name} sit on another variant chromosome")
    else:
        con.register("vi_sample", pa.table({"vidx": pa.array([int(v) for v in sample], pa.int64())}))
        con.execute(f"""CREATE OR REPLACE TABLE vi AS SELECT position, vidx FROM '{variant_idx_path(cfg)}'
            WHERE chr = '{chrom}' AND vidx IN (SELECT vidx FROM vi_sample)""")
    return con.execute(f"""
        SELECT vi.vidx, CASE WHEN t.qtl_type = 's' THEN 1 ELSE 0 END::UTINYINT AS kind, si.ord AS gene,
               t.pval::DOUBLE AS pval, t.beta::DOUBLE AS beta, NULL::DOUBLE AS slope, NULL::DOUBLE AS slope_se,
               NULL::DOUBLE AS pip, NULL::USMALLINT AS cs_id, false AS significant,
               CASE WHEN t.qtl_type = 's' THEN split_part(t.phenotype_id, ':', 2)::UINTEGER END AS intron_start,
               CASE WHEN t.qtl_type = 's' THEN split_part(t.phenotype_id, ':', 3)::UINTEGER END AS intron_end,
               CASE WHEN t.qtl_type = 's' THEN regexp_extract(split_part(t.phenotype_id, ':', 4), '^clu_([0-9]+)_[+-]$', 1)::UINTEGER END AS cluster,
               CASE WHEN t.qtl_type = 's' THEN right(split_part(t.phenotype_id, ':', 4), 1) END AS strand
        FROM read_parquet('{src}', hive_partitioning = false) t
        JOIN si USING (gene_id)
        JOIN vi ON vi.position = t.position
      UNION ALL
        SELECT vi.vidx, 2::UTINYINT, si.ord, g.pval_perm::DOUBLE, NULL::DOUBLE, g.slope::DOUBLE, g.slope_se::DOUBLE,
               NULL::DOUBLE, NULL::USMALLINT, g.is_egene,
               NULL::UINTEGER, NULL::UINTEGER, NULL::UINTEGER, NULL::VARCHAR
        FROM {genes_table(cfg)} g JOIN si USING (gene_id) JOIN vi ON vi.position = g.lead_position
        WHERE g.tested AND g.chr = '{chrom}'
      UNION ALL
        SELECT vi.vidx, 3::UTINYINT, si.ord, s.pval_perm::DOUBLE, NULL::DOUBLE, s.slope::DOUBLE, s.slope_se::DOUBLE,
               NULL::DOUBLE, NULL::USMALLINT, s.is_sqtl,
               {_INTRON.format(p='s.phenotype_id')}
        FROM {introns_table(cfg)} s JOIN si USING (gene_id) JOIN vi ON vi.position = s.lead_position
        WHERE s.chr = '{chrom}'
      UNION ALL
        SELECT vi.vidx, CASE WHEN c.qtl_type = 's' THEN 5 ELSE 4 END::UTINYINT, si.ord,
               NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, NULL::DOUBLE, c.pip::DOUBLE, c.cs_id::USMALLINT, false,
               CASE WHEN c.qtl_type = 's' THEN split_part(c.phenotype_id, ':', 2)::UINTEGER END,
               CASE WHEN c.qtl_type = 's' THEN split_part(c.phenotype_id, ':', 3)::UINTEGER END,
               CASE WHEN c.qtl_type = 's' THEN regexp_extract(split_part(c.phenotype_id, ':', 4), '^clu_([0-9]+)_[+-]$', 1)::UINTEGER END,
               CASE WHEN c.qtl_type = 's' THEN right(split_part(c.phenotype_id, ':', 4), 1) END
        FROM {credible_sets_table(cfg)} c JOIN si USING (gene_id) JOIN vi ON vi.position = c.position
        WHERE c.chr = '{chrom}'""").fetch_arrow_table()


def _source_counts(cfg: Config, con, chrom: str) -> dict[str, int]:
    """What each source holds for this chromosome, so the join can be shown to lose nothing."""
    tr = con.execute(f"""SELECT count(*) FILTER (WHERE qtl_type = 'e'), count(*) FILTER (WHERE qtl_type = 's')
        FROM read_parquet('{trans_source(cfg, chrom)}', hive_partitioning = false)""").fetchone()
    le = con.execute(f"SELECT count(*) FROM {genes_table(cfg)} WHERE tested AND chr = ?", [chrom]).fetchone()[0]
    ls = con.execute(f"SELECT count(*) FROM {introns_table(cfg)} WHERE chr = ?", [chrom]).fetchone()[0]
    cs = con.execute(f"""SELECT count(*) FILTER (WHERE qtl_type = 'e'), count(*) FILTER (WHERE qtl_type = 's')
        FROM {credible_sets_table(cfg)} WHERE chr = ?""", [chrom]).fetchone()
    return dict(zip(KIND_COUNT_KEY, (tr[0], tr[1], le, ls, cs[0], cs[1])))


def pack_hits(cfg: Config) -> None:
    d, pk = cfg.derived, cfg["packs"]
    level, F = int(pk["zstd_level"]), _frame_variants(cfg)
    con = connect(cfg, memory_limit=pk["duckdb_memory_limit"], threads=int(pk["duckdb_threads"]),
                  temp_dir=cfg.tmp / "pack_hits")
    idx = register_search_index(cfg, con, "si_all")
    _require("ord" in read_search_index(cfg).column_names,
             "pack_hits: search_index has no ord column; run `build --step search_index --force`")
    con.execute(f"CREATE TABLE si AS SELECT gene_id, ord FROM {idx}")
    n_ord, max_ord = con.execute("SELECT count(*), max(ord) FROM si").fetchone()
    _require(max_ord < 65536, f"pack_hits: search_index ord reaches {max_ord}, which no longer fits a u16")
    _variant_idx(cfg, con)
    pdir = pointer_dir(cfg)
    pdir.mkdir(parents=True, exist_ok=True)
    for chrom in CHROMS:
        t0 = time.time()
        key = f"hits/{chrom}"
        want = _source_counts(cfg, con, chrom)
        tb = _hit_rows(cfg, con, chrom)
        got = dict(zip(KIND_COUNT_KEY, np.bincount(tb["kind"].to_numpy(zero_copy_only=False).astype(np.int64),
                                                   minlength=6).tolist()))
        _require(got == want, f"pack_hits {chrom}: rows by kind after the joins {got} differ from the sources {want}; "
                              f"a source row did not join to a variant or to a gene")
        h = pt.read_header(pack_paths(cfg, chrom)[0])
        total = h["count"]
        tmp = stage(cfg, key, EXT["hits"])
        ptr = pt.write_hits_pack(tb, tmp, chrom, n_variants=total, frame_variants=F, level=level)
        out = publish_file(cfg, tmp, key, EXT["hits"])
        write_parquet(ptr, pdir / f"hits_{chrom}.parquet", 100_000)
        size = out.stat().st_size
        lens = np.array(ptr["hits_len"].to_pylist())
        rows = np.array(ptr["rows"].to_pylist())
        _require(size == packfmt.FILE_HEADER_LEN + int(lens.sum()) and size <= packfmt.U32_MAX,
                 f"pack_hits {chrom}: {size:,} bytes is not 32 + sum(hits_len) ({int(lens.sum()):,}), or is past 4 GiB")
        st = {"chrom": chrom, "frames": len(lens), "variants": total, "rows": int(rows.sum()),
              "by_kind": {k: int(got[k]) for k in KIND_COUNT_KEY}, "bytes": size,
              "source_bytes": trans_source(cfg, chrom).stat().st_size,
              "median_frame": int(np.median(lens)), "max_frame": int(lens.max()),
              "max_frame_first_vidx": int(ptr["first_vidx"][int(lens.argmax())].as_py()),
              "max_rows_in_a_frame": int(rows.max()), "empty_frames": int((rows == 0).sum()),
              "seconds": round(time.time() - t0, 1)}
        (pdir / f"hits_{chrom}.json").write_text(json.dumps(st, indent=2))
        log(f"pack_hits {chrom}: {st['frames']:,} frames over {total:,} variants, {st['rows']:,} rows, {size:,} B "
            f"({size / max(st['rows'], 1):.2f} B/row; trans parquet {st['source_bytes']:,} B); frame median "
            f"{st['median_frame']:,} B, max {st['max_frame']:,} B; {st['empty_frames']:,} empty; {st['seconds']} s")
        del tb
    stats = [json.loads((pdir / f"hits_{c}.json").read_text()) for c in CHROMS]
    tot = {k: sum(s[k] for s in stats) for k in ("frames", "variants", "rows", "bytes")}
    by_kind = {k: sum(s["by_kind"][k] for s in stats) for k in KIND_COUNT_KEY}
    for k in KIND_COUNT_KEY:
        log(f"pack_hits total {k}: {by_kind[k]:,} (expected {EXPECT[k]:,}{'' if by_kind[k] == EXPECT[k] else ', DIFFERS'})")
    _require(tot["rows"] == EXPECT["rows"] and tot["variants"] == EXPECT["variants"],
             f"pack_hits: {tot['rows']:,} rows over {tot['variants']:,} variants, expected {EXPECT['rows']:,} and {EXPECT['variants']:,}")
    big = max(stats, key=lambda s: s["max_frame"])
    log(f"pack_hits total: {tot['frames']:,} frames over {tot['variants']:,} variants in {len(stats)} files, "
        f"{tot['rows']:,} rows, {tot['bytes']:,} B ({tot['bytes'] / tot['rows']:.2f} B/row); largest frame "
        f"{big['max_frame']:,} B ({big['chrom']} from vidx {big['max_frame_first_vidx']:,})")


# ---- build: the rsID index and the startup file ------------------------------------------------------
def pack_variant_index(cfg: Config) -> None:
    d, pk = cfg.derived, cfg["packs"]
    level, B = int(pk["zstd_level"]), _block_records(cfg)
    con = connect(cfg, memory_limit=pk["duckdb_memory_limit"], threads=int(pk["duckdb_threads"]),
                  temp_dir=cfg.tmp / "pack_variant_index")
    vip = variant_idx_path(cfg)
    _require(vip.exists(), "pack_variant_index: _tmp/variant_idx.parquet is missing; run `build --step pack_hits`")
    pdir = pointer_dir(cfg)
    pdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # 1. rsid_index.qbr: one 8-byte record per rsID, sorted by rs number, uncompressed
    rs = con.execute(f"""SELECT rs_number::BIGINT AS rs_number, chr, vidx::BIGINT AS vidx
        FROM '{vip}' WHERE rs_number IS NOT NULL ORDER BY rs_number""").fetch_arrow_table()
    n_rs = rs.num_rows
    _require(n_rs == EXPECT["rsid"], f"rsid index: {n_rs:,} rsIDs, expected {EXPECT['rsid']:,}")
    a = rs["rs_number"].to_numpy()
    _require(bool(np.all(np.diff(a) > 0)), "rsid index: rs_number does not strictly increase, so an rsID repeats")
    tmp = stage(cfg, "rsid_index", EXT["rsid_index"])
    blocks = pt.write_rsid_index(rs, tmp, block_records=B)
    rsid_out = publish_file(cfg, tmp, "rsid_index", EXT["rsid_index"])
    rsid_bytes = rsid_out.stat().st_size
    _require(rsid_bytes == packfmt.FILE_HEADER_LEN + packfmt.RSID_RECORD_LEN * n_rs,
             f"rsid index: {rsid_bytes:,} bytes for {n_rs:,} records")
    log(f"pack_variant_index: {rsid_out.name} {n_rs:,} records in {blocks.num_rows:,} blocks of {B}, {rsid_bytes:,} B")
    del rs, a

    # 2. variant_index.qbx: page and frame offsets of every file, plus each rsID block's first number
    variants = {c: pack_paths(cfg, c)[0] for c in CHROMS}
    hits = {c: hits_path(cfg, c) for c in CHROMS}
    missing = [str(p) for p in list(variants.values()) + list(hits.values()) if not p.exists()]
    _require(not missing, f"pack_variant_index: missing inputs {missing[:3]}")
    tmp = stage(cfg, "variant_index", EXT["variant_index"])
    table = pt.write_variant_index(tmp, variants, hits, rsid_out, level=level)
    vx_out = publish_file(cfg, tmp, "variant_index", EXT["variant_index"])
    vx_bytes = vx_out.stat().st_size
    raw = len(packfmt.zstd_unframe(vx_out.read_bytes()[packfmt.FILE_HEADER_LEN:], None, "variant_index"))
    idx = pt.read_variant_index(vx_out)
    pages = sum(c["n_pages_cis"] + c["n_pages_trans"] for c in idx["chroms"].values())
    frames = sum(c["n_frames"] for c in idx["chroms"].values())
    variants_n = sum(c["n_cis"] + c["n_trans_only"] for c in idx["chroms"].values())
    _require(variants_n == EXPECT["variants"],
             f"variant_index: {variants_n:,} variants, expected {EXPECT['variants']:,}")
    st = {"rsid_records": n_rs, "rsid_blocks": int(blocks.num_rows), "rsid_block_records": B, "rsid_bytes": rsid_bytes,
          "variants": variants_n, "pages": pages, "frames": frames, "page_size": idx["page_size"],
          "frame_variants": idx["frame_variants"], "raw_bytes": raw, "bytes": vx_bytes,
          "seconds": round(time.time() - t0, 1)}
    (pdir / "variant_index.json").write_text(json.dumps(st, indent=2))
    write_parquet(table, pdir / "variant_index.parquet", 100_000)
    log(f"pack_variant_index: {vx_out.name} {pages:,} pages and {frames:,} frames over {variants_n:,} variants, "
        f"{raw:,} B raw -> {vx_bytes:,} B zstd-{level}; with the rsID index {rsid_bytes + vx_bytes:,} B; {st['seconds']} s")


# ---- validate: an independent reader of SPEC sections 13 to 15 --------------------------------------
SAMPLE_PER_CHROM = 440         # about 10,000 sampled variant indices over the 23 chromosomes
SEED = 20260915
TOL = {"nlp_half_step_divisor": 131066, "nlp_slack": 1e-9, "beta_half_step_divisor": 65534,
       "se_half_step_divisor": 131070, "slope_half_step_divisor": 65534, "pip_abs": 0.5 / 65535 + 1e-12}
EXACT_FIELDS = ("vidx", "kind", "gene", "intron_start", "intron_end", "cluster", "cs_id", "strand", "significant")
SORT_KEY = ("cs_id", "cluster", "intron_end", "intron_start", "gene", "kind", "vidx")   # lexsort: last is primary


def _source_arrays(tb: pa.Table) -> dict:
    """The source rows as plain arrays, with the same nulls-to-sentinels the codec uses."""
    col = lambda n: tb[n].to_numpy(zero_copy_only=False)
    out = {k: col(k).astype(np.float64) for k in ("pval", "beta", "slope", "slope_se", "pip")}
    out["vidx"] = col("vidx").astype(np.int64)
    out["kind"] = col("kind").astype(np.int64)
    out["gene"] = col("gene").astype(np.int64)
    for k in ("intron_start", "intron_end", "cluster"):
        out[k] = np.nan_to_num(col(k).astype(np.float64)).astype(np.int64)
    cs = col("cs_id").astype(np.float64)
    out["cs_id"] = np.where(np.isnan(cs), -1, np.nan_to_num(cs)).astype(np.int64)
    out["strand"] = np.array([x if x is not None else "" for x in tb["strand"].to_pylist()])
    out["significant"] = np.array([bool(x) for x in tb["significant"].to_pylist()])
    return out


def _decoded_arrays(path, vidx_sorted: np.ndarray, frame_variants: int) -> dict:
    """Decode only the frames holding `vidx_sorted`, and return those variants' rows with each row's
    frame scales, so a value can be checked against the scale it was quantized with."""
    buf = Path(path).read_bytes()
    frames = pt.walk_hits_frames(buf)
    parts, scales = [], []
    for g in sorted({int(v) // frame_variants for v in vidx_sorted.tolist()}):
        off, ln = frames[g]
        f = packfmt.decode_hits_frame(buf[off:off + ln], g * frame_variants, what=f"{path} frame {g}")
        want = vidx_sorted[(vidx_sorted >= g * frame_variants) & (vidx_sorted < (g + 1) * frame_variants)]
        for v in want.tolist():
            a, b = packfmt.hits_slice(f, int(v))
            if b == a:
                continue
            parts.append({k: f[k][a:b] for k in ("vidx", "kind", "gene", "intron_start", "intron_end", "cluster",
                                                 "cs_id", "pval", "beta", "slope", "slope_se", "pip", "significant")}
                         | {"strand": np.array([x if x is not None else "" for x in f["strand"][a:b]])})
            scales.append(np.tile([f["trans_nlp_max"], f["trans_beta_max"], f["perm_nlp_max"], f["se_max"],
                                   f["slope_max"]], (b - a, 1)))
    if not parts:
        return {k: np.zeros(0) for k in EXACT_FIELDS}, np.zeros((0, 5))
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    return out, np.concatenate(scales)


def _compare_hits(chrom: str, S: dict, D: dict, sc: np.ndarray, worst: dict) -> list[str]:
    """Exact fields exactly, values within SPEC section 13's limits. Returns the failures."""
    bad = []
    if S["vidx"].size != D["vidx"].size:
        return [f"{chrom}: {S['vidx'].size:,} source rows against {D['vidx'].size:,} decoded"]
    if not S["vidx"].size:
        return bad
    so = np.lexsort(tuple(S[k] for k in SORT_KEY))
    do = np.lexsort(tuple(D[k] for k in SORT_KEY))
    tup = np.stack([S[k][so] for k in SORT_KEY])
    if tup.shape[1] > 1 and bool((np.diff(tup, axis=1) == 0).all(axis=0).any()):
        bad.append(f"{chrom}: two rows share (vidx, kind, gene, intron, cs_id), so the comparison order is ambiguous")
    for k in EXACT_FIELDS:
        n = int((S[k][so] != D[k][do]).sum())
        if n:
            bad.append(f"{chrom}: {n:,} rows differ in {k}")
    kd, s2 = S["kind"][so], sc[do]
    is_tr, is_ld, is_cs = kd <= 1, (kd == 2) | (kd == 3), kd >= 4
    with np.errstate(divide="ignore", invalid="ignore"):
        snlp, dnlp = -np.log10(S["pval"][so]), -np.log10(D["pval"][do])
    for name, sel, a, b, lim in (
            ("trans -log10 p", is_tr, snlp, dnlp, s2[:, 0] / TOL["nlp_half_step_divisor"] + TOL["nlp_slack"]),
            ("trans beta", is_tr, S["beta"][so], D["beta"][do], s2[:, 1] / TOL["beta_half_step_divisor"] + 1e-12),
            ("lead -log10 pval_perm", is_ld, snlp, dnlp, s2[:, 2] / TOL["nlp_half_step_divisor"] + TOL["nlp_slack"]),
            ("lead slope_se", is_ld, S["slope_se"][so], D["slope_se"][do], s2[:, 3] / TOL["se_half_step_divisor"] + 1e-12),
            ("lead slope", is_ld, S["slope"][so], D["slope"][do], s2[:, 4] / TOL["slope_half_step_divisor"] + 1e-12),
            ("pip", is_cs, S["pip"][so], D["pip"][do], np.full(kd.size, TOL["pip_abs"]))):
        if not sel.any():
            continue
        err = np.abs(a[sel] - b[sel])
        over = int((err > lim[sel]).sum())
        worst[name] = max(worst.get(name, 0.0), float(err.max()))
        if over:
            bad.append(f"{chrom}: {over:,} {name} values outside SPEC section 13's limit")
    return bad


# ---- validate: the cis scan of SPEC section 7 ("The cis scan") --------------------------------------
def raw_eqtl_path(cfg: Config, chrom: str) -> Path:
    """One chromosome's raw Zenodo eQTL nominal file, the independent reference for the cis scan."""
    return cfg.raw_dir("cis_eQTL_nominal") / f"topchef_{chrom}_MaxPC70.cis_qtl_pairs.{chrom}.parquet"


def _scan_span(buf: bytes, base: int, named: dict[int, str], vidx: int, dof: int,
               details: bool) -> tuple[list[dict], list[tuple[str, int, int]], dict]:
    """Walk one span's 64-byte block headers, exactly as the browser does: a block is read only
    through its header, and only a block that covers `vidx` gives up its one pair and its
    credible-set records. `named` maps a block's absolute offset to the phenotype it belongs to;
    a block inside the span that is not named, or does not cover `vidx`, is skipped.

    `details` decompresses a named block's details frame for its `splice` entries, which is how the
    eQTL span names the intron blocks of the sQTL span. It happens for every named block, not only a
    covering one: a gene that is sQTL-tested but not eQTL-tested has an empty block (`n_rows` 0) that
    can never cover a variant, and its introns would otherwise be missed."""
    rows: list[dict] = []
    introns: list[tuple[str, int, int]] = []
    off = n_blocks = n_skipped = 0
    while off < len(buf):
        h = packfmt.parse_block_header(buf, off)
        n_blocks += 1
        n, vs, name = h["n_rows"], h["var_start"], named.get(off + base)
        if name is not None and details:
            body = packfmt.BLOCK_HEADER_LEN + 4 * n + packfmt.CS_RECORD_LEN * h["n_cs"]
            dj = packfmt.zstd_unframe(bytes(buf[off + body:off + body + h["details_zlen"]]),
                                      h["details_len"], "block details")
            for s in json.loads(dj.decode("utf-8")).get("splice", []):
                introns.append((s["phenotype_id"], int(s["blk_off"]), int(s["blk_len"])))
        if name is not None and n and vs <= vidx < vs + n:
            row = vidx - vs
            pair = np.frombuffer(buf, packfmt.PAIR_DTYPE, 1, off + packfmt.BLOCK_HEADER_LEN + 4 * row)
            nlp = packfmt.dequantize_nlp(pair["nlp"].copy(), h["nlp_max"])
            pval = packfmt.p_from_nlp(nlp)
            se, negative = packfmt.dequantize_se(pair["se"].copy(), h["lse_min"], h["lse_max"])
            cs = np.frombuffer(buf, packfmt.CS_DTYPE, h["n_cs"], off + packfmt.BLOCK_HEADER_LEN + 4 * n)
            mine = cs[cs["row"] == row]
            rows.append({"name": name, "pval_nominal": float(pval[0]), "slope_se": float(se[0]),
                         "slope": float(packfmt.slope_from_se(se, negative, pval, dof)[0]),
                         "nlp_max": h["nlp_max"], "lse_min": h["lse_min"], "lse_max": h["lse_max"],
                         "cs": sorted((int(r["cs_id"]), float(r["pip"])) for r in mine)})
        else:
            n_skipped += 1
        off += h["blk_len"]
    return rows, introns, {"bytes": len(buf), "blocks": n_blocks, "non_covering": n_skipped}


def cis_scan(cfg: Config, con, chrom: str, position: int, vidx: int, dof_e: int,
             dof_s: int) -> tuple[list[dict], list[dict], dict, dict]:
    """SPEC section 7's cis scan over the packs: two ranges, the eQTL span and the sQTL span."""
    genes = con.execute(f"""SELECT gene_id, blk_off, blk_len FROM {register_search_index(cfg, con, 'si_scan')}
        WHERE chr = ? AND w_lo <= ? AND ? <= w_hi AND blk_off IS NOT NULL ORDER BY blk_off""",
                        [chrom, position, position]).fetchall()
    _require(bool(genes), f"cis scan: no {chrom} gene window covers {position:,}")
    epath = pack_paths(cfg, chrom)[1]
    lo, hi = int(genes[0][1]), int(genes[-1][1]) + int(genes[-1][2])
    e_rows, introns, e_span = _scan_span(pt.read_range(epath, lo, hi - lo), lo,
                                         {int(g[1]): g[0] for g in genes}, vidx, dof_e, True)
    e_span |= {"byte_start": lo, "byte_end": hi, "genes_in_window": len(genes)}

    s_rows: list[dict] = []
    s_span = {"bytes": 0, "blocks": 0, "non_covering": 0, "byte_start": 0, "byte_end": 0, "introns_named": 0}
    if introns:
        slo = min(o for _, o, _ in introns)
        shi = max(o + n for _, o, n in introns)
        s_rows, _, s_span = _scan_span(pt.read_range(sqtl_path(cfg, chrom), slo, shi - slo), slo,
                                       {o: pid for pid, o, _ in introns}, vidx, dof_s, False)
        s_span |= {"byte_start": slo, "byte_end": shi, "introns_named": len(introns)}
    return e_rows, s_rows, e_span, s_span


def _scan_table(e_rows: list[dict], s_rows: list[dict], label: str) -> pa.Table:
    """One scan's decoded rows as a table: the eQTL genes then the sQTL introns, each with the block
    scales it was quantized with, so a browser decoder can be checked row by row."""
    kinds = ["eqtl"] * len(e_rows) + ["sqtl"] * len(s_rows)
    rows = e_rows + s_rows
    col = lambda k: [r[k] for r in rows]
    return pa.table({
        "scan": pa.array([label] * len(rows), pa.string()),
        "qtl_type": pa.array(kinds, pa.string()),
        "phenotype_id": pa.array(col("name"), pa.string()),
        "pval_nominal": pa.array(col("pval_nominal"), pa.float64()),
        "slope": pa.array(col("slope"), pa.float64()),
        "slope_se": pa.array(col("slope_se"), pa.float64()),
        "nlp_max": pa.array(col("nlp_max"), pa.float64()),
        "lse_min": pa.array(col("lse_min"), pa.float64()),
        "lse_max": pa.array(col("lse_max"), pa.float64()),
        "cs_id": pa.array([[c[0] for c in r["cs"]] for r in rows], pa.list_(pa.int8())),
        "pip": pa.array([[c[1] for c in r["cs"]] for r in rows], pa.list_(pa.float32())),
    })


def _scan_parity(cfg: Config, con, check, chrom: str, position: int, label: str) -> dict:
    """The scan's rows against the raw Zenodo nominal files at the same position: the same
    phenotypes, and p, slope_se and slope within SPEC section 9's limits. Returns the scan as a
    reference entry for `npm run pack-check`, with its decoded rows under `rows`."""
    m = json.loads((cfg.derived / "manifest.json").read_text())
    dof_e, dof_s = int(m["packs"]["dof"]["eqtl"]), int(m["packs"]["dof"]["sqtl"])
    vidx = con.execute(f"SELECT vidx FROM '{variant_idx_path(cfg)}' WHERE chr = ? AND position = ?",
                       [chrom, position]).fetchone()
    empty = {"id": label, "chrom": chrom, "position": position, "rows": _scan_table([], [], label)}
    if vidx is None:
        check(False, f"cis scan {label}: {chrom}:{position:,} is not a variant")
        return empty
    vidx = int(vidx[0])
    t0 = time.time()
    e_rows, s_rows, e_span, s_span = cis_scan(cfg, con, chrom, position, vidx, dof_e, dof_s)
    log(f"  cis scan {label} ({chrom}:{position:,}, vidx {vidx:,}): eQTL span {e_span['bytes']:,} B, "
        f"{e_span['blocks']} blocks ({e_span['non_covering']} non-covering), {len(e_rows)} genes of "
        f"{e_span['genes_in_window']} in the window; sQTL span {s_span['bytes']:,} B, {s_span['blocks']} blocks "
        f"({s_span['non_covering']} non-covering), {len(s_rows)} of {s_span['introns_named']} named introns; "
        f"{time.time() - t0:.1f} s")

    ref = {"id": label, "chrom": chrom, "position": position, "vidx": vidx,
           "eqtl": {"file": str(pack_paths(cfg, chrom)[1].relative_to(cfg.derived)), **e_span},
           "sqtl": {"file": str(sqtl_path(cfg, chrom).relative_to(cfg.derived)), **s_span},
           "rows": _scan_table(e_rows, s_rows, label)}
    for what, rows, raw, dof in (("eQTL", e_rows, raw_eqtl_path(cfg, chrom), dof_e),
                                 ("sQTL", s_rows, _raw_sqtl(cfg, chrom), dof_s)):
        src = con.execute(f"""SELECT phenotype_id, pval_nominal, slope, slope_se FROM '{raw}'
            WHERE position = ? ORDER BY phenotype_id""", [position]).fetch_arrow_table()
        got = {r["name"]: r for r in rows}
        same = sorted(got) == src["phenotype_id"].to_pylist()
        check(same, f"cis scan {label} {what}: {len(rows):,} phenotypes from the packs are exactly the "
                    f"{src.num_rows:,} in the raw file at {chrom}:{position:,}"
                    + ("" if same else f" (only in the packs {sorted(set(got) - set(src['phenotype_id'].to_pylist()))[:3]}, "
                                       f"only in the raw {sorted(set(src['phenotype_id'].to_pylist()) - set(got))[:3]})"))
        if not same:
            continue
        bad, worst = [], {"-log10 p": 0.0, "slope_se": 0.0, "slope": 0.0}
        for pid, p_src, sl_src, se_src in zip(src["phenotype_id"].to_pylist(), src["pval_nominal"].to_pylist(),
                                              src["slope"].to_pylist(), src["slope_se"].to_pylist()):
            r = got[pid]
            if p_src is None or not (p_src > 0):
                continue
            e_nlp = abs(-np.log10(r["pval_nominal"]) + np.log10(p_src))
            lim_nlp = r["nlp_max"] / TOL["nlp_half_step_divisor"] + TOL["nlp_slack"]
            se_b, sl_b = packfmt.error_bounds(np.array([p_src]), np.array([sl_src]), np.array([se_src]),
                                              r["nlp_max"], r["lse_min"], r["lse_max"], dof)
            e_se, e_sl = abs(r["slope_se"] - se_src), abs(r["slope"] - sl_src)
            worst["-log10 p"] = max(worst["-log10 p"], float(e_nlp))
            worst["slope_se"] = max(worst["slope_se"], float(e_se / se_b[0]) if se_b[0] > 0 else 0.0)
            worst["slope"] = max(worst["slope"], float(e_sl / sl_b[0]) if sl_b[0] > 0 else 0.0)
            if e_nlp > lim_nlp or e_se > packfmt.BOUND_FACTOR * se_b[0] or e_sl > packfmt.BOUND_FACTOR * sl_b[0]:
                bad.append(pid)
        check(not bad, f"cis scan {label} {what}: every phenotype's p, slope_se and slope are within SPEC section 9 "
                       f"(worst -log10 p error {worst['-log10 p']:.2e}, slope_se {worst['slope_se']:.4f} and slope "
                       f"{worst['slope']:.4f} of the per-row bound, limit {packfmt.BOUND_FACTOR}) ({bad[:3]})")
    return ref


def validate(cfg: Config, con, check) -> None:
    """SPEC sections 13 to 15: the hits packs against their sources, the rsID index's block math, the
    startup file against the files it indexes, and the reference files for `npm run pack-check`."""
    d, pk = cfg.derived, cfg["packs"]
    F, B = _frame_variants(cfg), _block_records(cfg)
    t_start = time.time()
    con.execute(f"SET memory_limit = '{pk['duckdb_memory_limit']}'")
    con.execute(f"SET threads = {int(pk['duckdb_threads'])}")
    out_dir = cfg.tmp / "pack_check" / "variant"
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(parents=True)

    # 1. ord: the row position in search_index's own (chr, tss, gene_id) order
    si = read_search_index(cfg, ["gene_id", "symbol", "chr", "tss", "gene_version", "ord"])
    if "ord" not in si.column_names:
        check(False, "search_index has an ord column (run `build --step search_index --force`)")
        return
    ordv = si["ord"].to_numpy()
    rows = si.to_pylist()
    in_order = rows == sorted(rows, key=lambda r: (r["chr"], r["tss"], r["gene_id"]))
    check(np.array_equal(ordv, np.arange(si.num_rows)) and si.num_rows <= 65536 and in_order,
          f"search_index ord is the row position for all {si.num_rows:,} genes, below 65,536, and the rows are in "
          f"(chr, tss, gene_id) order (max ord {int(ordv.max())}, sorted {in_order})")
    gene_of = {i: (r["gene_id"], r["gene_version"]) for i, r in enumerate(rows)}

    vip = variant_idx_path(cfg)
    if not vip.exists():
        check(False, f"{vip.name} exists (run `build --step pack_hits`)")
        return
    con.execute(f"CREATE OR REPLACE TABLE si AS SELECT gene_id, ord FROM {register_search_index(cfg, con, 'si_val')}")

    # 2. the hits packs: structure, totals, and a sample of rows against the source
    rng = np.random.default_rng(SEED)
    totals = {k: 0 for k in KIND_COUNT_KEY}
    n_rows = n_frames = n_bytes = 0
    struct_bad, row_bad, pid_bad, worst = [], [], [], {}
    sampled = 0
    ref = {"frames": [], "rows": []}
    for chrom in CHROMS:
        path = hits_path(cfg, chrom)
        try:
            h = pt.read_header(path)
            buf = path.read_bytes()
            frames = pt.walk_hits_frames(buf)
            total = pt.read_header(pack_paths(cfg, chrom)[0])["count"]
            want_frames = -(-total // F)
            off = packfmt.FILE_HEADER_LEN
            contiguous = True
            for o, ln in frames:
                contiguous &= o == off
                off += ln
            if not (h["kind"] == packfmt.KIND_HITS and h["chrom"] == chrom and h["page_size"] == F
                    and h["count"] == want_frames == len(frames) and contiguous and off == len(buf)):
                struct_bad.append(f"{chrom} ({h['count']} frames for {want_frames} expected, ends at {off:,} of {len(buf):,} B)")
            n_frames += len(frames)
            n_bytes += len(buf)
            c = pt.check_file(path)
            n_rows += c["rows"]
            for i, k in enumerate(KIND_COUNT_KEY):
                totals[k] += c["rows_by_kind"][packfmt.HITS_KIND_NAMES[i]]
            sample = np.sort(rng.choice(total, size=min(SAMPLE_PER_CHROM, total), replace=False))
            tb = _hit_rows(cfg, con, chrom, sample=sample)
            S = _source_arrays(tb)
            D, sc = _decoded_arrays(path, sample, F)
            sampled += len(sample)
            row_bad += _compare_hits(chrom, S, D, sc, worst)
            # sQTL phenotype ids rebuilt from the intron fields, on every sampled odd-kind row
            odd = (D["kind"] % 2) == 1
            if odd.any():
                built = [f"{chrom}:{int(a)}:{int(b)}:clu_{int(cl)}_{st}:{gene_of[int(g)][0]}.{gene_of[int(g)][1]}"
                         for a, b, cl, st, g in zip(D["intron_start"][odd], D["intron_end"][odd], D["cluster"][odd],
                                                    D["strand"][odd], D["gene"][odd])]
                src_pid = sorted(tb.filter(pa.compute.equal(pa.compute.bit_wise_and(tb["kind"], 1), 1))["phenotype_id"].to_pylist()) \
                    if "phenotype_id" in tb.column_names else None
                if src_pid is not None and sorted(built) != src_pid:
                    pid_bad.append(chrom)
        except (PackError, ValueError, OSError) as e:
            struct_bad.append(f"{chrom}: {e}")
    check(not struct_bad, f"every hits pack is kind 7 with ceil(variants / {F}) frames back to back from byte 32 to the "
                          f"file end ({struct_bad[:3]})")
    check(n_rows == EXPECT["rows"], f"hits packs hold {n_rows:,} rows in {n_frames:,} frames, {n_bytes:,} B "
                                    f"(expected {EXPECT['rows']:,} rows)")
    kind_ok = all(totals[k] == EXPECT[k] for k in KIND_COUNT_KEY)
    check(kind_ok, "hits rows by kind: " + ", ".join(f"{k} {totals[k]:,} (expected {EXPECT[k]:,})" for k in KIND_COUNT_KEY))
    check(not row_bad, f"{sampled:,} sampled variant indices decode to exactly their source rows, values within SPEC "
                       f"section 13 (max err " + ", ".join(f"{k} {v:.2e}" for k, v in sorted(worst.items())) + f") ({row_bad[:3]})")
    check(not pid_bad, f"every sampled kind 1, 3 and 5 row rebuilds its sQTL phenotype_id from the intron fields ({pid_bad})")

    # 3. the rsID index
    rp = rsid_index_path(cfg)
    rh = pt.read_header(rp)
    n_rec = rh["count"]
    prev, inc = None, True
    firsts = []
    for b in range(-(-n_rec // B)):
        o, ln = packfmt.rsid_block_range(b, n_rec, B)
        blk = packfmt.decode_rsid_block(pt.read_range(rp, o, ln), what=f"rsid block {b}")
        firsts.append(int(blk["rs_number"][0]))
        if prev is not None and int(blk["rs_number"][0]) <= prev:
            inc = False
        prev = int(blk["rs_number"][-1])
    check(rh["kind"] == packfmt.KIND_RSID and rh["chrom"] == "all" and n_rec == EXPECT["rsid"] and inc
          and rh["bytes"] == packfmt.FILE_HEADER_LEN + 8 * n_rec,
          f"rsid_index.qbr is kind 8 with {n_rec:,} strictly increasing 8-byte records in {len(firsts):,} blocks of {B}, "
          f"{rh['bytes']:,} B (expected {EXPECT['rsid']:,} records)")
    samp = con.execute(f"""SELECT rs_number, chr, position, vidx FROM '{vip}' WHERE rs_number IS NOT NULL
        USING SAMPLE 10000 ROWS (reservoir, {SEED})""").fetch_arrow_table()
    miss = []
    for rs, chrom, pos, vidx in zip(samp["rs_number"].to_pylist(), samp["chr"].to_pylist(),
                                    samp["position"].to_pylist(), samp["vidx"].to_pylist()):
        got = pt.rsid_lookup(rp, int(rs))
        if got is None or got["chr"] != chrom or got["vidx"] != vidx:
            miss.append(f"rs{rs} -> {got}")
    check(not miss, f"{samp.num_rows:,} sampled rs numbers resolve through the block math to their own (chromosome, vidx) "
                    f"({miss[:3]})")

    # 4. the startup file against the files it indexes
    vxp = variant_index_path(cfg)
    idx = pt.read_variant_index(vxp)
    vx_bad = []
    if idx["rsid_first"].tolist() != firsts:
        vx_bad.append("rsid_first differs from the blocks' first records")
    if idx["rsid_n_records"] != n_rec or idx["page_size"] != int(pk["variant_page_size"]) or idx["frame_variants"] != F:
        vx_bad.append("the header constants differ from the files")
    for chrom in CHROMS:
        c = idx["chroms"][chrom]
        vh, pages, first = packfmt.variants_page_firsts(pack_paths(cfg, chrom)[0].read_bytes())
        want_off = [p["offset"] for p in pages] + [pages[-1]["offset"] + pages[-1]["length"]]
        if c["page_off"].tolist() != want_off or c["page_first_position"].tolist() != first.tolist():
            vx_bad.append(f"{chrom}: page offsets or first positions differ from a fresh walk")
        if (c["n_cis"], c["n_trans_only"]) != (vh["n_cis"], vh["count"] - vh["n_cis"]):
            vx_bad.append(f"{chrom}: n_cis or n_trans_only differs from the variants file header")
        hb = hits_path(cfg, chrom).read_bytes()
        hf = pt.walk_hits_frames(hb)
        if c["hits_off"].tolist() != [o for o, _ in hf] + [len(hb)]:
            vx_bad.append(f"{chrom}: hits frame offsets differ from the hits pack")
    check(not vx_bad, f"variant_index.qbx holds every variants-file page offset and first position, every hits frame "
                      f"offset, and each rsID block's first number ({vx_bad[:3]})")

    # the paper's variants resolve both ways, by rsID and by chr:pos, as the browser does
    pv_bad = []
    for key, rsid in cfg["paper_variants"].items():
        chrom, pos = key.split(":")[0], int(key.split(":")[1])
        by_rs = pt.rsid_lookup(rp, int(rsid[2:]), variant_index=vxp)
        page, off, ln = packfmt.variant_index_position(idx, chrom, pos, "cis")
        recs = pt.variant_rows(pack_paths(cfg, chrom)[0], off=off, length=ln) if False else \
            pt.variant_rows(pack_paths(cfg, chrom)[0], vidx=by_rs["vidx"], n=1) if by_rs else None
        ok = by_rs is not None and by_rs["chr"] == chrom and recs is not None \
            and recs["position"][0].as_py() == pos and recs["rs_number"][0].as_py() == int(rsid[2:])
        if not ok:
            pv_bad.append(f"{key} {rsid} -> {by_rs}")
        else:
            ref["frames"].append({"variant": rsid, "chrom": chrom, "position": pos, "vidx": by_rs["vidx"],
                                  "rsid_block": by_rs["block"], "rsid_byte_start": by_rs["byte_start"],
                                  "rsid_byte_end": by_rs["byte_end"]})
    check(not pv_bad, f"the {len(cfg['paper_variants'])} paper variants resolve by rsID through the index and land on "
                      f"their own chr:pos ({pv_bad})")

    # 4.5 the cis scan of SPEC section 7 against the raw Zenodo nominal files
    scan_ref = []
    xlead = con.execute(f"""SELECT lead_position, lead_rsid, symbol FROM '{cfg.tables / 'genes.parquet'}'
        WHERE chr = 'chrX' AND is_egene ORDER BY pval_perm, gene_id LIMIT 1""").fetchone()
    for chrom, position, label in (("chr10", 73_661_450, "rs10824026"),
                                   ("chrX", int(xlead[0]), f"{xlead[1]} ({xlead[2]} lead)")):
        scan_ref.append(_scan_parity(cfg, con, check, chrom, position, label))

    # 5. reference files for `npm run pack-check`
    picks = []
    for chrom, why in (("chr10", "rs10824026's frame"), ("chrX", "a chrX frame"), ("chr21", "the short last frame")):
        h = pt.read_header(hits_path(cfg, chrom))
        c = idx["chroms"][chrom]
        g = c["n_frames"] - 1 if why.startswith("the short") else \
            (next((r["vidx"] for r in ref["frames"] if r["chrom"] == chrom), 0) // F)
        o, ln = int(c["hits_off"][g]), int(c["hits_off"][g + 1] - c["hits_off"][g])
        f = pt.read_hits_frame(hits_path(cfg, chrom), o, ln, g * F)
        rid = len(picks)
        picks.append({"id": rid, "chrom": chrom, "why": why, "file": str(hits_path(cfg, chrom).relative_to(d)),
                      "frame": g, "byte_start": o, "byte_end": o + ln, "first_vidx": g * F,
                      "n_variants": f["n_variants"], "rows": f["n_rows"],
                      "scales": {k: f[k] for k in ("trans_nlp_max", "trans_beta_max", "perm_nlp_max", "se_max", "slope_max")}})
        t = pt.hits_table(f, chrom, pt._gene_rows(search_index_path(cfg)))
        ref["rows"].append(t.append_column("ref", pa.array(np.full(t.num_rows, rid, dtype=np.int32))))
    (out_dir / "hits.json").write_text(json.dumps({
        "note": "Each frame is bytes byte_start..byte_end - 1 of its hits pack, one zstd frame. hits_rows.arrow holds "
                "the rows it must decode to, tagged by ref id, in file order; values carry the frame's own scales.",
        "frame_variants": F, "tolerances": TOL, "frames": picks, "variants": ref["frames"]}, indent=1))
    _write_arrow(pa.concat_tables(ref["rows"]), out_dir / "hits_rows.arrow")
    blocks = []
    for b in (0, -(-n_rec // B) - 1):
        o, ln = packfmt.rsid_block_range(b, n_rec, B)
        blocks.append({"block": b, "byte_start": o, "byte_end": o + ln, "records": ln // 8})
    (out_dir / "rsid.json").write_text(json.dumps({
        "note": "Uncompressed 8-byte records (u32 rs_number, u32 (chr_ordinal << 27) | vidx). rsid_rows.arrow holds the "
                "records each block must decode to, tagged by block.",
        "file": str(rp.relative_to(d)), "records": n_rec, "block_records": B, "blocks": blocks}, indent=1))
    _write_arrow(pa.concat_tables([pt.rsid_records(rp, block=b["block"]) for b in blocks]), out_dir / "rsid_rows.arrow")
    (out_dir / "variant_index.json").write_text(json.dumps({
        "note": "The whole startup file: one zstd frame after a 32-byte kind 9 header. Per chromosome the decoded "
                "counts, and the first and last page and frame offsets, so a decoder can be checked without the arrays.",
        "file": str(vxp.relative_to(d)), "bytes": vxp.stat().st_size, "page_size": idx["page_size"],
        "frame_variants": idx["frame_variants"], "rsid_block_records": idx["rsid_block_records"],
        "rsid_n_records": idx["rsid_n_records"], "rsid_n_blocks": idx["rsid_n_blocks"],
        "chroms": {c: {"n_cis": v["n_cis"], "n_trans_only": v["n_trans_only"], "n_pages_cis": v["n_pages_cis"],
                       "n_pages_trans": v["n_pages_trans"], "n_frames": v["n_frames"],
                       "page_off_first": int(v["page_off"][0]), "page_off_last": int(v["page_off"][-1]),
                       "hits_off_first": int(v["hits_off"][0]), "hits_off_last": int(v["hits_off"][-1])}
                   for c, v in idx["chroms"].items()}}, indent=1))
    _write_arrow(pt.variant_index_table(vxp), out_dir / "variant_index_chroms.arrow")
    scan_rows = [s.pop("rows") for s in scan_ref]
    (out_dir / "scan.json").write_text(json.dumps({
        "note": "SPEC section 7's cis scan. Each scan is two ranges, byte_start..byte_end - 1 of the eQTL and sQTL "
                "packs, walked by 64-byte block headers. scan_rows.arrow holds the rows each must decode to, tagged "
                "by scan id; p, slope_se and slope carry their block's own scales, so compare with the bounds below.",
        "bound_factor": packfmt.BOUND_FACTOR, "nlp_half_step_divisor": TOL["nlp_half_step_divisor"],
        "scans": scan_ref}, indent=1))
    _write_arrow(pa.concat_tables(scan_rows), out_dir / "scan_rows.arrow")
    log(f"pack_variant validate: reference files in {out_dir.relative_to(d)}/ ({time.time() - t_start:.0f} s)")
