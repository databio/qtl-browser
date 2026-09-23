"""Steps pack_sqtl and pack_eqtl: the binary packs of SPEC.md v0, and their checks.

`pack_sqtl` runs first. For every chromosome it streams the raw Zenodo sQTL nominal file (every tested
intron, not only the significant ones) into the published file for `sqtl/<chr>`, one
block per intron in the raw file's order, and writes `_tmp/pack_pointers/sqtl_<chr>.{parquet,json}`:
each intron's block (`blk_off`, `blk_len`) and variant run (`var_start`, `n_var`). One process per
chromosome (`packs.sqtl_workers`), each with its own small DuckDB.

`pack_eqtl` then writes, for every chromosome:

- `variants/<chr>`: the chromosome's cis variants (`in_cis` rows of `_tables/variants.parquet`,
  ordered position, A1, A2), `packs.variant_page_size` per page. af and counts come from the eQTL
  nominal rows, else from the raw sQTL nominal rows (the derived sQTL table keeps significant
  introns only, which would leave variants tested only in other introns without values);
- `eqtl/<chr>`: one block per binned gene, in TSS order, with details built from the source tables.
  Each `splice` entry carries its intron's `blk_off`/`blk_len` from the sQTL pointers, and the union
  variant range takes intron runs from the same files.

The separate `search_index` step joins per-chromosome pointer files after all packs are built.
Every one of these files is content-addressed: it publishes to `data/derived/immutable/` as
`<kind>.<chr>.<sha16>.<ext>` (SPEC section 3), so one byte changing changes its name and the
manifest is the only place the names live. The encoders are `packfmt`'s (the reference codec).

`validate` reads the files with its own decoder below, written from SPEC rather than from
`packfmt`, so an encoder bug cannot cancel out in the check. Only SPEC section 9's pass limit
(`packfmt.error_bounds`, `packfmt.BOUND_FACTOR`) is shared. It also writes reference files for the
browser decoder check (`cd ui && npm run pack-check`) to `data/derived/_tmp/pack_check/`. It also reads the GWAS
packs and index of `pack_gwas` (step in `steps_gwas`, SPEC section 11) against its own read of the source TSV.
"""
from __future__ import annotations

import json
import math
import mmap
import resource
import shutil
import struct
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import zstandard
from scipy.special import stdtrit

from . import packfmt_v0 as packfmt
from .common import (CHROMS, PACKS_METADATA_KEY, SEARCH_INDEX_EXT, Config, addressed_files, connect, log,
                     pack_file, phenotype_batches, publish_file, read_search_index, register_search_index,
                     search_index_path, stage, variants_path, variants_sql, write_parquet)

PACK_COLUMNS = ["blk_off", "blk_len", "var_start", "n_var", "var_off", "var_len", "w_lo", "w_hi"]
PACK_SCHEMA = pa.schema([("gene_id", pa.string())] + [(c, pa.uint32()) for c in PACK_COLUMNS])
BANDS = [0.0, 0.01, 0.1, 1.0, math.inf]      # |t| bands for reporting


# The extension each pack kind publishes under. Every one of these files is content-addressed:
# it lives flat under `immutable/` as `<kind>.<chr>.<sha16>.<ext>` (SPEC section 3).
EXT = {"variants": "qbv", "eqtl": "qbe", "sqtl": "qbs", "trans": "qbt", "gwas": "qbg",
       "gwas_index": "bin", "hits": "qbh", "rsid_index": "qbr", "variant_index": "qbx",
       "search_index": SEARCH_INDEX_EXT}
# the kinds whose bytes `search_index` points into, directly (variants, eqtl, trans) or through the
# gene-details frames of the eQTL pack (sqtl). `search_index` records their SHA-256s (SPEC section 3).
INDEX_PACK_KINDS = ("variants", "eqtl", "sqtl", "trans")


def pack_paths(cfg: Config, chrom: str) -> tuple[Path, Path]:
    """(variants file, eQTL pack) for one chromosome, as published under `immutable/`."""
    return pack_file(cfg, f"variants/{chrom}", EXT["variants"]), pack_file(cfg, f"eqtl/{chrom}", EXT["eqtl"])


def sqtl_path(cfg: Config, chrom: str) -> Path:
    return pack_file(cfg, f"sqtl/{chrom}", EXT["sqtl"])


def gwas_path(cfg: Config, chrom: str) -> Path:
    return pack_file(cfg, f"gwas/{chrom}", EXT["gwas"])


def gwas_index_path(cfg: Config) -> Path:
    return pack_file(cfg, "gwas_index", EXT["gwas_index"])


def pointer_dir(cfg: Config) -> Path:
    """Per-chromosome pointer and stats files that pack steps hand to later steps (never uploaded)."""
    return cfg.tmp / "pack_pointers"


def _nominal(cfg: Config, chrom: str, bin_: int | str = "*") -> str:
    return f"read_parquet('{cfg.tables}/cis_eqtl_nominal/chr={chrom}/bin={bin_}/data.parquet', hive_partitioning=false)"


def _raw_sqtl(cfg: Config, chrom: str) -> Path:
    return cfg.raw_dir("cis_sQTL_nominal") / f"topchefSplice_{chrom}_MaxPC25.cis_qtl_pairs.{chrom}.parquet"


def _require(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError(msg)


def _tree_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


# ---- build ----------------------------------------------------------------------------------------
def _variants(cfg: Config, con, chrom: str) -> pa.Table:
    """Table `v` (vidx, position, A1, A2) in DuckDB, and the per-variant values in vidx order, after
    asserting that af, ma_samples, and ma_count are one value per variant within and across types."""
    d = cfg.derived
    con.execute(f"""CREATE OR REPLACE TABLE v AS
        SELECT (row_number() OVER (ORDER BY position, A1, A2) - 1)::INTEGER AS vidx, position, A1, A2, rs_number, match
        FROM {variants_sql(cfg, chrom)} WHERE in_cis""")
    n, distinct, nulls = con.execute("""SELECT count(*), (SELECT count(*) FROM (SELECT DISTINCT position, A1, A2 FROM v)),
        count(*) FILTER (WHERE position IS NULL OR A1 IS NULL OR A2 IS NULL) FROM v""").fetchone()
    _require(distinct == n and nulls == 0, f"{chrom}: cis variant keys are not unique ({n - distinct} repeats) or null ({nulls})")
    unmatched = con.execute(f"SELECT count(*) FROM {_nominal(cfg, chrom)} n ANTI JOIN v USING (position, A1, A2)").fetchone()[0]
    _require(unmatched == 0, f"{chrom}: {unmatched} eQTL nominal rows do not join to the cis variant list")
    raw = _raw_sqtl(cfg, chrom)
    _require(raw.exists(), f"{chrom}: raw sQTL nominal file missing ({raw}); run `build --step extract`")
    unmatched = con.execute(f"SELECT count(*) FROM '{raw}' r ANTI JOIN v USING (position, A1, A2)").fetchone()[0]
    _require(unmatched == 0, f"{chrom}: {unmatched} raw sQTL nominal rows do not join to the cis variant list")

    def one_value(c: str) -> str:
        return f"count(DISTINCT {c}) > 1 OR count({c}) NOT IN (0, count(*))"

    con.execute(f"""CREATE OR REPLACE TABLE ve AS
        SELECT v.vidx, min(n.af) AS af, min(n.ma_samples) AS ms, min(n.ma_count) AS mc,
               {one_value('n.af')} OR {one_value('n.ma_samples')} OR {one_value('n.ma_count')} AS differs,
               bool_or(n.rs_number IS DISTINCT FROM v.rs_number) AS rs_differs
        FROM {_nominal(cfg, chrom)} n JOIN v USING (position, A1, A2) GROUP BY v.vidx""")
    con.execute(f"""CREATE OR REPLACE TABLE vs AS
        SELECT v.vidx, min(r.af::FLOAT) AS af, min(r.ma_samples::SMALLINT) AS ms, min(r.ma_count::SMALLINT) AS mc,
               {one_value('r.af::FLOAT')} OR {one_value('r.ma_samples')} OR {one_value('r.ma_count')} AS differs
        FROM '{raw}' r JOIN v USING (position, A1, A2) GROUP BY v.vidx""")
    de, drs = con.execute("SELECT count(*) FILTER (WHERE differs), count(*) FILTER (WHERE rs_differs) FROM ve").fetchone()
    ds = con.execute("SELECT count(*) FILTER (WHERE differs) FROM vs").fetchone()[0]
    dx = con.execute("""SELECT count(*) FROM ve JOIN vs USING (vidx)
        WHERE ve.af IS DISTINCT FROM vs.af OR ve.ms IS DISTINCT FROM vs.ms OR ve.mc IS DISTINCT FROM vs.mc""").fetchone()[0]
    _require(de == 0 and ds == 0 and dx == 0 and drs == 0,
             f"{chrom}: per-variant values differ: {de} within eQTL, {ds} within sQTL, {dx} between types, {drs} eQTL rs_number")
    t = con.execute("""SELECT v.vidx, v.position, v.A1, v.A2, v.rs_number, v.match,
            CASE WHEN ve.vidx IS NOT NULL THEN ve.af ELSE vs.af END AS af,
            CASE WHEN ve.vidx IS NOT NULL THEN ve.ms ELSE vs.ms END AS ma_samples,
            CASE WHEN ve.vidx IS NOT NULL THEN ve.mc ELSE vs.mc END AS ma_count,
            ve.vidx IS NOT NULL AS in_e, vs.vidx IS NOT NULL AS in_s
        FROM v LEFT JOIN ve USING (vidx) LEFT JOIN vs USING (vidx) ORDER BY vidx""").fetch_arrow_table()
    _require(np.array_equal(t["vidx"].to_numpy(), np.arange(n)), f"{chrom}: vidx is not 0..n-1")
    in_e = t["in_e"].to_numpy(zero_copy_only=False)
    in_s = t["in_s"].to_numpy(zero_copy_only=False)
    log(f"pack_eqtl {chrom}: {n:,} cis variants: {int(in_e.sum()):,} eQTL-tested, {int((~in_e & in_s).sum()):,} sQTL only, "
        f"{int((~in_e & ~in_s).sum()):,} tested by neither (af and counts stored as not known)")
    return t


def _pairs_for_bin(cfg: Config, con, chrom: str, b: int) -> dict[str, dict[str, np.ndarray]]:
    """One bin file's eQTL rows joined to vidx, grouped by gene (rows in vidx order)."""
    f = cfg.tables / "cis_eqtl_nominal" / f"chr={chrom}" / f"bin={b}" / "data.parquet"
    if not f.exists():
        return {}
    t = con.execute(f"""SELECT n.gene_id, v.vidx, n.position, n.tss_distance, n.pval_nominal, n.slope, n.slope_se, n.pip, n.cs_id
        FROM {_nominal(cfg, chrom, b)} n JOIN v USING (position, A1, A2) ORDER BY n.gene_id, v.vidx""").fetch_arrow_table()
    if t.num_rows == 0:
        return {}
    gid = t["gene_id"].to_numpy(zero_copy_only=False)
    cols = {k: t[k].to_numpy(zero_copy_only=False)
            for k in ("vidx", "position", "tss_distance", "pval_nominal", "slope", "slope_se", "pip", "cs_id")}
    starts = np.r_[0, np.flatnonzero(gid[1:] != gid[:-1]) + 1]
    ends = np.r_[starts[1:], len(gid)]
    out = {str(gid[s]): {k: a[s:e] for k, a in cols.items()} for s, e in zip(starts, ends)}
    _require(len(out) == len(starts), f"{chrom} bin {b}: a gene's rows are not contiguous")
    return out


def _details(cfg: Config, con, chrom: str) -> tuple[list[tuple], dict[str, dict]]:
    """Genes in block order and their compact details, built without gene_detail parquet."""
    d = cfg.derived
    genes_table = con.execute(f"SELECT * FROM '{cfg.tables / 'genes.parquet'}' WHERE chr = ? AND bin IS NOT NULL ORDER BY tss, gene_id", [chrom]).fetch_arrow_table()
    gene_rows = {r["gene_id"]: r for r in genes_table.to_pylist()}
    con.execute(f"""CREATE OR REPLACE TABLE ex AS
        WITH e AS (SELECT gene_id, start, \"end\" FROM '{cfg.tables / 'exons.parquet'}' WHERE gene_id IN (SELECT gene_id FROM '{cfg.tables / 'genes.parquet'}' WHERE chr = ? AND bin IS NOT NULL)),
        o AS (SELECT gene_id, start, \"end\", max(\"end\") OVER (PARTITION BY gene_id ORDER BY start, \"end\" ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_max FROM e),
        grp AS (SELECT *, sum(CASE WHEN prev_max IS NULL OR start > prev_max THEN 1 ELSE 0 END) OVER (PARTITION BY gene_id ORDER BY start, \"end\") AS k FROM o)
        SELECT gene_id, min(start)::INTEGER AS start, max(\"end\")::INTEGER AS \"end\" FROM grp GROUP BY gene_id, k ORDER BY gene_id, start""", [chrom])
    exons: dict[str, list[tuple[int, int]]] = {}
    for gid, start, end in con.execute("SELECT * FROM ex").fetchall():
        exons.setdefault(gid, []).append((start, end))
    splice: dict[str, list[dict]] = {}
    ptr = pointer_dir(cfg) / f"sqtl_{chrom}.parquet"
    _require(ptr.exists(), f"{chrom}: sQTL pointers missing ({ptr}); run `build --step pack_sqtl` first")
    rows = con.execute(f"""SELECT s.gene_id AS _gene, s.* EXCLUDE (gene_id, symbol, chr, tss), p.blk_off, p.blk_len
        FROM '{cfg.tables / 'splice_phenotypes.parquet'}' s LEFT JOIN '{ptr}' p USING (phenotype_id)
        WHERE s.chr = ? ORDER BY s.gene_id, s.cluster_id, s.intron_start, s.intron_end""", [chrom]).fetch_arrow_table()
    for row in rows.to_pylist():
        _require(row["blk_off"] is not None and row["blk_len"] is not None, f"{chrom}: intron {row['phenotype_id']} has no block in {ptr.name}")
        splice.setdefault(row.pop("_gene"), []).append(row)
    details = {gid: packfmt.details_from_tables(row, exons.get(gid, []), splice.get(gid, [])) for gid, row in gene_rows.items()}
    genes = [(r["gene_id"], r["tss"], r["bin"], r["tested"]) for r in gene_rows.values()]
    return genes, details


def _inputs(cfg: Config, chrom: str) -> list[Path]:
    d = cfg.derived
    return (list((cfg.tables / "cis_eqtl_nominal" / f"chr={chrom}").glob("bin=*/data.parquet")) +
            [variants_path(cfg), _raw_sqtl(cfg, chrom),
             cfg.tables / "genes.parquet", cfg.tables / "exons.parquet", cfg.tables / "splice_phenotypes.parquet", cfg.tables / "credible_sets.parquet",
             pointer_dir(cfg) / f"sqtl_{chrom}.parquet"])


def eqtl(cfg: Config, force: bool = False) -> None:
    d = cfg.derived
    pk = cfg["packs"]
    level, page_size, codec = pk["zstd_level"], pk["variant_page_size"], pk["variant_page_codec"]
    con = connect(cfg, memory_limit=pk["duckdb_memory_limit"], threads=pk["duckdb_threads"], temp_dir=cfg.tmp / "pack_eqtl")
    zc = zstandard.ZstdCompressor(level=level)
    totals = {"genes": 0, "rows": 0, "memberships": 0, "variants": 0, "eqtl_bytes": 0, "variant_bytes": 0}
    for chrom in CHROMS:
        var_key, eqtl_key = f"variants/{chrom}", f"eqtl/{chrom}"
        var_path, eqtl_path = pack_paths(cfg, chrom)

        inputs = _inputs(cfg, chrom)
        _require(all(p.exists() for p in inputs), f"{chrom}: missing pack input(s): {[str(p) for p in inputs if not p.exists()]}")
        if not force and var_path.exists() and eqtl_path.exists() and min(var_path.stat().st_mtime, eqtl_path.stat().st_mtime) > max(p.stat().st_mtime for p in inputs):
            log(f"pack_eqtl {chrom}: fresh, skipping")
            continue
        vt = _variants(cfg, con, chrom)
        vbytes, page_offsets = packfmt.encode_variants_file(
            chrom, vt["position"].to_numpy(), vt["rs_number"].to_numpy(zero_copy_only=False), vt["af"].to_numpy(zero_copy_only=False),
            vt["ma_samples"].to_numpy(zero_copy_only=False), vt["ma_count"].to_numpy(zero_copy_only=False),
            vt["A1"], vt["A2"], vt["match"], page_size, codec, level)
        var_tmp = stage(cfg, var_key, EXT["variants"])
        var_tmp.write_bytes(vbytes)

        genes, details = _details(cfg, con, chrom)
        # intron runs from pack_sqtl's pointers; pack_sqtl checked every run against the raw rows
        sp = pq.read_table(pointer_dir(cfg) / f"sqtl_{chrom}.parquet", columns=["gene_id", "var_start", "n_var"])
        intron_ranges: dict[str, tuple[int, int]] = {}
        for gid, i_start, i_n in zip(sp["gene_id"].to_pylist(), sp["var_start"].to_pylist(), sp["n_var"].to_pylist()):
            lo, hi = intron_ranges.get(gid, (i_start, i_start + i_n))
            intron_ranges[gid] = (min(lo, i_start), max(hi, i_start + i_n))
        _require(all(hi <= len(vt) for _, hi in intron_ranges.values()),
                 f"{chrom}: an intron run in the sQTL pointers ends past the {len(vt):,} cis variants")
        cs_rows: dict[str, list[tuple[int, float, int]]] = {}
        cst = con.execute(f"""SELECT c.phenotype_id, v.vidx, c.pip, c.cs_id FROM '{cfg.tables / 'credible_sets.parquet'}' c
            JOIN v USING (position, A1, A2) WHERE c.qtl_type = 'e' AND c.chr = ? ORDER BY c.phenotype_id, v.vidx, c.cs_id""", [chrom]).fetchall()
        for gid, vidx, pip, cid in cst:
            cs_rows.setdefault(gid, []).append((int(vidx), float(pip), int(cid)))

        # eQTL pack, one bin file of pairs at a time
        eqtl_tmp = stage(cfg, eqtl_key, EXT["eqtl"])
        n_pairs = pair_z = n_cs = n_tested = wider = more_pages = max_extra = 0
        index_rows: list[dict] = []
        with open(eqtl_tmp, "wb") as fh:
            fh.write(packfmt.file_header(packfmt.KIND_EQTL, chrom, len(genes), 0))
            cur_bin, pairs = None, {}
            for gene_id, _tss, b, tested in genes:
                if b != cur_bin:
                    _require(cur_bin is None or b > cur_bin, f"{chrom}: bins are not in TSS order at {gene_id}")
                    _require(not pairs, f"{chrom} bin {cur_bin}: nominal rows for genes not in search_index: {sorted(pairs)[:5]}")
                    cur_bin, pairs = b, _pairs_for_bin(cfg, con, chrom, b)
                p = pairs.pop(gene_id, None)
                if p is None:
                    _require(not tested, f"{gene_id}: tested for eQTL but has no nominal rows in {chrom} bin {b}")
                    blk = packfmt.encode_gene_block(details[gene_id], None, None, [], [], [], None, None, None, level)
                    var = (None, None, None, None)
                else:
                    _require(bool(tested), f"{gene_id}: has nominal rows but search_index says not tested")
                    vidx = p["vidx"].astype(np.int64)
                    n, var_start = len(vidx), int(vidx[0])
                    if not np.array_equal(vidx, np.arange(var_start, var_start + n)):
                        raise RuntimeError(f"{gene_id}: tested variants are not one contiguous run")
                    pos = p["position"].astype(np.int64)
                    off = pos - p["tss_distance"].astype(np.int64)
                    if np.any(off != off[0]):
                        raise RuntimeError(f"{gene_id}: position - tss_distance is not constant")
                    memberships = cs_rows.pop(gene_id, [])
                    csr = [v - var_start for v, _, _ in memberships]
                    _require(all(0 <= r < n for r in csr), f"{gene_id}: credible-set membership outside nominal run")
                    blk = packfmt.encode_gene_block(
                        details[gene_id], var_start, int(off[0]), p["pval_nominal"], p["slope"], p["slope_se"],
                        csr, [x[1] for x in memberships], [x[2] for x in memberships], level,
                        pos_first=int(pos[0]), pos_last=int(pos[-1]))
                    u_start, u_end = var_start, var_start + n
                    if gene_id in intron_ranges:
                        u_start = min(u_start, intron_ranges[gene_id][0]); u_end = max(u_end, intron_ranges[gene_id][1])
                    if (u_start, u_end) != (var_start, var_start + n): wider += 1
                    old_off, old_len = packfmt.variant_range(page_offsets, page_size, var_start, n)
                    var_off, var_len = packfmt.variant_range(page_offsets, page_size, u_start, u_end - u_start)
                    if (var_off, var_len) != (old_off, old_len): more_pages += 1
                    max_extra = max(max_extra, var_len - old_len)
                    var = (var_start, n, var_off, var_len, int(vt["position"][u_start].as_py()), int(vt["position"][u_end - 1].as_py()))
                    pair_z += len(zc.compress(blk[packfmt.BLOCK_HEADER_LEN:packfmt.BLOCK_HEADER_LEN + 4 * n]))
                    n_pairs += n
                    n_cs += len(memberships)
                    n_tested += 1
                blk_off = fh.tell()
                fh.write(blk)
                if p is None:
                    _require(gene_id in intron_ranges, f"{gene_id}: has a bin but no eQTL rows and no intron runs in {chrom}")
                    u_start, u_end = intron_ranges[gene_id]
                    var_off, var_len = packfmt.variant_range(page_offsets, page_size, u_start, u_end - u_start)
                    var = (None, None, var_off, var_len, int(vt["position"][u_start].as_py()), int(vt["position"][u_end - 1].as_py()))
                index_rows.append(dict(zip(["gene_id"] + PACK_COLUMNS, [gene_id, blk_off, len(blk), *var])))
            _require(not pairs, f"{chrom} bin {cur_bin}: nominal rows for genes not in search_index: {sorted(pairs)[:5]}")
            eqtl_bytes = fh.tell()
        _require(eqtl_bytes <= packfmt.U32_MAX, f"{chrom}: eQTL pack is {eqtl_bytes} bytes, over 4 GiB")
        _require(not cs_rows, f"{chrom}: credible-set genes absent from pack: {list(cs_rows)[:5]}")
        want_cs = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 'e' AND chr = ?", [chrom]).fetchone()[0]
        _require(n_cs == want_cs, f"{chrom}: {n_cs} credible-set records written, credible_sets has {want_cs} eQTL rows")
        _require(len(index_rows) == len(genes), f"{chrom}: {len(index_rows)} blocks for {len(genes)} binned genes")
        _require(n_pairs == pq.read_metadata(cfg.raw_dir("cis_eQTL_nominal") / f"topchef_{chrom}_MaxPC70.cis_qtl_pairs.{chrom}.parquet").num_rows,
                 f"{chrom}: eQTL row total differs from raw source")
        var_path = publish_file(cfg, var_tmp, var_key, EXT["variants"])
        eqtl_path = publish_file(cfg, eqtl_tmp, eqtl_key, EXT["eqtl"])

        today_nom = _tree_bytes(cfg.tables / "cis_eqtl_nominal" / f"chr={chrom}")
        today_gd = _tree_bytes(d / "gene_detail" / f"chr={chrom}")
        packs = eqtl_bytes + len(vbytes)
        log(f"pack_eqtl {chrom}: {len(genes):,} genes ({n_tested:,} eQTL-tested, {len(genes) - n_tested:,} sQTL only), "
            f"{n_pairs:,} pairs, {n_cs:,} credible-set records, {len(vt):,} variants")
        log(f"pack_eqtl {chrom}: eQTL pack {eqtl_bytes:,} B ({eqtl_path.relative_to(d)}), variants file {len(vbytes):,} B "
            f"({var_path.relative_to(d)}), together {packs:,} B")
        log(f"pack_eqtl {chrom}: pair arrays alone {4 * n_pairs:,} B raw (4.00 B/pair), zstd {level} per gene {pair_z:,} B "
            f"({pair_z / max(n_pairs, 1):.2f} B/pair)")
        log(f"pack_eqtl {chrom}: today cis_eqtl_nominal {today_nom:,} B + gene_detail {today_gd:,} B = {today_nom + today_gd:,} B; "
            f"packs are {packs / (today_nom + today_gd):.1%} of that")
        ptr = d / "_tmp" / "pack_pointers"; ptr.mkdir(parents=True, exist_ok=True)
        write_parquet(pa.Table.from_pylist(index_rows, schema=PACK_SCHEMA), ptr / f"eqtl_{chrom}.parquet", 100_000)
        stats = {"genes": len(genes), "eqtl_genes": n_tested, "rows": n_pairs, "memberships": n_cs,
                 "variants": len(vt), "pack_bytes": eqtl_bytes, "variants_bytes": len(vbytes),
                 "union_wider": wider, "union_more_pages": more_pages, "max_extra_bytes": max_extra}
        (ptr / f"eqtl_{chrom}.json").write_text(json.dumps(stats, indent=2))
        log(f"pack_eqtl {chrom}: union wider for {wider} genes, changed pages for {more_pages}; largest extra range {max_extra:,} B")
        con.execute("DROP TABLE IF EXISTS v; DROP TABLE IF EXISTS ve; DROP TABLE IF EXISTS vs; DROP TABLE IF EXISTS ex")
    ptr = d / "_tmp" / "pack_pointers"
    for chrom in CHROMS:
        f = ptr / f"eqtl_{chrom}.json"
        _require(f.exists(), f"{chrom}: pointer stats missing ({f}); re-run with --force")
        st = json.loads(f.read_text())
        for k in ("genes", "rows", "memberships", "variants", "union_wider", "union_more_pages"):
            totals[k] = totals.get(k, 0) + st[k]
        totals["eqtl_bytes"] += st["pack_bytes"]; totals["variant_bytes"] += st["variants_bytes"]
        totals["max_extra_bytes"] = max(totals.get("max_extra_bytes", 0), st["max_extra_bytes"])
    expect = {"genes": 20_628, "rows": 123_458_578, "memberships": 215_462, "variants": 8_872_723}
    for k, want in expect.items():
        log(f"pack_eqtl total {k}: {totals[k]:,} (expected {want:,}{'' if totals[k] == want else ', DIFFERS'})")
    today = _tree_bytes(cfg.tables / "cis_eqtl_nominal") + _tree_bytes(d / "gene_detail")
    packs = totals["eqtl_bytes"] + totals["variant_bytes"]
    log(f"pack_eqtl total bytes: eQTL packs {totals['eqtl_bytes']:,}, variants files {totals['variant_bytes']:,}, together {packs:,}; "
        f"today cis_eqtl_nominal + gene_detail {today:,} ({packs / max(today, 1):.1%})")
    log(f"pack_eqtl total union: wider than the eQTL run for {totals['union_wider']} genes, other pages for "
        f"{totals['union_more_pages']}, largest extra range {totals['max_extra_bytes']:,} B")


def index_pack_shas(cfg: Config) -> dict[str, str]:
    """Logical key -> full SHA-256 for every published file `search_index`'s offsets reach. Stored in
    the index's own Arrow metadata so one manifest fetch pins an index and the packs it matches."""
    from .common import digests
    files = addressed_files(cfg)
    keys = sorted(k for k in files if k.split("/")[0] in INDEX_PACK_KINDS)
    _require(bool(keys), "search_index: no packs are published; run the pack steps first")
    return {k: digests(files[k])[0] for k in keys}


def search_index(cfg: Config) -> None:
    d = cfg.derived
    idx_path = search_index_path(cfg)
    old = read_search_index(cfg) if idx_path.exists() else None
    old_bytes = idx_path.stat().st_size if old is not None else 0
    con = connect(cfg)
    # trans_off/trans_len from pack_trans (null for genes without trans rows); gene_version rebuilds sQTL phenotype ids
    out = con.execute(f"""SELECT g.gene_id, g.symbol, g.chr, g.tss, g.tested, g.is_egene,
        coalesce(s.n, 0)::SMALLINT n_sqtl_sig, g.bin, g.start, g.\"end\", g.strand, g.biotype,
        p.blk_off, p.blk_len, p.var_start, p.n_var, p.var_off, p.var_len, p.w_lo, p.w_hi,
        t.trans_off, t.trans_len, CAST(split_part(g.gene_id_version, '.', 2) AS UTINYINT) AS gene_version
        FROM '{cfg.tables / 'genes.parquet'}' g
        LEFT JOIN (SELECT gene_id, count(*) FILTER (WHERE is_sqtl) n FROM '{cfg.tables / 'splice_phenotypes.parquet'}' GROUP BY 1) s USING (gene_id)
        LEFT JOIN read_parquet('{d}/_tmp/pack_pointers/eqtl_*.parquet') p USING (gene_id)
        LEFT JOIN read_parquet('{d}/_tmp/pack_pointers/trans_*.parquet') t USING (gene_id)
        ORDER BY g.chr, g.tss, g.gene_id""").fetch_arrow_table()
    n_ptr = con.execute(f"SELECT count(*) FROM read_parquet('{d}/_tmp/pack_pointers/trans_*.parquet')").fetchone()[0]
    n_trans = out.num_rows - out["trans_off"].null_count
    _require(n_ptr == n_trans and out["gene_version"].null_count == 0,
             f"search_index: {n_ptr:,} trans pointer rows but {n_trans:,} genes with trans_off; {out['gene_version'].null_count} genes without a version")
    log(f"search_index: {n_trans:,} genes with trans_off/trans_len")
    # ord: the row's position in search_index order (chr, tss, gene_id). Hit rows point at it with a u16
    # instead of carrying a gene id (SPEC section 13), so it must fit one and equal the row position.
    _require(out.num_rows <= 65536, f"search_index: {out.num_rows:,} genes, more than a u16 ord can point at")
    out = out.append_column("ord", pa.array(np.arange(out.num_rows, dtype=np.uint16)))
    # `bin` names a nominal parquet partition, which only the pack builders read now; it stays in
    # _tables/genes.parquet and leaves the browser's index. "Has a gene block" is blk_off.
    if "bin" in out.column_names:
        out = out.drop_columns(["bin"])
    base = out.select([c for c in out.column_names[:12]])
    if old is not None:
        old_base = old.select(base.column_names)
        _require(old_base.schema.equals(base.schema) and old_base.equals(base), "search_index base columns changed")
    # one Arrow IPC stream in one zstd frame (SPEC sections 2 and 6): the browser inserts it
    # straight into DuckDB, so no parquet reader has to exist in the bundle
    #
    # the SHA-256 of every pack this index points into rides in the Arrow schema metadata, so
    # `manifest` can refuse an index built from packs that have since been rebuilt (SPEC section 3)
    shas = index_pack_shas(cfg)
    out = out.replace_schema_metadata({PACKS_METADATA_KEY: json.dumps(shas, sort_keys=True).encode()})
    sink = pa.BufferOutputStream()
    with pa.ipc.new_stream(sink, out.schema) as w:
        w.write_table(out)
    raw = sink.getvalue().to_pybytes()
    frame = packfmt.zstd_frame(raw, int(cfg["packs"]["zstd_level"]))
    tmp = stage(cfg, "search_index", EXT["search_index"])
    tmp.write_bytes(frame)
    idx_path = publish_file(cfg, tmp, "search_index", EXT["search_index"])
    log(f"search_index: {out.num_rows:,} rows, {len(out.column_names)} columns, {len(raw):,} B Arrow -> "
        f"{len(frame):,} B zstd (was {old_bytes:,} B); {idx_path.name}, "
        f"{len(shas)} pack SHA-256s in `{PACKS_METADATA_KEY.decode()}`")


# ---- pack_sqtl ------------------------------------------------------------------------------------
SQTL_POINTER_SCHEMA = pa.schema([("phenotype_id", pa.string()), ("gene_id", pa.string()), ("blk_off", pa.uint32()),
                                 ("blk_len", pa.uint32()), ("var_start", pa.uint32()), ("n_var", pa.uint32()),
                                 ("anchor", pa.int32()), ("n_cs", pa.uint32())])
SQTL_BATCH_ROWS = 1_000_000     # raw rows per streamed table
DOF_CHECK_ROWS = 16             # per intron: this many smallest-p rows, and this many seeded random rows
DOF_CHECK_MIN_T = 0.05
DOF_CHECK_REL = 1e-6


def _sqtl_inputs(cfg: Config, chrom: str) -> list[Path]:
    d = cfg.derived
    return [_raw_sqtl(cfg, chrom), cfg.tables / "credible_sets.parquet", cfg.tables / "splice_phenotypes.parquet",
            variants_path(cfg)]


def _first_bad_row(vpos: np.ndarray, vA1: pa.Array, vA2: pa.Array, vs: int, pos: np.ndarray, a1: pa.Array, a2: pa.Array) -> int:
    for i in range(len(pos)):
        j = vs + i
        if j >= len(vpos) or vpos[j] != pos[i] or vA1[j].as_py() != a1[i].as_py() or vA2[j].as_py() != a2[i].as_py():
            return i
    return len(pos) - 1


def _sqtl_chrom(cfg: Config, chrom: str) -> dict:
    """pack_sqtl worker for one chromosome, in its own process: `<chr>.qbs` and its pointer files."""
    t0 = time.time()
    d = cfg.derived
    pk = cfg["packs"]
    dof, level = int(pk["dof"]["sqtl"]), int(pk["zstd_level"])
    work = cfg.tmp / f"pack_sqtl-{chrom}"
    con = connect(cfg, memory_limit=pk["sqtl_duckdb_memory_limit"], threads=int(pk["sqtl_duckdb_threads"]), temp_dir=work)

    # small inputs: the cis variant list in vidx order, credible sets, intron -> gene
    vt = con.execute(f"""SELECT position, A1, A2 FROM {variants_sql(cfg, chrom)}
        WHERE in_cis ORDER BY position, A1, A2""").fetch_arrow_table()
    vpos = vt["position"].to_numpy().astype(np.int64)
    vA1 = vt["A1"].combine_chunks().cast(pa.string())
    vA2 = vt["A2"].combine_chunks().cast(pa.string())
    nv = len(vpos)
    var_path = pack_paths(cfg, chrom)[0]
    if var_path.exists():      # absent only on a first build, where pack_eqtl writes it next from the same list
        with open(var_path, "rb") as f:
            h = packfmt.parse_file_header(f.read(packfmt.FILE_HEADER_LEN))
        _require(h["kind"] == packfmt.KIND_VARIANTS and h["n_cis"] == nv,
                 f"{chrom}: {nv:,} cis variants, but {var_path.name} holds {h['n_cis']:,}")
    cs: dict[str, list[tuple[int, str, str, float, int]]] = {}
    for ph, cpos, ca1, ca2, pip, cid in con.execute(f"""SELECT phenotype_id, position, A1, A2, pip, cs_id
            FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 's' AND chr = ?""", [chrom]).fetchall():
        cs.setdefault(ph, []).append((int(cpos), ca1, ca2, float(pip), int(cid)))
    want_cs = sum(len(v) for v in cs.values())
    gene_of = dict(con.execute(f"SELECT phenotype_id, gene_id FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE chr = ?", [chrom]).fetchall())
    con.close()

    raw = _raw_sqtl(cfg, chrom)
    pfile = pq.ParquetFile(raw)
    n_raw = pfile.metadata.num_rows
    key = f"sqtl/{chrom}"
    tmp = stage(cfg, key, EXT["sqtl"])
    rng = np.random.default_rng(CHROMS.index(chrom))
    ptr: dict[str, list] = {k: [] for k in SQTL_POINTER_SCHEMA.names}
    seen: set[str] = set()
    n_rows = n_null = n_cs = dof_rows = 0
    dof_err, dof_worst = 0.0, None
    cols = ["phenotype_id", "position", "A1", "A2", "start_distance", "pval_nominal", "slope", "slope_se"]
    with open(tmp, "wb") as fh:
        fh.write(bytes(packfmt.FILE_HEADER_LEN))          # rewritten with the block count at the end
        for tb in phenotype_batches(pfile, cols, SQTL_BATCH_ROWS):
            nb = tb.num_rows
            if nb == 0:
                continue
            ids = tb["phenotype_id"].combine_chunks()
            starts = np.r_[0, pc.indices_nonzero(pc.not_equal(ids.slice(0, nb - 1), ids.slice(1))).to_numpy() + 1].astype(np.int64)
            ends = np.r_[starts[1:], nb]
            bpos = tb["position"].to_numpy().astype(np.int64)
            banchor = bpos - tb["start_distance"].to_numpy().astype(np.int64)
            bp, bsl, bse = (tb[c].to_numpy().astype(np.float64) for c in ("pval_nominal", "slope", "slope_se"))
            bA1 = tb["A1"].combine_chunks().cast(pa.string())
            bA2 = tb["A2"].combine_chunks().cast(pa.string())
            for s, e in zip(starts.tolist(), ends.tolist()):
                ph = ids[s].as_py()
                _require(ph not in seen, f"{chrom}: intron {ph} appears in two places in {raw.name}")
                seen.add(ph)
                gene_id = gene_of.get(ph)
                _require(gene_id is not None, f"{chrom}: intron {ph} is in {raw.name} but not in splice_phenotypes")
                n = e - s
                pos = bpos[s:e]
                a1, a2 = bA1.slice(s, n), bA2.slice(s, n)
                # the run starts at the first variant with this position and these alleles
                vs = int(np.searchsorted(vpos, pos[0]))
                key0 = (a1[0].as_py(), a2[0].as_py())
                while vs < nv and vpos[vs] == pos[0] and (vA1[vs].as_py(), vA2[vs].as_py()) != key0:
                    vs += 1
                ok = vs + n <= nv and np.array_equal(vpos[vs:vs + n], pos)
                if ok:
                    same = pc.and_(pc.equal(vA1.slice(vs, n), a1), pc.equal(vA2.slice(vs, n), a2))
                    ok = same.null_count == 0 and bool(pc.all(same).as_py())
                if not ok:
                    i = _first_bad_row(vpos, vA1, vA2, vs, pos, a1, a2)
                    raise RuntimeError(f"{chrom}: intron {ph} is not one run of the cis variant list from vidx {vs}: "
                                       f"first bad row {i} (position {pos[i]}, {a1[i].as_py()}/{a2[i].as_py()})")
                anchor = banchor[s:e]
                _require(bool(np.all(anchor == anchor[0])), f"{chrom}: intron {ph}: position - start_distance is not constant")
                memberships = []
                for cpos, ca1, ca2, pip, cid in cs.pop(ph, []):
                    lo, hi = int(np.searchsorted(pos, cpos)), int(np.searchsorted(pos, cpos, "right"))
                    row = next((r for r in range(lo, hi) if a1[r].as_py() == ca1 and a2[r].as_py() == ca2), None)
                    _require(row is not None and 0 <= row < n,
                             f"{chrom}: intron {ph}: credible-set variant {cpos} {ca1}/{ca2} is not among its rows")
                    memberships.append((row, cid, pip))
                memberships.sort()
                p, sl, se = bp[s:e], bsl[s:e], bse[s:e]
                blk = packfmt.encode_gene_block(None, vs, int(anchor[0]), p, sl, se, [m[0] for m in memberships],
                                                [m[2] for m in memberships], [m[1] for m in memberships], level,
                                                pos_first=int(pos[0]), pos_last=int(pos[-1]))
                # dof: the slope rebuilt from the unquantized p and SE matches the source
                with np.errstate(divide="ignore", invalid="ignore"):
                    cand = np.flatnonzero((p > 0) & (np.abs(sl) / se > DOF_CHECK_MIN_T))
                if cand.size:
                    pick = cand[np.argsort(p[cand], kind="stable")[:DOF_CHECK_ROWS]]
                    rest = np.setdiff1d(cand, pick)
                    if rest.size:
                        pick = np.r_[pick, rng.choice(rest, min(DOF_CHECK_ROWS, rest.size), replace=False)]
                    ref = np.sign(sl[pick]) * se[pick] * np.maximum(-stdtrit(dof, p[pick] / 2), 0.0)
                    err = float(np.max(np.abs(ref - sl[pick]) / np.abs(sl[pick])))
                    dof_rows += int(pick.size)
                    if err > dof_err:
                        dof_err, dof_worst = err, ph
                blk_off = fh.tell()
                fh.write(blk)
                for k, v in (("phenotype_id", ph), ("gene_id", gene_id), ("blk_off", blk_off), ("blk_len", len(blk)),
                             ("var_start", vs), ("n_var", n), ("anchor", int(anchor[0])), ("n_cs", len(memberships))):
                    ptr[k].append(v)
                n_rows += n
                n_null += int(np.isnan(p).sum())
                n_cs += len(memberships)
        end = fh.tell()
        fh.seek(0)
        fh.write(packfmt.file_header(packfmt.KIND_SQTL, chrom, len(seen), 0))

    missing = sorted(set(gene_of) - seen)
    _require(not missing, f"{chrom}: {len(missing)} splice_phenotypes introns have no raw rows, e.g. {missing[:3]}")
    _require(n_rows == n_raw, f"{chrom}: {n_rows:,} rows written, {raw.name} has {n_raw:,}")
    _require(n_cs == want_cs and not cs, f"{chrom}: {n_cs:,} credible-set records written, credible_sets has {want_cs:,} "
                                         f"sQTL rows ({len(cs)} introns with sets but no raw rows)")
    offs, lens = np.array(ptr["blk_off"], dtype=np.int64), np.array(ptr["blk_len"], dtype=np.int64)
    _require(offs.size > 0 and offs[0] == packfmt.FILE_HEADER_LEN and np.array_equal(offs[1:], offs[:-1] + lens[:-1])
             and end == offs[-1] + lens[-1] == tmp.stat().st_size and end <= packfmt.U32_MAX,
             f"{chrom}: blocks are not contiguous from byte 32 to the file end, or the file is over 4 GiB ({end:,} bytes)")
    _require(dof_err <= DOF_CHECK_REL, f"{chrom}: slope rebuilt from p and SE with dof {dof} misses the source by {dof_err:.3g} (intron {dof_worst})")
    publish_file(cfg, tmp, key, EXT["sqtl"])
    pdir = pointer_dir(cfg)
    pdir.mkdir(parents=True, exist_ok=True)
    write_parquet(pa.Table.from_pydict(ptr, schema=SQTL_POINTER_SCHEMA), pdir / f"sqtl_{chrom}.parquet", 100_000)
    stats = {"chrom": chrom, "introns": len(seen), "rows": n_rows, "null_rows": n_null, "memberships": n_cs, "bytes": end,
             "max_n_var": int(lens.size and max(ptr["n_var"])), "dof": dof, "dof_check_rows": dof_rows,
             "max_dof_rel_error": dof_err, "max_dof_error_intron": dof_worst,
             "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
             "seconds": round(time.time() - t0, 1)}
    (pdir / f"sqtl_{chrom}.json").write_text(json.dumps(stats, indent=2))
    shutil.rmtree(work, ignore_errors=True)
    return stats


def sqtl(cfg: Config, force: bool = False) -> None:
    pk = cfg["packs"]
    pdir = pointer_dir(cfg)
    jobs = []
    for chrom in CHROMS:
        inputs = _sqtl_inputs(cfg, chrom)
        _require(all(p.exists() for p in inputs), f"{chrom}: missing pack_sqtl input(s): {[str(p) for p in inputs if not p.exists()]}")
        outs = [sqtl_path(cfg, chrom), pdir / f"sqtl_{chrom}.parquet", pdir / f"sqtl_{chrom}.json"]
        if not force and all(p.exists() for p in outs) and min(p.stat().st_mtime for p in outs) > max(p.stat().st_mtime for p in inputs):
            log(f"pack_sqtl {chrom}: fresh, skipping")
            continue
        jobs.append(chrom)
    jobs.sort(key=lambda c: -_raw_sqtl(cfg, c).stat().st_size)
    workers = int(pk["sqtl_workers"])
    log(f"pack_sqtl: {len(jobs)} chromosome(s) with {workers} worker(s), DuckDB {pk['sqtl_duckdb_memory_limit']} and "
        f"{pk['sqtl_duckdb_threads']} threads each, biggest raw file first")
    if jobs:
        # a fresh process per chromosome, so its peak RSS is its own and memory is returned between chromosomes
        with ProcessPoolExecutor(max_workers=workers, max_tasks_per_child=1) as ex:
            futs = {ex.submit(_sqtl_chrom, cfg, c): c for c in jobs}
            try:
                for f in as_completed(futs):
                    st = f.result()
                    log(f"pack_sqtl {st['chrom']}: {st['introns']:,} introns, {st['rows']:,} rows ({st['null_rows']:,} null), "
                        f"{st['memberships']:,} credible-set records, {st['bytes']:,} B, largest run {st['max_n_var']:,}; "
                        f"dof check on {st['dof_check_rows']:,} rows, max relative error {st['max_dof_rel_error']:.3g}; "
                        f"peak RSS {st['peak_rss_mb']:,.0f} MB; {st['seconds']:,.0f} s")
            except BaseException:
                ex.shutdown(wait=False, cancel_futures=True)
                raise
    totals = {"introns": 0, "rows": 0, "null_rows": 0, "memberships": 0, "bytes": 0}
    dof_err, rss = 0.0, 0.0
    missing = [c for c in CHROMS if not (pdir / f"sqtl_{c}.json").exists()]
    if missing:
        log(f"pack_sqtl total: stats missing for {missing}; totals not reported")
        return
    for chrom in CHROMS:
        st = json.loads((pdir / f"sqtl_{chrom}.json").read_text())
        for k in totals:
            totals[k] += st[k]
        dof_err, rss = max(dof_err, st["max_dof_rel_error"]), max(rss, st["peak_rss_mb"])
    for k, want in {"introns": 80_750, "rows": 499_599_113, "memberships": 246_322}.items():
        log(f"pack_sqtl total {k}: {totals[k]:,} (expected {want:,}{'' if totals[k] == want else ', DIFFERS'})")
    log(f"pack_sqtl total: {totals['null_rows']:,} null rows; {totals['bytes']:,} B in {len(CHROMS)} files; "
        f"largest dof relative error {dof_err:.3g}; largest worker peak RSS {rss:,.0f} MB")


# ---- validate: an independent reader of SPEC sections 3 to 8 --------------------------------------
_FILE_H = struct.Struct("<4sBBH8sII8s")
_PAGE_H = struct.Struct("<IIHBB")
_BLOCK_H = struct.Struct("<4sIIIiIIIdddII")
_CS = np.dtype([("row", "<u4"), ("pip", "<f4"), ("cs_id", "u1"), ("pad", "u1", (3,))])
_SNP = [None, ("A", "C"), ("A", "G"), ("A", "T"), ("C", "A"), ("C", "G"), ("C", "T"),
        ("G", "A"), ("G", "C"), ("G", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
_ZD = zstandard.ZstdDecompressor()


class PackError(ValueError):
    pass


def _unframe(buf: bytes, expect: int | None, what: str) -> bytes:
    """One zstd frame: content size present (and equal to `expect`), checksum, no dictionary, no trailing bytes."""
    try:
        fp = zstandard.get_frame_parameters(buf)
    except zstandard.ZstdError as e:
        raise PackError(f"{what}: not a zstd frame ({e})") from None
    if fp.content_size in (zstandard.CONTENTSIZE_UNKNOWN, zstandard.CONTENTSIZE_ERROR):
        raise PackError(f"{what}: frame has no content size")
    if not fp.has_checksum or fp.dict_id != 0:
        raise PackError(f"{what}: frame needs a checksum and no dictionary")
    if expect is not None and fp.content_size != expect:
        raise PackError(f"{what}: content size {fp.content_size} != {expect}")
    try:
        out = _ZD.decompress(buf, max_output_size=fp.content_size, allow_extra_data=False)
    except zstandard.ZstdError as e:
        raise PackError(f"{what}: {e}") from None
    if len(out) != fp.content_size:
        raise PackError(f"{what}: decoded {len(out)} bytes, frame says {fp.content_size}")
    return out


def _reject(name):
    raise PackError(f"JSON constant {name}")


def _file_header(buf: bytes, kind: int, chrom: str, what: str) -> tuple[int, int, int]:
    """(count, page size, n_cis). n_cis is the u32 at offset 24 of a variants file (kind 1), at most count;
    in every other kind those 4 bytes are zero, and bytes 28 to 31 are zero in every kind."""
    if len(buf) < 32:
        raise PackError(f"{what}: shorter than the file header")
    magic, k, version, hlen, name, count, page_size, reserved = _FILE_H.unpack_from(buf, 0)
    n_cis, rest = struct.unpack("<II", reserved)
    if (magic, k, version, hlen, name.rstrip(b"\0"), rest) != (b"QTLB", kind, 0, 32, chrom.encode(), 0) \
            or (kind != 1 and n_cis != 0) or n_cis > count:
        raise PackError(f"{what}: file header {magic!r} kind {k} version {version} length {hlen} chrom {name!r} "
                        f"count {count} n_cis/reserved {reserved!r}")
    return count, page_size, n_cis


def _read_block(buf: bytes, dof: int, want_n: int | None, want_var_start: int | None, what: str,
                kind: int = 2, derive: bool = True) -> dict:
    """SPEC section 7 step 3 for a kind 2 (eQTL) or kind 3 (sQTL) block, then the section 5 values.
    With `derive=False` only the rules are checked and the header fields and credible sets returned."""
    if len(buf) < 64:
        raise PackError(f"{what}: shorter than the block header")
    magic, blen, n, var_start, anchor, pos_first, pos_last, n_cs, nlp_max, lse_min, lse_max, dz, dl = _BLOCK_H.unpack_from(buf, 0)
    if magic != b"QGB0" or blen != len(buf):
        raise PackError(f"{what}: magic {magic!r}, length field {blen}, range {len(buf)}")
    body = 64 + 4 * n + 12 * n_cs + dz
    if (body + 3) // 4 * 4 != blen or any(buf[body:]):
        raise PackError(f"{what}: body {body} padded to 4 != {blen}, or nonzero padding")
    if kind == 3:
        if dz or dl or n == 0:
            raise PackError(f"{what}: a kind 3 block has no details frame and at least one row (details_zlen {dz}, details_len {dl}, n_rows {n})")
    elif dz == 0:
        raise PackError(f"{what}: a kind 2 block needs a details frame")
    if n != (want_n or 0):
        raise PackError(f"{what}: n_rows {n} != index n_var {want_n}")
    if n == 0:
        if (var_start, anchor, pos_first, pos_last, n_cs, nlp_max, lse_min, lse_max) != (0xFFFFFFFF, 0, 0, 0, 0, 0.0, 0.0, 0.0) or want_var_start is not None:
            raise PackError(f"{what}: an empty block needs the SPEC empty values")
    elif var_start != want_var_start or not 1 <= pos_first <= pos_last or var_start + n - 1 >= 0xFFFFFFFF:
        raise PackError(f"{what}: var_start {var_start} (index {want_var_start}), positions {pos_first}..{pos_last}")
    if not (math.isfinite(nlp_max) and nlp_max >= 0 and math.isfinite(lse_min) and math.isfinite(lse_max) and lse_min <= lse_max):
        raise PackError(f"{what}: scales nlp_max {nlp_max}, lse {lse_min}..{lse_max}")
    codes = np.frombuffer(buf, "<u2", 2 * n, 64).reshape(n, 2)
    nq, sq = codes[:, 0].astype(np.int64), codes[:, 1].astype(np.int64)
    fin = nq[nq <= 65533]
    if (nlp_max > 0 and (fin.size == 0 or fin.max() != 65533)) or (nlp_max == 0 and np.any(fin != 0)):
        raise PackError(f"{what}: -log10 p codes break the scale rule")
    lq = sq[sq != 0xFFFF] & 0x7FFF
    if np.any(lq == 0x7FFF):
        raise PackError(f"{what}: SE code 0x7FFF")
    if (lq.size == 0 and (lse_min, lse_max) != (0.0, 0.0)) or (lq.size and lse_min == lse_max and np.any(lq != 0)) \
            or (lse_min < lse_max and (lq.size == 0 or lq.min() != 0 or lq.max() != 32766)):
        raise PackError(f"{what}: log(SE) codes break the scale rule")
    cs = np.frombuffer(buf, _CS, n_cs, 64 + 4 * n)
    ordered = all((int(cs["row"][i]), int(cs["cs_id"][i])) < (int(cs["row"][i + 1]), int(cs["cs_id"][i + 1]))
                  for i in range(n_cs - 1))
    if n_cs and (not ordered or cs["row"][-1] >= n or np.any(cs["pad"] != 0)
                 or not np.all((cs["pip"] >= 0) & (cs["pip"] <= 1)) or np.any(cs["cs_id"] > 127)):
        raise PackError(f"{what}: credible-set records break SPEC")
    details = None
    if kind == 2:
        raw = _unframe(buf[64 + 4 * n + 12 * n_cs:body], dl, f"{what} details")
        if raw.startswith(b"\xef\xbb\xbf"):
            raise PackError(f"{what}: details JSON has a BOM")
        try:
            details = json.loads(raw.decode("utf-8"), parse_constant=_reject)
        except (UnicodeDecodeError, json.JSONDecodeError) as e:
            raise PackError(f"{what}: details JSON: {e}") from None
        if not isinstance(details, dict) or details.get("v") != 0:
            raise PackError(f"{what}: details are not an object with v = 0")
    if not derive:
        return {"n": n, "var_start": var_start, "anchor": anchor, "pos_first": pos_first, "pos_last": pos_last,
                "n_cs": n_cs, "cs": cs}

    nlp = nq * (nlp_max / 65533)
    nlp[nq == 65534] = np.inf
    nlp[nq == 65535] = np.nan
    p = np.power(10.0, -nlp)
    se_null = sq == 0xFFFF
    se = np.exp(lse_min + (sq & 0x7FFF) * ((lse_max - lse_min) / 32766))
    se[se_null] = np.nan
    negative = ((sq & 0x8000) != 0) & ~se_null
    with np.errstate(invalid="ignore"):
        t = np.maximum(-stdtrit(dof, p / 2), 0.0)
    slope = np.where(negative, -1.0, 1.0) * se * t
    slope[np.isnan(p) | (p == 0) | se_null] = np.nan
    return {"n": n, "var_start": var_start, "anchor": anchor, "pos_first": pos_first, "pos_last": pos_last,
            "nlp_max": nlp_max, "lse_min": lse_min, "lse_max": lse_max, "nq": nq, "sq": sq, "nlp": nlp, "p": p,
            "se": se, "negative": negative, "slope": slope, "cs": cs, "details": details}


def _read_pages(buf: bytes, what: str, n_cis: int | None = None) -> dict:
    """SPEC section 7 step 4: consecutive pages from byte 0 of `buf`. Flags bit 2 marks a record whose alleles are
    not reported: allele code 0, no heap record, A1 and A2 returned as None. With `n_cis` (the variants file header's
    value, for a range that may hold both sections): no page straddles n_cis, bit 2 appears only at vidx >= n_cis, and
    positions may restart at vidx n_cis; without it positions never decrease over the whole range."""
    off, next_vidx, codec0 = 0, None, None
    firsts, sizes, starts = [], [], []
    cols = {k: [] for k in ("pos", "rs", "af", "ms", "mc", "flags")}
    a1: list[str | None] = []
    a2: list[str | None] = []
    while off < len(buf):
        w = f"{what}, page at byte {off}"
        if len(buf) - off < 12:
            raise PackError(f"{w}: header truncated")
        stored, first, n, codec, reserved = _PAGE_H.unpack_from(buf, off)
        if reserved or codec not in (0, 1) or (codec0 is not None and codec != codec0) or n == 0 \
                or (next_vidx is not None and first != next_vidx):
            raise PackError(f"{w}: header stored {stored} first {first} n {n} codec {codec} reserved {reserved}")
        if n_cis is not None and first < n_cis < first + n:
            raise PackError(f"{w}: page {first}..{first + n - 1} straddles n_cis {n_cis}")
        codec0 = codec
        end = off + 12 + stored
        stop = (end - off + 3) // 4 * 4 + off
        if stop > len(buf) or any(buf[end:stop]):
            raise PackError(f"{w}: runs past the range or nonzero padding")
        payload = buf[off + 12:end] if codec == 0 else _unframe(buf[off + 12:end], None, w)
        (heap_len,) = struct.unpack_from("<I", payload, 0)
        if len(payload) != 4 + 16 * n + heap_len:
            raise PackError(f"{w}: payload {len(payload)} != 4 + 16n + heap_len")
        allele = np.frombuffer(payload, "u1", n, 4 + 14 * n)
        flags = np.frombuffer(payload, "u1", n, 4 + 15 * n)
        no_al = (flags & 4) != 0
        if np.any(allele > 12) or np.any(flags & 0xF8) or np.any((flags & 3) == 3):
            raise PackError(f"{w}: reserved allele code or flags")
        if np.any(no_al & (allele != 0)) or (n_cis is not None and first < n_cis and np.any(no_al)):
            raise PackError(f"{w}: flags bit 2 (alleles not reported) on a record with an allele code, or below n_cis {n_cis}")
        heap = payload[4 + 16 * n:].decode("ascii")
        recs = heap.split("\n")
        if recs[-1] != "":
            raise PackError(f"{w}: heap does not end in a newline")
        recs = iter(recs[:-1])
        for c, missing in zip(allele.tolist(), no_al.tolist()):
            if missing:
                x = y = None
            elif c:
                x, y = _SNP[c]
            else:
                parts = next(recs, "").split("\t")
                if len(parts) != 2 or (parts[0], parts[1]) in _SNP:
                    raise PackError(f"{w}: bad heap record {parts!r}")
                x, y = parts
            a1.append(x)
            a2.append(y)
        if next(recs, None) is not None:
            raise PackError(f"{w}: heap holds more records than code-0 alleles")
        cols["pos"].append(np.cumsum(np.frombuffer(payload, "<u4", n, 4), dtype=np.int64))
        cols["rs"].append(np.frombuffer(payload, "<u4", n, 4 + 4 * n))
        cols["af"].append(np.frombuffer(payload, "<u2", n, 4 + 8 * n))
        cols["ms"].append(np.frombuffer(payload, "<u2", n, 4 + 10 * n))
        cols["mc"].append(np.frombuffer(payload, "<u2", n, 4 + 12 * n))
        cols["flags"].append(flags)
        firsts.append(first)
        sizes.append(n)
        starts.append(off)
        next_vidx = first + n
        off = stop
    if not firsts:
        raise PackError(f"{what}: no pages")
    c = {k: np.concatenate(v) for k, v in cols.items()}
    down = np.diff(c["pos"]) < 0
    if n_cis is not None:
        down &= firsts[0] + np.arange(1, len(c["pos"])) != n_cis      # positions restart at the first trans-only variant
    if c["pos"].min() < 1 or np.any(down):
        raise PackError(f"{what}: positions below 1 or decreasing within a section")
    return {"first": firsts[0], "page_firsts": firsts, "page_sizes": sizes, "page_starts": starts + [off],
            "position": c["pos"], "rs_number": c["rs"], "af_code": c["af"], "ms": c["ms"], "mc": c["mc"],
            "match": c["flags"] & 3, "no_alleles": (c["flags"] & 4) != 0, "A1": a1, "A2": a2}


def _nan_to_none(x):
    if isinstance(x, dict):
        return {k: _nan_to_none(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_nan_to_none(v) for v in x]
    if isinstance(x, float) and math.isnan(x):
        return None
    return x


def _fmt_p(x: float) -> str:
    """ui/src/lib/format.ts fmtP"""
    if x is None or math.isnan(x):
        return ""
    if x == 0:
        return "0"
    return f"{x:#.2g}" if x >= 0.001 else f"{x:.1e}"


def _fmt3(x: float) -> str:
    """fmtNum(x, 3) of a FLOAT column"""
    return "" if x is None or math.isnan(x) else f"{float(np.float32(x)):.3f}"


def _write_arrow(table: pa.Table, path: Path) -> None:
    with pa.OSFile(str(path), "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)


def _details_reference(cfg: Config, chrom: str) -> dict[str, dict]:
    """Details for every binned gene of `chrom`, rebuilt without packfmt or SQL: the genes.parquet row
    minus chr and bin, exons merged in Python, splice_phenotypes rows in (cluster_id, intron_start,
    intron_end) order; NaN -> None."""
    gt = pq.read_table(cfg.tables / "genes.parquet", filters=[("chr", "=", chrom)])
    gkeys = [c for c in gt.column_names if c not in ("chr", "bin")]
    genes = {r["gene_id"]: {k: r[k] for k in gkeys} for r in gt.to_pylist() if r["bin"] is not None}
    ivs: dict[str, list[tuple[int, int]]] = {}
    for r in pq.read_table(cfg.tables / "exons.parquet", columns=["gene_id", "start", "end"], filters=[("gene_id", "in", sorted(genes))]).to_pylist():
        ivs.setdefault(r["gene_id"], []).append((r["start"], r["end"]))
    exons: dict[str, list[list[int]]] = {}
    for gid, iv in ivs.items():
        merged: list[list[int]] = []
        for a, b in sorted(iv):
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        exons[gid] = merged
    st = pq.read_table(cfg.tables / "splice_phenotypes.parquet", filters=[("chr", "=", chrom)])
    skeys = [c for c in st.column_names if c not in ("gene_id", "symbol", "chr", "tss")]
    blocks = {r["phenotype_id"]: r for r in pq.read_table(cfg.tmp / "pack_pointers" / f"sqtl_{chrom}.parquet",
                                                            columns=["phenotype_id", "blk_off", "blk_len"]).to_pylist()}
    splice: dict[str, list[dict]] = {}
    for r in sorted(st.to_pylist(), key=lambda r: (r["cluster_id"], r["intron_start"], r["intron_end"])):
        b = blocks.get(r["phenotype_id"], {})
        splice.setdefault(r["gene_id"], []).append({**{k: r[k] for k in skeys}, "blk_off": b.get("blk_off"), "blk_len": b.get("blk_len")})
    return {gid: _nan_to_none({"v": 0, "gene": g, "exons": exons.get(gid, []), "splice": splice.get(gid, [])}) for gid, g in genes.items()}


def validate(cfg: Config, con, check) -> None:
    d = cfg.derived
    pk = cfg["packs"]
    dof = int(pk["dof"]["eqtl"])
    out_dir = cfg.tmp / "pack_check"
    out_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET memory_limit = '{pk['duckdb_memory_limit']}'")
    con.execute(f"SET threads = {pk['duckdb_threads']}")
    register_search_index(cfg, con, "si_all")
    cols = read_search_index(cfg).column_names
    if not all(c in cols for c in PACK_COLUMNS):
        check(False, f"search_index has the pack columns {PACK_COLUMNS} (run `build --step pack_eqtl --step search_index`)")
        return
    chroms = list(CHROMS)
    # a gene outside the pack chromosomes, or one the eQTL pack has no block for, must have every
    # pack column null. `blk_off IS NULL` is the "no gene block" test now that `bin` has left the index.
    stray = con.execute(f"""SELECT count(*) FROM si_all WHERE (chr NOT IN ({', '.join(repr(c) for c in chroms)}) OR blk_off IS NULL)
        AND NOT ({' AND '.join(f'{c} IS NULL' for c in PACK_COLUMNS)})""").fetchone()[0]
    check(stray == 0, f"pack columns are null for genes without a block ({stray} rows are not)")
    tol = {"bands": BANDS[:-1], "bound_factor": packfmt.BOUND_FACTOR, "nlp_half_step_divisor": 131066, "nlp_slack": 1e-9,
           "af_tol": 0.5 / 65534 + 1e-12, "formula_rel": 1e-6, "formula_min_t": 0.05, "flat_nlp": 0.0005,
           "t_grid": {"rel": 1e-9, "abs": 1e-12, "small_t": 1e-3}}
    (out_dir / "tolerances.json").write_text(json.dumps(tol, indent=1))
    t_grid_parts = []
    summary = {"sqtl_totals": {"blocks": 0, "rows": 0, "records": 0}, "sqtl_rt": [],
               "sqtl_ref": {"files": {}, "genes": [], "introns": [], "rows": [], "bounds": []}}
    for chrom in chroms:
        try:
            t_grid_parts += _validate_chrom(cfg, con, check, chrom, dof, out_dir, tol, summary)
        except PackError as e:
            check(False, f"{chrom} pack decodes under SPEC: {e}")
    _validate_sqtl_summary(cfg, check, summary, out_dir)
    _validate_gwas(cfg, check, out_dir)
    # inverse-t reference stays within the format's explicit -log10(p) ceiling.
    grid = np.r_[0.0, np.logspace(-8, math.log10(packfmt.NLP_LIMIT), 1500)]
    nl, df = [np.concatenate([grid, grid])], [np.r_[np.full(grid.size, 435), np.full(grid.size, 480)]]
    for extra in t_grid_parts:
        nl.append(extra)
        df.append(np.full(extra.size, dof))
    nlp_in, dofs = np.concatenate(nl), np.concatenate(df).astype(np.float64)
    half = np.power(10.0, -nlp_in) / 2
    p = 2 * half
    with np.errstate(divide="ignore"):
        nlp_col = np.where(p > 0.5, -np.log1p(p - 1) / math.log(10), -np.log10(p))
    nlp_col = np.maximum(nlp_col, 0.0)
    t = np.maximum(-stdtrit(dofs, half), 0.0)
    ok = np.isfinite(t) & np.isfinite(nlp_col)
    _write_arrow(pa.table({"nlp": nlp_col[ok], "dof": dofs[ok], "t": t[ok]}), out_dir / "t_grid.arrow")
    log(f"pack validate: reference files for `npm run pack-check` in {out_dir.relative_to(cfg.derived)}/")


def _validate_chrom(cfg: Config, con, check, chrom: str, dof: int, out_dir: Path, tol: dict, summary: dict) -> list[np.ndarray]:
    d = cfg.derived
    pk = cfg["packs"]
    var_path, eqtl_path = pack_paths(cfg, chrom)
    if not (var_path.exists() and eqtl_path.exists()):
        check(False, f"{chrom}: pack files exist ({var_path.relative_to(d)}, {eqtl_path.relative_to(d)})")
        return []
    pack = eqtl_path.read_bytes()
    vfile = var_path.read_bytes()
    idx = con.execute(f"""SELECT gene_id, symbol, chr, tss, tested, {', '.join(PACK_COLUMNS)}
        FROM si_all WHERE chr = ? AND blk_off IS NOT NULL ORDER BY tss, gene_id""", [chrom]).fetchall()
    names = ["gene_id", "symbol", "chr", "tss", "tested"] + PACK_COLUMNS
    genes = [dict(zip(names, r)) for r in idx]

    # ---- structure, every gene with a bin ----
    n_detail = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'genes.parquet'}' WHERE chr = ? AND bin IS NOT NULL", [chrom]).fetchone()[0]
    nonnull = sum(all(g[c] is not None for c in ("blk_off", "blk_len", "var_off", "var_len", "w_lo", "w_hi"))
                  and ((g["var_start"] is None) == (not g["tested"]))
                  and ((g["n_var"] is None) == (g["var_start"] is None)) for g in genes)
    count, zero, _ = _file_header(pack, 2, chrom, str(eqtl_path.name))
    check(nonnull == len(genes) == n_detail == count and zero == 0,
          f"{chrom}: search_index pack columns set for {nonnull} of {len(genes)} genes with a bin; genes.parquet binned genes {n_detail}; pack header blocks {count}")
    off, contiguous = 32, True
    for g in genes:
        contiguous &= g["blk_off"] == off and g["blk_off"] % 4 == 0
        off += g["blk_len"] or 0
    check(contiguous and off == len(pack), f"{chrom}: blocks are contiguous in TSS order from byte 32 to the file end ({off} vs {len(pack):,} bytes)")

    n_vars, page_size, hdr_n_cis = _file_header(vfile, 1, chrom, str(var_path.name))
    n_cis = con.execute(f"SELECT count(*) FROM {variants_sql(cfg, chrom)} WHERE in_cis").fetchone()[0]
    pages = _read_pages(vfile[32:], str(var_path.name), n_cis=hdr_n_cis)
    cis_pages = -(-hdr_n_cis // page_size)
    check(hdr_n_cis == n_cis and n_vars == len(pages["position"]) and page_size == pk["variant_page_size"] and pages["first"] == 0
          and all(s == page_size for s in pages["page_sizes"][:cis_pages - 1]),
          f"{chrom}: variants file holds {hdr_n_cis:,} cis variants (in_cis rows {n_cis:,}) in {cis_pages} pages of {page_size}, then "
          f"{n_vars - hdr_n_cis:,} trans-only variants (header count {n_vars:,}, decoded {len(pages['position']):,})")
    n_vars = hdr_n_cis      # every gene and intron run lies in the cis section
    page_starts = [s + 32 for s in pages["page_starts"]]
    raw_e = cfg.raw_dir("cis_eQTL_nominal") / f"topchef_{chrom}_MaxPC70.cis_qtl_pairs.{chrom}.parquet"
    con.execute(f"""CREATE OR REPLACE TABLE vv AS
        SELECT (row_number() OVER (ORDER BY position, A1, A2) - 1)::INTEGER vidx, position, A1, A2, rs_number
        FROM {variants_sql(cfg, chrom)} WHERE in_cis""")
    intron_runs: dict[str, list[tuple[str, int, int, int]]] = {}
    for gid, ph, lo, hi, cnt in con.execute(f"""SELECT s.gene_id, r.phenotype_id, min(v.vidx), max(v.vidx) + 1, count(*)
            FROM '{_raw_sqtl(cfg, chrom)}' r JOIN vv v ON v.position = r.position AND v.A1 = r.A1 AND v.A2 = r.A2
            JOIN '{cfg.tables / 'splice_phenotypes.parquet'}' s ON s.phenotype_id = r.phenotype_id WHERE s.chr = ? GROUP BY 1, 2""", [chrom]).fetchall():
        intron_runs.setdefault(gid, []).append((ph, int(lo), int(hi), int(cnt)))
    intron_ranges = {g: (min(r[1] for r in rs), max(r[2] for r in rs)) for g, rs in intron_runs.items()}

    nominal_n = dict(con.execute(f"SELECT phenotype_id, count(*) FROM '{raw_e}' GROUP BY 1").fetchall())
    num_var = dict(con.execute(f"SELECT gene_id, num_var FROM '{cfg.tables / 'genes.parquet'}' WHERE chr = ? AND bin IS NOT NULL", [chrom]).fetchall())
    want_details = _details_reference(cfg, chrom)
    blocks: dict[str, dict] = {}
    bad_n, bad_range, unknown, bad_detail, frames = [], [], [], [], 0
    steep = []
    for g in genes:
        gid = g["gene_id"]
        blk = _read_block(pack[g["blk_off"]:g["blk_off"] + g["blk_len"]], dof, g["n_var"], g["var_start"], f"{eqtl_path.name} {gid}")
        frames += 1
        blocks[gid] = blk
        if blk["n"] != nominal_n.get(gid, 0):
            bad_n.append(gid)
        if blk["nlp_max"] / 131066 > tol["flat_nlp"]:
            steep.append((blk["nlp_max"] / 131066, gid))
        run = (g["var_start"], g["var_start"] + g["n_var"]) if blk["n"] else None
        if gid in intron_ranges:
            lo, hi = intron_ranges[gid]
            run = (lo, hi) if run is None else (min(run[0], lo), max(run[1], hi))
        if run is None or run[1] > n_vars:
            bad_range.append(gid)
        else:
            u0, u1 = run
            span_ok = (g["var_off"] == page_starts[u0 // page_size]
                       and g["var_off"] + g["var_len"] == page_starts[(u1 - 1) // page_size + 1] <= len(vfile)
                       and pages["position"][u0] == g["w_lo"] and pages["position"][u1 - 1] == g["w_hi"])
            if blk["n"]:
                span_ok = (span_ok and pages["position"][g["var_start"]] == blk["pos_first"]
                           and pages["position"][g["var_start"] + g["n_var"] - 1] == blk["pos_last"])
            if not span_ok:
                bad_range.append(gid)
            if any(np.any(pages[k][u0:u1] == 65535) for k in ("af_code", "ms", "mc")):
                unknown.append(gid)
        want, got = want_details[gid], blk["details"]
        order_ok = (list(got) == ["v", "gene", "exons", "splice"] and list(got["gene"]) == list(want["gene"])
                    and all(list(a) == list(b) for a, b in zip(got["splice"], want["splice"])))
        if got != want or not order_ok:
            bad_detail.append(gid)
    differs_num_var = [gid for gid in blocks if blocks[gid]["n"] and blocks[gid]["n"] != num_var.get(gid)]
    check(not bad_n, f"{chrom}: block pair counts equal raw eQTL row counts for {len(genes) - len(bad_n)} of {len(genes)} genes {bad_n[:5]}")
    log(f"  {chrom}: {len(differs_num_var)} genes have one row more than permutation num_var (the null-p row), e.g. {differs_num_var[:3]}")
    check(frames == len(genes), f"{chrom}: every block passes SPEC section 7 (magic, length, padding, scale rule, credible sets, one zstd details frame, JSON without NaN, v = 0): {frames} blocks")
    check(not bad_range, f"{chrom}: every gene has a union run (eQTL and intron runs) within {n_vars:,} variants, var_off/var_len are exactly its covering pages, w_lo/w_hi are its first and last positions, and pos_first/pos_last match ({len(bad_range)} bad {bad_range[:5]})")
    check(not unknown, f"{chrom}: no variant inside a gene's union run has the not-known af or count code ({len(unknown)} genes {unknown[:5]})")
    check(not bad_detail, f"{chrom}: details JSON equals an independent rebuild from genes, exons (Python interval merge), and splice_phenotypes (NaN -> null, key order) for {len(genes) - len(bad_detail)} of {len(genes)} genes {bad_detail[:5]}")
    cs_total = sum(len(b["cs"]) for b in blocks.values())
    want_cs = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 'e' AND chr = ?", [chrom]).fetchone()[0]
    check(cs_total == want_cs, f"{chrom}: credible-set records in the pack ({cs_total:,}) equal the credible_sets eQTL rows ({want_cs:,})")
    steep.sort(reverse=True)
    log(f"  {chrom}: {len(steep)} genes have half a -log10 p step above {tol['flat_nlp']}: worst {steep[0][0]:.5f} ({steep[0][1]})" if steep
        else f"  {chrom}: no gene has half a -log10 p step above {tol['flat_nlp']}")

    by_id = {g["gene_id"]: g for g in genes}
    _validate_sqtl_chrom(cfg, con, check, chrom, {
        "pages": pages, "page_starts": page_starts, "page_size": page_size, "n_vars": n_vars, "by_id": by_id,
        "blocks": blocks, "intron_runs": intron_runs, "vfile": vfile, "var_path": var_path, "tol": tol, "summary": summary})

    # ---- 3b round trip on sample genes ----
    # `bin` has left the search index, so FLNC's bin-mates become its TSS neighbours: the same
    # intent (a run of adjacent genes, not scattered ones) from what the index still carries
    order = [g["gene_id"] for g in genes]           # already in (tss, gene_id) order
    neighbours: set[str] = set()
    if "ENSG00000128591" in order:
        i = order.index("ENSG00000128591")
        neighbours = set(order[max(0, i - 50):i + 50])
    tested_ids = [g["gene_id"] for g in genes if g["tested"]]
    random_ids = [str(x) for x in np.random.default_rng(3).choice(tested_ids, min(2, len(tested_ids)), replace=False)]
    sample = sorted({g["gene_id"] for g in genes if g["gene_id"] in set(pk["check_genes"])} | neighbours | set(random_ids))
    in_list = ', '.join(repr(x) for x in sample)
    # SPEC Q1: af from float32; pip/cs_id: the higher-PIP membership, lower cs_id on a tie (section 8)
    src = con.execute(f"""WITH cs AS (
            SELECT phenotype_id, position, A1, A2, pip, cs_id,
                   row_number() OVER (PARTITION BY phenotype_id, position, A1, A2 ORDER BY pip DESC, cs_id) AS rk
            FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 'e' AND chr = ? AND phenotype_id IN ({in_list}))
        SELECT r.phenotype_id AS gene_id, r.position, r.A1, r.A2, v.rs_number, r.start_distance::INTEGER AS tss_distance,
               r.af::FLOAT AS af, r.ma_samples::SMALLINT AS ma_samples, r.ma_count::SMALLINT AS ma_count,
               r.pval_nominal, r.slope::FLOAT AS slope, r.slope_se::FLOAT AS slope_se, cs.pip, cs.cs_id
        FROM '{raw_e}' r
        LEFT JOIN vv v ON v.position = r.position AND v.A1 = r.A1 AND v.A2 = r.A2
        LEFT JOIN cs ON cs.rk = 1 AND cs.phenotype_id = r.phenotype_id AND cs.position = r.position AND cs.A1 = r.A1 AND cs.A2 = r.A2
        WHERE r.phenotype_id IN ({in_list})
        ORDER BY gene_id, r.position, r.A1, r.A2""", [chrom]).fetch_arrow_table()
    sgid = src["gene_id"].to_numpy(zero_copy_only=False)
    S = {k: src[k].to_numpy(zero_copy_only=False) for k in ("position", "rs_number", "tss_distance", "af", "ma_samples", "ma_count",
                                                            "pval_nominal", "slope", "slope_se", "pip", "cs_id")}
    S["A1"], S["A2"] = src["A1"].to_pylist(), src["A2"].to_pylist()
    dec = {k: [] for k in ("position", "rs", "tss", "af", "ms", "mc", "nlp", "p", "se", "slope", "negative", "pip", "cs_id",
                           "nlp_max", "se_b", "sl_b", "nq", "sq")}
    dA1, dA2, count_bad = [], [], []
    got_cs: list[tuple] = []
    for gid in sample:
        g, blk = by_id[gid], blocks[gid]
        n = blk["n"]
        if (sgid == gid).sum() != n:
            count_bad.append(gid)
            continue
        if n == 0:
            continue
        rng = _read_pages(vfile[g["var_off"]:g["var_off"] + g["var_len"]], f"{var_path.name} {gid}")
        i0 = g["var_start"] - rng["first"]
        if i0 < 0 or i0 + n > len(rng["position"]):
            raise PackError(f"{var_path.name} {gid}: range pages do not cover the gene's rows")
        sl = slice(i0, i0 + n)
        dec["position"].append(rng["position"][sl])
        dA1 += rng["A1"][sl]
        dA2 += rng["A2"][sl]
        dec["rs"].append(rng["rs_number"][sl].astype(np.int64))
        dec["tss"].append(rng["position"][sl] - blk["anchor"])
        dec["af"].append(rng["af_code"][sl])
        dec["ms"].append(rng["ms"][sl])
        dec["mc"].append(rng["mc"][sl])
        for k in ("nlp", "p", "se", "slope", "negative", "nq", "sq"):
            dec[k].append(blk[k])
        pip = np.full(n, np.nan, dtype=np.float32)
        cs_id = np.full(n, -1, dtype=np.int64)
        for row, value, cid in zip(blk["cs"]["row"], blk["cs"]["pip"], blk["cs"]["cs_id"]):
            if np.isnan(pip[row]) or value > pip[row] or (value == pip[row] and cid < cs_id[row]):
                pip[row], cs_id[row] = value, cid
        got_cs += [(gid, int(rng["position"][i0 + r]), rng["A1"][i0 + r], rng["A2"][i0 + r], int(c), float(v))
                   for r, v, c in zip(blk["cs"]["row"].tolist(), blk["cs"]["pip"], blk["cs"]["cs_id"].tolist())]
        dec["pip"].append(pip)
        dec["cs_id"].append(cs_id)
        dec["nlp_max"].append(np.full(n, blk["nlp_max"]))
        m = sgid == gid
        se_b, sl_b = packfmt.error_bounds(S["pval_nominal"][m], S["slope"][m], S["slope_se"][m], blk["nlp_max"], blk["lse_min"], blk["lse_max"], dof)
        dec["se_b"].append(se_b)
        dec["sl_b"].append(sl_b)
    want_cs_rows = sorted((a, int(b), c, e, int(f), float(g)) for a, b, c, e, f, g in con.execute(
        f"""SELECT phenotype_id, position, A1, A2, cs_id, pip FROM '{cfg.tables / 'credible_sets.parquet'}'
        WHERE qtl_type = 'e' AND chr = ? AND phenotype_id IN ({in_list})""", [chrom]).fetchall())
    check(sorted(got_cs) == want_cs_rows, f"{chrom}: credible-set (position, A1, A2, cs_id, pip) of the {len(sample)} sample genes equal credible_sets.parquet "
          f"({len(got_cs)} decoded, {len(want_cs_rows)} in the table)")
    s_ph = [(gid, *r) for gid in sample for r in intron_runs.get(gid, [])]
    cov_bad: list[str] = []
    if s_ph:
        rt = con.execute(f"""SELECT phenotype_id, position, A1, A2 FROM '{_raw_sqtl(cfg, chrom)}'
            WHERE phenotype_id IN ({', '.join(repr(x) for x in sorted({r[1] for r in s_ph}))}) ORDER BY phenotype_id, position, A1, A2""").fetch_arrow_table()
        rph = rt["phenotype_id"].to_numpy(zero_copy_only=False)
        rpos, rA1, rA2 = rt["position"].to_numpy(), rt["A1"].to_pylist(), rt["A2"].to_pylist()
        ranges: dict[str, dict] = {}
        for gid, ph, lo, hi, cnt in s_ph:
            g = by_id[gid]
            if gid not in ranges:
                ranges[gid] = _read_pages(vfile[g["var_off"]:g["var_off"] + g["var_len"]], f"{var_path.name} {gid}")
            rng = ranges[gid]
            m = np.flatnonzero(rph == ph)
            a, b = lo - rng["first"], hi - rng["first"]
            ok = (cnt == hi - lo == len(m) and a >= 0 and b <= len(rng["position"])
                  and np.array_equal(rng["position"][a:b], rpos[m]) and rng["A1"][a:b] == [rA1[i] for i in m]
                  and rng["A2"][a:b] == [rA2[i] for i in m])
            if not ok:
                cov_bad.append(ph)
    check(not cov_bad, f"{chrom}: intron coverage, {len(s_ph)} introns of the sample genes: every raw sQTL row matches the decoded variant "
          f"at its run, inside the gene's range ({len(cov_bad)} bad {cov_bad[:3]})")
    check(not count_bad, f"{chrom}: round-trip sample of {len(sample)} genes: block row counts equal raw eQTL row counts ({count_bad[:5]})")
    D = {k: np.concatenate(v) if v else np.array([]) for k, v in dec.items()}
    keep = np.isin(sgid, [x for x in sample if x not in count_bad])
    S = {k: (np.asarray(v, dtype=object)[keep].tolist() if k in ("A1", "A2") else v[keep]) for k, v in S.items()}
    sgid = sgid[keep]
    nrow = len(sgid)
    if _compare_rows(check, chrom, f"{nrow:,} rows of {len(sample)} genes", D, dA1, dA2, S, dof, tol) is None:
        return []

    # ---- 3c reference files ----
    _write_arrow(src.filter(pa.array(keep)), out_dir / f"{chrom}_rows.arrow")
    _write_arrow(pa.table({"gene_id": pa.array(sgid.tolist()), "position": S["position"].astype(np.int32), "A1": S["A1"], "A2": S["A2"],
                           "se_bound": D["se_b"], "slope_bound": D["sl_b"]}), out_dir / f"{chrom}_bounds.arrow")
    (out_dir / f"{chrom}_index.json").write_text(json.dumps({
        "chrom": chrom, "dof": dof, "variant_page_size": pk["variant_page_size"],
        "files": {"eqtl": str(eqtl_path.relative_to(d)), "variants": str(var_path.relative_to(d))},
        "rows": [by_id[gid] for gid in sample]}, indent=1))
    (out_dir / f"{chrom}_details.json").write_text(json.dumps({gid: want_details[gid] for gid in sample}))
    con.execute("DROP TABLE IF EXISTS vv")
    flnc = blocks.get("ENSG00000128591")
    return [np.unique(flnc["nlp"][np.isfinite(flnc["nlp"])])] if flnc and flnc["n"] else []


def _compare_rows(check, label: str, scope: str, D: dict, dA1: list, dA2: list, S: dict, dof: int, tol: dict,
                  af_exact: bool = False) -> dict | None:
    """Decoded rows `D` against source rows `S` (same order): exact fields, then af, -log10 p, SE, and
    slope within SPEC section 9's limits. Returns the worst errors, or None when an exact field differs."""

    def isnull(a):
        return np.isnan(a.astype(np.float64))

    rs_src = np.where(isnull(S["rs_number"]), 0, np.nan_to_num(S["rs_number"].astype(np.float64))).astype(np.int64)
    ms_null, mc_null, af_null = isnull(S["ma_samples"]), isnull(S["ma_count"]), isnull(S["af"])
    p_src, sl_src, se_src = S["pval_nominal"].astype(np.float64), S["slope"].astype(np.float64), S["slope_se"].astype(np.float64)
    pip_src = S["pip"].astype(np.float32)
    cs_src = np.where(isnull(S["cs_id"]), -1, np.nan_to_num(S["cs_id"].astype(np.float64))).astype(np.int64)
    exact = {
        "position, A1, A2": np.array_equal(D["position"], S["position"].astype(np.int64)) and dA1 == S["A1"] and dA2 == S["A2"],
        "rs_number": np.array_equal(D["rs"], rs_src),
        "tss_distance": np.array_equal(D["tss"], S["tss_distance"].astype(np.int64)),
        "ma_samples": np.array_equal(D["ms"] == 65535, ms_null) and np.array_equal(D["ms"][~ms_null], S["ma_samples"][~ms_null].astype(np.int64)),
        "ma_count": np.array_equal(D["mc"] == 65535, mc_null) and np.array_equal(D["mc"][~mc_null], S["ma_count"][~mc_null].astype(np.int64)),
        "pip": np.array_equal(np.isnan(D["pip"]), np.isnan(pip_src)) and np.array_equal(D["pip"][~np.isnan(pip_src)], pip_src[~np.isnan(pip_src)]),
        "cs_id": np.array_equal(D["cs_id"], cs_src),
        "null p": np.array_equal(D["nq"] == 65535, np.isnan(p_src)),
        "null SE code": np.array_equal(D["sq"] == 0xFFFF, np.isnan(se_src) | np.isnan(sl_src)),
        "null slope": np.array_equal(np.isnan(D["slope"]), np.isnan(sl_src) | np.isnan(p_src) | (p_src == 0) | np.isnan(se_src)),
        "null af": np.array_equal(D["af"] == 65535, af_null),
    }
    if af_exact:
        exact["af code"] = np.array_equal(D["af"][~af_null], np.rint(S["af"][~af_null].astype(np.float64) * 65534).astype(np.int64))
    bad = [k for k, ok in exact.items() if not ok]
    check(not bad, f"{label}: round trip exact against the raw Zenodo rows on {scope}: row keys, rs_number, tss_distance, counts, pip, cs_id, null pattern{', af code' if af_exact else ''} ({bad})")
    if bad:
        return None
    af_err = np.abs(D["af"][~af_null] / 65534 - S["af"][~af_null].astype(np.float64))
    ok_p = ~np.isnan(p_src) & (p_src > 0)
    nlp_err = np.abs(D["nlp"][ok_p] - (-np.log10(p_src[ok_p])))
    nlp_lim = D["nlp_max"][ok_p] / tol["nlp_half_step_divisor"] + tol["nlp_slack"]
    nn = ~np.isnan(D["se_b"])                                    # rows with p, slope, SE present and p > 0
    se_err = np.abs(D["se"] - se_src)
    sl_err = np.abs(D["slope"] - sl_src)
    r_se = np.where(nn, se_err / D["se_b"], 0.0)
    r_sl = np.where(nn, sl_err / D["sl_b"], 0.0)
    check(bool(af_err.max(initial=0) <= tol["af_tol"]), f"{label}: af within half a code step: max error {af_err.max(initial=0):.3g} (limit {tol['af_tol']:.3g})")
    check(bool(np.all(nlp_err <= nlp_lim)), f"{label}: -log10 p within half the gene's step: max error {nlp_err.max(initial=0):.3g}, max error / limit {np.max(nlp_err / nlp_lim, initial=0):.4f}")
    check(bool(r_se.max(initial=0) <= tol["bound_factor"] and r_sl.max(initial=0) <= tol["bound_factor"]),
          f"{label}: SE and slope within {tol['bound_factor']} x the per-row bound: max error / bound {r_se.max(initial=0):.4f} (SE), {r_sl.max(initial=0):.4f} (slope); "
          f"worst SE error {se_err[nn].max(initial=0):.3g}, worst slope error {sl_err[nn].max(initial=0):.3g}")
    present = ~np.isnan(se_src) & ~np.isnan(sl_src) & ~np.isnan(p_src) & (p_src > 0)
    nz = present & (sl_src != 0)
    sane = (np.all(np.isfinite(D["se"][present])) and np.all(D["se"][present] > 0) and np.all(np.isfinite(D["slope"][present]))
            and np.array_equal(np.signbit(D["slope"][nz]), sl_src[nz] < 0))
    check(bool(sane), f"{label}: every non-null row decodes to a finite SE > 0 and a finite slope with the source's sign ({int(present.sum()):,} rows)")
    t_src = np.abs(sl_src) / se_src
    ref = np.sign(sl_src) * se_src * np.maximum(-stdtrit(dof, p_src / 2), 0.0)
    fm = present & (t_src > tol["formula_min_t"])
    rel = np.abs(ref[fm] - sl_src[fm]) / np.abs(sl_src[fm])
    check(bool(rel.max(initial=0) <= tol["formula_rel"]), f"{label}: slope from the unquantized p and SE with dof {dof} matches the source: max relative error {rel.max(initial=0):.3g} where |t| > {tol['formula_min_t']}")

    log(f"  {label} worst errors: af {af_err.max(initial=0):.3g}, -log10 p {nlp_err.max(initial=0):.3g}, SE {se_err[nn].max(initial=0):.3g}, slope {sl_err[nn].max(initial=0):.3g}; other fields exact")
    for lo, hi in zip(BANDS[:-1], BANDS[1:]):
        m = nn & (t_src >= lo) & (t_src < hi)
        if m.any():
            log(f"  {label} |t| [{lo}, {hi}): {int(m.sum()):,} rows, SE max error {se_err[m].max():.3g} (error / bound {r_se[m].max():.4f}), "
                f"slope max error {sl_err[m].max():.3g} (error / bound {r_sl[m].max():.4f})")
    one = nn & (D["nq"] == 0)
    if one.any():
        log(f"  {label}: {int(one.sum())} rows read back as p = 1 (slope 0); worst slope error {sl_err[one].max():.3g}")
    ps = [_fmt_p(a) != _fmt_p(b) for a, b in zip(D["p"][present].tolist(), p_src[present].tolist())]
    ss = [_fmt3(a) != _fmt3(b) for a, b in zip(D["slope"][present].tolist(), sl_src[present].tolist())]
    es = [_fmt3(a) != _fmt3(b) for a, b in zip(D["se"][present].tolist(), se_src[present].tolist())]
    log(f"  {label}: printed differently from the source rows: p (fmtP) {np.mean(ps):.2%}, slope (3 decimals) {np.mean(ss):.2%}, SE (3 decimals) {np.mean(es):.2%}")
    return {"label": label, "rows": int(nn.sum()), "af": float(af_err.max(initial=0)), "nlp": float(nlp_err.max(initial=0)),
            "nlp_ratio": float(np.max(nlp_err / nlp_lim, initial=0)), "se": float(se_err[nn].max(initial=0)),
            "slope": float(sl_err[nn].max(initial=0)), "r_se": float(r_se.max(initial=0)), "r_sl": float(r_sl.max(initial=0)),
            "formula_rel": float(rel.max(initial=0))}


SQTL_REF_CHROMS = ("chr1", "chr2", "chr4", "chr10", "chrX")     # reference files for the browser decoder check (cd ui && npm run pack-check)


def _validate_sqtl_chrom(cfg: Config, con, check, chrom: str, ctx: dict) -> None:
    """The sQTL pack of one chromosome (SPEC section 5 "sQTL results pack", section 7 step 5): every
    block's structure, pointers and gene details against the file, runs against the raw rows, then a
    round trip of sample introns against the raw Zenodo rows."""
    d = cfg.derived
    pk = cfg["packs"]
    dof = int(pk["dof"]["sqtl"])
    tol, summary = ctx["tol"], ctx["summary"]
    path, ptr_path = sqtl_path(cfg, chrom), pointer_dir(cfg) / f"sqtl_{chrom}.parquet"
    if not (path.exists() and ptr_path.exists()):
        check(False, f"{chrom}: sQTL pack and pointer file exist ({path.relative_to(d)}, {ptr_path.name})")
        return
    ptrs = sorted(pq.read_table(ptr_path).to_pylist(), key=lambda r: r["blk_off"])
    by_ph = {r["phenotype_id"]: r for r in ptrs}
    raw = _raw_sqtl(cfg, chrom)
    n_raw = pq.read_metadata(raw).num_rows
    want_cs = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 's' AND chr = ?", [chrom]).fetchone()[0]
    splice_ids = {r[0] for r in con.execute(f"SELECT phenotype_id FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE chr = ?", [chrom]).fetchall()}
    pages, page_starts, P, by_id = ctx["pages"], ctx["page_starts"], ctx["page_size"], ctx["by_id"]
    with open(path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        try:
            # ---- structure, every block ----
            size = len(mm)
            count, zero, _ = _file_header(mm[:32], 3, chrom, path.name)
            walk, off = [], 32
            while off < size:
                if size - off < 64 or mm[off:off + 4] != b"QGB0":
                    raise PackError(f"{path.name}: no block magic at byte {off}")
                (ln,) = struct.unpack_from("<I", mm, off + 4)
                if ln < 64 or ln % 4 or off + ln > size:
                    raise PackError(f"{path.name}: block length {ln} at byte {off}")
                walk.append((off, ln))
                off += ln
            check(count == len(walk) == len(ptrs) == len(splice_ids) and zero == 0 and set(by_ph) == splice_ids,
                  f"{chrom}: sQTL pack header is kind 3 with {count:,} blocks; its length fields walk {len(walk):,} contiguous blocks from byte 32 "
                  f"to the file end ({size:,} bytes); pointer rows {len(ptrs):,}; splice_phenotypes introns {len(splice_ids):,}")
            check([(r["blk_off"], r["blk_len"]) for r in ptrs] == walk, f"{chrom}: sQTL pointer blk_off/blk_len equal the blocks walked in the file")
            raw_runs = {ph: (gid, lo, hi, cnt) for gid, rs in ctx["intron_runs"].items() for ph, lo, hi, cnt in rs}
            n_sum = cs_sum = 0
            hdr_bad, raw_bad, range_bad = [], [], []
            for r in ptrs:
                ph = r["phenotype_id"]
                blk = _read_block(mm[r["blk_off"]:r["blk_off"] + r["blk_len"]], dof, r["n_var"], r["var_start"], f"{path.name} {ph}",
                                  kind=3, derive=False)
                n, vs = blk["n"], blk["var_start"]
                n_sum += n
                cs_sum += blk["n_cs"]
                if (blk["anchor"], blk["n_cs"]) != (r["anchor"], r["n_cs"]):
                    hdr_bad.append(ph)
                if raw_runs.get(ph) != (r["gene_id"], vs, vs + n, n):
                    raw_bad.append(ph)
                g = by_id.get(r["gene_id"])
                ok = g is not None and g["var_off"] is not None and vs + n <= ctx["n_vars"]
                if ok:
                    ok = (g["var_off"] <= page_starts[vs // P] and page_starts[(vs + n - 1) // P + 1] <= g["var_off"] + g["var_len"]
                          and pages["position"][vs] == blk["pos_first"] and pages["position"][vs + n - 1] == blk["pos_last"])
                if not ok:
                    range_bad.append(ph)
            check(True, f"{chrom}: every sQTL block passes SPEC section 7 as kind 3 (magic, length, padding, no details frame, "
                        f"at least one row, scale rule, credible sets): {len(ptrs):,} blocks")
            check(n_sum == n_raw and cs_sum == want_cs, f"{chrom}: sQTL block rows {n_sum:,} equal the raw file's {n_raw:,}; "
                                                        f"credible-set records {cs_sum:,} equal the credible_sets sQTL rows {want_cs:,}")
            check(not hdr_bad, f"{chrom}: sQTL pointer anchor and n_cs equal the block headers ({len(hdr_bad)} bad {hdr_bad[:3]})")
            check(not raw_bad, f"{chrom}: every intron's block var_start and n_rows, and its gene, equal its raw rows' run in the cis variant list ({len(raw_bad)} bad {raw_bad[:3]})")
            check(not range_bad, f"{chrom}: every intron run lies inside its gene's var_off/var_len pages, and pos_first/pos_last equal the decoded positions ({len(range_bad)} bad {range_bad[:3]})")
            got = {(sp["phenotype_id"], sp.get("blk_off"), sp.get("blk_len"), gid) for gid, b in ctx["blocks"].items() for sp in b["details"]["splice"]}
            want = {(r["phenotype_id"], r["blk_off"], r["blk_len"], r["gene_id"]) for r in ptrs}
            check(got == want, f"{chrom}: gene details carry every intron's blk_off/blk_len, equal to the pointer file, under its gene "
                               f"({len(want):,} introns, {len(got ^ want)} differ)")
            T = summary["sqtl_totals"]
            T["blocks"] += len(walk)
            T["rows"] += n_sum
            T["records"] += cs_sum

            # ---- round trip on sample introns ----
            reasons: dict[str, str] = {}
            for r in ptrs:
                if r["gene_id"] in set(pk["sqtl_check_genes"]):
                    reasons[r["phenotype_id"]] = "gene"
            reasons.setdefault(max(ptrs, key=lambda r: (r["n_var"], r["phenotype_id"]))["phenotype_id"], "largest")
            for (ph,) in con.execute(f"""SELECT DISTINCT phenotype_id FROM (SELECT phenotype_id FROM '{cfg.tables / 'credible_sets.parquet'}'
                    WHERE qtl_type = 's' AND chr = ? GROUP BY phenotype_id, position, A1, A2 HAVING count(*) > 1)""", [chrom]).fetchall():
                reasons[ph] = "two sets"
            for ph in np.random.default_rng(100 + CHROMS.index(chrom)).choice(sorted(by_ph), min(5, len(by_ph)), replace=False).tolist():
                reasons.setdefault(ph, "random")
            sample = sorted(reasons)
            in_list = ", ".join(repr(x) for x in sample)
            src = con.execute(f"""WITH cs AS (
                    SELECT phenotype_id, position, A1, A2, pip, cs_id,
                           row_number() OVER (PARTITION BY phenotype_id, position, A1, A2 ORDER BY pip DESC, cs_id) AS rk
                    FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 's' AND chr = ? AND phenotype_id IN ({in_list}))
                SELECT r.phenotype_id, r.position, r.A1, r.A2, v.rs_number, r.start_distance::INTEGER AS tss_distance,
                       r.af::FLOAT AS af, r.ma_samples::SMALLINT AS ma_samples, r.ma_count::SMALLINT AS ma_count,
                       r.pval_nominal, r.slope::FLOAT AS slope, r.slope_se::FLOAT AS slope_se, cs.pip, cs.cs_id
                FROM '{raw}' r
                LEFT JOIN vv v ON v.position = r.position AND v.A1 = r.A1 AND v.A2 = r.A2
                LEFT JOIN cs ON cs.rk = 1 AND cs.phenotype_id = r.phenotype_id AND cs.position = r.position AND cs.A1 = r.A1 AND cs.A2 = r.A2
                WHERE r.phenotype_id IN ({in_list})
                ORDER BY r.phenotype_id, r.position, r.A1, r.A2""", [chrom]).fetch_arrow_table()
            sph = src["phenotype_id"].to_numpy(zero_copy_only=False)
            S = {k: src[k].to_numpy(zero_copy_only=False) for k in ("position", "rs_number", "tss_distance", "af", "ma_samples", "ma_count",
                                                                    "pval_nominal", "slope", "slope_se", "pip", "cs_id")}
            S["A1"], S["A2"] = src["A1"].to_pylist(), src["A2"].to_pylist()
            dec = {k: [] for k in ("position", "rs", "tss", "af", "ms", "mc", "nlp", "p", "se", "slope", "negative", "pip", "cs_id",
                                   "nlp_max", "se_b", "sl_b", "nq", "sq")}
            dA1, dA2, count_bad, got_cs = [], [], [], []
            ranges: dict[str, dict] = {}
            for ph in sample:
                r = by_ph[ph]
                g = by_id[r["gene_id"]]
                blk = _read_block(mm[r["blk_off"]:r["blk_off"] + r["blk_len"]], dof, r["n_var"], r["var_start"], f"{path.name} {ph}", kind=3)
                n = blk["n"]
                m = sph == ph
                if int(m.sum()) != n:
                    count_bad.append(ph)
                    continue
                if g["gene_id"] not in ranges:
                    ranges[g["gene_id"]] = _read_pages(ctx["vfile"][g["var_off"]:g["var_off"] + g["var_len"]], f"{ctx['var_path'].name} {g['gene_id']}")
                rng = ranges[g["gene_id"]]
                i0 = blk["var_start"] - rng["first"]
                if i0 < 0 or i0 + n > len(rng["position"]) or rng["position"][i0] != blk["pos_first"] or rng["position"][i0 + n - 1] != blk["pos_last"]:
                    raise PackError(f"{path.name} {ph}: the gene's variants range does not hold the intron run at pos_first/pos_last")
                sl = slice(i0, i0 + n)
                dec["position"].append(rng["position"][sl])
                dA1 += rng["A1"][sl]
                dA2 += rng["A2"][sl]
                dec["rs"].append(rng["rs_number"][sl].astype(np.int64))
                dec["tss"].append(rng["position"][sl] - blk["anchor"])
                dec["af"].append(rng["af_code"][sl])
                dec["ms"].append(rng["ms"][sl])
                dec["mc"].append(rng["mc"][sl])
                for k in ("nlp", "p", "se", "slope", "negative", "nq", "sq"):
                    dec[k].append(blk[k])
                pip = np.full(n, np.nan, dtype=np.float32)
                cs_id = np.full(n, -1, dtype=np.int64)
                for row, value, cid in zip(blk["cs"]["row"], blk["cs"]["pip"], blk["cs"]["cs_id"]):
                    if np.isnan(pip[row]) or value > pip[row] or (value == pip[row] and cid < cs_id[row]):
                        pip[row], cs_id[row] = value, cid
                got_cs += [(ph, int(rng["position"][i0 + rr]), rng["A1"][i0 + rr], rng["A2"][i0 + rr], int(c), float(v))
                           for rr, v, c in zip(blk["cs"]["row"].tolist(), blk["cs"]["pip"], blk["cs"]["cs_id"].tolist())]
                dec["pip"].append(pip)
                dec["cs_id"].append(cs_id)
                dec["nlp_max"].append(np.full(n, blk["nlp_max"]))
                se_b, sl_b = packfmt.error_bounds(S["pval_nominal"][m], S["slope"][m], S["slope_se"][m], blk["nlp_max"], blk["lse_min"], blk["lse_max"], dof)
                dec["se_b"].append(se_b)
                dec["sl_b"].append(sl_b)
            ranges.clear()
            want_cs_rows = sorted((a, int(b), c, e, int(f), float(g)) for a, b, c, e, f, g in con.execute(
                f"""SELECT phenotype_id, position, A1, A2, cs_id, pip FROM '{cfg.tables / 'credible_sets.parquet'}'
                WHERE qtl_type = 's' AND chr = ? AND phenotype_id IN ({in_list})""", [chrom]).fetchall())
            check(sorted(got_cs) == want_cs_rows, f"{chrom}: sQTL credible-set (position, A1, A2, cs_id, pip) of the {len(sample)} sample introns equal "
                                                  f"credible_sets.parquet ({len(got_cs)} decoded, {len(want_cs_rows)} in the table)")
            check(not count_bad, f"{chrom}: sQTL round-trip sample of {len(sample)} introns: block row counts equal raw row counts ({count_bad[:5]})")
            D = {k: np.concatenate(v) if v else np.array([]) for k, v in dec.items()}
            keep = ~np.isin(sph, count_bad)
            S = {k: (np.asarray(v, dtype=object)[keep].tolist() if k in ("A1", "A2") else v[keep]) for k, v in S.items()}
            sph = sph[keep]
            why = ", ".join(f"{sum(1 for x in reasons.values() if x == k)} {k}" for k in ("gene", "largest", "two sets", "random"))
            worst = _compare_rows(check, f"{chrom} sQTL", f"{len(sph):,} rows of {len(sample)} introns ({why})", D, dA1, dA2, S, dof, tol, af_exact=True)
            if worst is None:
                return
            worst["chrom"] = chrom
            worst["largest"] = max((by_ph[ph]["n_var"], ph) for ph in sample)
            summary["sqtl_rt"].append(worst)

            # ---- reference files for the browser decoder check ----
            if chrom in SQTL_REF_CHROMS:
                ref_ids = [ph for ph in sample if reasons[ph] in ("gene", "two sets")] if chrom != "chrX" else \
                    [next(ph for ph in sample if reasons[ph] == "random")]
                R = summary["sqtl_ref"]
                R["files"][chrom] = {"sqtl": str(path.relative_to(d)), "variants": str(ctx["var_path"].relative_to(d))}
                gids = sorted({by_ph[ph]["gene_id"] for ph in ref_ids})
                R["genes"] += [by_id[gid] for gid in gids]
                R["introns"] += [{"chrom": chrom, "gene_id": by_ph[ph]["gene_id"], "phenotype_id": ph, "blk_off": by_ph[ph]["blk_off"],
                                  "blk_len": by_ph[ph]["blk_len"]} for ph in ref_ids]
                mask = np.isin(sph, ref_ids)
                rows = src.filter(pa.array(keep))
                rows = rows.append_column("gene_id", pa.array([by_ph[ph]["gene_id"] for ph in sph.tolist()], type=pa.string()))
                R["rows"].append(rows.filter(pa.array(mask)))
                R["bounds"].append(pa.table({"phenotype_id": pa.array(sph[mask].tolist(), type=pa.string()),
                                             "position": S["position"][mask].astype(np.int32),
                                             "A1": pa.array([a for a, k in zip(S["A1"], mask) if k], type=pa.string()),
                                             "A2": pa.array([a for a, k in zip(S["A2"], mask) if k], type=pa.string()),
                                             "se_bound": D["se_b"][mask], "slope_bound": D["sl_b"][mask]}))
        except PackError as e:
            check(False, f"{chrom} sQTL pack decodes under SPEC: {e}")


def _validate_sqtl_summary(cfg: Config, check, summary: dict, out_dir: Path) -> None:
    T = summary["sqtl_totals"]
    check(T["blocks"] == 80_750 and T["rows"] == 499_599_113 and T["records"] == 246_322,
          f"sQTL packs hold {T['blocks']:,} blocks, {T['rows']:,} rows, and {T['records']:,} credible-set records (expected 80,750, 499,599,113, 246,322)")
    rt = summary["sqtl_rt"]
    if rt:
        def worst(k):
            w = max(rt, key=lambda x: x[k])
            return f"{w[k]:.3g} ({w['chrom']})"
        log(f"sQTL round trip over {sum(x['rows'] for x in rt):,} sample rows on {len(rt)} chromosomes: worst SE error {worst('se')} "
            f"(SPEC A.2 over every row: 1.93e-05), slope error {worst('slope')} (1.45e-03), -log10 p error {worst('nlp')} (1.56e-03), af error {worst('af')}")
        log(f"sQTL round trip: largest error / bound {max(x['r_se'] for x in rt):.4f} SE, {max(x['r_sl'] for x in rt):.4f} slope "
            f"(SPEC A.2: 1.0000, 1.0008; limit {packfmt.BOUND_FACTOR}); -log10 p error / half step {max(x['nlp_ratio'] for x in rt):.4f}; "
            f"slope from unquantized p and SE {max(x['formula_rel'] for x in rt):.3g}; largest sample intron {max(x['largest'] for x in rt)}")
    R = summary["sqtl_ref"]
    if R["rows"]:
        sq = out_dir / "sqtl"
        sq.mkdir(parents=True, exist_ok=True)
        pk = cfg["packs"]
        (sq / "introns.json").write_text(json.dumps({"dof": int(pk["dof"]["sqtl"]), "variant_page_size": pk["variant_page_size"],
                                                     "files": R["files"], "genes": R["genes"], "introns": R["introns"]}, indent=1))
        _write_arrow(pa.concat_tables(R["rows"]), sq / "rows.arrow")
        _write_arrow(pa.concat_tables(R["bounds"]), sq / "bounds.arrow")
        log(f"pack validate: sQTL reference files in {sq.relative_to(cfg.derived)}/: {len(R['introns'])} introns of {len(R['genes'])} genes "
            f"on {', '.join(R['files'])}, {sum(t.num_rows for t in R['rows']):,} rows")


# ---- validate: GWAS pack and index (SPEC section 11) ----------------------------------------------
GWAS_EXPECT_ROWS = 12_504_079
GWAS_WINDOW_GENES = {"ENSG00000128591": "FLNC", "ENSG00000177791": "MYOZ1", "ENSG00000145349": "CAMK2D"}
GWAS_REF_CHROMS = ("chr1", "chr22")      # reference windows for the browser decoder check
_GWAS_COLS = (("position_delta", "<u4"), ("beta", "<i4"), ("rs_number", "<u4"), ("se", "<u2"), ("eaf", "<u2"),
              ("p_mant", "<u2"), ("p_exp", "i1"), ("n_code", "u1"), ("allele", "u1"))
_P10 = np.array([float(f"1e{k}") for k in range(129)])


def _gwas_index(buf: bytes) -> tuple[int, np.ndarray, dict]:
    """gwas_index.bin under SPEC section 11: (rows per block, n table, {chrom: (first_position, end_offset)})."""
    count, block_rows, _ = _file_header(buf[:32], 5, "all", "gwas_index.bin")
    p = _unframe(buf[32:], None, "gwas_index.bin frame")
    off = 0

    def take(k: int) -> np.ndarray:
        nonlocal off
        if off + 4 * k > len(p):
            raise PackError("gwas_index.bin: payload ends early")
        a = np.frombuffer(p, "<u4", k, off).astype(np.int64)
        off += 4 * k
        return a

    nv = take(int(take(1)[0]))
    n_chroms = int(take(1)[0])
    if n_chroms != count or not 1 <= block_rows <= 65535 or not 1 <= nv.size <= 255 or np.any(np.diff(nv) <= 0):
        raise PackError(f"gwas_index.bin: {n_chroms} chromosomes (header {count}), page size {block_rows}, {nv.size} n values")
    chroms = {}
    for _ in range(n_chroms):
        if off + 8 > len(p):
            raise PackError("gwas_index.bin: payload ends early")
        name = p[off:off + 8].rstrip(b"\0").decode("ascii")
        off += 8
        nb = int(take(1)[0])
        fp, eo = take(nb), take(nb)
        if nb < 1 or name in chroms or fp[0] < 1 or np.any(np.diff(fp) < 0) or eo[0] <= 32 or np.any(np.diff(eo) <= 0):
            raise PackError(f"gwas_index.bin {name}: {nb} blocks; positions must not decrease and offsets must increase from past byte 32")
        chroms[name] = (fp, eo)
    if off != len(p):
        raise PackError(f"gwas_index.bin: {len(p) - off} bytes after the last chromosome")
    return block_rows, nv, chroms


def _gwas_block(frame: bytes, nv: np.ndarray, what: str) -> dict:
    """One GWAS block under SPEC section 11: its frame, payload length, codes, and heap; values as codes, p as float64."""
    p = _unframe(frame, None, what)
    if len(p) < 8:
        raise PackError(f"{what}: payload shorter than 8 bytes")
    n, heap_len = struct.unpack_from("<II", p, 0)
    if n < 1 or len(p) != 8 + 21 * n + heap_len:
        raise PackError(f"{what}: payload {len(p)} bytes for n {n} and heap_len {heap_len} (need 8 + 21n + heap_len)")
    c, off = {}, 8
    for name, t in _GWAS_COLS:
        c[name] = np.frombuffer(p, t, n, off)
        off += np.dtype(t).itemsize * n
    exp = c["p_exp"].astype(np.int64)
    mant = c["p_mant"].astype(np.int64)
    if (c["allele"] > 12).any() or (c["eaf"] > 10000).any() or ((mant < 1000) | (mant > 9999)).any() \
            or ((exp > -3) | ((exp == -3) & (mant != 1000))).any() or (c["n_code"] >= nv.size).any():
        raise PackError(f"{what}: an allele code above 12, eaf code above 10000, p mantissa outside 1000..9999, p above 1, or n code outside the table")
    codes = c["allele"].tolist()
    ea = [_SNP[x][0] if x else None for x in codes]
    nea = [_SNP[x][1] if x else None for x in codes]
    zero = [i for i, x in enumerate(codes) if x == 0]
    heap = p[off:]
    if not zero:
        if heap:
            raise PackError(f"{what}: heap of {len(heap)} bytes with no code-0 rows")
    else:
        if max(heap) >= 0x80 or not heap.endswith(b"\n"):
            raise PackError(f"{what}: heap is not ASCII records ending in a newline")
        recs = heap[:-1].decode("ascii").split("\n")
        if len(recs) != len(zero):
            raise PackError(f"{what}: {len(recs)} heap records for {len(zero)} code-0 rows")
        for i, rec in zip(zero, recs):
            f = rec.split("\t")
            if len(f) != 2 or tuple(f) in _SNP:
                raise PackError(f"{what}: heap record {rec!r} for row {i} needs one tab and must not be a coded SNP")
            ea[i], nea[i] = f
    pos = np.cumsum(c["position_delta"], dtype=np.int64)
    return {"n": n, "position": pos, "ea": ea, "nea": nea, "rs_number": c["rs_number"].astype(np.int64), "n_value": nv[c["n_code"]],
            "beta_q": c["beta"].astype(np.int64), "se_q": c["se"].astype(np.int64), "eaf_q": c["eaf"].astype(np.int64),
            "p": mant.astype(np.float64) / _P10[-exp]}


def _gwas_concat(parts: list[dict]) -> dict:
    out = {k: np.concatenate([x[k] for x in parts]) for k in ("position", "rs_number", "n_value", "beta_q", "se_q", "eaf_q", "p")}
    out["ea"] = [a for x in parts for a in x["ea"]]
    out["nea"] = [a for x in parts for a in x["nea"]]
    return out


def _gwas_rows_equal(D: dict, S: dict, i0: int, i1: int) -> tuple[bool, float]:
    """Decoded rows D against source rows S[i0:i1]: exact for position, alleles, rs_number, n and the 1e4 codes; p within 1e-12."""
    if len(D["position"]) != i1 - i0:
        return False, math.inf
    sp = S["p"][i0:i1]
    rel = float((np.abs(D["p"] - sp) / sp).max()) if i1 > i0 else 0.0
    ok = (np.array_equal(D["position"], S["position"][i0:i1]) and np.array_equal(D["rs_number"], S["rs_number"][i0:i1])
          and np.array_equal(D["n_value"], S["n"][i0:i1]) and D["ea"] == S["ea"][i0:i1] and D["nea"] == S["nea"][i0:i1]
          and all(np.array_equal(D[f"{k}_q"], S[f"{k}_q"][i0:i1]) for k in ("beta", "se", "eaf")) and rel <= 1e-12)
    return ok, rel


def _validate_gwas(cfg: Config, check, out_dir: Path) -> None:
    """SPEC section 11: every GWAS row decodes to the source TSV row (its own read of the TSV, not pack_gwas's scratch copy),
    and windows decode to exactly the source rows inside them. Writes reference files for the browser decoder check under pack_check/gwas/."""
    from .steps_gwas import COLUMNS, SOURCE_FILTER
    d, pk = cfg.derived, cfg["packs"]
    ipath = gwas_index_path(cfg)
    if not ipath.exists():
        check(False, f"GWAS index {ipath.relative_to(d)} exists (run `build --step pack_gwas`)")
        return
    tmpdir = cfg.tmp / "validate_gwas"
    shutil.rmtree(tmpdir, ignore_errors=True)
    tmpdir.mkdir(parents=True)
    spq = tmpdir / "source.parquet"
    gcon = connect(cfg, memory_limit=pk["gwas_duckdb_memory_limit"], threads=int(pk["gwas_duckdb_threads"]), temp_dir=tmpdir / "duckdb")
    gcon.execute(f"""COPY (
        SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr, CAST(split_part(CHRBP_B38, ':', 2) AS BIGINT) AS position,
               EA AS ea, NEA AS nea, BETA AS beta, SE AS se, EAFREQ AS eaf, P AS p,
               CASE WHEN regexp_matches(rsID, '^rs[0-9]+$') THEN CAST(substr(rsID, 3) AS BIGINT) ELSE 0 END AS rs_number, N AS n
        FROM read_csv('{cfg.raw / cfg["dcm_gwas"]}', delim='\\t', header=true, columns={COLUMNS}) WHERE {SOURCE_FILTER}
    ) TO '{spq}' (FORMAT parquet)""")
    per_chr = dict(gcon.execute(f"SELECT chr, count(*) FROM '{spq}' GROUP BY 1").fetchall())
    try:
        block_rows, nv, index = _gwas_index(ipath.read_bytes())
    except PackError as e:
        check(False, f"GWAS index decodes under SPEC: {e}")
        return
    with_rows = [c for c in CHROMS if per_chr.get(c)]
    check(list(index) == with_rows and block_rows == int(pk["gwas_block_rows"]) and set(per_chr) <= set(CHROMS),
          f"gwas_index.bin: kind 5, {len(index)} chromosomes in order ({', '.join(index)}) equal the source's chromosomes with rows; "
          f"{block_rows} rows per block; n table {nv.tolist()}; {ipath.stat().st_size:,} B")
    check("chrX" not in per_chr and "chrX" not in index and not gwas_path(cfg, "chrX").exists(), "no chrX GWAS rows, index entry, or file")
    genes = {r["gene_id"]: r for r in read_search_index(cfg, ["gene_id", "chr", "w_lo", "w_hi"]).to_pylist()
             if r["gene_id"] in GWAS_WINDOW_GENES}
    rng = np.random.default_rng(20260915)
    seeded = {c: 0 for c in index}
    for c in rng.choice(list(index), 50):
        seeded[str(c)] += 1
    total, worst_p, results = 0, 0.0, []
    ref = {"windows": [], "rows": []}
    for chrom in index:
        path = gwas_path(cfg, chrom)
        fp, eo = index[chrom]
        t = gcon.execute(f"SELECT position, ea, nea, beta, se, eaf, p, rs_number, n FROM '{spq}' WHERE chr = ? ORDER BY position, ea, nea, rs_number, p",
                         [chrom]).fetch_arrow_table()
        S = {k: t.column(k).to_numpy() for k in ("position", "p", "rs_number", "n")}
        S["ea"], S["nea"] = t.column("ea").to_pylist(), t.column("nea").to_pylist()
        for k in ("beta", "se", "eaf"):
            S[f"{k}_q"] = np.rint(t.column(k).to_numpy() * 1e4).astype(np.int64)
        try:
            with open(path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                count, page, _ = _file_header(mm[:32], 4, chrom, path.name)
                if eo[-1] != len(mm) or page != block_rows or count != (fp.size - 1) * block_rows + (count - (fp.size - 1) * block_rows) \
                        or not 1 <= count - (fp.size - 1) * block_rows <= block_rows:
                    raise PackError(f"{path.name}: {count:,} rows, page size {page}, {fp.size} index blocks ending at {eo[-1]:,}, file {len(mm):,} B")
                starts = np.r_[32, eo[:-1]]
                parts, last = [], 0
                for k in range(fp.size):
                    b = _gwas_block(mm[starts[k]:eo[k]], nv, f"{path.name} block {k}")
                    want_n = block_rows if k < fp.size - 1 else count - (fp.size - 1) * block_rows
                    if b["n"] != want_n or b["position"][0] != fp[k] or b["position"][0] < last or np.any(np.diff(b["position"]) < 0):
                        raise PackError(f"{path.name} block {k}: {b['n']} rows (want {want_n}), first position {b['position'][0]} "
                                        f"(index {fp[k]}), previous block ends at {last}")
                    last = int(b["position"][-1])
                    parts.append(b)
                D = _gwas_concat(parts)
                ok, rel = _gwas_rows_equal(D, S, 0, t.num_rows)
                worst_p = max(worst_p, rel)
                check(ok and count == t.num_rows, f"{chrom}: GWAS pack is kind 4 with {count:,} rows in {fp.size:,} blocks of {block_rows} "
                      f"back to back from byte 32 to the file end ({len(mm):,} B), each passing the payload, code, and heap rules and starting at "
                      f"its index first_position; the rows equal the source's {t.num_rows:,} in order (position, ea, nea, rs_number, n exact; "
                      f"rint(beta, se, eaf x 1e4) equal; p relative error {rel:.2e})")
                total += count
                del parts, D

                # windows: seeded, a position repeated across a block boundary, first and last rows, empty, gene windows
                pos = S["position"]
                wins = []
                for _ in range(seeded[chrom]):
                    lo = int(rng.integers(max(1, int(pos[0]) - 1_000_000), int(pos[-1]) + 1_000_000))
                    wins.append(("seeded", lo, lo + int(rng.integers(0, 3_000_000))))
                dup = next((k for k in range(1, fp.size) if pos[k * block_rows - 1] == fp[k]), None)
                if dup is not None:
                    wins.append(("boundary duplicate", int(fp[dup]), int(fp[dup]) + 20_000))
                wins += [("first row", int(pos[0]), int(pos[0]) + 50_000), ("last row", int(pos[-1]) - 50_000, int(pos[-1]))]
                gap = int(np.argmax(np.diff(pos)))
                wins += [("empty, before the first row", 1, int(pos[0]) - 1), ("empty, inside the largest gap", int(pos[gap]) + 1, int(pos[gap + 1]) - 1)]
                wins += [(f"{GWAS_WINDOW_GENES[g]} [w_lo, w_hi]", int(r["w_lo"]), int(r["w_hi"])) for g, r in genes.items() if r["chr"] == chrom]
                for name, lo, hi in wins:
                    i0, i1 = int(np.searchsorted(pos, lo, "left")), int(np.searchsorted(pos, hi, "right"))
                    end = int(np.searchsorted(fp, hi, "right")) - 1
                    if end < 0:
                        rng_b, ok, nb = None, i1 == i0, 0
                    else:
                        start = max(int(np.searchsorted(fp, lo, "left")) - 1, 0)
                        b0, b1 = (32 if start == 0 else int(eo[start - 1])), int(eo[end])
                        rng_b, nb = (b0, b1), end - start + 1
                        chunk = mm[b0:b1]
                        wparts = [_gwas_block(chunk[int(starts[k]) - b0:int(eo[k]) - b0], nv, f"{path.name} window block {k}") for k in range(start, end + 1)]
                        W = _gwas_concat(wparts)
                        m = (W["position"] >= lo) & (W["position"] <= hi)
                        sel = np.flatnonzero(m)
                        F = {k: (v[sel] if isinstance(v, np.ndarray) else [v[i] for i in sel]) for k, v in W.items()}
                        ok = _gwas_rows_equal(F, S, i0, i1)[0]
                    results.append({"chrom": chrom, "name": name, "lo": lo, "hi": hi, "range": rng_b, "blocks": nb, "rows": i1 - i0, "ok": ok})
                    if chrom in GWAS_REF_CHROMS and ok:
                        wid = len(ref["windows"])
                        ref["windows"].append({"id": wid, "name": name, "chrom": chrom, "lo": lo, "hi": hi, "rows": i1 - i0,
                                               "byte_start": None if rng_b is None else rng_b[0], "byte_end": None if rng_b is None else rng_b[1]})
                        ref["rows"].append(t.slice(i0, i1 - i0).append_column("window", pa.array(np.full(i1 - i0, wid, dtype=np.int32))))
        except PackError as e:
            check(False, f"{chrom} GWAS pack decodes under SPEC: {e}")
    gcon.close()
    check(total == GWAS_EXPECT_ROWS == sum(per_chr.values()), f"GWAS packs hold {total:,} rows (source {sum(per_chr.values()):,}, expected {GWAS_EXPECT_ROWS:,}); "
          f"largest p relative error {worst_p:.2e}")
    for label, pick in (("50 seeded windows", lambda r: r["name"] == "seeded"),
                        ("windows starting at a position repeated across a block boundary", lambda r: r["name"] == "boundary duplicate"),
                        ("windows touching each chromosome's first and last rows", lambda r: r["name"] in ("first row", "last row")),
                        ("empty windows", lambda r: r["name"].startswith("empty")),
                        ("gene [w_lo, w_hi] windows", lambda r: "w_lo" in r["name"])):
        rs = [r for r in results if pick(r)]
        bad = [f"{r['chrom']}:{r['lo']}-{r['hi']}" for r in rs if not r["ok"]]
        eg = "; ".join(f"{r['name']} {r['chrom']}:{r['lo']:,}-{r['hi']:,} {r['rows']:,} rows, "
                       + ("no range" if r["range"] is None else f"{r['blocks']} blocks, {r['range'][1] - r['range'][0]:,} B") for r in rs[:3])
        check(bool(rs) and not bad, f"GWAS {label}: {len(rs)} decode, through the window rule, to exactly the source rows with lo <= position <= hi"
              + (f" (failed: {bad[:5]})" if bad else "") + f"; e.g. {eg}")
    empty = [r for r in results if r["name"].startswith("empty")]
    check(all(r["rows"] == 0 for r in empty) and any(r["range"] is None for r in empty) and any(r["range"] is not None for r in empty),
          "GWAS empty windows hold no source rows; one has no byte range and one has a range whose rows all fall outside")
    if ref["windows"]:
        gd = out_dir / "gwas"
        gd.mkdir(parents=True, exist_ok=True)
        (gd / "windows.json").write_text(json.dumps({
            "index": str(ipath.relative_to(d)), "block_rows": block_rows, "n_values": nv.tolist(),
            "files": {c: str(gwas_path(cfg, c).relative_to(d)) for c in GWAS_REF_CHROMS},
            "note": "byte_start is the first byte and byte_end one past the last byte of the window's blocks; null when the range is empty. "
                    "rows.arrow holds the source rows with lo <= position <= hi, in file order, tagged by window id.",
            "windows": ref["windows"]}, indent=1))
        rows = pa.concat_tables(ref["rows"])
        rows = rows.select(["window", "position", "ea", "nea", "rs_number", "beta", "se", "eaf", "p", "n"]).cast(pa.schema([
            ("window", pa.int32()), ("position", pa.int32()), ("ea", pa.string()), ("nea", pa.string()), ("rs_number", pa.int64()),
            ("beta", pa.float64()), ("se", pa.float64()), ("eaf", pa.float64()), ("p", pa.float64()), ("n", pa.int32())]))
        _write_arrow(rows, gd / "rows.arrow")
        log(f"pack validate: GWAS reference files in {gd.relative_to(d)}/: {len(ref['windows'])} windows on {', '.join(GWAS_REF_CHROMS)}, {rows.num_rows:,} rows")
    shutil.rmtree(tmpdir)

