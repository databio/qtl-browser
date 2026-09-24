"""Step 7: re-encode the cis nominal pair files, one row group per phenotype: per gene for
eQTL, per intron for sQTL. Every UI query filters on the phenotype, and DuckDB reads whole
row groups, so a per-gene sQTL group would pull every significant intron of the gene to show
one (CAMK2D: 1.4 MB for nine introns instead of ~150 KB).

Memory-bounded: DuckDB writes the joined, sorted rows to a temporary parquet under its own
memory limit (spilling to disk if needed); pyarrow then streams that file batch by batch and
re-emits it with one row group per phenotype. Nothing holds a whole chromosome in memory.

Encodings (measured on chr22, 2026-09-03): per-phenotype row groups lose the cross-phenotype
repetition that dictionary encoding exploits in the raw million-row groups, so positions use
DELTA_BINARY_PACKED, floats use BYTE_STREAM_SPLIT, and rsIDs are stored as delta-encoded
integer rs_number. Result: eQTL files ~10% smaller than raw. sQTL files grew ~15% (1.99 to
2.28 GB) when the groups went from per gene to per intron (2026-09-05), since the variant
columns (rs_number, af, ma_count, alleles) repeat per intron and are no longer compressed
across the gene block.

`sqtl_nominal: significant` in the config keeps sQTL rows only for introns flagged is_sqtl in
splice_phenotypes (84.5M of 499.6M rows, 1.8 GB instead of 10.8), which is what fits the R2
free tier. eQTL rows are always kept for every tested gene.
"""
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .common import CHROMS, Config, connect, log, variants_sql
from .steps_tables import SPLICE_PARSE

STATS = {"e": ["gene_id", "position"], "s": ["gene_id", "phenotype_id", "position"]}
FLOATS = ("af", "slope", "slope_se", "pval_nominal", "pip")
ENCODING = {
    "e": {"position": "DELTA_BINARY_PACKED", "tss_distance": "DELTA_BINARY_PACKED", "rs_number": "DELTA_BINARY_PACKED",
          **{c: "BYTE_STREAM_SPLIT" for c in FLOATS}},
    "s": {"position": "DELTA_BINARY_PACKED", "tss_distance": "DELTA_BINARY_PACKED", "rs_number": "DELTA_BINARY_PACKED",
          **{c: "BYTE_STREAM_SPLIT" for c in FLOATS}},
}


def _regroup(src: Path, out_dir: Path, group_col: str, stats_columns: list[str], encoding: dict[str, str]) -> tuple[int, int]:
    """Stream `src` (sorted by bin, then group_col) into `out_dir/bin=<bin>/data.parquet` files,
    one row group per group value. The `bin` column is consumed by the path and not written.
    Returns (row groups, files)."""
    pf = pq.ParquetFile(src)
    schema = pf.schema_arrow.remove(pf.schema_arrow.get_field_index("bin"))
    writer: pq.ParquetWriter | None = None
    buf: list[pa.Table] = []
    cur = None
    cur_bin = None
    n = files = 0

    def open_writer(b: int):
        nonlocal writer, files
        path = out_dir / f"bin={b}" / "data.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(path, schema, compression="zstd", compression_level=9,
                                  use_dictionary=[c for c in schema.names if c not in encoding],
                                  column_encoding=encoding, write_statistics=stats_columns)
        files += 1

    def flush():
        nonlocal buf, n
        if buf:
            t = pa.concat_tables(buf).drop_columns(["bin"])
            writer.write_table(t, row_group_size=t.num_rows)
            n += 1
            buf = []

    try:
        for batch in pf.iter_batches(batch_size=250_000):
            t = pa.Table.from_batches([batch])
            col = t[group_col].combine_chunks()
            if len(col) == 0:
                continue
            bins = t["bin"].combine_chunks()
            shifted = col.slice(1)
            idx = pc.indices_nonzero(pc.not_equal(col.slice(0, len(col) - 1), shifted)).to_pylist()
            starts = [0] + [i + 1 for i in idx]
            ends = starts[1:] + [t.num_rows]
            for s, e in zip(starts, ends):
                g = col[s].as_py()
                if g != cur:
                    flush()
                    cur = g
                    b = bins[s].as_py()
                    if b != cur_bin:
                        if writer is not None:
                            writer.close()
                        open_writer(b)
                        cur_bin = b
                buf.append(t.slice(s, e - s))
        flush()
    finally:
        if writer is not None:
            writer.close()
    return n, files


def _one(args) -> tuple[str, str, int, int]:
    qtl_type, chrom, src, out = args
    out = Path(out)
    cfg = Config()
    work_tmp = cfg.tmp / f"nominal-{qtl_type}-{chrom}"
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"], temp_dir=work_tmp)
    con.execute("SET preserve_insertion_order = true")
    genes = cfg.tables / "genes.parquet"
    vpos = variants_sql(cfg, chrom)
    susie = cfg.raw_dir("cis_eQTL_SuSiE" if qtl_type == "e" else "cis_sQTL_SuSiE") / (
        f"topchef_{chrom}_MaxPC70.SuSiE_summary.parquet" if qtl_type == "e" else f"topchefSplice_{chrom}_MaxPC25.SuSiE_summary.parquet")
    con.execute(f"CREATE TABLE v AS SELECT position, A1, A2, rs_number FROM {vpos}")
    # a variant can belong to two credible sets of one phenotype: keep the higher-PIP membership
    con.execute(f"""CREATE TABLE s AS
        SELECT phenotype_id, position, A1, A2, max(pip)::FLOAT AS pip, arg_max(cs_id, pip)::TINYINT AS cs_id
        FROM '{susie}' GROUP BY 1, 2, 3, 4""")
    # tested genes carry their partition bin (TSS rank / nominal_bin_genes) in genes.parquet
    con.execute(f"CREATE TABLE g AS SELECT gene_id, tss, bin FROM '{genes}' WHERE chr = '{chrom}' AND bin IS NOT NULL")
    if qtl_type == "e":
        con.execute(f"CREATE VIEW n AS SELECT *, phenotype_id AS gene_id FROM '{src}'")
        select_extra, order = "", "g.bin, g.tss, n.gene_id, n.position"
    elif cfg["sqtl_nominal"] == "significant":
        sp = cfg.tables / "splice_phenotypes.parquet"
        con.execute(f"""CREATE VIEW n AS SELECT r.*, {SPLICE_PARSE.replace('phenotype_id', 'r.phenotype_id')}
            FROM '{src}' r SEMI JOIN (SELECT phenotype_id FROM '{sp}' WHERE is_sqtl) k USING (phenotype_id)""")
        select_extra, order = "n.phenotype_id,", "g.bin, g.tss, n.gene_id, n.phenotype_id, n.position"
    else:
        con.execute(f"CREATE VIEW n AS SELECT *, {SPLICE_PARSE} FROM '{src}'")
        select_extra, order = "n.phenotype_id,", "g.bin, g.tss, n.gene_id, n.phenotype_id, n.position"
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "sorted.tmp.parquet"
    con.execute(f"""
        COPY (
            SELECT g.bin, {select_extra} n.gene_id, n.position::INTEGER AS position, n.A1, n.A2, v.rs_number,
                   n.start_distance::INTEGER AS tss_distance, n.af::FLOAT AS af,
                   n.ma_samples::SMALLINT AS ma_samples, n.ma_count::SMALLINT AS ma_count,
                   n.pval_nominal, n.slope::FLOAT AS slope, n.slope_se::FLOAT AS slope_se,
                   s.pip, s.cs_id
            FROM n
            JOIN g ON g.gene_id = n.gene_id
            LEFT JOIN v ON v.position = n.position AND v.A1 = n.A1 AND v.A2 = n.A2
            LEFT JOIN s ON s.phenotype_id = n.phenotype_id AND s.position = n.position AND s.A1 = n.A1 AND s.A2 = n.A2
            ORDER BY {order}
        ) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000)
    """)
    raw_n = con.execute("SELECT count(*) FROM n").fetchone()[0]
    tmp_n = pq.read_metadata(tmp).num_rows
    if tmp_n != raw_n:
        tmp.unlink()
        raise RuntimeError(f"{src}: row count changed {raw_n} -> {tmp_n} (join duplicated rows)")
    con.close()
    shutil.rmtree(work_tmp, ignore_errors=True)
    groups, files = _regroup(tmp, out, "gene_id" if qtl_type == "e" else "phenotype_id", STATS[qtl_type], ENCODING[qtl_type])
    tmp.unlink()
    return qtl_type, chrom, raw_n, groups, files


def run(cfg: Config, force: bool = False) -> None:
    jobs = []
    for qtl_type, src_dir, out_dir, pat in [
        ("e", "cis_eQTL_nominal", "cis_eqtl_nominal", "topchef_{c}_MaxPC70.cis_qtl_pairs.{c}.parquet"),
        ("s", "cis_sQTL_nominal", "cis_sqtl_nominal", "topchefSplice_{c}_MaxPC25.cis_qtl_pairs.{c}.parquet"),
    ]:
        for c in CHROMS:
            src = cfg.raw_dir(src_dir) / pat.format(c=c)
            # under _tables/, where pack_eqtl, manifest, and packcheck read it; nothing the browser
            # fetches lives here
            out = cfg.tables / out_dir / f"chr={c}"           # bin=<n>/data.parquet files go inside
            if not src.exists():
                log(f"nominal: missing raw file {src.name}, skipping")
                continue
            done = list(out.glob("bin=*/data.parquet"))
            if not force and done and min(p.stat().st_mtime for p in done) > src.stat().st_mtime:
                continue
            jobs.append((qtl_type, c, str(src), str(out)))
    jobs.sort(key=lambda j: -Path(j[2]).stat().st_size)  # biggest first
    log(f"nominal: {len(jobs)} chromosomes to build with {cfg['workers']} workers")
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as ex:
        futs = {ex.submit(_one, j): j for j in jobs}
        for f in as_completed(futs):
            qtl_type, c, n, groups, files = f.result()
            log(f"nominal: {qtl_type} {c}: {n:,} rows in {groups} row groups across {files} bin files")
