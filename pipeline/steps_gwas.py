"""DCM GWAS (Jurgens et al. 2024): the landing-track bins and the GWAS pack.

`gwas_bins`: one row per fixed genomic window, with the strongest variant (min p, its position and
rsID), how many variants in the window pass genome-wide significance, and how many were tested.

`pack_gwas`: every row, lossless at the source's precision, in the binary format of SPEC.md section 11:
one content-addressed `gwas/<chr>` file (chr1-22; the source has no chrX rows) and `gwas_index`, the
small index the gene page loads at startup to turn its `[w_lo, w_hi]` window into one byte range.
Both publish under `immutable/` (SPEC section 3). It also writes
`gwas_dcm.json` (file, cases, controls, variants) for the manifest, and `_tmp/pack_pointers/gwas.json`
with counts, sizes, and the byte range of every gene's window.

Both read straight from the raw meta-analysis TSV; GRCh38 coordinates come from its CHRBP_B38 column.
"""
import json
import os
import shutil
import time

import numpy as np
import pyarrow.parquet as pq

from . import packfmt
from .common import CHROMS, Config, connect, log, publish_file, read_search_index, stage, write_parquet

COLUMNS = ("{'CHRBP_B37':'VARCHAR','CHRBP_B38':'VARCHAR','ID_B38':'VARCHAR','CHR':'VARCHAR','POS':'BIGINT','EA':'VARCHAR',"
           "'NEA':'VARCHAR','BETA':'DOUBLE','SE':'DOUBLE','P':'DOUBLE','EAFREQ':'DOUBLE','HetDf':'INTEGER','HetPVal':'DOUBLE',"
           "'N':'BIGINT','N_cases':'BIGINT','N_controls':'BIGINT','rsID':'VARCHAR','ID_B37':'VARCHAR','INDEL':'VARCHAR'}")
# the rows the browser shows: a GRCh38 position and p above 0 (pack_gwas and validate)
SOURCE_FILTER = "CHRBP_B38 IS NOT NULL AND CHRBP_B38 <> '' AND P > 0"
# row order within a chromosome; rs_number and p only break ties, so the order is deterministic (the file has no two rows
# with one position and allele pair)
PACK_ORDER = "position, ea, nea, rs_number, p"


def pack(cfg: Config) -> None:
    from .steps_pack import EXT, gwas_index_path, gwas_path, pointer_dir
    pk = cfg["packs"]
    src = cfg.raw / cfg["dcm_gwas"]
    if not src.exists():
        raise FileNotFoundError(f"pack_gwas: {src} is missing (data/raw/download.py --only dcm_gwas, then unzip it)")
    block_rows, level = int(pk["gwas_block_rows"]), int(pk["zstd_level"])
    t0 = time.time()
    scratch = cfg.tmp / "pack_gwas"
    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True)
    gpq = scratch / "g.parquet"
    con = connect(cfg, memory_limit=pk["gwas_duckdb_memory_limit"], threads=int(pk["gwas_duckdb_threads"]), temp_dir=scratch / "duckdb")
    # the one read of the TSV; rs_number is the digits of an `rs` rsID (null when they do not parse, which fails below), else 0
    con.execute(f"""COPY (
        SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr, TRY_CAST(split_part(CHRBP_B38, ':', 2) AS BIGINT) AS position,
               EA AS ea, NEA AS nea, BETA AS beta, SE AS se, EAFREQ AS eaf, P AS p,
               CASE WHEN rsID LIKE 'rs%' THEN nullif(TRY_CAST(substr(rsID, 3) AS BIGINT), 0) ELSE 0 END AS rs_number,
               N AS n, N_cases AS n_cases, N_controls AS n_controls
        FROM read_csv('{src}', delim='\\t', header=true, columns={COLUMNS})
        WHERE {SOURCE_FILTER}
        ORDER BY chr, {PACK_ORDER}
    ) TO '{gpq}' (FORMAT parquet, ROW_GROUP_SIZE 1000000)""")
    log(f"pack_gwas: read {src.name} to {gpq.relative_to(cfg.derived)} in {time.time() - t0:.0f} s")
    g = f"'{gpq}'"
    nulls = con.execute(f"""SELECT count(*) FILTER (WHERE position IS NULL), count(*) FILTER (WHERE rs_number IS NULL),
        count(*) FILTER (WHERE n IS NULL), count(*) FILTER (WHERE ea IS NULL OR nea IS NULL) FROM {g}""").fetchone()
    if any(nulls):
        raise RuntimeError(f"pack_gwas: rows with an unparsed position, an unparsed rsID, a null n, or a null allele: {nulls}")
    n_cases, n_controls = con.execute(f"SELECT max(n_cases), max(n_controls) FROM {g}").fetchone()
    per_chr = dict(con.execute(f"SELECT chr, count(*) FROM {g} GROUP BY 1").fetchall())
    stray = {c: k for c, k in sorted(per_chr.items()) if c not in CHROMS}
    if stray:
        raise RuntimeError(f"pack_gwas: rows on chromosomes outside {CHROMS[0]}..{CHROMS[-1]}: {stray}")
    n_values = [int(r[0]) for r in con.execute(f"SELECT DISTINCT n FROM {g} ORDER BY 1").fetchall()]

    index_chroms, stats = [], {"block_rows": block_rows, "n_values": n_values, "rows_by_chrom": {}, "blocks_by_chrom": {}, "bytes_by_chrom": {}}
    heap_rows = 0
    for chrom in CHROMS:
        key, path = f"gwas/{chrom}", gwas_path(cfg, chrom)
        t = con.execute(f"SELECT position, ea, nea, beta, se, eaf, p, rs_number, n FROM {g} WHERE chr = ? ORDER BY {PACK_ORDER}", [chrom]).fetch_arrow_table()
        if t.num_rows == 0:
            path.unlink(missing_ok=True)      # a published file from an earlier build; the placeholder never exists
            continue
        col = {c: t.column(c).to_numpy() for c in ("position", "beta", "se", "eaf", "p", "rs_number", "n")}
        try:
            codes = packfmt.gwas_codes(col["position"], col["beta"], col["se"], col["eaf"], col["p"], col["rs_number"], col["n"], n_values)
        except ValueError as e:
            raise RuntimeError(f"pack_gwas {chrom}: {e}") from None
        ea, nea = t.column("ea").to_pylist(), t.column("nea").to_pylist()
        heap_rows += sum(1 for a, b in zip(ea, nea) if (a, b) not in packfmt.SNP_CODES)
        tmp = stage(cfg, key, EXT["gwas"])
        first_position, end_offset = [], []
        with open(tmp, "wb") as fh:
            fh.write(packfmt.file_header(packfmt.KIND_GWAS, chrom, t.num_rows, block_rows))
            off = packfmt.FILE_HEADER_LEN
            for s in range(0, t.num_rows, block_rows):
                e = min(s + block_rows, t.num_rows)
                try:
                    frame = packfmt.encode_gwas_block({k: v[s:e] for k, v in codes.items()}, ea[s:e], nea[s:e], level)
                except ValueError as err:
                    raise RuntimeError(f"pack_gwas {chrom}: block of rows {s}..{e - 1}: {err}") from None
                fh.write(frame)
                off += len(frame)
                first_position.append(int(codes["position"][s]))
                end_offset.append(off)
        if tmp.stat().st_size != off or off > packfmt.U32_MAX:
            raise RuntimeError(f"pack_gwas {chrom}: file is {tmp.stat().st_size:,} bytes, blocks end at {off:,} (limit 4 GiB)")
        publish_file(cfg, tmp, key, EXT["gwas"])
        index_chroms.append((chrom, np.array(first_position), np.array(end_offset)))
        stats["rows_by_chrom"][chrom], stats["blocks_by_chrom"][chrom], stats["bytes_by_chrom"][chrom] = t.num_rows, len(end_offset), off
        log(f"pack_gwas {chrom}: {t.num_rows:,} rows in {len(end_offset):,} blocks, {off:,} B")
    idx = packfmt.encode_gwas_index(n_values, index_chroms, block_rows, level)
    tmp = stage(cfg, "gwas_index", EXT["gwas_index"])
    tmp.write_bytes(idx)
    ipath = publish_file(cfg, tmp, "gwas_index", EXT["gwas_index"])

    # the gene page's third request: the byte range of every gene's [w_lo, w_hi] window
    fp_eo = {c: (fp, eo) for c, fp, eo in index_chroms}
    si = read_search_index(cfg, ["symbol", "chr", "w_lo", "w_hi"]).to_pylist()
    lens, no_file = [], 0
    for r in si:
        if r["w_lo"] is None:
            continue
        if r["chr"] not in fp_eo:
            no_file += 1
            continue
        w = packfmt.gwas_window(*fp_eo[r["chr"]], r["w_lo"], r["w_hi"])
        lens.append((0 if w is None else w[3] - w[2], r["symbol"], r["chr"]))
    L = np.array([x[0] for x in lens])
    top = max(lens)
    stats.update({
        "rows": sum(stats["rows_by_chrom"].values()), "blocks": sum(stats["blocks_by_chrom"].values()),
        "bytes": sum(stats["bytes_by_chrom"].values()), "index_bytes": len(idx), "heap_rows": heap_rows,
        "window_bytes": {"genes": len(lens), "genes_without_gwas_file": no_file, "mean": float(L.mean()), "median": float(np.median(L)),
                         "p90": float(np.percentile(L, 90)), "max": int(top[0]), "max_gene": f"{top[1]} ({top[2]})"},
        "seconds": round(time.time() - t0, 1),
    })
    pointer_dir(cfg).mkdir(parents=True, exist_ok=True)
    (pointer_dir(cfg) / "gwas.json").write_text(json.dumps(stats, indent=1))
    # which of the Jurgens sets this is, for the manifest and the About page
    (cfg.derived / "gwas_dcm.json").write_text(json.dumps({
        "file": src.name, "n_cases": n_cases, "n_controls": n_controls, "variants": stats["rows"],
    }, indent=2))
    shutil.rmtree(scratch)
    W = stats["window_bytes"]
    log(f"pack_gwas: {stats['rows']:,} rows (expected 12,504,079) in {stats['blocks']:,} blocks of {block_rows}, {stats['bytes']:,} B in "
        f"{len(index_chroms)} files; {heap_rows:,} rows with heap alleles; n values {n_values}; {n_cases:,} cases / {n_controls:,} controls")
    log(f"pack_gwas: index {len(idx):,} B" + (" (over 100 KB: double packs.gwas_block_rows and rebuild)" if len(idx) > 100_000 else ""))
    log(f"pack_gwas: window byte range over {W['genes']:,} genes' [w_lo, w_hi] ({no_file} genes on chromosomes without a GWAS file): "
        f"mean {W['mean']:,.0f} B, median {W['median']:,.0f}, 90th percentile {W['p90']:,.0f}, max {W['max']:,} ({W['max_gene']}); "
        f"{stats['seconds']:.0f} s")


def run(cfg: Config) -> None:
    src = cfg.raw / cfg["dcm_gwas"]
    bin_bp = int(cfg["gwas_bin_bp"])
    con = connect(cfg)
    t = con.execute(f"""
        WITH v AS (
            SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr,
                   split_part(CHRBP_B38, ':', 2)::BIGINT AS position,
                   P AS p, rsID AS rsid, EA, NEA, BETA AS beta
            FROM read_csv('{src}', delim='\\t', header=true, columns={COLUMNS})
            WHERE CHRBP_B38 IS NOT NULL AND CHRBP_B38 <> '' AND P > 0
        )
        SELECT chr, (position // {bin_bp}) * {bin_bp} AS bin_start, (position // {bin_bp}) * {bin_bp} + {bin_bp} AS bin_end,
               min(p) AS min_p, arg_min(position, p) AS lead_position, arg_min(rsid, p) AS lead_rsid,
               arg_min(beta, p) AS lead_beta, arg_min(EA, p) AS lead_ea,
               count(*) FILTER (WHERE p < 5e-8)::INTEGER AS n_gws, count(*)::INTEGER AS n_variants
        FROM v GROUP BY 1, 2, 3 ORDER BY 1, 2
    """).fetch_arrow_table()
    # JSON only: the landing track fetches this with one plain request, before the engine has
    # booted, and nothing ever read the parquet copy. Columnar (one array per column) because
    # r2.dev serves JSON uncompressed and row objects repeat every key 569 times.
    cols = {name: t.column(name).to_pylist() for name in t.column_names}
    cols["min_p"] = [float(f"{p:.3g}") for p in cols["min_p"]]
    cols["lead_beta"] = [round(b, 3) for b in cols["lead_beta"]]
    (cfg.derived / "gwas_dcm_bins.json").write_text(json.dumps({"n": t.num_rows, "columns": cols}, separators=(",", ":")))
    n_gws = sum(1 for x in t.column("n_gws").to_pylist() if x > 0)
    log(f"gwas_bins: {t.num_rows} windows of {bin_bp // 1_000_000} Mb, {n_gws} with genome-wide significant variants -> gwas_dcm_bins.json")
