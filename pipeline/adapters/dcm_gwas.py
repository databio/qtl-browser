"""The DCM GWAS adapter: Jurgens et al. 2024 (biobanks-only meta-analysis) -> the contract `gwas` table.

Writes `gwas.parquet`, `gwas.json` and `gwas_bins.parquet` into an experiment's tables directory (default the TOPCHeF one,
`_tables/topchef/`: the experiment the GWAS is shown next to, CONTRACT.md "gwas"). The browser uses
it for the gene page's GWAS panel and LocusCompare.

    uv run python -m pipeline.adapters.dcm_gwas [--tables DIR] [--bins-only]

Rows are every source row with a GRCh38 position (`CHRBP_B38`) and p > 0, as the v0 GWAS pack took
them. Orientation is the store's one rule (`qtlstore.orient_to_ref`): the reference is read at every
row from the refgetstore; where the source's effect allele `EA` is the reference base(s), the alleles
swap, `beta` negates and `af` becomes `1 - EAFREQ`; where it is the other allele nothing changes; where
the reference reads neither, the row is dropped and counted by class. Negation and `1 - x` are exact
at the source's 4 decimals, so the GWAS object stays lossless.

`gwas_bins.parquet` is the landing track's summary: per `gwas_bin_bp` window (5 Mb) the strongest p, its
row, and the row counts, over the **source** rows (v0's `gwas_bins` step, unchanged), so it matches v0's
`gwas_dcm_bins.json` bin for bin. `--bins-only` writes just that table (one pass over the TSV).
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import qtlstore as qs
from ..common import CHROMS, Config, connect, log, write_parquet
from ..steps_refget import A1, A2, BOTH, NEITHER, UNCHECKED, classify, load_sequence, open_store, reference_json

GWAS_ID = "dcm_jurgens2024_biobanks"
TITLE = "Dilated cardiomyopathy GWAS, Jurgens et al. 2024 (biobanks-only meta-analysis)"
COLUMNS = ("{'CHRBP_B37':'VARCHAR','CHRBP_B38':'VARCHAR','ID_B38':'VARCHAR','CHR':'VARCHAR','POS':'BIGINT','EA':'VARCHAR',"
           "'NEA':'VARCHAR','BETA':'DOUBLE','SE':'DOUBLE','P':'DOUBLE','EAFREQ':'DOUBLE','HetDf':'INTEGER','HetPVal':'DOUBLE',"
           "'N':'BIGINT','N_cases':'BIGINT','N_controls':'BIGINT','rsID':'VARCHAR','ID_B37':'VARCHAR','INDEL':'VARCHAR'}")
SOURCE_FILTER = "CHRBP_B38 IS NOT NULL AND CHRBP_B38 <> '' AND P > 0"
CLASS_NAMES = {A1: "ea_is_ref", A2: "nea_is_ref", BOTH: "both", NEITHER: "neither", UNCHECKED: "unchecked"}


def run(cfg: Config, out: Path) -> dict:
    src = cfg.raw / cfg["dcm_gwas"]
    if not src.exists():
        raise FileNotFoundError(f"dcm_gwas: {src} is missing (data/raw/download.py --only dcm_gwas, then unzip it)")
    t0 = time.time()
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    con.execute(f"""CREATE TABLE g AS
        SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr, TRY_CAST(split_part(CHRBP_B38, ':', 2) AS BIGINT) AS pos,
               upper(EA) AS ea, upper(NEA) AS nea, BETA AS beta, SE AS se, EAFREQ AS eaf, P AS pvalue, N AS n,
               CASE WHEN rsID LIKE 'rs%' THEN nullif(TRY_CAST(substr(rsID, 3) AS BIGINT), 0) ELSE 0 END AS rs_number,
               N_cases AS n_cases, N_controls AS n_controls
        FROM read_csv('{src}', delim='\\t', header=true, columns={COLUMNS})
        WHERE {SOURCE_FILTER}""")
    n_src = con.execute("SELECT count(*) FROM g").fetchone()[0]
    bad = con.execute("""SELECT count(*) FILTER (WHERE pos IS NULL), count(*) FILTER (WHERE rs_number IS NULL),
        count(*) FILTER (WHERE ea IS NULL OR nea IS NULL) FROM g""").fetchone()
    if any(bad):
        raise ValueError(f"dcm_gwas: rows with an unparsed position, an unparsed rsID or a null allele: {bad}")
    stray = con.execute(f"SELECT chr, count(*) FROM g WHERE chr NOT IN ({', '.join(repr(c) for c in CHROMS)}) "
                        "GROUP BY 1").fetchall()
    n_cases, n_controls = con.execute("SELECT max(n_cases), max(n_controls) FROM g").fetchone()
    seqs = json.loads(reference_json(cfg).read_text())["sequences"]
    store = open_store(cfg)
    counts = {"as_is": 0, "swapped": 0, "dropped": 0}
    by_class: dict[str, int] = {}
    parts = []
    for chrom in CHROMS:
        t = con.execute("SELECT pos, ea, nea, beta, se, eaf, pvalue, n, rs_number FROM g WHERE chr = ? "
                        "ORDER BY pos, ea, nea", [chrom]).fetch_arrow_table()
        if t.num_rows == 0:
            continue
        seq = load_sequence(store, seqs[chrom]["digest"], seqs[chrom]["length"])
        pos = t["pos"].to_numpy().astype(np.int64)
        ea, nea = t["ea"].to_pylist(), t["nea"].to_pylist()
        cls = classify(seq, pos, ea, nea)
        del seq
        for k, v in zip(*np.unique(cls, return_counts=True)):
            by_class[CLASS_NAMES[int(k)]] = by_class.get(CLASS_NAMES[int(k)], 0) + int(v)
        ref_base = np.where(cls == A1, np.array(ea, dtype=object), np.where(cls == A2, np.array(nea, dtype=object), None))
        r = qs.orient_to_ref(ref_base, ea, nea, t["beta"].to_numpy(), t["eaf"].to_numpy())
        keep = r["keep"]
        for k in counts:
            counts[k] += r["counts"][k]
        # -beta is exact; 1 - af can carry a float residue (1 - 0.3333 = 0.6667000000000001), which the
        # builder's 4-decimal rule absorbs: the stored code is rint(af * 1e4), the same integer either way
        parts.append(pa.table({
            "chr": pa.array([chrom] * int(keep.sum()), pa.string()),
            "pos": pa.array(pos[keep], pa.int32()), "ref": pa.array(list(r["ref"]), pa.string()),
            "alt": pa.array(list(r["alt"]), pa.string()),
            "beta": pa.array(r["beta"], pa.float64()), "se": t["se"].filter(pa.array(keep)),
            "af": pa.array(r["af"], pa.float64()), "pvalue": t["pvalue"].filter(pa.array(keep)),
            "n": t["n"].filter(pa.array(keep)).cast(pa.int64()),
            "rs_number": t["rs_number"].filter(pa.array(keep)).cast(pa.int64())}))
        log(f"dcm_gwas {chrom}: {t.num_rows:,} rows, {r['counts']}")
    tab = pa.concat_tables(parts)
    # the source lists some indels twice, once per allele order, with different N and statistics; after
    # orientation both rows name one site. Two measurements, both kept (CONTRACT.md `gwas`), counted here.
    # Two source records are the same row written once per allele order with mirrored values; oriented,
    # they are identical, so one copy of each goes (counted)
    con.register("tab0", tab)
    n0 = tab.num_rows
    tab = con.execute("SELECT DISTINCT * FROM tab0 ORDER BY chr, pos, ref, alt, rs_number, pvalue, n, beta").arrow()
    if not isinstance(tab, pa.Table):
        tab = tab.read_all()
    con.unregister("tab0")
    identical = n0 - tab.num_rows
    con.register("tab", tab)
    sharing = con.execute("SELECT count(*) FROM (SELECT count(*) OVER (PARTITION BY chr, pos, ref, alt) k FROM tab) "
                          "WHERE k > 1").fetchone()[0]
    con.unregister("tab")
    out.mkdir(parents=True, exist_ok=True)
    write_parquet(tab, out / "gwas.parquet", 500_000, stats_columns=["chr", "pos"])
    meta = {"id": GWAS_ID, "title": TITLE,
            "source": {"file": src.name, "config_key": "dcm_gwas", "n_cases": n_cases, "n_controls": n_controls,
                       "rows_with_grch38_position_and_p_above_0": int(n_src),
                       "rows_on_other_chromosomes": {c: int(k) for c, k in stray}},
            "orientation": {**counts, "by_reference_class": by_class, "rule": "qtlstore.orient_to_ref"},
            "rows": tab.num_rows, "rows_sharing_a_site": int(sharing), "identical_rows_dropped": int(identical)}
    tmp = out / "gwas.json.tmp"
    tmp.write_text(json.dumps(meta, indent=2) + "\n")
    os.replace(tmp, out / "gwas.json")
    log(f"dcm_gwas: {tab.num_rows:,} of {n_src:,} rows -> {out / 'gwas.parquet'}; {counts}; {by_class}; "
        f"{identical} identical rows dropped, {sharing:,} rows share a site; "
        f"{time.time() - t0:.0f} s")
    return meta


BINS_SQL = """
    WITH v AS (
        SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr, split_part(CHRBP_B38, ':', 2)::BIGINT AS position,
               P AS p, rsID AS rsid, EA, BETA AS beta
        FROM read_csv('{src}', delim='\\t', header=true, columns={columns})
        WHERE {where}
    )
    SELECT chr, (position // {bin_bp}) * {bin_bp} AS bin_start, (position // {bin_bp}) * {bin_bp} + {bin_bp} AS bin_end,
           min(p) AS min_p, arg_min(position, p) AS lead_position, arg_min(rsid, p) AS lead_rsid,
           arg_min(beta, p) AS lead_beta, arg_min(EA, p) AS lead_ea,
           count(*) FILTER (WHERE p < 5e-8)::INTEGER AS n_gws, count(*)::INTEGER AS n_variants
    FROM v GROUP BY 1, 2, 3 ORDER BY 1, 2
"""


def bins(cfg: Config, out: Path, con=None) -> int:
    """`gwas_bins.parquet`: v0's landing-track summary (steps_gwas.run), from the source rows as given.

    Source rows, not the oriented table: v0 binned every row with a GRCh38 position and p > 0, the 9,346
    rows the reference does not read included, with the lead's `EA` and its beta as published. Values
    unrounded here; the builder (`gwas.build`) rounds them as v0 did. Returns the number of bins."""
    src = cfg.raw / cfg["dcm_gwas"]
    con = con or connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    t = con.execute(BINS_SQL.format(src=src, columns=COLUMNS, where=SOURCE_FILTER,
                                    bin_bp=int(cfg["gwas_bin_bp"]))).fetch_arrow_table()
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "gwas_bins.parquet.tmp"
    pq.write_table(t, tmp)
    os.replace(tmp, out / "gwas_bins.parquet")
    log(f"dcm_gwas: {t.num_rows} bins of {int(cfg['gwas_bin_bp']) // 1_000_000} Mb -> {out / 'gwas_bins.parquet'}")
    return t.num_rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tables", type=Path, help="experiment tables dir (default: the TOPCHeF one)")
    ap.add_argument("--bins-only", action="store_true", help="write only gwas_bins.parquet")
    a = ap.parse_args(argv)
    cfg = Config()
    from .topchef import tables
    out = a.tables or tables(cfg)
    if not a.bins_only:
        run(cfg, out)
    bins(cfg, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
