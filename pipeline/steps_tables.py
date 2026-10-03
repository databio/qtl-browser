"""The v0 small tables: genes, splice_phenotypes, credible_sets.

These are the v0 tables, not the contract tables; `adapters/topchef.py` writes those. What the two
share is where the rows come from and how a leafcutter phenotype id is spelled, so both read the
archive names and `SPLICE_PARSE` from the adapter. The TOPCHeF acceptance gate
(`adapters/verify_topchef.py`) counts against them.
"""
from .adapters import topchef
from .adapters.topchef import SPLICE_PARSE
from .common import Config, connect, log, write_parquet


def _setup(cfg: Config):
    con = connect(cfg)
    con.execute(f"CREATE VIEW ann AS SELECT * FROM '{cfg.tables / 'gene_annotation.parquet'}'")
    con.execute(f"CREATE VIEW vpos AS SELECT * FROM '{cfg.tables / 'variants.parquet'}'")
    return con


def _bh(con, table: str, pcol: str) -> None:
    """Add a Benjamini-Hochberg qval column to `table` based on `pcol`."""
    con.execute(f"ALTER TABLE {table} ADD COLUMN qval DOUBLE")
    con.execute(f"""
        WITH r AS (
            SELECT phenotype_id, {pcol} AS p,
                   row_number() OVER (ORDER BY {pcol}) AS rk, count(*) OVER () AS n
            FROM {table} WHERE {pcol} IS NOT NULL
        ), q AS (
            SELECT phenotype_id,
                   least(1.0, min(p * n / rk) OVER (ORDER BY rk DESC ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)) AS qval
            FROM r
        )
        UPDATE {table} SET qval = q.qval FROM q WHERE {table}.phenotype_id = q.phenotype_id
    """)


def permutation_tables(cfg: Config) -> None:
    con = _setup(cfg)
    sig_col, thr = cfg["sig_column"], cfg["sig_threshold"]

    con.execute(f"CREATE TABLE sperm AS SELECT *, {SPLICE_PARSE} FROM read_parquet('{topchef.archive_glob(cfg, 's', 'permutation')}')")

    # ---- genes ----
    con.execute(f"CREATE TABLE perm AS SELECT * FROM read_parquet('{topchef.archive_glob(cfg, 'e', 'permutation')}')")
    _bh(con, "perm", "pval_beta")
    con.execute(f"""
        CREATE TABLE ncs AS SELECT phenotype_id, count(DISTINCT cs_id) AS n_credible_sets
        FROM read_parquet('{topchef.archive_glob(cfg, 'e', 'susie')}') GROUP BY 1
    """)
    con.execute(f"""
        CREATE TABLE ntrans AS SELECT phenotype_id, count(*) AS n_trans_pairs
        FROM read_parquet('{topchef.archive_glob(cfg, 'e', 'trans')}') GROUP BY 1
    """)
    genes = con.execute(f"""
        SELECT a.gene_id, a.gene_id_version, a.symbol, a.chr, a.start, a.end, a.strand, a.tss, a.biotype,
               p.phenotype_id IS NOT NULL AS tested,
               p.num_var, p.position AS lead_position, p.A1 AS lead_A1, p.A2 AS lead_A2, v.rsid AS lead_rsid,
               p.af AS lead_af, p.start_distance AS lead_tss_distance,
               p.slope, p.slope_se, p.pval_nominal, p.pval_perm, p.pval_beta, p.qval,
               CASE WHEN p.phenotype_id IS NULL THEN NULL ELSE p.{sig_col} < {thr} END AS is_egene,
               coalesce(c.n_credible_sets, 0)::INTEGER AS n_credible_sets,
               coalesce(t.n_trans_pairs, 0)::INTEGER AS n_trans_pairs,
               -- partition file of the gene's nominal rows: TSS rank among
               -- genes tested for eQTL or sQTL on the chromosome, {cfg['nominal_bin_genes']} genes per bin
               CASE WHEN p.phenotype_id IS NULL AND sg.gene_id IS NULL THEN NULL ELSE
                 ((row_number() OVER (PARTITION BY a.chr, p.phenotype_id IS NOT NULL OR sg.gene_id IS NOT NULL ORDER BY a.tss, a.gene_id) - 1)
                  // {int(cfg['nominal_bin_genes'])})::INTEGER END AS bin
        FROM ann a
        LEFT JOIN (SELECT DISTINCT gene_id FROM sperm) sg ON sg.gene_id = a.gene_id
        LEFT JOIN perm p ON p.phenotype_id = a.gene_id
        LEFT JOIN vpos v ON v.chr = p.chr AND v.position = p.position AND v.A1 = p.A1 AND v.A2 = p.A2
        LEFT JOIN ncs c ON c.phenotype_id = p.phenotype_id
        LEFT JOIN ntrans t ON t.phenotype_id = p.phenotype_id
        ORDER BY a.chr, a.tss, a.gene_id
    """).fetch_arrow_table()
    # small row groups + chr/tss statistics so a page can range-read one gene's row:
    # the app knows chr and tss from search_index and filters on them, not on gene_id alone
    write_parquet(genes, cfg.tables / "genes.parquet", 1_000, stats_columns=["chr", "tss", "gene_id"])
    missing = con.execute("SELECT count(*) FROM perm p LEFT JOIN ann a ON a.gene_id = p.phenotype_id WHERE a.gene_id IS NULL").fetchone()[0]
    n_e = con.execute(f"SELECT count(*) FROM perm WHERE {sig_col} < {thr}").fetchone()[0]
    log(f"genes: {genes.num_rows} rows, {n_e} eGenes by {sig_col} < {thr}, {missing} tested genes missing from GTF")

    # ---- splice phenotypes ----
    _bh(con, "sperm", "pval_beta")
    con.execute(f"""
        CREATE TABLE sncs AS SELECT phenotype_id, count(DISTINCT cs_id) AS n_credible_sets
        FROM read_parquet('{topchef.archive_glob(cfg, 's', 'susie')}') GROUP BY 1
    """)
    sp = con.execute(f"""
        SELECT p.phenotype_id, p.gene_id, a.symbol, p.chr, p.intron_start, p.intron_end, p.cluster_id, p.strand, a.tss,
               p.num_var, p.position AS lead_position, p.A1 AS lead_A1, p.A2 AS lead_A2, v.rsid AS lead_rsid,
               p.af AS lead_af, p.start_distance AS lead_tss_distance,
               p.slope, p.slope_se, p.pval_nominal, p.pval_perm, p.pval_beta, p.qval,
               p.{sig_col} < {thr} AS is_sqtl,
               coalesce(c.n_credible_sets, 0)::INTEGER AS n_credible_sets
        FROM sperm p
        LEFT JOIN ann a ON a.gene_id = p.gene_id
        LEFT JOIN vpos v ON v.chr = p.chr AND v.position = p.position AND v.A1 = p.A1 AND v.A2 = p.A2
        LEFT JOIN sncs c ON c.phenotype_id = p.phenotype_id
        ORDER BY p.chr, a.tss, p.gene_id, p.phenotype_id
    """).fetch_arrow_table()
    write_parquet(sp, cfg.tables / "splice_phenotypes.parquet", 1_000, stats_columns=["chr", "tss", "gene_id", "lead_position"])
    n_s = con.execute(f"SELECT count(*) FROM sperm WHERE {sig_col} < {thr}").fetchone()[0]
    nogene = con.execute("SELECT count(*) FROM sperm p LEFT JOIN ann a ON a.gene_id = p.gene_id WHERE a.gene_id IS NULL").fetchone()[0]
    log(f"splice_phenotypes: {sp.num_rows} rows, {n_s} significant, {nogene} with gene not in GTF")


def credible_sets(cfg: Config) -> None:
    con = _setup(cfg)
    t = con.execute(f"""
        WITH e AS (
            SELECT 'e' AS qtl_type, phenotype_id, phenotype_id AS gene_id, chr, position, A1, A2, af, cs_id::TINYINT AS cs_id, pip
            FROM read_parquet('{topchef.archive_glob(cfg, 'e', 'susie')}')
        ), s AS (
            SELECT 's' AS qtl_type, phenotype_id, split_part(split_part(phenotype_id, ':', 5), '.', 1) AS gene_id,
                   chr, position, A1, A2, af, cs_id::TINYINT AS cs_id, pip
            FROM read_parquet('{topchef.archive_glob(cfg, 's', 'susie')}')
        ), u AS (SELECT * FROM e UNION ALL SELECT * FROM s)
        SELECT u.qtl_type, u.phenotype_id, u.gene_id, a.symbol, u.chr, a.tss, u.position, u.A1, u.A2, v.rsid,
               u.af::FLOAT AS af, u.cs_id, u.pip::FLOAT AS pip
        FROM u LEFT JOIN ann a ON a.gene_id = u.gene_id
               LEFT JOIN vpos v ON v.chr = u.chr AND v.position = u.position AND v.A1 = u.A1 AND v.A2 = u.A2
        ORDER BY u.chr, a.tss, u.gene_id, u.phenotype_id, u.cs_id, u.pip DESC
    """).fetch_arrow_table()
    write_parquet(t, cfg.tables / "credible_sets.parquet", 2_000, stats_columns=["chr", "tss", "gene_id", "position"])
    log(f"credible_sets: {t.num_rows} rows")
