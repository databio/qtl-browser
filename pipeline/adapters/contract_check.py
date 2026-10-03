"""Check an adapter's contract tables against the rules in pipeline/CONTRACT.md.

    uv run python -m pipeline.adapters.contract_check <tables_dir>      # e.g. $QTLB_DERIVED/_tables/topchef

Works for any adapter (TOPCHeF, eQTL Catalogue, ...): it reads only the contract tables and
`ingestion.json`, never a source file. Read-only. Prints one PASS/FAIL line per rule with the count
of offending rows, then `RESULT PASS` or the failed names; exits 1 on any failure.

Rules checked: column types of every table; sites unique, on known chromosomes, sorted by
(chromosome, pos), with non-empty upper-case alleles and af in [0, 1]; phenotypes unique, `extra` a
JSON object, `gene_id` unversioned; every nominal row on its partition's chromosome, naming a site
and a phenotype, unique per (type, phenotype, site); `has_nominal` agrees with the nominal rows; one
permuted row per group, its lead phenotype in the group and its lead variant a site and among the
lead phenotype's nominal rows; credible-set variants are sites, rows unique, phenotypes known, pip in
[0, 1] and `cs_size` equal to the rows per set; the optional `trans` table (sites and phenotypes known,
unique, p and beta valid) and `gwas` table (no identical rows, known chromosomes, alleles, values in
range; rows sharing a site are counted);
`ingestion.json` sections and row counts.

Nominal tables are read one chromosome partition at a time. DuckDB sizing: QTLB_DUCKDB_MEMORY
(default 24GB), QTLB_DUCKDB_THREADS (default 4). Joins are equi-joins only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import duckdb

from ..common import CHROMS

SITES = {"chr": "VARCHAR", "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR", "af": "FLOAT",
         "ma_samples": "INTEGER", "in_cis": "BOOLEAN"}
PHENOTYPES = {"phenotype_type": "VARCHAR", "phenotype_id": "VARCHAR", "phenotype_object_id": "VARCHAR",
              "gene_id": "VARCHAR", "has_nominal": "BOOLEAN", "extra": "VARCHAR"}
NOMINAL = {"phenotype_type": "VARCHAR", "phenotype_id": "VARCHAR", "gene_id": "VARCHAR", "chr": "VARCHAR",
           "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR", "beta": "FLOAT", "se": "FLOAT", "pvalue": "DOUBLE"}
PERMUTED = {"phenotype_type": "VARCHAR", "phenotype_object_id": "VARCHAR", "phenotype_id": "VARCHAR",
            "gene_id": "VARCHAR", "n_variants": "INTEGER", "p_perm": "DOUBLE", "p_beta": "DOUBLE",
            "lead_chr": "VARCHAR", "lead_pos": "INTEGER", "lead_ref": "VARCHAR", "lead_alt": "VARCHAR"}
CREDIBLE_SETS = {"phenotype_type": "VARCHAR", "phenotype_object_id": "VARCHAR", "phenotype_id": "VARCHAR",
                 "cs_id": "SMALLINT", "chr": "VARCHAR", "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR",
                 "pip": "FLOAT", "z": "FLOAT", "cs_size": "INTEGER", "cs_min_r2": "FLOAT"}
TRANS = {"phenotype_type": "VARCHAR", "phenotype_id": "VARCHAR", "gene_id": "VARCHAR", "chr": "VARCHAR",
         "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR", "beta": "FLOAT", "se": "FLOAT", "pvalue": "DOUBLE"}
GWAS = {"chr": "VARCHAR", "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR", "beta": "DOUBLE", "se": "DOUBLE",
        "af": "DOUBLE", "pvalue": "DOUBLE", "n": "BIGINT", "rs_number": "BIGINT"}
# the keys every adapter's report shares; CONTRACT.md "The ingestion report" asks for per-type coverage
# and drop counts but not their key names, which differ between adapters today
INGESTION_SECTIONS = {"experiment_id", "allele_orientation_source", "phenotype_types", "chromosomes", "dof",
                      "significance", "source", "rows"}


class Checker:
    def __init__(self, tables: Path):
        self.T = Path(tables)
        self.con = duckdb.connect()
        self.con.execute(f"SET memory_limit = '{os.environ.get('QTLB_DUCKDB_MEMORY', '24GB')}'; "
                         f"SET threads = {int(os.environ.get('QTLB_DUCKDB_THREADS', '4'))}; "
                         "SET preserve_insertion_order = false")
        self.fails: list[str] = []

    def check(self, name: str, bad: int, detail: str = "") -> None:
        print(f"{'PASS' if bad == 0 else 'FAIL'} {name}: {bad}{' ' + detail if detail else ''}", flush=True)
        if bad:
            self.fails.append(name)

    def q(self, sql: str):
        return self.con.execute(sql).fetchone()

    def src(self, name: str) -> str:
        return f"read_parquet('{self.T / name}.parquet')"

    def types(self, src: str, want: dict) -> tuple[int, str]:
        got = {r[0]: r[1] for r in self.con.execute(f"DESCRIBE SELECT * FROM {src}").fetchall()}
        bad = {k: (got.get(k), v) for k, v in want.items() if got.get(k) != v}
        return len(bad), str(bad) if bad else ""

    def run(self) -> int:
        print(f"tables {self.T}")
        self.sites()
        self.phenotypes()
        n_nom = self.nominal()
        self.permuted()
        self.credible_sets()
        self.trans()
        self.gwas()
        self.ingestion(n_nom)
        print("RESULT", "PASS" if not self.fails else f"{len(self.fails)} FAILED: {self.fails}")
        return 1 if self.fails else 0

    def sites(self) -> None:
        c, q = self.check, self.q
        S = self.src("sites")
        c("sites types", *self.types(S, SITES))
        self.con.execute(f"CREATE TABLE s AS SELECT chr, pos, ref, alt, af, row_number() OVER () AS rn FROM {S}")
        n = q("SELECT count(*) FROM s")[0]
        print(f"sites rows {n:,}")
        c("sites unique (chr,pos,ref,alt)", n - q("SELECT count(*) FROM (SELECT DISTINCT chr, pos, ref, alt FROM s)")[0])
        order = ", ".join(f"('{x}', {i})" for i, x in enumerate(CHROMS))
        self.con.execute(f"CREATE TABLE co AS SELECT * FROM (VALUES {order}) t(chr, k)")
        c("sites chr known", q("SELECT count(*) FROM s ANTI JOIN co USING (chr)")[0])
        c("sites sorted by chr then pos", q("""SELECT count(*) FROM (SELECT co.k, s.pos, lag(co.k) OVER (ORDER BY rn) pk,
            lag(s.pos) OVER (ORDER BY rn) pp FROM s JOIN co USING (chr)) WHERE k < pk OR (k = pk AND pos < pp)""")[0])
        c("sites alleles upper case, non-empty", q("""SELECT count(*) FROM s WHERE ref IS NULL OR alt IS NULL
            OR ref = '' OR alt = '' OR ref <> upper(ref) OR alt <> upper(alt) OR ref = alt""")[0])
        c("sites af in [0,1] or NaN", q("SELECT count(*) FROM s WHERE NOT isnan(af) AND (af < 0 OR af > 1)")[0])

    def phenotypes(self) -> None:
        c, q = self.check, self.q
        P = self.src("phenotypes")
        c("phenotypes types", *self.types(P, PHENOTYPES))
        self.con.execute(f"CREATE TABLE ph AS SELECT * FROM {P}")
        c("phenotypes unique (type,id)", q("SELECT count(*) - count(DISTINCT (phenotype_type, phenotype_id)) FROM ph")[0])
        c("phenotypes extra is JSON object", q("SELECT count(*) FROM ph WHERE json_type(extra) <> 'OBJECT'")[0])
        c("phenotypes gene_id unversioned", q("SELECT count(*) FROM ph WHERE gene_id LIKE '%.%'")[0])

    def nominal(self) -> int:
        """Per chromosome partition; returns the total row count."""
        c, q, con = self.check, self.q, self.con
        files = sorted((self.T / "nominal").glob("chr=*/*.parquet"))
        if not files:
            c("nominal partitions present", 1)
            return 0
        c("nominal types", *self.types(f"read_parquet('{files[0]}', hive_partitioning = false)", NOMINAL))
        orphan = dup = n_nom = badchr = 0
        con.execute("CREATE TABLE nomids (phenotype_type VARCHAR, phenotype_id VARCHAR)")
        con.execute("CREATE TABLE lead_hit (phenotype_type VARCHAR, phenotype_object_id VARCHAR)")
        con.execute(f"CREATE TABLE pm AS SELECT * FROM {self.src('permuted')}")
        for part in sorted({f.parent for f in files}):
            ch = part.name.split("=", 1)[1]
            N = f"read_parquet('{part}/*.parquet', hive_partitioning = false)"
            n = q(f"SELECT count(*) FROM {N}")[0]
            n_nom += n
            badchr += q(f"SELECT count(*) FROM {N} WHERE chr <> '{ch}'")[0]
            orphan += q(f"""SELECT count(*) FROM {N} n ANTI JOIN (SELECT pos, ref, alt FROM s WHERE chr = '{ch}') s
                USING (pos, ref, alt)""")[0]
            dup += n - q(f"SELECT count(*) FROM (SELECT DISTINCT phenotype_type, phenotype_id, pos, ref, alt FROM {N})")[0]
            con.execute(f"INSERT INTO nomids SELECT DISTINCT phenotype_type, phenotype_id FROM {N}")
            con.execute(f"""INSERT INTO lead_hit SELECT DISTINCT p.phenotype_type, p.phenotype_object_id FROM pm p
                JOIN {N} n ON n.phenotype_type = p.phenotype_type AND n.phenotype_id = p.phenotype_id
                           AND n.chr = p.lead_chr AND n.pos = p.lead_pos AND n.ref = p.lead_ref AND n.alt = p.lead_alt
                WHERE p.lead_chr = '{ch}'""")
            print(f"  nominal {ch}: {n:,}", flush=True)
        print(f"nominal rows {n_nom:,}")
        c("nominal chr matches partition", badchr)
        c("nominal every site in sites", orphan)
        c("nominal unique (type,phenotype,site)", dup)
        c("nominal phenotypes are phenotypes", q("""SELECT count(*) FROM (SELECT DISTINCT * FROM nomids)
            ANTI JOIN ph USING (phenotype_type, phenotype_id)""")[0])
        c("phenotypes.has_nominal agrees with nominal", q("""SELECT count(*) FROM ph LEFT JOIN
            (SELECT DISTINCT *, true AS h FROM nomids) n USING (phenotype_type, phenotype_id)
            WHERE has_nominal <> coalesce(h, false)""")[0])
        return n_nom

    def permuted(self) -> None:
        c, q = self.check, self.q
        c("permuted types", *self.types(self.src("permuted"), PERMUTED))
        if "pm" not in {r[0] for r in self.con.execute("SHOW TABLES").fetchall()}:
            self.con.execute(f"CREATE TABLE pm AS SELECT * FROM {self.src('permuted')}")
        c("permuted one row per group",
          q("SELECT count(*) - count(DISTINCT (phenotype_type, phenotype_object_id)) FROM pm")[0])
        c("permuted lead phenotype is a phenotype of its group", q("""SELECT count(*) FROM pm ANTI JOIN ph
            USING (phenotype_type, phenotype_id, phenotype_object_id)""")[0])
        c("permuted lead variant in sites", q("""SELECT count(*) FROM pm p ANTI JOIN s
            ON s.chr = p.lead_chr AND s.pos = p.lead_pos AND s.ref = p.lead_ref AND s.alt = p.lead_alt""")[0])
        if "lead_hit" in {r[0] for r in self.con.execute("SHOW TABLES").fetchall()}:
            c("permuted lead variant among lead phenotype's nominal rows (where it has any)", q("""SELECT count(*)
                FROM pm p SEMI JOIN (SELECT DISTINCT * FROM nomids) n USING (phenotype_type, phenotype_id)
                ANTI JOIN lead_hit h USING (phenotype_type, phenotype_object_id)""")[0])

    def credible_sets(self) -> None:
        c, q = self.check, self.q
        CS = self.src("credible_sets")
        c("credible_sets types", *self.types(CS, CREDIBLE_SETS))
        self.con.execute(f"CREATE TABLE cs AS SELECT * FROM {CS}")
        print(f"credible_sets rows {q('SELECT count(*) FROM cs')[0]:,}")
        c("credible_sets variants in sites", q("SELECT count(*) FROM cs ANTI JOIN s USING (chr, pos, ref, alt)")[0])
        c("credible_sets unique (type,phenotype,cs,site)", q("""SELECT count(*) - count(DISTINCT
            (phenotype_type, phenotype_id, cs_id, chr, pos, ref, alt)) FROM cs""")[0])
        c("credible_sets phenotype is a phenotype of its group",
          q("SELECT count(*) FROM cs ANTI JOIN ph USING (phenotype_type, phenotype_id, phenotype_object_id)")[0])
        c("credible_sets pip in [0,1]", q("SELECT count(*) FROM cs WHERE pip < 0 OR pip > 1 OR pip IS NULL")[0])
        c("credible_sets cs_size equals rows per set (sets)", q("""SELECT count(*) FROM (SELECT phenotype_type,
            phenotype_id, cs_id, any_value(cs_size) sz, count(DISTINCT cs_size) nd, count(*) n FROM cs
            GROUP BY 1, 2, 3) WHERE sz <> n OR nd > 1""")[0])

    def trans(self) -> None:
        """Optional table: absent means no trans results."""
        if not (self.T / "trans.parquet").exists():
            print("trans: no table (no trans results)")
            return
        c, q = self.check, self.q
        T = self.src("trans")
        c("trans types", *self.types(T, TRANS))
        print(f"trans rows {q(f'SELECT count(*) FROM {T}')[0]:,}")
        c("trans every site in sites", q(f"SELECT count(*) FROM {T} t ANTI JOIN s USING (chr, pos, ref, alt)")[0])
        c("trans phenotypes are phenotypes", q(f"""SELECT count(*) FROM (SELECT DISTINCT phenotype_type, phenotype_id
            FROM {T}) ANTI JOIN ph USING (phenotype_type, phenotype_id)""")[0])
        c("trans unique (type,phenotype,site)", q(f"""SELECT count(*) - count(DISTINCT
            (phenotype_type, phenotype_id, chr, pos, ref, alt)) FROM {T}""")[0])
        c("trans pvalue in [0,1], beta finite", q(f"""SELECT count(*) FROM {T} WHERE pvalue IS NULL OR isnan(pvalue)
            OR pvalue < 0 OR pvalue > 1 OR beta IS NULL OR NOT isfinite(beta)""")[0])

    def gwas(self) -> None:
        """Optional table (an experiment's GWAS): absent means none."""
        if not (self.T / "gwas.parquet").exists():
            print("gwas: no table")
            return
        c, q = self.check, self.q
        G = self.src("gwas")
        c("gwas types", *self.types(G, GWAS))
        print(f"gwas rows {q(f'SELECT count(*) FROM {G}')[0]:,}")
        c("gwas no identical rows", q(f"SELECT count(*) - count(DISTINCT (chr, pos, ref, alt, beta, se, af, pvalue, n, rs_number)) FROM {G}")[0])
        # a site may appear twice (a source reporting it once per allele order): counted, not a failure
        print(f"gwas rows sharing a site with another row: "
              f"{q(f'SELECT count(*) FROM (SELECT count(*) OVER (PARTITION BY chr, pos, ref, alt) k FROM {G}) WHERE k > 1')[0]:,}")
        c("gwas chr known", q(f"SELECT count(*) FROM {G} ANTI JOIN co USING (chr)")[0])
        c("gwas alleles upper case, non-empty, different", q(f"""SELECT count(*) FROM {G} WHERE ref IS NULL
            OR alt IS NULL OR ref = '' OR alt = '' OR ref <> upper(ref) OR alt <> upper(alt) OR ref = alt""")[0])
        c("gwas values present, p in (0,1], af in [0,1]", q(f"""SELECT count(*) FROM {G} WHERE beta IS NULL OR se IS NULL
            OR af IS NULL OR pvalue IS NULL OR n IS NULL OR pvalue <= 0 OR pvalue > 1 OR af < 0 OR af > 1""")[0])
        c("gwas.json present", 0 if (self.T / "gwas.json").exists() else 1)
        if (self.T / "gwas_bins.parquet").exists():   # optional: the landing track's bin summary
            B = f"read_parquet('{self.T / 'gwas_bins.parquet'}')"
            c("gwas_bins chr known", q(f"SELECT count(*) FROM {B} ANTI JOIN co USING (chr)")[0])
            c("gwas_bins one row per bin, lead inside, 0 <= n_gws <= n_variants, p in (0,1]", q(f"""SELECT
                (SELECT count(*) - count(DISTINCT (chr, bin_start)) FROM {B}) + (SELECT count(*) FROM {B}
                WHERE lead_position < bin_start OR lead_position >= bin_end OR n_gws < 0 OR n_gws > n_variants
                OR n_variants < 1 OR min_p <= 0 OR min_p > 1 OR lead_ea IS NULL OR lead_beta IS NULL)""")[0])

    def ingestion(self, n_nom: int) -> None:
        c, q = self.check, self.q
        rep = json.loads((self.T / "ingestion.json").read_text())
        c("ingestion.json has required sections", len(INGESTION_SECTIONS - set(rep)),
          str(sorted(INGESTION_SECTIONS - set(rep))))
        rows = rep.get("rows") or {}
        for t in ("sites", "phenotypes", "permuted", "credible_sets"):
            c(f"ingestion rows.{t} matches table", abs(rows.get(t, -1) - q(f"SELECT count(*) FROM {self.src(t)}")[0]))
        c("ingestion rows.nominal matches tables", abs(rows.get("nominal", -1) - n_nom))
        if (self.T / "trans.parquet").exists():
            c("ingestion rows.trans matches table", abs(rows.get("trans", -1) - q(f"SELECT count(*) FROM {self.src('trans')}")[0]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("tables", type=Path, help="an adapter's contract tables dir, e.g. $QTLB_DERIVED/_tables/topchef")
    args = ap.parse_args(argv)
    return Checker(args.tables).run()


if __name__ == "__main__":
    sys.exit(main())
