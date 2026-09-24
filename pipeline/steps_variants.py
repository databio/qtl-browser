"""Steps 3-4: collect distinct tested variants, then assign rsIDs from dbSNP.

Which archives hold the tested variants, and how each one names its alleles, is TOPCHeF knowledge
and lives in `adapters/topchef.py`; this file only knows how to take a distinct set of (chr,
position, A1, A2) and hang a dbSNP rsID on it.
"""
import gzip
import shutil
import subprocess

import pyarrow as pa

from .adapters import topchef
from .common import CHROMS, Config, connect, log, write_parquet


def collect(cfg: Config) -> None:
    """Distinct tested variants with an `in_cis` flag. Cis files carry alleles and define the
    cis set. The trans scan is genome-wide, so trans files add positions outside every cis
    window: trans_sQTL has A1/A2, trans_eQTL has only a chr:pos variant_id, so its positions
    get null alleles unless the same position appears in trans_sQTL. A position never gets
    both an allele-bearing and a null-allele row."""
    con = connect(cfg)
    union = " UNION ALL ".join(
        f"SELECT chr, position, A1, A2 FROM read_parquet('{cfg.raw_glob(s)}')" for s in topchef.cis_sources()
    )
    out = cfg.tmp / "variants_raw.parquet"
    log("variants_collect: scanning cis files for distinct (chr, position, A1, A2), then trans files for positions outside cis")
    con.execute(f"""
        COPY (
            WITH cis AS (SELECT DISTINCT chr, position::INTEGER AS position, A1, A2 FROM ({union})),
            cispos AS (SELECT DISTINCT chr, position FROM cis),
            ts AS (SELECT DISTINCT chr, position::INTEGER AS position, A1, A2
                   FROM read_parquet('{topchef.archive_glob(cfg, 's', 'trans')}')),
            te AS (SELECT DISTINCT split_part(variant_id, ':', 1) AS chr, split_part(variant_id, ':', 2)::INTEGER AS position
                   FROM read_parquet('{topchef.archive_glob(cfg, 'e', 'trans')}')),
            trans_alleles AS (SELECT * FROM ts ANTI JOIN cispos USING (chr, position)),
            trans_noallele AS (
                SELECT chr, position, NULL::VARCHAR AS A1, NULL::VARCHAR AS A2
                FROM te ANTI JOIN cispos USING (chr, position)
                ANTI JOIN (SELECT DISTINCT chr, position FROM ts) tsp USING (chr, position)
            )
            SELECT *, true AS in_cis FROM cis
            UNION ALL SELECT *, false FROM trans_alleles
            UNION ALL SELECT *, false FROM trans_noallele
            ORDER BY chr, position, A1, A2
        )
        TO '{out}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    for in_cis, alleles, n in con.execute(
        f"SELECT in_cis, A1 IS NOT NULL, count(*) FROM '{out}' GROUP BY 1, 2 ORDER BY 1 DESC, 2 DESC"
    ).fetchall():
        log(f"variants_collect: {'cis' if in_cis else 'trans-only'}, {'alleles' if alleles else 'no alleles'}: {n:>12,}")
    n = con.execute(f"SELECT count(*) FROM '{out}'").fetchone()[0]
    log(f"variants_collect: {n:,} distinct variants -> {out}")


def _accession_map(cfg: Config) -> dict[str, str]:
    """RefSeq accession (NC_000001.11) -> chr name, from the assembly report."""
    m = {}
    for line in cfg.assembly_report.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        # columns: Sequence-Name, Sequence-Role, Assigned-Molecule, ..., RefSeq-Accn(6), ..., UCSC-style-name(9)
        if f[1] == "assembled-molecule" and f[6] != "na":
            m[f[6]] = f[9] if f[9] != "na" else f"chr{f[2]}"
    return m


def rsid(cfg: Config) -> None:
    if shutil.which("bcftools") is None:
        raise SystemExit("bcftools not found. It comes from a bulker crate, not a package manager: "
                         "`bulker activate bulker/qtlb-format.yaml` in the qtlb-format analysis "
                         "project, which pins bcftools 1.24. On Rivanna an sbatch script must also "
                         "`module load apptainer` for bulker to find its container runtime.")
    con = connect(cfg)
    acc = _accession_map(cfg)
    chr_to_acc = {v: k for k, v in acc.items()}
    raw = cfg.tmp / "variants_raw.parquet"

    # 1. targets file for bcftools: RefSeq accession + position, sorted in VCF order. Both
    #    cached files are rebuilt when the variant list is newer than they are; an existence
    #    check alone would silently reuse a cache built from a shorter list.
    targets = cfg.tmp / "dbsnp_targets.tsv"
    matched = cfg.tmp / "dbsnp_matched.tsv.gz"
    for cache in (targets, matched):
        if cache.exists() and cache.stat().st_mtime < raw.stat().st_mtime:
            log(f"variants_rsid: {cache.name} is older than {raw.name}, rebuilding it")
            cache.unlink()
    if not targets.exists():
        log("variants_rsid: writing targets file")
        with open(targets, "w") as fh:
            for c in CHROMS:
                rows = con.execute(f"SELECT DISTINCT position FROM '{raw}' WHERE chr = ? ORDER BY position", [c]).fetchall()
                a = chr_to_acc[c]
                for (p,) in rows:
                    fh.write(f"{a}\t{p}\n")

    # 2. stream dbSNP once, keep only records at tested positions.
    #
    #    -T streams the whole 29.5 GB VCF; -R seeks with the tabix index instead. For a genome-wide
    #    build -T wins, because the targets are dense enough that seeking costs more than reading.
    #    For a CHROMS subset (a smoke build) the targets cover a few percent of the genome and -R
    #    turns nine minutes into seconds, which is the difference between a usable debug loop and a
    #    useless one. Both flags select the same records; only the access pattern differs.
    subset = len(CHROMS) < 23
    flag = "-R" if subset else "-T"
    if not matched.exists():
        log(f"variants_rsid: {'seeking' if subset else 'streaming'} the dbSNP VCF through bcftools "
            f"({flag}, {len(CHROMS)} chromosome(s))" + ("" if subset else "; this is the long step"))
        with gzip.open(matched, "wt") as out:
            p = subprocess.Popen(
                ["bcftools", "query", flag, str(targets), "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\n", str(cfg.dbsnp_vcf)],
                stdout=subprocess.PIPE, text=True,
            )
            n = 0
            for line in p.stdout:
                out.write(line)
                n += 1
            if p.wait() != 0:
                matched.unlink(missing_ok=True)
                raise SystemExit("bcftools failed")
        log(f"variants_rsid: {n:,} dbSNP records at tested positions")

    # 3. allele-aware match in DuckDB
    log("variants_rsid: matching alleles")
    acc_rows = pa.table({"acc": list(acc.keys()), "chr": list(acc.values())})
    con.register("acc_map", acc_rows)
    con.execute(f"""
        CREATE TABLE dbsnp AS
        SELECT m.chr, d.pos::INTEGER AS position,
               split_part(d.id, ';', 1) AS rsid,
               try_cast(substr(split_part(d.id, ';', 1), 3) AS BIGINT) AS rs_number,
               d.ref, unnest(string_split(d.alt, ',')) AS alt
        FROM read_csv('{matched}', delim='\t', header=false, columns={{'acc':'VARCHAR','pos':'BIGINT','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR'}}) d
        JOIN acc_map m ON m.acc = d.acc
    """)
    con.execute(f"CREATE TABLE v AS SELECT * FROM '{raw}'")
    # arg_min, not min: two separate min()s take the text of one dbSNP record and the number of
    # another, so the rsID text stops matching its own rs_number (14,803 variants before this fix).
    # `bypos` below already used arg_min.
    con.execute("""
        CREATE TABLE exact AS
        SELECT v.chr, v.position, v.A1, v.A2, arg_min(d.rsid, d.rs_number) AS rsid, min(d.rs_number) AS rs_number
        FROM v JOIN dbsnp d ON d.chr = v.chr AND d.position = v.position
         AND ((v.A2 = d.ref AND v.A1 = d.alt) OR (v.A1 = d.ref AND v.A2 = d.alt))
        GROUP BY 1,2,3,4
    """)
    con.execute("""
        CREATE TABLE bypos AS
        SELECT chr, position, arg_min(rsid, rs_number) AS rsid, min(rs_number) AS rs_number
        FROM dbsnp GROUP BY 1,2
    """)
    # null alleles (trans eQTL-only positions) never satisfy the exact join and fall through
    # to the position match
    con.execute("""
        CREATE TABLE variants AS
        SELECT v.chr, v.position, v.A1, v.A2,
               coalesce(e.rsid, b.rsid) AS rsid,
               coalesce(e.rs_number, b.rs_number) AS rs_number,
               CASE WHEN e.rsid IS NOT NULL THEN 'exact' WHEN b.rsid IS NOT NULL THEN 'position' ELSE 'none' END AS match,
               v.in_cis
        FROM v LEFT JOIN exact e USING (chr, position, A1, A2)
               LEFT JOIN bypos b USING (chr, position)
    """)
    for in_cis in (True, False):
        stats = con.execute("SELECT match, count(*) FROM variants WHERE in_cis = ? GROUP BY 1 ORDER BY 1", [in_cis]).fetchall()
        total = sum(n for _, n in stats)
        for m, n in stats:
            log(f"variants_rsid: {'cis       ' if in_cis else 'trans-only'} {m:9s} {n:>12,} ({100 * n / total:.2f}%)")

    # one build intermediate, in the order the packs read it. The browser gets the variants files
    # and the rsID index instead (SPEC sections 4 and 14), so no per-chromosome or rsID-keyed copy
    # is written any more.
    t = con.execute("SELECT * FROM variants ORDER BY chr, position, A1, A2").fetch_arrow_table()
    write_parquet(t, cfg.tables / "variants.parquet", 200_000, stats_columns=["chr", "position"])
    log(f"variants_rsid: wrote _tables/variants.parquet, {t.num_rows:,} variants")
