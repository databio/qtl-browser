"""Steps 3-4: collect distinct tested variants, then assign rsIDs from dbSNP."""
import gzip
import shutil
import subprocess

import pyarrow as pa

from .common import CHROMS, Config, connect, log, write_parquet


def collect(cfg: Config) -> None:
    """Distinct tested variants with an `in_cis` flag. Cis files carry alleles and define the
    cis set. The trans scan is genome-wide, so trans files add positions outside every cis
    window: trans_sQTL has A1/A2, trans_eQTL has only a chr:pos variant_id, so its positions
    get null alleles unless the same position appears in trans_sQTL. A position never gets
    both an allele-bearing and a null-allele row."""
    con = connect(cfg)
    srcs = ["cis_eQTL_nominal", "cis_sQTL_nominal", "cis_eQTL_permutation", "cis_sQTL_permutation",
            "cis_eQTL_SuSiE", "cis_sQTL_SuSiE"]
    union = " UNION ALL ".join(
        f"SELECT chr, position, A1, A2 FROM read_parquet('{cfg.raw_glob(s)}')" for s in srcs
    )
    out = cfg.tmp / "variants_raw.parquet"
    log("variants_collect: scanning cis files for distinct (chr, position, A1, A2), then trans files for positions outside cis")
    con.execute(f"""
        COPY (
            WITH cis AS (SELECT DISTINCT chr, position::INTEGER AS position, A1, A2 FROM ({union})),
            cispos AS (SELECT DISTINCT chr, position FROM cis),
            ts AS (SELECT DISTINCT chr, position::INTEGER AS position, A1, A2
                   FROM read_parquet('{cfg.raw_glob('trans_sQTL')}')),
            te AS (SELECT DISTINCT split_part(variant_id, ':', 1) AS chr, split_part(variant_id, ':', 2)::INTEGER AS position
                   FROM read_parquet('{cfg.raw_glob('trans_eQTL')}')),
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
        raise SystemExit("bcftools not found; install with `brew install bcftools` (plan D6)")
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

    # 2. stream dbSNP once, keep only records at tested positions
    if not matched.exists():
        log("variants_rsid: streaming dbSNP VCF through bcftools (this is the long step)")
        with gzip.open(matched, "wt") as out:
            p = subprocess.Popen(
                ["bcftools", "query", "-T", str(targets), "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\n", str(cfg.dbsnp_vcf)],
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
    con.execute("""
        CREATE TABLE exact AS
        SELECT v.chr, v.position, v.A1, v.A2, min(d.rsid) AS rsid, min(d.rs_number) AS rs_number
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

    rg = cfg["row_group_sizes"]["variants"]
    for c in CHROMS:
        t = con.execute("SELECT * FROM variants WHERE chr = ? ORDER BY position, A1, A2", [c]).fetch_arrow_table()
        write_parquet(t, cfg.derived / "variants_by_position" / f"chr={c}" / "data.parquet", rg, stats_columns=["position"])
    t = con.execute("SELECT * FROM variants WHERE rs_number IS NOT NULL ORDER BY rs_number").fetch_arrow_table()
    write_parquet(t, cfg.derived / "variants_by_rsid.parquet", rg, stats_columns=["rs_number", "rsid"])
    log("variants_rsid: wrote variants_by_position/ and variants_by_rsid.parquet")
