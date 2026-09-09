"""Steps 10-11: manifest and validation."""
import json
import subprocess

import yaml

from .common import ROOT, Config, connect, log

TABLES = {
    "search_index": "search_index.parquet",
    "genes": "genes.parquet", "splice_phenotypes": "splice_phenotypes.parquet",
    "credible_sets": "credible_sets.parquet", "coloc": "coloc.parquet",
    "gwas_dcm_bins": "gwas_dcm_bins.parquet",
    "gwas_dcm": "gwas_dcm/*/*.parquet",
    "gene_annotation": "gene_annotation.parquet", "exons": "exons.parquet",
    "variants_by_rsid": "variants_by_rsid.parquet",
    "variants_by_position": "variants_by_position/*/*.parquet",
    "gene_detail": "gene_detail/*/*/*.parquet",
    "cis_eqtl_nominal": "cis_eqtl_nominal/*/*/*.parquet",
    "cis_sqtl_nominal": "cis_sqtl_nominal/*/*/*.parquet",
    "trans_pairs": "trans_pairs/*/*.parquet",
    "trans_by_variant": "trans_by_variant/*/*.parquet",
}
WHOLE = {"search_index", "gwas_dcm_bins"}


def manifest(cfg: Config) -> None:
    import datetime as dt
    con = connect(cfg)
    sources = yaml.safe_load((cfg.raw / "sources.yaml").read_text())
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip() or None
    except FileNotFoundError:
        commit = None
    tables = {}
    for name, rel in TABLES.items():
        pat = cfg.derived / rel
        files = sorted(cfg.derived.glob(rel)) if "*" in rel else [pat]
        rows = con.execute(f"SELECT count(*) FROM read_parquet('{pat}', hive_partitioning={'true' if '*' in rel else 'false'})").fetchone()[0]
        cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{pat}', hive_partitioning={'true' if '*' in rel else 'false'})").fetchall()]
        tables[name] = {
            "path": rel, "rows": rows, "bytes": sum(f.stat().st_size for f in files), "files": len(files),
            "columns": cols, "load": "whole" if name in WHOLE else "range",
        }
    tables["cis_sqtl_nominal"]["scope"] = "significant introns" if cfg["sqtl_nominal"] == "significant" else "all introns"
    sig_col, thr = cfg["sig_column"], cfg["sig_threshold"]
    counts = {
        "egenes": con.execute(f"SELECT count(*) FROM '{cfg.derived / 'genes.parquet'}' WHERE is_egene").fetchone()[0],
        "genes_tested": con.execute(f"SELECT count(*) FROM '{cfg.derived / 'genes.parquet'}' WHERE tested").fetchone()[0],
        "sqtl_sig_phenotypes": con.execute(f"SELECT count(*) FROM '{cfg.derived / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0],
        "sqtl_sig_genes": con.execute(f"SELECT count(DISTINCT gene_id) FROM '{cfg.derived / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0],
        # match breakdown for cis variants; trans-only rows include allele-less positions that can only match by position
        "rsid_match": dict(con.execute(f"SELECT match, count(*) FROM read_parquet('{cfg.derived}/variants_by_position/*/*.parquet') WHERE in_cis GROUP BY 1").fetchall()),
        "variants_cis": con.execute(f"SELECT count(*) FROM read_parquet('{cfg.derived}/variants_by_position/*/*.parquet') WHERE in_cis").fetchone()[0],
        "variants_trans_only": con.execute(f"SELECT count(*) FROM read_parquet('{cfg.derived}/variants_by_position/*/*.parquet') WHERE NOT in_cis").fetchone()[0],
    }
    out = {
        "built": dt.datetime.now().isoformat(timespec="seconds"),
        "pipeline_commit": commit,
        "significance_rule": f"{sig_col} < {thr}",
        "paper_counts": cfg["paper_counts"],
        "counts": counts,
        "sources": {s["name"]: {"version": s.get("version"), "description": s.get("description")} for s in sources["sources"]},
        "tables": tables,
        # small JSON files the app fetches directly, without the query engine
        "assets": {n: {"path": n, "bytes": (cfg.derived / n).stat().st_size} for n in ("coloc_loci.json", "gwas_dcm_bins.json") if (cfg.derived / n).exists()},
        "gwas_dcm": json.loads((cfg.derived / "gwas_dcm.json").read_text()) if (cfg.derived / "gwas_dcm.json").exists() else None,
    }
    (cfg.derived / "manifest.json").write_text(json.dumps(out, indent=2))
    log(f"manifest: {len(tables)} tables, eGenes={counts['egenes']}, sQTL sig={counts['sqtl_sig_phenotypes']}")


def validate(cfg: Config) -> None:
    con = connect(cfg)
    fails = []

    def check(ok: bool, msg: str):
        log(("PASS " if ok else "FAIL ") + msg)
        if not ok:
            fails.append(msg)

    # 1. row counts per nominal partition equal raw (sQTL: raw restricted to significant introns when configured)
    sig_only = cfg["sqtl_nominal"] == "significant"
    sp = cfg.derived / "splice_phenotypes.parquet"
    for qtl, raw_dir, out_dir, pat in [
        ("e", "cis_eQTL_nominal", "cis_eqtl_nominal", "topchef_{c}_MaxPC70.cis_qtl_pairs.{c}.parquet"),
        ("s", "cis_sQTL_nominal", "cis_sqtl_nominal", "topchefSplice_{c}_MaxPC25.cis_qtl_pairs.{c}.parquet"),
    ]:
        bad = []
        for c in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
            r = cfg.raw_dir(raw_dir) / pat.format(c=c)
            d = cfg.derived / out_dir / f"chr={c}" / "bin=*" / "data.parquet"
            if not (r.exists() and list(d.parent.parent.glob("bin=*/data.parquet"))):
                bad.append(f"{c}:missing")
                continue
            if qtl == "s" and sig_only:
                nr = con.execute(f"SELECT count(*) FROM '{r}' SEMI JOIN (SELECT phenotype_id FROM '{sp}' WHERE is_sqtl) USING (phenotype_id)").fetchone()[0]
            else:
                nr = con.execute(f"SELECT count(*) FROM '{r}'").fetchone()[0]
            nd = con.execute(f"SELECT count(*) FROM '{d}'").fetchone()[0]
            if nr != nd:
                bad.append(f"{c}:{nr}!={nd}")
        what = "raw significant introns" if qtl == "s" and sig_only else "raw"
        check(not bad, f"{out_dir} row counts match {what}" + (f" ({', '.join(bad)})" if bad else ""))

    # 2. every nominal gene appears in genes as tested
    n = con.execute(f"""
        SELECT count(*) FROM (SELECT DISTINCT gene_id FROM read_parquet('{cfg.derived}/cis_eqtl_nominal/*/*/*.parquet')) x
        LEFT JOIN '{cfg.derived / 'genes.parquet'}' g USING (gene_id) WHERE g.gene_id IS NULL OR NOT g.tested
    """).fetchone()[0]
    check(n == 0, f"all nominal eQTL genes present and tested in genes ({n} missing)")

    # 3. counts vs paper
    pc, tol = cfg["paper_counts"], cfg["paper_counts"]["tolerance"]
    ne = con.execute(f"SELECT count(*) FROM '{cfg.derived / 'genes.parquet'}' WHERE is_egene").fetchone()[0]
    ns = con.execute(f"SELECT count(*) FROM '{cfg.derived / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0]
    check(abs(ne - pc["egenes"]) / pc["egenes"] <= tol, f"eGenes {ne} within {tol:.0%} of paper {pc['egenes']}")
    check(abs(ns - pc["sqtl_sig_phenotypes"]) / pc["sqtl_sig_phenotypes"] <= tol, f"sQTL sig {ns} within {tol:.0%} of paper {pc['sqtl_sig_phenotypes']}")

    # 4. paper variants resolve
    for key, rs in cfg["paper_variants"].items():
        c, p = key.split(":")
        got = [r[0] for r in con.execute(f"SELECT DISTINCT rsid FROM read_parquet('{cfg.derived}/variants_by_position/chr={c}/data.parquet') WHERE position = {p}").fetchall()]
        check(rs in got, f"{key} -> {rs} (got {got})")

    # 5. row-group pruning on FLNC, in its bin file; and that every tested gene has a gene_detail row
    import pyarrow.parquet as pq
    flnc_bin = con.execute(f"SELECT bin FROM '{cfg.derived / 'genes.parquet'}' WHERE gene_id = 'ENSG00000128591'").fetchone()[0]
    flnc = cfg.derived / "cis_eqtl_nominal" / "chr=chr7" / f"bin={flnc_bin}" / "data.parquet"
    md = pq.read_metadata(flnc)
    hits = 0
    for i in range(md.num_row_groups):
        st = md.row_group(i).column(0).statistics
        if st and st.min <= "ENSG00000128591" <= st.max:
            hits += 1
    check(hits == 1, f"FLNC query touches {hits} row group(s) of {md.num_row_groups} in chr7 bin {flnc_bin} eQTL file ({md.serialized_size / 1e3:.0f} KB footer)")
    # sQTL rows are grouped per intron: one intron of CAMK2D (nine significant) touches one row group
    camk2d = "ENSG00000145349"
    camk2d_bin = con.execute(f"SELECT bin FROM '{cfg.derived / 'genes.parquet'}' WHERE gene_id = '{camk2d}'").fetchone()[0]
    introns = [r[0] for r in con.execute(f"SELECT phenotype_id FROM '{cfg.derived / 'splice_phenotypes.parquet'}' WHERE gene_id = '{camk2d}' AND is_sqtl").fetchall()]
    md = pq.read_metadata(cfg.derived / "cis_sqtl_nominal" / "chr=chr4" / f"bin={camk2d_bin}" / "data.parquet")
    pcol = [md.row_group(0).column(j).path_in_schema for j in range(md.row_group(0).num_columns)].index("phenotype_id")
    touched = {p: 0 for p in introns}
    for i in range(md.num_row_groups):
        st = md.row_group(i).column(pcol).statistics
        for p in introns:
            if st and st.min <= p <= st.max:
                touched[p] += 1
    check(all(n == 1 for n in touched.values()), f"each of CAMK2D's {len(introns)} significant introns touches one row group of {md.num_row_groups} in chr4 bin {camk2d_bin} sQTL file ({md.serialized_size / 1e3:.0f} KB footer): {sorted(set(touched.values()))}")
    n_binned = con.execute(f"SELECT count(*) FROM '{cfg.derived / 'genes.parquet'}' WHERE bin IS NOT NULL").fetchone()[0]
    n_detail = con.execute(f"SELECT count(*) FROM read_parquet('{cfg.derived}/gene_detail/*/*/*.parquet', hive_partitioning=false)").fetchone()[0]
    check(n_detail == n_binned, f"gene_detail has one row per binned (eQTL- or sQTL-tested) gene ({n_detail} vs {n_binned})")

    # 6. rsID exact rate among cis variants (trans eQTL-only positions have no alleles and can only match by position)
    vpos = f"read_parquet('{cfg.derived}/variants_by_position/*/*.parquet')"
    tot, ex = con.execute(f"SELECT count(*), sum(match = 'exact') FROM {vpos} WHERE in_cis").fetchone()
    check(ex / tot >= 0.90, f"rsID exact match rate among cis variants {ex / tot:.1%}")

    # 7. the cis side of the variant table is unchanged by adding trans positions: same row count as
    #    the cis-only build (2026-09-03), and no position carries both allele-bearing and allele-less rows
    check(tot == 8_872_723, f"cis variant rows {tot:,} (expected 8,872,723)")
    mixed = con.execute(f"SELECT count(*) FROM (SELECT chr, position FROM {vpos} GROUP BY 1, 2 HAVING bool_or(A1 IS NULL) AND bool_or(A1 IS NOT NULL))").fetchone()[0]
    check(mixed == 0, f"positions with both allele-bearing and allele-less rows: {mixed}")

    # 8. trans rows resolve to an rsID except where dbSNP has no record at the position
    tn, tnull = con.execute(f"SELECT count(*), sum(rsid IS NULL) FROM read_parquet('{cfg.derived}/trans_pairs/*/*.parquet', hive_partitioning=false)").fetchone()
    check(tnull / tn < 0.001, f"trans rows without rsID {tnull:,} of {tn:,} ({tnull / tn:.3%})")

    if fails:
        raise SystemExit(f"validate: {len(fails)} check(s) failed")
    log("validate: all checks passed")
