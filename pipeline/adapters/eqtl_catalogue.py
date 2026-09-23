"""The eQTL Catalogue adapter: release 6 TSVs -> the five contract tables (pipeline/CONTRACT.md).

    QTLB_CHROMS=chr21,chr22 QTLB_DERIVED=/scratch/$USER/qtl-browser/derived-eqtlcat \
        uv run python -m pipeline.adapters.eqtl_catalogue [--experiment gtex_v8_heart_lv]

This module is the only place that knows what an eQTL Catalogue dataset looks like on disk. The
experiment is named in `config.yaml` under `eqtl_catalogue:`; everything below is keyed by that
name, so a second tissue is a second config entry, not a second adapter.

What the source looks like (measured in analysis/qtlb-format/docs/eqtl-catalogue-survey.md):

- `ref`/`alt` with ALT the effect allele. No reorientation is expected; the refget check still runs,
  and a site whose `ref` the reference does not read is swapped (if `alt` reads) or dropped.
- ALT frequency is `ac / an`. `maf` mirrors whenever ALT is the major allele, so it is not used.
- Chromosomes are spelled unprefixed in `chromosome` and prefixed in `variant`
  (`chr1_13550_G_A`). The adapter reads the site from `variant` and counts disagreements.
- `ge` has a genome-wide nominal file (`.all.tsv.gz`). `leafcutter` does not, anywhere: its nominal
  rows are `.cc.tsv.gz`, which keeps only the introns that have a credible set (2,243 of 176,603).
  `phenotypes.has_nominal` says which, and `permuted.n_variants` is recounted for `ge` only.
- The permutation pass is per `molecular_trait_object_id` (a leafcutter **cluster**), with
  `molecular_trait_id` its lead trait (an **intron**). No results file enumerates the introns, so
  the Catalogue's leafcutter phenotype metadata (Zenodo 7850746,
  `leafcutter_<dataset>_Ensembl_105_phenotype_metadata.tsv.gz`) supplies them, and their genes: it
  lists every intron of every tested cluster, `n_traits` agreeing on all 46,249 clusters, with one
  gene per intron. So `gene_id` is published for every leafcutter phenotype after all; the survey's
  "44,226 clusters without a gene" was a statement about the results files only.

Outputs go to `_tables/<experiment>/`: `sites.parquet`, `nominal/chr=<c>/data.parquet`,
`permuted.parquet`, `credible_sets.parquet`, `phenotypes.parquet`, `ingestion.json`. Staged
per-chromosome copies of the big nominal TSVs live in `_source/` next to them, so iterating on the
tables does not re-read 3.5 GB of gzip each time.
"""
from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .. import dof as doffit
from .. import qtlstore as qs
from ..common import CHROMS, Config, connect, log, write_parquet
from ..steps_refget import _matches, load_sequence, open_store, reference_json
from .topchef import _regroup_by_phenotype

DEFAULT_EXPERIMENT = "gtex_v8_heart_lv"
ALLELE_ORIENTATION_SOURCE = "eqtl_catalogue_ref_alt"
TYPES = ("ge", "leafcutter")
DOF_SAMPLE = 1_000_000

# The summary-statistics columns, typed. `r2` is NA on every row of r6 and is not read.
NOMINAL_COLUMNS = {
    "molecular_trait_id": "VARCHAR", "chromosome": "VARCHAR", "position": "BIGINT", "ref": "VARCHAR",
    "alt": "VARCHAR", "variant": "VARCHAR", "ma_samples": "INTEGER", "maf": "DOUBLE", "pvalue": "DOUBLE",
    "beta": "DOUBLE", "se": "DOUBLE", "type": "VARCHAR", "ac": "INTEGER", "an": "INTEGER", "r2": "VARCHAR",
    "molecular_trait_object_id": "VARCHAR", "gene_id": "VARCHAR", "median_tpm": "DOUBLE", "rsid": "VARCHAR",
}
PERMUTED_COLUMNS = {
    "molecular_trait_object_id": "VARCHAR", "molecular_trait_id": "VARCHAR", "n_traits": "INTEGER",
    "n_variants": "INTEGER", "variant": "VARCHAR", "chromosome": "VARCHAR", "position": "BIGINT",
    "pvalue": "DOUBLE", "beta": "DOUBLE", "p_perm": "DOUBLE", "p_beta": "DOUBLE",
}
CS_COLUMNS = {
    "molecular_trait_id": "VARCHAR", "gene_id": "VARCHAR", "cs_id": "VARCHAR", "variant": "VARCHAR",
    "rsid": "VARCHAR", "cs_size": "INTEGER", "pip": "DOUBLE", "pvalue": "DOUBLE", "beta": "DOUBLE",
    "se": "DOUBLE", "z": "DOUBLE", "cs_min_r2": "DOUBLE", "region": "VARCHAR",
}


# ---- config and paths --------------------------------------------------------------------------
class Experiment:
    """One `eqtl_catalogue:` entry of config.yaml, with its source files resolved."""

    def __init__(self, cfg: Config, name: str = DEFAULT_EXPERIMENT, chroms: list[str] | None = None):
        block = cfg["eqtl_catalogue"][name]
        self.cfg, self.name, self.block = cfg, name, block
        self.root = Path(os.environ.get("QTLB_EQTLCAT_ROOT") or block["root"])
        self.study = block["study"]
        self.datasets = dict(block["datasets"])
        self.nominal_kind = dict(block["nominal"])
        self.n_samples = int(block["n_samples"])
        self.chroms = list(chroms or CHROMS)
        self.tables = cfg.tables / name

    # sources
    def sumstats(self, ptype: str, kind: str) -> Path:
        ds = self.datasets[ptype]
        return self.root / "sumstats" / self.study / ds / f"{ds}.{kind}.tsv.gz"

    def nominal_source(self, ptype: str) -> Path:
        return self.sumstats(ptype, self.nominal_kind[ptype])

    def permuted_source(self, ptype: str) -> Path:
        return self.sumstats(ptype, "permuted")

    def cs_source(self, ptype: str) -> Path:
        ds = self.datasets[ptype]
        return self.root / "susie" / self.study / ds / f"{ds}.credible_sets.tsv.gz"

    def leafcutter_metadata(self) -> Path:
        return Path(self.block["leafcutter_metadata"])

    # outputs
    def stage_dir(self, ptype: str) -> Path:
        return self.tables / "_source" / f"{ptype}_nominal"

    def stage_glob(self, ptype: str, chrom: str | None = None) -> str:
        return str(self.stage_dir(ptype) / f"chr={chrom or '*'}" / "*.parquet")

    def stage_sql(self, ptype: str, chrom: str | None = None) -> str:
        """The staged rows; `chr` comes back from the partition path."""
        return f"read_parquet('{self.stage_glob(ptype, chrom)}', hive_partitioning = true)"

    def orientation_path(self, chrom: str) -> Path:
        return self.tables / "_orientation" / f"chr={chrom}" / "data.parquet"

    def orientation_sql(self) -> str:
        return f"read_parquet('{self.tables / '_orientation' / 'chr=*' / 'data.parquet'}', hive_partitioning = false)"

    def nominal_path(self, chrom: str) -> Path:
        return self.tables / "nominal" / f"chr={chrom}" / "data.parquet"

    def nominal_sql(self) -> str:
        files = [self.nominal_path(c) for c in self.chroms if self.nominal_path(c).exists()]
        if not files:
            return "(SELECT NULL::VARCHAR AS phenotype_type, NULL::VARCHAR AS phenotype_id, 0 AS pos WHERE false)"
        return "read_parquet([" + ", ".join(f"'{f}'" for f in files) + "], hive_partitioning = false)"

    def path(self, table: str) -> Path:
        return self.tables / f"{table}.parquet"

    def report_path(self) -> Path:
        return self.tables / "ingestion.json"

    def significance(self) -> dict:
        s = self.block["significance"]
        return {"column": s["column"], "op": s["op"], "threshold": float(s["threshold"])}

    def chrom_list(self, column: str) -> str:
        return f"{column} IN (" + ", ".join(f"'{c}'" for c in self.chroms) + ")"


def _tsv(path: Path, columns: dict) -> str:
    cols = ", ".join(f"'{k}': '{v}'" for k, v in columns.items())
    return (f"read_csv('{path}', delim = '\t', header = true, quote = '', nullstr = 'NA', "
            f"columns = {{{cols}}})")


# `variant` is chr_pos_ref_alt with a chr prefix; alleles never contain `_`.
VARIANT_PARSE = """split_part({v}, '_', 1) AS chr, split_part({v}, '_', 2)::INTEGER AS pos,
                   split_part({v}, '_', 3) AS ref, split_part({v}, '_', 4) AS alt"""


def _variant(col: str = "variant") -> str:
    return VARIANT_PARSE.format(v=col)


# ---- staging -----------------------------------------------------------------------------------
def extract(exp: Experiment, force: bool = False) -> dict:
    """Stage each nominal TSV as `_source/<type>_nominal/chr=<c>/*.parquet`, one pass over the gzip.

    Cached: reused when a marker says it covers the requested chromosomes and the source has not
    changed since. The marker also carries the naming check: rows whose `chromosome`/`position`/
    `ref`/`alt` columns disagree with `variant` (which the adapter trusts) are counted there.
    """
    out = {}
    for ptype in TYPES:
        src = exp.nominal_source(ptype)
        marker = exp.stage_dir(ptype).with_suffix(".json")
        if marker.exists() and not force:
            m = json.loads(marker.read_text())
            if set(exp.chroms) <= set(m["chroms"]) and m["source_mtime"] == src.stat().st_mtime:
                out[ptype] = m
                continue
        t0 = time.time()
        log(f"eqtl_catalogue extract {ptype}: {src.name} ({src.stat().st_size / 1e9:.2f} GB)")
        con = connect(exp.cfg, memory_limit=exp.cfg["duckdb_memory_limit"], threads=exp.cfg["duckdb_threads"],
                      temp_dir=exp.cfg.tmp / f"eqtlcat-extract-{ptype}")
        import shutil
        shutil.rmtree(exp.stage_dir(ptype), ignore_errors=True)
        exp.stage_dir(ptype).parent.mkdir(parents=True, exist_ok=True)
        # Streamed straight to parquet: the `ge` file is 158M rows, too many to hold. The naming
        # check rides along as a column and is counted from the staged files afterwards.
        con.execute(f"""COPY (
            SELECT molecular_trait_id, molecular_trait_object_id, gene_id, variant, ma_samples,
                   pvalue, beta, se, ac, an, rsid, {_variant()},
                   ('chr' || chromosome = split_part(variant, '_', 1)
                    AND position = split_part(variant, '_', 2)::BIGINT
                    AND ref = split_part(variant, '_', 3) AND alt = split_part(variant, '_', 4)) AS naming_ok
            FROM {_tsv(src, NOMINAL_COLUMNS)}
            WHERE {exp.chrom_list("'chr' || chromosome")}
        ) TO '{exp.stage_dir(ptype)}' (FORMAT PARQUET, COMPRESSION ZSTD, PARTITION_BY (chr))""")
        naming = con.execute(f"SELECT count(*), count(*) FILTER (WHERE NOT naming_ok) "
                             f"FROM {exp.stage_sql(ptype)}").fetchone()
        con.close()
        m = {"source": str(src), "source_bytes": src.stat().st_size, "source_mtime": src.stat().st_mtime,
             "chroms": exp.chroms, "rows": int(naming[0]), "rows_variant_disagrees_with_columns": int(naming[1])}
        marker.write_text(json.dumps(m, indent=2) + "\n")
        log(f"eqtl_catalogue extract {ptype}: {m['rows']:,} rows in {(time.time() - t0) / 60:.1f} min")
        out[ptype] = m
    return out


def _load_small(exp: Experiment, con) -> None:
    """Permuted, credible-set and leafcutter metadata tables into DuckDB. Small; read whole."""
    for ptype in TYPES:
        con.execute(f"""CREATE OR REPLACE TABLE perm_{ptype} AS
            SELECT *, {_variant()}, 'chr' || chromosome AS lead_chromosome
            FROM {_tsv(exp.permuted_source(ptype), PERMUTED_COLUMNS)}""")
        # Like the nominal files, the credible-set file repeats a (set, variant) row once per dbSNP id
        # the variant carries, identical but for `rsid` (GTEx heart: 4,871 ge and 2,002 leafcutter
        # extra rows). One row each here; `n_src` keeps the count for the ingestion report, and
        # `credible_sets()` stops the run if a (set, variant) still repeats with different values.
        con.execute(f"""CREATE OR REPLACE TABLE cs_{ptype} AS
            SELECT * EXCLUDE (cs_id, rsid), cs_id AS cs_label, {_variant()}, count(*) AS n_src
            FROM {_tsv(exp.cs_source(ptype), CS_COLUMNS)} GROUP BY ALL""")
    con.execute(f"""CREATE OR REPLACE TABLE meta_lc AS
        SELECT * FROM read_csv('{exp.leafcutter_metadata()}', delim = '\t', header = true, quote = '',
                               nullstr = 'NA', all_varchar = true)""")


# ---- the reference check -----------------------------------------------------------------------
def allele_reads(seq: np.ndarray, pos: np.ndarray, ref: list, alt: list) -> tuple[np.ndarray, np.ndarray]:
    """Does the reference read `ref`, and does it read `alt`, at each 1-based `pos`?

    Asked separately rather than through `steps_refget.classify`, because `classify` breaks an
    indel tie (both alleles prefix-match) toward the longer allele. That tie-break is right for
    TOPCHeF's unordered A1/A2 and wrong here: the Catalogue names its REF, and `A>AT` where the
    genome reads `AT` is correctly anchored as given. So the source `ref` stands wherever it reads.
    """
    n = len(pos)
    ref_ok = np.zeros(n, dtype=bool)
    alt_ok = np.zeros(n, dtype=bool)
    snp = np.array([len(r) == 1 and len(a) == 1 for r, a in zip(ref, alt)], dtype=bool)
    inside = (pos >= 1) & (pos <= len(seq))
    idx = np.flatnonzero(snp & inside)
    if len(idx):
        base = seq[pos[idx] - 1]
        ref_ok[idx] = base == np.frombuffer("".join(ref[i] for i in idx).upper().encode(), dtype=np.uint8)
        alt_ok[idx] = base == np.frombuffer("".join(alt[i] for i in idx).upper().encode(), dtype=np.uint8)
    for i in np.flatnonzero(~snp):
        p = int(pos[i])
        ref_ok[i] = _matches(seq, p, ref[i])
        alt_ok[i] = _matches(seq, p, alt[i])
    return ref_ok, alt_ok


def reference_base(ref: np.ndarray, alt: np.ndarray, ref_ok: np.ndarray, alt_ok: np.ndarray) -> np.ndarray:
    """The allele the reference reads, preferring the source's `ref`; None where neither reads."""
    out = np.where(ref_ok, ref, np.where(alt_ok, alt, None)).astype(object)
    return out


def sequence_loader(cfg: Config) -> Callable[[str], np.ndarray]:
    """chrom -> upper-case uint8 sequence from the refgetstore. Digests come from
    `_tables/reference.json` when the tree has one, else straight from the configured collection."""
    store = open_store(cfg)
    if reference_json(cfg).exists():
        seqs = json.loads(reference_json(cfg).read_text())["sequences"]
        table = {c: (s["digest"], s["length"]) for c, s in seqs.items()}
    else:
        coll = cfg["reference"]["collection"]
        table = {r.metadata.name: (r.metadata.sha512t24u, r.metadata.length) for r in store.get_collection(coll)}

    def load(chrom: str) -> np.ndarray:
        digest, length = table[chrom]
        return load_sequence(store, digest, length)
    return load


# ---- sites -------------------------------------------------------------------------------------
def sites(exp: Experiment, con, sequence: Callable[[str], np.ndarray]) -> dict:
    """`sites.parquet` and the internal `_orientation/` table, with the orientation counts.

    A site is any variant in either nominal table, any lead variant, any credible-set variant.
    `af = ac/an` and `ma_samples` from the `ge` rows where the variant has any, else `leafcutter`;
    a variant only a lead or credible set names has none (`af` NaN, counts -1), which the report
    counts. `ma_count = min(ac, an - ac)`. Everything here is a cis result, so `in_cis` is true.
    """
    report = {"as_is": 0, "swapped": 0, "dropped_neither_reads": 0, "by_origin": {},
              "af_disagrees_between_types": 0, "stats_disagree_within_type": 0,
              # the source repeats a row per dbSNP id; `rsid` keeps the smallest string
              "sites_with_several_rsids": 0}
    parts = []
    for chrom in exp.chroms:
        per_type = []
        for ptype in TYPES:
            if not list(exp.stage_dir(ptype).glob(f"chr={chrom}/*.parquet")):
                con.execute(f"""CREATE OR REPLACE TABLE v_{ptype} AS SELECT NULL::INTEGER AS pos,
                    NULL::VARCHAR AS ref, NULL::VARCHAR AS alt, NULL::INTEGER AS ac, NULL::INTEGER AS an,
                    NULL::INTEGER AS ma_samples, NULL::VARCHAR AS rsid, 0 AS n_rsid, false AS differs WHERE false""")
            else:
                con.execute(f"""CREATE OR REPLACE TABLE v_{ptype} AS
                    SELECT pos, ref, alt, min(ac) AS ac, min(an) AS an, min(ma_samples) AS ma_samples,
                           min(rsid) AS rsid, count(DISTINCT rsid) AS n_rsid,
                           count(DISTINCT ac) > 1 OR count(DISTINCT an) > 1 OR count(DISTINCT ma_samples) > 1 AS differs
                    FROM {exp.stage_sql(ptype, chrom)} GROUP BY 1, 2, 3""")
            per_type.append(f"v_{ptype}")
        report["stats_disagree_within_type"] += con.execute(
            "SELECT (SELECT count(*) FILTER (WHERE differs) FROM v_ge) + (SELECT count(*) FILTER (WHERE differs) FROM v_leafcutter)").fetchone()[0]
        report["sites_with_several_rsids"] += con.execute(
            "SELECT count(*) FROM (SELECT pos, ref, alt FROM v_ge WHERE n_rsid > 1 UNION SELECT pos, ref, alt FROM v_leafcutter WHERE n_rsid > 1)").fetchone()[0]
        extra = " UNION ".join(
            [f"SELECT pos, ref, alt, 'lead' AS origin FROM perm_{t} WHERE chr = '{chrom}'" for t in TYPES]
            + [f"SELECT pos, ref, alt, 'credible_set' AS origin FROM cs_{t} WHERE chr = '{chrom}'" for t in TYPES])
        t = con.execute(f"""
            WITH other AS (SELECT pos, ref, alt, min(origin) AS origin FROM ({extra}) GROUP BY 1, 2, 3),
            allv AS (SELECT pos, ref, alt FROM v_ge UNION SELECT pos, ref, alt FROM v_leafcutter
                     UNION SELECT pos, ref, alt FROM other)
            SELECT a.pos, a.ref, a.alt,
                   coalesce(g.ac::DOUBLE / g.an, l.ac::DOUBLE / l.an) AS af,
                   coalesce(g.ma_samples, l.ma_samples, -1)::INTEGER AS ma_samples,
                   coalesce(least(g.ac, g.an - g.ac), least(l.ac, l.an - l.ac), -1)::INTEGER AS ma_count,
                   coalesce(g.rsid, l.rsid) AS rsid,
                   CASE WHEN g.pos IS NOT NULL THEN 'ge_nominal' WHEN l.pos IS NOT NULL THEN 'leafcutter_nominal'
                        ELSE o.origin || '_only' END AS origin,
                   (g.pos IS NOT NULL AND l.pos IS NOT NULL AND g.ac::DOUBLE / g.an <> l.ac::DOUBLE / l.an) AS af_differs
            FROM allv a
            LEFT JOIN v_ge g USING (pos, ref, alt)
            LEFT JOIN v_leafcutter l USING (pos, ref, alt)
            LEFT JOIN other o USING (pos, ref, alt)
            ORDER BY a.pos, a.ref, a.alt""").fetch_arrow_table()
        report["af_disagrees_between_types"] += int(pc.sum(pc.cast(t["af_differs"], pa.int64())).as_py() or 0)
        for k, v in zip(*np.unique(np.array(t["origin"].to_pylist(), dtype=object), return_counts=True)):
            report["by_origin"][k] = report["by_origin"].get(k, 0) + int(v)
        pos = t["pos"].to_numpy().astype(np.int64)
        ref = t["ref"].to_pylist()
        alt = t["alt"].to_pylist()
        ref_ok, alt_ok = allele_reads(sequence(chrom), pos, ref, alt)
        af = t["af"].to_numpy(zero_copy_only=False).astype(np.float64)
        r = qs.orient_to_ref(reference_base(np.array(ref, dtype=object), np.array(alt, dtype=object), ref_ok, alt_ok),
                             np.array(alt, dtype=object), np.array(ref, dtype=object), np.ones(len(pos)), af)
        keep = r["keep"]
        report["as_is"] += r["counts"]["as_is"]
        report["swapped"] += r["counts"]["swapped"]
        report["dropped_neither_reads"] += r["counts"]["dropped"]
        swapped = r["beta"] < 0
        n = int(keep.sum())
        write_parquet(pa.table({
            "chr": pa.array([chrom] * n, pa.string()), "pos": pa.array(pos[keep], pa.int32()),
            "src_ref": pa.array([x for x, k in zip(ref, keep) if k], pa.string()),
            "src_alt": pa.array([x for x, k in zip(alt, keep) if k], pa.string()),
            "ref": pa.array(list(r["ref"]), pa.string()), "alt": pa.array(list(r["alt"]), pa.string()),
            "swapped": pa.array(swapped, pa.bool_()),
        }), exp.orientation_path(chrom), 200_000, stats_columns=["chr", "pos"])
        rsid = t["rsid"].filter(pa.array(keep))
        rs_num = pc.if_else(pc.is_null(rsid), -1,
                            pc.cast(pc.utf8_slice_codeunits(pc.fill_null(rsid, "rs-1"), 2), pa.int64()))
        part = pa.table({
            "chr": pa.array([chrom] * n, pa.string()), "pos": pa.array(pos[keep], pa.int32()),
            "ref": pa.array(list(r["ref"]), pa.string()), "alt": pa.array(list(r["alt"]), pa.string()),
            "af": pa.array(r["af"], pa.float32()),
            "ma_samples": t["ma_samples"].filter(pa.array(keep)),
            "in_cis": pa.array(np.ones(n, dtype=bool), pa.bool_()),
            "ma_count": t["ma_count"].filter(pa.array(keep)),
            "rsid": rsid, "rs_number": rs_num,
        })
        # (chr, pos) order within the chromosome; a swap can reorder alleles at one position
        parts.append(part.sort_by([("pos", "ascending"), ("ref", "ascending"), ("alt", "ascending")]))
        log(f"eqtl_catalogue sites {chrom}: {n:,} kept of {len(pos):,}")
    tbl = pa.concat_tables(parts)
    write_parquet(tbl, exp.path("sites"), 200_000, stats_columns=["chr", "pos"])
    report["rows"] = tbl.num_rows
    return report


# ---- phenotypes --------------------------------------------------------------------------------
def phenotypes(exp: Experiment, con) -> dict:
    """`phenotypes.parquet`, plus the `pheno` DuckDB table the other tables join to for gene ids.

    `ge`: one row per permuted gene; object = trait = gene. `leafcutter`: every intron the
    metadata lists for a cluster the permutation pass tested. Introns of clusters it did not test
    are counted, not written. `has_nominal` is read from the staged nominal rows.
    """
    rep: dict = {}
    for ptype in TYPES:
        con.execute(f"""CREATE OR REPLACE TABLE nomids_{ptype} AS
            SELECT DISTINCT molecular_trait_id AS phenotype_id
            FROM {exp.stage_sql(ptype)}
            WHERE {exp.chrom_list('chr')}""")
    # The gene list of an intron: the metadata's `gene_id`, split in case one names several.
    con.execute(f"""CREATE OR REPLACE TABLE pheno AS
        SELECT 'ge' AS phenotype_type, molecular_trait_id AS phenotype_id,
               molecular_trait_object_id AS phenotype_object_id,
               split_part(molecular_trait_id, '.', 1) AS gene_id, [split_part(molecular_trait_id, '.', 1)] AS gene_ids,
               '{{}}' AS extra, lead_chromosome AS chr
        FROM perm_ge WHERE {exp.chrom_list('lead_chromosome')}
        UNION ALL
        SELECT 'leafcutter', m.phenotype_id, m.group_id,
               nullif(split_part(trim(string_split(regexp_replace(m.gene_id, ';', ',', 'g'), ',')[1]), '.', 1), ''),
               list_transform(string_split(regexp_replace(coalesce(m.gene_id, ''), ';', ',', 'g'), ','),
                              x -> split_part(trim(x), '.', 1)),
               to_json({{'intron_start': m.intron_start::INTEGER, 'intron_end': m.intron_end::INTEGER,
                         'cluster_id': m.group_id, 'strand': right(m.group_id, 1)}})::VARCHAR,
               'chr' || m.chromosome
        FROM meta_lc m SEMI JOIN perm_leafcutter p ON p.molecular_trait_object_id = m.group_id
        WHERE {exp.chrom_list("'chr' || m.chromosome")}""")
    multi = con.execute("SELECT count(*) FROM pheno WHERE len(list_filter(gene_ids, x -> x <> '')) > 1").fetchone()[0]
    t = con.execute("""
        SELECT p.phenotype_type, p.phenotype_id, p.phenotype_object_id, p.gene_id,
               (n1.phenotype_id IS NOT NULL OR n2.phenotype_id IS NOT NULL) AS has_nominal,
               CASE WHEN len(list_filter(p.gene_ids, x -> x <> '')) > 1
                    THEN json_merge_patch(p.extra, to_json({'gene_ids': p.gene_ids}))::VARCHAR ELSE p.extra END AS extra
        FROM pheno p
        LEFT JOIN nomids_ge n1 ON p.phenotype_type = 'ge' AND n1.phenotype_id = p.phenotype_id
        LEFT JOIN nomids_leafcutter n2 ON p.phenotype_type = 'leafcutter' AND n2.phenotype_id = p.phenotype_id
        ORDER BY 1, 2""").fetch_arrow_table()
    write_parquet(t, exp.path("phenotypes"), 10_000, stats_columns=["phenotype_type", "phenotype_id", "gene_id"])
    for ptype in TYPES:
        row = con.execute(f"""SELECT count(*), count(DISTINCT phenotype_object_id), count(gene_id)
            FROM pheno WHERE phenotype_type = '{ptype}'""").fetchone()
        has = con.execute(f"SELECT count(*) FROM pheno p SEMI JOIN nomids_{ptype} n USING (phenotype_id) "
                          f"WHERE p.phenotype_type = '{ptype}'").fetchone()[0]
        orphan = con.execute(f"SELECT count(*) FROM nomids_{ptype} n ANTI JOIN "
                             f"(SELECT phenotype_id FROM pheno WHERE phenotype_type = '{ptype}') USING (phenotype_id)").fetchone()[0]
        if orphan:
            raise ValueError(f"eqtl_catalogue phenotypes {ptype}: {orphan:,} nominal traits are not phenotypes")
        rep[ptype] = {"phenotypes": row[0], "groups": row[1], "with_gene_id": row[2], "has_nominal": has,
                      "nominal_coverage": round(has / row[0], 6) if row[0] else None}
    # the leafcutter enumeration: metadata against the permutation pass
    rep["leafcutter"]["metadata"] = dict(zip(
        ("metadata_introns", "metadata_clusters", "tested_clusters_missing_from_metadata",
         "n_traits_disagree_with_metadata", "untested_clusters_in_metadata", "untested_introns_in_metadata"),
        con.execute(f"""WITH mc AS (SELECT group_id, count(*) AS n FROM meta_lc
                                    WHERE {exp.chrom_list("'chr' || chromosome")} GROUP BY 1),
                         pc AS (SELECT * FROM perm_leafcutter WHERE {exp.chrom_list('lead_chromosome')})
            SELECT (SELECT sum(n) FROM mc)::BIGINT, (SELECT count(*) FROM mc),
                   (SELECT count(*) FROM pc ANTI JOIN mc ON mc.group_id = pc.molecular_trait_object_id),
                   (SELECT count(*) FROM pc JOIN mc ON mc.group_id = pc.molecular_trait_object_id WHERE mc.n <> pc.n_traits),
                   (SELECT count(*) FROM mc ANTI JOIN perm_leafcutter p ON p.molecular_trait_object_id = mc.group_id),
                   (SELECT coalesce(sum(n), 0) FROM mc ANTI JOIN perm_leafcutter p ON p.molecular_trait_object_id = mc.group_id)::BIGINT
        """).fetchone()))
    rep["leafcutter"]["metadata"]["source"] = str(exp.leafcutter_metadata())
    rep["multi_gene_phenotypes"] = int(multi)
    log(f"eqtl_catalogue phenotypes: {t.num_rows:,} rows; {rep}")
    return rep


# ---- permuted ----------------------------------------------------------------------------------
def permuted(exp: Experiment, con) -> dict:
    """`permuted.parquet`: one row per group, `phenotype_id` its lead trait.

    `n_variants` is recounted from the contract nominal table for `ge`, whose `.all` file is the
    complete window, and copied from the source for `leafcutter`, whose nominal rows are per intron
    and exist only for fine-mapped introns, so no cluster's window can be counted from them.
    A lead variant missing from `sites` (the reference check dropped it) stops the run.
    """
    rep = {}
    parts = []
    for ptype in TYPES:
        recount = exp.nominal_kind[ptype] == "all"
        n, distinct = con.execute(f"""SELECT count(*), count(DISTINCT molecular_trait_object_id)
            FROM perm_{ptype} WHERE {exp.chrom_list('lead_chromosome')}""").fetchone()
        if n != distinct:
            raise ValueError(f"eqtl_catalogue permuted {ptype}: {n:,} rows for {distinct:,} groups")
        counts = (f"(SELECT phenotype_id, count(*)::INTEGER AS n FROM {exp.nominal_sql()} "
                  f"WHERE phenotype_type = '{ptype}' GROUP BY 1)")
        t = con.execute(f"""
            SELECT '{ptype}' AS phenotype_type, p.molecular_trait_object_id AS phenotype_object_id,
                   p.molecular_trait_id AS phenotype_id, ph.gene_id,
                   {'coalesce(c.n, p.n_variants)' if recount else 'p.n_variants'}::INTEGER AS n_variants,
                   p.p_perm, p.p_beta,
                   p.chr AS lead_chr, p.pos::INTEGER AS lead_pos, o.ref AS lead_ref, o.alt AS lead_alt,
                   c.n AS n_nominal, p.n_variants AS n_source,
                   ph.phenotype_object_id AS ph_object
            FROM perm_{ptype} p
            LEFT JOIN {exp.orientation_sql()} o ON o.chr = p.chr AND o.pos = p.pos AND o.src_ref = p.ref AND o.src_alt = p.alt
            LEFT JOIN pheno ph ON ph.phenotype_type = '{ptype}' AND ph.phenotype_id = p.molecular_trait_id
            LEFT JOIN {counts} c ON c.phenotype_id = p.molecular_trait_id
            WHERE {exp.chrom_list('p.lead_chromosome')}
            ORDER BY p.molecular_trait_object_id""").fetch_arrow_table()
        bad_obj = int(pc.sum(pc.cast(pc.invert(pc.fill_null(pc.equal(t["ph_object"], t["phenotype_object_id"]), False)), pa.int64())).as_py() or 0)
        if bad_obj:
            raise ValueError(f"eqtl_catalogue permuted {ptype}: {bad_obj:,} lead traits are not phenotypes of their group")
        n_nom, n_src = t["n_nominal"], t["n_source"]
        has = pc.is_valid(n_nom)
        rep[ptype] = {
            "groups": t.num_rows,
            "n_variants": "recounted_from_nominal" if recount else "copied_from_source",
            "groups_with_nominal_for_lead": int(pc.sum(pc.cast(has, pa.int64())).as_py() or 0),
            # lead-trait nominal rows against the group's source count; only meaningful where a
            # group is one trait (ge), so left out otherwise
            "recount_differs_from_source": (int(pc.sum(pc.cast(pc.and_(has, pc.not_equal(n_nom, n_src)), pa.int64())).as_py() or 0)
                                            if recount else None),
            "lead_variant_missing_from_sites": int(pc.sum(pc.cast(pc.is_null(t["lead_ref"]), pa.int64())).as_py() or 0),
        }
        if recount and rep[ptype]["groups_with_nominal_for_lead"] != t.num_rows:
            rep[ptype]["copied_where_no_nominal"] = t.num_rows - rep[ptype]["groups_with_nominal_for_lead"]
        parts.append(t.drop_columns(["n_nominal", "n_source", "ph_object"]))
    t = pa.concat_tables(parts)
    missing = sum(r["lead_variant_missing_from_sites"] for r in rep.values())
    if missing:
        # CONTRACT.md: the lead variant must exist in sites. A lead the reference check dropped
        # would need a rule (re-pick the lead, or drop the group) that has not been decided; none
        # occurs in GTEx heart (0 of 68,122 genome-wide), so stop rather than write a broken row.
        raise ValueError(f"eqtl_catalogue permuted: {missing:,} lead variants are not in sites")
    write_parquet(t, exp.path("permuted"), 10_000, stats_columns=["phenotype_type", "phenotype_object_id", "phenotype_id"])
    log(f"eqtl_catalogue permuted: {t.num_rows:,} rows; {rep}")
    return rep


# ---- credible sets -----------------------------------------------------------------------------
CS_LABEL = re.compile(r"^(?P<trait>.+)_L(?P<k>\d+)$")


def credible_sets(exp: Experiment, con) -> dict:
    """`credible_sets.parquet`: `cs_id` is the `_L<k>` suffix of the source's string id. The
    Catalogue fine-maps per trait (per intron for leafcutter), so `phenotype_id` is the trait and
    `phenotype_object_id` its group from `phenotypes`."""
    rep = {}
    parts = []
    for ptype in TYPES:
        t = con.execute(f"""
            SELECT '{ptype}' AS phenotype_type, ph.phenotype_object_id, c.molecular_trait_id AS phenotype_id,
                   regexp_extract(c.cs_label, '_L([0-9]+)$', 1) AS k,
                   c.cs_label,
                   c.chr, c.pos::INTEGER AS pos, o.ref, o.alt, c.pip::FLOAT AS pip,
                   CASE WHEN o.swapped THEN -c.z ELSE c.z END::FLOAT AS z,
                   c.cs_size::INTEGER AS cs_size, coalesce(c.cs_min_r2, 'nan'::DOUBLE)::FLOAT AS cs_min_r2,
                   count(*) OVER (PARTITION BY c.cs_label) AS n_rows,
                   count(*) OVER (PARTITION BY c.cs_label, c.variant) AS n_same_site,
                   c.n_src
            FROM cs_{ptype} c
            LEFT JOIN {exp.orientation_sql()} o ON o.chr = c.chr AND o.pos = c.pos AND o.src_ref = c.ref AND o.src_alt = c.alt
            LEFT JOIN pheno ph ON ph.phenotype_type = '{ptype}' AND ph.phenotype_id = c.molecular_trait_id
            WHERE {exp.chrom_list('c.chr')}
            ORDER BY c.molecular_trait_id, TRY_CAST(k AS INTEGER), c.pip DESC""").fetch_arrow_table()
        labels = t["cs_label"].to_pylist()
        pids = t["phenotype_id"].to_pylist()
        bad_label = sum(1 for lab, pid in zip(labels, pids)
                        if not (m := CS_LABEL.match(lab)) or m["trait"] != pid)
        if bad_label:
            raise ValueError(f"eqtl_catalogue credible_sets {ptype}: {bad_label:,} cs_id values are not <trait>_L<k>")
        no_pheno = int(pc.sum(pc.cast(pc.is_null(t["phenotype_object_id"]), pa.int64())).as_py() or 0)
        if no_pheno:
            raise ValueError(f"eqtl_catalogue credible_sets {ptype}: {no_pheno:,} rows name a trait that is not a phenotype")
        conflict = int(pc.sum(pc.cast(pc.greater(t["n_same_site"], 1), pa.int64())).as_py() or 0)
        if conflict:
            raise ValueError(f"eqtl_catalogue credible_sets {ptype}: {conflict:,} rows repeat a (set, variant) "
                             f"with different values")
        missing = pc.is_null(t["ref"])
        size_off = pc.not_equal(t["cs_size"], pc.cast(t["n_rows"], pa.int32()))
        rep[ptype] = {
            "rows": t.num_rows,
            "sets": len(set(labels)),
            "traits": len(set(pids)),
            "duplicate_rows_rsid_only": int(pc.sum(pc.subtract(t["n_src"], 1)).as_py() or 0),
            "cs_size_differs_from_rows": int(pc.sum(pc.cast(size_off, pa.int64())).as_py() or 0),
            "sets_cs_size_differs_from_rows": len(set(t["cs_label"].filter(size_off).to_pylist())),
            "rows_dropped_site_not_in_sites": int(pc.sum(pc.cast(missing, pa.int64())).as_py() or 0),
        }
        t = t.filter(pc.invert(missing))
        t = t.set_column(t.schema.get_field_index("k"), "cs_id", pc.cast(t["k"], pa.int16()))
        parts.append(t.drop_columns(["cs_label", "n_rows", "n_same_site", "n_src"]).select(
            ["phenotype_type", "phenotype_object_id", "phenotype_id", "cs_id", "chr", "pos", "ref", "alt",
             "pip", "z", "cs_size", "cs_min_r2"]))
    t = pa.concat_tables(parts)
    write_parquet(t, exp.path("credible_sets"), 2_000, stats_columns=["phenotype_type", "phenotype_id", "chr", "pos"])
    log(f"eqtl_catalogue credible_sets: {t.num_rows:,} rows; {rep}")
    return rep


# ---- nominal -----------------------------------------------------------------------------------
def nominal(exp: Experiment, con) -> dict:
    """`nominal/chr=<c>/data.parquet`, one row group per phenotype (topchef's layout, reused).
    Rows whose site the reference check dropped are counted, not written."""
    rep = {"rows": 0, "rows_dropped_site": 0, "duplicate_rows_rsid_only": 0, "by_type": {t: 0 for t in TYPES}}
    for chrom in exp.chroms:
        selects, raw = [], 0
        for ptype in TYPES:
            if not list(exp.stage_dir(ptype).glob(f"chr={chrom}/*.parquet")):
                continue
            # The source repeats a (trait, variant) row once per dbSNP id the variant carries: same
            # statistics, different `rsid` (184,289 such rows on chr21-22 alone). One row each here;
            # rows that repeat with *different* statistics would be a real conflict, and stop the run.
            con.execute(f"""CREATE OR REPLACE TABLE nsrc AS
                SELECT molecular_trait_id, chr, pos, ref, alt, any_value(beta) AS beta, any_value(se) AS se,
                       any_value(pvalue) AS pvalue, count(*) AS k, count(DISTINCT (beta, se, pvalue)) AS nd
                FROM {exp.stage_sql(ptype, chrom)} GROUP BY 1, 2, 3, 4, 5""")
            r0, dup, conflict = con.execute("SELECT sum(k), sum(k - 1), count(*) FILTER (WHERE nd > 1) FROM nsrc").fetchone()
            if conflict:
                raise ValueError(f"eqtl_catalogue nominal {chrom} {ptype}: {conflict:,} (trait, variant) pairs repeat with different statistics")
            rep["duplicate_rows_rsid_only"] += int(dup)
            con.execute(f"CREATE OR REPLACE TABLE nsrc_{ptype} AS SELECT * EXCLUDE (k, nd) FROM nsrc")
            src = f"nsrc_{ptype}"
            raw += int(r0) - int(dup)
            selects.append(f"""
                SELECT '{ptype}' AS phenotype_type, n.molecular_trait_id AS phenotype_id, ph.gene_id,
                       n.chr, n.pos, o.ref, o.alt,
                       (CASE WHEN o.swapped THEN -n.beta ELSE n.beta END)::FLOAT AS beta,
                       n.se::FLOAT AS se, n.pvalue
                FROM {src} n
                JOIN {exp.orientation_sql()} o ON o.chr = n.chr AND o.pos = n.pos AND o.src_ref = n.ref AND o.src_alt = n.alt
                LEFT JOIN pheno ph ON ph.phenotype_type = '{ptype}' AND ph.phenotype_id = n.molecular_trait_id""")
        if not selects:
            continue
        out = exp.nominal_path(chrom)
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_name("sorted.tmp.parquet")
        con.execute("SET preserve_insertion_order = true")
        con.execute(f"""COPY ({' UNION ALL '.join(selects)} ORDER BY phenotype_type, gene_id, phenotype_id, pos)
            TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000)""")
        con.execute("SET preserve_insertion_order = false")
        n = pq.read_metadata(tmp).num_rows
        groups = _regroup_by_phenotype(tmp, out)
        tmp.unlink()
        rep["rows"] += n
        rep["rows_dropped_site"] += raw - n
        for ptype, k in con.execute(f"SELECT phenotype_type, count(*) FROM read_parquet('{out}') GROUP BY 1").fetchall():
            rep["by_type"][ptype] += k
        log(f"eqtl_catalogue nominal {chrom}: {n:,} rows ({raw - n:,} dropped) in {groups:,} row groups")
    return rep


# ---- degrees of freedom ------------------------------------------------------------------------
def fit_dof(exp: Experiment, con) -> dict:
    """Fit dof per phenotype type from the contract nominal rows (`pipeline/dof.py`)."""
    out = {}
    for ptype in TYPES:
        t = con.execute(f"""SELECT beta::DOUBLE AS beta, se::DOUBLE AS se, pvalue FROM {exp.nominal_sql()}
            WHERE phenotype_type = '{ptype}' USING SAMPLE reservoir({DOF_SAMPLE} ROWS) REPEATABLE (0)""").fetch_arrow_table()
        if t.num_rows == 0:
            out[ptype] = {"dof": None, "usable": False, "reason": "no nominal rows"}
            continue
        r = doffit.fit(t["beta"].to_numpy(), t["se"].to_numpy(), t["pvalue"].to_numpy(), n_samples=exp.n_samples)
        out[ptype] = {k: r.get(k) for k in ("dof", "usable", "identified", "residual_log10p", "margin",
                                            "rows", "rows_usable", "n_samples", "implied_covariates", "reason")}
        out[ptype]["dof_for_manifest"] = doffit.dof_for_manifest(r)
    log(f"eqtl_catalogue dof: {out}")
    return out


# ---- the whole adapter -------------------------------------------------------------------------
def run(exp: Experiment, sequence: Callable[[str], np.ndarray] | None = None, force_extract: bool = False) -> dict:
    t0 = time.time()
    staged = extract(exp, force=force_extract)
    sequence = sequence or sequence_loader(exp.cfg)
    con = connect(exp.cfg, memory_limit=exp.cfg["duckdb_memory_limit"], threads=exp.cfg["duckdb_threads"],
                  temp_dir=exp.cfg.tmp / f"eqtlcat-{exp.name}")
    _load_small(exp, con)
    report: dict = {"experiment_id": exp.name, "allele_orientation_source": ALLELE_ORIENTATION_SOURCE,
                    "chromosomes": exp.chroms, "phenotype_types": list(TYPES)}
    report["orientation"] = sites(exp, con, sequence)
    report["phenotypes"] = phenotypes(exp, con)
    report["nominal"] = nominal(exp, con)
    report["permuted"] = permuted(exp, con)
    report["credible_sets"] = credible_sets(exp, con)
    report["dof"] = fit_dof(exp, con)
    report["significance"] = exp.significance()
    report["source"] = {
        "root": str(exp.root), "study": exp.study, "datasets": exp.datasets,
        "nominal_files": {t: exp.nominal_kind[t] for t in TYPES},
        "staged": {t: {k: staged[t][k] for k in ("source", "source_bytes", "rows", "rows_variant_disagrees_with_columns")}
                   for t in TYPES},
        "leafcutter_metadata": str(exp.leafcutter_metadata()),
        "gene_annotation": exp.block.get("gene_annotation"),
    }
    report["rows"] = {k: pq.read_metadata(exp.path(k)).num_rows
                      for k in ("sites", "phenotypes", "permuted", "credible_sets")}
    report["rows"]["nominal"] = report["nominal"]["rows"]
    exp.report_path().write_text(json.dumps(report, indent=2, default=str) + "\n")
    log(f"eqtl_catalogue {exp.name}: {report['rows']} in {(time.time() - t0) / 60:.1f} min -> {exp.report_path()}")
    return report


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--experiment", default=DEFAULT_EXPERIMENT)
    ap.add_argument("--force-extract", action="store_true", help="re-stage the nominal TSVs")
    a = ap.parse_args()
    run(Experiment(Config(), a.experiment), force_extract=a.force_extract)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
