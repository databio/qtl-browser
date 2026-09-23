"""The TOPCHeF adapter: Zenodo 21382723 -> the five contract tables (pipeline/CONTRACT.md).

This module is the only place that knows what the TOPCHeF release looks like on disk: which archive
holds which result, that a nominal file is named `topchef_<chr>_MaxPC70.cis_qtl_pairs.<chr>.parquet`,
that the effect column is `slope` and not `beta`, that a splice phenotype id is
`chr:start:end:clu_<n>_<strand>:<gene>.<version>`, and that the alleles come as `A1`/`A2` with no
stated orientation. The v0 steps still call in here for those facts rather than spelling them a
second time, so there is one definition of "where the eQTL nominal rows live" in the pipeline.

Outputs go to `_tables/topchef/`: `sites.parquet`, `nominal/chr=<c>/data.parquet`,
`permuted.parquet`, `credible_sets.parquet`, `phenotypes.parquet`, and `ingestion.json` with the
orientation counts the contract asks an adapter to report.

Orientation
-----------
A2 is the reference allele and `af`/`slope` describe A1, so A1 is already ALT and the reorientation
is a relabelling: `A2 -> ref`, `A1 -> alt`, no number changes. Both halves are measured, not
assumed; `qtlstore.orient_to_ref` carries the evidence and is the one implementation of the rule.
This module supplies it with the reference allele per variant and applies its verdict.

The reference allele comes from `_tables/refcheck/<chr>.parquet`, not from re-reading the genome.
At SNPs that table is unambiguous: A2 at all 8,419,594 cis SNPs, A1 at none. At indels it is a
tie-break, not a reading -- both alleles share a leading base, so prefix matching cannot separate
them and `classify` takes the longer allele. `INDEL_REFERENCE` below says which of the two to
believe there and why, and `orient` settles each one by asking the refgetstore whether the
reference actually reads A2 at that position rather than by preferring an answer.
"""
from __future__ import annotations

import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .. import qtlstore as qs
from ..common import CHROMS, Config, connect, log, variants_sql, write_parquet
from ..steps_refget import refcheck_table

EXPERIMENT_ID = "topchef"
ALLELE_ORIENTATION_SOURCE = "topchef_refcheck_a2_is_ref"

# The two phenotype types this experiment has. The single-letter keys are the v0 `qtl_type`, which
# the raw archives, the pack kinds and the v0 tables are all keyed by; the values are the contract's
# open `phenotype_type` set. Nothing below `_tables/topchef/` sees the letters.
PHENOTYPE_TYPE = {"e": "ge", "s": "leafcutter"}
TYPES = tuple(PHENOTYPE_TYPE)

# Zenodo archive directory names. `cis_sources` is what the variant collection scans.
ARCHIVE = {
    "e": {"nominal": "cis_eQTL_nominal", "permutation": "cis_eQTL_permutation",
          "susie": "cis_eQTL_SuSiE", "trans": "trans_eQTL"},
    "s": {"nominal": "cis_sQTL_nominal", "permutation": "cis_sQTL_permutation",
          "susie": "cis_sQTL_SuSiE", "trans": "trans_sQTL"},
}
# Per-chromosome file names inside the nominal and SuSiE archives. The permutation archives hold one
# file per chromosome too, but nothing needs to name them individually, so those are globbed.
FILENAME = {
    "e": {"nominal": "topchef_{c}_MaxPC70.cis_qtl_pairs.{c}.parquet",
          "susie": "topchef_{c}_MaxPC70.SuSiE_summary.parquet"},
    "s": {"nominal": "topchefSplice_{c}_MaxPC25.cis_qtl_pairs.{c}.parquet",
          "susie": "topchefSplice_{c}_MaxPC25.SuSiE_summary.parquet"},
}

# The four intron fields and the gene a leafcutter phenotype id carries, as DuckDB expressions over
# a `phenotype_id` column. `chr:start:end:clu_<n>_<strand>:<gene>.<version>`; the strand is the last
# character of the cluster field, and the gene id is unversioned.
SPLICE_PARSE = """
    split_part(phenotype_id, ':', 1)                          AS s_chr,
    split_part(phenotype_id, ':', 2)::INTEGER                 AS intron_start,
    split_part(phenotype_id, ':', 3)::INTEGER                 AS intron_end,
    split_part(phenotype_id, ':', 4)                          AS cluster_id,
    right(split_part(phenotype_id, ':', 4), 1)                AS strand,
    split_part(split_part(phenotype_id, ':', 5), '.', 1)      AS gene_id
"""

# Which allele to call the reference at an indel the refcheck class calls `a1`. That class is a
# tie-break there, not a reading: an indel's alleles share a leading base, so prefix matching cannot
# separate them and `classify` takes the longer one.
#
#   "a2"       -- A2, the same answer as at every SNP, but only where the reference really does read
#                 A2 at the position; `orient` checks that against the refgetstore variant by
#                 variant and keeps the refcheck call where it does not.
#   "refcheck" -- the refcheck call as stored, so those indels swap.
#
# "a2", because `af` says A1 is the minor allele at those sites and a minor allele is not what
# GRCh38 reads. dbSNP does not settle it: over the 97,433 cis indels called `a1`, dbSNP puts REF on
# A1 at 41,560, on A2 at 38,034, carries both orientations at 12,894 and has no matching record at
# 4,945 -- the even split a position with both a deletion and an insertion recorded at it produces.
# `af` does settle it. `af` is A1's frequency (EVIDENCE.md A.11, analysis repo, +0.989 against the
# Jurgens DCM GWAS over 818,435 indels), and on the 41,560 where dbSNP puts REF on A1 it averages
# 0.1504 with median 0.0833 and exceeds 0.5 at only 5.3% of sites -- indistinguishable from the
# 11.7% of the undisputed `a2` indels, and nothing like the ~90% a genuine reference allele would
# show. So A1 is ALT at those sites too, the dbSNP REF=A1 records are a different variant at the
# same position, and calling A1 the reference would negate 97,433 slopes and mirror 97,433
# frequencies on the strength of a prefix tie-break. Both `G>GC` and `GC>G` are correctly anchored
# where the reference reads `GC`, which is why the check below asks only whether A2 reads, not
# which form is canonical.
INDEL_REFERENCE = "a2"

# Classes with no reference allele, so no site. `no_source_alleles` is a trans eQTL variant the
# release names only as chr:pos: excluded until the authors supply its alleles, never inferred.
NO_REFERENCE = ("neither", "unchecked", "no_source_alleles")


# ---- source paths ------------------------------------------------------------------------------
def archive(cfg: Config, qtl_type: str, kind: str) -> Path:
    return cfg.raw_dir(ARCHIVE[qtl_type][kind])


def archive_glob(cfg: Config, qtl_type: str, kind: str) -> str:
    return cfg.raw_glob(ARCHIVE[qtl_type][kind])


def source_file(cfg: Config, qtl_type: str, kind: str, chrom: str) -> Path:
    """One chromosome's raw nominal or SuSiE file."""
    return archive(cfg, qtl_type, kind) / FILENAME[qtl_type][kind].format(c=chrom)


def cis_sources() -> list[str]:
    """Every cis archive that reports alleles. These six define the cis variant set."""
    return [ARCHIVE[t][k] for k in ("nominal", "permutation", "susie") for t in TYPES]


def significance(cfg: Config) -> dict:
    """The experiment's significance rule. The contract keeps the rule out of `permuted` and in the
    experiment JSON, so the adapter reports it rather than applying it to the tables it writes."""
    col = {"pval_perm": "p_perm", "pval_beta": "p_beta"}[cfg["sig_column"]]
    return {"column": col, "op": "<", "threshold": float(cfg["sig_threshold"])}


def dof(cfg: Config) -> dict[str, int]:
    """Student-t degrees of freedom, published by this study, by contract phenotype type."""
    d = cfg["packs"]["dof"]
    return {"ge": int(d["eqtl"]), "leafcutter": int(d["sqtl"])}


# ---- output paths ------------------------------------------------------------------------------
def tables(cfg: Config) -> Path:
    return cfg.tables / EXPERIMENT_ID


def sites_path(cfg: Config) -> Path:
    return tables(cfg) / "sites.parquet"


def nominal_path(cfg: Config, chrom: str) -> Path:
    return tables(cfg) / "nominal" / f"chr={chrom}" / "data.parquet"


def permuted_path(cfg: Config) -> Path:
    return tables(cfg) / "permuted.parquet"


def credible_sets_path(cfg: Config) -> Path:
    return tables(cfg) / "credible_sets.parquet"


def phenotypes_path(cfg: Config) -> Path:
    return tables(cfg) / "phenotypes.parquet"


def trans_path(cfg: Config) -> Path:
    return tables(cfg) / "trans.parquet"


def orientation_path(cfg: Config, chrom: str) -> Path:
    """Internal: the per-variant ref/alt call the other four tables join to. Not a contract table;
    it exists so the orientation is decided once, counted once, and reported once."""
    return tables(cfg) / "_orientation" / f"chr={chrom}" / "data.parquet"


def report_path(cfg: Config) -> Path:
    return tables(cfg) / "ingestion.json"


# ---- orientation -------------------------------------------------------------------------------
def disputed(match: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """The indels the refcheck class calls `a1`: the only variants where the class is a tie-break
    and not a reading, and so the only ones `INDEL_REFERENCE` governs."""
    indel = np.array([len(x) != 1 or len(y) != 1 for x, y in zip(a1, a2)], dtype=bool)
    return (match == "a1") & indel


def reference_allele(match: np.ndarray, a1: np.ndarray, a2: np.ndarray,
                     a2_reads: np.ndarray | None = None) -> np.ndarray:
    """The reference allele string per variant, or None where there is none.

    `match` is the refcheck class (`a1`, `a2`, `both`, `neither`, `unchecked`), or
    `no_source_alleles` for a variant the release names without alleles. `both` only happens when
    A1 equals A2, so either is the reference. The `NO_REFERENCE` classes have no reference allele
    and the site is dropped by `orient_to_ref`.

    `a2_reads` says, per variant, whether the reference actually reads A2 at the position. It is
    consulted only for `disputed` variants, and only when `INDEL_REFERENCE` is `a2`: those take A2
    where the reference really does read it and keep the refcheck call where it does not. Passing
    None means the check has not been run, and then every disputed variant keeps the refcheck call.
    """
    called_a1 = match == "a1"
    if INDEL_REFERENCE == "a2" and a2_reads is not None:
        called_a1 = called_a1 & ~(disputed(match, a1, a2) & a2_reads)
    out = np.where(called_a1, a1, a2).astype(object)
    out[np.isin(match, list(NO_REFERENCE))] = None
    return out


def _a2_reads(cfg: Config, chrom: str, pos: np.ndarray, a2: np.ndarray, which: np.ndarray) -> np.ndarray:
    """Does the reference read A2 at these positions? Asked of the refgetstore, one substring per
    variant, for the `which` rows only -- about 100,000 genome-wide, against 9.2 million variants."""
    from ..steps_refget import open_store, reference_json
    out = np.zeros(len(pos), dtype=bool)
    idx = np.flatnonzero(which)
    if not len(idx):
        return out
    seqs = json.loads(reference_json(cfg).read_text())["sequences"]
    digest = seqs[chrom]["digest"]
    store = open_store(cfg)
    for i in idx:
        p, allele = int(pos[i]), a2[i].upper()
        out[i] = store.get_substring(digest, p - 1, p - 1 + len(allele)).upper() == allele
    return out


def orient(cfg: Config) -> dict:
    """Write `_orientation/chr=<c>/data.parquet` and return the ingestion counts.

    One row per (chr, position, A1, A2) with `ref`, `alt` and `swapped`. `qtlstore.orient_to_ref`
    is the one implementation of the rule; this only feeds it the reference allele and records what
    it said. The counts it returns are the `swapped` and `dropped` totals the contract's acceptance
    gate asks for, split cis against trans-only because only the cis half is claimed to be zero.

    A trans-only variant the release gives no alleles for (a trans eQTL `chr:pos` found in no cis
    and no trans sQTL file) has no `ref`/`alt` to anchor, so it gets no site. It is dropped here
    under the class `no_source_alleles`, never inferred (CONTRACT.md, `sites`).
    """
    counts = {k: {"as_is": 0, "swapped": 0, "dropped": 0} for k in ("cis", "trans_only")}
    dropped_by_class: dict[str, int] = {}
    tie = {"disputed": 0, "a2_reads": 0, "kept_refcheck_call": 0}
    for chrom in CHROMS:
        rc = refcheck_table(cfg, chrom)
        if not rc.exists():
            log(f"topchef orient: no refcheck table for {chrom}, skipping")
            continue
        out = orientation_path(cfg, chrom)
        t = pq.read_table(rc, columns=["position", "A1", "A2", "in_cis", "match"])
        a1 = np.array(t["A1"].to_pylist(), dtype=object)
        a2 = np.array(t["A2"].to_pylist(), dtype=object)
        match = np.array(t["match"].cast(pa.string()).to_pylist(), dtype=object)
        # a trans-only variant may carry no alleles at all; it has no reference allele, so no site
        missing = np.array([x is None or y is None for x, y in zip(a1, a2)], dtype=bool)
        a1[missing], a2[missing], match[missing] = "", "", "no_source_alleles"
        pos = t["position"].to_numpy().astype(np.int64)
        d = disputed(match, a1, a2)
        reads = _a2_reads(cfg, chrom, pos, a2, d) if INDEL_REFERENCE == "a2" else None
        if reads is not None:
            tie["disputed"] += int(d.sum())
            tie["a2_reads"] += int((d & reads).sum())
            tie["kept_refcheck_call"] += int((d & ~reads).sum())
        # beta = +1 going in, so a negative beta coming out is `orient_to_ref` reporting a swap.
        # Reading the verdict back beats restating the rule here as a second implementation.
        r = qs.orient_to_ref(reference_allele(match, a1, a2, reads), a1, a2, np.ones(len(a1)), np.zeros(len(a1)))
        keep = r["keep"]
        swapped = np.zeros(len(a1), dtype=bool)
        swapped[keep] = r["beta"] < 0
        in_cis = t["in_cis"].to_numpy(zero_copy_only=False)
        for label, mask in (("cis", in_cis), ("trans_only", ~in_cis)):
            counts[label]["as_is"] += int(np.sum(mask & keep & ~swapped))
            counts[label]["swapped"] += int(np.sum(mask & keep & swapped))
            counts[label]["dropped"] += int(np.sum(mask & ~keep))
        for cls in np.unique(match[~keep]):
            dropped_by_class[str(cls)] = dropped_by_class.get(str(cls), 0) + int(np.sum(~keep & (match == cls)))
        write_parquet(pa.table({
            "chr": pa.array([chrom] * int(keep.sum()), pa.string()),
            "position": t["position"].filter(pa.array(keep)).cast(pa.int32()),
            "A1": pa.array(list(a1[keep]), pa.string()), "A2": pa.array(list(a2[keep]), pa.string()),
            "ref": pa.array(list(r["ref"]), pa.string()), "alt": pa.array(list(r["alt"]), pa.string()),
            "swapped": pa.array(swapped[keep], pa.bool_()),
            "in_cis": pa.array(in_cis[keep], pa.bool_()),
        }), out, 200_000, stats_columns=["chr", "position"])
    log(f"topchef orient: cis {counts['cis']}, trans-only {counts['trans_only']}; "
        f"dropped by class {dropped_by_class}; indel tie-breaks {tie}")
    return {"counts": counts, "dropped_by_refcheck_class": dropped_by_class,
            "indel_reference": INDEL_REFERENCE, "indel_tie_break": tie,
            "allele_orientation_source": ALLELE_ORIENTATION_SOURCE}


def orientation_sql(cfg: Config, chrom: str | None = None) -> str:
    """The orientation table as a SQL source."""
    if chrom is None:
        # hive_partitioning off: the files carry their own `chr` column, and a hive column of the
        # same name would collide with it.
        return f"read_parquet('{tables(cfg) / '_orientation' / 'chr=*' / 'data.parquet'}', hive_partitioning = false)"
    return f"read_parquet('{orientation_path(cfg, chrom)}')"


def chrom_filter(column: str = "chr") -> str:
    """`chr IN (...)` over the build's chromosome list. A genome-wide build names all of them and
    this changes nothing; a `QTLB_CHROMS` subset needs it, because the permutation and SuSiE
    archives are read whole and their other chromosomes have no orientation table to join to."""
    return f"{column} IN (" + ", ".join(f"'{c}'" for c in CHROMS) + ")"


def _swap(col: str, swapped: str = "o.swapped") -> str:
    """Apply an orientation swap to a value that describes the effect allele. Kept next to
    `orient_to_ref`'s rule in one place: `beta` negates, `af` mirrors, `se` does neither."""
    return f"CASE WHEN {swapped} THEN -{col} ELSE {col} END"


def _swap_af(col: str, swapped: str = "o.swapped") -> str:
    return f"CASE WHEN {swapped} THEN 1.0 - {col} ELSE {col} END"


# ---- sites -------------------------------------------------------------------------------------
def sites(cfg: Config) -> None:
    """`sites.parquet`: every distinct variant this experiment tested, once, in (chr, pos) order.

    `af`, `ma_samples` and `ma_count` come from the raw nominal rows, one value per variant, the
    eQTL value where the variant was tested for eQTL and the sQTL value otherwise (SPEC section 4).
    Trans-only variants (all from trans sQTL, the only trans file with alleles) take their `af` from
    the trans rows and have no minor-allele counts. `rsid`, `rs_number` and `match` ride along as
    variant catalog attributes, not contract columns.
    """
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    parts = []
    for chrom in CHROMS:
        if not orientation_path(cfg, chrom).exists():
            continue
        e = source_file(cfg, "e", "nominal", chrom)
        s = source_file(cfg, "s", "nominal", chrom)
        if not (e.exists() and s.exists()):
            log(f"topchef sites: missing raw nominal file for {chrom}, skipping")
            continue
        for name, src in (("ve", e), ("vs", s)):
            con.execute(f"""CREATE OR REPLACE TABLE {name} AS
                SELECT position, A1, A2, min(af)::FLOAT AS af,
                       min(ma_samples)::INTEGER AS ma_samples, min(ma_count)::INTEGER AS ma_count,
                       count(DISTINCT af) > 1 OR count(DISTINCT ma_samples) > 1 OR count(DISTINCT ma_count) > 1 AS differs
                FROM '{src}' GROUP BY 1, 2, 3""")
        bad = con.execute("SELECT (SELECT count(*) FILTER (WHERE differs) FROM ve), (SELECT count(*) FILTER (WHERE differs) FROM vs)").fetchone()
        if any(bad):
            raise ValueError(f"topchef sites {chrom}: {bad[0]} eQTL and {bad[1]} sQTL variants report more than one af or count")
        parts.append(con.execute(f"""
            SELECT o.chr, o.position AS pos, o.ref, o.alt,
                   {_swap_af('coalesce(ve.af, vs.af, tr.af)')}::FLOAT AS af,
                   coalesce(ve.ma_samples, vs.ma_samples, -1)::INTEGER AS ma_samples,
                   coalesce(ve.ma_count, vs.ma_count, -1)::INTEGER AS ma_count,
                   o.in_cis, v.rsid, v.rs_number, v.match
            FROM {orientation_sql(cfg, chrom)} o
            LEFT JOIN ve ON ve.position = o.position AND ve.A1 = o.A1 AND ve.A2 = o.A2
            LEFT JOIN vs ON vs.position = o.position AND vs.A1 = o.A1 AND vs.A2 = o.A2
            LEFT JOIN '{_trans_af_table(cfg, con)}' tr ON tr.variant_chr = o.chr AND tr.position = o.position
            LEFT JOIN {variants_sql(cfg, chrom)} v ON v.position = o.position
              AND v.A1 IS NOT DISTINCT FROM o.A1 AND v.A2 IS NOT DISTINCT FROM o.A2
            ORDER BY o.position, o.ref, o.alt
        """).fetch_arrow_table())
    t = pa.concat_tables(parts) if parts else _empty_sites()
    write_parquet(t, sites_path(cfg), 200_000, stats_columns=["chr", "pos"])
    n_cis = int(pc.sum(pc.cast(t["in_cis"], pa.int64())).as_py() or 0)
    log(f"topchef sites: {t.num_rows:,} sites ({n_cis:,} cis) -> {sites_path(cfg).relative_to(cfg.derived)}")


def _empty_sites() -> pa.Table:
    return pa.table({"chr": pa.array([], pa.string()), "pos": pa.array([], pa.int32()),
                     "ref": pa.array([], pa.string()), "alt": pa.array([], pa.string()),
                     "af": pa.array([], pa.float32()), "ma_samples": pa.array([], pa.int32()),
                     "ma_count": pa.array([], pa.int32()), "in_cis": pa.array([], pa.bool_()),
                     "rsid": pa.array([], pa.string()), "rs_number": pa.array([], pa.int64()),
                     "match": pa.array([], pa.string())})


def _trans_af_table(cfg: Config, con) -> Path:
    """One af per trans variant, from the trans rows. Built once and cached under `_tables/topchef/`
    because both `sites` and the trans results set want it and it costs a full scan of both trans
    archives. The trans eQTL file names its variant only as `chr:pos`, so this is keyed by position."""
    out = tables(cfg) / "_trans_variant_af.parquet"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name("_trans_variant_af.tmp.parquet")   # keep .parquet: DuckDB reads by extension
    union = " UNION ALL ".join(
        f"SELECT split_part(variant_id, ':', 1) AS variant_chr, split_part(variant_id, ':', 2)::INTEGER AS position, af "
        f"FROM read_parquet('{archive_glob(cfg, t, 'trans')}')" for t in TYPES)
    con.execute(f"""COPY (
        SELECT variant_chr, position, min(af)::FLOAT AS af, count(DISTINCT af) AS n_af
        FROM ({union}) GROUP BY 1, 2
    ) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD)""")
    multi = con.execute(f"SELECT count(*) FILTER (WHERE n_af <> 1) FROM '{tmp}'").fetchone()[0]
    if multi:
        tmp.unlink()
        raise ValueError(f"topchef sites: {multi:,} trans variants report more than one af over their trans rows")
    os.replace(tmp, out)
    return out


# ---- phenotypes and permuted -------------------------------------------------------------------
def _permutation_view(cfg: Config, con, qtl_type: str, name: str) -> None:
    """The permutation rows of one phenotype type, with the leafcutter fields parsed out."""
    extra = f", {SPLICE_PARSE}" if qtl_type == "s" else ", phenotype_id AS gene_id"
    con.execute(f"CREATE OR REPLACE TABLE {name} AS SELECT *{extra} "
                f"FROM read_parquet('{archive_glob(cfg, qtl_type, 'permutation')}') WHERE {chrom_filter()}")


def _nominal_counts(cfg: Config, con) -> Path:
    """Nominal rows per (phenotype_type, phenotype_id), read back from the contract nominal files.
    `has_nominal` and the `n_variants` recount both come from here, so they describe the table a
    reader actually gets and not the raw archives. Cached next to the tables and rebuilt whenever a
    nominal file is newer than the cache.

    `n_nominal` is every row; `n_tested` leaves out rows with no p-value. Those are variants with no
    genotype variance in the samples (29 genome-wide, every sample heterozygous: af 0.5,
    ma_samples 516): tensorQTL's nominal pass writes them with null statistics, but its permutation
    pass drops monomorphic variants before counting `num_var`. They were never tested, so
    `n_variants` counts `n_tested`, which then equals the source `num_var` everywhere."""
    out = tables(cfg) / "_nominal_counts.parquet"
    srcs = [nominal_path(cfg, c) for c in CHROMS if nominal_path(cfg, c).exists()]
    if not srcs:
        raise FileNotFoundError("topchef: no nominal tables yet; run the nominal step first")
    if (out.exists() and out.stat().st_mtime > max(p.stat().st_mtime for p in srcs)
            and "n_tested" in pq.read_schema(out).names):
        return out
    files = ", ".join(f"'{p}'" for p in srcs)
    tmp = out.with_name("_nominal_counts.tmp.parquet")
    con.execute(f"""COPY (
        SELECT phenotype_type, phenotype_id, count(*)::INTEGER AS n_nominal,
               count(*) FILTER (WHERE pvalue IS NOT NULL AND NOT isnan(pvalue))::INTEGER AS n_tested
        FROM read_parquet([{files}], hive_partitioning = false) GROUP BY 1, 2
    ) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD)""")
    os.replace(tmp, out)
    return out


def phenotypes(cfg: Config) -> None:
    """`phenotypes.parquet`: what each phenotype is, with no annotation columns. The leafcutter
    intron fields go to `extra` as JSON, which is why a third phenotype type needs no new file.

    tensorQTL permutes every phenotype on its own, so each is its own group:
    `phenotype_object_id = phenotype_id` throughout. `has_nominal` is read from the contract
    nominal tables, so the nominal step must have run. Every tested TOPCHeF phenotype has nominal
    rows; the column is there for sources that publish nominal rows for only some phenotypes.
    """
    con = connect(cfg)
    counts = _nominal_counts(cfg, con)
    parts = []
    for qtl_type in TYPES:
        _permutation_view(cfg, con, qtl_type, "p")
        extra = ("'{}'" if qtl_type == "e" else
                 "to_json({'intron_start': intron_start, 'intron_end': intron_end, "
                 "'cluster_id': cluster_id, 'strand': strand})::VARCHAR")
        pt = PHENOTYPE_TYPE[qtl_type]
        parts.append(con.execute(f"""
            SELECT '{pt}' AS phenotype_type, p.phenotype_id, p.phenotype_id AS phenotype_object_id,
                   p.gene_id, n.phenotype_id IS NOT NULL AS has_nominal, {extra} AS extra
            FROM p LEFT JOIN (SELECT * FROM '{counts}' WHERE phenotype_type = '{pt}') n USING (phenotype_id)
            ORDER BY p.phenotype_id""").fetch_arrow_table())
        # phenotypes the trans file names that have no cis result (17 genes, most on chrM, and 36 introns
        # genome-wide): the source describes them, so they get a row, with no nominal rows and no group
        gene = "phenotype_id" if qtl_type == "e" else "split_part(split_part(phenotype_id, ':', 5), '.', 1)"
        tx = ("'{}'" if qtl_type == "e" else
              "to_json({'intron_start': split_part(phenotype_id, ':', 2)::INTEGER, "
              "'intron_end': split_part(phenotype_id, ':', 3)::INTEGER, "
              "'cluster_id': split_part(phenotype_id, ':', 4), "
              "'strand': right(split_part(phenotype_id, ':', 4), 1)})::VARCHAR")
        parts.append(con.execute(f"""
            SELECT DISTINCT '{pt}' AS phenotype_type, phenotype_id, phenotype_id AS phenotype_object_id,
                   {gene} AS gene_id, false AS has_nominal, {tx} AS extra
            FROM read_parquet('{archive_glob(cfg, qtl_type, 'trans')}')
            WHERE phenotype_id NOT IN (SELECT phenotype_id FROM p)
            ORDER BY phenotype_id""").fetch_arrow_table().cast(parts[-1].schema))
    t = pa.concat_tables(parts)
    write_parquet(t, phenotypes_path(cfg), 10_000, stats_columns=["phenotype_type", "phenotype_id", "gene_id"])
    log(f"topchef phenotypes: {t.num_rows:,} rows -> {phenotypes_path(cfg).relative_to(cfg.derived)}")


def permuted(cfg: Config) -> None:
    """`permuted.parquet`: one row per phenotype group, with the lead variant in ref/alt orientation.

    tensorQTL permutes each phenotype on its own, so a TOPCHeF group is one phenotype
    (`phenotype_object_id = phenotype_id`) and the key `(phenotype_type, phenotype_object_id)` is
    also one row per phenotype. `n_variants` is recounted from the contract nominal table wherever
    the phenotype has nominal rows (complete for that phenotype), counting only rows with a p-value
    (see `_nominal_counts`), and copied from the source `num_var` otherwise; `coverage` reports the
    split and the disagreements.
    """
    con = connect(cfg)
    counts = _nominal_counts(cfg, con)
    parts = []
    for qtl_type in TYPES:
        _permutation_view(cfg, con, qtl_type, "p")
        n, distinct = con.execute("SELECT count(*), count(DISTINCT phenotype_id) FROM p").fetchone()
        if n != distinct:
            raise ValueError(f"topchef permuted {qtl_type}: {n:,} rows for {distinct:,} phenotypes; "
                             "the contract requires exactly one row per phenotype")
        parts.append(con.execute(f"""
            SELECT '{PHENOTYPE_TYPE[qtl_type]}' AS phenotype_type,
                   p.phenotype_id AS phenotype_object_id, p.phenotype_id, p.gene_id,
                   coalesce(n.n_tested, p.num_var)::INTEGER AS n_variants,
                   p.pval_perm AS p_perm, p.pval_beta AS p_beta,
                   p.chr AS lead_chr, p.position::INTEGER AS lead_pos, o.ref AS lead_ref, o.alt AS lead_alt
            FROM p LEFT JOIN {orientation_sql(cfg)} o
              ON o.chr = p.chr AND o.position = p.position AND o.A1 = p.A1 AND o.A2 = p.A2
            LEFT JOIN (SELECT * FROM '{counts}' WHERE phenotype_type = '{PHENOTYPE_TYPE[qtl_type]}') n
              ON n.phenotype_id = p.phenotype_id
            ORDER BY p.phenotype_id""").fetch_arrow_table())
    t = pa.concat_tables(parts)
    orphan = int(pc.sum(pc.cast(pc.is_null(t["lead_ref"]), pa.int64())).as_py() or 0)
    if orphan:
        raise ValueError(f"topchef permuted: {orphan:,} lead variants are not in sites")
    write_parquet(t, permuted_path(cfg), 10_000, stats_columns=["phenotype_type", "phenotype_id", "gene_id"])
    log(f"topchef permuted: {t.num_rows:,} rows -> {permuted_path(cfg).relative_to(cfg.derived)}")


# ---- credible sets -----------------------------------------------------------------------------
def credible_sets(cfg: Config) -> None:
    """`credible_sets.parquet`: one row per (phenotype, credible set, variant).

    A variant in two sets of one phenotype keeps both rows; that is real and the hits pack shows it.
    SuSiE here reports neither a z-score nor a within-set r2, so `z` and `cs_min_r2` are NaN;
    `cs_size` is counted from the memberships.
    """
    con = connect(cfg)
    parts = []
    for qtl_type in TYPES:
        gene = ("phenotype_id" if qtl_type == "e"
                else "split_part(split_part(phenotype_id, ':', 5), '.', 1)")
        con.execute(f"""CREATE OR REPLACE TABLE cs AS
            SELECT phenotype_id, {gene} AS gene_id, chr, position, A1, A2, cs_id::SMALLINT AS cs_id, pip
            FROM read_parquet('{archive_glob(cfg, qtl_type, 'susie')}') WHERE {chrom_filter()}""")
        parts.append(con.execute(f"""
            SELECT '{PHENOTYPE_TYPE[qtl_type]}' AS phenotype_type,
                   cs.phenotype_id AS phenotype_object_id, cs.phenotype_id, cs.cs_id,
                   cs.chr, cs.position::INTEGER AS pos, o.ref, o.alt, cs.pip::FLOAT AS pip,
                   'nan'::FLOAT AS z,
                   count(*) OVER (PARTITION BY cs.phenotype_id, cs.cs_id)::INTEGER AS cs_size,
                   'nan'::FLOAT AS cs_min_r2
            FROM cs LEFT JOIN {orientation_sql(cfg)} o
              ON o.chr = cs.chr AND o.position = cs.position AND o.A1 = cs.A1 AND o.A2 = cs.A2
            ORDER BY cs.phenotype_id, cs.cs_id, cs.pip DESC""").fetch_arrow_table())
    t = pa.concat_tables(parts)
    orphan = int(pc.sum(pc.cast(pc.is_null(t["ref"]), pa.int64())).as_py() or 0)
    if orphan:
        raise ValueError(f"topchef credible_sets: {orphan:,} variants are not in sites")
    write_parquet(t, credible_sets_path(cfg), 2_000, stats_columns=["phenotype_type", "phenotype_id", "chr", "pos"])
    log(f"topchef credible_sets: {t.num_rows:,} rows -> {credible_sets_path(cfg).relative_to(cfg.derived)}")


# ---- nominal -----------------------------------------------------------------------------------
# One row group per phenotype: every reader of this table filters on the phenotype, and a reader
# that pulls a whole row group should not get its neighbours too. Same reasoning, and the same
# encodings, as the v0 build intermediates (`steps_nominal`).
NOMINAL_ENCODING = {"pos": "DELTA_BINARY_PACKED",
                    **{c: "BYTE_STREAM_SPLIT" for c in ("beta", "se", "pvalue")}}
NOMINAL_STATS = ["phenotype_type", "phenotype_id", "gene_id", "pos"]


# Rows per sorted chunk. A chromosome's nominal rows are sorted a run of genes at a time and the
# runs written back to back, which is the same order as one global sort (a phenotype belongs to one
# gene, and the runs are in gene order). chr1 holds 58 million rows once every intron keeps its
# rows; sorted in one piece that overran a 6 GB DuckDB.
NOMINAL_CHUNK_ROWS = 8_000_000


def _gene_chunks(con, src: Path, gene: str, limit: int) -> list[tuple[str, str]]:
    """Consecutive (first, last) gene id ranges of `src`, each holding at most `limit` rows unless
    one gene alone holds more."""
    chunks: list[tuple[str, str]] = []
    lo = prev = None
    n = 0
    for g, c in con.execute(f"SELECT {gene} AS g, count(*) FROM '{src}' n GROUP BY 1 ORDER BY 1").fetchall():
        if lo is not None and n + c > limit:
            chunks.append((lo, prev))
            lo, n = None, 0
        lo = g if lo is None else lo
        n += c
        prev = g
    if lo is not None:
        chunks.append((lo, prev))
    return chunks


def _nominal_chrom(args) -> tuple[str, int, int]:
    """One chromosome's contract nominal file: both phenotype types, one row group per phenotype."""
    chrom, out = args[0], Path(args[1])
    cfg = Config()
    work = cfg.tmp / f"topchef-nominal-{chrom}"
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"], temp_dir=work)
    con.execute("SET preserve_insertion_order = true")
    con.execute(f"CREATE TABLE o AS SELECT position, A1, A2, ref, alt, swapped FROM {orientation_sql(cfg, chrom)} WHERE in_cis")
    out.parent.mkdir(parents=True, exist_ok=True)
    parts: list[Path] = []
    raw_total = 0
    # `TYPES` is in `phenotype_type` order (ge < leafcutter), so the parts come out in the contract
    # table's order: phenotype_type, gene_id, phenotype_id, pos
    for qtl_type in TYPES:
        src = source_file(cfg, qtl_type, "nominal", chrom)
        if not src.exists():
            continue
        gene = ("n.phenotype_id" if qtl_type == "e"
                else "split_part(split_part(n.phenotype_id, ':', 5), '.', 1)")
        # Every tested phenotype keeps every row, leafcutter included. `sqtl_nominal` governs only the
        # v0 `_tables/cis_sqtl_nominal` intermediate; the v0 sQTL packs stream this same raw file
        # (`steps_pack._raw_sqtl`) for every tested intron, so a contract table trimmed to the
        # significant introns would leave most introns with `has_nominal = false` that v0 serves.
        for lo, hi in _gene_chunks(con, src, gene, NOMINAL_CHUNK_ROWS):
            part = out.with_name(f"sorted.{len(parts):03d}.tmp.parquet")
            con.execute(f"""COPY (
                SELECT '{PHENOTYPE_TYPE[qtl_type]}' AS phenotype_type, n.phenotype_id, {gene} AS gene_id,
                       n.chr, n.position::INTEGER AS pos, o.ref, o.alt,
                       {_swap('n.slope')}::FLOAT AS beta, n.slope_se::FLOAT AS se, n.pval_nominal AS pvalue
                FROM (SELECT * FROM '{src}' n WHERE {gene} BETWEEN ? AND ?) n
                JOIN o ON o.position = n.position AND o.A1 = n.A1 AND o.A2 = n.A2
                ORDER BY gene_id, n.phenotype_id, pos
            ) TO '{part}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000)""", [lo, hi])
            parts.append(part)
        raw_total += pq.read_metadata(src).num_rows
    con.close()
    shutil.rmtree(work, ignore_errors=True)
    if not parts:
        return chrom, 0, 0
    n = sum(pq.read_metadata(p).num_rows for p in parts)
    if n != raw_total:
        for p in parts:
            p.unlink()
        raise RuntimeError(f"topchef nominal {chrom}: {raw_total:,} raw rows became {n:,}; "
                           "a row either lost its site or the orientation join duplicated it")
    groups = _regroup_by_phenotype(parts, out)
    for p in parts:
        p.unlink()
    return chrom, n, groups


def _regroup_by_phenotype(srcs: Path | list[Path], out: Path) -> int:
    """Stream `srcs` (one file, or several read one after another, already ordered by phenotype)
    into `out` with one row group per phenotype. Nothing holds a whole chromosome: batches come in,
    whole phenotypes go out."""
    srcs = [srcs] if isinstance(srcs, Path) else list(srcs)
    schema = pq.read_schema(srcs[0])
    writer = pq.ParquetWriter(out, schema, compression="zstd", compression_level=9,
                              use_dictionary=[c for c in schema.names if c not in NOMINAL_ENCODING],
                              column_encoding=NOMINAL_ENCODING, write_statistics=NOMINAL_STATS)
    buf: list[pa.Table] = []
    cur = None
    n = 0

    def flush():
        nonlocal buf, n
        if buf:
            t = pa.concat_tables(buf)
            writer.write_table(t, row_group_size=t.num_rows)
            n += 1
            buf = []

    try:
        for batch in (b for src in srcs for b in pq.ParquetFile(src).iter_batches(batch_size=250_000)):
            t = pa.Table.from_batches([batch])
            col = t["phenotype_id"].combine_chunks()
            if len(col) == 0:
                continue
            idx = pc.indices_nonzero(pc.not_equal(col.slice(0, len(col) - 1), col.slice(1))).to_pylist()
            starts = [0] + [i + 1 for i in idx]
            for s, e in zip(starts, starts[1:] + [t.num_rows]):
                g = col[s].as_py()
                if g != cur:
                    flush()
                    cur = g
                buf.append(t.slice(s, e - s))
        flush()
    finally:
        writer.close()
    return n


def nominal(cfg: Config, force: bool = False) -> None:
    """`nominal/chr=<c>/data.parquet`: one row per (phenotype, variant) tested, in ref/alt
    orientation. `beta` is the ALT effect, `se` is untouched, `pvalue` is the source's."""
    jobs = []
    for chrom in CHROMS:
        if not orientation_path(cfg, chrom).exists():
            continue
        srcs = [source_file(cfg, t, "nominal", chrom) for t in TYPES]
        srcs = [p for p in srcs if p.exists()]
        if not srcs:
            log(f"topchef nominal: no raw nominal file for {chrom}, skipping")
            continue
        out = nominal_path(cfg, chrom)
        if not force and out.exists() and out.stat().st_mtime > max(p.stat().st_mtime for p in srcs):
            continue
        jobs.append((chrom, str(out), sum(p.stat().st_size for p in srcs)))
    jobs.sort(key=lambda j: -j[2])                     # biggest first
    log(f"topchef nominal: {len(jobs)} chromosomes with {cfg['workers']} workers")
    total = groups = 0
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as ex:
        futs = {ex.submit(_nominal_chrom, j): j for j in jobs}
        for f in as_completed(futs):
            chrom, n, g = f.result()
            total, groups = total + n, groups + g
            log(f"topchef nominal: {chrom}: {n:,} rows in {g:,} row groups")
    log(f"topchef nominal: {total:,} rows in {groups:,} row groups")


# ---- trans -------------------------------------------------------------------------------------
def trans(cfg: Config) -> dict:
    """`trans.parquet`: one row per trans (phenotype, variant) pair whose variant has source alleles, in
    ref/alt orientation, sorted by (phenotype_type, phenotype_id).

    Trans sQTL rows carry A1/A2 and join the orientation table on (chr, position, A1, A2). Trans eQTL
    rows name the variant as `chr:pos` only; they join on (chr, position), which is one variant per
    position in this release (the authors' plink2 `--rm-dup force-first` kept one tested allele per
    position, and the orientation table is checked to hold one row per position here). The rows whose
    position has no orientation row are the allele-less variants: left out, never inferred, and their
    count must equal `trans_eqtl_excluded.trans_eqtl_rows`. `beta` negates where the orientation
    swapped; `se` and `pvalue` are the source's."""
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    o = orientation_sql(cfg)
    multi = con.execute(f"SELECT count(*) FROM (SELECT chr, position FROM {o} GROUP BY 1, 2 HAVING count(*) > 1)").fetchone()[0]
    if multi:
        raise ValueError(f"topchef trans: {multi:,} positions carry more than one variant; a chr:pos trans eQTL row "
                         "would be ambiguous")
    e = f"""SELECT 'ge' AS phenotype_type, t.phenotype_id, t.phenotype_id AS gene_id, o.chr, o.position::INTEGER AS pos,
               o.ref, o.alt, ({_swap('t.b')})::FLOAT AS beta, t.b_se::FLOAT AS se, t.pval AS pvalue
        FROM read_parquet('{archive_glob(cfg, 'e', 'trans')}') t
        JOIN {o} o ON o.chr = split_part(t.variant_id, ':', 1) AND o.position = split_part(t.variant_id, ':', 2)::INTEGER
        WHERE {chrom_filter('split_part(t.variant_id, chr(58), 1)')}"""
    s_ = f"""SELECT 'leafcutter' AS phenotype_type, t.phenotype_id,
               split_part(split_part(t.phenotype_id, ':', 5), '.', 1) AS gene_id, o.chr, o.position::INTEGER AS pos,
               o.ref, o.alt, ({_swap('t.b')})::FLOAT AS beta, t.b_se::FLOAT AS se, t.pval AS pvalue
        FROM read_parquet('{archive_glob(cfg, 's', 'trans')}') t
        LEFT JOIN {o} o ON o.chr = t.chr AND o.position = t.position AND o.A1 = t.A1 AND o.A2 = t.A2
        WHERE {chrom_filter('t.chr')}"""
    lost_s = con.execute(f"SELECT count(*) FROM ({s_}) WHERE ref IS NULL").fetchone()[0]
    if lost_s:
        raise ValueError(f"topchef trans: {lost_s:,} trans sQTL rows have no oriented variant")
    tmp = trans_path(cfg).with_name("trans.tmp.parquet")
    trans_path(cfg).parent.mkdir(parents=True, exist_ok=True)
    con.execute("SET preserve_insertion_order = true")
    con.execute(f"""COPY (SELECT * FROM ({e} UNION ALL {s_}) ORDER BY phenotype_type, phenotype_id, chr, pos)
        TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 250000)""")
    os.replace(tmp, trans_path(cfg))
    by = dict(con.execute(f"SELECT phenotype_type, count(*) FROM '{trans_path(cfg)}' GROUP BY 1").fetchall())
    src = {t_: con.execute(f"SELECT count(*) FROM read_parquet('{archive_glob(cfg, t_, 'trans')}') t "
                           f"WHERE {chrom_filter('split_part(t.variant_id, chr(58), 1)')}").fetchone()[0] for t_ in TYPES}
    rep = {"rows": int(sum(by.values())), "by_type": {k: int(v) for k, v in by.items()},
           "source_rows": {PHENOTYPE_TYPE[k]: int(v) for k, v in src.items()},
           "ge_rows_without_source_alleles": int(src["e"] - by.get("ge", 0))}
    log(f"topchef trans: {rep}")
    return rep


# ---- the whole adapter -------------------------------------------------------------------------
def run(cfg: Config, force: bool = False) -> None:
    """The five contract tables plus the ingestion report."""
    t0 = time.time()
    report = orient(cfg)
    sites(cfg)
    nominal(cfg, force=force)          # before phenotypes/permuted: has_nominal and n_variants read it
    phenotypes(cfg)
    permuted(cfg)
    credible_sets(cfg)
    trans_rep = trans(cfg)
    excluded = trans_excluded(cfg)
    if trans_rep["ge_rows_without_source_alleles"] != excluded["trans_eqtl_rows"]:
        raise ValueError(f"topchef trans: {trans_rep['ge_rows_without_source_alleles']:,} trans eQTL rows left out, "
                         f"but {excluded['trans_eqtl_rows']:,} name an allele-less variant")
    st = pq.read_metadata(sites_path(cfg))
    report |= {
        "trans": trans_rep,
        "trans_eqtl_excluded": excluded,
        "phenotype_coverage": coverage(cfg),
        "experiment_id": EXPERIMENT_ID,
        "source": {"zenodo": cfg["zenodo_dir"], "archives": sorted({a for t in TYPES for a in ARCHIVE[t].values()})},
        "phenotype_types": sorted(PHENOTYPE_TYPE.values()),
        "dof": dof(cfg),
        "significance": significance(cfg),
        "sqtl_nominal": "all",          # every tested intron; `cfg["sqtl_nominal"]` is the v0 intermediate's
        "rows": {"sites": st.num_rows,
                 "permuted": pq.read_metadata(permuted_path(cfg)).num_rows,
                 "credible_sets": pq.read_metadata(credible_sets_path(cfg)).num_rows,
                 "phenotypes": pq.read_metadata(phenotypes_path(cfg)).num_rows,
                 "trans": pq.read_metadata(trans_path(cfg)).num_rows,
                 "nominal": sum(pq.read_metadata(nominal_path(cfg, c)).num_rows
                                for c in CHROMS if nominal_path(cfg, c).exists())},
        "chromosomes": [c for c in CHROMS if orientation_path(cfg, c).exists()],
    }
    report_path(cfg).parent.mkdir(parents=True, exist_ok=True)
    report_path(cfg).write_text(json.dumps(report, indent=2) + "\n")
    log(f"topchef: {report['rows']} in {(time.time() - t0) / 60:.1f} min -> {report_path(cfg).relative_to(cfg.derived)}")


TRANS_EXCLUSION_REASON = (
    "the trans eQTL file names these variants only as chr:pos, and no cis or trans sQTL file gives "
    "their alleles; without REF/ALT a variant cannot be anchored, so they are left out until the "
    "authors supply the alleles. Alleles are not inferred (CONTRACT.md, sites).")


def trans_excluded(cfg: Config) -> dict:
    """What the no-allele rule leaves out: the trans-only variants with no source alleles, and the
    trans eQTL rows (variant, gene pairs) that name them. Trans sQTL rows carry A1/A2 and lose
    nothing. Counted from the source, so the report says how much trans eQTL is missing, not just
    how many sites."""
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    con.execute(f"""CREATE TABLE na AS SELECT DISTINCT chr, position FROM {variants_sql(cfg)}
                    WHERE (A1 IS NULL OR A2 IS NULL) AND NOT in_cis AND {chrom_filter()}""")
    con.execute(f"""CREATE TABLE te AS
        SELECT split_part(variant_id, ':', 1) AS chr, split_part(variant_id, ':', 2)::INTEGER AS position,
               phenotype_id
        FROM read_parquet('{archive_glob(cfg, 'e', 'trans')}')""")
    variants, rows, genes, total = con.execute(f"""
        SELECT (SELECT count(*) FROM na),
               count(*) FILTER (WHERE na.position IS NOT NULL),
               count(DISTINCT phenotype_id) FILTER (WHERE na.position IS NOT NULL),
               count(*)
        FROM te LEFT JOIN na ON na.chr = te.chr AND na.position = te.position
        WHERE {chrom_filter('te.chr')}""").fetchone()
    log(f"topchef trans: excluded {variants:,} allele-less variants and {rows:,} of {total:,} trans eQTL rows")
    return {"variants": int(variants), "trans_eqtl_rows": int(rows), "trans_eqtl_rows_total": int(total),
            "genes_with_excluded_rows": int(genes), "trans_sqtl_rows": 0,
            "reason": TRANS_EXCLUSION_REASON}


def coverage(cfg: Config) -> dict:
    """Per phenotype type: phenotypes, groups, how many have nominal rows, and where each
    `permuted.n_variants` came from (CONTRACT.md, the ingestion report)."""
    con = connect(cfg)
    counts = _nominal_counts(cfg, con)
    out = {}
    for pt, n, groups, with_nom in con.execute(f"""
            SELECT phenotype_type, count(*), count(DISTINCT phenotype_object_id), count(*) FILTER (WHERE has_nominal)
            FROM '{phenotypes_path(cfg)}' GROUP BY 1 ORDER BY 1""").fetchall():
        out[pt] = {"phenotypes": n, "groups": groups, "with_nominal": with_nom}
    for pt, rec, cop, differs, untested in con.execute(f"""
            SELECT p.phenotype_type, count(n.n_tested), count(*) - count(n.n_tested),
                   count(*) FILTER (WHERE n.n_tested IS NOT NULL AND n.n_tested <> s.num_var),
                   coalesce(sum(n.n_nominal - n.n_tested), 0)
            FROM '{permuted_path(cfg)}' p
            LEFT JOIN '{counts}' n USING (phenotype_type, phenotype_id)
            LEFT JOIN ({" UNION ALL ".join(
                f"SELECT '{PHENOTYPE_TYPE[t]}' AS phenotype_type, phenotype_id, num_var "
                f"FROM read_parquet('{archive_glob(cfg, t, 'permutation')}')" for t in TYPES)}) s
              USING (phenotype_type, phenotype_id)
            GROUP BY 1 ORDER BY 1""").fetchall():
        out.setdefault(pt, {})["n_variants"] = {"recounted_from_nominal": rec, "copied_from_source": cop,
                                                "recount_differs_from_source": differs,
                                                "nominal_rows_without_pvalue": int(untested)}
    return out


def _trans_step(cfg: Config) -> None:
    """`--step trans`: the trans table alone, with its counts merged into an existing ingestion report."""
    rep = trans(cfg)
    excluded = trans_excluded(cfg)
    if rep["ge_rows_without_source_alleles"] != excluded["trans_eqtl_rows"]:
        raise ValueError(f"topchef trans: {rep['ge_rows_without_source_alleles']:,} trans eQTL rows left out, "
                         f"but {excluded['trans_eqtl_rows']:,} name an allele-less variant")
    if report_path(cfg).exists():
        r = json.loads(report_path(cfg).read_text())
        r["trans"], r["trans_eqtl_excluded"] = rep, excluded
        r.setdefault("rows", {})["trans"] = pq.read_metadata(trans_path(cfg)).num_rows
        # the trans-only phenotypes join `phenotypes` (run `--step phenotypes` first), so its count moves
        r["rows"]["phenotypes"] = pq.read_metadata(phenotypes_path(cfg)).num_rows
        r["phenotype_coverage"] = coverage(cfg)
        report_path(cfg).write_text(json.dumps(r, indent=2) + "\n")


def main() -> int:
    """`uv run python -m pipeline.adapters.topchef [--force] [--step orient|sites|...]`.

    Not in the `pipeline build` step list: the contract tables are the v1 ingestion, and the v0
    build has its own steps that write the v0 tables. Both read the same raw archives.
    """
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--step", action="append", choices=["orient", "sites", "phenotypes", "permuted",
                                                        "credible_sets", "nominal", "trans"],
                    help="run only this table (repeatable; default all, plus the ingestion report)")
    ap.add_argument("--force", action="store_true", help="rebuild files that are already current")
    a = ap.parse_args()
    cfg = Config()
    if not a.step:
        run(cfg, force=a.force)
        return 0
    fns = {"orient": lambda c: orient(c), "sites": sites, "phenotypes": phenotypes,
           "permuted": permuted, "credible_sets": credible_sets,
           "nominal": lambda c: nominal(c, force=a.force), "trans": _trans_step}
    for name in a.step:
        t0 = time.time()
        log(f"== topchef {name}")
        fns[name](cfg)
        log(f"== topchef {name} finished in {(time.time() - t0) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
