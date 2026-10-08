"""The ARIC plasma cis-pQTL adapter: PLINK2 `--glm` output -> the five contract tables.

    QTLB_CHROMS=chr21,chr22 QTLB_DERIVED=/scratch/$USER/qtl-browser/derived-aric \
        uv run python -m pipeline.adapters.aric_pqtl [--experiment aric_plasma_ea]

Zhang, Dutta, Chatterjee et al. 2022 Nat Genet (10.1038/s41588-022-01975-5 for the DCM GWAS; this
study is 10.1038/s41588-022-01051-w), from nilanjanchatterjeelab.org/pwas. One experiment per
cohort, named in `config.yaml` under `aric_pqtl:`, so a second cohort is a config entry.

What the source looks like (measured over all 4,657 proteins of both cohorts, 2026-10-07):

- Two zips, `EA.zip` and `AA.zip`, each holding 4,659 members: one PLINK2 `.glm.linear` per SOMAmer
  named `<COHORT>/SeqId_<n>_<m>.PHENO1.glm.linear`, plus the directory entry. The members are read
  straight out of the zip; nothing is unpacked.
- Columns `#CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT BETA SE T_STAT P ERRCODE`. `TEST` is `ADD`
  and `ERRCODE` is `.` on every row of both cohorts, and no cell is NA, so there is nothing to
  filter. `OBS_CT` is exactly 7,213 on all 12,850,239 EA rows and 1,871 on all 22,983,529 AA rows,
  so the cohort's N is a fact rather than a fit.
- **`BETA` is relative to `A1`, and `A1` is the minor allele, not `ALT`.** A1 is REF on 18.5% of EA
  rows and 11.8% of AA rows. `qtlstore.orient_to_ref` handles it: the effect allele is A1, the other
  allele is whichever of REF/ALT is not A1, and it negates `beta` and mirrors `af` wherever the
  effect allele turns out to be the reference. Skipping that step would sign-flip a sixth of the
  corpus, and nothing downstream could detect it.
- `REF` is the GRCh38 reference base: over one member checked against the local refgetstore it reads
  REF on 2,814 of 2,814 rows, and the full refcheck here reports the same per chromosome. `ID` is an
  rsID on all but 93 (EA) and 114 (AA) rows; the rest get a null `rs_number`, and examples go in the
  ingestion report.
- The cis window is +/-500 kb of the protein-coding gene's TSS, and **no row sits off its protein's
  own chromosome**, so staging by chromosome needs no cross-chromosome shuffle.
- `seqid.txt` maps every SeqId to `uniprot_id`, a gene symbol, a chromosome and a TSS. It publishes
  **no ENSG**, so `gene_id` comes from `pipeline/genemap.py`, which resolves the UniProt accession
  before the symbol -- the symbol alone maps `PACAP` to ADCYAP1 when the accession says MZB1.

Two things the source does not have, and what this adapter does about them:

- **No permutation pass.** PLINK2 `--glm` publishes no `p_perm` or `p_beta`, so `permuted` carries
  neither column (CONTRACT.md: they are conventional, not required) and the experiment declares a
  rule on `p_bonferroni`, which this adapter computes as `min(1, p_lead * n_variants)` from two
  published numbers. On TOPCHeF, where both exist, that rule agrees with `p_perm < 0.05` on 93.3% of
  genes. A looser calibrated correction was rejected: the factor would have to be borrowed from
  TOPCHeF's heart tissue and European-ancestry LD, and three of these cohorts have shorter LD
  blocks, so it would inflate significance in exactly the experiments it was least valid for.
- **No fine-mapping.** `credible_sets.parquet` is written with zero rows, which the contract asks
  for explicitly rather than omitting the table.

Outputs go to `_tables/<experiment>/`: `sites.parquet`, `nominal/chr=<c>/data.parquet`,
`permuted.parquet`, `credible_sets.parquet`, `phenotypes.parquet`, `ingestion.json`. Staged
per-chromosome copies of the members live in `_source/` next to them, so iterating on the contract
tables does not re-read 3.1 GB out of the zips each time.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import time
import zipfile
from pathlib import Path
from typing import Callable

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .. import dof as doffit
from .. import genemap
from .. import qtlstore as qs
from ..common import CHROMS, Config, connect, log, write_parquet
from ..steps_refget import reference_json
from .eqtl_catalogue import allele_reads, reference_base, sequence_loader

DEFAULT_EXPERIMENT = "aric_plasma_ea"
PHENOTYPE_TYPE = "protein_somalogic"
ALLELE_ORIENTATION_SOURCE = "aric_plink2_a1_is_minor"
DOF_SAMPLE = 1_000_000
ROW_GROUP = 1 << 16
RSID = re.compile(r"^rs(\d+)$")
MEMBER = re.compile(r"(SeqId_[0-9_]+)\.PHENO1\.glm\.linear$")

# the member columns, in the order PLINK2 writes them; `#CHROM` loses its hash
COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "OBS_CT",
           "BETA", "SE", "T_STAT", "P", "ERRCODE")

SITES_SCHEMA = pa.schema([
    ("chr", pa.string()), ("pos", pa.int32()), ("ref", pa.string()), ("alt", pa.string()),
    ("af", pa.float32()), ("ma_samples", pa.int32()), ("ma_count", pa.int32()),
    ("rs_number", pa.int64()), ("rsid", pa.string()), ("match", pa.string()), ("in_cis", pa.bool_()),
])

STAGE_SCHEMA = pa.schema([
    ("phenotype_id", pa.string()), ("pos", pa.int32()), ("rs_number", pa.int64()),
    ("ref", pa.string()), ("alt", pa.string()), ("a1", pa.string()), ("a1_freq", pa.float64()),
    ("beta", pa.float64()), ("se", pa.float64()), ("pvalue", pa.float64()), ("obs_ct", pa.int32()),
])


class Experiment:
    """One `aric_pqtl:` entry of config.yaml, with its source files resolved."""

    def __init__(self, cfg: Config, name: str = DEFAULT_EXPERIMENT, chroms: list[str] | None = None):
        if name not in (cfg["aric_pqtl"] or {}):
            raise SystemExit(f"no aric_pqtl entry {name!r} in config.yaml; have {sorted(cfg['aric_pqtl'])}")
        block = cfg["aric_pqtl"][name]
        self.cfg, self.name, self.block = cfg, name, block
        root = Path(os.environ.get("QTLB_ARIC_ROOT") or cfg.raw)
        self.zip_path = root / block["zip"]
        self.seqid_path = root / block["seqid"]
        self.n_samples = int(block["n_samples"])
        self.chroms = list(chroms or CHROMS)
        self.tables = cfg.tables / name

    # sources
    def members(self) -> dict[str, str]:
        """SeqId -> member name, for every `.glm.linear` in the zip."""
        with zipfile.ZipFile(self.zip_path) as z:
            out = {}
            for n in z.namelist():
                m = MEMBER.search(n)
                if m:
                    out[m.group(1)] = n
            return out

    def seqids(self) -> dict[str, dict]:
        """SeqId -> its `seqid.txt` row: uniprot id, symbol, chromosome, TSS."""
        with open(self.seqid_path, newline="", encoding="utf-8") as fh:
            return {r["seqid_in_sample"]: {"uniprot_id": r["uniprot_id"].strip(),
                                           "symbol": r["entrezgenesymbol"].strip(),
                                           "chr": f"chr{r['chromosome_name'].strip()}",
                                           "tss": int(r["transcription_start_site"])}
                    for r in csv.DictReader(fh, delimiter="\t")}

    def significance(self) -> dict:
        return dict(self.block["significance"])

    # outputs
    def stage_dir(self) -> Path:
        return self.tables / "_source" / "nominal"

    def stage_glob(self, chrom: str | None = None) -> str:
        return str(self.stage_dir() / f"chr={chrom or '*'}" / "*.parquet")

    def stage_sql(self, chrom: str | None = None) -> str:
        return f"read_parquet('{self.stage_glob(chrom)}', hive_partitioning = true)"

    def path(self, name: str) -> Path:
        return self.tables / f"{name}.parquet"

    def nominal_dir(self) -> Path:
        return self.tables / "nominal"

    def report_path(self) -> Path:
        return self.tables / "ingestion.json"


# ---- staging -----------------------------------------------------------------------------------
def _parse_member(raw: bytes, seqid: str, want_chrom: str, cols: dict[str, list], bad_ids: list) -> int:
    """Append one member's rows to `cols`; returns how many were dropped for sitting on another
    chromosome. The header is checked rather than assumed: a silent column reorder would mis-assign
    every number, and these files carry no other integrity marker.

    `want_chrom` is the protein's chromosome from `seqid.txt`. The staged rows carry no chromosome
    of their own -- the partition is the chromosome -- so a row that disagrees cannot be stored
    truthfully and is dropped and counted instead of being filed under the wrong one.
    """
    lines = raw.split(b"\n")
    head = tuple(lines[0].decode().lstrip("#").rstrip("\r").split("\t"))
    if head != COLUMNS:
        raise ValueError(f"{seqid}: unexpected columns {head}")
    off = 0
    for line in lines[1:]:
        if not line or line == b"\r":
            continue
        f = line.decode().rstrip("\r").split("\t")
        if f"chr{f[0]}" != want_chrom:
            off += 1
            continue
        m = RSID.match(f[2])
        if not m and len(bad_ids) < 20 and f[2] not in bad_ids:
            bad_ids.append(f[2])        # examples, deduped: one id repeats across every protein testing it
        cols["phenotype_id"].append(seqid)
        cols["pos"].append(int(f[1]))
        cols["rs_number"].append(int(m.group(1)) if m else -1)
        cols["ref"].append(f[3])
        cols["alt"].append(f[4])
        cols["a1"].append(f[5])
        cols["a1_freq"].append(float(f[6]))
        cols["obs_ct"].append(int(f[8]))
        cols["beta"].append(float(f[9]))
        cols["se"].append(float(f[10]))
        cols["pvalue"].append(float(f[12]))
    return off


def extract(exp: Experiment, force: bool = False) -> dict:
    """Zip members -> `_source/nominal/chr=<c>/data.parquet`, one file per chromosome.

    Grouped by the protein's own chromosome from `seqid.txt`, which the survey showed is also every
    one of its rows' chromosome; a row that disagrees is counted and dropped rather than written to
    the wrong partition, so the partition always means what it says.
    """
    rep = {"source": str(exp.zip_path), "members": 0, "rows": 0, "rows_off_chromosome": 0,
           "non_rsid_examples": [], "chromosomes": {}, "seqids_without_member": [],
           "members_without_seqid": []}
    members, meta = exp.members(), exp.seqids()
    rep["members_without_seqid"] = sorted(set(members) - set(meta))[:20]
    rep["seqids_without_member"] = sorted(set(meta) - set(members))[:20]
    by_chrom: dict[str, list[str]] = {}
    for sid in members:
        if sid in meta:
            by_chrom.setdefault(meta[sid]["chr"], []).append(sid)
    bad_ids: list[str] = []
    with zipfile.ZipFile(exp.zip_path) as z:
        for chrom in exp.chroms:
            out = exp.stage_dir() / f"chr={chrom}" / "data.parquet"
            sids = sorted(by_chrom.get(chrom, []))
            if out.exists() and not force:
                rep["chromosomes"][chrom] = {"staged": pq.read_metadata(out).num_rows, "cached": True}
                rep["rows"] += rep["chromosomes"][chrom]["staged"]
                rep["members"] += len(sids)
                continue
            if not sids:
                rep["chromosomes"][chrom] = {"staged": 0, "members": 0}
                continue
            cols: dict[str, list] = {k: [] for k in STAGE_SCHEMA.names}
            off = 0
            for sid in sids:
                off += _parse_member(z.read(members[sid]), sid, chrom, cols, bad_ids)
                rep["members"] += 1
            t = pa.table({k: pa.array(v, type=STAGE_SCHEMA.field(k).type) for k, v in cols.items()},
                         schema=STAGE_SCHEMA)
            write_parquet(t, out, ROW_GROUP)
            rep["chromosomes"][chrom] = {"staged": t.num_rows, "members": len(sids), "rows_off_chromosome": off}
            rep["rows"] += t.num_rows
            rep["rows_off_chromosome"] += off
            log(f"aric_pqtl {exp.name}: staged {chrom} {t.num_rows:,} rows from {len(sids)} proteins")
    rep["non_rsid_examples"] = bad_ids
    return rep


# ---- sites -------------------------------------------------------------------------------------
def sites(exp: Experiment, con, sequence: Callable[[str], np.ndarray]) -> dict:
    """`sites.parquet` and the internal `_orientation/` table, with the orientation counts.

    Every row of the source is a cis result, so `in_cis` is true throughout. `af` is the **ALT**
    frequency, which is why the orientation runs here and not only over `nominal`: `A1_FREQ`
    describes A1, and A1 is the reference allele on a sixth of the rows.

    `ma_samples` and `ma_count` are -1: PLINK2 `--glm` reports `OBS_CT`, the samples tested, which
    is not the minor-allele count, and this release carries no `MA_CT` column. Writing `OBS_CT`
    there would publish a number that means something else under a contract-defined name.

    The per-site orientation is decided once, here, and written to `_orientation/` for `nominal` to
    join against -- so one site cannot be oriented one way in `sites` and another in `nominal`.
    """
    rep = {"as_is": 0, "swapped": 0, "dropped_neither_reads": 0, "sites": 0,
           "rows_a1_is_ref": 0, "rows_a1_is_alt": 0, "a1_neither_allele": 0,
           "af_disagrees_within_site": 0, "sites_without_rsid": 0, "by_chrom": {}}
    parts = []
    for chrom in exp.chroms:
        if not list(exp.stage_dir().glob(f"chr={chrom}/*.parquet")):
            continue
        # one row per site: the source repeats a site across proteins, and every protein that tests
        # it reports the same A1, A1_FREQ and rs_number. `count(DISTINCT a1_freq) > 1` checks that
        # rather than trusting it, because a disagreement would make `af` arbitrary.
        t = con.execute(f"""
            SELECT pos, ref, alt, min(a1) AS a1, min(a1_freq) AS a1_freq, max(rs_number) AS rs_number,
                   count(DISTINCT a1_freq) > 1 AS af_differs,
                   count(*) FILTER (WHERE a1 = ref) AS n_a1_ref,
                   count(*) FILTER (WHERE a1 = alt) AS n_a1_alt
            FROM {exp.stage_sql(chrom)} GROUP BY 1, 2, 3 ORDER BY 1, 2, 3""").fetch_arrow_table()
        if not t.num_rows:
            continue
        pos = t["pos"].to_numpy().astype(np.int64)
        ref = t["ref"].to_pylist()
        alt = t["alt"].to_pylist()
        a1 = np.array(t["a1"].to_pylist(), dtype=object)
        ref_a = np.array(ref, dtype=object)
        alt_a = np.array(alt, dtype=object)
        ref_ok, alt_ok = allele_reads(sequence(chrom), pos, ref, alt)
        rbase = reference_base(ref_a, alt_a, ref_ok, alt_ok)
        # the effect allele is A1; the other allele is whichever of REF/ALT it is not. A1 that is
        # neither cannot be oriented at all, and `orient_to_ref` drops it as `neither`.
        other = np.where(a1 == alt_a, ref_a, np.where(a1 == ref_a, alt_a, None)).astype(object)
        rep["a1_neither_allele"] += int(((a1 != alt_a) & (a1 != ref_a)).sum())
        r = qs.orient_to_ref(rbase, a1, other, np.zeros(len(pos)), t["a1_freq"].to_numpy())
        keep = r["keep"]
        # `orient_to_ref`'s swap, recomputed per row because `nominal` needs it per row and the
        # helper returns only totals. Note what it does **not** mean here: ARIC already publishes
        # REF and ALT correctly, so a swap leaves the allele pair alone and only negates the effect.
        # Comparing the oriented ref against the source's would therefore never fire.
        swap = (other != rbase) & (a1 == rbase)
        for k in ("as_is", "swapped"):
            rep[k] += r["counts"][k]
        rep["dropped_neither_reads"] += r["counts"]["dropped"]
        rep["af_disagrees_within_site"] += int(np.count_nonzero(t["af_differs"].to_numpy(zero_copy_only=False)))
        rep["rows_a1_is_ref"] += int(t["n_a1_ref"].to_numpy().sum())
        rep["rows_a1_is_alt"] += int(t["n_a1_alt"].to_numpy().sum())
        rs = t["rs_number"].to_numpy()[keep]
        rep["sites_without_rsid"] += int((rs < 0).sum())
        chrom_sites = pa.table({
            "chr": pa.array([chrom] * int(keep.sum()), pa.string()),
            "pos": pa.array(pos[keep], pa.int32()),
            "ref": pa.array(r["ref"], pa.string()),
            "alt": pa.array(r["alt"], pa.string()),
            "af": pa.array(r["af"], pa.float32()),
            "ma_samples": pa.array(np.full(int(keep.sum()), -1), pa.int32()),
            "ma_count": pa.array(np.full(int(keep.sum()), -1), pa.int32()),
            "rs_number": pa.array(rs, pa.int64()),
            "rsid": pa.array([None if n < 0 else f"rs{n}" for n in rs], pa.string()),
            "match": pa.array([None if n < 0 else "exact" for n in rs], pa.string()),
            "in_cis": pa.array(np.ones(int(keep.sum()), dtype=bool), pa.bool_()),
        })
        parts.append(chrom_sites)
        # what `nominal` joins to: the source's (pos, a1) -> the oriented site and whether beta flips
        orient = pa.table({
            "pos": pa.array(pos[keep], pa.int32()),
            "src_ref": pa.array(ref_a[keep].tolist(), pa.string()),
            "src_alt": pa.array(alt_a[keep].tolist(), pa.string()),
            "ref": pa.array(r["ref"], pa.string()),
            "alt": pa.array(r["alt"], pa.string()),
            "flip": pa.array(swap[keep].tolist(), pa.bool_()),
        })
        write_parquet(orient, exp.tables / "_orientation" / f"chr={chrom}" / "data.parquet", ROW_GROUP)
        rep["by_chrom"][chrom] = {"sites": chrom_sites.num_rows, "dropped": int((~keep).sum())}
        log(f"aric_pqtl {exp.name}: {chrom} {chrom_sites.num_rows:,} sites, {int((~keep).sum())} dropped")
    table = pa.concat_tables(parts) if parts else pa.table({f.name: pa.array([], f.type) for f in SITES_SCHEMA})
    write_parquet(table, exp.path("sites"), ROW_GROUP, stats_columns=["chr", "pos"])
    rep["sites"] = table.num_rows
    return rep


# ---- phenotypes, nominal, permuted, credible sets -----------------------------------------------
def phenotypes(exp: Experiment, con, resolver: genemap.Resolver) -> dict:
    """`phenotypes.parquet`: one row per SOMAmer.

    `phenotype_object_id` equals `phenotype_id`: PLINK2 tested each protein on its own, so a group is
    one phenotype, as it is for TOPCHeF. `gene_id` comes from `genemap`, accession first.

    `extra` carries what is true of the assay rather than of the gene: the SomaLogic SeqId's UniProt
    accession, the symbol the source published, the TSS its cis window was drawn around, and how the
    gene was resolved. The TSS is **not** annotation -- it is the number this study centred its
    window on, which may differ from the annotation's TSS for the same gene, and a reader comparing
    window extents needs the study's own.
    """
    rep = {"phenotypes": 0, "with_nominal": 0, "without_gene": 0, "genemap": {}}
    meta = exp.seqids()
    staged = con.execute(f"""SELECT phenotype_id, count(*) AS n FROM {exp.stage_sql()} GROUP BY 1""").fetchall() \
        if list(exp.stage_dir().glob("chr=*/*.parquet")) else []
    n_rows = dict(staged)
    gm = resolver.report([(m["symbol"], m["uniprot_id"]) for m in meta.values()])
    rep["genemap"] = {k: v for k, v in gm.items() if k != "map"}
    gene_of = gm["map"]
    ids = sorted(meta)
    rows = []
    for sid in ids:
        m = meta[sid]
        gid = gene_of.get((m["symbol"], m["uniprot_id"]))
        _, how = resolver.resolve(m["symbol"], m["uniprot_id"])
        rows.append((PHENOTYPE_TYPE, sid, sid, gid, sid in n_rows,
                     json.dumps({"uniprot_id": m["uniprot_id"] or None, "source_symbol": m["symbol"] or None,
                                 "source_tss": m["tss"], "source_chr": m["chr"], "gene_id_from": how},
                                separators=(",", ":"), sort_keys=True)))
        rep["without_gene"] += gid is None
    t = pa.table({
        "phenotype_type": pa.array([r[0] for r in rows], pa.string()),
        "phenotype_id": pa.array([r[1] for r in rows], pa.string()),
        "phenotype_object_id": pa.array([r[2] for r in rows], pa.string()),
        "gene_id": pa.array([r[3] for r in rows], pa.string()),
        "has_nominal": pa.array([r[4] for r in rows], pa.bool_()),
        "extra": pa.array([r[5] for r in rows], pa.string()),
    })
    write_parquet(t, exp.path("phenotypes"), ROW_GROUP)
    rep["phenotypes"] = t.num_rows
    rep["with_nominal"] = int(sum(1 for r in rows if r[4]))
    return rep


def _oriented(exp: Experiment, chrom: str) -> str:
    return f"read_parquet('{exp.tables / '_orientation' / f'chr={chrom}' / 'data.parquet'}')"


def nominal(exp: Experiment, con) -> dict:
    """`nominal/chr=<c>/data.parquet`: one row per (protein, variant), `beta` ALT-relative.

    The sign comes from the per-site `flip` that `sites` decided, joined on the source's own
    `(pos, ref, alt)`. Deciding it again here from the row's A1 would be the same arithmetic twice,
    and two copies of a sign rule is how one of them ends up wrong.

    A row whose site `sites` dropped is dropped here too -- the contract fails the build on a
    nominal row with no site, so the two must agree by construction rather than by luck.
    """
    rep = {"rows": 0, "rows_dropped_site": 0, "rows_beta_negated": 0, "by_chrom": {}}
    out_dir = exp.nominal_dir()
    for chrom in exp.chroms:
        if not list(exp.stage_dir().glob(f"chr={chrom}/*.parquet")):
            continue
        t = con.execute(f"""
            SELECT n.phenotype_id, p.gene_id, n.pos, o.ref, o.alt,
                   CASE WHEN o.flip THEN -n.beta ELSE n.beta END AS beta, n.se, n.pvalue,
                   o.flip
            FROM {exp.stage_sql(chrom)} n
            JOIN {_oriented(exp, chrom)} o ON o.pos = n.pos AND o.src_ref = n.ref AND o.src_alt = n.alt
            LEFT JOIN read_parquet('{exp.path("phenotypes")}') p
                   ON p.phenotype_id = n.phenotype_id
            ORDER BY n.pos, o.ref, o.alt, n.phenotype_id""").fetch_arrow_table()
        staged = con.execute(f"SELECT count(*) FROM {exp.stage_sql(chrom)}").fetchone()[0]
        rep["rows_dropped_site"] += staged - t.num_rows
        rep["rows_beta_negated"] += int(sum(t["flip"].to_pylist()))
        out = pa.table({
            "phenotype_type": pa.array([PHENOTYPE_TYPE] * t.num_rows, pa.string()),
            "phenotype_id": t["phenotype_id"].cast(pa.string()),
            "gene_id": t["gene_id"].cast(pa.string()),
            "chr": pa.array([chrom] * t.num_rows, pa.string()),
            "pos": t["pos"].cast(pa.int32()),
            "ref": t["ref"].cast(pa.string()),
            "alt": t["alt"].cast(pa.string()),
            "beta": t["beta"].cast(pa.float32()),
            "se": t["se"].cast(pa.float32()),
            "pvalue": t["pvalue"].cast(pa.float64()),
        })
        write_parquet(out, out_dir / f"chr={chrom}" / "data.parquet", ROW_GROUP, stats_columns=["pos"])
        rep["by_chrom"][chrom] = out.num_rows
        rep["rows"] += out.num_rows
        log(f"aric_pqtl {exp.name}: nominal {chrom} {out.num_rows:,} rows")
    return rep


def permuted(exp: Experiment, con) -> dict:
    """`permuted.parquet`: one row per protein, with its lead variant and a Bonferroni p.

    There is no permutation pass to copy, so `p_perm` and `p_beta` are absent (CONTRACT.md: they are
    conventional, not required) and `p_bonferroni = min(1, p_lead * n_variants)` stands in, which the
    experiment's significance rule names. `n_variants` is the recount from this experiment's own
    nominal rows, which is both what the contract prefers and what the Bonferroni needs.

    The lead is the smallest `pvalue`, tie-broken by position then by the allele pair, so a protein
    with two equally small p-values gets the same lead on every rebuild.
    """
    rep = {"groups": 0, "lead_ties_broken": 0, "p_bonferroni_at_one": 0}
    t = con.execute(f"""
        WITH n AS (SELECT * FROM read_parquet('{exp.nominal_dir()}/chr=*/data.parquet', hive_partitioning = true)),
        counted AS (SELECT phenotype_id, count(*) AS n_variants, min(pvalue) AS p_lead FROM n GROUP BY 1),
        ranked AS (
            SELECT n.phenotype_id, n.gene_id, n.chr, n.pos, n.ref, n.alt, n.pvalue,
                   row_number() OVER (PARTITION BY n.phenotype_id ORDER BY n.pvalue, n.pos, n.ref, n.alt) AS rk,
                   count(*) OVER (PARTITION BY n.phenotype_id, n.pvalue) AS at_this_p
            FROM n JOIN counted c USING (phenotype_id) WHERE n.pvalue = c.p_lead)
        SELECT r.phenotype_id, r.gene_id, c.n_variants, c.p_lead,
               least(1.0, c.p_lead * c.n_variants) AS p_bonferroni,
               r.chr AS lead_chr, r.pos AS lead_pos, r.ref AS lead_ref, r.alt AS lead_alt,
               r.at_this_p > 1 AS tied
        FROM ranked r JOIN counted c USING (phenotype_id) WHERE r.rk = 1
        ORDER BY r.phenotype_id""").fetch_arrow_table()
    rep["lead_ties_broken"] = int(sum(t["tied"].to_pylist()))
    rep["p_bonferroni_at_one"] = int(sum(1 for v in t["p_bonferroni"].to_pylist() if v is not None and v >= 1.0))
    out = pa.table({
        "phenotype_type": pa.array([PHENOTYPE_TYPE] * t.num_rows, pa.string()),
        "phenotype_object_id": t["phenotype_id"].cast(pa.string()),
        "phenotype_id": t["phenotype_id"].cast(pa.string()),
        "gene_id": t["gene_id"].cast(pa.string()),
        "n_variants": t["n_variants"].cast(pa.int32()),
        "p_bonferroni": t["p_bonferroni"].cast(pa.float64()),
        "p_nominal_lead": t["p_lead"].cast(pa.float64()),
        "lead_chr": t["lead_chr"].cast(pa.string()),
        "lead_pos": t["lead_pos"].cast(pa.int32()),
        "lead_ref": t["lead_ref"].cast(pa.string()),
        "lead_alt": t["lead_alt"].cast(pa.string()),
    })
    write_parquet(out, exp.path("permuted"), ROW_GROUP)
    rep["groups"] = out.num_rows
    return rep


def credible_sets(exp: Experiment) -> dict:
    """`credible_sets.parquet` with zero rows: this release publishes no fine-mapping. The contract
    asks for the empty table rather than no table, so a reader can tell "none" from "not built"."""
    t = pa.table({f.name: pa.array([], f.type) for f in pa.schema([
        ("phenotype_type", pa.string()), ("phenotype_object_id", pa.string()), ("phenotype_id", pa.string()),
        ("cs_id", pa.int16()), ("chr", pa.string()), ("pos", pa.int32()), ("ref", pa.string()),
        ("alt", pa.string()), ("pip", pa.float32()), ("z", pa.float32()), ("cs_size", pa.int32()),
        ("cs_min_r2", pa.float32())])})
    write_parquet(t, exp.path("credible_sets"), ROW_GROUP)
    return {"rows": 0, "reason": "the release publishes no fine-mapping"}


def fit_dof(exp: Experiment, con) -> dict:
    """Degrees of freedom from the contract nominal rows (`pipeline/dof.py`).

    `OBS_CT` is constant across the whole cohort, so `n_samples` is the published N rather than an
    estimate, and the fit's `implied_covariates` is then a real check on it: a wildly wrong value
    there means the N or the model is not what the file says.
    """
    t = con.execute(f"""SELECT beta, se, pvalue FROM
        read_parquet('{exp.nominal_dir()}/chr=*/data.parquet', hive_partitioning = true)
        USING SAMPLE reservoir({DOF_SAMPLE} ROWS) REPEATABLE (0)""").fetch_arrow_table()
    if not t.num_rows:
        return {PHENOTYPE_TYPE: {"dof": None, "usable": False, "reason": "no nominal rows"}}
    r = doffit.fit(t["beta"].to_numpy(), t["se"].to_numpy(), t["pvalue"].to_numpy(), n_samples=exp.n_samples)
    out = {k: r.get(k) for k in ("dof", "usable", "identified", "residual_log10p", "margin",
                                 "rows", "rows_usable", "n_samples", "implied_covariates", "reason")}
    out["dof_for_manifest"] = doffit.dof_for_manifest(r)
    return {PHENOTYPE_TYPE: out}


# ---- the whole adapter -------------------------------------------------------------------------
def run(exp: Experiment, sequence: Callable[[str], np.ndarray] | None = None,
        resolver: genemap.Resolver | None = None, force_extract: bool = False) -> dict:
    t0 = time.time()
    sequence = sequence or sequence_loader(exp.cfg)
    resolver = resolver or genemap.load(exp.cfg)
    con = connect(exp.cfg, memory_limit=exp.cfg["duckdb_memory_limit"], threads=exp.cfg["duckdb_threads"],
                  temp_dir=exp.cfg.tmp / f"aric-{exp.name}")
    report: dict = {"experiment_id": exp.name, "allele_orientation_source": ALLELE_ORIENTATION_SOURCE,
                    "chromosomes": exp.chroms, "phenotype_types": [PHENOTYPE_TYPE]}
    report["staging"] = extract(exp, force=force_extract)
    report["orientation"] = sites(exp, con, sequence)
    report["phenotypes"] = phenotypes(exp, con, resolver)
    report["nominal"] = nominal(exp, con)
    report["permuted"] = permuted(exp, con)
    report["credible_sets"] = credible_sets(exp)
    report["dof"] = fit_dof(exp, con)
    report["significance"] = exp.significance()
    report["source"] = {
        "zip": str(exp.zip_path), "seqid": str(exp.seqid_path), "member_prefix": exp.block.get("member_prefix"),
        "n_samples": exp.n_samples, "cis_window_bp": exp.block.get("cis_window_bp"),
        "gene_annotation": exp.block.get("gene_annotation"),
        # provenance the store can carry: without it a published derivative of an S3 bucket with no
        # DOI says only that some bytes at some path produced it (see the plan's blast radius)
        "files": _provenance(exp),
    }
    report["rows"] = {k: pq.read_metadata(exp.path(k)).num_rows
                      for k in ("sites", "phenotypes", "permuted", "credible_sets")}
    report["rows"]["nominal"] = report["nominal"]["rows"]
    exp.report_path().write_text(json.dumps(report, indent=2, default=str) + "\n")
    log(f"aric_pqtl {exp.name}: {report['rows']} in {(time.time() - t0) / 60:.1f} min -> {exp.report_path()}")
    return report


def _provenance(exp: Experiment) -> list[dict]:
    """The `sources.yaml` entries for this experiment's inputs, matched the way
    `annotation.gtf_source` does, so the experiment pointer can name a URL and a checksum."""
    import yaml
    p = exp.cfg.raw / "sources.yaml"
    if not p.exists():
        return []
    doc = yaml.safe_load(p.read_text())
    want = {exp.zip_path.name, exp.seqid_path.name}
    out = []
    for src in doc.get("sources", []):
        if src.get("dir") != exp.zip_path.parent.name:
            continue
        for f in src.get("files", []):
            url = f.get("url") or (src.get("url_base") or "").format(file=f.get("file", ""))
            name = Path(f.get("file") or url).name
            if name in want:
                out.append({"file": name, "url": url, "md5": f.get("md5"), "size": f.get("size"),
                            "source": src.get("name"), "version": src.get("version")})
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--experiment", default=DEFAULT_EXPERIMENT)
    ap.add_argument("--force-extract", action="store_true", help="restage the zip members even when cached")
    args = ap.parse_args(argv)
    cfg = Config()
    chroms = None
    env = os.environ.get("QTLB_CHROMS")
    if env:
        chroms = [c.strip() for c in env.split(",") if c.strip()]
    exp = Experiment(cfg, args.experiment, chroms)
    if not reference_json(cfg).exists():
        log(f"note: no {reference_json(cfg)}; sequence digests come straight from the configured collection")
    run(exp, force_extract=args.force_extract)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
