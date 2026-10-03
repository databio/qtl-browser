"""Tests for the eQTL Catalogue adapter (pipeline/adapters/eqtl_catalogue.py).

    uv run python -m pipeline.adapters.test_eqtl_catalogue

Plain asserts, synthetic data only: a ten-base-repeat chr21, a handful of TSVs laid out like the
Catalogue's FTP tree, and a fake sequence loader in place of the refgetstore. The rules pinned here
are the ones a real run cannot show: a site whose `ref` the genome does not read swaps or drops, an
indel whose two alleles both prefix-match keeps the source's REF, every intron of a tested cluster
becomes a phenotype, and `permuted` is one row per group.
"""
from __future__ import annotations

import gzip
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from scipy.special import stdtr

from ..common import Config
from . import eqtl_catalogue as ec

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


# 1-based: 1 A, 2 C, 3 G, 4 T, 5 T, 6 G, 7 C, 8 A, 9 A, 10 C, then repeat
SEQ = np.frombuffer(b"ACGTTGCAAC" * 10, dtype=np.uint8)

NOM_HEADER = list(ec.NOMINAL_COLUMNS)
PERM_HEADER = list(ec.PERMUTED_COLUMNS)
CS_HEADER = list(ec.CS_COLUMNS)
META_HEADER = ["phenotype_id", "quant_id", "group_id", "gene_id", "chromosome", "gene_start", "gene_end",
               "strand", "gene_name", "gene_type", "gene_version", "phenotype_pos", "gene_gc_content",
               "intron_start", "intron_end", "group_start", "group_end", "gene_count", "gene_pos"]

# (variant, ac, an, ma_samples)
V1 = ("chr21_1_A_G", 100, 764, 90)      # ref reads
V2 = ("chr21_2_T_C", 600, 764, 150)     # the genome reads C = alt: swaps
V3 = ("chr21_3_A_T", 50, 764, 45)       # the genome reads G: neither, dropped
V4 = ("chr21_4_T_TT", 30, 764, 28)      # both prefix-match (genome TT); source REF stands
V5 = "chr21_6_G_A"                      # only a leafcutter lead names it


def _p(beta, se, dof=367):
    return float(2 * stdtr(dof, -abs(beta / se)))


def _nom(trait, obj, gene, v, beta, se):
    variant, ac, an, ma = v
    c, pos, ref, alt = variant.split("_")
    return [trait, c[3:], pos, ref, alt, variant, ma, min(ac, an - ac) / an, f"{_p(beta, se):.6g}", beta, se,
            "SNP", ac, an, "NA", obj, gene, 1.0, "rs" + pos]


def _write(path: Path, header: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("NA" if x is None else str(x) for x in r) + "\n")


def build(tmp: Path) -> tuple[ec.Experiment, dict]:
    root = tmp / "r6"
    ss, su = root / "sumstats" / "S1", root / "susie" / "S1"
    ge_rows = [_nom("ENSG1", "ENSG1", "ENSG1", V1, 0.8, 0.1), _nom("ENSG1", "ENSG1", "ENSG1", V2, 0.5, 0.1),
               _nom("ENSG1", "ENSG1", "ENSG1", V3, 0.2, 0.1), _nom("ENSG1", "ENSG1", "ENSG1", V4, -0.3, 0.1),
               _nom("ENSG2", "ENSG2", "ENSG2", V1, 0.05, 0.1)]
    ge_rows.append(ge_rows[0][:-1] + ["rs999"])    # the source repeats a row per dbSNP id
    lc_rows = [_nom("21:10:20:clu_1_+", "clu_1_+", "ENSG1", V1, 0.9, 0.1),
               _nom("21:10:20:clu_1_+", "clu_1_+", "ENSG1", V4, 0.4, 0.1)]
    _write(ss / "D1" / "D1.all.tsv.gz", NOM_HEADER, ge_rows)
    _write(ss / "D2" / "D2.cc.tsv.gz", NOM_HEADER, lc_rows)
    _write(ss / "D1" / "D1.permuted.tsv.gz", PERM_HEADER, [
        ["ENSG1", "ENSG1", 1, 10, "chr21_1_A_G", "21", 1, 1e-10, 0.8, 0.001, 1e-8],
        ["ENSG2", "ENSG2", 1, 9, "chr21_1_A_G", "21", 1, 0.6, 0.05, 0.9, 0.9]])
    _write(ss / "D2" / "D2.permuted.tsv.gz", PERM_HEADER, [
        ["clu_1_+", "21:10:20:clu_1_+", 2, 50, "chr21_1_A_G", "21", 1, 1e-12, 0.9, 0.001, 1e-9],
        ["clu_2_-", "21:40:45:clu_2_-", 3, 60, V5, "21", 6, 1e-3, 0.2, 0.4, 0.4]])
    cs = lambda trait, gene, k, v, size, z: [trait, gene, f"{trait}_L{k}", v, "rs1", size, 0.5, 1e-8, 0.3, 0.1, z, 0.9, "chr21:-10-100"]
    _write(su / "D1" / "D1.credible_sets.tsv.gz", CS_HEADER, [
        cs("ENSG1", "ENSG1", 1, V1[0], 2, 8.0), cs("ENSG1", "ENSG1", 1, V2[0], 2, 5.0), cs("ENSG1", "ENSG1", 2, V4[0], 1, -3.0),
        cs("ENSG1", "ENSG1", 1, V2[0], 2, 5.0)[:4] + ["rs999"] + cs("ENSG1", "ENSG1", 1, V2[0], 2, 5.0)[5:]])  # repeat per dbSNP id
    _write(su / "D2" / "D2.credible_sets.tsv.gz", CS_HEADER, [cs("21:10:20:clu_1_+", "ENSG1", 1, V1[0], 1, 9.0)])
    meta = tmp / "meta.tsv.gz"
    mrow = lambda pid, clu, gene, s, e: [pid, clu, clu, gene, "21", 1, 100, 1, "G", "protein_coding", 1, 5, 50, s, e, 1, 100, 1, 1]
    _write(meta, META_HEADER, [
        mrow("21:10:20:clu_1_+", "clu_1_+", "ENSG1", 10, 20), mrow("21:10:30:clu_1_+", "clu_1_+", "ENSG1", 10, 30),
        mrow("21:40:50:clu_2_-", "clu_2_-", "ENSG3", 40, 50), mrow("21:45:50:clu_2_-", "clu_2_-", "ENSG3", 45, 50),
        mrow("21:40:45:clu_2_-", "clu_2_-", "ENSG3,ENSG4", 40, 45),
        mrow("21:70:80:clu_3_+", "clu_3_+", "ENSG5", 70, 80), mrow("21:70:90:clu_3_+", "clu_3_+", "ENSG5", 70, 90)])
    cfg = Config()
    cfg.cfg = dict(cfg.cfg)
    cfg.cfg["eqtl_catalogue"] = {"t": {
        "root": str(root), "study": "S1", "datasets": {"ge": "D1", "leafcutter": "D2"},
        "nominal": {"ge": "all", "leafcutter": "cc"}, "leafcutter_metadata": str(meta), "n_samples": 382,
        "significance": {"column": "p_perm", "op": "<", "threshold": 0.05}}}
    cfg.cfg["duckdb_memory_limit"], cfg.cfg["duckdb_threads"] = "1GB", 1
    cfg.derived = tmp / "derived"
    cfg.tables, cfg.tmp = cfg.derived / "_tables", cfg.derived / "_tmp"
    exp = ec.Experiment(cfg, "t", chroms=["chr21"])
    report = ec.run(exp, sequence=lambda chrom: SEQ)
    return exp, report


_BUILT: dict = {}


def built():
    if not _BUILT:
        tmp = Path(tempfile.mkdtemp(prefix="eqtlcat-test-"))
        _BUILT["exp"], _BUILT["report"] = build(tmp)
    return _BUILT["exp"], _BUILT["report"]


def _rows(exp, table):
    return pq.read_table(exp.path(table)).to_pylist()


@case
def allele_reads_keeps_source_ref_at_ambiguous_indel():
    ref_ok, alt_ok = ec.allele_reads(SEQ, np.array([1, 2, 3, 4]), ["A", "T", "A", "T"], ["G", "C", "T", "TT"])
    assert list(ref_ok) == [True, False, False, True]
    assert list(alt_ok) == [False, True, False, True]
    got = ec.reference_base(np.array(["A", "T", "A", "T"], dtype=object), np.array(["G", "C", "T", "TT"], dtype=object),
                            ref_ok, alt_ok)
    assert list(got) == ["A", "C", None, "T"], list(got)


@case
def sites_swap_drop_and_af():
    exp, rep = built()
    s = {(r["pos"], r["ref"], r["alt"]): r for r in _rows(exp, "sites")}
    assert set(s) == {(1, "A", "G"), (2, "C", "T"), (4, "T", "TT"), (6, "G", "A")}, set(s)
    assert abs(s[(1, "A", "G")]["af"] - 100 / 764) < 1e-6
    assert abs(s[(2, "C", "T")]["af"] - (1 - 600 / 764)) < 1e-6          # mirrored with the swap
    assert s[(2, "C", "T")]["ma_count"] == 164 and s[(2, "C", "T")]["ma_samples"] == 150
    assert np.isnan(s[(6, "G", "A")]["af"]) and s[(6, "G", "A")]["ma_samples"] == -1 and s[(6, "G", "A")]["ma_count"] == -1
    assert s[(1, "A", "G")]["rsid"] == "rs1" and s[(1, "A", "G")]["rs_number"] == 1
    assert all(r["in_cis"] for r in s.values())
    o = rep["orientation"]
    assert (o["as_is"], o["swapped"], o["dropped_neither_reads"]) == (3, 1, 1), o
    assert o["by_origin"]["lead_only"] == 1


@case
def nominal_negates_swapped_beta_and_drops_orphans():
    exp, rep = built()
    rows = pq.read_table(exp.nominal_path("chr21")).to_pylist()
    assert len(rows) == 6 and rep["nominal"]["rows_dropped_site"] == 1, (len(rows), rep["nominal"])
    by = {(r["phenotype_id"], r["pos"]): r for r in rows}
    assert by[("ENSG1", 2)]["beta"] == -0.5 and by[("ENSG1", 2)]["ref"] == "C"
    assert by[("ENSG1", 1)]["beta"] == np.float32(0.8)
    assert by[("21:10:20:clu_1_+", 4)]["gene_id"] == "ENSG1"
    assert ("ENSG1", 3) not in by
    assert rep["nominal"]["duplicate_rows_rsid_only"] == 1


@case
def phenotypes_enumerate_every_intron_of_tested_clusters():
    exp, rep = built()
    p = {(r["phenotype_type"], r["phenotype_id"]): r for r in _rows(exp, "phenotypes")}
    assert len(p) == 7, sorted(p)
    assert ("leafcutter", "21:70:80:clu_3_+") not in p                 # untested cluster
    assert rep["phenotypes"]["leafcutter"]["metadata"]["untested_introns_in_metadata"] == 2
    lc = p[("leafcutter", "21:45:50:clu_2_-")]
    assert lc["phenotype_object_id"] == "clu_2_-" and lc["gene_id"] == "ENSG3" and not lc["has_nominal"]
    assert json.loads(lc["extra"]) == {"intron_start": 45, "intron_end": 50, "cluster_id": "clu_2_-", "strand": "-"}
    multi = json.loads(p[("leafcutter", "21:40:45:clu_2_-")]["extra"])
    assert multi["gene_ids"] == ["ENSG3", "ENSG4"] and p[("leafcutter", "21:40:45:clu_2_-")]["gene_id"] == "ENSG3"
    assert p[("ge", "ENSG2")]["has_nominal"] and p[("leafcutter", "21:10:20:clu_1_+")]["has_nominal"]
    assert p[("ge", "ENSG1")]["phenotype_object_id"] == "ENSG1" and p[("ge", "ENSG1")]["extra"] == "{}"
    assert rep["phenotypes"]["leafcutter"]["has_nominal"] == 1 and rep["phenotypes"]["leafcutter"]["phenotypes"] == 5


@case
def permuted_one_row_per_group_with_recount_for_ge_only():
    exp, rep = built()
    p = {(r["phenotype_type"], r["phenotype_object_id"]): r for r in _rows(exp, "permuted")}
    assert len(p) == 4
    assert p[("ge", "ENSG1")]["n_variants"] == 3                        # recounted: V3 was dropped
    assert p[("ge", "ENSG2")]["n_variants"] == 1
    assert p[("leafcutter", "clu_1_+")]["n_variants"] == 50               # copied
    lead = p[("leafcutter", "clu_2_-")]
    assert lead["phenotype_id"] == "21:40:45:clu_2_-" and lead["gene_id"] == "ENSG3"
    assert (lead["lead_chr"], lead["lead_pos"], lead["lead_ref"], lead["lead_alt"]) == ("chr21", 6, "G", "A")
    assert rep["permuted"]["ge"]["n_variants"] == "recounted_from_nominal"
    assert rep["permuted"]["leafcutter"]["n_variants"] == "copied_from_source"


@case
def credible_sets_parse_cs_id_and_flip_z():
    exp, rep = built()
    t = pq.read_table(exp.path("credible_sets"))
    assert str(t.schema.field("cs_id").type) == "int16"
    rows = t.to_pylist()
    assert sorted((r["phenotype_id"], r["cs_id"], r["pos"]) for r in rows) == [
        ("21:10:20:clu_1_+", 1, 1), ("ENSG1", 1, 1), ("ENSG1", 1, 2), ("ENSG1", 2, 4)]
    swapped = [r for r in rows if r["pos"] == 2][0]
    assert swapped["z"] == -5.0 and swapped["ref"] == "C"
    assert [r for r in rows if r["phenotype_type"] == "leafcutter"][0]["phenotype_object_id"] == "clu_1_+"
    assert rep["credible_sets"]["ge"]["duplicate_rows_rsid_only"] == 1
    assert rep["credible_sets"]["ge"]["cs_size_differs_from_rows"] == 0


@case
def ingestion_report_is_complete():
    exp, rep = built()
    on_disk = json.loads(exp.report_path().read_text())
    assert on_disk["allele_orientation_source"] == "eqtl_catalogue_ref_alt"
    assert set(on_disk["dof"]) == {"ge", "leafcutter"}
    assert on_disk["rows"] == {"sites": 4, "phenotypes": 7, "permuted": 4, "credible_sets": 4, "nominal": 6}
    assert on_disk["source"]["staged"]["ge"]["rows_variant_disagrees_with_columns"] == 0


def main() -> int:
    bad = 0
    for fn in CASES:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001 - report every case
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(CASES) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
