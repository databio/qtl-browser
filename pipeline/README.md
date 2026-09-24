# pipeline

Builds the qtlb v1 store (`SPEC.md` at the repo root) that the browser reads. Each study goes
through an adapter into contract tables (`CONTRACT.md`), then into a refget-anchored store, then
through checks and a benchmark. Every step runs on Rivanna through Slurm from the synced checkout
(`~/scratch/qtl-browser-live`).

```bash
# v1: adapters, store, checks, benchmark (Slurm; details below)
QTLB_DERIVED=<tree> sbatch adapter.sbatch                          # TOPCHeF adapter + gate + contract check
QTLB_STORE=<store> EXPERIMENTS="topchef:<tree>/_tables/topchef" sbatch store.sbatch   # store build + validate + verify_v0
QTLB_STORE=<store> sbatch bench_store.sbatch                       # v0 vs v1 sizes and read speed
uv run python -m pipeline.qtlstore validate | remove-experiment ID | gc [--dry-run] | crosscat A B --store <store>

# the shared input tables the TOPCHeF adapter and its gate read (build.sbatch)
uv run python -m pipeline steps                                    # list steps in order
uv run python -m pipeline build                                    # run what is not done yet
uv run python -m pipeline build --step nominal --force

# tests (no data needed)
uv run python -m pipeline.test_qtlstore      # also test_catalog, test_annotation, test_results, test_dof,
                                             # test_gtf, test_refcheck, adapters.test_topchef, adapters.test_eqtl_catalogue
uv run python -m pipeline.test_packfmt --synthetic && uv run python -m pipeline.test_packtool   # the v0 reader

uv run python -m pipeline.packtool header | blocks | block | variants | frames | trans | gwas | gwas-index | check ...   # inspect frozen v0 packs
uv run python -m pipeline.figures manhattan | density 5 10 | gwas 5 | themes                                          # quick-look PNGs from the v0 tables
```

## The v1 build: adapters, store, checks, benchmark

The v1 path turns each study into contract tables (`CONTRACT.md`), then builds a refget-anchored
store (`qtlstore.py`) holding one variant catalog and one results set per experiment (cis results,
search index, paged hits, trans results, and a GWAS when the tables have one), sharing annotations.
The v1 modules (`qtlstore`, `catalog`, `annotation`, `results`, `gwas`) use the codec
`packfmt_v1.py` and never `packfmt_v0.py`; `verify_v0.py` and `bench_store.py` compare the two formats
and are the only modules that import both. Each adapter writes into its own `QTLB_DERIVED` tree; nothing here
writes into the frozen v0 tree (`adapter.sbatch` refuses it). Iterate with
`QTLB_CHROMS=chr21,chr22` first.

| Step | Command | Module |
|---|---|---|
| Adapter + TOPCHeF gate + contract check | `QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/topchef.py`, `adapters/verify_topchef.py`, `adapters/contract_check.py` |
| eQTL Catalogue adapter + contract check | `ADAPTER=eqtl_catalogue QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/eqtl_catalogue.py` |
| One step alone | `STEPS=gate` or `STEPS=contract` (or `"adapter contract"`) with the same variables | |
| One TOPCHeF table | `STEPS=adapter QTLB_DERIVED=<tree> sbatch adapter.sbatch --step phenotypes --step trans` (trans needs the trans-only phenotypes, so both) | `adapters/topchef.py` |
| DCM GWAS into the TOPCHeF tables | `ADAPTER=dcm_gwas STEPS="adapter contract" QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/dcm_gwas.py` |
| Only the GWAS bin table (`gwas_bins.parquet`, v0's `gwas_dcm_bins.json` bins; seconds) | `ADAPTER=dcm_gwas STEPS="adapter contract" QTLB_DERIVED=<tree> sbatch adapter.sbatch --bins-only` | `adapters/dcm_gwas.py` |
| Store build | `QTLB_CHROMS=all QTLB_STORE=<store> EXPERIMENTS="topchef:<tree>/_tables/topchef:gencode_v34 gtex_v8_heart_lv:<tree>/_tables/gtex_v8_heart_lv" CROSSCAT="topchef_grch38 gtex_v8_heart_lv_grch38" sbatch --time=6:00:00 --mem=64G store.sbatch` | `annotation.py`, `catalog.py`, `results.py`, `gwas.py`, `verify_v0.py` |
| Store maintenance | `uv run python -m pipeline.qtlstore validate \| remove-experiment ID \| gc [--dry-run] \| crosscat A B --store <store>` | `qtlstore.py` |
| Per-chromosome objects for a store built before them | `uv run python -m pipeline.annotation add-split --store <store> --id <annotation>`, then `uv run python -m pipeline.results add-split --store <store> --id <experiment>` (rewrites the pointer; no other object changes) | `annotation.py`, `results.py` |
| Benchmark, v0 vs v1 | `QTLB_STORE=<store> sbatch bench_store.sbatch`; `EXPERIMENT=gtex_v8_heart_lv ... sbatch bench_store.sbatch --no-reads` for a study with no v0 twin | `bench_store.py` |

Outside Slurm the same modules run directly, e.g.
`uv run python -m pipeline.adapters.contract_check $QTLB_DERIVED/_tables/topchef` (prints one
PASS/FAIL line per `CONTRACT.md` rule, exits 1 on any failure).

- **Gate** (TOPCHeF only): every contract number against the frozen v0 build, bit for bit, slope sign
  included. `GENES` sets the per-chromosome pack sample (default 20).
- **Store build** (`store.sbatch`): per experiment an annotation (named in `EXPERIMENTS`, else the
  adapter's `ingestion.json` `source.gene_annotation`, else `gencode_v34`; GTFs from
  `ANNOTATION_GTFS`, with its per-chromosome genes and exon models and its gene lookup, SPEC.md
  section 6), a variant catalog `<id>_grch38` and the results (with the search index split by
  chromosome and the page counts, SPEC.md section 8, and the trans objects when
  the tables hold `trans.parquet`, and the GWAS object when they hold `gwas.parquet`, with its bin
  summary when they also hold `gwas_bins.parquet`); then
  `Store.validate`, a per-experiment count of genes missing from the annotation, the cross-catalog
  lookup check when `CROSSCAT` names two variant catalogs, and the v0 comparison for `topchef`
  (`python -m pipeline.verify_v0`, which also compares the GWAS bin summary to v0's
  `gwas_dcm_bins.json` field by field; `SKIP_V0=1` skips it). Each phase logs its wall time and peak RSS. `results.py` holds one chromosome's nominal rows
  at a time as compact arrays (DuckDB join, `QTLB_DUCKDB_MEMORY`, default 12GB).
- **Benchmark**: bytes per table type (variants, eQTL, sQTL, hits, trans, GWAS), warm random-access
  read medians/p95 for gene and intron blocks, the gene page's trans table and GWAS window, and startup cost;
  writes `bench_store.{json,md}` under `/scratch/ns5bc/qtl-browser/bench/`. Results are kept in the
  `results_analysis/qtlb_format` brick under `data/` (genome-wide: `data/format_v0_v1_2026-09-24/`).
  The browser side has its own smoke suites and gene-page benchmark in `ui/bench/`.

### Experiment modularity

Every object is named by its bytes, so experiments can come and go without touching each other:

```bash
uv run python -m pipeline.qtlstore remove-experiment gtex_v8_heart_lv --store <store>  # drop the pointer, rewrite store.json
uv run python -m pipeline.qtlstore gc --store <store> [--dry-run]                      # delete objects no pointer names
uv run python -m pipeline.qtlstore validate --store <store>
```

`remove-experiment` deletes only `experiments/<id>.json` and rewrites `store.json`; the experiment's
variant catalog and annotation pointers stay (another experiment may use them; delete those pointer
files by hand first if they should go too), and its objects stay until `gc`. `gc` deletes every
`immutable/` object that no pointer file on disk names, plus leftover `*.tmp` files.
`test_results.py::test_remove_experiment_and_gc` covers the cycle.

Checked on chr21/22 scratch stores on 2026-09-24 (`store-c2122-mod`, Slurm job 20443347): build
TOPCHeF alone (18 objects), add GTEx (31), remove GTEx (`validate` PASS), `gc` (7 GTEx experiment
objects deleted, 28.8 MB; the GTEx variant catalog and annotation pointers stay and keep theirs),
`validate` PASS again, and every TOPCHeF object (its experiment, variant catalog and annotation) has
the same SHA-256 after each step. A second `gc` deletes nothing.

Paths, the significance rule, worker counts, the eQTL Catalogue experiments and the gate's sample
genes live in `config.yaml`. Nothing is hard-coded in the steps.

The `reference:` block in `config.yaml` names the reference FASTA, the local refgetstore, and
`reference.collection`, the seqcol digest that pins the assembly the positions sit on.
The `refget_store` step and the adapters need the `gtars` package, which `uv sync` installs. On a
machine with no store yet `refget_store` builds one from the FASTA and prints the collection digest to paste into `reference.collection`;
after that the digest is required to match.

## Input tables (`python -m pipeline build`)

The TOPCHeF adapter reads the variant table and the reference allele check from `_tables/`, and its
acceptance gate compares against the v0 tables there. These steps build them. Each leaves a marker
in `data/derived/.done/` and is skipped on the next run unless `--force` is given. On Rivanna the
frozen v0 tree (`/scratch/ns5bc/qtl-browser/derived/`) already holds them, and an adapter tree
symlinks its `_tables` inputs.

| Step | Module | Reads | Writes |
|---|---|---|---|
| extract | `steps_extract` | Zenodo `*.tar.gz` | per-chromosome parquet unpacked next to the archives |
| gtf | `steps_gtf` | GENCODE v34 GTF | `_tables/gene_annotation.parquet`, `_tables/exons.parquet` (sorted by gene, small row groups) |
| variants_collect | `steps_variants` | every cis file, then both trans files | `_tmp/variants_raw.parquet`: distinct (chr, position, A1, A2) with `in_cis`. 8.87M cis variants plus 343k positions seen only in the genome-wide trans scan; trans_eQTL has no allele columns, so its positions get null alleles unless trans_sQTL has them |
| variants_rsid | `steps_variants` | dbSNP b157 VCF via `bcftools query -T` on those positions | `_tables/variants.parquet`, sorted by (chr, position, A1, A2), with an exact / position / none match flag (allele-less rows can only match by position) and `rsid` always `'rs' || rs_number` (`arg_min`, not two separate `min`s). The two `_tmp/dbsnp_*` caches are rebuilt when older than `variants_raw.parquet` |
| refget_store | `steps_refget` | the reference FASTA named by `reference.fasta`, or an existing refgetstore | `_tables/reference.json`: the seqcol collection digest, its three attribute digests, and every chromosome's sha512t24u sequence digest, length and md5. Builds the local store from the FASTA when there is none, and fails when a name in `reference.chromosomes` is not in the collection |
| variants_refcheck | `steps_refget` | `_tables/variants.parquet`, `_tables/reference.json`, the refgetstore | `_tables/refcheck/<chr>.parquet`, one row per variant saying which allele the reference carries (`a1`, `a2`, `both`, `neither`, `unchecked`), and `_checks/refcheck/summary.json`, the counts per class and chromosome split cis against trans-only, with strand-flip and on-N diagnostics and example mismatches (EVIDENCE.md A.10). Fails when the cis match fraction is below `reference.allele_check.min_match_fraction` |
| permutation_tables | `steps_tables` | cis permutation files, SuSiE, trans, annotation | `_tables/genes.parquet`, `_tables/splice_phenotypes.parquet` |
| credible_sets | `steps_tables` | SuSiE files | `_tables/credible_sets.parquet` |
| nominal | `steps_nominal` | cis nominal files, one process per chromosome | `_tables/cis_eqtl_nominal/chr=*/bin=*/`, `_tables/cis_sqtl_nominal/chr=*/bin=*/`: one file per `nominal_bin_genes` tested genes (by TSS rank, the `bin` column of `genes`), one row group per gene, delta/byte-stream-split encodings, rsIDs as `rs_number`. With `sqtl_nominal: significant` the sQTL side keeps only introns flagged `is_sqtl` |

## Frozen v0 and its reader

The v0 pack build is gone (it is in git history at commit `c62bca3`). The frozen v0 build on Rivanna
(`/scratch/ns5bc/qtl-browser/derived/`: `immutable/`, `manifest.json`, `_tables/`) stays as the
reference the v1 checks compare against. The v0 format is described in the `analysis` repo at
`qtlb-format/docs/` (`PACKS.md`, historical; the v0 `SPEC.md` at commit `fe5a606`; `EVIDENCE.md`).
What reads it here:

- `packfmt_v0.py`: the v0 codec. `verify_v0.py`, `bench_store.py` and the TOPCHeF gate decode v0
  packs with it. No v1 builder imports it (`test_qtlstore.py` checks).
- `packtool.py`: a CLI and API over `packfmt_v0` for looking inside v0 packs. The checks use its
  readers (`read_block`, `find_gene_block`, `list_blocks`, `variant_rows`).
- `test_packfmt.py`, `test_packtool.py`: round-trip tests on synthetic data; `test_packfmt` without
  `--synthetic` also runs a real-data case that needs the v0 tables.

Subset runs: `QTLB_CHROMS=chr21,chr22` narrows the chromosome list and `QTLB_DERIVED` moves the
whole output tree, step markers included. Always set both for a subset run, and never point one at
the frozen v0 tree.

## Conventions

- `gene_id` is an unversioned ENSG; `symbol` is the GENCODE name; `chr` is `chr1`..`chrX`.
- `A1` is the effect (minor) allele, `A2` the reference, as in the Zenodo release. Verified
  against GRCh38: all 8,419,594 cis SNPs read A2 at their position, none read A1.
- sQTL `phenotype_id` is the leafcutter string `chr:start:end:clu_N_strand:ENSG.v`.
- `pval_perm < 0.05` is the significance flag; a BH `qval` on `pval_beta` is stored alongside.
