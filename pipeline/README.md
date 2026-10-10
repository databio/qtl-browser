# pipeline

Builds the qtlstore (`SPEC.md` at the repo root, format version 2) that the browser reads. Each
study goes through an adapter into contract tables (`CONTRACT.md`), then into a refget-anchored
store, then through checks, then into the bucket. The adapters run on Rivanna through Slurm from a
checkout (`QTLB_REPO`, default `/scratch/$USER/qtl-browser`; logs go to `/scratch/<user>/`); the
store builds anywhere the tables are.

```bash
# adapters (Slurm; details below)
QTLB_DERIVED=<tree> sbatch adapter.sbatch                          # TOPCHeF adapter + contract check

# the store: one command from the configured TOPCHeF tables to a validated store (about 10 min
# genome-wide; without --keep it deletes the directory first, so never aim it at a combined store)
uv run python -m pipeline.build_store --store data/derived/store-topchef-v2 --validate
# a second cohort into an existing store (how store-release-v2 got GTEx; log build_15_store_release_v2.log)
uv run python -m pipeline.annotation build --gtf data/raw/gencode_v39/gencode.v39.annotation.gtf.gz --store <store> --id gencode_v39
uv run python -m pipeline.catalog --sites <tables>/sites.parquet --store <store> --id gtex_v8_heart_lv_grch38
uv run python -m pipeline.results build --tables <tables> --store <store> --id gtex_v8_heart_lv --catalog gtex_v8_heart_lv_grch38 --annotation gencode_v39
uv run python -m pipeline.overlap build --store <store> --id grch38 && uv run python -m pipeline.coloc add --store <store>
# then rewrite store.json and validate, as store.sbatch's validate step does
uv run python -m pipeline.which_codec --store <store>              # which codec wrote each object, from its bytes
uv run python -m pipeline.scan_codes --store <store> --experiment topchef_ge_cis [--chrom chr22]
uv run python -m pipeline.qtlstore validate | remove-experiment ID | gc [--dry-run] | crosscat A B --store <store>
uv run python -m pipeline.overlap build --store <store> --id <id> --catalogs A B ...   # catalog overlap index
QTLB_STORE=<store> EXPERIMENTS="topchef:<tree>/_tables/topchef" sbatch store.sbatch   # the same on Rivanna, several cohorts

# the bucket (dry runs unless --yes; credentials from .env or the environment)
uv run python -m pipeline.upload plan --store <store>
uv run python -m pipeline.upload push --store <store> --yes && uv run python -m pipeline.upload check --store <store>
uv run python -m pipeline.upload prune --store <store> [--yes]

# the shared input tables the TOPCHeF adapter reads (build.sbatch)
uv run python -m pipeline steps                                    # list steps in order
uv run python -m pipeline build                                    # run what is not done yet
uv run python -m pipeline build --step variants_refcheck --force

# tests (no data needed)
uv run python -m pipeline.test_qtlstore      # also test_packfmt_v2, test_catalog, test_annotation, test_results, test_upload,
                                             # test_dof, test_genemap, test_refcheck, test_overlap,
                                             # adapters.test_topchef, adapters.test_eqtl_catalogue
```

## The build: adapters, store, checks, upload

Each study goes through an adapter into contract tables (`CONTRACT.md`), then into a
refget-anchored store (`qtlstore.py`): one variant catalog per cohort, one experiment per scan (a
cis scan with its results blocks, search index and hits; a trans scan with its trans object, index
and hits; a GWAS with its blocks, index and bin summary), every experiment naming the collection it
was published in, sharing annotations. The builders (`qtlstore`, `catalog`, `annotation`,
`results`, `gwas`) write version 2 unless told `--version 1`, take the codec as a parameter and
import `packfmt_v1` only for the pieces every version shares. Each adapter writes into its own `QTLB_DERIVED` tree; nothing here writes into the
frozen v0 tree (`adapter.sbatch` refuses it). Measure on a genome-wide build: a chr21/22 subset
changes the shape of the trans data, not just its size.

| Step | Command | Module |
|---|---|---|
| TOPCHeF adapter + contract check | `QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/topchef.py`, `adapters/contract_check.py` |
| eQTL Catalogue adapter + contract check | `ADAPTER=eqtl_catalogue QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/eqtl_catalogue.py` |
| One step alone | `STEPS=adapter` or `STEPS=contract` with the same variables | |
| One TOPCHeF table | `STEPS=adapter QTLB_DERIVED=<tree> sbatch adapter.sbatch --step phenotypes --step trans` (trans needs the trans-only phenotypes, so both) | `adapters/topchef.py` |
| DCM GWAS into the TOPCHeF tables | `ADAPTER=dcm_gwas STEPS="adapter contract" QTLB_DERIVED=<tree> sbatch adapter.sbatch` | `adapters/dcm_gwas.py` |
| Store build, one cohort, locally | `uv run python -m pipeline.build_store --store <store> --validate` (`--tables`, `--sites`, `--gtf` override the configured paths; `--keep` adds to an existing store) | `build_store.py` |
| Store build, several cohorts, Rivanna | `QTLB_CHROMS=all QTLB_STORE=<store> EXPERIMENTS="topchef:<tree>/_tables/topchef:gencode_v34 gtex_v8_heart_lv:<tree>/_tables/gtex_v8_heart_lv" CROSSCAT="topchef_grch38 gtex_v8_heart_lv_grch38" sbatch --time=6:00:00 --mem=64G store.sbatch` | `annotation.py`, `catalog.py`, `results.py`, `gwas.py`, `which_codec.py` |
| Store maintenance | `uv run python -m pipeline.qtlstore validate \| remove-experiment ID \| gc [--dry-run] \| crosscat A B --store <store>` | `qtlstore.py` |
| Which codec wrote each object | `uv run python -m pipeline.which_codec --store <store>`: reads the block and frame magics back out of the bytes and judges them against `store.json`'s declared version | `which_codec.py` |
| Upload | `uv run python -m pipeline.upload inventory \| plan \| push [--yes] \| check \| prune [--yes] --store <store>` | `upload.py` |
| Catalog overlap index | `uv run python -m pipeline.overlap build --store <store> --id <id> [--catalogs A B ...]` (run by `store.sbatch`; `OVERLAP=` skips it) | `overlap.py` |
| Gene symbols -> ENSG, for a source that publishes no gene id | `uv run python -m pipeline.genemap --symbols aric` (the report for ARIC's `seqid.txt`; any file of one symbol per line otherwise) | `genemap.py` |
| Search index column rename, in place | `uv run python -m pipeline.results migrate-index --store <store> --id <experiment>` (25 objects per experiment; no results, variant, hits or GWAS object changes) | `results.py` |
| Results-block code sweep | `uv run python -m pipeline.scan_codes --store <store> --experiment <id> [--chrom chr22 ...]`, or `--base <url>` to read a published store over HTTP | `scan_codes.py` |
| Per-chromosome objects for a store built before them | `uv run python -m pipeline.annotation add-split --store <store> --id <annotation>`, then `uv run python -m pipeline.results add-split --store <store> --id <experiment>` (rewrites the pointer; no other object changes) | `annotation.py`, `results.py` |

Outside Slurm the same modules run directly, e.g.
`uv run python -m pipeline.adapters.contract_check $QTLB_DERIVED/_tables/topchef` (prints one
PASS/FAIL line per `CONTRACT.md` rule, exits 1 on any failure).

- **Store build** (`build_store.py`, or `store.sbatch` per cohort): an annotation (named in
  `EXPERIMENTS`, else the adapter's `ingestion.json` `source.gene_annotation`, else `gencode_v34`;
  GTFs from `ANNOTATION_GTFS`, with its per-chromosome genes and exon models and its gene lookup,
  SPEC.md section 6), a variant catalog `<id>_grch38`, and one experiment per scan: a cis and a
  trans experiment per phenotype type (the trans ones when the tables hold `trans.parquet`), each
  with its search index split by chromosome and its hits, and a GWAS experiment when they hold
  `gwas.parquet`, with the bin summary the builder computes from its rows (SPEC.md sections 8
  to 10); then `Store.validate`, a per-experiment count of genes missing from the annotation, the
  cross-catalog lookup check when `CROSSCAT` names two variant catalogs, and `which_codec`. Each
  phase logs its wall time and peak RSS. `results.py` holds one chromosome's nominal rows at a time
  as compact arrays (DuckDB join, `QTLB_DUCKDB_MEMORY`, default 12GB).
- **Upload** (`upload.py`): objects first, then pointers, then `store.json` last, so the live store
  switches in one PUT and a reader mid-session keeps working; `prune` is separate and refuses to
  run until the public `store.json` is the local one. The browser side has its own smoke suites and
  gene-page benchmark in `ui-topchef/bench/`.

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

Checked end to end on a chr21/22 store: build TOPCHeF alone (18 objects), add GTEx (31), remove
GTEx, `gc` (7 objects, 28.8 MB; GTEx's variant catalog and annotation pointers stay and keep
theirs), with `validate` passing at each step and every TOPCHeF object keeping the same SHA-256
throughout. A second `gc` deletes nothing.

### Catalog overlap index

`overlap.py` (SPEC.md section 19) records, for every site any of the store's variant catalogs holds,
which of them hold it: membership in one read, a catalog's vidx for a shared site, and exact
pairwise overlap counts written into the pointer, so the statistics need no object read at all. It is
derived: rebuild it any time, and a store without one is still complete. `store.sbatch` rebuilds it
before `validate`, which fails on an index whose catalogs no longer carry the identity it was
built against.

The pointer also counts positions two catalogs hold with no shared allele pair, and indels that may
be one event written two ways. The second kind are candidates to check rather than
conclusions: confirming one needs the reference sequence, and `CONTRACT.md` does not require
left-aligned indels.

### Results-block code sweep

`validate` checks object digests, headers and pointers; it never reads the codes inside a block.
`scan_codes.py` sweeps those directly, for one cis scan -- fast, because a block stores them raw --
picking the codec from each object's own header, and counts the reserved values: `p_null` (65535,
the phenotype was never tested against that variant), `p_zero` (65534, p underflowed to 0 in the
source, so no z exists and a reader must drop the row), `se_null`, `beta_null`, and `both`. Those
are legitimate data. The failures it catches are the three scale rules (a block's largest finite
code must sit exactly at the top of its ruler, or the codes and the header's scale disagree and
every value in the block decodes wrong) and a `measured_worst` above its `theoretical_limit` in the
pointer. It exits 1 on any.

It is not part of the build: the full sweep reads every results object (3.8 GB for TOPCHeF), which
is too much for a routine validate. Run it after a format change, or to answer a question about a
published store.

TOPCHeF, genome-wide, both cis scans: 623,057,691 rows, **no underflowed p at all**, 2,649 untested
rows (1 in 235,000), null in all three codes together, and every scale rule holding in all 100,173
blocks. The row count reconciles to the pointers: 623,057,691 - 2,649 equals `rows_compared`
summed over the two cis scans.

### Configuration

Paths, the significance rule, worker counts, the experiments of each adapter, the collections and
the coloc lists live in `config.yaml`. Nothing is hard-coded in the steps.

Subset runs: `QTLB_CHROMS=chr21,chr22` narrows the chromosome list and `QTLB_DERIVED` moves the
whole output tree, step markers included. Always set both, and never point either at the frozen v0
tree.

The `reference:` block in `config.yaml` names the reference FASTA, the local refgetstore, and
`reference.collection`, the seqcol digest that pins the assembly the positions sit on.
The `refget_store` step and the adapters need the `gtars` package, which `uv sync` installs. On a
machine with no store yet `refget_store` builds one from the FASTA and prints the collection digest to paste into `reference.collection`;
after that the digest is required to match.

## Input tables (`python -m pipeline build`)

The TOPCHeF adapter reads the variant table and the reference allele check from `_tables/`. These
steps build them. Each leaves a marker in `data/derived/.done/` and is skipped on the next run unless
`--force` is given. On Rivanna the frozen v0 tree (`/scratch/ns5bc/qtl-browser/derived/`) already
holds them.

| Step | Module | Reads | Writes |
|---|---|---|---|
| extract | `steps_extract` | Zenodo `*.tar.gz` | per-chromosome parquet unpacked next to the archives |
| variants_collect | `steps_variants` | every cis file, then both trans files | `_tmp/variants_raw.parquet`: distinct (chr, position, A1, A2) with `in_cis`. 8.87M cis variants plus 343k positions seen only in the genome-wide trans scan; trans_eQTL has no allele columns, so its positions get null alleles unless trans_sQTL has them |
| variants_rsid | `steps_variants` | dbSNP b157 VCF via `bcftools query -T` on those positions | `_tables/variants.parquet`, sorted by (chr, position, A1, A2), with an exact / position / none match flag (allele-less rows can only match by position) and `rsid` always `'rs' || rs_number` (`arg_min`, not two separate `min`s). The two `_tmp/dbsnp_*` caches are rebuilt when older than `variants_raw.parquet` |
| refget_store | `steps_refget` | the reference FASTA named by `reference.fasta`, or an existing refgetstore | `_tables/reference.json`: the seqcol collection digest, its three attribute digests, and every chromosome's sha512t24u sequence digest, length and md5. Builds the local store from the FASTA when there is none, and fails when a name in `reference.chromosomes` is not in the collection |
| variants_refcheck | `steps_refget` | `_tables/variants.parquet`, `_tables/reference.json`, the refgetstore | `_tables/refcheck/<chr>.parquet`, one row per variant saying which allele the reference carries (`a1`, `a2`, `both`, `neither`, `unchecked`), and `_checks/refcheck/summary.json`, the counts per class and chromosome split cis against trans-only, with strand-flip and on-N diagnostics and example mismatches (EVIDENCE.md A.10). Fails when the cis match fraction is below `reference.allele_check.min_match_fraction` |

## Earlier formats

Version 1 is in git history at commit `be0bb44` (its spec, `git show be0bb44:SPEC.md`) and in the
bucket until the version 2 release. `packfmt_v1.py` still reads it: `results.read_block`,
`results.read_trans` and `gwas.read_window` pick the codec from each object's header, and the
builders write it with `--version 1` for an object-for-object comparison.

v0 is retired. The pack build is in git history at commit `c62bca3`; its codec (`packfmt_v0.py`),
reader (`packtool.py`), the TOPCHeF acceptance gate that compared the contract tables against it
(`adapters/verify_topchef.py`), the steps that wrote its tables (`steps_gtf`, `steps_tables`,
`steps_nominal`) and `figures.py` are in history at `8e1fb44`. The frozen v0 build on Rivanna
(`/scratch/ns5bc/qtl-browser/derived/`) still holds the shared input tables; nothing reads its packs.
The v0 format is described in the `analysis` repo at `qtlb-format/docs/` (`PACKS.md`, historical;
the v0 `SPEC.md` at commit `fe5a606`; `EVIDENCE.md`).

The v0-against-v1 comparison and benchmark (`verify_v0.py`, `bench_store.py`, `bench_store.sbatch`)
are in git history before 2026-10-09; they read the version 1 single-document store and nothing
current.

## Conventions in the input tables

These are the `_tables/` conventions, which follow the Zenodo release. The store itself is
ref/alt with ALT-relative effects (`CONTRACT.md`, SPEC.md section 7), so an adapter converts.

- `gene_id` is an unversioned ENSG; `symbol` is the GENCODE name; `chr` is `chr1`..`chrX`.
- `A1` is the effect (minor) allele, `A2` the reference. Verified against GRCh38: all 8,419,594 cis
  SNPs read A2 at their position, none read A1, so the TOPCHeF adapter maps A2 to `ref` and keeps
  the sign.
- sQTL `phenotype_id` is the leafcutter string `chr:start:end:clu_N_strand:ENSG.v`.
- `pval_perm < 0.05` is the significance flag; a BH `qval` on `pval_beta` is stored alongside.
