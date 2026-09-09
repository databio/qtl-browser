# pipeline

Turns `data/raw/` into the parquet tables the browser reads, under `data/derived/`, and writes
`manifest.json` describing them. Idempotent: each step leaves a marker in `data/derived/.done/`
and is skipped on the next run unless `--force` is given.

```bash
uv run python -m pipeline steps                       # list steps in order
uv run python -m pipeline build                       # run what is not done yet
uv run python -m pipeline build --step nominal --force
uv run python -m pipeline validate
uv run python -m pipeline.figures manhattan | density 5 10 | gwas 5 | themes
uv run python -m pipeline.upload cors | sync [--dryrun] [--delete] | check   # R2, see config `r2:`
```

Paths, the significance rule, window sizes, worker counts, and the paper's reference counts
live in `config.yaml`. Nothing is hard-coded in the steps.

## Steps

| Step | Module | Reads | Writes |
|---|---|---|---|
| extract | `steps_extract` | Zenodo `*.tar.gz` | per-chromosome parquet unpacked next to the archives |
| gtf | `steps_gtf` | GENCODE v34 GTF | `gene_annotation.parquet`, `exons.parquet` (sorted by gene, small row groups) |
| variants_collect | `steps_variants` | every cis file, then both trans files | `_tmp/variants_raw.parquet`: distinct (chr, position, A1, A2) with `in_cis`. 8.87M cis variants plus 343k positions seen only in the genome-wide trans scan; trans_eQTL has no allele columns, so its positions get null alleles unless trans_sQTL has them |
| variants_rsid | `steps_variants` | dbSNP b157 VCF via `bcftools query -T` on those positions | `variants_by_position/chr=*/`, `variants_by_rsid.parquet`, with an exact / position / none match flag (allele-less rows can only match by position). The two `_tmp/dbsnp_*` caches are rebuilt when older than `variants_raw.parquet` |
| permutation_tables | `steps_tables` | cis permutation files, SuSiE, trans, annotation | `genes.parquet`, `splice_phenotypes.parquet` (sorted chr, tss; stats for range reads), `search_index.parquet` (the one table the browser loads whole) |
| credible_sets | `steps_tables` | SuSiE files | `credible_sets.parquet` |
| nominal | `steps_nominal` | cis nominal files, one process per chromosome | `cis_eqtl_nominal/chr=*/bin=*/`, `cis_sqtl_nominal/chr=*/bin=*/`: one file per `nominal_bin_genes` tested genes (by TSS rank, the `bin` column of `genes`), one row group per gene, delta/byte-stream-split encodings, rsIDs as `rs_number`. With `sqtl_nominal: significant` the sQTL side keeps only introns flagged `is_sqtl` |
| gene_detail | `steps_tables` | genes, exons, splice_phenotypes | `gene_detail/chr=*/bin=*/`: one row per tested gene with the genes row, the collapsed exon model (`exons` list) and every tested intron (`splice` list); one row group per gene |
| trans | `steps_tables` | trans files | `trans_pairs/chr=<gene chr>/` sorted by gene, and `trans_by_variant/chr=<variant chr>/` sorted by position |
| coloc_stub | `steps_tables` | `coloc_genes` in the config, genes | empty `coloc.parquet` with the intended schema; `coloc_loci.json`, the landing track's 25 loci (one plain fetch, no engine) |
| gwas_bins | `steps_gwas` | the Jurgens 2024 file named by `dcm_gwas` in the config (biobanks-only; the CVDKP zip has five) | `gwas_dcm_bins.parquet` and the same rows as `gwas_dcm_bins.json`: strongest p per 5 Mb window |
| gwas_full | `steps_gwas` | same | `gwas_dcm/chr=*/` position-sorted in 10k-row groups (about 3 Mb each), for LocusCompare; `gwas_dcm.json` with the file name, cases, controls, and variant count for the manifest |
| manifest | `steps_finish` | everything above | `manifest.json`: tables, rows, sizes, columns, counts, source versions |

Underscore directories in `data/derived/` (`_tmp`, `_full`) are scratch and never uploaded;
`_full/cis_sqtl_nominal` is the all-introns build parked for local use.

`validate` checks row counts against the raw files, eGene and sQTL counts against the preprint,
that six paper-named variants resolve to their rsIDs, that a FLNC query touches one row group,
and the rsID exact-match rate.

## Layout rules the browser depends on

- Big tables are hive-partitioned by chromosome, and the per-gene ones (nominal, gene_detail)
  by chromosome and `bin`, so a footer covers about 100 genes instead of a chromosome. The app
  always names the exact partition file, since there is no directory listing over HTTP; it
  takes `bin` from the search index.
- Row-group statistics are written only for the columns queries filter on. A per-gene read
  needs `gene_id` (nominal, trans, exons) or `chr` + `tss` (genes, splice phenotypes, credible
  sets: the app knows both from the search index). Position lookups need `position`.
- Sort order decides what a row group's statistics cover. Sorting by TSS instead of gene id made
  a gene's trans read the whole chromosome; sorting exons in GTF order made one gene touch
  sixteen row groups. Both were fixed by sorting on the filtered column.
- Concurrent DuckDB workers must not share a `temp_directory`; `connect()` takes one per worker.

## Conventions

- `gene_id` is an unversioned ENSG; `symbol` is the GENCODE name; `chr` is `chr1`..`chrX`.
- `A1` is the effect (minor) allele, `A2` the reference, as in the Zenodo release.
- sQTL `phenotype_id` is the leafcutter string `chr:start:end:clu_N_strand:ENSG.v`.
- `pval_perm < 0.05` is the significance flag; a BH `qval` on `pval_beta` is stored alongside.
