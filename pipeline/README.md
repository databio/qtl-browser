# pipeline

Turns `data/raw/` into the binary pack files the browser reads. Every one of them is published to
`data/derived/immutable/` under a content-addressed name, and `manifest.json` is what names them.
Idempotent: each step leaves a marker in `data/derived/.done/` and is skipped on the next run unless
`--force` is given.

```bash
uv run python -m pipeline steps                       # list steps in order
uv run python -m pipeline build                       # run what is not done yet
uv run python -m pipeline build --step nominal --force
uv run python -m pipeline validate
uv run python -m pipeline packcheck dof                       # pack format checks (SPEC.md); not a build step, see below
uv run python -m pipeline packcheck run [--type e|s] [--chrom chr7 ...] [--workers 3]
uv run python -m pipeline packcheck report
uv run python -m pipeline packcheck measure --chrom chr7 chr10 chr22
uv run python -m pipeline packcheck roundtrip [--chrom chr7 ...] [--steepest 50] [--workers 2]
uv run python -m pipeline.packtool header | blocks | block | variants | frames | trans | gwas | gwas-index | check ...   # inspect packs, decode rows
uv run python -m pipeline.packtool pack-variants | pack-results | pack-trans | pack-gwas ...                           # build packs from tables
uv run python -m pipeline.test_packfmt --synthetic && uv run python -m pipeline.test_packtool          # codec and tool round-trip tests
uv run python -m pipeline.packtool_example                                                             # worked example on synthetic tables
uv run python -m pipeline.figures manhattan | density 5 10 | gwas 5 | themes
uv run python -m pipeline.upload inventory [--public] -o FILE          # R2 (config `r2:`); read-only listing, reusable as --listing
uv run python -m pipeline.upload budget [--listing FILE] [--harness JSON] [--overlap-days N]
uv run python -m pipeline.upload stage [--dryrun] [--listing FILE] | check [--staged|--retired] | release [--dryrun] [--listing FILE]
uv run python -m pipeline.upload pointer FILE | prune [--stale] [--yes] [--listing FILE] | cors [--print]
uv run python -m pipeline.test_upload                                  # upload.py against a fake bucket; no network
```

Paths, the significance rule, window sizes, worker counts, and the paper's reference counts
live in `config.yaml`. Nothing is hard-coded in the steps.

## Steps

| Step | Module | Reads | Writes |
|---|---|---|---|
| extract | `steps_extract` | Zenodo `*.tar.gz` | per-chromosome parquet unpacked next to the archives |
| gtf | `steps_gtf` | GENCODE v34 GTF | `_tables/gene_annotation.parquet`, `_tables/exons.parquet` (sorted by gene, small row groups) |
| variants_collect | `steps_variants` | every cis file, then both trans files | `_tmp/variants_raw.parquet`: distinct (chr, position, A1, A2) with `in_cis`. 8.87M cis variants plus 343k positions seen only in the genome-wide trans scan; trans_eQTL has no allele columns, so its positions get null alleles unless trans_sQTL has them |
| variants_rsid | `steps_variants` | dbSNP b157 VCF via `bcftools query -T` on those positions | `_tables/variants.parquet`, sorted by (chr, position, A1, A2), with an exact / position / none match flag (allele-less rows can only match by position) and `rsid` always `'rs' || rs_number` (`arg_min`, not two separate `min`s). The two `_tmp/dbsnp_*` caches are rebuilt when older than `variants_raw.parquet` |
| permutation_tables | `steps_tables` | cis permutation files, SuSiE, trans, annotation | `_tables/genes.parquet`, `_tables/splice_phenotypes.parquet` |
| credible_sets | `steps_tables` | SuSiE files | `_tables/credible_sets.parquet` |
| nominal | `steps_nominal` | cis nominal files, one process per chromosome | `_tables/cis_eqtl_nominal/chr=*/bin=*/`, `_tables/cis_sqtl_nominal/chr=*/bin=*/`: one file per `nominal_bin_genes` tested genes (by TSS rank, the `bin` column of `genes`), one row group per gene, delta/byte-stream-split encodings, rsIDs as `rs_number`. With `sqtl_nominal: significant` the sQTL side keeps only introns flagged `is_sqtl` |
| pack_sqtl | `steps_pack` | raw cis sQTL nominal (every intron), `_tables/variants.parquet`, splice_phenotypes, credible_sets | all 23 sQTL packs, published to `immutable/sqtl.<chr>.<sha16>.qbs`: one block per tested intron (80,750), streamed in the raw file's order without sorting, one process per chromosome (`packs.sqtl_workers`, each with a small DuckDB). Plus `_tmp/pack_pointers/sqtl_<chr>.{parquet,json}`: each intron's block and variant run, and per-chromosome counts, the dof check, and peak RSS. Fails when an intron's rows are not one run of the variants list, or a count differs from its source |
| pack_eqtl | `steps_pack` | `_tables/variants.parquet`, cis eQTL nominal, raw cis sQTL nominal (af and counts for variants tested only for sQTL), the sQTL pointer files, genes, exons, splice_phenotypes, credible_sets | all 23 variants files and 23 eQTL packs, published to `immutable/variants.<chr>.<sha16>.qbv` and `immutable/eqtl.<chr>.<sha16>.qbe`, plus `_tmp/pack_pointers/eqtl_<chr>.{parquet,json}`. Each gene's details give every intron's `blk_off`/`blk_len` in the sQTL pack, and the union variant range takes intron runs from the sQTL pointers, so it runs after `pack_sqtl` |
| trans | `steps_tables` | trans files | `_tables/trans/chr=<gene chr>/` sorted by gene (read by `pack_trans`), and `_tables/trans_by_variant/chr=<variant chr>/` sorted by position (read by `pack_hits`; one file per chromosome there instead of filtering all 24, which is why that step takes 1.4 minutes). Fails if any trans phenotype has no gene chromosome |
| pack_variants_trans | `steps_pack_trans` | `_tables/trans`, `_tables/variants.parquet`, the variants files | republishes every variants file with a trans-only section after the cis pages (SPEC.md section 4): the `NOT in_cis` variants by position, af from the trans rows (one value per variant, checked), no sample counts, flags bit 2 where the source reports no alleles. The bytes changed, so the file gets a new `immutable/variants.<chr>.<sha16>.qbv` name; it refuses to publish one whose cis bytes would change, so `var_off`/`var_len` and every eQTL and sQTL run stay valid and only `search_index` has to be rebuilt after it. Also `_tmp/trans_variant_af.parquet` and `_tmp/pack_pointers/variants_trans_<chr>.json` |
| pack_trans | `steps_pack_trans` | `_tables/trans`, `_tables/variants.parquet`, genes | `immutable/trans.<gene chr>.<sha16>.qbt` for chr1-22, chrX, chrM (SPEC.md section 12): one zstd frame per gene with trans rows, in `search_index` order, with variant chromosome, position, rs_number, af, p, and beta inline. Plus `_tmp/pack_pointers/trans_<chr>.{parquet,json}` (each gene's `trans_off`/`trans_len`, counts, frame sizes) |
| search_index | `steps_pack` | genes, splice_phenotypes, pack pointer files | `immutable/search_index.<sha16>.arrow.zst` (SPEC.md section 6): an Arrow IPC stream in one zstd frame, with block, union-range, GWAS-window and trans-frame pointers, `gene_version` and `ord`. Its schema metadata records the SHA-256 of every pack its offsets reach, so a stale pairing cannot reach the bucket. No `bin`: that named a nominal partition, which only the builders read |
| pack_hits | `steps_pack_variant` | `_tables/trans_by_variant`, genes, splice_phenotypes, credible_sets, `_tables/variants.parquet`, the variants files, search_index | all 23 hits packs, published to `immutable/hits.<chr>.<sha16>.qbh` (SPEC.md section 13): one zstd frame per `packs.hits_frame_variants` variant indices, over both sections of the variants file, holding each variant's trans eQTL and sQTL rows, the phenotypes it is the lead of, and its credible-set memberships, each carrying a `search_index` `ord` instead of a gene id. Builds `_tmp/variant_idx.parquet` first, the (chr, position) -> `vidx` map, and checks every position in it against the decoded variants file. Fails when a source row does not join to a variant or to a gene, or when the rows by kind differ from the sources. Plus `_tmp/pack_pointers/hits_<chr>.{parquet,json}` (each frame's `hits_off`/`hits_len`, counts by kind, frame sizes) |
| pack_variant_index | `steps_pack_variant` | `_tmp/variant_idx.parquet`, the variants files, the hits packs | `immutable/rsid_index.<sha16>.qbr` (SPEC.md section 14): every rsID as an uncompressed 8-byte record sorted by rs number, in blocks of `packs.rsid_block_records`, so a lookup is one range request at a computed offset. And `immutable/variant_index.<sha16>.qbx` (section 15), the startup file the variant page needs: every variants-file page's byte offset and first position, every hits frame's offset, and each rsID block's first number. Rebuild it whenever a variants file, a hits pack, or the rsID index changes. Plus `_tmp/pack_pointers/variant_index.{parquet,json}` |
| coloc_stub | `steps_tables` | `coloc_genes` in the config, genes | `coloc_loci.json`, the landing track's 25 loci (one plain fetch, no engine) |
| gwas_bins | `steps_gwas` | the Jurgens 2024 file named by `dcm_gwas` in the config (biobanks-only; the CVDKP zip has five) | `gwas_dcm_bins.json`: strongest p per 5 Mb window, columnar, one plain fetch |
| pack_gwas | `steps_gwas` | same | `immutable/gwas.<chr>.<sha16>.qbg` (chr1-22) and `immutable/gwas_index.<sha16>.bin` (SPEC.md section 11): every row, lossless at the source's 4 decimals and 4 significant digits, in zstd blocks of `packs.gwas_block_rows` rows, and the startup index that turns a gene's `[w_lo, w_hi]` into one byte range. One DuckDB read of the TSV (`packs.gwas_duckdb_*`); fails on any row that breaks the lossless rules. Also `gwas_dcm.json` (file, cases, controls, variants) for the manifest and `_tmp/pack_pointers/gwas.json` (counts, sizes, window byte ranges). Reads `search_index` for the window sizes, so it runs after that step |
| manifest | `steps_finish` | everything above | `manifest.json`: counts, source versions, the `packs` block (the published path, bytes and counts of every pack kind, the search index, the GWAS index, the rsID index, the variant index, dof), the `immutable` block (per file: logical key, bytes, sha256, md5, and the old bucket prefixes it `replaces`), and the `precision` block of rounding maximums. Refuses to write when `search_index` was built from packs that have since been rebuilt. No `tables` block: every table is a build intermediate |

Underscore directories in `data/derived/` are never uploaded: `_tables`, `_tmp`, `_full`, `_old`,
`_checks`, `_deploy` (rollback material) and `_retired` (local copies of pruned bucket paths).
Nothing keeps them out by name, because **uploads follow `manifest.json`, never a directory walk**:
the bucket is a copy of what the manifest names and nothing else. **`_tables/` holds every parquet
the build produces**: the browser reads none of them, only the `immutable/` files the manifest names
and two JSON assets. `_full/cis_sqtl_nominal` is the all-introns build parked for local use.

`validate` checks row counts against the raw files, eGene and sQTL counts against the preprint,
that six paper-named variants resolve to their rsIDs, that a FLNC query touches one row group,
and the rsID exact-match rate. For every chromosome it also reads the packs with its own decoder
(`steps_pack.validate`): block and page structure against SPEC.md section 7, gene details against
an independent rebuild from genes, exons, and splice phenotypes, credible-set record counts, and, for
sample genes, credible-set memberships, intron-run coverage against the raw sQTL rows, and a round
trip against the raw Zenodo eQTL rows within SPEC's section 9 limits. For the sQTL packs it checks
every block as kind 3, block, row, and credible-set totals against the raw files and credible_sets,
every intron's run against its raw rows and inside its gene's variants range, and the block pointers
in the gene details. It then runs a round trip against the raw Zenodo sQTL rows on sample introns:
every intron of `packs.sqtl_check_genes`, each chromosome's largest intron, introns with a variant in
two credible sets, and 5 random introns per chromosome. It writes reference files to
`data/derived/_tmp/pack_check/` for the browser decoder check, `cd ui && npm run pack-check`; the
sQTL ones go in `pack_check/sqtl/`, because the check reads every top-level `*_index.json` as an
eQTL chromosome. For GWAS it reads the TSV itself (its own DuckDB copy, not `pack_gwas`'s) and checks
every row of every GWAS pack against it in order, the index structure, and windows through the index:
50 seeded ones, a position repeated across a block boundary, each chromosome's first and last rows,
empty windows, and the FLNC, MYOZ1, and CAMK2D `[w_lo, w_hi]`. Its reference files, for chr1 and
chr22 windows, go in `pack_check/gwas/`. For the trans side (`steps_pack_trans.validate`, its own
frame reader) it checks every variants file's two sections against the variant table (counts,
records, af codes, that no run reaches `n_cis`), every trans frame of every file against the
source rows (order, rebuilt phenotype ids, positions, rsIDs exact; p, beta, `beta_se`, `r2`, and af
within SPEC section 9), the `search_index` pointers, and writes reference files to
`pack_check/trans/`.

For the variant page (`steps_pack_variant.validate`) it checks that `ord` is the row position, every
hits pack's structure and per-kind totals, about 10,000 sampled variant indices against their source
rows, the rsID index's block maths with 10,000 sampled lookups, `variant_index.qbx` against a fresh
walk of every file it indexes, and the six paper variants both ways. It then runs SPEC section 7's
**cis scan** in Python over the packs and compares it to the raw Zenodo nominal files at the same
position: rs10824026 (42 genes, 239 introns) and the lead of the first chrX eGene (24 and 75). Its
reference files go in `pack_check/variant/`.

Two whole-build checks close it out: every stored rsID text equals `'rs' || rs_number` of its
variant (the text columns, not just the variant table), and no `.parquet` exists anywhere in
`data/derived/` outside an `_*` folder.

`packcheck` checks, genome-wide, the facts the binary pack format in `SPEC.md` rests on: every
phenotype's tested variants are one contiguous run of the chromosome's cis variant list; `af`,
`ma_samples`, and `ma_count` never differ for a variant, within or across QTL types; one Student-t
degrees-of-freedom value per type rebuilds `slope_se` from `slope` and `pval_nominal`; and the window
start (`position - tss_distance`) is constant per phenotype. `roundtrip` (check 6 of SPEC Appendix
A.1) streams the raw nominal rows through the 16-bit codes and back, and checks `slope_se` and the
derived slope against SPEC's per-row error bound; it is cheap (every row genome-wide in minutes, about 1 GB). `measure`
settles the variant page size and codec. It reads the extracted Zenodo nominal files (falling back to the derived tables when a raw
file is missing, which for sQTL means significant introns only) and writes `report.md`,
`report.json`, `measure.md`, `roundtrip.md`, `roundtrip_<scope>.json`, `se_reference.json`, and
per-phenotype and per-variant parquet to
`data/derived/_checks/packcheck/` (`packs.checks_dir` in the config), which is never uploaded. It is
not a build step and leaves no `.done` marker: re-run it when the inputs change. `dof` and `measure`
update `packs:` in `config.yaml` when the data disagree with it.

## Pack tooling

`../PACKS.md` is the plain-language overview of the pack files (why, which file answers which
question, glossary); `../SPEC.md` is the byte layout. The code:

- `packfmt.py`: the reference codec. Every pack byte is encoded and decoded here; the build steps
  import its encoders and `validate` compares against its error bounds. No config, no file paths.
- `packtool.py`: a CLI and API over `packfmt` for anyone outside the pipeline. `header`, `blocks`,
  `block` (a gene's or an intron's rows as the SPEC section 8 table, its details, or its
  credible sets), `variants` (a run, a byte range, or a section as rows), `frames` and `trans` (a
  gene's trans frame as rows, through `search_index`), `gwas` (a window through the index),
  `gwas-index`, and `check` (decode a whole file); `pack-variants` (with `--trans-only`),
  `pack-results`, `pack-trans`, and `pack-gwas` build packs from parquet, Arrow, TSV, CSV, or
  JSON tables. Takes plain paths; the derived values' `dof` comes from `--dof` (or
  `--dof-eqtl`/`--dof-sqtl`), `--manifest`, or a `manifest.json` above the pack.
- `test_packfmt.py`, `test_packtool.py`: round-trip tests on synthetic data (plus FLNC and
  SYNPO2L, and chr21, when the build exists). `packtool_example.py` is the runnable walkthrough.
- `packcheck.py`: the genome-wide checks above; `steps_pack.validate` is the independent reader.

## Layout rules the browser depends on

The browser reads no parquet. These are the invariants its readers assume; breaking one silently
produces wrong pages, so `validate` checks every one of them.

**Publishing.**

1. **Every file the browser reads at byte offsets lives under `immutable/`** as
   `<stem>.<sha256[:16]>.<ext>`, flat, with the stem being the logical key and `/` written as `.`
   (`eqtl/chr1` -> `immutable/eqtl.chr1.<sha16>.qbe`). Change one byte and the name changes, so a
   key under `immutable/` is never overwritten and never has two builds on disk.
2. **`manifest.json` is the only file whose content changes under the same name.** The UI fetches
   it with `cache: 'no-cache'` and finds every content-addressed file through `manifest.packs`. No
   hashed name is written down in UI code.
3. **`search_index` records the SHA-256 of every pack it indexes,** in its Arrow schema metadata.
   The `manifest` step fails when they disagree, so manifest, index and packs always move together.
4. **Headers are set at upload** (`pipeline/upload.py`), never by a later server-side copy:
   `public, max-age=31536000, immutable` under `immutable/`, `no-cache` for everything else it
   writes. Packs are `application/octet-stream` with no `Content-Encoding`, because their byte
   offsets address the stored bytes.
5. **Old files leave the bucket only through `upload.py prune`,** only under the prefixes the
   `replaces:` map names, and only after `release` has made the new manifest live.
6. **`_tables/` is never uploaded.** Every parquet the build writes lives there, and only the
   builders and `validate` read it. What ships is the `immutable/` files the manifest names,
   `manifest.json`, and three small JSON files at the root: `coloc_loci.json` and
   `gwas_dcm_bins.json`, which the app fetches whole, plus `gwas_dcm.json`.

**Pack contents.**

- **A position is unique within a chromosome.** `(chr, position)` identifies one variant across
  all 9,215,026 of them, which is what lets a `chr:pos` lookup end after one page read.
- **A variants file is the cis section, then the trans-only section.** `vidx` is a variant's row
  number in that order, and `n_cis` (header byte 24) is the boundary. Every gene and intron run
  lies inside the cis section.
- **Hits rows name a gene by `search_index.ord`,** its row number in (chr, tss, gene_id) order,
  never by gene id. So `ord` must equal the row position, and the index must stay under 65,536 rows.
- **A hits frame covers 1,024 variant indices** (`packs.hits_frame_variants`), an rsID block holds
  **4,096 records** (`packs.rsid_block_records`), and a variants page holds **512** variants
  (`packs.variant_page_size`). All three are in the manifest and in the variant index's own header,
  and a reader checks that they agree.
- **An intron's window is anchored at its gene's TSS,** while a gene's window start is not the
  GENCODE TSS for 2,773 genes, so each block stores its own anchor.
- **Sort order decides what a row group's statistics cover** in `_tables/`. No browser reads those
  files any more, but the pack builders do, and the rule still costs build time: sorting by TSS
  instead of gene id made a gene's trans read the whole chromosome; sorting exons in GTF order made
  one gene touch sixteen row groups. Both were fixed by sorting on the filtered column.
- Concurrent DuckDB workers must not share a `temp_directory`; `connect()` takes one per worker.

## Conventions

- `gene_id` is an unversioned ENSG; `symbol` is the GENCODE name; `chr` is `chr1`..`chrX`.
- `A1` is the effect (minor) allele, `A2` the reference, as in the Zenodo release.
- sQTL `phenotype_id` is the leafcutter string `chr:start:end:clu_N_strand:ENSG.v`.
- `pval_perm < 0.05` is the significance flag; a BH `qval` on `pval_beta` is stored alongside.
