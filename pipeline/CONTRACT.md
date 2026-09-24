# The ingestion contract

`_tables/` is the seam in the qtlstore build. Above it sits one **adapter** per data source, which
knows that source's file names, column spellings and allele convention and nothing else. Below it
sit the **builders** (variant catalog, annotation, results, store), which know the qtlb format and nothing
about where the numbers came from. This file defines the seam: the five Parquet tables an adapter
must produce, the three it may produce (`trans`, `gwas`, `gwas_bins`), and the rules those tables obey.

Adding a second experiment should mean writing one adapter. If it means touching a builder, either
the contract is wrong or the adapter is cheating; fix whichever it is rather than special-casing
downstream.

Everything here is qtlb **v1**. The v0 tables the TOPCHeF build writes today are described at the
end under [What TOPCHeF's current tables look like](#what-topchefs-current-tables-look-like), with
the delta from this contract, because that delta is the work the TOPCHeF adapter has to do.

## Scope

An adapter produces these tables for **one experiment**: one cohort, one tissue, one set of
molecular phenotypes. An experiment may have several *phenotype types* (`ge` and `leafcutter` are
the two we have), and those share a `sites` table but get their own `nominal`, `permuted`,
`credible_sets` and `phenotypes` rows, distinguished by the `phenotype_type` column.

An adapter does **not**: fetch the reference, assign variant ordinals, compute digests, decide page
sizes, read or write pack bytes, or touch the gene annotation. Those are builder jobs.

## The five tables

Paths are under `data/derived/_tables/<experiment_id>/`. Every table is Parquet. Column order is
not significant; names and types are.

### `sites.parquet`

Every distinct variant this experiment tested, once. The variant catalog builder reads only this.

| column | type | meaning |
|---|---|---|
| `chr` | `string` | sequence name, exactly as it appears in the refget collection (`chr1`, not `1`) |
| `pos` | `int32` | 1-based, on the sequence named by `chr` |
| `ref` | `string` | reference allele: the base(s) the anchored sequence actually reads at `pos`, upper case |
| `alt` | `string` | alternate allele, upper case |
| `af` | `float` | **ALT** allele frequency in this cohort, in [0, 1]; `NaN` when the source does not report it |
| `ma_samples` | `int32` | samples carrying the minor allele; `-1` when the source does not report it |
| `in_cis` | `bool` | true if the site appears in any cis result for this experiment |

Rules:

- `(chr, pos, ref, alt)` is unique and is the site's identity. Row order does not matter: the
  variant catalog builder sorts, and the identity digest runs over the set of sites in a canonical order
  (SPEC.md section 5). Sorting by `chr` then `pos` is still good practice for the parquet reader.
- `ref` is the anchored reference. It is not "the first allele in the source file" and not "the
  major allele". The adapter must verify it against the refgetstore (reuse `steps_refget.classify`)
  and **drop** sites where neither source allele matches, counting them in the ingestion report.
  A site whose `ref` is not the reference cannot be joined across studies or given a VRS id, which
  is the whole point of the table.

  **Decided: variants that report no alleles are left out.** *Decided 2026-09-24; replaces the
  dbSNP rsID recovery of 2026-09-22.* TOPCHeF has 32,765 positions that only the trans eQTL file
  reports, as a bare `chr:pos` with no A1/A2 (v0 SPEC section 4, flags bit 2). The release does
  know them: the authors' pipeline kept exactly one tested allele per position, but it joined A1/A2
  onto the cis and trans sQTL files and not onto the trans eQTL files. The rule:

  - A variant with no **source-provided** alleles gets **no `sites` row**, and no row that names it
    (its trans eQTL pairs) is ingested. Alleles are **not inferred** -- not from dbSNP by rsID, not
    from population frequency, not from any other study. A guessed allele that is wrong looks
    exactly like a right one once it is anchored, and it would carry a real VRS id.
  - The exclusion is counted, not silent: `ingestion.json` records `trans_eqtl_excluded` with the
    number of variants, the number of trans eQTL rows left out (and the total), the genes those rows
    belong to, and the reason. The orientation counts list the same variants under the class
    `no_source_alleles`.
  - Trans sQTL rows carry A1/A2 in the source and are unaffected.
  - This holds **until the authors supply the alleles** (asked for: A1/A2 on the trans eQTL files,
    or their `.bim` files). When they do, the variants enter as ordinary sites from the source's
    own alleles, through the same refcheck as every other variant.
- `af` is the ALT frequency, always, even when the source reports the minor-allele frequency. If the
  source's effect allele is REF, the adapter swaps (see [Orientation](#orientation)).
- Indels follow VCF convention: both alleles carry the shared leading base, so `AT>A` is a one-base
  deletion recorded at the position of the `A`.
- The seven columns above are required. An adapter may carry **declared variant catalog attributes**
  alongside them: extra per-site columns the variant catalog object lists and the identity digest ignores.
  TOPCHeF carries `ma_count`, `rsid`, `rs_number` and `match`, which the variant page shows and
  which the rsID index is built from. A builder reads an attribute only because the variant catalog JSON
  named it, never because it happened to be in the file.
- Declared attributes with a fixed meaning, so two adapters that carry them agree:

  | attribute | type | meaning |
  |---|---|---|
  | `ma_count` | `int32` | minor-allele count in this cohort; `-1` when unknown. Optional. |
  | `rsid` | `string` | dbSNP id (`rs123`), or null. Must be exactly `'rs' || rs_number` (null where `rs_number` is `-1`); the builder fails otherwise. Not stored: a reader rebuilds it from `rs_number`, so the variant catalog lists `rs_number`, not `rsid` |
  | `rs_number` | `int64` | numeric part of `rsid`, or `-1` |
  | `match` | `string` | how the rsID was matched: `exact`, `position` or `none` (null = `none`); no other value |

### `nominal.parquet`

One row per (phenotype, variant) association tested. This is the big one: tens of millions of rows
for a single experiment, so it is partitioned by `chr` (`nominal/chr=chr1/data.parquet`).

| column | type | meaning |
|---|---|---|
| `phenotype_type` | `string` | `ge`, `leafcutter`, ... |
| `phenotype_id` | `string` | the source's phenotype identifier |
| `gene_id` | `string`, nullable | unversioned gene id (`ENSG00000128591`) the phenotype belongs to; null when the source publishes none. For a multi-gene phenotype, the primary gene |
| `chr` | `string` | variant sequence name |
| `pos` | `int32` | variant position |
| `ref` | `string` | as in `sites` |
| `alt` | `string` | as in `sites` |
| `beta` | `float` | effect of the **ALT** allele |
| `se` | `float` | standard error of `beta` |
| `pvalue` | `double` | nominal p-value as the source reports it |

Rules:

- Every `(chr, pos, ref, alt)` must exist in `sites`. The builder fails loudly on an orphan rather
  than dropping it.
- `beta` is ALT-relative. When the adapter swaps alleles it negates `beta` and leaves `se` alone.
- `pvalue` is kept as given. The builder derives `-log10(p)` from it and, separately, checks whether
  `beta/se` reproduces it under some degrees of freedom (see `dof` below). Those must agree, but the
  stored p is the source's.

### `permuted.parquet`

One row per **phenotype group**: the permutation pass, and the lead variant. A group is the unit the
source permuted (see `phenotype_object_id` under [`phenotypes`](#phenotypesparquet)).

| column | type | meaning |
|---|---|---|
| `phenotype_type` | `string` | |
| `phenotype_object_id` | `string` | the group that was permuted |
| `phenotype_id` | `string` | the group's **lead** phenotype (e.g. the lead intron of a leafcutter cluster) |
| `gene_id` | `string`, nullable | primary gene of the group; null when the source publishes none |
| `n_variants` | `int32` | variants tested in this group's window |
| `p_perm` | `double` | permutation p-value (the significance rule reads this by default) |
| `p_beta` | `double` | beta-approximated p-value (read instead when the rule's `column` is `p_beta`) |
| `lead_chr` | `string` | lead variant |
| `lead_pos` | `int32` | |
| `lead_ref` | `string` | |
| `lead_alt` | `string` | |

Rules:

- Exactly one row per `(phenotype_type, phenotype_object_id)`. `phenotype_id` is the lead
  phenotype, which must exist in `phenotypes` with the same `phenotype_object_id`.
- The lead variant must exist in `sites`. Where the lead phenotype has `nominal` rows, the lead
  variant must be one of them.
- `n_variants` is **recounted from `nominal`** wherever nominal is complete for the group, and
  otherwise **copied from the source**; the ingestion report says, per `phenotype_type`, which of
  the two it did. For GTEx `ge` the two agree once the source's duplicate rows (one per dbSNP id
  at a site) are removed: 153,275,168 nominal rows vs 153,275,183 summed `n_variants`, 15 genes
  differing.
- The recount counts **tested** rows: nominal rows whose `pvalue` is null or NaN are kept in
  `nominal` (the source publishes them) but not counted. *Measured 2026-09-24.* TOPCHeF's nominal
  files carry 29 variants genome-wide (e.g. `chr18:77490468`, `chr11:1438027`) with `af` 0.5,
  `ma_samples` = `ma_count` = 516, and null slope, SE and p: every sample is heterozygous, so the
  genotype has no variance and nothing can be tested. tensorQTL's nominal pass writes such rows;
  its permutation pass drops monomorphic variants before counting `num_var`. Counting every row
  put the recount one above `num_var` for exactly the 525 genes and 2,124 introns whose window
  holds one of them (each holds exactly one). The source is right; counting tested rows makes the
  recount equal `num_var` for every TOPCHeF phenotype. The ingestion report gives
  `nominal_rows_without_pvalue` per `phenotype_type`.
- The eQTL Catalogue permuted file's lead p-value and beta differ slightly from the nominal row
  for the same variant. `permuted` carries the permuted file's values; do not mix in nominal ones.
- The significance rule itself (`p_perm < 0.05` for TOPCHeF) is **not** here: it is a field in the
  experiment JSON, so that two experiments can use different rules and a reader can see which.
  The adapter reports the numbers; the experiment declares the rule: `column` (`p_perm` or
  `p_beta`), `op` (`<` or `<=`) and `threshold`. The results builder tests exactly that column.

**Decided: `permuted` is keyed by the group, not the phenotype.** *Decided 2026-09-22.* The eQTL
Catalogue permutes on `molecular_trait_object_id`, not `molecular_trait_id`. For `ge` the two are
the same. For `leafcutter` the object is the **cluster** and the trait is the **intron**: GTEx
leafcutter is 46,249 clusters over 176,603 introns, and keying `permuted` or `phenotypes` on the
cluster would silently drop 130,354 introns. So `phenotype_object_id` is carried in `phenotypes`,
`permuted` and `credible_sets`; `permuted` has one row per group; `phenotypes` keeps one row per
phenotype (every intron), each naming its group. Fanning the cluster's row out to its introns was
rejected because it writes the cluster's lead variant and p-values onto every intron, a number that
is wrong without looking wrong. For `ge`, and for TOPCHeF everywhere (tensorQTL permutes each
phenotype separately), `phenotype_object_id = phenotype_id` and nothing about the tables changes.

### `credible_sets.parquet`

Fine-mapping output, one row per (phenotype group, credible set, variant).

| column | type | meaning |
|---|---|---|
| `phenotype_type` | `string` | |
| `phenotype_object_id` | `string` | the group that was fine-mapped |
| `phenotype_id` | `string` | the phenotype the source attaches the set to (the group's lead where the source fine-maps groups) |
| `cs_id` | `int16` | credible-set number within the group, 1-based, parsed from the source (GTEx `ENSG..._L3` -> 3) |
| `chr` | `string` | variant |
| `pos` | `int32` | |
| `ref` | `string` | |
| `alt` | `string` | |
| `pip` | `float` | posterior inclusion probability, in [0, 1] |
| `z` | `float` | z-score; `nan` when the source does not report it |
| `cs_size` | `int32` | variants in this credible set |
| `cs_min_r2` | `float` | minimum pairwise r2 within the set; `nan` when not reported |

Rules:

- A variant may appear in more than one credible set of the same phenotype. That is real, and the
  hits pack represents it; do not deduplicate.
- Every variant must exist in `sites`.
- An experiment with no fine-mapping writes the table with zero rows, not no table.

### `phenotypes.parquet`

One row per phenotype: what it is, independent of any result. Every phenotype the source
describes gets a row, including those (every non-lead intron) that never appear in `permuted`, and
those the source reports only in trans (TOPCHeF: 17 genes, most on chrM, and 36 introns with trans
rows and no cis result), which get `has_nominal` false and no group.

| column | type | meaning |
|---|---|---|
| `phenotype_type` | `string` | |
| `phenotype_id` | `string` | |
| `phenotype_object_id` | `string` | the group this phenotype was permuted/fine-mapped in: eQTL Catalogue's `molecular_trait_object_id`; for leafcutter, the cluster. Equals `phenotype_id` for `ge` and for all of TOPCHeF |
| `gene_id` | `string`, nullable | unversioned ENSG; null when the source publishes none. The sumstats name no gene for most GTEx leafcutter clusters, but the eQTL Catalogue leafcutter metadata (Zenodo 7850746) gives every intron one, so the adapter fills it from there. For a multi-gene phenotype, the primary gene |
| `has_nominal` | `bool` | whether any `nominal` rows exist for this phenotype |
| `extra` | `string` | JSON object, source-specific fields |

Rules:

- Exactly one row per `(phenotype_type, phenotype_id)`.
- A phenotype that belongs to several genes keeps the primary in `gene_id` and the full list in
  `extra` as `{"gene_ids": [...]}`.
- `has_nominal` exists because nominal coverage can be sparse: the eQTL Catalogue keeps nominal
  rows for only 2,243 of 176,603 GTEx leafcutter introns. A reader uses it to say "no per-variant
  data published" instead of showing an empty result as if nothing were associated. The ingestion
  report records the counts per `phenotype_type` (phenotypes, phenotypes with nominal rows).

- `extra` is where phenotype-type-specific structure goes, and it is the reason a third phenotype
  type does not need a new file kind. For `leafcutter`: `{"intron_start": ..., "intron_end": ...,
  "cluster_id": "clu_1234", "strand": "+"}`. For `ge` it is usually `{}`.
- **No annotation columns.** No `symbol`, no `tss`, no `biotype`, no gene `start`/`end`/`strand`.
  Those belong to the annotation object, which is keyed by `gene_id` and shared across experiments.
  This is the single most common way to break the contract, because the v0 tables are full of them.

## The optional tables

### `trans.parquet`

Trans associations: one row per (phenotype, variant) pair the source tested in trans and reports.
Absent when the experiment has no trans results (the eQTL Catalogue GTEx datasets have none).

| column | type | meaning |
|---|---|---|
| `phenotype_type` | `string` | |
| `phenotype_id` | `string` | as in `phenotypes` |
| `gene_id` | `string`, nullable | as in `phenotypes` |
| `chr` | `string` | variant sequence name |
| `pos` | `int32` | variant position |
| `ref` | `string` | as in `sites` |
| `alt` | `string` | as in `sites` |
| `beta` | `float` | effect of the **ALT** allele; finite |
| `se` | `float` | standard error of `beta` as the source reports it; not stored (a reader derives it from p and dof) |
| `pvalue` | `double` | in [0, 1], not null |

Rules:

- Every `(chr, pos, ref, alt)` exists in `sites` and every `(phenotype_type, phenotype_id)` in
  `phenotypes`; the builder fails otherwise. A trans-only variant is in `sites` with `in_cis` false.
- `(phenotype_type, phenotype_id, chr, pos, ref, alt)` is unique.
- **Only variants with source alleles.** A trans row whose variant the source names without alleles
  is left out, as in `sites` (TOPCHeF: 37,951 of 2,680,117 trans eQTL rows, counted in
  `trans_eqtl_excluded`). A trans row that names the variant only as `chr:pos` may take the alleles
  another file of the same release gives that position, but only when the release holds one variant
  per position (TOPCHeF: the authors' plink2 `--rm-dup force-first` kept one tested allele per
  position; the adapter checks the orientation table has one row per position and fails otherwise).
  This is the release's own allele, not an inferred one.
- `beta` is ALT-relative: negated where the adapter swapped the alleles, as in `nominal`.
- Sort by `(phenotype_type, phenotype_id)` for the reader; the builder sorts anyway.
- The ingestion report gives `rows.trans` and a `trans` section with rows per type and the source's
  row counts.

### `gwas.parquet` and `gwas.json`

A GWAS shown next to this experiment (the design plan attaches a GWAS to an experiment; for TOPCHeF
it is the DCM GWAS, written by `pipeline/adapters/dcm_gwas.py` into `_tables/topchef/`). Absent when
the experiment has none.

| column | type | meaning |
|---|---|---|
| `chr` | `string` | sequence name, in the variant catalog's chromosome list |
| `pos` | `int32` | 1-based |
| `ref` | `string` | the reference base(s) at `pos` |
| `alt` | `string` | the other allele |
| `beta` | `double` | effect of **ALT**, at most 4 decimals |
| `se` | `double` | at most 4 decimals |
| `af` | `double` | **ALT** frequency, at most 4 decimals |
| `pvalue` | `double` | in (0, 1], at most 4 significant digits |
| `n` | `int64` | sample size (at most 255 distinct values) |
| `rs_number` | `int64` | dbSNP number, 0 or -1 for none |

Rules:

- Orientation as everywhere: the adapter reads the reference at every row and applies
  `qtlstore.orient_to_ref` with the source's effect allele; rows whose alleles the reference does not
  read are dropped and counted in `gwas.json`. The GWAS rows do not have to be sites of the
  variant catalog: the GWAS object stores its own positions and alleles.
- Values are lossless at the source's printed precision; the builder refuses anything finer.
- **A site may appear more than once** when the source reports it more than once with different
  numbers (the DCM GWAS lists 796,531 indels once per allele order, with different N and statistics:
  two measurements, both kept). Identical rows are not allowed; the adapter drops the second copy and
  counts it. `gwas.json` reports `rows_sharing_a_site` and `identical_rows_dropped`.
- `gwas.json`: `id`, `title`, `source` (file, cases, controls, source row counts), `orientation`
  (`as_is`, `swapped`, `dropped`, counts by reference class), `rows`.

### `gwas_bins.parquet`

The landing track's summary of the GWAS (SPEC.md section 10, "Bin summary"), optional next to
`gwas.parquet`. One row per window of `bin_bp` (5 Mb, `gwas_bin_bp` in the config) that has source rows,
from **the source rows as published**, not the oriented `gwas` rows: rows the reference does not read
count, and the lead's allele and beta are the source's effect allele and its beta. The DCM GWAS adapter
writes it with v0's `gwas_bins` query unchanged (`dcm_gwas.BINS_SQL`; `--bins-only` writes just this
table), so it equals v0's `gwas_dcm_bins.json`.

| column | type | meaning |
|---|---|---|
| `chr` | `string` | sequence name, in the variant catalog's chromosome list |
| `bin_start`, `bin_end` | integer | `[k * bin_bp, (k + 1) * bin_bp)`; a position `pos` is in bin `pos // bin_bp` |
| `min_p` | `double` | smallest p in the bin, unrounded (the builder rounds to 3 significant digits) |
| `lead_position`, `lead_rsid`, `lead_beta`, `lead_ea` | | the row holding `min_p`: position, rsID text (null when none), beta (unrounded; the builder rounds to 3 decimals), effect allele |
| `n_gws`, `n_variants` | integer | rows with p < 5e-8; rows in the bin |

## The ingestion report

Each adapter writes `ingestion.json` next to the tables. Beyond what individual rules above
require, it records, per `phenotype_type`: the number of phenotypes, groups, and phenotypes with
`nominal` rows (`has_nominal`), and whether `permuted.n_variants` was recounted from nominal or
copied from the source. Every site dropped or held out of `sites` is counted there by reason.

## Orientation

The one convention everything else rests on: **`ref` is the anchored reference base, `alt` is the
other allele, `af` and `beta` are ALT-relative.**

Sources do not agree on this. The eQTL Catalogue already gives `ref`/`alt` with the effect allele
being ALT, so its adapter passes them through. TOPCHeF gives `A1`/`A2` with no stated orientation.
Two separate questions had to be answered for it, and both are now measured:

1. **Which allele is REF.** The refcheck reads GRCh38 at every variant: **A2 is the reference
   allele**, at all 8,419,594 cis SNPs without exception, with `neither` equal to 0 across all
   8,872,723 cis variants and a match fraction of exactly 1.000000. The 97,433 apparent A1 matches
   are all indels, where prefix matching cannot break the tie and `classify` takes the longer
   allele. EVIDENCE.md A.10, in the analysis repo at `qtlb-format/docs/`.
2. **Which allele the values describe.** `af` and `beta` both describe **A1**. Measured against the
   Jurgens 2024 DCM GWAS: over the 6,793,566 biallelic SNPs this release shares with it, our `af`
   correlates **+0.9932** with that study's `EAFREQ` oriented to A1 (mean absolute difference
   0.022), and -0.9932 against the mirror `1 - eaf` (mean absolute difference 0.668). The 35,444
   rows where the GWAS effect allele is our A2, so the frequency had to be mirrored to compare,
   agree at +0.9975 -- the join is not answering its own question. No palindrome could confound it:
   that GWAS carries zero strand-ambiguous SNPs. Cis `af` averages 0.2420 with median 0.1318 over
   SNPs, an ALT-side distribution. `beta` follows from `af` by construction: `steps_nominal` copies
   `n.af` and `n.slope` from one tensorQTL row unchanged, and tensorQTL codes both against the same
   dosage vector. EVIDENCE.md A.11.

Put together: A2 is REF, A1 is ALT, `af` is already the ALT frequency and `beta` is already the ALT
effect. **So the TOPCHeF reorientation is a relabelling, not a sign flip.** The adapter maps
`A2 -> ref` and `A1 -> alt`, keeps `af` and keeps `beta`, both untouched.

> **Corrected: what to do at the 97,433 indels the refcheck table calls `a1`.** *Measured by the
> TOPCHeF adapter, 2026-09-22. This paragraph used to say "take the call from
> `_tables/refcheck/<chr>.parquet` rather than re-deriving it", which cannot be done at the same
> time as the acceptance gate's "zero swapped sites": taking that call literally makes ref = A1 at
> 97,433 cis sites, which is 97,433 swaps and 97,433 negated slopes.*
>
> The refcheck class is a **reading** at SNPs and a **tie-break** at indels. An indel's two alleles
> share a leading base, so prefix matching cannot separate them, and `steps_refget.classify` takes
> the longer allele (EVIDENCE.md A.10 says so outright). The 97,433 are where that tie-break fired.
>
> dbSNP does not settle which allele GRCh38 reads there. Joining those 97,433 to the dbSNP records
> the rsID step already pulled at exactly those positions: REF is A1 at 41,560, A2 at 38,034, both
> orientations are present at 12,894, and 4,945 have no matching record. That split is what a
> position carrying both a deletion and an insertion looks like, not a verdict.
>
> `af` settles it. `af` is A1's frequency, and on the 41,560 where dbSNP puts REF on A1 it averages
> **0.1504** with median **0.0833**, above 0.5 at only **5.3%** of sites -- statistically the same
> as the **11.7%** of the undisputed `a2` indels, and nothing like the near-90% a genuine reference
> allele would show. A1 is the minor allele at those sites, and a minor allele is not what the
> reference reads. So A1 is ALT there too.
>
> The rule, and what the adapter implements: **take A2 as `ref` at a disputed indel wherever the
> reference actually reads A2 at that position**, asked of the refgetstore variant by variant, and
> keep the refcheck call where it does not. Nothing is assumed: both `G>GC` and `GC>G` are correctly
> anchored where the reference reads `GC`, so the only question that needs asking is whether A2
> reads, and it is asked rather than argued. `pipeline/adapters/topchef.py:INDEL_REFERENCE` carries
> this reasoning next to the code that acts on it.

Adapters do not implement this by hand. `pipeline/qtlstore.orient_to_ref(ref_base, effect_allele,
other_allele, beta, af)` is the one implementation, and its parameters name the convention rather
than any source's column letters: `beta` and `af` must describe `effect_allele`. For TOPCHeF the
plain call `orient_to_ref(ref, A1, A2, slope, af)` is the correct one and is a no-op on 98.9% of
sites. `pipeline/test_qtlstore.py::orientation_topchef_is_a_relabelling_not_a_sign_flip` is the
regression guard.

Each adapter records how it decided, as `allele_orientation_source` in the experiment JSON
(`topchef_refcheck_a2_is_ref`, `eqtl_catalogue_ref_alt`). A reader should be able to tell a
convention that was verified from one that was assumed.

### Acceptance gate for a TOPCHeF re-ingestion

This is the authoritative wording; it replaces the phrase "bit-identical up to the sign flip on
swapped sites" in the qtlstore design plan, step 2. That phrasing pre-authorized the exact failure
it should have caught. SPEC.md section 11 stores a slope's sign in one bit of the SE u16, so a
genome-wide inversion re-encodes to perfectly valid bytes, survives every precision and round-trip
check, and is indistinguishable from a legitimate result under a gate that expects sign flips.

A rebuild of TOPCHeF through the contract tables is accepted only when, against the v0 packs:

- Every decoded `-log10 p` and SE array is **bit-identical**.
- Every decoded slope is **bit-identical, sign included**. Zero sign flips is the expected count,
  not an allowance: A2 is REF and `beta` is A1's, so no cis site swaps.
- Every `af` is **bit-identical**. No site takes `1 - af`.
- `ref` equals the old `A2` and `alt` equals the old `A1` at every SNP. At every indel the pair is
  one of those two orderings and the choice follows the corrected rule above: A2 wherever the
  refgetstore confirms the reference reads A2 there, the `_tables/refcheck/<chr>.parquet` call
  otherwise. The ingestion report gives the split, so a reader can see how many were decided which
  way rather than taking the count on trust.
- The count of swapped sites and the count of dropped sites are both reported in the ingestion
  report, and both are **0** for cis variants.

If any slope changes sign, the rebuild is wrong. Do not accept it as orientation.

## Degrees of freedom

Not a contract table. The results builder needs a `dof` per results set so a reader can rebuild
`beta` from `se` and p. Where the source publishes it (TOPCHeF: 435 eQTL, 480 sQTL) the adapter
passes it through. Where it does not (eQTL Catalogue), `pipeline/dof.py` fits it from `nominal`
and records the fit residual. A fit that does not converge is not fatal: the experiment stores
`dof: null` and the reader uses the stored `-log10(p)` directly, which is what it does today
anyway.

## What TOPCHeF's current tables look like

The v0 build writes these, and they are *not* the contract. The delta is the TOPCHeF adapter's job:

| v0 table | rows | what has to change |
|---|---|---|
| `variants.parquet` | 9,215,026 | `A1`/`A2` -> `ref`/`alt` by the refcheck call; `position` -> `pos`; keep `rsid`/`rs_number` as variant catalog attributes, not contract columns |
| `genes.parquet` | 60,624 | split three ways: annotation columns to the annotation object, result columns to `permuted`, `bin` dropped (it named a build partition) |
| `splice_phenotypes.parquet` | 80,750 | same split; `intron_start`, `intron_end`, `cluster_id`, `strand` go into `phenotypes.extra` |
| `credible_sets.parquet` | 461,784 | drop `symbol`, `tss`, `gene_id`, `rsid`, `af` (all derivable); `qtl_type` -> `phenotype_type`; add `z`, `cs_size`, `cs_min_r2` |
| `trans/`, `trans_by_variant/` | 15,862,525 pairs | `trans.parquet` (15,824,574 rows with source alleles); `gene_chr`/`gene_tss`/`symbol` come out |
| nominal build intermediates | | become `nominal/chr=*/` with `ref`/`alt` |

Two things to notice. First, almost every deletion above is an annotation column that v0 copied into
three different tables; the annotation object exists so that a second experiment on GENCODE v39
cannot silently disagree with this one on where a gene starts. Second, `qtl_type` becoming
`phenotype_type` is not a rename: v0 has exactly two values baked into the file kinds, and the
contract makes it an open set carried in data.
