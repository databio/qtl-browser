---
date: 2026-09-08
status: complete
model: Claude Fable 5.1
description: Give trans-only variants rsIDs and a variant page by adding trans file positions to the dbSNP lookup
---

# Trans-only variants: rsIDs and variant pages

## Why

The variant table (`variants_by_position/`, `variants_by_rsid.parquet`) is built from the cis
files only, and rsIDs are assigned by a dbSNP lookup over those positions. The trans step
left-joins that table by position. TensorQTL's trans scan is genome-wide, so any trans variant
outside every cis window falls through the join with a null rsID, and its variant page says it
was never tested.

Measured on the current derived tables:

| | Trans rows | Rows without rsID | Distinct positions without rsID |
|---|---|---|---|
| eQTL | 2,680,117 | 133,497 (5.0%) | 107,071 |
| sQTL | 13,182,408 | 684,320 (5.2%) | 310,002 |
| Both | 15,862,525 | 817,817 (5.2%) | 342,817 |

Of the 342,817 positions, 342,284 are more than 1 Mb from the TSS of every gene tested in cis;
the remaining 533 sit at window edges where tensorQTL measured the window from gene start or
end. Allele frequency is not a factor (same AF distribution as rows with an rsID). Example:
chr1:82572797 lies in a 550 kb gap between two cis windows; the eight annotated genes within
1 Mb are untested lncRNAs and pseudogenes.

Allele availability differs by file. `trans_sQTL` carries `chr`, `position`, `A1`, `A2`;
`trans_eQTL` carries only `variant_id` as `chr:pos`. Of the 342,817 positions, 310,002 appear
in trans sQTL rows (alleles known) and 32,815 appear only in trans eQTL rows (alleles unknown).

A trans-only variant page can show the rsID, position, AF, the DCM GWAS panel where the GWAS
covers the position, and the full trans association table, which is the main content. It
cannot show cis rows, credible-set membership, or lead status, and should say why.

## Change

### Pipeline

**`steps_variants.collect`** unions the trans files into the distinct-variant scan and adds an
`in_cis` boolean:

- `cis`: distinct `(chr, position, A1, A2)` from the six cis files, as now. `in_cis = true`.
- `trans_sqtl`: distinct `(chr, position, A1, A2)` from `trans_sQTL`, anti-joined on
  `(chr, position)` against the cis positions. `in_cis = false`.
- `trans_eqtl`: distinct `(chr, position)` from `trans_eQTL` (`variant_id` split), anti-joined
  against both the cis positions and the trans sQTL positions. `A1`, `A2` null. `in_cis = false`.

A position that appears in cis gets only its cis rows. A position that appears in trans sQTL
with alleles never also gets a null-allele row from trans eQTL.

**`steps_variants.rsid`**:

- Rebuilds `dbsnp_targets.tsv` and `dbsnp_matched.tsv.gz` when `variants_raw.parquet` is newer
  than they are, instead of trusting their existence. Without this the step silently reuses the
  cis-only cache and the new positions never reach bcftools.
- The exact-allele join is unchanged. Null alleles never satisfy it, so eQTL-only trans
  variants fall through to the position match, as `match = 'position'`.
- `in_cis` is carried into both output tables. The match-rate log line is split by `in_cis`.
- Row groups, sort order, and statistics columns are unchanged.

**`steps_tables.trans`** needs no code change. Its per-position rsID lookup reads the variant
table, so a forced re-run fills the `rsid` column for the trans-only positions in both
`trans_pairs/` and `trans_by_variant/`.

**`steps_finish`**:

- `manifest`: `counts.rsid_match` becomes the match breakdown for `in_cis` rows only, and a new
  `counts.variants_trans_only` records how many rows have `in_cis = false`.
- `validate` check 6 (exact-match rate at or above 90%) runs on `in_cis` rows only, so adding
  position-only rows cannot move it.
- New check: the share of `trans_pairs/` rows with a null `rsid` is below 0.1%. The remaining
  nulls are positions dbSNP b157 has no record for.
- New check: the number of `in_cis` rows equals the previous variant-table row count
  (8,872,723), so the cis side is provably untouched.

**`pipeline/README.md`**: the `variants_collect` row of the step table says it scans the cis
and trans files and records `in_cis`.

### UI

**`lib/queries.ts`**: `VariantRow` gains `in_cis: boolean`. `A1` and `A2` become nullable.

**`routes/Variant.tsx`**:

- Not-found message: "… is not among the variants tested in TOPCHeF (MAF ≥ 0.01; cis windows
  within 1 Mb of a tested gene, or genome-wide trans)." The dbSNP fallback link stays.
- Alleles row: when `A1` is null, shows "not reported (trans eQTL file has no allele columns)".
- rsID match row: when `A1` is null and `match = 'position'`, says "by position (alleles not
  reported)" instead of "alleles differ from dbSNP record".
- gnomAD and Open Targets links need alleles; hidden when `A1` is null. UCSC, dbSNP, and
  Ensembl stay.
- When `in_cis` is false, the "Lead variant for", "Credible-set membership", and "All cis
  associations" sections each show one line: "Outside every cis window (more than 1 Mb from any
  gene tested in cis)". The scan button is not rendered, since the scan would read every
  overlapping window and find none. The trans section is unchanged and is the page's content.

**`components/TransTable.tsx`**: no change. Its dimmed chr:pos fallback now appears only for
positions dbSNP lacks.

### Run

```
rm data/derived/.done/{variants_collect,variants_rsid,trans,manifest}
uv run python -m pipeline build --step variants_collect --step variants_rsid --step trans --step manifest \
  2>&1 | tee data/derived/build_12_trans_variants.log
uv run python -m pipeline validate
```

Expected time: bcftools re-streams dbSNP, about five minutes last time; the trans step took
under a minute; the rest is seconds. Expected size changes are in the next section.

No R2 upload until local testing and Sam's explicit go-ahead.

## Decisions & ownership

| # | Decision | Owner | Note |
|---|---|---|---|
| D1 | Include trans-only variants in the variant table and give them a variant page | user-owned, surfaced | Sam confirmed 2026-09-08 after seeing the counts |
| D2 | Add an `in_cis` column rather than inferring cis status in the UI | AI-owned, defended | inferring needs the cis window scan, which is the slow on-request query; a boolean is free |
| D3 | eQTL-only trans variants (32,815 positions) get position-only rsID matches and null alleles | silently implied, now surfaced | the trans eQTL file has no allele columns; nothing else in the release carries them. The UI must handle null alleles (D6) |
| D4 | A position never gets both an allele-bearing and a null-allele row | AI-owned, defended | avoids duplicate variant pages for one position |
| D5 | Invalidate the dbSNP cache files by mtime instead of deleting them by hand | AI-owned, defended | the existing `if not exists` check would silently skip the new positions; a manual delete is easy to forget on the next rebuild |
| D6 | Hide gnomAD and Open Targets links when alleles are null; UCSC, dbSNP, Ensembl stay | AI-owned, defended | both URLs embed ref and alt; a null would build a broken link |
| D7 | `permutation_tables`, `credible_sets`, `nominal`, `gene_detail` are not re-run | AI-owned, defended | they join the variant table only on cis positions, whose rows do not change; the validate check on the `in_cis` row count proves it |
| D8 | Cis sections on a trans-only variant page show a one-line explanation and no scan button | AI-owned, default | alternative is to drop the three sections entirely; keeping them with a reason reads clearer than a page that is silently shorter |
| D9 | `match` keeps its three values; null alleles are distinguished in the UI text | AI-owned, default | a fourth value like `position_no_alleles` would be cleaner but touches the validate check and manifest; can be added later |
| D10 | Not-found message rewritten to name both cis and trans coverage | AI-owned, defended | the current text becomes false once trans-tested variants resolve |
| D11 | No R2 upload without Sam's go-ahead | user-owned, surfaced | free tier cap; same rule as the prior plans |
| D12 | Search is unchanged | AI-owned, defended | rsID and chr:pos searches already route to the variant page, which now resolves these positions |

## What this changes elsewhere

- **Variant tables grow.** 342,817 rows added to 8,872,723 (3.9%), two new columns' worth of
  bytes for `in_cis`. Roughly +6 MB on `variants_by_position/` (145 MB) and +4 MB on
  `variants_by_rsid.parquet` (109 MB). Row-group count grows in step, so the per-lookup read
  size is unchanged. The manifest's `columns` list for both tables changes.
- **Trans tables are rewritten in full.** `trans_pairs/` (438 MB) and `trans_by_variant/`
  (540 MB) get the same rows with the `rsid` column filled. Local rewrite, and about 980 MB of
  R2 re-upload when Sam approves. Bucket total is roughly unchanged since the rows are the same.
- **R2 upload delta**: about 1.25 GB across the four table dirs. Class A operations: 24 + 23
  + 23 + 1 files. Free-tier storage cap is 10 GB; the bucket size should move by under 20 MB.
- **`.done` markers**: four removed and re-created. The other nine steps stay done.
- **`_tmp` caches**: `dbsnp_targets.tsv` (197 MB) and `dbsnp_matched.tsv.gz` (93 MB) are
  regenerated and grow slightly.
- **`VariantRow` type change** ripples to every consumer of `A1`/`A2`: the Variant page header
  and KV tables, and the gnomAD and Open Targets URL builders.
- **Validate check semantics change**: the exact-match rate is now a cis-only statistic. The
  manifest's `rsid_match` likewise.
- **Position-only matches appear for a new reason.** Until now `match = 'position'` meant the
  alleles disagreed with dbSNP. For null-allele rows it means the alleles were never reported.
  The UI text distinguishes the two; the parquet column does not (D9).
- **The memory note on the variant table** (8.87M cis variants) stays true for cis; the table
  itself is no longer cis-only. Update the note after the build.
- **Gene page trans tables** show rsIDs for about 5% more rows and the chr:pos fallback almost
  never.

## Implementation log (2026-09-08)

Decision review: Sam chose the recommended option on all three surfaced items (include the
eQTL-only positions with null alleles; keep the three cis sections with a one-line reason;
rebuild only the four steps).

Code: `steps_variants.collect` and `.rsid`, `steps_finish.manifest` and `.validate`,
`pipeline/README.md`, `ui/src/lib/queries.ts` (`VariantRow`), `ui/src/routes/Variant.tsx`,
`ui/src/routes/About.tsx` (counts split into cis and trans-only). All as specified above.

Build: `.done` markers for the four steps removed, then `build --step variants_collect --step
variants_rsid --step trans --step manifest`, then `validate`. Log in
`data/derived/build_12_trans_variants.log`. The mtime check fired for both dbSNP caches.

| Step | Time |
|---|---|
| variants_collect | 0.1 min |
| variants_rsid (bcftools stream 6.5 min) | 7.0 min |
| trans | 0.9 min |
| manifest | seconds |

Variant table after the build:

| Group | Rows | exact | position | none |
|---|---|---|---|---|
| cis | 8,872,723 | 8,826,684 (99.48%) | 45,334 (0.51%) | 705 (0.01%) |
| trans-only, alleles | 309,538 | 307,851 | 1,499 | 188 |
| trans-only, no alleles | 32,765 | 0 | 32,765 | 0 |

The cis row is identical to the 2026-09-03 build. Trans-only rows total 342,303, not the
342,817 estimated from null rsIDs: the other 514 positions are cis variants dbSNP has no
record for. Trans rows without an rsID fell from 817,817 (5.2%) to 2,070 (0.013%).
chr1:82572797 resolves to rs141809548 (exact match, alleles ATGTCT/A, from the trans sQTL file).

Validate: all 18 checks pass, including the three new ones (cis row count, no mixed-allele
positions, trans rsID coverage).

Sizes on disk (before → after): `variants_by_position/` 145 → 165 MB, `variants_by_rsid.parquet`
109 → 128 MB, `trans_pairs/` 438 → 452 MB, `trans_by_variant/` 540 → 554 MB. Net +67 MB,
larger than the +10 MB estimate because the rsID strings for 342k new rows and the filled
`rsid` column in both trans tables compress less than assumed. R2 re-upload when approved is
about 1.3 GB across 71 files. Not uploaded.
