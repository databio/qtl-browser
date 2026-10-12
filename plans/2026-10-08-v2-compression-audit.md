---
date: 2026-10-08
status: draft
model: Claude Fable 5.1
description: Shelved audit of the v2 schema for design-level size reductions, with every figure measured on the genome-wide store; nothing here is scheduled
---

# v2 compression audit (shelved)

Shelved 2026-10-08 (Sam: "we can just shelve the audit findings for now"). Kept so the measurements
are not redone. Everything below was measured on `data/derived/store-genome-v2` as built that day;
projections are labelled.

## Where the bytes are

| object | MB | share | stored per unit |
|---|---:|---:|---|
| cis blocks `.qbe` | 3,780 | 83.5% | 6.07 B/row, raw |
| trans frames `.qbt` | 223 | 4.9% | 14.1 B/row |
| GWAS blocks `.qbg` | 171 | 3.8% | 13.7 B/row |
| hits `.qbh` | 120 | 2.7% | 7.3 B/record |
| rsID index `.qbr` | 110 | 2.4% | 12.0 B/record, uncompressed |
| variant catalog `.qbv` | 95 | 2.1% | 10.2 B/site |
| Arrow indexes and the rest | 29 | 0.6% | |

The v2 plan's size tables are MiB labelled MB; 3,604.7 there is 3,780 here.

## Levers, ranked

| lever | measured on | projected saving | cost |
|---|---|---:|---|
| cis `se` to 8 bits **with the effect coded as z**, per-column zstd | 507 blocks, 3.1M rows: 4.68 B/row | 830 MB | the budget is restated: z absolute in SE units, se and beta relative (0.3% half-step on a median block). Coloc is insensitive because z is exact and se enters only through `r = W/(W+se²)` |
| per-column zstd of cis blocks alone | same sample: 5.44 B/row | 380 MB | one inflate per gene page; partial reads survive with per-column lengths in the header |
| rsID index blocked (4,096 records, already defined), columnar, rs delta-coded, zstd | every 50th block: 3.55 B/record | 78 MB | a block offset table in the variant index |
| hits record 20 -> 13 B: vidx as u16 frame offset, value/beta as 16-bit codes per frame, no pad | 1.46M records: 4.22 vs 7.28 B/record | 50 MB | the trans `value` is already a u16 code in the frame; storing it as f32 in hits is more precise than the primary object |
| cis details frame dropped (repeats the search-index row) | index arithmetic: 34 MB incl. cs records and pad | up to 30 MB | blocks stop being self-describing |
| catalog `ma_samples`/`ma_count` dropped | chr22: 2.3 of 10.2 B/site | 21 MB | only if the UI never shows them |
| trans rows drop `rs` and `af` (display-only, recoverable by lookup) | 6 of 18 raw bytes | about 60 MB | a lookup per hover |
| GWAS on a per-source decimal ruler (lossless for fixed-decimal text) | chr21/22, earlier | 4 MB | a second ruler kind |

About 1.0 GB in total, to roughly 3.5 GB, with the precision restatement; about 560 MB without
touching `se`.

## Measured and rejected

- **z instead of beta is size-neutral** (1.79 vs 1.74 B/row after zstd, 127 blocks) but tightens the
  error bound in SE units from a median of 2.4e-04 to 5.8e-05, on 100% of blocks, and makes the
  bound per row (`z_max / 65532`), which closes G1 structurally. Worth doing for precision alone.
- **Delta-coding rsIDs is worse than absolute** everywhere it was tried: GWAS chr22 2.85 -> 3.42
  B/row, catalog chr22 3.00 -> 3.51, trans frames 3.49 -> 3.60. rs numbers are not monotone in
  position and carry about 30 bits of entropy per site.
- **8-bit nlp**: a step of 0.12 in −log10 p on a block whose maximum is 30, a 30% error on a
  displayed p near a threshold. Keep 16.
- **`se` as a residual against the t-predicted value**: 4.43 B/row, the ceiling, and it puts dof
  and the inverse t back on the decode path. Not recommended.

## The number that is not a primitive

Leafcutter is 499.6M of the 623M cis rows, about 3.0 GB of 4.5. Storing nominal rows for all
80,750 introns is a product decision that outweighs every lever above.

## Found on the way

- `results.read_block` is still hardcoded to the v1 block magic.
- The ge trans index's unplaced residue is 16 phenotypes whose genes are absent from GENCODE v34.
