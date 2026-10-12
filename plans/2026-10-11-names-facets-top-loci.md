---
date: 2026-10-11
status: complete
model: Claude Opus 5.5
description: Short collection names, optional experiment facets (with condition), and lead columns in the search index for a top-loci list
---

# Short names, facets, and top loci

Three store changes that land in one pass, then the UI that reads them.

1. **Short collection titles.** The collection `title` becomes the name a list prints ("TOPCHeF
   LV"); the citation and `source` already carry the long description.
2. **Optional facets on each experiment document**, as plain labels: `organ`, `tissue`,
   `cell_type`, `ancestry`, `n_samples`, `assay`, `sex`, `condition`, and `trait` for a GWAS
   (brainstorm section 7; `condition` added 2026-10-11).
3. **`lead_nlp` and `lead_pos` in the cis search index**, so a study's strongest loci can be ranked
   from the object the search box already reads whole. Today the index's only p-like column is
   `sig_value`, the permutation p, which floors: 6,242 TOPCHeF eQTL genes tie at 1e-4 and 3,457
   GTEx eQTL genes at 1e-3, so there is no top 20 to take. The continuous values sit in each
   phenotype's details frame, one range read per gene.

The UI then shows the short names everywhere, the facets on the landing rows, and a **Top loci**
list on the collection page in place of the GWAS track (off since 2026-10-11): a study's 20
strongest genes, or a GWAS's 20 strongest windows from its bin summary.

## Proposed values

Names (Sam, 2026-10-11: the proposed ones):

| Collection | Candidates | Proposed |
|---|---|---|
| TOPCHeF (Murray et al. 2026) | TOPCHeF · TOPCHeF heart · TOPCHeF LV | TOPCHeF LV |
| GTEx v8 via eQTL Catalogue r6 | GTEx v8 heart LV · GTEx heart LV · GTEx LV | GTEx v8 LV |
| Jurgens 2024 DCM, biobanks-only | DCM · Jurgens 2024 · Jurgens 2024 DCM | Jurgens 2024 DCM |

Facets (every one optional; blank means not stated):

| | TOPCHeF (4 experiments) | GTEx v8 LV (2) | Jurgens 2024 DCM |
|---|---|---|---|
| organ | heart | heart | |
| tissue | left ventricle | left ventricle | |
| cell_type | | | |
| ancestry | multi-ancestry | mostly European | European |
| n_samples | blank (Sam). 516 is the likely mapping n by the dof arithmetic, but the paper states no per-scan n | 382 | 937,963 (5,022 cases / 932,941 controls) |
| assay | bulk RNA-seq | bulk RNA-seq | blank: assay describes a molecular phenotype's measurement, which a GWAS has none of |
| sex | both | both | both |
| condition | heart failure (DCM, ICM) + non-failing | postmortem, non-diseased | |
| trait | | | dilated cardiomyopathy |

## Steps

1. **SPEC.**
   - Section 8 Collections: `title` is the collection's display name, short.
   - Section 8 `experiments/<id>.json`: an optional `facets` object holding the nine keys. Values
     are strings, except `n_samples`, which is an integer. Unknown keys are refused.
   - Section 8 Search index: two cis columns.
     - `lead_nlp` (float64): the phenotype's largest −log₁₀ nominal p as its block decodes it.
     - `lead_pos` (uint32): that row's position, the first such row in run order on a tie.
     - Both are null when `n_rows` is 0.
   - Section 14: validate checks the facet keys and types, and checks both columns against a sample
     of decoded blocks.
2. **Builder.**
   - `results._block` returns the run's best `nlp` code and its vidx. The index row carries the
     decoded value, so the index and the block agree by construction.
   - `build_store` reads the facets from config per experiment and writes them.
   - `config.yaml` gets the short titles, and a `facets:` block keyed by experiment id.
3. **Backfill the existing store, with no rebuild.**
   - `results.py add-leads EXP`: decodes every block's `nlp` column (the first column, at offset
     72). It rewrites the index with the two columns, and rewrites its parts and pointer. This
     follows the `migrate-index` / `add-split` precedent.
   - `qtlstore set-collection KEY`: writes the new collection document from config, repoints that
     key's experiments, and rewrites `store.json`.
   - `qtlstore set-facets`: writes the facets into the pointers from config.
   - Test: on the fixture store, build-with-columns and build-then-backfill give identical
     index bytes.
4. **UI (`ui/`).**
   - `shortName` returns the collection title.
   - Landing rows show the facet union: tissue, ancestry, n, condition.
   - **Top loci** on the collection page:
     - **Study:** genes ranked by their best `lead_nlp` over eQTL and sQTL, one per 1 Mb window, top
       20. Each row shows gene, type, lead position, −log₁₀ p, and links to the gene page.
     - **GWAS:** bins ranked by `min_p`, top 20. Each row shows the window, lead rsID and p, and
       links to the region or variant page.
5. Update the smoke suites' digests, the README and the brainstorm doc. Validate the backfilled
   store and run `store-check`.

## Decisions & ownership

| # | Decision | Tag | Reason / consequence |
|---|---|---|---|
| D1 | Collection names: TOPCHeF LV, GTEx v8 LV, Jurgens 2024 DCM | user-owned, surfaced | Sam, 2026-10-11. Each title change makes a new collection digest |
| D2 | Jurgens `assay` left blank | user-owned, surfaced | Sam, 2026-10-11: `assay` reads oddly for a GWAS |
| D11 | TOPCHeF `n_samples` left blank | user-owned, surfaced | Sam, 2026-10-11: the paper states no per-scan n; 516 is inferred, not stated |
| D3 | `condition` is an optional facet | user-owned, surfaced | Sam, 2026-10-11 |
| D4 | Rank by lead nominal p, not `p_beta` | AI-owned, defended | Computable from nominal rows for any source; `p_beta` exists only where the source ran permutations |
| D5 | `lead_nlp` is the block's decoded value, not the source float | AI-owned, defended | Index and gene page then never disagree in the last digit |
| D6 | Backfill commands instead of a full rebuild | AI-owned, defended | Decodes stored blocks locally (minutes) rather than a Rivanna rerun; config still changes so a later rebuild yields the same bytes, which the fixture test checks |
| D7 | No format version bump for the two index columns | AI-owned, defended | Additive: the UI reads missing columns as null. Validate will require them, so an un-backfilled v2 store fails validate until `add-leads` runs |
| D8 | Top loci: one gene per 1 Mb, 20 rows, eQTL and sQTL pooled | AI-owned, default | Keeps one locus (e.g. HLA) from filling the list; 1 Mb matches the cis window |
| D9 | Top loci replaces the track's spot on the collection page | **user-owned, surfaced** | Placement is Sam's; the list is the track's "see the important loci" job without the genome plot |
| D10 | Facets stay off collections; a landing row shows the union of its experiments' values | user-owned, surfaced | Brainstorm section 7 |

## What this changes elsewhere

- **URLs:**
  - New titles mean new collection digests, so every `/c/<digest>` link changes, and old links get
    the not-found page (Sam's policy).
  - The smoke suites and `locus-bench.mjs` hardcode the TOPCHeF and GTEx digests.
  - The deployed site is unaffected until upload.
- **Object store:**
  - The old collection documents and old index objects become unreferenced; `gc` removes them
    locally.
  - On B2 they stay until the separate prune step, which needs its own go-ahead. **R2 stays frozen
    and untouched.**
- **Index size:**
  - Two columns over ~19k–80k rows per cis experiment, about 0.3–1 MB uncompressed and less after
    zstd.
  - The search box reads every cis index whole, so the first search gets a little slower.
  - This needs measuring on the TOPCHeF sQTL index (1.7 MB today).
- **ui-topchef:**
  - It reads the same store, but it pins its study by an experiment id
    (`VITE_COLLECTION_OF=topchef_ge_cis`), not by a digest. So the rename doesn't break it. It
    would only show the new title wherever it prints one.
  - It ignores the new index columns.
- **Validate:**
  - Stricter. A store built before this plan fails validate until it is backfilled.
- **Coloc lists:**
  - These name experiments, not collections, so they are unaffected.

## Implementation log

**Pipeline.**
- `results.py`:
  - `lead_nlp`/`lead_pos` are in `INDEX_SCHEMA` and `CIS_INDEX_FIELDS`.
  - `_block` takes them from the codec's own quantize/dequantize of the run's p (`run_lead`).
  - `block_leads` reads them back out of stored blocks.
  - `add_leads` (CLI `add-leads`) rewrites the index through `_rewrite_index`, which `migrate_index`
    now shares.
  - `check_leads` is validate's sample check.
  - `build(..., facets=)` applies `qtlstore.with_facets` per experiment.
- `qtlstore.py`:
  - `FACET_KEYS`, `check_facets`, `with_facets`.
  - Validate checks facets and leads (SPEC section 14, items 14 and 15).
  - `set_collection` / `set_facets` (CLI `set-collection KEY`, `set-facets`). Members are a GWAS
    whose id is the key, or QTL scans whose `source.experiment_id` is.
- `config.yaml`: short titles, and a `facets:` block (YAML anchors per study).
- SPEC sections 8, 14 and 15. `pipeline/README.md` documents the commands.
- Tests: `test_lead_columns_match_blocks_and_backfill` (the backfill lands on the builder's index
  bytes) and `test_facets_and_set_collection`. All 12 suites pass.

**Release store backfilled, no rebuild.** The pointers from before are in
`data/derived/_pointers_before_2026-10-11/`.
- `add-leads` on the four cis scans took 5–16 s each. Index sizes:

  | Scan | Before | After | Change |
  |---|---|---|---|
  | TOPCHeF eQTL | 0.56 MB | 0.76 MB | +37% |
  | TOPCHeF sQTL | 1.72 MB | 2.53 MB | +47% |
  | GTEx eQTL | 0.60 MB | 0.84 MB | +41% |
  | GTEx sQTL | 2.62 MB | 2.67 MB | +2% |

  GTEx sQTL barely changes because its nominal file is fine-mapped introns only, so most introns
  have no lead.
- Then `set-collection` for the three keys, then `set-facets`. New digests:

  | Collection | Digest |
  |---|---|
  | TOPCHeF LV | `mm4EHuZR-o3yuq37E-LKhkWIS1s6kOzr` |
  | GTEx v8 LV | `0XIkm9wOG1L6qTyHou7h5WljxlJW6ITt` |
  | Jurgens 2024 DCM | `xv3pLGNS6VlRt2M5wouPF-IPrU6c3aYM` |

- Validate passes with no notes, and `store-check` passes for TOPCHeF + DCM and for GTEx.

**UI.**
- `shortName` is the collection title.
- `facetValues` / `facetLine` drive the landing rows (citation on line 1, facets on line 2) and the
  collection header.
- `gene-index.ts` reads the lead columns when present (`ifPresent`, so an older store still opens)
  and `topLoci` ranks them.
- `components/TopLoci.tsx`: both lists are drawn as the whole-genome track, with markers only (no
  bins) and no table (Sam, 2026-10-12).
  - `StudyTopLoci`: at most 20 genes, one per 1 Mb, eQTL and sQTL pooled. Each marker is labelled
    with its symbol, colored eQTL or sQTL, and opens the gene.
  - `GwasTopLoci`: the lead variants of the 20 strongest bins. Each marker is labelled with its rsID
    only and opens UCSC.
- The smoke suite checks Top loci instead of the removed picker text. Type-check, lint and build are
  clean.

**Measured on the store.**
- TOPCHeF's top 20 is 18 sQTL. PQBP1 (−log₁₀ p 204) leads; the strongest eQTL is DDX11 (166).
- GTEx's top 20 mixes the two, DDX11 eQTL first (152).
- Pooled, splicing dominates a study whose intron p run stronger. Showing eQTL and sQTL separately
  would change that.

**The collection page now reads the whole-genome genes table and every cis index on load**, the
same reads the search box makes on first use: about 3.3 MB of index for TOPCHeF, which was 2.3 MB
before.
