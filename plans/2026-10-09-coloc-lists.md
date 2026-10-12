---
date: 2026-10-09
status: in-progress
model: Claude Opus 5.5
description: Optional per-pair coloc lists in the store (membership only), replacing the TOPCHeF UI's hardcoded coloc gene list
---

# Coloc lists in the store

Decided with Sam on 2026-10-09 (`plans/2026-10-09-general-ui-brainstorm.md`, "Coloc lists"): no coloc
list is hardcoded in either UI; the pipeline never runs coloc; a list is stored only for a pair
someone declares; a list is membership only, since the browser computes PP.H4 and the lead variant.
Goes in PR #5.

## Today

`ui-topchef/src/lib/coloc.ts` (identical in `ui/`): `COLOC_EQTL_GENES` (21 symbols),
`COLOC_SQTL_GENES` (4), and `COLOC_LOCI`, each symbol's gene id, chromosome and TSS pinned to GENCODE
v34's identity digest. Read by `ColocLoci.tsx` (Home track markers, panel copy), `Home.tsx` (count
sentence), `Gene.tsx` (the "DCM coloc" chips), `About.tsx` (prose), and `scripts/store-check.ts`
(the pin check). The authors will not supply their table; this list is all there is.

## Design

**Input.** Coloc results as commonly written, not produced by the pipeline (Sam, 2026-10-09): a
CSV or TSV in `pipeline/colocs/` with whatever of these columns it has, matched by name:

| meaning | accepted names |
|---|---|
| phenotype | `phenotype_id`, `molecular_trait_id`, `trait_id` |
| gene | `gene_id`, `ensembl_gene_id`, `ensembl_id` (versions stripped) |
| symbol | `symbol`, `gene_name`, `gene_symbol`, `hgnc_symbol` |
| PP.H4 | `PP.H4.abf`, `PP.H4`, `pp_h4`, `PPH4`, `H4` |

A row needs a phenotype, a gene or a symbol. A symbol is resolved through the QTL experiment's
annotation and must name exactly one gene there; a phenotype without a gene takes its gene from the
experiment's search index. With a PP.H4 column the list keeps the rows at or above the entry's
`min_pp_h4`; without one, every row is a result as the source reported it. Only membership is
stored. A `colocs:` entry in `config.yaml` per list names the QTL experiment, the GWAS experiment,
the file, `source` (free-form), `threshold` (text, e.g. `PP.H4 > 0.8`), optional `min_pp_h4`, and
`note`. TOPCHeF gets two lists, from the authors' symbols as given:

| id | QTL experiment | GWAS | rows |
|---|---|---|---|
| `murray2026_ge_dcm` | `topchef_ge_cis` | `dcm_jurgens2024_biobanks` | 21 symbols, gene level |
| `murray2026_leafcutter_dcm` | `topchef_leafcutter_cis` | `dcm_jurgens2024_biobanks` | 4 symbols, gene level (which 17 introns is not recoverable) |

**Store.** A new pointer level, `colocs/`, in `qtlstore.POINTER_DIRS`, so `store.json` lists it and
`gc`, `validate` and the upload order pick it up unchanged. The pointer:

```json
{"id": "murray2026_ge_dcm", "qtl": "topchef_ge_cis", "gwas": "dcm_jurgens2024_biobanks",
 "source": {...}, "threshold": "PP.H4 > 0.8", "note": "...",
 "object": "<digest>.json", "n": 21}
```

The object is a content-addressed JSON array of `{gene_id, phenotype_id}` rows, sorted. Written by
`python -m pipeline.coloc add --store S [--id ID]` from the config entries (all whose two experiments
are in the store when no id is given), and by `build_store` and `store.sbatch` the same way.

**Validate** (new check 15): `qtl` names a `qtl_cis` or `qtl_trans` experiment in the store and
`gwas` a `gwas` one; the object decodes to rows with unversioned Ensembl `gene_id`s; every
`phenotype_id` is a phenotype of the QTL experiment, and every `gene_id` the gene of at least one of
its phenotypes (so a rebuilt experiment that lost a phenotype fails the list rather than leaving it
silently wrong); rows are unique; `n` equals the row count.

**UI** (`ui-topchef`, then copied to `ui/` as before, since `ui/src` is still the same code until the
redesign):

- `lib/store.ts`: `colocLists()` fetches the `colocs` pointers `store.json` lists whose `qtl` is in
  the open collection and whose `gwas` is the open GWAS, then their objects; memoized.
- `ColocLoci.tsx`: markers from the lists' gene ids, placed through the annotation's gene lookup
  (one small read per gene, as the non-pinned path does today); trait from the list's QTL
  experiment's phenotype type; panel title and description from the list (`threshold`, GWAS
  collection title, `source`). No list for the pair: the panel shows the GWAS bars with no markers
  and says no coloc list is stored for this pair.
- `Gene.tsx`: the chips from list membership by gene id, tooltip from the list.
- `Home.tsx`: the count sentence from the lists' distinct gene counts.
- `About.tsx`: the "Colocalized loci" item describes the list from its pointer (`source`,
  `threshold`, `note`); the PJVK/CDKN1A remark moves into the list's `note`.
- `lib/coloc.ts` deleted. `store-check.ts` loses the pin check and gains a check that every list
  object decodes and its gene ids resolve through the annotation lookup; `store_reference.py` emits
  the lists.
- Smoke suite: `PP\.H4 > 0\.8, 22 loci` follows the new panel copy.

## Decisions & ownership

| # | Decision | Owner | Note |
|---|---|---|---|
| K1 | Lists live in the store, optional, per declared pair, membership only | user-owned, taken (Sam, 2026-10-09) | |
| K2 | The source of truth is a coloc results file as commonly written, in `pipeline/colocs/` (tracked) | user-owned, taken (Sam, 2026-10-09) | The TOPCHeF list exists nowhere downloadable, so it is tracked; a larger list can sit under `data/raw/` with the config naming it |
| K3 | **No identity pinning of the pair; validate checks the rows against the experiments instead** | AI-owned, defended | The brainstorm said the pointer carries the experiments' catalog identity digests. A catalog digest does not change when a scan is rebuilt on the same sites, so it would not catch the stale case anyway; checking every row's phenotype and gene against the QTL experiment's index does |
| K4 | No version number changes | user-owned, taken (Sam, 2026-10-09: "v2 is not published or released anywhere") | |
| K5 | sQTL rows at gene level | forced by the source | The track and the chips are gene-level today, so nothing is lost relative to now |
| K6 | `ui/` gets the same change in this PR | AI-owned, defended | Sam: "drop the hardcoded list in the ui"; `ui/src` is the same code until the redesign |
| K7 | The GWAS bin summary is computed by `gwas.build` for every GWAS, from the oriented rows | user-owned, taken (Sam, 2026-10-09: "fully generic") | `gwas_bins.parquet`, `dcm_gwas.BINS_SQL`, `--bins-only` and `gwas_bin_bp` go; `lead_ea` becomes `lead_alt`; ties go to the first row in stored order. DCM's bin object changes digest; `min_p` and `n_gws` are unchanged in all 569 bins |

## What this changes elsewhere

- **Store**: no object or experiment pointer changes digest. `store-release-v2` gains two pointers
  and two small objects; the upload dry run grows by four keys.
- **SPEC**: section 20 for `colocs/`; section 2's pointer list; section 14's check 15.
- **`qtlstore`**: `POINTER_DIRS` gains `colocs`, so `remove-experiment` leaves a list naming a
  removed experiment, which then fails `validate` (by design: a list about a missing experiment is
  stale; remove it by hand like a catalog pointer).
- **TOPCHeF Home**: markers now need one gene lookup each instead of the pinned table (22 small
  reads, parallel, after the annotation lookup header); the pinned path existed to avoid that.
- **About copy** changes wording; the smoke suites' panel text assertion changes.

## Steps

1. Pipeline: `POINTER_DIRS`, `add-coloc`, validate check, `build_store` hook, config entries and the
   two TSVs, tests in `test_qtlstore.py`. Add both lists to `store-release-v2`; validate.
2. SPEC section 20 and the touched sections.
3. UI (`ui-topchef`): `colocLists()`, the four components, `coloc.ts` deleted, store-check and
   reference, smoke text. `npx tsc -b`, `vite build`, `store-check` on the release store. Copy to `ui/`.
4. Sam looks at Home and a coloc gene page in the dev server.

## Implementation log (2026-10-09)

**Generic GWAS bins (K7).** `gwas.bin_rows` / `build_bins` compute the summary from the rows
`gwas.build` already holds; `bin_bp` is a parameter (default `BIN_BP` 5 Mb). Removed:
`dcm_gwas.BINS_SQL`, `bins()`, `--bins-only`, `gwas_bin_bp` in config, the contract table and its
`contract_check` rule, the README row. `BINS_SCHEMA`'s `lead_ea` is `lead_alt` (and in
`store-decode.ts`). `test_gwas_bins` rewritten (fixture rows, narrow bins, a tie). TOPCHeF rebuilt
into `store-release-v2` (`results build`, 11 min, log `data/derived/build_16_release_generic_bins.log`):
every TOPCHeF experiment pointer byte-identical to `store-genome-v2`; the DCM pointer differs only in
`bins` (569 bins, 15,650 bytes against 16,507); the old bin object removed by `gc`.

**Coloc lists.** `pipeline/coloc.py` (`resolve`, `add`, `add_configured`, `check`, CLI `add`/`show`);
`qtlstore.COLOCS` in `POINTER_DIRS` and validate check 15; hooks in `build_store` and `store.sbatch`;
SPEC sections 2 and 20; `config.yaml` `colocs:` with the two lists from the authors' symbols in
`pipeline/colocs/`. Every symbol resolved to one GENCODE v34 gene; validate PASS on the release store.
`test_coloc_lists` added; all 12 pipeline suites pass.

**UI.** `colocLists()` in `lib/store.ts`; `lib/coloc.ts` holds helpers and `useColocLists()`, no
genes. The TOPCHeF wording is unchanged (Sam: restore the original wording everywhere); only what
the list decided now comes from the store: `ColocLoci.tsx` places markers by gene id through the
lookup and counts the loci, `Gene.tsx` decides the "DCM coloc" chips and the ColocSection's
"Reported as colocalized in the preprint" note by list membership, and `Home.tsx` and `About.tsx`
take the 21 and 4 from the lists. `store-check`'s pin
check replaced by a coloc list round trip; `store_reference.py` emits the lists with their genes'
lookup rows; the Home smoke assertion follows the new copy. Mirrored to `ui/` (`src` identical;
scripts and bench with `ui/` paths). `npx tsc -b` and `vite build` pass in both; `store-check` passes
on both collections of `store-release-v2`.

Remaining: Sam's look in the dev server (step 4); the smoke suites are unrun.
