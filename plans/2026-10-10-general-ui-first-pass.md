---
date: 2026-10-10
status: in-progress
model: Claude Opus 5.5 (1M context)
description: First pass of the general UI (ui/): a collection browser around today's pages, unpinned from TOPCHeF and DCM, so Sam can explore the store before the design is refined
---

# General UI, first pass

The aim is to see what the store holds and what each page can show for every collection, not a
finished design. Sam will refine routing, layout and workflow after using it, and expects that to
change much of this. So: the smallest change that turns `ui/` from a pinned TOPCHeF copy into a
browser of every collection, reusing today's pages wherever they already work.

Decisions this rests on are in `plans/2026-10-09-general-ui-brainstorm.md`: pages belong to a
collection, collection digests in URLs, experiments by id, any GWAS for any study joined by
`(pos, ref, alt)`, one sequence collection per store, no default experiment, public audience, no
shared code with `ui-topchef/`.

## What the visitor sees

**Landing `/`.** Root collections as row cards in two groups, Studies and GWAS (the GWAS are a layer
any study is read against, not peers of the studies). Line 1: short name ("Murray 2026", from first
author and year), venue, and a right-aligned summary ("eQTL, sQTL · cis, trans" or "12,504,079
variants · 5,022 cases"); line 2: the full title, truncated. Nothing is opened until a row is
clicked. Built from `store.json` and the collection documents alone.

**Names.** Lists, crumbs and the GWAS picker use the short name; the full title appears only as the
collection page's description. Crumbs are links: Collections leads back to the landing page, and
gene, variant and region pages lead with their collection.

**GWAS picker.** A dropdown of short names, each with its size and a "coloc list" tag where a list
pairs it with this collection, and the full title dimmed below. A filter box appears past 6 GWAS.

**Collection `/c/<digest>`.**

- Header: title, citation, source.
- Experiments table: kind, phenotype type, phenotypes tested, significant, the significance rule.
- Genome track with a GWAS picker. A collection's own loci when it is a GWAS
  collection; for a QTL collection the chosen GWAS's loci, marked from a coloc list when one names
  the pair, otherwise with the collection's significant genes near each locus.
- Link to the gene table.
- A parent collection lists its children instead.

**Inside a QTL collection**, today's pages under the prefix:

- `/c/<digest>/genes`: the gene table, a tab per experiment.
- `/c/<digest>/gene/<id>`: a tab per experiment that tested the gene; coloc and LocusCompare
  against `?gwas=`.
- `/c/<digest>/variant/<id>`, `/c/<digest>/region/<loc>`.
- A breadcrumb (store > collection > page) and a search box scoped to the collection.

**Elsewhere:** `/about` (the browser and the store, no TOPCHeF story); a digest not in `store.json`
gets a "not in this store" page with the collection list.

**Not in this pass:** facet labels and filters (SPEC and config work first), store-wide search,
presence across collections, TOPCHeF's colocalized-loci prose.

## Steps

1. **Store index.** `getIndex()`: `store.json` plus every collection document its experiments name
   and their `parent` chains, one fetch each. Feeds the landing page, collection pages and the
   breadcrumb.
2. **Store per collection.** `openStore()` takes a collection digest instead of
   `VITE_COLLECTION_OF`, memoized per digest; every memo that calls `getStore()` today
   (`wholeGenes`, `chromGenes`, `experimentPart`, `gwasBins`, `colocLists`, ...) is keyed by it.
   A collection must hold QTL experiments on one catalog and one annotation, as today; one that
   does not gets an error page naming why. The remembered-ids guess in `localStorage` is keyed by
   digest. `VITE_COLLECTION_OF` and `VITE_GWAS` leave `ui/` (`store.ts`, `.env.production`,
   README).
3. **GWAS by choice.** `/c/<digest>/gwas/<experiment id>` opens that GWAS's document and
   collection. No GWAS is chosen until the visitor picks one,
   however many the store holds; until then the GWAS panels show the picker and nothing else.
   `openStore` stops refusing a GWAS on another catalog.
4. **Window join.** Where coloc, LocusCompare and the gene page's GWAS panel read GWAS rows by
   vidx, a GWAS on another catalog is matched by `(pos, ref, alt)`: the GWAS window rows, the GWAS
   catalog's pages for the window, a map from site to row. The same catalog keeps the vidx path.
   Checked against TOPCHeF × DCM (both paths must agree there) and run on GTEx × DCM.
5. **Routes and links.** The routes above; a `useCollectionPath()` helper so every internal link
   carries the prefix; the not-found page for an unknown digest.
6. **Pages.** New: landing, collection, parent collection, GWAS-only collection (track only).
   Reused: Genes, Gene, Variant, Region, with TOPCHeF-only text removed. Experiment tabs are
   labelled from `phenotype_type` (`ge` as eQTL, `leafcutter` as sQTL, anything else as itself).
   Coloc chips and markers come from the lists that name the shown experiment and the chosen GWAS.
7. **Checks.** `npx tsc -b`; `store-check` over every collection rather than one; the smoke suites
   in `ui/bench` moved to the new routes and run for TOPCHeF and GTEx against Sam's dev server; the
   cold-load bench for one TOPCHeF gene page against today's numbers (29 requests, 0.5 MB).

## Decisions & ownership

| # | Decision | Ownership | Reason / consequence |
|---|---|---|---|
| D1 | A collection opens only if its QTL experiments share one catalog and one annotation | AI-owned, defended | True of every collection in the store, and the reused pages assume it. A future collection mixing catalogs (D24 allows it) gets an error page until pages handle it |
| D2 | No GWAS is chosen until the visitor picks one, even when the store holds one | user-owned, surfaced (Sam, 2026-10-10: "no preselection of any experiment anywhere") | The track, coloc and LocusCompare panels start as a picker; a shared link carries `?gwas=` |
| D3 | TOPCHeF-only text removed from `ui/`, and generic pages get short new copy (card text, collection header, About) | **user-owned, surfaced** | Sam objected to unrequested wording changes before. New copy is minimal and data-driven, listed in the report for Sam to rewrite |
| D4 | Facets left out of this pass | user-owned, surfaced | Sam asked for a simple first pass. Cards show only what the store already holds |
| D5 | Tab labels: `ge` as eQTL, `leafcutter` as sQTL, others raw | AI-owned, default | Needed to label tabs at all; a label map is one line to change |
| D6 | GWAS-only collection page shows the track and nothing else | AI-owned, default | A GWAS has no gene pages; anything more is a design question for after the first look |
| D7 | `ui/bench` smoke suites rewritten for the new routes; `ui-topchef/` untouched | AI-owned, defended | Its suites pin routes that no longer exist in `ui/`; no shared code (brainstorm section 8) |
| D8 | The window join is written in `ui/` only | silently implied | Code sharing was skipped, so `ui-topchef` keeps its same-catalog path and never gets the join |
| D9 | The GWAS is a path segment, `/c/<digest>/gwas/<id>/...`, not `?gwas=` | **AI-owned, changed from the agreed scheme** | The collection and GWAS prefix is the router's basename, so all ~50 existing links keep it without edits; a query parameter would be dropped by every link. Brainstorm section 5 said `?gwas=` |
| D10 | One collection and one GWAS per page load; changing either is a full navigation | AI-owned, defended | The readers cache per page load by chromosome, gene and vidx; scoping the load is what keeps two collections from mixing without rekeying every cache |
| D12 | Landing as row cards grouped Studies / GWAS, short names derived from the citation | **user-owned, surfaced** | Sam chose rows (B) grouped by role (C) over cards and a table. The short name is derived, not stored: a `label` field would change every collection digest |
| D13 | GWAS genome tracks off: the collection page keeps only the GWAS picker; the track components stay in `GwasTrack.tsx` unrendered | **user-owned, surfaced** | Sam: the track reads as tacky and is hard to make general. D11 is moot while it is off |
| D14 | The GWAS picker leaves the collection page; it stays only where it is used, the gene page's coloc and LocusCompare panels. Its final home is open | **user-owned, surfaced** | Sam: not on the collection page; placement undecided |
| D15 | One nav everywhere: Collections, About. Context inside a collection is the crumbs (collection › Genes › gene) | **user-owned, surfaced** | Sam: the in-collection nav (Collection, Genes, Collections, About) was neither universal nor clear |
| D16 | A section or tab exists only for scans the study has: no trans section without that trans scan, no eQTL or sQTL tab without that cis scan | **user-owned, surfaced** | Sam. Replaces the "not available" notices |
| D11 | Without a coloc list, the track marks each genome-wide significant bin's nearest significant gene within 1 Mb of the bin's lead | AI-owned, default | Reads the whole-genome genes and cis indexes the search box reads. One marker per bin, so a bin with two loci shows one |

## What this changes elsewhere

- **`ui-topchef/`**: nothing. Its pinned collection, GWAS, routes and wording stay.
- **The store**: nothing. Same `store-release-v2`, same B2 upload; no SPEC or pipeline change.
- **Deploy**: `deploy-ui.yml` unchanged; `ui/.env.production` loses `VITE_COLLECTION_OF` and
  `VITE_GWAS`. The general Worker serves the new app once PR #5 merges.
- **Links**: anything pointing at the general Worker's old `/gene/:id` routes breaks. Nothing public
  links there (its Worker has been broken since 2026-10-07).
- **Requests**: the landing page costs `store.json` plus one small JSON per collection (3 today).
  A gene page against a GWAS on another catalog adds the GWAS catalog's window pages, a few small
  reads. A study page with no GWAS chosen reads no GWAS objects.
- **Release plan** (`plans/2026-10-09-v2-release.md`): the merge waits for this pass and Sam's look
  at it; refinements after that may land in PR #5 or a follow-up, Sam's call.
- **Dev server**: Sam runs `QTL_DATA_DIR=../data/derived/store-release-v2 VITE_DATA_BASE= npm run dev`
  in `ui/` as before; no new settings.

## Implementation log

2026-10-10, uncommitted.

- **Scope** (`src/lib/scope.ts`): the collection digest and chosen GWAS are parsed from the URL
  once; their prefix is `BrowserRouter`'s basename (`main.tsx`). `chooseGwas` reopens the page.
- **Store** (`src/lib/store.ts`): `VITE_COLLECTION_OF` and `VITE_GWAS` are gone. `storeDoc()`
  (one `store.json` fetch), `getIndex()` (every experiment, collection with parents, and coloc
  pointer), and `openStore()` opening the scope's collection. The refusal of a GWAS on another
  catalog is gone: GWAS blocks carry positions and alleles, and coloc and the locus SQL already
  aligned on them, so no new join code was needed; the catalog's `collection_digest` must equal
  `store.json`'s. Pointers are fetched once per page load. `gwasBinsOf(g)` replaces `gwasBins()`.
- **Context** (`src/contexts/store-context.tsx`): every page gets the index; a QTL collection's
  study opens beside it, not after it.
- **Pages**: `routes/Landing.tsx` (grouped rows), `routes/Collection.tsx` (study, GWAS-only and parent
  collections, the not-in-store page), `routes/About.tsx` rewritten for the store,
  `components/GwasTrack.tsx` (picker, study track, GWAS-only track) replacing `ColocLoci.tsx`;
  `routes/Home.tsx` removed. `lib/describe.ts` holds the short-name, summary and table wording helpers.
- **Wording made data-driven**: gene chips ("coloc · eQTL", tooltip from the list), coloc section,
  LocusCompare axis and empty states, the not-tested and variant-not-found messages, the Genes
  header, the nav wordmark ("qtl.browser").
- **Chromosome sizes** by `store.json`'s `collection_digest` rather than a build name.
- **Checks**: `npx tsc -b`, `oxlint` (no new warnings beyond `only-export-components` on the
  context's hooks), `vite build` clean. `store-check` passes on the release store for TOPCHeF +
  DCM, GTEx alone and GTEx + DCM; its seqcol check now compares `store.json`'s collection with the
  catalog's. Smoke suites and the bench moved onto `/c/<digest>/gwas/dcm_jurgens2024_biobanks/...`;
  the gene suite gains landing, collection, About and GTEx (no GWAS, then against DCM) checks.
  **Not run**: the smoke suites and bench need Sam's preview server.
