---
date: 2026-10-09
status: abandoned
model: Claude Opus 5.5
description: The general app (ui/) opens on store.json alone and lists the store's collections; a choice puts the study in the URL as a router basename, so the reader underneath is unchanged
---

# The general app's chooser

Step 3 of `plans/2026-10-09-v2-release.md`, under D23 of `plans/2026-10-07-experiment-model.md`
(Sam, 2026-10-07: "nothing should be a default entry point. no particular experiment or collection
loaded on initial load for general ui").

## Where it starts

`ui/src` is a copy of `ui-topchef/src` (mirrored 2026-10-09). The study is fixed at build time:
`VITE_COLLECTION_OF` names an anchor experiment, `openStore` opens that experiment's collection,
`VITE_GWAS` pins the GWAS beside it. Every cache in `lib/store.ts` (`getStore`, `getStoreInfo`,
`wholeGenes`, the `memoBy` per-chromosome and per-part caches, `gwasBins`, about a dozen) is a
module singleton built on that one study, and `StoreProvider` opens it at the app root. Ten literal
internal paths (`/gene/…`, `/variant/…`, `/genes`, `/about`) are written in five files, plus links
built from template strings.

The store after the release holds two QTL collections and one GWAS:

| collection | members | catalog | annotation |
|---|---|---|---|
| TOPCHeF (Murray 2026) | `topchef_{ge,leafcutter}_{cis,trans}` | `topchef_grch38` | `gencode_v34` |
| GTEx v8 heart LV | `gtex_v8_heart_lv_{ge,leafcutter}_cis` | `gtex_v8_heart_lv_grch38` | `gencode_v39` |
| DCM GWAS (Jurgens 2024) | `dcm_jurgens2024_biobanks` | `topchef_grch38` | none |

## Design

1. **The study is a path prefix, `/s/<experiment id>/…`, used as React Router's `basename`.**
   `main.tsx` reads the first two path segments once, before React renders. With `/s/<id>` it sets
   `COLLECTION_OF` to `<id>` and renders the app as it is today under `basename="/s/<id>"`; every
   existing link, `navigate()` and route resolves inside the study with no change. With anything
   else it renders the chooser. Changing study is a full document navigation (`<a href>`, not a
   router link), so every module cache starts empty for the new study and nothing needs a study key.
2. **The key is an experiment id, not the collection name.** The collection name is a content
   digest that changes whenever its metadata does; an experiment id is stable and readable.
   `/s/gtex_v8_heart_lv_ge_cis/gene/ENSG…` opens that experiment's collection, as
   `VITE_COLLECTION_OF` does now. The chooser links to each collection's cis eQTL experiment (or its
   first cis experiment).
3. **The chooser** (`routes/Chooser.tsx`) fetches `store.json`, groups `experiments` by
   `collection`, fetches each collection document (one small JSON object each), and lists per
   collection: title, citation (first author, year, venue, DOI link), source, and the member scans
   by kind and phenotype type. GWAS-only collections are listed under their own heading as data
   available to studies on the same catalog, not as something to open. It is the site's `/`.
4. **GWAS by catalog, not by name.** `VITE_GWAS` goes in `ui/`. A study shows a GWAS when exactly
   one GWAS experiment in `store.json` is on its catalog; with none (GTEx today) the GWAS panels say
   so, as the reader already does when `hasGwas` is false; with several, the first by id and a
   console note (no picker until there is a second GWAS).
5. **`openStore` without the speculative pointers.** `remembered`, `remember`, `speculative` and
   `REMEMBER_KEY` go; the experiment, catalog and annotation pointers are fetched after
   `store.json` names them.
6. **TOPCHeF-specific content out of `ui/`.** The pinned coloc loci (`lib/coloc.ts`, `ColocLoci`,
   Home's coloc track), the DCM wording in `Gene.tsx`, `LocusCompare.tsx` and `LocusPlot.tsx`, and
   About's TOPCHeF text are removed or written from the collection document. The brand in the top
   nav becomes the store's name plus the chosen collection's title, linking back to `/`.
7. **`ui-topchef/` does not change.** It keeps `VITE_COLLECTION_OF`, `VITE_GWAS`, the speculative
   pointers and its content. The two copies stop being byte-identical here.

## Decisions & ownership

| # | Decision | Owner | Note |
|---|---|---|---|
| C1 | The study lives in the URL path as a router basename, and switching study reloads the page | AI-owned, defended | The alternative, a study key on every cache, touches every reader in `lib/` and buys an in-app switch nobody asked for. A shared link carries the study, and local storage would not |
| C2 | The URL key is an experiment id, not the collection name | AI-owned, defended | Collection names are digests and change with metadata edits (R10 in the release plan), which would break every saved link |
| C3 | The GWAS is picked by catalog: the one GWAS on the study's catalog | AI-owned, default | **Silently implied: GTEx has no GWAS panels**, since the DCM GWAS is on `topchef_grch38`. Putting it beside GTEx needs a cross-catalog join through the overlap index, which D24 leaves to the reader and nothing implements |
| C4 | Speculative pointers go | user-owned, taken (D23, 2026-10-07) | **Silently implied: a deep link costs two more round trips** before the first experiment-scoped read (store.json, then the experiment documents, then catalog and annotation). Measurable with the bench; the TOPCHeF app keeps the guess |
| C5 | The chooser lists TOPCHeF as well as GTEx | open / vague | It is in the same store and nothing hides it. Hiding it needs a list in the build or a flag in the store |
| C6 | The coloc track and DCM copy are removed from `ui/`, not generalized | AI-owned, default | Coloc loci are a TOPCHeF publication table, not store data; a store-level coloc object is a later format change |
| C7 | The general app keeps the fault page and the error displays | AI-owned, default | They came with the mirror and cost nothing in production |

## What this changes elsewhere

- **Deploy**: Cloudflare's SPA fallback already serves `index.html` for unknown paths, so `/s/…`
  deep links work without a Worker change; to be confirmed against `ui/wrangler.jsonc`
  (`not_found_handling`).
- **Bench and smoke suites in `ui/`** assume the TOPCHeF paths and numbers. They get the
  `/s/topchef_ge_cis` prefix or stay TOPCHeF-only; the general app's own check is the chooser and
  one GTEx gene page by hand.
- **`store-check` and `store_reference.py` in `ui/`** take `--anchor`, so they run on either
  collection unchanged; the GTEx round trip runs with `--anchor gtex_v8_heart_lv_ge_cis --gwas ''`.
- **C3** makes GTEx gene pages narrower than TOPCHeF's: no GWAS panel, no LocusCompare, no trans tab.
- **C4** costs the general app's first gene page two serial round trips on a cold load.
- **The live general app** (Worker `qtl-browser`, no custom domain) has been broken since 2026-10-07
  and is fixed by this; its URLs change shape, which matters to nobody since it has no users.

## Steps

1. `main.tsx` routing split and `COLLECTION_OF` from the path; `Chooser.tsx`.
2. `openStore` without speculation; GWAS by catalog.
3. Remove the TOPCHeF content; nav brand and About from the collection document.
4. `npx tsc -b`, `vite build`, `store-check` on the GTEx collection of the release store; Sam looks
   at `/`, a TOPCHeF gene page and a GTEx gene page in the dev server.
