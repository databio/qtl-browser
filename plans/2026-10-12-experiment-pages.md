---
date: 2026-10-12
status: complete
model: Claude Opus 5.5
description: The general UI's working unit becomes one experiment: /c/<digest>/e/<id> pages, with the collection as a grouping and citation page
---

# Experiment pages

The general UI's search, gene, variant and region pages currently open a whole collection and join
its cis scans by gene. That join assumes at most one cis scan per phenotype type in a collection,
and the store reader throws otherwise. A collection that spans tissues or ancestries breaks it:
GTEx with several tissues, ARIC EA/AA, MESA.

Sam, 2026-10-12: the working unit becomes the experiment, and the collection keeps its digest in
the URL.

## Routes

```
/                                    landing: studies and GWAS
/about
/c/<digest>                          collection: citation, facets, its experiments (links)
/c/<digest>/e/<id>                   experiment: its landing (search, top loci, summary)
/c/<digest>/e/<id>/genes             its genes
/c/<digest>/e/<id>/gene/<g>, /variant/<v>, /region/<loc>
/c/<digest>/e/<id>/gwas/<gwas id>/…  the same, with a GWAS for coloc and LocusCompare
```

The router's basename is `/c/<digest>` on a collection page and `/c/<digest>/e/<id>[/gwas/<g>]`
inside an experiment. Moving between experiments, collections or GWAS is a full navigation, as now.

## What each page holds

- **Collection:**
  - shows the header (name, facets, citation, DOI) and an experiments table whose rows link to the
    experiment pages;
  - loses its search and Top loci.
- **QTL cis experiment:**
  - the landing page has search, Top loci from this scan's lead columns, and a summary line (size,
    significance rule, facets);
  - gene, variant and region pages show this scan only. The gene page's eQTL/sQTL tab switch
    disappears on its own (D16 of the first-pass plan), since there is one cis scan;
  - the trans section is hidden, because the experiment holds no trans scan.
- **QTL trans experiment:**
  - search runs over the annotation's genes;
  - the gene page shows the gene's trans associations of this scan;
  - the variant page shows only its trans associations;
  - no Top loci (there are no cis leads) and no locus plot.
- **GWAS experiment:** the landing page has its Top loci track and summary. It has no gene, variant
  or region pages.

## Decisions & ownership

| # | Decision | Tag | Reason / consequence |
|---|---|---|---|
| D1 | Pages are scoped to one experiment | user-owned, surfaced | Sam. Removes the one-cis-scan-per-type constraint |
| D2 | The collection digest stays in the URL | user-owned, surfaced | Sam. A collection edit (title, citation) breaks experiment links too; old links get the not-found page |
| D3 | The GWAS stays a path segment, after the experiment | AI-owned, defended | Same reason as before: it is part of the basename, so in-app links keep it |
| D4 | The combined gene page goes: an eQTL page no longer shows the gene's sQTL introns, and the reverse | **silently implied, surfaced** | The cost of D1. Seeing both means switching experiments. A "sibling scan" link could come back later |
| D5 | A trans experiment searches the annotation's genes, and its gene and variant pages show trans associations only | AI-owned, default | The smallest page a trans scan supports. No Top loci: trans scans have no lead columns |
| D6 | A GWAS experiment page is its Top loci track and summary | AI-owned, default | It has no gene-level data of its own |
| D7 | Old `/c/<digest>/gene/…` URLs get the not-found page | AI-owned, defended | The redirect-free policy for stale links (2026-10-10) |
| D8 | The smoke suites move to per-experiment URLs: sQTL checks go to the sQTL experiment | AI-owned, default | Not run here: they need Sam's server |

## What this changes elsewhere

- **`lib/scope.ts`:**
  - parses the experiment and builds the basename;
  - `chooseGwas` keeps the experiment.
- **`lib/store.ts` `openStore`:**
  - the study is one experiment;
  - its collection must match the URL's;
  - the speculative-open key includes it.
- **`StoreProvider`:** opens the study only inside a QTL experiment.
- **Gene page:**
  - with no cis scan (a trans experiment), shows the trans section instead of "not tested";
  - the load effect fetches the trans table for such a gene.
- **Variant page:** with no cis scan, hides the cis sections.
- **Crumbs:** Collections › collection › experiment › … (the collection is now a full navigation).
- **Bench:** `locus-bench.mjs` and the three smoke suites rebuild their prefixes.
- **Unchanged:** ui-topchef, which has its own app, and the store.

## Implementation log

- **Scope and store.**
  - `lib/scope.ts` parses `/c/<digest>[/e/<id>[/gwas/<g>]]`; `collectionPath(digest, experiment, gwas)`.
  - `openStore` opens the one experiment the URL names, refuses one outside the URL's collection, and
    keys the speculative open by it.
  - `StoreProvider` opens the study only for a QTL experiment of the collection. New hook:
    `useExperiment`.
- **Pages.**
  - `routes/Collection.tsx` lists the experiments as links (a full navigation) and loses search and
    Top loci.
  - New `routes/Experiment.tsx`: the header shows the collection name, the scan, the experiment's
    facets, its size and its rule. A QTL scan gets search, a cis scan Top loci, a GWAS Top loci.
  - `App.tsx` routes the experiment's pages behind `Study`, which now asks about the experiment.
  - Crumbs: Collections › collection › experiment › …
- **Trans experiment.**
  - The gene page loads and shows the gene's trans table(s) instead of "not tested".
  - The variant page hides its cis sections.
- **Genes page.** Filters follow the scan: eGenes for an eQTL scan, sQTL genes for an sQTL scan,
  "All" for a trans scan.
- **Fixes the split exposed.**
  - The "no eQTL test" chip shows only in an eQTL experiment; it would otherwise mark every gene of
    an sQTL scan.
  - The GWAS picker's coloc tag keys on the experiment.
- **Bench.**
  - The three smoke suites split their checks across the four TOPCHeF experiments. Combined counts
    are summed over the pages: rs34599497's 50 credible sets over the two cis scans, rs4819361's 10
    trans rows as 2 + 8.
  - `locus-bench.mjs` opens `tab=sqtl` pages in the sQTL experiment. Its `pair` scenario between an
    eQTL page and an sQTL page is now a full navigation, so it measures something different from
    before.
- **Smoke suites against Sam's preview (2026-10-12):** all pass. Gene 49/49, trans 15/15 (with
  `SMOKE_REF`), variant 33/33. rs34599497's credible sets split 2 + 48 between the eQTL and sQTL
  experiments. The bench was not run. Type-check, lint and build are clean.
