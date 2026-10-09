---
date: 2026-10-09
status: draft
model: Claude Opus 5.5
description: Brainstorm for the general UI (ui/): collection- and experiment-based navigation, gene/variant/locus pages across experiments, what is shared with ui-topchef, and what the store would need
---

# The general UI: brainstorm

Not a build plan. This lays out the design space for `ui/` once it stops being a copy of the TOPCHeF
app, so the calls can be made before any code is written. Each section ends in the questions it
leaves; the ones only Sam can answer are collected at the end.

## What is already fixed

- **D23** (2026-10-07): no default entry point. Nothing experiment-scoped loads until the visitor
  has chosen something.
- **D24**: a collection is provenance, not a compatibility claim. The reader checks catalogs itself
  before putting two experiments on one axis.
- **D25**: `kind` is closed (`qtl_cis | qtl_trans | gwas`); `phenotype_type` is an open label.
- **Sam, 2026-10-09**: `ui/` is not a mirror of `ui-topchef/`. What both need (the v2 reader, the
  error page and toasts, the fault page) lives in `ui/`; the UX is redone for collection- or
  experiment-based navigation, with gene, variant and locus pages not tied to one experiment.
- `ui-topchef/` stays the pinned single-study site at topchef.databio.org. The general UI leaves
  out what only applies to it (the published coloc table, the DCM wording, the pinned GWAS pairing,
  the coloc pin in `store-check`).
- **Sam, 2026-10-09: presence first.** Gene and variant pages do not fan out every experiment's
  results up front. They show where the gene or variant is present, and the visitor expands the
  experiments that have it.
- **Sam, 2026-10-09: scope for the v2 release.** No cross-study aggregation yet. Gene and variant
  pages stay experiment-focused. Coloc, LocusCompare and the GWAS genome track work against any
  GWAS experiment the visitor picks, not one hardcoded to DCM.
- **Sam, 2026-10-09: the GWAS genome track is universal.** The bin summary is computed by
  `gwas.build` for every GWAS from the contract `gwas.parquet` rows, rather than by the DCM adapter
  from its source TSV, so any GWAS experiment has a track and no adapter does anything for it.
  Measured on DCM's 569 bins: `min_p` and `n_gws` unchanged in every bin, `n_variants` lower in
  167 (9,348 rows the reference does not read); the lead's beta and allele become alt-relative and
  its rsID the catalog's. The bin object's digest changes, with SPEC section 10's "equals v0"
  wording and the store-check reference. Pairing a GWAS with a QTL experiment on another catalog joins the two windows by
  `(pos, ref, alt)`; the track itself needs no join.
- **Sam, 2026-10-09: collections are addressed by their digest.** No separate key: one inside the
  document would be in the digest anyway, so it buys nothing.

## What the store gives a reader today

| Question | How | Cost |
|---|---|---|
| Which experiments, collections, kinds exist | `store.json` `experiments` map: `{kind, collection}` per id | 1 fetch |
| A collection's title, citation, source | its content-addressed JSON document | 1 small fetch per collection, cached forever |
| An experiment's catalog, annotation, phenotype type, significance rule, dof | `experiments/<id>.json` | 1 fetch per experiment |
| Gene name or id to chromosome | annotation's `lookup` (`.qgl`): header plus one bucket | 2 ranges per annotation |
| A gene's phenotypes and results location in one experiment | that experiment's search index part for the chromosome | 1 object per experiment (60 to 200 kB each, genome-wide TOPCHeF) |
| A variant's significant and credible-set records in one experiment | that experiment's hits frame | 1 range per experiment |
| Whether a site is in catalog X, and its vidx there | the store's overlap index (SPEC 19): directory, one union page, one mask chunk | directory once per session (~170 kB at two catalogs), then 2 ranges per site, for all catalogs at once. No browser reader yet |
| rsID to site | each catalog's rsID index | 1 to 2 ranges per catalog |
| Nominal stats for a gene's window | the results block | 1 range per experiment per phenotype |

Everything below the first three rows is **per experiment**. A page that shows "everything the store
knows about gene X" fans out to every experiment that might have tested it. That fan-out is the
central performance problem of an experiment-independent page.

## What the store will hold within the next year

From `pipeline/config.yaml` and `plans/2026-10-06-pqtl-sources.md`:

| Cohort | Experiments | Phenotype types | Tissue | Catalog / assembly | Annotation |
|---|---|---|---|---|---|
| TOPCHeF | cis and trans, `ge` and `leafcutter` (4) | gene expression, intron | heart LV | `topchef_grch38` | GENCODE v34 |
| GTEx v8 heart LV | cis `ge`, cis `leafcutter` (2) | gene expression, intron | heart LV | `gtex_v8_heart_lv_grch38` | GENCODE v39 |
| DCM GWAS (Jurgens 2024) | 1 GWAS | disease | none | `topchef_grch38` | none |
| ARIC plasma pQTL | EA, AA (2+) | protein (SomaScan aptamer) | plasma | own catalogs | GENCODE v34 |
| MESA | eQTL GRCh37, pQTL GRCh38 | expression, protein | blood | GRCh37 and GRCh38 | ? |

So about a dozen experiments, three or more phenotype types, at least three annotations, five or
more catalogs, two ancestries, and possibly two assemblies. The design should hold at 50
experiments without rework, though a dozen is the honest near term.

## 1. Navigation model

Three ways in, not exclusive:

1. **Collection-first.** Browse publications; open one; inside it, the TOPCHeF-style app. Easy,
   because it is what exists, but it is single-study navigation with a menu in front. It does not
   give experiment-independent gene pages.
2. **Experiment-set-first.** The visitor assembles a set of experiments (a "workspace", from
   collections, filters, or picking scans), and every entity page shows that set. The set lives in
   the URL, so a link reproduces the view. Explicit and fast, since the fan-out is bounded by the
   set, but it asks a first-time visitor to understand experiments before seeing anything.
3. **Entity-first.** Search a gene, variant or region; the page shows every experiment with
   something to say, grouped by collection; facets (collection, phenotype type, tissue, kind,
   ancestry) narrow it. This is the GTEx portal / Open Targets / gnomAD shape. It satisfies D23
   naturally, since nothing loads until a search, but it has the full fan-out cost.

A likely shape is entity-first pages with facets, plus collection and experiment pages as browse
and landing pages that link into the entity pages pre-filtered (`/gene/FLNC?c=topchef`). The facet
state is the experiment set of model 2, made implicit.

**Questions:** which way in does a typical visitor take: a gene name, an rsID, a paper? Should
facets default to "everything", or to something narrower once the store is large?

## 2. Pages

**Home.** A search box (gene, rsID, `chr:pos`, region) and the store's collections as cards: title,
citation, member scans, tissue, n. No study loaded. Store-level counts come from `store.json` and
collection docs alone.

**Collection page** (`/collection/<key>`). The publication: citation, source, member experiments
with kind, phenotype type, significance rule, n, and counts (tested, significant). Links to each
experiment page and to entity pages filtered to this collection. The natural place for a study's
own extras (TOPCHeF's coloc table), if the general app carries them at all.

**Experiment page** (`/experiment/<id>`). One scan: its counts, its rule, dof, precision, a
genome-wide overview of lead signals (from the search index: the lead p per phenotype, as a
Manhattan of phenotypes), and its gene table (today's Genes page, per experiment).

**Gene page** (`/gene/<ensembl id>`). The heart of the redesign.

- Header: symbol, id, coordinates per annotation (they differ between GENCODE releases).
- **Presence**: the experiments where the gene is present (tested), grouped by collection, and
  nothing else loaded. What "present" can cheaply say is section 6's question.
- **Expanding an experiment** loads its own results: the phenotypes the gene has there (1
  gene-level `ge`, n introns, n aptamers), best p, lead variant, significant by that experiment's
  rule, credible sets. Only expanded experiments cost requests.
- **Detail**: selecting expanded experiments draws locus plots. Options: one experiment at a time (today's page);
  stacked plots on a shared x axis (same assembly only); a LocusCompare-style scatter between two
  experiments joined by site through the overlap index; a GWAS track from any GWAS whose catalog
  joins through the overlap index.
- Trans: experiments where the gene is a trans target, and its trans hits.

**Variant page** (`/variant/<rsid | site>`). Presence from the overlap index: which catalogs hold
the site, so which experiments tested it, in one read. Expanding an experiment reads its hits frame
(cis and trans records); a forest plot of the variant's effect on the same gene or
protein across experiments; optionally a PheWAS-style plot across every phenotype it touches.

**Locus or region page** (`/region/<chr:start-end>`). A track view: genes, then one row per
experiment showing its lead signals in the window (from hits or the search index), then GWAS tracks.
Full nominal stats per experiment in a window cost one block per phenotype per experiment, so the
region view would show leads by default and expand to full stats per row.

**Compare page** (later). Two experiments genome-wide: shared significant genes, sign concordance of
lead effects, replication rate. Expensive to compute in the browser; a precomputed object if wanted.

**Question:** which detail views matter first: stacked locus plots, the cross-experiment
LocusCompare, or the effect forest plot?

## 3. Identity across experiments

**Genes across annotations.** The key is the unversioned Ensembl gene id; symbols drift between
GENCODE releases, coordinates shift, and some genes exist in one release only. A gene page needs to
tell three states apart per experiment: *not in its annotation*, *in the annotation but not tested*,
*tested and not significant*. Today's UI conflates them.

**Proteins.** A SomaScan aptamer maps to one or more genes (complexes, paralogues), and a gene can
have several aptamers. So "the gene's protein QTLs" is a many-to-many join through the phenotype
table's `gene_id`. Whether ARIC's phenotypes carry a usable `gene_id` per aptamer, and what to do
with multi-gene aptamers, needs checking against its tables.

**Variants across catalogs.** Same assembly: the overlap index answers membership and vidx exactly.
A variant page has to tell *not in this catalog* from *in the catalog with no significant record*
(hits only hold significant and credible-set records; nominal stats need the block of a specific
phenotype). Different assembly (MESA's GRCh37 eQTL): a different seqcol collection, so no overlap.
Either liftover at build time (a new catalog, a pipeline decision), or the UI shows GRCh37
experiments only on their own pages.

**Questions:** is GRCh37 data worth supporting, or is it lifted at build time and never seen as
GRCh37 by the UI?

## 4. Comparing effects

- **Orientation**: every catalog is `ref_alt`, so betas from different experiments are on the same
  allele. Signs compare directly.
- **Scale**: magnitudes do not compare across pipelines (TensorQTL on TOPCHeF's normalisation
  against the eQTL Catalogue's, and protein levels in their own units). Showing z or signed
  -log10 p alongside beta, and labelling beta as experiment-scaled, avoids an implied equivalence.
- **Significance**: each experiment keeps its own rule (`p_perm < 0.05`, Bonferroni for ARIC), and
  every badge says whose rule. No store-wide threshold.
- **Colocalization across experiments**: `coloc-abf.ts` exists. Coloc between two QTL experiments
  (TOPCHeF and GTEx for one gene), or a QTL and any GWAS, needs both sides' nominal stats on
  matched sites (the overlap index again) and per-pair priors. A strong feature, and the most
  expensive one in the browser.

For the v2 release, coloc is computed live per gene or per locus (below); concordance and other
cross-experiment statistics wait.

### Coloc lists

Decided with Sam, 2026-10-09:

- **No coloc list is hardcoded in either UI.** `ui-topchef/src/lib/coloc.ts` (21 eQTL and 4 sQTL
  symbols with GENCODE v34 positions pinned) goes. The authors will not supply their table.
- **The pipeline never runs coloc, and no pair is mandatory.** No per-pair precomputation, no
  Worker, no genome-wide on-demand screen (it reads every candidate's block: 357 eGenes and 618
  introns near DCM's 26 loci, 14 MB and 23 MB).
- **Optional stored lists, one per declared pair.** `colocs/<id>.json`, a store-level pointer like
  `overlaps/`, naming the QTL experiment and the GWAS experiment by id with their catalog identity
  digests (a rebuilt experiment fails `validate` rather than leaving the list stale), a `source`
  saying who produced the list and how (reported by a publication, with its collection digest, or
  computed, with the method), the threshold, and a small content-addressed table: `phenotype_id`
  (the QTL experiment's own id; null when the source names only the gene) and `gene_id`
  (unversioned Ensembl). **Membership only**: no PP.H4, no lead variant, since the browser computes
  both when a locus is opened. The pipeline only validates and stores a supplied list.
- **The track**: the selected GWAS's bins and loci; with a list for the selected pair, its loci
  marked and labelled with the list's source; without one, the GWAS loci with the QTL experiment's
  significant phenotypes near them (from its search index, already loaded). Opening a locus
  computes coloc for its candidates against that one window.
- **TOPCHeF's list**: 21 eQTL rows with gene ids, and 4 sQTL rows at gene level (`phenotype_id`
  null; which 17 introns is not recoverable), source "reported in Murray et al. 2026, PP.H4 > 0.8,
  against `dcm_jurgens2024_biobanks`".

Mapped from the pinned positions: all 22 reported genes sit within 1 Mb of a DCM locus, in 16 of
the 26.

## 5. URL and state

- Entity routes are experiment-free: `/gene/:id`, `/variant/:id`, `/region/:loc`.
- Facets in the query string: `?c=<collection>` or `?e=<experiment,...>`, `?type=ge,leafcutter`.
  Shareable; empty means everything.
- **Collections by digest**: `/collection/<digest>`. No format change. A metadata edit gives the
  collection a new digest, so a link to the old one finds nothing in the current `store.json`; the
  page can say so and offer the collection list.
- Selected experiment rows on a gene page, locus window, active tab: in the URL too, so a link
  reproduces the view (today's TOPCHeF app already keeps the window in the URL).

## 6. Performance at N experiments

Presence first moves the cost to expansion, which the visitor controls. What remains is the cost of
presence itself:

- **Variants: already cheap.** A site's presence in a catalog is presence in every experiment
  on it, and the overlap index answers all catalogs with one union page and one mask chunk after
  its directory. The browser has no reader for it yet (Python's `overlap.find` and `membership`). Being present in a catalog
  means tested, not significant; the expansion says which.
- **Genes: not answerable cheaply today.** Whether an experiment tested a gene lives only in that
  experiment's search index part, so presence would itself fan out (a dozen index parts, 1 to 2 MB
  at a dozen experiments), which is the thing presence first is meant to avoid. The annotation is
  no proxy: GTEx's leafcutter scan names 14,124 genes of GENCODE v39's 61,487.
- **A gene presence index** would close the gap: the overlap index's idea applied to genes. A
  derived object keyed by unversioned Ensembl gene id (the union across the store's annotations),
  one mask bit per experiment, in hashed buckets like the gene lookup so a gene costs one small
  range. Optional extra columns per experiment (significant or not, phenotype count) would let
  presence show "significant in 3 of 5" without expansion. New derived object kind, `validate`
  checks, rebuilt when an experiment comes or goes, ignored by `ui-topchef`.

Reader side, regardless: the module singletons in `lib/store.ts` become caches keyed by experiment,
catalog and annotation id (`openExperiment(id)`, `openCatalog(id)`, `openAnnotation(id)`), each
fetched once and shared across pages. DuckDB holds the union of what a page materialises, tagged by
experiment.

**Question:** is a gene presence index acceptable as format work, and should it carry significance
or presence alone?

## 7. Metadata the facets need

Facets on tissue, ancestry, sample size, assay or platform are not in the experiment documents
today; the collection `source` is free-form. Options: structured fields on the experiment document
(`tissue`, `ancestry`, `n_samples`, `assay`), each a SPEC addition, or a free-form `labels` map the
UI facets on without the format knowing the vocabulary. Ontology terms (UBERON for tissue, EFO for
traits) would make them comparable with other resources.

**Question:** which facets matter: tissue, ancestry, phenotype type, assay, study? Ontology ids, or
plain labels?

## 8. Code shared with ui-topchef

Sam's direction: the shared parts live in `ui/`. Ways to do that:

| Option | How | Consequence |
|---|---|---|
| A. Path alias | `ui-topchef` imports `@core/*` from `../ui/src/core` via the Vite and TS config | One copy, no package machinery; ui-topchef's deploy must also trigger on `ui/src/core/**`; a change made for the general app can break TOPCHeF silently unless its checks run |
| B. Workspace package | `packages/qtl-core` (or `ui/core`) as an npm workspace both apps depend on | The cleanest boundary; one more `package.json`; both deploy workflows watch it |
| C. Keep copying | As now | Already drifting; the reason Sam asked |

What would be shared: `store-decode.ts`, a per-experiment store reader, `significance.ts`,
`rounding.ts`, `coloc-abf.ts`, `db.ts`, the error page, toasts and faults, theming, and probably
the plot primitives (`LocusPlot`, `GeneTrack`). TOPCHeF would then read its pinned collection
through the same per-experiment reader, which changes how it opens the store (see section 6).

**Question:** A or B, and does ui-topchef move onto the shared reader now or stay on its copy until
the general UI settles?

## 9. Prior art worth looking at before designing

GTEx portal (gene page with per-tissue rows, multi-tissue eQTL views), the eQTL Catalogue's own
browser, Open Targets Platform (entity-first with evidence rows), PheWeb and LocusZoom (association
plots, PheWAS), gnomAD (variant and gene pages, search-first), eQTLGen, QTLbase. Worth a session of
screenshots and notes on what each does with "many studies, one gene".

## 10. What this means for the release

The release plan (`plans/2026-10-09-v2-release.md`) assumed the general app's changes go in PR #5.
With the UX redo, that no longer fits one PR. Options for `ui/` in the release:

- deploy the v2 mirror as it stands (fixes the general Worker, which has been broken since
  2026-10-07; it shows TOPCHeF only, since it opens `topchef_ge_cis`);
- or leave `ui/` out of the release, and its Worker broken, until the redesign lands.

The GTEx rebuild, the upload and the TOPCHeF deploy do not depend on either.

## Questions for Sam

1. **Audience**: who is the general UI for: the lab, collaborators, the public? It decides how much
   the UI explains experiments, and whether search or browse comes first.
2. **Way in** (section 1): entity-first with facets, an explicit experiment set, or both?
3. **First cross-experiment views** (sections 2, 4): stacked locus plots, cross-experiment
   LocusCompare, effect forest plots, coloc between QTL experiments: which first?
4. **Format work** (sections 6, 7): a gene presence index (presence alone or with significance),
   structured facet metadata.
5. **GRCh37** (section 3): supported in the UI, or lifted at build time?
6. **Code sharing** (section 8): path alias or workspace package, and when TOPCHeF moves onto it.
7. **The release** (section 10): ship the v2 mirror to the general Worker now as a stopgap, or hold
   `ui/` out of PR #5?
8. **Prior art** (section 9): any of these you like or dislike, or others you have in mind?
