---
date: 2026-10-07
status: in-progress (D1-D26 settled; pipeline steps 1-3 built, F2 reversed; reader and apps not started)
model: Claude Opus 5 (1M context)
description: The experiment level is a cohort bundle holding other people's data — split it into one-run experiments, durable collections for provenance, and ephemeral tags for curation
---

# The experiment model

## The defect

`experiments/topchef.json` holds, in one document:

- two phenotype types (`ge`, `leafcutter`), each with a cis scope and a trans scope — four separate
  tensorQTL invocations;
- the **DCM GWAS of Jurgens et al. 2024**, a different publication's meta-analysis, as a field.

The second is not a naming quibble. TOPCHeF did not produce that GWAS; it colocalized against it.
The store records it as TOPCHeF's, and that has a concrete consequence:

```
remove_experiment(store, "topchef")   # unlinks the pointer
gc(store)                             # collects every object it named
```

The only thing referencing the 22 `.qbg` objects, the `.qgi` index and the bin summary is the
TOPCHeF pointer. **Removing a Virginia heart-failure cohort deletes a European DCM
meta-analysis.** The reverse costs too: a second cohort colocalizing against the same GWAS would
dedupe the objects by content address but copy the metadata — `source`, `orientation`,
cases/controls, title — into a second pointer, leaving two independent claims about one analysis
with nothing keeping them in agreement.

The same container shape makes four of the five files in the CVDKP zip unrepresentable. `dcm_gwas`
in `pipeline/config.yaml` picks one, and the full META build sits parked at
`data/derived/_full/gwas_meta/` because there is nowhere in the store to put a second GWAS of the
same trait — exactly the comparison a reader of the paper would want, given the figures were drawn
from BiobanksOnly and the Methods cite META (`plans/2026-09-04-r2-deploy.md` D9).

## The model

Two levels, replacing one. A third — tags — was intended and is deferred (see below), so nothing in
the store records "shown together" and the browser's deployment config carries it meanwhile.

### 1. Experiment — one analysis run

One experiment is **one modality, one scope, one run**: `ge` cis, `ge` trans, `leafcutter` cis,
`leafcutter` trans, `dcm_biobanks_only`. Standalone and self-describing: its own source and
citation, its own `dof`, its own `significance` rule, its own `precision` block, its own variant
catalog and annotation.

This is a flattening, not new machinery. `dof`, `significance`, `precision` and `counts` are already
per-`results[]`-entry today; they become fields of the document. A GWAS stops being a field on
another experiment and becomes an experiment, which removes the one-GWAS-per-cohort limit outright.

### 2. Collection — published together, durable

A dedicated grouping that records **provenance**: these runs are one publication's output. The 5
Jurgens files are one collection; TOPCHeF's four runs are another. A property of the data, so it
does not change once written, and it is what a citation attaches to.

`overlaps/` is the precedent — a pointer level whose documents name other pointers and carry their
own identity checks (SPEC section 19). `collections/` follows it, and because `store.json` lists
pointer levels, adding one is additive. **No object bytes change.**

Collections may name collections, so a release containing sub-analyses nests. `validate` must
therefore reject a cycle and a collection naming a missing member, the way `_check_overlap` already
rejects a stale `identity_digest`.

### 3. Tags — shown together — deferred, not designed

"TOPCHeF's QTLs beside the DCM GWAS it colocalized against" spans two publications and is a curation
decision someone made for the browser, not a fact about the data. The intent (Sam, 2026-10-07) was a
lightweight tag **at the experiment level**, carrying no identity digests and no validation
guarantees beyond "the things named exist" — deliberately not durable, and deliberately kept apart
from collections, since folding curation into provenance would recreate today's conflation one level
up.

**Out of scope for now** (Sam, 2026-10-07). Two consequences to keep straight:

- Tagging at the experiment level bounds the cardinality to tens, because there are only ever a
  handful of experiments. Any design needing a queryable table — parquet over DuckDB httpfs, a SQL
  service — is permanently out of proportion to the problem. If tags return they are a small
  document, and the only real question is whether it ships in the repo (reviewed, versioned,
  per-deployment) or sits beside the store in the bucket (changeable without a deploy).
- **Something must still express "shown together" in the interim**, or removing the `gwas` field
  leaves the browser unable to pair them. **Deployment config already does this.** `EXPERIMENT_ID`
  is `import.meta.env.VITE_EXPERIMENT` (`ui/src/lib/store.ts:33`) and `VITE_DATA_BASE` lives in
  `ui/.env.production`, so which experiment an app shows is already a per-deployment build-time
  decision, with two deployments from one repo (`qtl-browser`, `qtl-browser-topchef`). A second
  variable naming the GWAS experiment — or one naming a collection — puts the pairing exactly where
  the other per-app curation already lives. No new mechanism, no store mutation.

## Decisions & ownership

| # | Decision | Ownership | Why |
|---|---|---|---|
| D1 | An experiment is one modality and one scope: one instance of a cis scan, a trans scan, or a GWAS | **user-owned, surfaced** — answered 2026-10-07 | Sam: "each experiment gets its own single modality attempt... the 5 sumstats is 5 separate gwas experiments", then "cis trans separate. a qtl can be cis or trans, but be one instance of a cis or trans". So TOPCHeF is 4 experiments in 1 collection, and the DCM GWAS is a 5th in its own |
| D2 | "Published together" and "shown together" are **separate types**, not two flavours of one grouping | **user-owned, surfaced** — answered 2026-10-07 | Sam: provenance "should get dedicated collection label or grouping", curation "much more flexible and more ephemeral like tagging two recommended results together". Typed separately so curation cannot silently become provenance. Only the provenance half is being built (D7), and it is a field rather than a level (D16) |
| D3 | The DCM GWAS becomes its own experiment, not a field on TOPCHeF | **user-owned, surfaced** — answered 2026-10-07 | Sam: "the dcm gwas wasn't even part of topchef, so putting it under topchef experiment here doesn't make much sense" |
| D4 | ~~Collections nest as a pointer level with a cycle check~~ → **superseded**. Nesting survives only if the collection field is a digest (B) | AI-owned, defended | With no entity there is no member list to corrupt. A `parent` digest in the collection object still needs a terminating check, but it is one chain, not a graph |
| D5 | Settled before anything is uploaded under v2 | **user-owned, surfaced** | v2 already forces a full re-upload and nothing is published under it, so pointer restructuring is free now and a second migration later |
| D6 | ~~Object bytes do not change; pointer-level JSON only~~ → **no longer true**, superseded by D18 and D19 | AI-owned, defended | It held while the format was assumed fixed. D18 opens the format, so F1 and F2 change object bytes; D19 changes `store.json`'s shape. The *weaker* claim still holds and is what D13 relies on: no object is re-encoded **except** those a chosen F-item touches |
| D7 | Tags are **out of scope**; "shown together" is deployment config for now | **user-owned, surfaced** — answered 2026-10-07 | Sam: "we can just skip tags for now". The relation does not disappear, it moves to the place that already carries per-app curation (`VITE_EXPERIMENT`), which means the rework can proceed on two levels instead of three. Revisit when a second pairing is wanted that config cannot express |
| D8 | ~~The grouping level is called `collections/`~~ → **superseded by D16**: there is no level, only a field named `collection` | **user-owned, surfaced** — answered 2026-10-07 | The word is still "collection"; it is no longer a pointer level |
| D16 | **A collection is a field on an experiment, not an entity.** Membership is derived by matching | **user-owned, surfaced** — answered 2026-10-07 | Sam: "just a field an experiment can have and if they either hash match or string match between experiments they are in the same collection". Removes the member list, the two-way agreement check and the cycle guard; matches how `catalog` and `annotation` are already referenced; and `source.zenodo` is already such a key |
| D17 | The field is a **digest** of a small immutable metadata object (B) | **user-owned, surfaced** — answered 2026-10-07 | Sam: "digest then". Gives the citation one home so sibling experiments cannot disagree, makes digest equality the validation, nests via a `parent` digest, and caches forever where a pointer document would be `no-cache` |
| D19 | Collections are tabulated from `store.json`, whose `experiments` array becomes a map carrying `kind` and `collection` | AI-owned, defended | Zero extra requests: it is the first fetch on every page load, 195 bytes today, and already generated by `write_store` from the pointer directory — a derived index, not the hand-maintained duplicate D16 removed. A range-readable or key-value object in the bucket is the right shape only at thousands of experiments, and the format already has that pattern in the kind 9 variant index |
| D18 | A format change is in scope alongside the pointer change | **user-owned, surfaced** — answered 2026-10-07 | Sam: "a format change can be made as well". Opens two things the pointer rework cannot fix on its own — see below |
| D9 | **A cis scan and a trans scan are separate experiments.** TOPCHeF is 4, plus the GWAS as a 5th | **user-owned, surfaced** — answered 2026-10-07 | Sam: "cis trans separate. a qtl can be cis or trans, but be one instance of a cis or trans". An earlier draft of this plan assumed the opposite on the grounds that one document with two scopes was the smaller change. It was the wrong instinct: the split **removes** machinery, because the search index and hits currently serve both scopes and carry `search_index_trans_only` solely to hold the residue |
| D10 | **Each experiment owns its own search index and hits**, with its own ordinal space | **user-owned, surfaced** — answered 2026-10-07 | Sam: "a". Decided on build coupling rather than read cost: `ord` is a row number in the table, so one shared table means one shared numbering, which cannot be allocated without building all four experiments together in a fixed order — rebuilding one would renumber the others and rewrite all four pointers, leaving documents only nominally separate. Bytes are about neutral (chr1's 323 KB part splits by rows and columns) and the four reads are parallel. The cost accepted: a variant that is a cis lead for one gene and a trans signal for another is two records in two objects for the reader to join |
| D11 | Experiment id convention | AI-owned, default | `topchef_ge`, `dcm_jurgens2024_biobanks` — flat, lowercase, underscore. It is URL-visible and it is the value of `VITE_EXPERIMENT`, so it is hard to change later |
| D12 | An experiment carries a `collection` field | **folded into D16** | It is no longer a *back*-pointer, because there is no forward one. One direction, so nothing can disagree |
| D13 | ~~Migration is a converter~~ → **Migration is a rebuild** | **user-owned, surfaced** — answered 2026-10-09 | A converter existed to migrate the live store without the contract tables. The tables are local and a genome-wide rebuild is nine minutes, so the migration is `build_store` plus `upload`; content addressing still reuses every unchanged object in the bucket |
| D14 | Coloc's priors and `sdY` live in **deployment config** | **user-owned, surfaced** — answered 2026-10-07 | Sam: "coloc priors and sdy deployment conflig". They are a per-app analysis choice, not a property of either experiment, so they sit beside the other build-time values rather than in the store. v2 gap G8 closes this way: the assumption stays documented in `coloc-abf.ts`, but the number stops being hardcoded |
| D15 | **A GWAS experiment owns no `hits` and no search index** | **settled by F1** | `hits` was the only place a GWAS credible set fitted, and a hit record is keyed by `(vidx, ord)` where `ord` indexes a search index — so a one-phenotype experiment would have carried a degenerate ordinal space purely to host its credible sets. With F1 the GWAS block carries them, and the question dissolves |
| D20 | **F1 only.** The trans frame and the GWAS block stay separate file types; each gains a credible-set section | **user-owned, surfaced** — "f1 and f2 combined" 2026-10-07, **reversed** "reverse the f2" 2026-10-07 | The merge was taken on advice that it closed G4 for free, which was false (no QTL source has a per-row n; D-cis-n below), and that it fixed G12, which is real but does not need a merge. Built and measured, it cost more than it bought: cutting a phenotype's rows at every chromosome turned 100,226 frames into a projected 1.7 M blocks of ~6 rows, because a trans phenotype's hits lie on a median of 16 chromosomes. A two-chromosome subset could not show this, since the subset itself capped every phenotype at two. Reversed before anything was uploaded, at the cost of one local rebuild. See "F2" below for the numbers |
| D26 | **No per-row n on the cis row** | **user-owned, surfaced** — answered 2026-10-07 | Sam: "no n for cis then". No QTL source in hand carries one — tensorQTL writes none, ARIC's PLINK2 `OBS_CT` is a single constant per cohort — and the one reader that would consume it, coloc's `sdY.est`, is one G8 decided must not be used. An unconditional byte on 623 M rows is 623 MB of copies of `dof`. If a source with a varying n arrives, the path is a column gated by the cis block's `flags` bit, appended after `se` so the `nlp`-first prefix invariant holds, with the n table in the experiment document |
| D21 | Collections **nest**, via a `parent` digest | **user-owned, surfaced** — answered 2026-10-07 | Sam: "digest collections should support nesting". `validate` checks the chain resolves. A cycle is impossible among well-formed objects, since each name depends on the other's bytes, so the loop and depth guards exist only to terminate on a store whose bytes no longer match their names — which validate's digest check reports on its own |
| D22 | **Folded into v2.** This work ships in the v2 PR, under `format_version` 2 | **user-owned, surfaced** — answered 2026-10-07 | Sam: "fold into v2. this is the v2 pr". v2 is built and unpublished, so one version and one upload. The consequence to hold onto: `format_version` 2 now means the pointer shape *and* the codec, so the stores already built locally under v2 are not v2 as finally defined and will be rebuilt |
| D23 | **No default entry point.** The general UI loads no experiment or collection on first paint | **user-owned, surfaced** — answered 2026-10-07 | Sam: "nothing should be a default entry point. no particular experiment or collection loaded on initial load for general ui". So startup is `store.json` alone — the map is enough to render a chooser — and every experiment-scoped read waits for a selection. This deletes the speculative-pointer machinery in `openStore` (`ui/src/lib/store.ts:230-250`), which exists to guess a catalog and annotation id before the experiment document arrives — there is nothing to guess when nothing loads. **The TOPCHeF app stays a separate single-study app** and pins its own collection; that is what distinguishes it from the general one, where today the two differ only by `VITE_EXPERIMENT` and branding |
| D24 | **Collection members need not agree on anything** — no `validate` check | **user-owned, surfaced** — answered 2026-10-07 | Sam: "collections do not have to agree on anything". A collection is provenance, not a compatibility claim, and MESA is the case that proves it: one publication whose eQTL release is GRCh37 and whose pQTL release would be GRCh38. The burden moves to the reader — a UI must never put two experiments on one axis without checking their catalogs itself |
| D25 | `kind` is a **closed, versioned** enum | **user-owned, surfaced** — answered 2026-10-07 | Sam: "closed and versioned enum then". The rule this follows, worth stating because it decides where a future axis goes: **open where the value is a label, closed where the value selects a layout.** `phenotype_type` is open and always was — SPEC section 8, "data, not a file kind" — which is why `protein_somalogic` cost no format change. `kind` picks which object kind, which fields exist and whether there is an annotation, so a reader cannot shrug at one it does not know |

## What this changes elsewhere

- **Search index and hits objects currently span phenotype types.** One `.arrow.zst` index and one
  `.qbh` per chromosome cover `ge` and `leafcutter` together. Split per experiment and a gene page
  goes from one index read to N, and from one hits read to N. Measured genome-wide on chr1 after
  the split: the cis parts are 60 kB (`ge`) and 199 kB (`leafcutter`) against 322 kB before, since
  each is narrower. The trans index is split the same way, by the phenotype's chromosome; an
  earlier build left it whole on the mistaken ground that a trans phenotype has no chromosome, and
  the gene page's trans lookup was then a 307 kB plus 1,706 kB read.
- **Coloc gains a home.** coloc.abf runs in the browser, so nothing stored is in the wrong place
  today. But the posteriors are a property of the *pair* and belong to neither experiment — under
  this model a collection or tag is where its priors and the `sdY = 1` assumption would sit, rather
  than a hardcoded constant (v2 plan gap G8).
- **`remove_experiment` and `gc` become safe.** Removing a cohort stops being able to delete another
  publication's data, because that data is no longer reachable only through the cohort's pointer.
- **The parked META build becomes storable**, so BiobanksOnly and META can sit side by side and the
  open question from D9 of `plans/2026-09-04-r2-deploy.md` — which set the coloc itself used — can
  be answered by comparison rather than correspondence.
- **SPEC sections 8, 10 and 14 are rewritten**, and section 8's `precision` example is already stale
  against v2 (it still shows `slope_max_error_over_se`).
- **Scale**: 160 references to `experiment` in `pipeline/`, 95 in `ui/src`. The UI's gene page
  assembles one view from several runs today by reading one experiment; it would assemble it from
  several experiments named by a collection.

## Proposed shape

Pointer levels are **unchanged** — no `collections/` level, because a collection is a field (below):

```
store.json                     name, format_version, refget, the ids at each level
experiments/<id>.json          one analysis run; carries a `collection` field
variant_catalogs/<id>.json     unchanged
annotations/<id>.json          unchanged
overlaps/<id>.json             unchanged
```

### One experiment is one scan: cis, trans, or GWAS

**A QTL experiment is a cis scan or a trans scan, never both** (Sam, 2026-10-07). So TOPCHeF is four
experiments and the DCM GWAS is a fifth:

```
collection <digest-A>.json   (TOPCHeF)      <- topchef_ge_cis, topchef_ge_trans,
                                               topchef_leafcutter_cis, topchef_leafcutter_trans
collection <digest-B>.json   (Jurgens 2024) <- dcm_jurgens2024_biobanks
```

Three kinds, one per object kind, so a document's shape is known from its `kind` alone:

| `kind` | object | phenotypes | window |
|---|---|---|---|
| `qtl_cis` | kind 2 `.qbe` | many | a window per phenotype |
| `qtl_trans` | kind 6 `.qbt` | many | genome-wide |
| `gwas` | kind 4 `.qbg` | one | genome-wide |

Written out, `qtl_trans` and `gwas` are visibly adjacent — both genome-wide, differing only in
phenotype count. Merging their file types was tried as F2 and reversed; the two are adjacent in
the row and opposite in shape, and the shape is what the framing has to fit.

**The split removes machinery rather than adding it**, which is the opposite of what an earlier
draft of this plan assumed when it defaulted to two scopes in one document:

- `INDEX_SCHEMA` is 21 columns because it serves both scopes: 8 are cis-only (`blk_off`, `blk_len`,
  `var_start`, `n_var`, `var_off`, `var_len`, `w_lo`, `w_hi`) and 3 are trans-only (`trans_off`,
  `trans_len`, `n_trans`). Split, each index carries only its own.
- **`search_index_trans_only` disappears**, and a smaller residue takes its place. It existed
  because a phenotype with trans rows and no cis result had a null `chr` and could not be filed
  under a chromosome part. In a trans experiment the index is split by the phenotype's chromosome
  like the cis one — its cis placement, else its gene's chromosome from the annotation — and only
  a phenotype with neither lands in `search_index_unplaced`: genome-wide, 17 of 19,440 `ge` and 36
  of 80,786 `leafcutter` phenotypes before the annotation fill.
- Hits stop mixing: `HIT_LEAD` and `HIT_CS` come from the cis path (`results.py:473,478`),
  `HIT_TRANS` from the trans path (`:838`), and today one `.qbh` holds all three. Split, a cis
  experiment's hits hold leads and credible sets, a trans experiment's hold trans rows.

The cost is read fan-out, which is D10.

### Shapes

What actually differs between a QTL run and a GWAS run is narrower than it first looks, and getting
the list right matters because it is the list the discriminator exists to separate:

| | QTL run | GWAS run |
|---|---|---|
| phenotypes | many (19,423 genes) | **one** — the trait |
| window | a cis window per phenotype, optionally genome-wide trans | genome-wide |
| annotation | needed, phenotypes are genes | not needed |
| `dof` / `N` | `dof`, provenance only in v2 | `N` per row, in the n table |
| credible sets | yes | **yes** — fine-mapped GWAS loci are routine |

**A GWAS can have credible sets** (Sam, 2026-10-07). An earlier draft of this plan said it could
not, and used that as part of the argument for a `kind` discriminator. It is wrong: SuSiE or FINEMAP
over GWAS summary statistics is standard practice, and a fine-mapped GWAS locus is precisely the
input colocalization wants. v2 cannot currently store them for a GWAS — only kind 2 `.qbe` has a
credible-set section — which is now recorded as gap **G11** in the v2 plan, and is a layout change,
so free before upload and a re-encode after.

With that corrected, a GWAS reads as **a single-phenotype genome-wide scan**: structurally a trans
scan with one phenotype rather than a different sort of thing. Acting on that by merging the two
file types was tried and reversed (F2 below). The discriminator stays, for the two differences that
are real: a GWAS has no annotation, and it has one phenotype rather than an index of them. An
explicit `kind` lets a reader reject a document it does not understand instead of probing for
absent fields.

The fields every kind shares, then the ones its kind adds. `catalog`, `catalog_identity` and
`allele_orientation_source` are common to all three; `annotation`, `phenotype_type`, `dof` and
`significance` to the two QTL kinds.

```jsonc
// experiments/topchef_ge_cis.json
{"id": "topchef_ge_cis", "kind": "qtl_cis",
 "collection": "<digest-A>.json",       // the object name of the collection document (D17)
 "phenotype_type": "ge",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 "annotation": "gencode_v34", "annotation_source": null,
 "allele_orientation_source": "topchef_refcheck_a2_is_ref",
 "dof": 435, "dof_fit": null, "dof_source": "ingestion.json",
 "significance": {"column": "p_perm", "op": "<", "threshold": 0.05, "label": "Permutation p"},
 "counts": {"phenotypes": 19423, "with_rows": 19423, "significant": 10220, "significant_genes": 10220},
 "unplaced": {"count": 0, "examples": []},
 "files": {"chr1": "<digest>.qbe"},            // one per chromosome built
 "hits": {"chr1": "<digest>.qbh"},             // leads and credible sets only
 "search_index": "<digest>.arrow.zst",
 "search_index_parts": {"chr1": {"file": "<digest>.arrow.zst", "rows": 986, "ords": [[0, 189]]}},
 "precision": {"theoretical_limit": {}, "measured_worst": {}},
 "source": {}}

// experiments/topchef_ge_trans.json
{"id": "topchef_ge_trans", "kind": "qtl_trans",
 "collection": "<digest-A>.json",
 "phenotype_type": "ge",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 "annotation": "gencode_v34",
 "allele_orientation_source": "topchef_refcheck_a2_is_ref",
 "dof": 435, "dof_source": "ingestion.json",
 "file": "<digest>.qbt",                        // one object, chromosome `all`, one frame per phenotype
 "hits": {"chr1": "<digest>.qbh"},              // trans records only
 "search_index": "<digest>.arrow.zst",
 "search_index_parts": {"chr1": {"file": "<digest>.arrow.zst", "rows": 1959, "ords": [[0, 1958]]}},
 "search_index_unplaced": {"file": "<digest>.arrow.zst", "rows": 17, "ords": [[19423, 19439]]},
 "n_rows": 2642166, "n_phenotypes": 19440,
 "excluded": {"variants": 32765, "rows": 37951},
 "rows_in_source": 2680117,
 "precision": {"theoretical_limit": {}, "measured_worst": {}},
 "source": {}}

// experiments/dcm_jurgens2024_biobanks.json      kind: "gwas"
{"id": "dcm_jurgens2024_biobanks", "kind": "gwas",
 "collection": "<digest-B>.json",
 "title": "Dilated cardiomyopathy GWAS, Jurgens et al. 2024 (biobanks-only meta-analysis)",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 // credible sets, when a source has them, live in the .qbg blocks themselves (F1), so a GWAS
 // carries no hits object and no search index
 "files": {"chr1": "<digest>.qbg"}, "index": "<digest>.qgi",
 "n_rows": 12494731, "rows_by_chrom": {}, "rows_sharing_a_site": 1593062,
 "block_rows": 2048, "n_values": [],
 "bins": {"file": "<digest>.arrow.zst", "bin_bp": 5000000, "n_bins": 569},
 "precision": {"theoretical_limit": {}, "measured_worst": {}},
 "orientation": {}, "source": {}}
```

What collapses relative to today: `results[]` disappears, because the experiment *is* the result.
`dof`, `significance`, `precision`, `counts` and `unplaced` move up from the list entry to the
document. The top-level `trans` row counts and `trans_excluded` merge into the one `trans` block.
The `gwas` field is gone.

### A collection is a field, not an entity

**Not a pointer level** (Sam, 2026-10-07): "collections should probably not even be an entity. i
think just a field an experiment can have and if they either hash match or string match between
experiments they are in the same collection."

Membership is then *derived* rather than declared, and that removes three things the entity version
needed: a member list to keep in agreement with the back-pointers (so D12's two-way check goes), a
cycle guard (D4), and a document to edit when an experiment joins. It also matches how the store
already treats `catalog` and `annotation` — an experiment names them by id, two experiments naming
the same id share them, and no catalog carries a list of its experiments.

It is also already half-present: `source.zenodo` is `"zenodo_21382723"` for TOPCHeF, which is
exactly a release key that all four of its experiments would match on.

**The field holds a digest** (D17), not a name, because that is what decides whether the citation has
one home:

| | grouping | metadata home | nesting |
|---|---|---|---|
| a name, `"collection": "topchef_2026"` | string equality | none — title and DOI repeat per experiment, with nothing keeping the copies in agreement | no |
| **a digest**, `"collection": "<sha512t24u>"` of a small immutable object | digest equality, which **proves** identical metadata | the object | yes, via a `parent` digest |
| no new field, group by `source.zenodo` | existing data, zero additions | none | no |

The digest keeps the one property the entity version existed for — a single home for the citation, so
four sibling experiments cannot disagree about their publication — while giving up the entity.
Equality of digests *is* the validation, so nothing new has to be checked, and a content-addressed
object caches forever where a pointer document is served `no-cache` with no validator
(`ui/src/lib/store.ts:230`).

Grouping by `source.zenodo` would cost nothing at all, but it breaks on a release with no Zenodo id,
on two publications in one record, and on anything not from Zenodo — ARIC is a supplement download,
MESA is Dropbox.

### What a collection object holds

One immutable object per publication, named by its own digest (D17). **No member list** — the members
are whichever experiments carry the digest:

```jsonc
// immutable/<digest>.json, referenced as "collection": "<digest>" by each member
{"title": "Jurgens et al. 2024, dilated cardiomyopathy GWAS meta-analyses",
 "citation": {"doi": "10.1038/s41588-024-01975-5", "year": 2024},
 "source": {"release": "CVD Knowledge Portal Jurgens2024_DCM_EU", "fetched": "2026-09-04"},
 "parent": null}                                // a digest, for nesting
```

What `validate` checks shrinks to almost nothing: the object exists, and if `parent` is set it
resolves and the chain terminates. There is no agreement to verify, because there is only one copy of
the metadata and one direction of reference.

### Tabulating collections

Membership is derived, so something has to answer "what collections exist, and who is in each"
without fetching every experiment document — those are 14 KB each today and served `no-cache`.

**It goes in `store.json`, which costs nothing.** That file is already the first fetch on every page
load (`ui/src/lib/store.ts:244`), is 195 bytes, and already enumerates experiment ids. It is also
**generated** — `write_store` globs the pointer directory — so it is a derived index in the same
category as the search index, not a hand-maintained duplicate of the kind D16 removed. Turn the
experiment array into a map and carry the two fields a reader needs for routing:

```jsonc
{"name": "topchef", "format_version": 2, "refget": [],
 "variant_catalogs": ["topchef_grch38"],
 "annotations": ["gencode_v34"],
 "experiments": {
   "topchef_ge_cis":            {"kind": "qtl_cis",   "collection": "<digest-A>"},
   "topchef_ge_trans":          {"kind": "qtl_trans", "collection": "<digest-A>"},
   "topchef_leafcutter_cis":    {"kind": "qtl_cis",   "collection": "<digest-A>"},
   "topchef_leafcutter_trans":  {"kind": "qtl_trans", "collection": "<digest-A>"},
   "dcm_jurgens2024_biobanks":  {"kind": "gwas",      "collection": "<digest-B>"}},
 "overlaps": []}
```

Grouping is then a `GROUP BY` over an object already in hand: **zero extra requests**. A collections
page costs `1 + K` fetches where K is the number of distinct digests — and those K are
content-addressed, so they cache forever, unlike the pointer documents.

Carrying `kind` as well is what lets a reader route without opening anything: "which experiment here
is the GWAS" is the question the deployment config has to answer now that tags are deferred (D7), and
this answers it from the file the app already fetches first.

`validate` gains one check — each entry agrees with its experiment document's own `collection` and
`kind`. It cannot drift in practice because `write_store` generates it, and the check costs nothing
because validate already reads every pointer.

**Why not a range-readable or key-value structure in the bucket.** It is the right shape at a
different scale. At tens of experiments, one 1 KB JSON blob beats any index: no offsets to maintain,
no second request, and the whole thing fits in the fetch already being made. A range-readable object
only earns its complexity when the list stops fitting — thousands of experiments, a GWAS-Catalog-size
store — and the format already has that pattern if it ever comes to that: the variant index (kind 9
`.qbx`) is exactly a paged, offset-indexed directory, and a `.qbx`-shaped experiment index would
follow it. Until then it is the tags-cardinality mistake again: machinery out of proportion to tens
of rows.

### The one hard coupling: the ordinal space

This is where the split is not free, and it is the part to get right rather than to decide quickly.

A hit record is `(vidx, kind, ord, cs_id, value, beta, flags)`, and **`ord` is the phenotype's
ordinal in the experiment's search index**. That index spans phenotype types today:
`search_index_parts.chr21.ords` is `[[18097, 18286], [94836, 95631]]` — two ranges, because `ge` and
`leafcutter` occupy different stretches of one experiment-wide ordering. So the search index, the
hits objects and the ordinal space are one unit, and D9 splits that unit four ways.

#### The question in one paragraph

TOPCHeF becomes four experiments. Each needs two lookup tables: **"where are gene X's numbers?"**
(the search index) and **"which phenotypes care about variant V?"** (hits). Today there is one of
each, covering all of TOPCHeF. After the split: do we build **four of each**, one per experiment, or
**one of each** that all four reference?

#### What the two objects are, and who reads them

**The search index** is the experiment's table of contents: one row per phenotype, carrying where its
block is (`blk_off`, `blk_len`), which variants it covers (`var_start`, `n_var`, `var_off`,
`var_len`), its window (`w_lo`, `w_hi`), its significance, and its trans frame (`trans_off`,
`trans_len`). It is split per chromosome so a gene page fetches only its own. Measured on the
genome-wide v2 store: the whole index is 2.90 MB, a chromosome part is **143 KB median and 324 KB at
chr1** (10,182 rows), and all 23 parts sum to 3.31 MB.

**Hits** answer the reverse question — given a *variant*, which phenotypes say something about it.
Records are `(vidx, ord, value, beta, kind, cs_id)`, sorted by variant, framed 1,024 variants at a
time so a variant page range-reads one frame rather than the object. The objects are large —
**5.4 MB median per chromosome, 9.5 MB at chr1, 119 MB in total** — but a read is an offset table
plus one frame, a few kB.

**The ordinal space** is what ties them together. `ord` in a hit record is a row number in the search
index, so the two objects are one data structure in two files, and `ord` is only meaningful within
the experiment that minted it.

#### What the split actually costs

Not bytes. The rows are the same rows either way, and each one gets *narrower* — a cis index drops
the three trans columns, a trans index drops the eight cis ones — so the total may fall. What
multiplies is **requests**: a gene page goes from one index part to four, and a variant page from one
hits frame to four.

And those four are **parallel, not sequential**, because D19 puts every experiment id and kind in
`store.json`, which is already the first fetch. That matters more than the count: the existing
comment in `openStore` (`ui/src/lib/store.ts:230`) is about *sequential waves* — "the reader needs
two sequential waves of pointer fetches before it can request a single object" — and that is the cost
worth avoiding. Four parallel range reads in one wave is close to free on HTTP/2 against a CDN.

#### The two options, concretely

Measured on chr1, whose index part is 323,527 B over 10,182 rows — 1,959 of them `ge` (19%) and
8,223 `leafcutter` (81%):

**(a) Four tables, one per experiment.** chr1's one part becomes four, split by rows and by columns:
a cis index needs about 13 of the 21 columns, a trans index about 8. Same rows, narrower, so **the
total is roughly the 323 KB it is now** — but it arrives in four parallel requests rather than one.
Each experiment builds independently.

**(b) One table, referenced by all four.** Each experiment's pointer carries the *same* index digest,
so no owner is needed and nothing has to be an entity — content addressing makes four references to
one object provable rather than duplicative. A gene page reads one table, exactly as today.

**The objection to (b) is not read count, it is build coupling.** `ord` is a row number in a shared
table, so a shared table means a shared ordinal space, and a shared ordinal space cannot be allocated
without building all four experiments together and in a fixed order. Rebuilding one — a corrected
`dof`, a re-run trans scan — renumbers the others and rewrites all four pointers. That is precisely
the coupling D9 and D16 set out to remove: it would leave four documents that are only nominally
separate, since none of them can be produced alone.

(b) does keep one thing (a) gives up. A whole ordinal space lets a **single** hits record say "this
variant is the lead for gene X's cis scan *and* carries a trans signal to gene Y". Under (a) that is
two records in two objects for the reader to join. Whether that join is awkward is a UI question, and
D23 lightens it: with no default entry point nothing loads until a user has chosen, so the fan-out
follows a deliberate selection rather than landing on first paint.

**Recommendation: (a).** Same bytes, four parallel requests, independent builds. Step 2 is to confirm
the parallel-fetch assumption against the bucket rather than reason about it.

## What a format change opens (D18)

The pointer rework alone cannot reach these. Both were deferred earlier only because the format was
assumed fixed, and both are free now and a re-encode after upload.

**The kinds, since the numbers are shorthand** (SPEC section 3): every object's header carries a
`kind` byte naming its layout. The three that hold statistics are **kind 2 `.qbe`**, cis results —
one block per phenotype; **kind 6 `.qbt`**, trans results — one zstd frame per phenotype; and
**kind 4 `.qbg`**, GWAS — fixed blocks of 2,048 rows plus a separate `.qgi` index (kind 5) that turns
a position window into byte ranges.

### F1 — a credible-set section for the genome-wide kinds

**Recommended.** Only kind 2 can hold credible sets today: its block header has an `n_cs` count and
its payload a section of 12-byte `(row, pip, cs_id)` records. Kind 6's columns end at `allele` and
kind 4's payload ends at the allele heap, so a **fine-mapped GWAS has nowhere to put its credible
sets** — gap G11 in the v2 plan, and the gap sits on the feature the DCM GWAS exists to serve, since
a fine-mapped locus is what colocalization wants.

The fix mirrors kind 2: an `n_cs` field and a record section before the allele heap, in both the
trans frame and the GWAS block, since trans fine-mapping has the same gap. The trans frame already
has a reserved `u32` at payload offset 44 that becomes `n_cs`; the GWAS block's payload header grows
for it. Small, and it dissolves D15 — with its own credible sets a GWAS needs no hits object and no
degenerate one-row ordinal space.

### F2 — merge the GWAS file type into the trans file type: built, measured, reversed

The idea: a GWAS is structurally a trans scan with one phenotype, so the two file types could be
one. What actually differs between them:

| | trans, `.qbt` | GWAS, `.qbg` |
|---|---|---|
| the `u8` column | `ordinal`, the variant's chromosome | `n_code`, an index into the N table |
| framing | one frame per phenotype | fixed runs of 2,048 rows |
| windowed access | none — a frame is read whole | the `.qgi` index gives byte ranges per locus |
| phenotypes | many | one |

The merge was adopted (D20) on two arguments: that it fixed G12, and that freeing byte 16 for
`n_code` closed G4 on the trans side. The second was false — no QTL source carries a per-row n, so
the byte would hold zeros (D26). The first is real, but is a property of the framing and not of the
kind: a trans frame can be subdivided and indexed without becoming a GWAS block.

It was built for the trans side and measured, and the cut rule is where it failed. To free byte 16
the chromosome had to become a block property, so a block could never span a chromosome. That rule
suits a GWAS, whose rows are dense along one chromosome. Trans rows are the opposite shape:

| genome-wide v2, 4,010 phenotypes sampled | rows per phenotype, median | chromosomes per phenotype, median | blocks under the per-chromosome cut | frames today |
|---|---:|---:|---:|---:|
| `ge` | 98 | 16 | 318,000 | 19,440 |
| `leafcutter` | 105 | 17 | 1,358,000 | 80,786 |

About 1.7 million blocks of six rows each, every one paying a 56-byte header and zstd's per-frame
cost, roughly 114 MB of overhead on a 213 MB object before counting how badly six-row blocks
compress, plus a 22 MB index the reader would hold to window anything. The chr21/22 build showed
+17% and could not show more, because a two-chromosome subset caps every phenotype at two
chromosomes; on that store the median block was already 5 rows and headers were 39% of the object.

So the scan block and its index are deleted, the trans frame keeps the chromosome per row, and the
two file types stay two. What survives of F2's motive is G12's path, which costs nothing now: a
frame starts with an absolute position, so a builder can later split a large phenotype into several
frames inside one span and add an index object beside it, with no change to the frame layout. The
one obligation on readers is to inflate a span frame by frame rather than assume it is one frame.

Two things the attempt also found, and which stand regardless: the per-block ruler does not
compress (zstd takes 8–12% off a quantized column, where v1's fixed-point GWAS integers lost far
more), and estimates made on a subset build can be wrong in kind rather than in degree when the
subset changes the shape of the data.

## Relation to qtlstore v2

`plans/2026-10-07-qtlstore-v2.md` steps 1–4 are done and committed: the codec, the version
threading, the build driver, the code sweep, and a verified genome-wide build. **Its steps 5–7 —
the reader, coloc, both apps, the gap re-check, upload — are on hold pending this plan** (Sam,
2026-10-07), because the restructure changes what a reader opens and how many objects a gene page
assembles. Writing the v2 reader first would mean writing it twice.

Two v2 questions therefore land in this plan's territory and must not be lost:

- **The credible-set section's position in a cis block.** v2's columnar layout makes a one-request
  partial read possible (12,802 B against 38,632 B, and below v1's 25,892 B), but the credible-set
  records sit past all three code columns while the locus plot needs them. Moving them ahead of the
  code arrays would make the plot's actual read a single prefix. It is a byte-layout change: free
  before upload, a re-encode after. Full numbers and the decision criterion are in the v2 plan's
  "Open question for the UI step".
- **Where coloc's priors and the `sdY = 1` assumption live** (v2 gap G8). Under this plan a
  collection or tag is the natural home, since the posteriors belong to the pair rather than to
  either experiment.

## Steps

**Every decision is settled** (D1–D25, 2026-10-07). What follows is implementation order, and it is
ordered by what blocks what rather than by size.

1. **The format, since everything downstream encodes against it.** F1 only: `n_cs` and a
   credible-set record section in the trans frame (taking its reserved `u32`, so a frame without
   any is byte-identical) and in the GWAS block (payload header 40 -> 48). Done 2026-10-08, with
   round-trip tests for a fine-mapped locus in each and for the no-credible-set case costing
   nothing. No source in hand supplies any yet; the adapters' job when one does.
2. **`qtlstore.py`**: `kind` as a closed versioned enum, the `collection` digest field, the
   `parent` chain check, `store.json`'s `experiments` map, and `validate` over all of it — including
   that each map entry agrees with its experiment document, and **not** checking that collection
   members agree on anything (D24).
3. **`results.py` and `gwas.py`** emit one document per scan instead of one per cohort. The search
   index splits into a cis schema and a trans schema, `search_index_trans_only` is deleted, and hits
   split by record kind. This is the largest pipeline change.
4. ~~The converter~~ **Dropped** (Sam, 2026-10-09). It existed to migrate the live store without the
   contract tables; the tables are local and a genome-wide rebuild is nine minutes, so the
   migration is a rebuild and an upload.
5. **Rebuild and verify**: `build_store.py --version 2` genome-wide, `validate`, `which_codec`,
   `scan_codes`, and the v2 precision comparison re-run against the new objects.
6. **The reader**, then the two apps — which now diverge: the general app gains a chooser and a
   no-selection home, the TOPCHeF app pins its collection. `openStore`'s speculative-pointer
   machinery goes. Coloc priors and `sdY` move to deployment config.
7. **Re-check every gap** in the v2 plan, G1–G12, against the result. G11 is closed by step 1. G4
   stays open by decision (D26) and G12 stays open with a zero-cost path recorded.
8. Only then upload.

Steps 1–3 are built and tested; step 4 is dropped; step 5's genome-wide build is run after each
format change; steps 6–8 are not started. F2 was built, measured and reversed (D20). Collection
metadata is `collections:` in `pipeline/config.yaml`, with the papers in `data/raw/papers/`.
