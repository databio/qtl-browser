---
date: 2026-10-07
status: draft
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

Three levels, replacing one.

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

### 3. Tags — shown together, ephemeral

Deliberately **not** a durable structure. "TOPCHeF's QTLs beside the DCM GWAS it colocalized
against" spans two publications and is a curation decision someone made for the browser, not a fact
about the data. Same for "these two results are worth looking at together" as a recommendation.

So this is lightweight tagging: cheap to add, cheap to remove, carrying no identity digests and no
validation guarantees beyond "the things named exist". Keeping it separate from collections is the
point — folding curation into provenance would recreate today's conflation one level up, which is
the whole defect this plan exists to remove.

## Decisions & ownership

| # | Decision | Ownership | Why |
|---|---|---|---|
| D1 | An experiment is one modality and one scope, so TOPCHeF becomes 4 experiments | **user-owned, surfaced** — answered 2026-10-07 | Sam: "each experiment gets its own single modality attempt... the 5 sumstats is 5 separate gwas experiments". tensorQTL's cis and trans are separate invocations, which is the line this follows |
| D2 | "Published together" gets a dedicated durable collection level; "shown together" is ephemeral tagging | **user-owned, surfaced** — answered 2026-10-07 | Sam: provenance "should get dedicated collection label or grouping", curation "much more flexible and more ephemeral like tagging two recommended results together". Typed separately so curation cannot silently become provenance |
| D3 | The DCM GWAS becomes its own experiment, not a field on TOPCHeF | **user-owned, surfaced** — answered 2026-10-07 | Sam: "the dcm gwas wasn't even part of topchef, so putting it under topchef experiment here doesn't make much sense" |
| D4 | Collections nest; `validate` gains a cycle and missing-member check | AI-owned, defended | Nesting is wanted (D2), and a pointer level that can reference its own kind can corrupt itself. The check is the difference between a working level and a trap |
| D5 | Settled before anything is uploaded under v2 | **user-owned, surfaced** | v2 already forces a full re-upload and nothing is published under it, so pointer restructuring is free now and a second migration later |
| D6 | Object bytes do not change; this is pointer-level JSON only | AI-owned, defended | Objects are content-addressed and carry no pointer identity. A restructured store reuses every existing object |
| D7 | Whether tags live in the store or outside it | **open / vague** | D2 says ephemeral, which argues for outside the content-addressed store entirely — a sidecar the browser reads. Not decided. If they go in `store.json`, "ephemeral" and "part of the validated store" are in tension |
| D8 | What the grouping level is called | AI-owned, default | `collections/` here, because Sam used "collection". `study`, `dataset` and `release` are all defensible and the word will end up in the UI and in citations |

## What this changes elsewhere

- **Search index and hits objects currently span phenotype types.** One `.arrow.zst` index and one
  `.qbh` per chromosome cover `ge` and `leafcutter` together. Split per experiment and a gene page
  goes from one index read to N, and from one hits read to N. **This is the only place the new model
  is more expensive at read time, and it should be measured before committing** — a gene page is
  already 98 kB in v1 / 111 kB in v2, 72% of it variant-catalog pages, so the question is how much
  N-way splitting adds to the fixed part.
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

## Steps

1. Settle D7 (tags in the store or beside it) and D8 (the name). Nothing else can be written down
   without them.
2. Measure the search-index and hits read cost of an N-way split on the genome-wide v2 store, since
   it is the one real objection to the model.
3. Write the pointer schemas — experiment, collection, tag — into SPEC, with the validate rules
   including D4's cycle check.
4. Migrate: a converter from a v1/v2 experiment document to the new set, reusing every object.
   TOPCHeF -> 4 experiments + 1 collection; the DCM GWAS -> 1 experiment + 1 collection; one tag
   joining them.
5. Reader and both apps.
6. Only then upload.

Nothing is implemented. This plan is a record of the design and the open questions, written before
the v2 store is uploaded so the restructure is not a second migration.
