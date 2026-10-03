---
date: 2026-09-30
status: in-progress
model: Claude Opus 5 (1M context)
description: A catalog-overlap index object (kind 11, .qbo) — union site order plus a per-site membership mask — with the overlap statistics in its pointer
---

# Catalog overlap index

## What it is for

Three questions, none of which the store can answer today without merging catalogs by hand:

1. **Which catalogs hold this site?** Membership, in one read rather than a page lookup per catalog.
2. **How much do two catalogs overlap?** Exact counts, recomputed on every build, so the union and
   dedupe strategies stop being guesses. This is the "continuously characterize" part.
3. **This locus in study B.** Given a window, the other catalog's vidx for every shared site —
   what a cross-study browser view needs, without N per-site lookups.

`crosscat` samples and checks; it does not index and does not count. This object does both.

## Shape

A fourth pointer level, symmetric with the three that exist:

```
overlaps/<id>.json          mutable pointer: which catalogs, the statistics, the object name
immutable/<digest>.qbo      kind 11: union site order + per-catalog presence
```

### The object

```
[64-byte header]        kind 11, chrom "all", count = catalogs, seq_digest = collection
[directory]             chromosome table, union page offsets and first positions,
                        mask chunk offsets, per-chunk per-catalog prefix counts,
                        per-catalog cis/trans breakpoints
[union pages]           the union site list, keys only, paged like .qbv
[mask chunks]           one membership mask per union site, 65,536 sites a chunk
```

**Union order** is the canonical order of SPEC §5: chromosomes by `seq_digest` ASCII, then `pos`,
`ref`, `alt`. The same order `catalog_identity` digests, so it is already specified and already
implemented — and it is deliberately *not* any catalog's vidx order.

**Union pages carry the site key only** — position delta, allele code, heap — and none of the
per-catalog attributes (`af`, `ma_samples`, `ma_count`, `rs_number`, `match`). Those are properties
of a study's cohort, not of a site, and they live in the catalogs. Cost estimate below turns on
this: keys only is roughly 3–4 B/site against `.qbv`'s ~10.

**Presence** is one **membership mask per union site**, not one bitmap per catalog: `u8` for up to
8 catalogs, widening to `u16`/`u32`/`u64`, bit *i* set when catalog *i* holds the site. Masks are
stored per chunk of 65,536 union indices and zstd-framed.

This is the change the six-catalog horizon buys. With a handful of catalogs the masks are extremely
repetitive — most sites are in all of them, or in one recurring subset — so they compress far below
a byte a site, membership is a single byte to read with no rank at all, and the three-container
roaring machinery disappears. Alongside the masks, a per-chunk per-catalog **prefix count** table
(`k × n_chunks × u32`; about 5 KB at six catalogs) makes rank a prefix lookup plus a scan bounded by
one chunk.

The ceiling is 64 catalogs. Past that the masks should become per-catalog bitmaps, and SPEC should
say so rather than pretending the format scales forever.

### Recovering a vidx is not rank-select

The obvious move — catalog B's vidx is `rank(bitmap_B, u)` — **is wrong**, and this is the part of
the design most likely to be got wrong later. A catalog's vidx order is, per chromosome, the cis
section sorted by `(pos, ref, alt)` and *then* the trans-only section sorted the same way. Union
order interleaves them. So a rank over the presence bitmap gives the site's index among B's sites
in canonical order, which is not its vidx wherever a chromosome has any trans-only site.

Two fixes, and the second is cheaper:

- a second bitmap per catalog marking which of its sites are cis, so
  `vidx = rank_cis(u)` or `n_cis + rank_trans(u)`; or
- a per-catalog, per-chromosome table of `(union_index_first, vidx_first)` breakpoints, one entry
  per contiguous run — which for catalogs whose trans-only section is empty collapses to one entry
  per chromosome.

TOPCHeF has 19 trans-only sites on chr22 out of 121,157; GTEx has none at all. So the breakpoint
table is tiny in practice and degrades gracefully. Take that, and store the cis/trans split as
breakpoints.

### The pointer carries the statistics

```json
{"id": "grch38_heart", "object": "<digest>.qbo",
 "collection_digest": "EiFob05aCWgVU_B_Ae0cypnQut3cxUP1",
 "catalogs": [{"id": "topchef_grch38", "identity_digest": "...", "n_sites": 9182261},
              {"id": "gtex_v8_heart_lv_grch38", "identity_digest": "...", "n_sites": 9613073}],
 "union": 0,
 "pairs": [{"a": "topchef_grch38", "b": "gtex_v8_heart_lv_grch38",
            "shared": 0, "jaccard": 0.0, "a_only": 0, "b_only": 0}],
 "by_chrom": {"chr22": {"union": 148485, "shared": {"topchef_grch38|gtex_v8_heart_lv_grch38": 110598}}},
 "conflicts": {"shared_positions_no_shared_allele_pair": 0, "examples": []},
 "normalisation_suspects": {"count": 0, "examples": []}}
```

`identity_digest` per catalog is what makes a stale index detectable: if a catalog is rebuilt, its
digest changes and the index no longer describes the store.

## Steps

1. **SPEC §19** — the object, kind 11, union order, page and chunk layouts, mask widths, prefix
   counts, breakpoint tables, the pointer schema, and what `validate` checks. Update §3's kind
   table, §2's store layout, and §18's open list (this is a step toward the reserved VRS index, not
   a replacement for it).
2. **Docs and references** — every statement that the store has three mutable pointer levels:
   SPEC §2, `README.md`, `pipeline/README.md`. The README layout table and the maintenance command
   list gain the new level.
3. **`pipeline/overlap.py`** — `build` (merge catalogs into union order, emit pages, bitmaps,
   breakpoints, statistics), `decode`, `membership(site)`, `vidx_in(catalog, union_index)`,
   `counts()`. Reuses `catalog.read_sites` and `packfmt_v1`'s page codec.
4. **`qtlstore.py`** — `overlaps/` as a pointer level in `write_store` and `object_names`; `gc`
   protects it; `validate` checks kind, collection digest, that every named catalog is in the store
   with a matching `identity_digest`, that the union is exactly the merge of their sites, that each
   catalog's set-bit count equals its `n_sites`, and that the prefix counts agree with the masks.
5. **`pipeline/test_overlap.py`** — round trip on synthetic catalogs; the cis/trans vidx recovery
   case explicitly, since that is the trap; mask width selection at k = 8 and k = 9; rank across a
   chunk boundary; a stale-catalog digest failing validation.
6. **`store.sbatch`** — a phase after the catalogs, before or after results, rebuilding the index
   whenever the set of catalogs changed.
7. **First run** — genome-wide numbers on the two existing catalogs, recorded here.

## Decisions & ownership

| decision | tag | note |
|---|---|---|
| Build now, at two catalogs | **user-owned, decided 2026-09-30** | Sam is adding at least four more QTL studies. Building now means the cis/trans vidx trap is verified while the data is small enough to check by hand, and the design targets k ≈ 6 rather than k = 2. |
| Union pages carry site keys only, no attributes | AI-owned, defended | Attributes are cohort properties and already live in the catalogs. Keys-only is roughly 40 MB against roughly 110 MB for the two-catalog union. |
| Per-site membership mask, not per-catalog bitmaps | AI-owned, defended | Chosen once Sam said four more studies minimum. At k ≈ 6 the masks compress better than separate bitmaps, membership needs no rank, and roaring's three container types are not needed at all. Documented ceiling: 64 catalogs. |
| No `pyroaring` dependency | AI-owned, defended | Every other structure here is byte-specified and hand-decoded, and the browser will need a decoder too. |
| Merge per chromosome, never genome-wide | AI-owned, defended | Six catalogs at roughly 9.5M sites each is 57M keys; held as Python objects that is gigabytes. Per chromosome it is a few million, and union order is defined chromosome by chromosome anyway. |
| vidx recovery by per-chromosome breakpoints, not rank-select | AI-owned, defended | **The trap in this design.** Rank-select silently returns wrong vidx for any chromosome with trans-only sites. Named here so it is not rediscovered as a bug. |
| A fourth pointer level rather than a `store.json` key | AI-owned, default | Symmetric with `variant_catalogs/`, `annotations/`, `experiments/`, and gives the statistics somewhere to live. A bare key would be less code and less symmetry. |
| Pipeline only; no browser decoder in this change | AI-owned, defended | The browser use lands with the general UI work, and the object should be proven on real data first. |
| Normalisation suspects are a heuristic, not a verdict | AI-owned, default | Proper left-alignment needs the reference, so it only runs where a refgetstore is configured. Without one, report indel pairs with equal length delta within a small window as *suspects*, never as identity. |
| Kind number 11 | AI-owned, defended | 3 is free (v0's detail-less sQTL kind) but reusing it would make old and new stores ambiguous. |

## What this changes elsewhere

- **SPEC becomes a five-level format.** Three mutable pointer levels become four. Every statement
  of the form "three mutable levels" in SPEC §2, the README and `pipeline/README.md` is now wrong
  and has to be updated with it.
- **The index is invalidated by any catalog change.** Rebuild a catalog and the index silently
  describes a store that no longer exists. Mitigated by recording each catalog's `identity_digest`
  and failing `validate` on a mismatch — but it means `store.sbatch` ordering matters, and a
  partial build can leave a store that validates only because the index is absent.
- **`gc` must learn about it**, or it deletes the object the first time it runs.
- **Store size** grows by roughly 40–45 MB at two catalogs. The union grows sub-linearly — on chr22
  two catalogs union to 1.08x the larger of them — so six studies on comparable panels project to
  roughly 13–15M sites, or 55–65 MB of union pages plus a few MB of masks.
- **The statistics are O(k²) in the pointer**, not in the object: six catalogs is 15 pairs, and 15
  pairs across 23 chromosomes is 345 entries in `by_chrom`. Fine at six, worth revisiting at twenty.
- **Two experiments may share one catalog.** The index is over catalogs, not experiments, so studies
  whose `identity_digest` matches appear once. That is the dedupe story working, and the pointer
  should make it visible rather than looking like a missing study.
- **`crosscat` overlaps in purpose.** It stays as the sampled correctness check; the index is the
  exhaustive one. Worth deciding later whether `crosscat` becomes a consumer of the index rather
  than its own code path.
- **Build time.** Merging two 9M-site catalogs is a few minutes; N catalogs is an N-way merge and
  the memory floor is the union's keys.
- **Estimates in this plan are extrapolated from chr22**, the only measurement taken so far: union
  148,485 against 121,157 and 137,926 sites, 91.3%/80.2% overlap. Genome-wide union is projected at
  roughly 11M sites. The first run replaces every number here.

## Not in this plan

The browser decoder and any cross-study UI. The VRS-id index (SPEC §18), which this is a step
toward but does not deliver: VRS ids would let sites be matched across stores, not just across
catalogs in one store. Left-alignment normalisation of indels at ingestion, which is the actual fix
for the suspects this index will report.

## Implementation log

**2026-09-30.** Steps 1-6 done, uncommitted. Full pipeline suite green, no regressions from the
fourth pointer level: 14 / 8 / 11 / 19 / 12 / 5 / 4 / **11 new** / 8 / 7 / 28.

- **SPEC §19** written, with §2 (store layout, four levels), §3 (kind table, `qbo` = 11), §16 and
  §18 updated. Two gaps in my own spec surfaced while implementing against it and were fixed in it:
  the directory had no `n_chrom` field, and it did not say whether a union page may span
  chromosomes. It may not — each chromosome starts a new page, so `page_first_position` is
  unambiguous.
- **`pipeline/overlap.py`**: `build`, `decode`, `membership`, `rank`, `vidx_in`, `find`,
  `union_sites`, `encode/decode_union_page`, and a CLI.
- **`pipeline/qtlstore.py`**: `OVERLAPS` in `POINTER_DIRS`, `KIND_OVERLAP`, and `_check_overlap` as
  validate check 14. Adding the level to `POINTER_DIRS` gives `gc` protection for free, since
  `referenced()` already walks every level.
- **`store.sbatch`**: an `overlap` phase before `validate`, skippable with `OVERLAP=`.
- **Docs**: `pipeline/README.md` gains the command, the test in its list, and a paragraph on what
  the index is; the repo README and SPEC no longer claim three mutable levels.

### Two things the implementation taught the design

**The vidx trap is not hypothetical, and it hides.** The fixture gives `cat_a` a trans-only site on
chr1, so its vidx order is `1 A>C (0), 5 T>G (1), 3 G>A (2)` while union order is `1, 2, 3, 5`. At
`5 T>G` the rank is 3 and the vidx is 1. `cat_b` has no trans-only site, and for it rank and vidx
agree — which is exactly why an implementation returning the rank would pass a naive test suite.
`test_vidx_is_not_a_rank` asserts both halves, and round-trips every site of both catalogs against
the catalogs' own pages.

**Chromosomes sort by digest, not by name.** The first test run failed on an assumption that chr1
comes first; in the fixture chr2's `seq_digest` sorts first, so its sites lead the union. That is
§5's canonical order behaving correctly. The test now pins the digest order *and* asserts the
fixture still exercises it, so it cannot quietly stop testing the interesting case.

**One ordering bug, caught by wiring it into `store.sbatch`:** the CLI read `store.json` in order to
rewrite it, but on a fresh build no `store.json` exists when the overlap phase runs. The build now
writes only the object and its pointer and leaves `store.json` to the validate step, which is the
write order SPEC §2 mandates anyway.

### Step 7: first run, genome-wide (2026-09-30)

Both live catalogs pulled from B2 into a local store (199 MB of `.qbv`), built with
`python -m pipeline.overlap build --store … --id grch38_heart`. **83 s wall, 2.5 GB peak RSS.**

| | sites |
|---|---:|
| `gtex_v8_heart_lv_grch38` | 9,613,073 |
| `topchef_grch38` | 9,182,261 |
| **union** | **10,751,935** |
| shared | 8,043,399 — 87.6% of TOPCHeF, 83.7% of GTEx, Jaccard 0.748 |
| GTEx only | 1,569,674 |
| TOPCHeF only | 1,138,862 |

**Object: 30.9 MB, 2.878 B/site** — against the 40–45 MB the plan estimated.

| part | bytes | B/site |
|---|---:|---:|
| union pages (site keys) | 29,334,934 | 2.728 |
| mask chunks | 1,439,867 | 0.134 |
| directory | 173,728 | 0.016 |

The masks are the vindication of the k≈6 redesign: 10.75M raw bytes compress to 1.44 MB, **7.5x**,
because membership is overwhelmingly "in both" or one of two other patterns. Per-catalog roaring
containers would have carried three container types to beat a byte a site that zstd already beats.

**Breakpoints confirm the vidx design exactly as predicted.** GTEx has **23** — one per chromosome,
because it has no trans-only sites anywhere. TOPCHeF has **310**: 23 plus 287 more wherever a
trans-only site interleaves into union order. An implementation using rank as a vidx would be
correct on GTEx and wrong on TOPCHeF at 287 places.

`validate` check 14 passes on the built index.

**Conflicts: 3,964** positions where both catalogs have a site and share no allele pair — 0.049% of
the 8.0M shared positions, and the examples are ordinary different variants (`T>TCG` against `T>G`).
The chr22 hand measurement gave 48 of 110,646, or 0.043%: the same rate.

**Normalisation suspects: 2,897** within 50 bp, and they look real:

```
gtex 8,238,865 C>CG   vs  topchef 8,238,871 G>GT
gtex 9,014,809 C>CT   vs  topchef 9,014,810 T>TA
topchef 9,515,246 A>AT vs  gtex 9,515,250 T>TA
```

Single-base insertions a few bp apart in homopolymer or short-repeat context — exactly the shape of
one event written at two positions. 2,897 of 10.75M is 0.027%, so it changes no headline number, but
it is the first direct evidence that the contract's silence on left-alignment has a real cost. It
argues for requiring normalised indels at ingestion (CONTRACT.md), which is the actual fix; this
index only surfaces them.

**The chr22 extrapolation held.** The plan projected a ~11M union from chr22's ratio; the answer is
10,751,935, within 0.5%.
