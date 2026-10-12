---
date: 2026-10-07
status: draft
model: Claude Opus 5 (1M context)
description: One encoding for p, beta and se across cis, trans and GWAS, prototyped locally on TOPCHeF and measured against the current codec and the source values
---

# A unified statistics encoding for cis, trans and GWAS

## Why

The store holds the same three numbers three different ways, and derives a different one of them each
time:

| | p | beta | se |
|---|---|---|---|
| cis `.qbe` | `u16` code / block `nlp_max` | **derived** `sign x se x t(nlp, dof)` | 15-bit log code + sign bit |
| trans `.qbt` | `u16` code / frame `nlp_max` | `i16` linear / frame `beta_max` | **derived** `|beta| / t(nlp, dof)` |
| GWAS `.qbg` | `u16` mantissa + `i8` exponent, lossless | `i32` = `rint(b x 1e4)`, lossless | `u16` = `rint(se x 1e4)`, lossless |

Every problem hit while building the ARIC adapter traced to the derivation, not to the encoding:

- **`p = 0` loses the effect size.** 240 of 495,450 ARIC rows on chr21/22 alone, concentrated in the
  strongest associations because that is what underflows. The row survives (`store-decode.ts:555`
  keeps it with `slope = NaN`), but its beta is gone -- not because the source withheld it, but
  because the format discards beta and rebuilds it through p.
- **`NLP_LIMIT = 300` fails the build outright.** 4 of ARIC's 170 chr21/22 blocks exceed it, worst at
  317.86. The limit exists only to keep the inverse-t inside double range; nothing about a `u16`
  code against an `f64 nlp_max` cares how large nlp gets.
- **A null `dof` costs an experiment every effect size.** ARIC's fit was being refused by an
  `at_edge` guard written for a lower-edge failure, and that alone would have cost it all
  colocalization. Fixed separately (see `2026-10-06-pqtl-sources.md`), but the fragility is the
  point: one fitted integer gates every beta in the experiment.
- **Trans has the mirror problem.** It stores beta and derives se, so a null dof costs it its SE and
  r2 instead (`store-decode.ts:768`).

`nlp`, `beta` and `se` are mutually determined: store two, derive one. A row with one of the stored
two unavailable is down to a single number and nothing can be derived. Swapping *which* two are kept
relocates the loss without removing it (measured: see "Rejected alternatives").

So: store all three, everywhere, in one encoding. Nothing is derived.

## The schema

### The trio, 6 bytes per row

| Field | Width | Code | Reserved |
|---|---|---|---|
| `nlp` | `u16` | `rint(-log10(p) / nlp_max * 65533)` | `65535` not tested, `65534` p underflowed to 0 in the source |
| `beta` | `i16` | `rint(beta / beta_max * 32766)` | `-32768` null |
| `se` | `u16` | `rint((ln se - lse_min) / (lse_max - lse_min) * 65534)` | `65535` null |

The rulers are each chosen for the quantity, and the reasons are measured (ARIC chr21/22, 170
blocks, within-block dynamic range `max/min`):

- `se` spans a median factor of **5.2** inside one block -- it is set by n, which is fixed, and MAF,
  which is a narrow band. A **log** ruler over that range gives a relative error of 1.3e-05 with 16
  bits. `se` is strictly positive, so the log is always defined.
- `|beta|` spans a median factor of **28,264** (95th pct 481,074), because it tracks z from ~0 up to
  the lead variant. It crosses zero, so a log ruler is unavailable; a **linear** ruler against
  `beta_max` gives a worst-case error of ~1.5e-03 of the row's own SE at `i16`, inside the codec's
  4.47e-03 slope budget. The bound generalises: the error in SE units is about
  `z_max / 65,532`, since `beta ~ z * se` and `se` is near-constant within a block. It is the
  geometry of a cis window, not a property of this dataset.

**The sign moves to beta**, where it is intrinsic to the `i16`. That frees bit 15 of the SE field,
so `se` gets 16 bits instead of 15 -- half the error it has today, for free.

### Scales, 32 bytes per scope

Four `f64`: `nlp_max`, `beta_max`, `lse_min`, `lse_max`. A scope is one cis block (a phenotype), one
trans frame (a phenotype), one GWAS block (a position range). Quantizing against the scope's own
range is what keeps the rulers tight.

### Per-row totals

Identity encoding is unchanged; only the statistics change.

| | identity | stats | total | today |
|---|---|---|---|---|
| cis | 0 -- rows are a contiguous `vidx` run from `var_start` | 6 | **6** | 4 |
| trans | 12 (`u32` pos delta, `u32` rs, `u16` af, `u8` chrom ordinal, `u8` allele) | 6 | **18** | 16 |
| GWAS | 12 (same, with `u8` n_code in place of the ordinal) | 6 | **18** | 21 |

### What is deleted

- `NLP_LIMIT = 300` and the inverse-t slope rebuild, with its `slope_err` measurement.
- `SE_SIGN`, bit 15 of the SE field.
- GWAS `p_mant`/`p_exp`, its `i32` beta, and `gwas_scaled`'s lossless enforcement.
- Trans's derived SE (`transSe` in the reader).
- The README's "GWAS values are not rounded" sentence, which exists to document the exception.

### What is added

- A `precision` block on GWAS objects, recording measured error per build, the way results sets
  already carry `neglog10p_max_error` and `slope_se_max_rel_error`. This is the replacement for
  losslessness: the error becomes a stated number rather than an unstated one.
- A consistency check: three dependent values quantized independently can disagree. See D4.

## Decisions & ownership

| # | Decision | Tag | Why / what it forces |
|---|---|---|---|
| D1 | GWAS joins the quantized trio; it stops being lossless | user-owned, decided 2026-10-07 | The README already states the principle -- "QTL values are rounded... Exact values are in the source releases" -- and GWAS was the exception to it. Where the two sides meet, in `logAbf`'s `r = W/(W+se^2)`, a 1.3e-05 relative error on se is invisible against a QTL side carrying ~1e-04. **Forces:** the store stops round-tripping Jurgens' published numbers bit for bit, so it is a rounded view like everything else; and `gwas_scaled`'s build-time refusal of over-precise input is lost, which D1a replaces. |
| D1a | GWAS objects gain a measured `precision` block | AI-owned, defended | `gwas_scaled` currently *fails the build* when a source carries more precision than `rint(x*1e4)` holds. That is a guard, not just a property, and dropping it would mean silently rounding an unexpected source. A stated per-build error is strictly better than today's QTL situation and no worse than losslessness, since the number is published either way. |
| D2 | `beta` is `i16`, not `i32`/`f32` | AI-owned, defended | Measured: `u16` linear gives a worst case of 7.6e-04 of an SE, `i16` about 1.5e-03, against a 4.47e-03 budget -- 0.00% of ARIC rows exceed it. `i32` would double the per-row cost of the trio's largest field to buy precision nothing needs. Trans already chose `i16`. |
| D3 | `se` keeps the log ruler and takes the freed sign bit | AI-owned, defended | Log because the within-block range is 5.2x and se is strictly positive; 16 bits because beta now carries the sign. Halves the SE error at no cost. |
| D4 | *Silently implied*: **the store can now be internally inconsistent** | silently implied -- surfaced | `nlp`, `beta` and `se` are mathematically dependent, and quantizing all three independently means a reader can find `beta/se` disagreeing with `t(nlp, dof)` by more than quantization error, with no single source of truth. Today consistency is guaranteed by construction, because beta is computed. This is what SPEC decision D13 was avoiding. **Mitigation:** a sweep in `scan_codes.py` asserting the triple agrees within the summed per-field bounds, and `precision` recording the worst observed disagreement. Without it the format trades a hard guarantee for an unchecked assumption. |
| D5 | `FORMAT_VERSION` bump | AI-owned, defended | The layout changes in all three object kinds. That byte sits in every object header (`qtlstore.py:78`) and readers reject a mismatch, so the bump re-digests variant catalogs, annotations, hits and indexes that do not themselves change -- a full ~5.5 GB re-upload. Trying to avoid it by flagging the layout in the dead `anchor` `i32` only saves the 18% of bytes that are not `.qbe`/`.qbt`/`.qbg`, and leaves old readers silently misreading the pairs array at the wrong stride. Take the bump. |
| D6 | +30% store, +50% on the cis read path | AI-owned, defended | 4.24 GB -> ~5.5 GB; a median TOPCHeF phenotype's block goes 25 KB -> 37 KB, and that block is what a gene page range-reads. Absolute cost is tens of KB per page, which compounds with the known `cloud2` ETag problem where range reads are never cached. Accepted because it buys the removal of every derivation failure. |
| D7 | `af` is unified to `rint(af * 65534)` everywhere | AI-owned, default | The variant catalog codes af `x65534` and GWAS codes it `x1e4` -- the same inconsistency one level down, free to fix while the layout changes. Challenge this if 4-decimal af matters somewhere. |
| D8 | The probe wires into nothing -- no UI, no B2, no live store | user-owned, decided 2026-10-07 | It writes its own objects under a scratch tree inside `data/` and reports numbers. Nothing in `ui/`, `ui-topchef/`, the bucket, or `store-genome-v1f` is touched. |
| D9 | The probe measures the codec; a **local store build** then measures the objects | AI-owned, corrected 2026-10-07 | An earlier version of this row claimed a build needs Rivanna. It does not. Only TOPCHeF's *already-built* contract tables are Rivanna-only; every raw input is local -- the unpacked Zenodo release (28 GB), dbSNP b157 (28 GB), the GENCODE GTF, the refgetstore built this session, `bcftools`, and crucially the September run's dbSNP join cache (`_tmp/dbsnp_matched.tsv.gz`, 96 MB) and `_tmp/variants_raw.parquet`. Nothing *looks* built because that run predates the `_tables/` reorganisation: its step markers all sit in `data/derived/.done/` while the current code reads `_tables/`, so `pipeline build` skips every step as done and finds no outputs. The adapter needs only `_tables/refcheck/<chr>.parquet` from `_tables` (`topchef.py:246`), reading statistics straight from the raw archives, so `gtf`, `permutation_tables`, `credible_sets` and the expensive `nominal` step are all unnecessary for a store build. Minimal chain: `variants_rsid --force`, `variants_refcheck`, `adapters.topchef`, then annotation + catalog + results. |
| D10 | Supersedes two queued items | AI-owned, defended | The sparse-beta-in-`anchor` idea and the t-statistic p recovery both become dead: a stored beta makes the first unnecessary and the second pointless. The `at_model_limit` dof fix stays useful on its own -- a fitted dof is still recorded as provenance and still wanted for `verify_v0`. |

## What this changes elsewhere

- **`packfmt_v1.py` -> v2.** `quantize_nlp` keeps its reserved codes; `quantize_se` loses the sign
  and gains a bit; a new `quantize_beta`. `encode_gene_block`, `encode_trans_frame`,
  `gwas_codes`/`encode_gwas_block` and all three decoders change. Block header grows by 8 bytes for
  `beta_max` (64 -> 72) or puts it in the details JSON; trans and GWAS headers likewise.
- **Builders.** `results.py` stops computing `slope_err` and starts encoding beta; `gwas.py` drops
  `gwas_scaled`'s enforcement and gains `precision`.
- **Reader.** `store-decode.ts` loses `tFromNlp` on the hot path, `transSe` entirely, and the
  `dof === null` branches; `coloc-abf` takes `z = beta/se` on both sides, one path.
- **`dof.py` becomes provenance.** The fit still runs and is still recorded, and `verify_v0` still
  wants it, but nothing decodes through it. The experiment pointer keeps `dof` and `dof_fit`.
- **`scan_codes.py`** gains the reserved-code sweep for the beta field and the D4 consistency check.
- **`bench_store.py`** is the natural place for the size comparison once a real build exists.
- **`verify_v0.py`** compares v0 packs against v1; it will need the v2 decode path.
- **SPEC.md** sections 8, 9, 10 and 13, plus the v0/v1 delta table; `CONTRACT.md` is **unaffected** --
  the adapter tables already carry beta, se and pvalue for every row, which is the whole reason this
  is possible without touching any adapter.
- **README** loses the GWAS rounding exception.

## The probe

A single module, `pipeline/schema_probe.py`, that reads the local source tables and reports. It
builds nothing into a store and imports no UI or upload code.

For each sampled scope (a phenotype for cis, a position range for GWAS):

1. Read the source `pval_nominal`, `slope`, `slope_se` (cis) or `p`, `beta`, `se` (GWAS).
2. **Current**: encode with today's scheme, decode, and rebuild beta through `t(nlp, dof)`.
3. **Unified**: encode with the trio, decode.
4. Compare both against the source -- the reference -- and report per-field error distributions in
   the format's own currency (error / se for beta, relative for se, absolute in -log10 p for p).
5. Report bytes per row for each, and the implied object and store totals.

It also counts what each scheme cannot represent: rows where the current scheme yields no beta
(`p = 0`, null dof, nlp > 300) against rows where the unified scheme yields no beta (null beta in
the source only).

Outputs a JSON blob and a markdown table, and the numbers go into this plan's implementation log --
not into a temp directory.

## Steps

1. `pipeline/schema_probe.py`: the two encoders, the decoders, the comparison, the CLI.
2. Run on a chromosome of eQTL, a chromosome of sQTL, and a chromosome of GWAS; then widen if the
   numbers look stable.
3. Record the results here: bytes/row, store projection, per-field accuracy against the reference,
   and the unrepresentable-row counts.
4. Add the D4 consistency check to the probe, so the inconsistency risk is quantified before it is
   designed around rather than after.
5. Decide, on those numbers, whether to proceed to a real v2 and in what order relative to the pQTL
   builds.

Nothing past step 5 is in scope here.

## Implementation log

**2026-10-07.** `pipeline/schema_probe.py` written and run. Nothing wired into `ui/`,
`ui-topchef/`, the bucket or the live store; it writes only JSON under `data/derived/_probe/`.

### Accuracy, against the source values

Error in the format's own currency: beta as a fraction of the row's own SE (the denomination of
SPEC section 9's **4.47e-03** budget), se as a relative error. `nlp` is unchanged by the proposal --
same quantizer -- and came out identical in every run, so it is omitted below.

| Scope | rows | scheme | beta err / se, median | **worst** | se rel err, median |
|---|---|---|---|---|---|
| TOPCHeF eQTL chr21 | 1,288,640 | current | 3.72e-05 | 2.28e-03 | 1.42e-05 |
| | | unified | 6.05e-05 | **6.35e-04** | **7.13e-06** |
| TOPCHeF sQTL chr21 | 902,323 | current | 4.75e-05 | 3.90e-03 | 1.39e-05 |
| | | unified | 6.12e-05 | **7.89e-04** | **6.99e-06** |
| DCM GWAS chr21 | 172,318 | current | **0.0** | **0.0** | **0.0** |
| | | unified | 1.10e-04 | 1.81e-03 | 1.25e-05 |
| ARIC EA chr22 | 331,128 | current | 2.55e-05 | **6.01e-03** | 1.26e-05 |
| | | unified | 5.60e-05 | **1.50e-03** | **6.32e-06** |

Three things to read off it:

- **The trade on cis is median against tail, and the tail is what the budget guards.** Unified's
  median beta error is 1.3-2.2x worse, because a derived beta inherits only se's very tight log error
  plus nlp's. But its worst case is **3.6x to 4.9x better**, because the inverse t amplifies error as
  |t| grows while a linear code's error is bounded by `beta_max / 65,532` no matter what.
- **se improves by exactly 2x everywhere**, which is the freed sign bit, as predicted.
- **GWAS is where unification actually costs something**: exactly zero error becomes ~1.1e-04 of an
  SE. That is the whole price of D1, stated as a number.

### The finding that was not part of the proposal

**On ARIC the current scheme exceeds its own budget.** Worst-case beta error 6.01e-03 against a
4.47e-03 ceiling -- because ARIC's |t| reaches 85.7 and the inverse t is steepest exactly there.
This is a defect in the format as it stands today, not an argument for the new one, and it means the
slope error bound in SPEC section 9 is not a bound for studies with ARIC's dynamic range.

### What each scheme cannot represent

| Scope | current: rows with no beta | current: scopes whose **build fails** | unified: rows with no beta |
|---|---|---|---|
| TOPCHeF eQTL chr21 | 0 | 0 | 0 |
| TOPCHeF sQTL chr21 | 0 | 0 | 0 |
| ARIC EA chr22 | **35** | **2 of 118** | **0** |

TOPCHeF exercises neither failure mode, which is why the probe had to be pointed at ARIC as well:
on TOPCHeF the proposal is pure cost. On ARIC the current scheme **fails the build** on 2 of 118
blocks (`NLP_LIMIT`) and silently drops 35 betas (`p = 0`); the unified scheme builds everything and
drops none.

### D4, quantified

Worst observed disagreement between the two now-independent routes to |t| -- `|beta|/se` against
`t(nlp, dof)`:

| Scope | worst | median scope's worst |
|---|---|---|
| TOPCHeF eQTL chr21 | 2.42e-03 | 3.43e-04 |
| TOPCHeF sQTL chr21 | 4.03e-03 | 3.92e-04 |
| ARIC EA chr22 | **2.17e-02** | — |

In absolute |t|. It grows with |t|, which is why ARIC is an order of magnitude worse than TOPCHeF --
at |t| = 85 a 0.022 disagreement is 0.026% relative, small but no longer negligible. This is the
cost of over-determination and the reason D4's consistency check belongs in `scan_codes.py` rather
than being assumed away.

### Storage, from a real local build

A validated v1 store for TOPCHeF chr21/22 was built locally (see below), which replaces the earlier
projection with measured numbers.

| | objects | MB |
|---|---|---|
| `.qbe` results | 4 | **88.85** |
| `.arrow.zst` index + annotation | 56 | 18.07 |
| `.qbt` trans | 2 | 8.44 |
| `.qbh` hits | 2 | 3.19 |
| `.qbr` rsID index | 1 | 2.87 |
| `.qgl` gene lookup | 1 | 2.79 |
| `.qbv` variant catalog | 23 | 2.44 |
| total | 90 | **126.65** |

`.qbe` is **4.064 bytes per row** over 21,863,503 rows, which matches SPEC's documented 4.07 and
confirms that essentially all of a results object is the pair array -- block headers and details
frames are rounding error. So +2 B/row is **+49% on `.qbe` bytes**, not less.

| | rows | now | unified | delta |
|---|---|---|---|---|
| local chr21/22 store | 21.9M | 126.65 MB | 170.4 MB | +35% |
| **live full store** | **~808M** | **4,240 MB** | **5,856 MB** | **+1,616 MB, +38%** |
| of which `.qbe` | | 3,284 MB | 4,900 MB | +1,616 MB |
| `.qbg` (21 -> 18 B/row) | 12.5M | 167 MB | ~143 MB | -24 MB |

**This corrects an earlier +29% in this plan.** That figure multiplied 2 bytes by TOPCHeF's
623,055,042 rows alone, but the live `.qbe` total covers *both* experiments: 3,284 MB at the measured
4.064 B/row implies ~808M rows, so ~185M belong to GTEx -- consistent with the ~153M `ge` nominal
rows `CONTRACT.md` records for it. The right multiplier is every experiment's rows, and the growth
is **38%**.

The local store over-states growth slightly on its own (35% vs 38% the other way) because a two-
chromosome subset carries proportionally more index and annotation bytes: 30% non-`.qbe` against the
full store's 23%.

### The local v1 build (D9, resolved)

Built end to end on this machine in about five minutes of compute, no cluster:

| Step | Result |
|---|---|
| `variants_rsid --force` | 9,215,026 variants; the September dbSNP join cache was reusable |
| `variants_refcheck` | 0.2 min, and **reproduces the Rivanna build exactly** -- A1 is the reference for 97,433 cis variants, A2 for 8,775,290, both 0, neither 0, strand flips 0, match fraction 1.00000, the same counts `qtlstore.orient_to_ref`'s docstring records from the cluster |
| `adapters.topchef` chr21/22 | 0.3 min; 239,260 sites, 21,863,503 nominal rows, 100,226 phenotypes, 399,828 trans rows; cis orientation `swapped: 0, dropped: 0` |
| `contract_check` | **RESULT PASS**, every rule |
| annotation + catalog + results | 90 objects, 126.65 MB |
| `Store.validate` | **PASS with empty `notes`** -- a full pass, so the refget anchor check and the site-reader identity check both ran rather than being skipped |

Two things this build settles beyond the sizes:

- **It independently confirms the probe's error measurement.** The builder's own
  `slope_max_error_over_se` came out at **3.9008e-03** for leafcutter and **2.4512e-03** for ge. The
  probe, measuring the same thing by a completely different route on chr21 alone, reported
  **3.901e-03** (sQTL) and **2.276e-03** (eQTL). The sQTL figures agree to four significant digits.
  So the probe's "current scheme" numbers are not an artifact of how it reimplements the codec.
- **It is the first store ever built with the metric-agnostic significance change** (step 0 of
  `2026-10-06-pqtl-sources.md`), and it validates clean, with `sig_value` in the index and real
  counts (`ge` 372 of 686 significant, `leafcutter` 512 of 2,862).

### A bug this exposed, unrelated to the proposal

`results._vidx` crashed on any chromosome with no credible sets: an empty pandas slice does not keep
its column dtypes, so `ref` and `alt` came back as `float64` and the merge against the catalog's
string columns raised `"trying to merge on str and float64"` instead of yielding nothing.

**It blocked the ARIC build entirely.** ARIC's `credible_sets.parquet` has zero rows by design --
`CONTRACT.md` requires the empty table rather than no table for a source with no fine-mapping -- so
every chromosome would have hit it. Fixed with an early return, with
`test_builds_with_no_credible_sets` covering both the empty-chromosome and empty-experiment cases.

### Probe limitations (D9)

It measures the codec, not a build. Not covered: post-zstd sizes for trans, hits and index objects,
real gene-page read cost, and the trans frame's own re-encoding (local TOPCHeF trans tables were not
exercised). Those need a Rivanna build. The GWAS reference is the source's 4-decimal value
recovered by rounding, because the local table stores beta and se as **float32** and
`gwas_codes` correctly refuses that as over-precise -- an error measured against float32 noise would
not be an error in the codec.

### Where this leaves the decision

The proposal is a clear win on ARIC and a clear cost on TOPCHeF and the GWAS. It removes a failure
mode that currently **fails the build**, it tightens the worst case by ~4x where it matters most,
and it is the only option that leaves `p = 0` rows usable. Against that: +29% store, +50% on the cis
read path, a `FORMAT_VERSION` bump with a full re-upload, and a new consistency obligation.

Not yet decided: whether to proceed, and whether before or after the pQTL builds. Step 5 of the plan
stands.
