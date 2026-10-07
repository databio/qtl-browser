---
date: 2026-10-06
status: draft
model: Claude Opus 5 (1M context)
description: Plasma pQTL studies (ARIC, MESA) into the qtlstore, built locally rather than on Rivanna
---

# Plasma pQTL sources into the store

Four candidate sources, all plasma proteome. Three are usable, one is not. Everything below was
measured on 2026-10-06 against the downloaded files, not read off a paper.

The store gains a phenotype that is not a transcript: a protein abundance measured by an affinity
assay. The contract already has the shape for it -- `phenotype_id` is the assay, `gene_id` is the
gene it belongs to, several phenotypes per gene -- which is the leafcutter-intron relationship.

## The sources

### ARIC plasma cis-pQTL (Zhang, Dutta, Chatterjee 2022, Nat Genet 10.1038/s41588-022-01051-w)

`data/raw/aric_pqtl_zhang2022/`, 929 MB. `EA.zip` and `AA.zip` each hold 4,659 members, one PLINK2
`.glm.linear` per SOMAmer, cis window +/-500 kb of the protein-coding gene's TSS. `seqid.txt` maps
SeqId to `uniprot_id` and a gene symbol.

| | EA | AA |
|---|---|---|
| proteins | 4,657 | 4,657 |
| rows | 12,850,239 | 22,983,529 |
| `OBS_CT` | 7,213 | 1,871 |

Columns: `#CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT BETA SE T_STAT P ERRCODE`. Everything the
contract wants is published: explicit REF/ALT, a real `SE`, per-variant N, rsIDs in `ID`, and
`ERRCODE` (`.` on every row of the member read). `TEST` is `ADD` throughout. Chromosomes are
unprefixed.

Two things settled against the local refgetstore, on all 2,814 rows of
`EA/SeqId_9201_13.PHENO1.glm.linear`:

- **`REF` is the reference base on 2,814 of 2,814 rows.** No allele swapping is needed, and the
  contract's `min_match_fraction: 0.999` is met with room. (Verify over the whole corpus during
  ingestion, not just this member.)
- **`A1` equals `REF` on 530 rows, 18.8%.** PLINK2 reported the effect against the minor allele, so
  `BETA` is REF-relative on nearly a fifth of rows. The adapter must negate `BETA` and mirror
  `A1_FREQ` to `1 - A1_FREQ` wherever `A1 == REF`. Nothing downstream can detect this being skipped;
  it would present as quietly degraded colocalization.

Build is GRCh38: 6 of 6 sampled REF alleles match hg38 via the UCSC API, and `classify` against the
local refgetstore agrees on every row of the member above.

Gene identity: `seqid.txt` publishes no ENSG. 4,657 SOMAmers carry 4,435 distinct symbols; against
GENCODE v34, 4,270 resolve to exactly one ENSG (96.3%), 161 match no v34 symbol, and 4 are
ambiguous across two ENSGs. The 161 are stale symbols, not ambiguity -- `C10orf54` (VSIR), `PVRL3`
(NECTIN3), `UFD1L` (UFD1), `FAM19A1` (TAFA1), `FDX1L` (FDX2), `ALPPL2` (ALPG), `HIST3H2A` -- so an
HGNC previous-symbol table resolves nearly all of them mechanically.

### MESA TOPMed plasma cis-pQTL (Schubert et al., Zenodo 6536687 v3)

`data/raw/zenodo_6536687/`, 1.2 GB, CC-BY-4.0. Five populations. Only the `allele_info` files are
used: their `all_cis` twins hold the same rows with `snps` as a bare `chr:pos` and no alleles at
all, and a row with no source alleles cannot become a site.

| population | rows | proteins | distinct sites | `statistic = 0` |
|---|---|---|---|---|
| AFA | 13,736,693 | 1,303 | 6,391,742 | 0 |
| HIS | 10,183,210 | 1,303 | 4,744,245 | 0 |
| CHN | 6,661,632 | 1,303 | 3,073,824 | **126** |
| EUR | 6,272,050 | 1,302 | 2,930,570 | 0 |
| ALL | 4,527,724 | 1,219 | 2,201,883 | 0 |

Columns: `snps gene statistic pvalue FDR beta ref_allele0 alt_allele1`. `gene` is
`SOMAmerID_ENSG.version` and parses on 100% of rows in every file, so unlike ARIC this source
publishes a gene id. There is no SE column; `statistic` is the t, so `se = beta / statistic` exactly
-- except on CHN's 126 zero-statistic rows, which have no derivable SE and must be dropped and
counted.

No `pvalue = 0` in any file, so the underflow case that `scan_codes.py` looks for does not arise
here. Minimum p is 3.9e-125 (EUR).

Build is GRCh38: 8 of 8 sampled `ref_allele0` values match hg38.

**Which allele `beta` counts is not stated in the files.** The `allele0`/`allele1` naming implies
`alt_allele1`, and that is an assumption, not a fact. See D9.

### MESA Olink PWAS fine-mapping (Krueger et al., Zenodo 15483845)

`data/raw/zenodo_15483845/`, 17 MB, CC-BY-4.0. Credible sets and PIPs only
(`{cis,trans}_finemap_<anc>.tsv`, `gene_id` an Olink OID, `variant_id` as `chr:pos:ref:alt`), plus a
five-method fine-mapping comparison table and protein-coding gene boundaries.

The record holds **no nominal summary statistics**. `CONTRACT.md` requires `nominal.parquet`, so this
cannot be an experiment. The other 48 files in the record are PrediXcan prediction models and are
not QTL results; they were deliberately not downloaded.

### MESA cis-eQTL (Dropbox folder) -- deferred to phase 2

`data/raw/mesa_eqtl_2020/MESA_cis-eQTLs.zip`, 11.1 GB holding 11.14 GB uncompressed (members are
stored, not deflated, so an adapter streams them straight out). Five population sets: AFA (3.07 GB),
HIS (2.46 GB), AFHI (2.12 GB), CAU (1.90 GB), ALL (1.58 GB) -- **two** of them pooled, `ALL` being
everything and `AFHI` African plus Hispanic. Fetched by hand 2026-10-07; the folder link serves an
on-demand zip that ignores Range, so `download.py` cannot get it and the enclosing zip's bytes are
not stable. Per-member md5s are in `members.md5` and are the identity.

Same cohort as `mesa_pqtl`, so plasma protein and transcript abundance are measurable in the same
people -- whether a pQTL acts through transcript level becomes answerable within one cohort.

Its bundled README settles two questions for the pQTL side as well, the two sources coming from the
same group and the same MatrixEQTL pipeline with the same 0/1-dosage naming: `alt` is "1 allele in
dosage (beta effect allele)", and `statistic` is "t-statistic from modelLINEAR", so
`se = beta / statistic`. `FDR` is Benjamini-Hochberg.

**It is GRCh37**, which is why it is deferred. The README says "rsid hg19" and "gene start hg19",
and 8 of 8 sampled `ref` alleles read hg19 rather than hg38.

A folder-level `dl=1` link streams an on-demand zip that ignores `Range`, so it cannot be listed,
sized, or resumed -- a probe for the first 400 bytes returned 5.7 GB before being killed. The URL is
in `sources.yaml` for provenance only; `download.py` reports the file `present` and never calls it.

## Local build, not Rivanna

Measured, not assumed:

- `pipeline/config.yaml` paths are already relative to the repo root, and `duckdb_memory_limit: 6GB`
  / `workers: 3` are laptop-sized. Rivanna's `store.sbatch` is what overrides upward.
- `gtars` 0.10.0 is installed in `.venv` and imports; Rust is present.
- 36 GB RAM, 14 cores, 165 GB free. `aws` installed, B2 keys in `.env`.
- The live store is 215 objects, 4.19 GB -- a download, not an obstacle.
- The refgetstore builds here in 11 s and prints the digest `config.yaml` demands,
  `EiFob05aCWgVU_B_Ae0cypnQut3cxUP1`, 195 sequences. 1.3 GB at `data/derived/_refget`.

What Rivanna has that this machine does not: the frozen v0 tree (used only by `verify_v0` and the
TOPCHeF acceptance gate, neither of which applies to a new study -- `SKIP_V0=1`) and the eQTL
Catalogue brick. Neither is in the way.

`data/derived` here is the September v0 parquet build, 24 GB, and no v1 `_tables/` ever existed
locally. Nothing needs the old build.

## Decisions & ownership

| # | Decision | Tag | Why / what it forces |
|---|---|---|---|
| D1 | Build locally; Rivanna not used for these sources | user-owned, confirmed 2026-10-06 | Measured above. Forces the B2 upload to happen from this machine. |
| D2 | Add to the live store `store-genome-v1f`, not a separate one | user-owned, decided 2026-10-06 | One store is what makes the overlap index and cross-catalog lookup span heart eQTL and plasma pQTL, which is the reason to add these at all. Forces: `store.json` and `overlaps/grch38.json` get rewritten in the bucket `topchef.databio.org` reads. See blast radius. |
| D3 | All 7 populations: ARIC EA, ARIC AA, MESA AFA/HIS/CHN/EUR/ALL | user-owned, decided 2026-10-06 | 77.2M rows, 7 variant catalogs, 9-catalog pairwise overlap, ~0.8-1.0 GB bucket growth. ARIC EA/AA are disjoint (`OBS_CT` exactly 7,213 and 1,871, summing to the site's stated N~9,000). **MESA `ALL` re-analyses the same participants as AFA+HIS+CHN+EUR pooled, and nothing in the store records that**, so any pairwise comparison will read it as independent replication. Unresolved sub-question: whether to note the pooling in its `ingestion.json`. |
| D3a | Sequencing: both adapters, then the general UI | user-owned, decided 2026-10-06 | All 7 experiments land before any UI work. Forces the MESA blockers (D9, D10) to be resolved up front, and means nothing is visible in a browser until the UI work that follows. |
| D4 | `phenotype_type` is `protein_somalogic` (ARIC, MESA) and `protein_olink` (if the 2026 record is ever used) | AI-owned, defended | Names the assay platform, because SOMAmer and Olink disagree on the same proteins and that is a real analytic caveat. `results.py` is generic over the value; the UI is not (D13). |
| D5 | ARIC: negate `BETA` and mirror `A1_FREQ` wherever `A1 == REF` | AI-owned, defended | Measured at 18.8%. Not doing it sign-flips a fifth of the betas, undetectably. |
| D6 | ARIC: keep `REF`/`ALT` as published | AI-owned, defended | Verified on 2,814/2,814 rows against the local refgetstore. Still re-verify corpus-wide at ingestion. |
| D7 | `gene_id` from the **UniProt accession** first, then the symbol: `pipeline/genemap.py` | AI-owned, decided 2026-10-07 | Resolved, and better than planned: **4,438 of 4,438 (symbol, accession) keys map**, so every ARIC SOMAmer gets a gene. 4,272 by a unique GENCODE v34 symbol, 164 by accession, 1 by an HGNC previous symbol, 1 by override. The accession rule is the one that matters: ARIC's symbol `PACAP` resolves by name to ADCYAP1 but its accession Q8WU39 is **MZB1**, whose previous symbol was also PACAP -- a name lookup would have attached a plasma protein to the wrong gene on a page that looked entirely normal. It also handles three Excel-mangled symbols (`3-Sep`, `10-Sep`, `11-Sep` -> the septins) and four immunoglobulin assays named by a gene list. One override remains: see D23. |
| D23 | SOD2 is pinned to `ENSG00000112096` by a recorded override | AI-owned, defended | Not resolvable from outside: GENCODE v34 carries **two** overlapping protein-coding genes named SOD2 (`ENSG00000112096` chr6:159,669,069-159,745,186 and `ENSG00000285441` chr6:159,679,119-159,762,529), and HGNC maps both the symbol and P04179 to `ENSG00000291237`, an id v34 does not have. So no authority can pick between the annotation's two candidates. `ENSG00000112096` is the long-standing id, the one GTEx and the literature use. The override lives in `genemap.OVERRIDES` with its reason, and a target the annotation lacks **raises at construction**, so a stale entry fails the build rather than going quiet. Annotation-version specific by nature -- a GENCODE v39 store may not need it. |
| D8 | MESA: `se = beta / statistic`; drop and count the 126 CHN rows where `statistic = 0` | AI-owned, defended | No SE is derivable there. Dropping is what the contract does with unusable rows, and the count goes in `ingestion.json`. |
| D9 | **MESA orientation: does `beta` count `alt_allele1`?** | **open / vague -- blocks MESA** | The files do not say. Resolve by sign concordance against ARIC on shared strong cis signals (both are GRCh38, ref/alt anchored, so sites join directly), or from the paper. Shipping the assumption unverified risks a wholesale sign flip on one of two studies -- worse than ARIC's 18.8%, because nothing would look wrong. |
| D10 | **MESA `n_samples` per population** | **open** | Needed for the dof fit. `dof.fit` reports `usable` and `implied_covariates`, so a wrong N shows up as an unusable fit rather than silently. A null dof means no slope can be rebuilt (`results.py:948`) **and no browser-side coloc**, since z comes from the nlp code plus dof. |
| D11 | Zenodo 15483845 is not an experiment | AI-owned, defended | No nominal rows to build from. Keep the files; revisit only if the contract ever grows a credible-sets-only experiment kind, which is a contract change, not an adapter. |
| D12 | MESA rsIDs: leave `rs_number` null | AI-owned, default | MESA publishes none. The local dbSNP b157 (28 GB, already here) could supply them via the `variants_rsid` path, at real cost. Challenge this if cross-study rsID search matters. ARIC publishes rsIDs in `ID` and keeps them. |
| D17 | Neither source ran a permutation pass; significance becomes metric-agnostic rather than permutation-shaped | user-owned, decided 2026-10-07 | See [Step 0](#step-0-significance-becomes-metric-agnostic). The store stops presuming the metric; each experiment declares its own. |
| D18 | The pQTL significance metric is a per-phenotype Bonferroni on the lead nominal p: `p_bonferroni = min(1, p_nominal_lead * n_variants)`, tested at 0.05 | AI-owned, defended | Both inputs are published (the lead is the minimum nominal p by definition; `n_variants` is recounted from nominal, which the contract already prefers). The variable correction folds into the column so a fixed threshold means the same thing for every phenotype. **Measured on TOPCHeF**, where both a permutation p and a lead nominal p exist for 19,423 genes: Bonferroni agrees with `p_perm < 0.05` on 93.3% of genes (8,910 significant vs 10,220). Rejected: a calibrated looser factor (median `Me/num_var` = 0.306, so Bonferroni over-corrects ~3.3x) -- that 0.306 is measured in heart tissue in a mostly European-ancestry cohort over ~6,400-variant windows, and ARIC AA / MESA AFA / MESA HIS have shorter LD blocks, so importing it would inflate significance for exactly those experiments in a direction nothing downstream would reveal. Also rejected: MESA's published `FDR`, a BH-across-tests quantity that is not comparable to a family-wise-within-phenotype one and would make MESA incomparable to ARIC. MESA's `FDR` is still carried in `permuted` as provenance. |
| D19 | No UI read fallback (`sig_value ?? p_perm`); clean cut | user-owned, decided 2026-10-07 | Neither app has been shared publicly, and store plus UI ship the same day. The transitional two-character read is not worth carrying. Consequence: between the store upload and the UI deploy, the gene page's significance column reads blank. Nobody is watching. |
| D20 | TOPCHeF's existing index objects are migrated in place by a new command, not rebuilt | AI-owned, defended | Renaming the index column would otherwise mean re-running `results.build` for TOPCHeF, which needs its contract tables -- those live on Rivanna (`/scratch/ns5bc/qtl-browser/derived/_tables/topchef`) and are tens of GB. The search index is a self-contained Arrow file: decode, rename the column, re-encode, new digest, rewrite the pointer. **25 objects per experiment** (1 whole index, 23 per-chromosome parts, 1 trans-only), 50 across both. No `.qbe`, `.qbv`, `.qbh` or `.qbg` object changes, so no re-upload of the 3.28 GB of results objects. |
| D13 | *Silently implied*: **nothing built here is visible in either app.** | silently implied | `ui/src/lib/store.ts:32` reads one experiment from `VITE_EXPERIMENT` (default `topchef`) and throws if `store.json` does not list it; `store.ts:35-36` hardcode `EQTL_TYPE = 'ge'` and `SQTL_TYPE = 'leafcutter'` as the gene page's two tabs. A `protein_somalogic` experiment is therefore unreachable and untabbed until the general UI makes both dynamic. |
| D14 | *Silently implied*: the coloc story becomes cross-tissue | silently implied | ARIC and MESA are plasma; TOPCHeF is heart. A DCM-GWAS colocalization against a plasma pQTL is a different claim than against a heart eQTL. `coloc.abf`'s `W1 = 0.15^2` prior was chosen for expression effect sizes and is not obviously right for an affinity-assay protein level. |
| D22 | MESA cis-eQTL is phase 2; the 7 GRCh38 experiments ship first | user-owned, decided 2026-10-07 | It is GRCh37, and the store's premise is refget-anchored GRCh38 -- `reference.collection = EiFob05aCWgVU_B_Ae0cypnQut3cxUP1`. A second refget collection would give two disjoint coordinate universes with no cross-catalog overlap, defeating the point. The route in is **not** a chain-file liftover: every row carries an rsID and dbSNP b157 is already local (28 GB, and `steps_variants` already does this lookup via `bcftools query -T`), so re-anchor by variant identity, then refcheck the published ref/alt at the new GRCh38 position and swap or drop per the contract. Open for phase 2: merged and deprecated rsIDs, and multi-allelic sites where one rsID maps to several allele pairs. |
| D21 | ARIC redistribution is fine; attribution is the obligation, not a license string | user-owned, decided 2026-10-07 | No license is stated anywhere, but the paper's data availability releases the sumstats unconditionally ("irrespective of significance level"), and its access restrictions attach to individual-level cohort data, not these files. Redistributing public sumstats is the field norm (eQTL Catalogue rehosts GTEx; OpenGWAS and GWAS Catalog rehost at scale). What this *does* oblige: the store carries the DOI, the URL and the checksum, so a derivative -- and ours is quantized and reoriented -- is traceable to its source. That is the `ingestion.json` provenance gap in the blast radius, which makes it a real work item rather than a legal one. |
| D15 | B2 upload | user-owned, surfaced, deferred | No uploader in the repo; the shelved B2 rewrite is on local `wip/b2-upload-profile` (`9b8b410`, unpushed). Upload needs its own explicit go-ahead and must write objects first, pointers and `store.json` last. |
| D16 | No `sdY` / trait-scale field; the assumption is documented instead | user-owned, decided 2026-10-07 | **Resolved, and nothing is mis-scaled today.** `coloc.abf`'s `W1 = 0.15^2` is a prior per *phenotype SD*, so it needs `sdY = 1` -- and all three results sets are, measured: tensorQTL inverse-normal transforms expression and splice ratios, and ARIC's phenotypes were standardised before PLINK2 (`sdY` estimate 0.9990 on a 2-parameter model, where the estimator is unbiased). A `phenotype_scale`/`sdY` field was designed and dropped: without a per-phenotype escape hatch every representable scale means `sdY = 1`, so the field could only hold `1.0`. The assumption is instead written in `coloc-abf.ts`'s `AbfPriors` docstring and the About page, with the trap recorded -- coloc's `sdY.est` returns residual SD, giving 0.5639 for TOPCHeF eQTL against a true 1, so "fixing" the prior with it would make the posteriors 3x wrong. See `2026-10-07-qtlstore-v2.md` G8. |

## What this changes elsewhere

- **The live store, if D2 says live.** `overlap.py` rebuilds across *every* catalog in the store, so
  adding 7 catalogs rewrites the overlap pointer that TOPCHeF's variant page reads, and
  `Store.validate` fails on an index whose catalogs no longer carry the identity it was built
  against. Object digests are content-addressed, so existing objects are untouched -- but
  `store.json` and `overlaps/grch38.json` are not objects, they are pointers, and they do change.
  The TOPCHeF app only asserts `store.experiments.includes('topchef')`, so it survives; this is a
  risk to the overlap index, not to TOPCHeF's results.
- **Bucket growth.** 77.2M rows at ~4.07 B/row is roughly 315 MB of results objects, plus 7 variant
  catalogs over 2.2-6.4M sites each. Order 0.8-1.0 GB on a 4.19 GB store. Known issue unchanged:
  `cloud2.databio.org` sends no ETag, so Chromium re-fetches every range read.
- **Local disk.** 165 GB free; raw pQTL sources are 2.2 GB and the store copy is 4.2 GB. The 24 GB
  v0 `data/derived` build is reclaimable if ever needed.
- **`sources.yaml` provenance does not reach the store.** `annotation.gtf_source`
  (`annotation.py:422`) gives an annotation real provenance from `sources.yaml`; an experiment gets
  only whatever its adapter puts in `ingestion.json["source"]`, which for the eQTL Catalogue is
  paths and accessions with no URL and no checksum (`eqtl_catalogue.py:657`). For ARIC -- whose
  host is an S3 bucket rather than an archive with a DOI -- that is not reconstructible, and D21
  makes traceability the thing we owe the source.
  Fix: have the adapter match its inputs in `sources.yaml` the way `gtf_source` does and copy
  name/version/url/md5/size into `ingestion.json["source"]["files"]`. Small, and it pays off exactly
  for non-archival sources.
- **Two adapters, not three.** ARIC (per-SOMAmer PLINK2 files inside a zip) and MESA (one long
  gzip per population, config block per population) are genuinely different source shapes. The 2026
  record is not an experiment.
- **`TYPES` in the eQTL Catalogue adapter** stays broken for single-phenotype-type datasets
  (`eqtl_catalogue.py:55`, hardcoded `("ge", "leafcutter")`). Unrelated to this work, but the pQTL
  adapters must not copy the pattern: both have exactly one phenotype type, so derive it from config.

## Step 0: significance becomes metric-agnostic

The reason this comes first: both pQTL adapters depend on it, and doing it now means TOPCHeF's store
is touched once rather than twice.

### Where the assumption sat

The rule was already data, not code -- `SPEC.md:1068` records that v0 had significance "fixed in
code" and v1 moved it into the experiment document so two experiments could differ. The leak was
narrower than it looked, in four places:

1. `results.py:82`, `SIG_COLUMNS = ("p_perm", "p_beta")`, a literal whitelist that raises. Its real
   job is guarding `getattr(g, col)` against an unvalidated config string -- `column: "count"` would
   return a bound method and then be compared with `<`, and `column: "pval_perm"` (the spelling
   `config.yaml` actually uses) would silently mark every phenotype non-significant. Input
   validation, not a claim about which metrics are legitimate.
2. `OPS = {"<", "<="}`: no metric where larger is better (a posterior, a Bayes factor, a PIP).
3. `contract_check.py:42` requires `p_perm` and `p_beta` to exist as `DOUBLE`, so a study with
   neither must still write both as null columns.
4. The **index column and the hit value are named `p_perm`** -- the one leak into the format itself
   (`SPEC.md:640,697`). And `results.py:378` writes `g.p_perm` unconditionally rather than the
   column the rule tested, so a `p_beta` rule already stores a boolean from one column beside a
   p-value from another. `test_results.py:491-504` pins that on purpose.

Nothing carries a *label*, which is why `Gene.tsx:250` hard-codes "Permutation p" and
`Gene.tsx:103-104` hard-code "permutation p < 0.05".

### The fix, by layer

The indirection resolves at **build** time, so no read-path component names a metric.

| Layer | After |
|---|---|
| `permuted.parquet` | `p_perm` and `p_beta` stay as **conventional, optional** names -- source vocabulary, so two tensorQTL-family studies agree on spelling, exactly as `CONTRACT.md`'s declared-attribute table already does for `ma_count` / `rsid` / `rs_number` / `match`. An adapter may add others (`p_bonferroni`, `fdr`). |
| `significance` in the experiment pointer | `{column, op, threshold, label}`. `column` names any numeric column of `permuted`; `label` is the display name. Gains a `null` meaning: the source assessed no significance. |
| Search index | `p_perm` -> **`sig_value`** (float64): one fixed generic slot holding the value the rule tested. `significant` becomes **nullable**; null means not assessed. |
| Hits, kind 0 | The 20-byte record's field is already named `value`; `SPEC.md:697` just redefines it from "`p_perm`" to "the value the rule tested". **Zero bytes change.** |
| Details `group` | Unchanged: keeps `p_perm` / `p_beta` as the source's own published numbers. Provenance, not semantics. Avoids a `DETAILS_VERSION` bump, which would rewrite all 92 `.qbe` objects (3.28 GB, 77% of the store). |

No `FORMAT_VERSION` bump. That byte sits in every binary header (`qtlstore.py:78`) and readers
reject a mismatch (`:88`, `store-decode.ts:142`), so bumping it rewrites every object's bytes and
every digest -- a full 4.24 GB re-upload for a column rename.

Once the builder writes the tested value into `sig_value`, **a reader never needs
`significance.column`**: it needs `label`, `op` and `threshold` to render "Bonferroni p < 0.05", and
`column` is pure provenance. That is what collapses the UI's ten call sites into one module.

### Changes

- `results.py`: `SIG_COLUMNS` -> validate that the rule names a numeric column present in
  `permuted` (swap lines 238-240, which compile the rule before the tables are read); `OPS` gains
  `>` and `>=`; `label` accepted; a null rule yields `significant = None`; `INDEX_SCHEMA` renames
  `p_perm` to `sig_value`; `:378`, `:400` and `:428` write the tested column's value.
- `contract_check.py:42`: require the rule's column rather than the two permutation names.
- `SPEC.md:527,640,697` and `CONTRACT.md`'s permuted section: the rule grammar, `sig_value`, the
  kind-0 meaning, `p_perm`/`p_beta` as conventional-and-optional.
- `results.py` gains `migrate-index --store --id`: decode each index object, rename the column,
  re-encode, rewrite the pointer (D20). 25 objects per experiment.
- `test_results.py:491-504` is rewritten: the index now reports the tested column, which is the
  opposite of what it currently asserts.
- **UI, both apps**: new `lib/significance.ts` exposing `{assessed, label, short, ruleText}` from the
  experiment pointer, and it becomes the only place those strings exist. Call sites:
  `Gene.tsx:250` (`sig.label`), `:317` (`sig.short`), `:94,103,104` (render iff `assessed`),
  `Search.tsx:86,87`, `Region.tsx:40,41`, `Variant.tsx:193` (badges hidden iff not assessed),
  `Genes.tsx:19` (default filter `assessed ? 'egenes' : 'tested'`, and drop the egenes/sQTL filter
  options), `Home.tsx:23` and `About.tsx:77` (omit the eGene clause when unassessed),
  `gene-index.ts:25,34` and `queries.ts:27,48` and `gene.ts:130,149` (field rename).
  Not covered by this and still TOPCHeF-specific in a shared component: `About.tsx:41`'s prose
  describing the permutation rule and `:53`'s sentence naming PJVK and CDKN1A.

## Steps

0. Step 0 above: the significance change, the index migration, and the UI consolidation.
1. ARIC adapter (`pipeline/adapters/aric_pqtl.py`), config block per cohort. Reads zip members
   without extracting, writes the five contract tables. Orientation per D5/D6, verified
   corpus-wide. `gene_id` per D7. Record the A1-flip count, the refcheck classes and the dropped
   rows in `ingestion.json`.
2. ~~Symbol -> ENSG lookup, built once and shared.~~ **Done**: `pipeline/genemap.py`, accession
   first, 4,438/4,438 (D7, D23).
3. `contract_check` on both ARIC cohorts.
4. Resolve D9 and D10, then the MESA adapter, config block per population.
5. Download the live store (368 objects, 4.24 GB), migrate its index objects (D20), store build of
   the 7 new experiments into it, `SKIP_V0=1`, chr21/22 first.
6. `Store.validate`, then `scan_codes.py` over the new experiments.
7. Adapter provenance from `sources.yaml` into `ingestion.json` (blast radius above).
8. Upload: separate, explicit, D15. Objects first, pointers and `store.json` last. Then `gc` the
   superseded index objects.

## Schedule

Target is end of day 2026-10-07. In likely-duration order, step 0 and the ARIC path are the parts
most likely to land; the store build and upload are where the time goes and are the most likely to
slip, because the 4.24 GB download, 77.2M rows across 7 experiments and 7 variant catalogs, an
overlap index over 9 catalogs, and a ~1 GB upload are each measured in tens of minutes to hours on
this machine. If the day runs short, the order above degrades gracefully: step 0 plus ARIC in a
local store is a complete, checkable result, and MESA is config blocks on a proven adapter.

## Implementation log

**Prerequisites, 2026-10-06.** Sources recorded in `data/raw/sources.yaml` with checksums, 2.9 GB
fetched and verified, refgetstore built and its collection digest confirmed
(`EiFob05aCWgVU_B_Ae0cypnQut3cxUP1`, 195 sequences, 11 s). ARIC's orientation checked end to end
against it: `classify` reads `REF` on 2,814 of 2,814 rows of one member, and `A1 == REF` on 530
(18.8%).

**Step 0, 2026-10-07.** Done. 92 pipeline tests green across 8 suites; both apps type-check and
build; `diff -rq ui/src ui-topchef/src` empty.

- `results.py`: the whitelist became a shape check against `permuted`'s real columns, which also
  catches the `pval_perm` spelling `config.yaml` uses; `OPS` gained `>` and `>=`; a null rule yields
  `significant = None` through the index, the details `group` and the type counts (null, not 0);
  `INDEX_SCHEMA`'s `p_perm` became `sig_value`, written from the tested column via a new
  `significance_value`; the hits kind-0 `value` likewise; the details `group` reads `p_perm` and
  `p_beta` with `getattr(..., None)` so a source with neither column does not crash; new
  `migrate-index` command.
- `contract_check.py`: `p_perm`/`p_beta` moved from `PERMUTED` into `PERMUTED_CONVENTIONAL`, checked
  only when present, plus a check that the rule's column exists and is numeric. `ingestion.json` is
  read once in the constructor.
- `SPEC.md`, `CONTRACT.md`: the rule grammar, `label`, the null meaning, `sig_value`, kind 0's
  meaning, and a dated decision note on the two conventional columns.
- Tests: `test_significance_uses_the_rules_column` rewritten (it had asserted the opposite),
  `test_significance_not_assessed` and `test_migrate_index_renames_p_perm` added.
- Both apps: new `lib/significance.ts` (pure) plus `useSignificance()` on the existing
  `StoreProvider`, which already held the experiment doc. Wired through `Gene.tsx`, `Genes.tsx`,
  `Search.tsx`, `Region.tsx`, `Variant.tsx`, `Home.tsx`, `About.tsx`, and the nullable `significant`
  through `cis-scan.ts`, `gene-index.ts`, `queries.ts`, `store-decode.ts`, `variant.ts`.

Worth knowing for the rest of this work: **`Row = Record<string, unknown>` (`lib/db.ts:13`), so
renaming a field on `Gene`, `SplicePhenotype` or `SearchHit` raises no type error.** `tsc` passed
with three stale `g.pval_perm` reads still live in `Gene.tsx`; only a grep found them.

**Step 2, 2026-10-07.** Done. `pipeline/genemap.py` plus `test_genemap.py` (10 tests, synthetic, no
data needed). HGNC's complete set added to `sources.yaml` (17 MB, rolled forward in place so its
md5 is its identity). Resolution order: override, unique GENCODE symbol, UniProt accession, HGNC
current symbol as a tie-break, HGNC previous symbol, HGNC alias. An HGNC record naming a gene the
annotation does not carry resolves to nothing, because a `gene_id` the annotation lacks would
dangle. Withdrawn HGNC records are skipped. ARIC: 4,438/4,438 (D7, D23).
