# qtl-browser

A web browser for heart QTL results: the TOPCHeF cis and trans eQTL and sQTL summary statistics
(Murray et al. 2026, medRxiv 10.64898/2026.01.12.26343934), the GTEx v8 heart left ventricle
results from the eQTL Catalogue, and the Jurgens et al. 2024 DCM GWAS for colocalization views.

Live site: https://topchef.databio.org

The site is a static React app with no backend. Its data is a **qtlstore**: a directory of files,
each named by the hash of its own bytes, with every position tied to a refget reference sequence.
One store can hold several studies over several variant catalogs. The app fetches the parts of those
files it needs with HTTP range requests and decodes them in the browser. `SPEC.md` has the byte
layout.

## Layout

| Path | What |
|---|---|
| `SPEC.md` | The qtlstore format (version 2): store layout and object names, the file header, variant catalogs, annotations, allele orientation, experiments and collections, results, search index, hits, trans, GWAS, precision, and what `Store.validate` checks. |
| `pipeline/` | Python code that builds the store. Adapters turn each study into standard input tables (`pipeline/CONTRACT.md`); `qtlstore.py`, `catalog.py`, `annotation.py`, `results.py` and `gwas.py` write the store, `build_store.py` runs them in one command, `upload.py` puts the result in the bucket. `pipeline/README.md` explains every step. |
| `*.sbatch` | Slurm jobs for Rivanna: `adapter.sbatch` (adapter, contract check), `store.sbatch` (store build, validate, codec check), `build.sbatch` (shared input tables). |
| `ui-topchef/` | The TOPCHeF browser app at https://topchef.databio.org: Vite, React, TypeScript, Tailwind 4 and DaisyUI 5, DuckDB-WASM, Mosaic plots. |
| `ui-topchef/bench/` | Browser smoke suites for the gene, variant and trans pages, and a Playwright harness that measures what a gene page costs. |
| `ui/` | The general multi-study browser, forked from `ui-topchef/` and being generalised: the experiment becomes part of the route rather than a build-time constant. Deployed to its own Worker, with no custom domain yet. |
| `.github/workflows/` | `deploy-topchef.yml` and `deploy-ui.yml`: build and deploy each app to Cloudflare Workers on push to main. |
| `data/raw/` | `sources.yaml` lists every input with URLs, versions and checksums; `download.py` fetches them. The data itself is not in git. |
| `plans/` | Dated plans from earlier work. |

## Build the store

The adapters run on Rivanna, where the inputs live; the store builds anywhere the contract tables
are, in about ten minutes genome-wide:

1. **Adapters.** `sbatch adapter.sbatch` turns each study into contract tables (TOPCHeF, the eQTL
   Catalogue, the DCM GWAS).
2. **Store.** `uv run python -m pipeline.build_store --store DIR --validate` builds the annotation,
   variant catalog, one experiment per scan and the overlap index from the configured tables, then
   validates; `sbatch store.sbatch` does the same on Rivanna for several cohorts.
3. **Checks.** `Store.validate` and `which_codec` run inside both; `scan_codes` sweeps every code
   in a cis scan's blocks on demand; the contract check runs inside `adapter.sbatch`.
4. **Upload.** `uv run python -m pipeline.upload plan --store DIR`, then `push --yes`, then
   `check`; see Deploy below.

`pipeline/README.md` has the exact commands and variables. The current whole-genome store is
`data/derived/store-release-v2` (TOPCHeF, the DCM GWAS, GTEx v8 heart LV and the coloc lists),
built locally; the bucket serves the version 1 store until the version 2 release.

Local setup:

```bash
uv sync                                   # Python dependencies (includes gtars for refget)
uv run python -m pipeline.test_qtlstore   # one of the test suites; pipeline/README.md lists them all
```

## Run the site locally

```bash
cd ui-topchef && npm install     # or `cd ui` for the general browser
QTL_DATA_DIR=/path/to/store VITE_DATA_BASE= npm run dev     # serves the store at /data with Range support
```

`ui-topchef/bench/README.md` shows how to copy a small store from Rivanna and run the smoke suites.

## Deploy

The data lives in the Backblaze B2 bucket `cloud-databio` under `qtl-browser/`, served at
`https://cloud2.databio.org/qtl-browser`. Each app is a Cloudflare Workers static-assets project in
the databio account:

| App | Worker | URL | Workflow |
|---|---|---|---|
| `ui-topchef/` | `qtl-browser-topchef` | https://topchef.databio.org | `.github/workflows/deploy-topchef.yml` |
| `ui/` | `qtl-browser` | workers.dev only | `.github/workflows/deploy-ui.yml` |

A push to main that touches an app's directory builds it (`npm ci && npm run build`) and runs
`wrangler deploy` there; `workflow_dispatch` reruns either by hand. Both workflows read the repo
secrets `CLOUDFLARE_API_TOKEN` (Workers edit permission, plus Workers Routes and DNS edit on the
databio.org zone for the custom domain) and `CLOUDFLARE_ACCOUNT_ID`. `topchef.databio.org` is
declared in `ui-topchef/wrangler.jsonc`, so a deploy of that app attaches it to
`qtl-browser-topchef`.

**Known issue:** `cloud2.databio.org` sends no `ETag` or `Last-Modified`, and Chromium
stores a 206 (range) response only with a strong validator, so range reads are re-fetched on every
visit instead of coming from the browser cache. Whole-file reads are cached normally.

The uploader is `pipeline/upload.py`: `plan` says what would change, `push --yes` sends new
`immutable/` objects first, then the pointer files, then `store.json` last, `check` reads the
result back over the public URL, and `prune --yes` removes what the store no longer names. It
takes `B2_ACCESS_KEY_ID` and `B2_SECRET_ACCESS_KEY` from the environment or the gitignored `.env`
(the key is in the `databio/secrets` pass store). A store upload needs no app redeploy. To deploy
an app by hand instead of through Actions:

```bash
cd ui-topchef            # or ui
npm run build
npx wrangler deploy
```

Each app's `.env.production` sets `VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser`, so a
plain `npm run build` reads the B2 bucket. `VITE_DATA_BASE= npm run build` (empty) makes a bundle that reads `/data` on its own origin.

## Data notes

- Coordinates are GRCh38, tied to the reference by refget sequence digests. TOPCHeF genes are
  GENCODE v34; GTEx genes are GENCODE v39.
- Alleles are stored as `ref`/`alt` against the reference, with allele frequency and effect size
  for the ALT allele.
- eGene / sQTL intron: permutation p < 0.05, the preprint's wording.
- Every stored statistic, QTL and GWAS alike, is quantized to 16-bit codes; each experiment records
  the worst error that introduced (`SPEC.md` section 13). Exact values are in the source releases.
- One experiment is one scan: TOPCHeF is four experiments (eQTL and sQTL, cis and trans) in one
  collection, and the DCM GWAS is a fifth in its own, since TOPCHeF did not produce it.
- DCM GWAS: Jurgens 2024 biobanks-only meta-analysis (5,022 cases / 932,941 controls), the set the
  preprint's figures were drawn from.
- Earlier formats: version 1 is in git history at `be0bb44` and in the bucket until the release;
  v0 (one study, `manifest.json`) is in the `analysis` repo at `qtlb-format/docs/`. `SPEC.md`
  section 17 says where both live.
