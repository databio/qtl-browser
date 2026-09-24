# qtl-browser

A web browser for heart QTL results: the TOPCHeF cis and trans eQTL and sQTL summary statistics
(Murray et al. 2026, medRxiv 10.64898/2026.01.12.26343934), the GTEx v8 heart left ventricle
results from the eQTL Catalogue, and the Jurgens et al. 2024 DCM GWAS for colocalization views.

There is no server. The data is a **qtlstore**: a folder of files, each named by the hash of its
own bytes, with every position tied to a refget reference sequence. The site is a static React app
that reads those files with plain HTTP range requests and decodes them in the browser. `SPEC.md`
describes the store byte by byte.

Live site: https://topchef.databio.org

## Layout

| Path | What |
|---|---|
| `SPEC.md` | The qtlstore format (v1): store layout and object names, the file header, variant catalogs, annotations, allele orientation, results, search index, hits, trans, GWAS, and what `Store.validate` checks. |
| `pipeline/` | Python code that builds the store. Adapters turn each study into standard input tables (`pipeline/CONTRACT.md`); `qtlstore.py`, `catalog.py`, `annotation.py`, `results.py` and `gwas.py` write the store. `pipeline/README.md` explains every step. |
| `*.sbatch` | Slurm jobs for Rivanna: `adapter.sbatch` (adapter, TOPCHeF gate, contract check), `store.sbatch` (store build, validate, v0 comparison), `bench_store.sbatch` (format benchmark), `build.sbatch` (shared input tables). |
| `ui/` | The browser app: Vite, React, TypeScript, Tailwind 4 and DaisyUI 5, DuckDB-WASM, Mosaic plots. |
| `ui/bench/` | Browser smoke suites for the gene, variant and trans pages, and a Playwright harness that measures what a gene page costs. |
| `data/raw/` | `sources.yaml` lists every input with URLs, versions and checksums; `download.py` fetches them. The data itself is not in git. |
| `plans/` | Dated plans from earlier work. |

## Build the store

The build runs on Rivanna, where the inputs live. It goes in four stages:

1. **Adapters** turn each study into contract tables: `sbatch adapter.sbatch` (TOPCHeF, the eQTL
   Catalogue, the DCM GWAS).
2. **Store**: `sbatch store.sbatch` builds the annotation, variant catalog and results for each
   experiment into one store.
3. **Checks**: `Store.validate` (run by `store.sbatch`), `verify_v0` (the TOPCHeF store against the
   frozen v0 build), and the TOPCHeF acceptance gate plus the contract check (run by
   `adapter.sbatch`).
4. **Benchmark**: `sbatch bench_store.sbatch` compares v0 and v1 sizes and read speed; the smoke
   suites in `ui/bench/` check the pages in a browser.

`pipeline/README.md` has the exact commands and variables. The current whole-genome store is
`/scratch/ns5bc/qtl-browser/store-genome-v1e` on Rivanna.

Local setup:

```bash
uv sync                                   # Python dependencies (includes gtars for refget)
uv run python -m pipeline.test_qtlstore   # one of the test suites; pipeline/README.md lists them all
```

## Run the site locally

```bash
cd ui && npm install
QTL_DATA_DIR=/path/to/store VITE_DATA_BASE= npm run dev     # serves the store at /data with Range support
```

`ui/bench/README.md` shows how to copy a small store from Rivanna and run the smoke suites.

## Deploy

The data lives in the Backblaze B2 bucket `cloud-databio` under `qtl-browser/`, served at
`https://cloud2.databio.org/qtl-browser`. The site is a Cloudflare Workers static-assets project
(`ui/wrangler.jsonc`) in the databio account, at https://topchef.databio.org.

The full steps, with credentials and checks, are in the cloud-management repo:
`~/workspaces/assistant/cloud-management/backblaze/qtl-browser.md`. In short: upload the store
(`immutable/` files first, the pointer files and `store.json` last), then build and deploy the app:

```bash
cd ui
VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser npm run build
npx wrangler deploy
```

`ui/.env.production` already holds that `VITE_DATA_BASE`, so a plain `npm run build` reads the same
bucket. `VITE_DATA_BASE= npm run build` (empty) makes a bundle that reads `/data` on its own origin.

## Data notes

- Coordinates are GRCh38, tied to the reference by refget sequence digests. TOPCHeF genes are
  GENCODE v34; GTEx genes are GENCODE v39.
- Alleles are stored as `ref`/`alt` against the reference, with allele frequency and effect size
  for the ALT allele.
- eGene / sQTL intron: permutation p < 0.05, the preprint's wording.
- QTL values are rounded to keep files small; `SPEC.md` gives the error bounds. Exact values are in
  the source releases. GWAS values are not rounded.
- DCM GWAS: Jurgens 2024 biobanks-only meta-analysis (5,022 cases / 932,941 controls), the set the
  preprint's figures were drawn from.
- The v0 format (one study, `manifest.json`) is retired. Its spec and measurements are in the
  `analysis` repo at `qtlb-format/docs/`; `SPEC.md` section 17 says where the frozen v0 build lives.
