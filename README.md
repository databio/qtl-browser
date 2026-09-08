# qtl-browser

Static browser for the TOPCHeF cis/trans eQTL and sQTL summary statistics (Murray et al. 2026,
medRxiv 10.64898/2026.01.12.26343934), with the Jurgens et al. 2024 DCM GWAS alongside for
colocalization views. All queries run in the browser with DuckDB-WASM over range-read parquet;
there is no server.

## Layout

| Path | What |
|---|---|
| `data/raw/` | `sources.yaml` lists every input (Zenodo QTL archives, GENCODE v34, dbSNP b157, DCM GWAS) with URLs, versions, and checksums; `download.py` fetches and verifies them. Everything else here is gitignored. |
| `pipeline/` | Python build that turns `data/raw/` into browser-ready parquet in `data/derived/`, plus `manifest.json`. `config.yaml` holds paths, the significance rule, window sizes. `figures.py` makes quick-look PNGs. |
| `ui/` | Vite + React + TypeScript app: Tailwind 4 and DaisyUI 5, DuckDB-WASM, Mosaic/vgplot plots. |
| `plans/` | Dated plans with decision ledgers and implementation logs. |

## Setup

```bash
uv sync                              # python deps
data/raw/download.py                 # ~45 GB of inputs; resumable, md5-verified
uv run python -m pipeline build      # -> data/derived/ (~7 GB), skips finished steps
uv run python -m pipeline validate   # counts vs the preprint, rsID checks, row-group pruning

cd ui && npm install
npm run dev                          # serves ../data/derived at /data with Range support
npm run build && npm run preview     # production bundle, same data plugin
```

`ui/.env.production` points production builds at the R2 bucket through `VITE_DATA_BASE`;
`VITE_DATA_BASE= npm run build` (empty) makes a bundle that reads `/data` on its own origin.

## Deploy

```bash
# .env at the repo root (gitignored): R2_ACCESS_KEY_ID and R2_SECRET_ACCESS_KEY,
# an R2 API token with Object Read & Write on the bucket. Run from the repo root.
uv run python -m pipeline.upload cors                # CORS rule from pipeline/config.yaml
uv run python -m pipeline.upload sync --dryrun       # then without --dryrun; re-runnable
uv run python -m pipeline.upload check               # HEAD + range GET on the public URL
```

Bucket, endpoint, public URL, allowed origins, and excludes are the `r2:` block in
`pipeline/config.yaml`. Setting the CORS policy needs an Admin token or the dashboard
(`upload.py cors --print` gives the JSON to paste).

The UI is a Cloudflare Workers static-assets project built from the repo: root directory
`ui`, build `npm run build`, deploy `npx wrangler deploy`, with `ui/wrangler.jsonc` naming
`dist` and the single-page-application fallback. Workers caps assets at 25 MiB, so the DuckDB
wasm modules are not bundled; the app loads DuckDB-WASM's jsDelivr bundles, as drumbeat-viewer
does.

## Data notes

- `manifest.json`, written last by the pipeline and published at the bucket root with the
  parquet, is the provenance record: build time, pipeline commit, significance rule, source
  versions, and for every table its path, row count, bytes, file count, columns, and whether
  the app loads it whole at startup or by range reads. The About page shows only the counts,
  source versions, and build date from it.
- Coordinates GRCh38; genes GENCODE v34; rsIDs from dbSNP by position and alleles.
- eGene / sQTL intron: permutation p < 0.05, the preprint's wording (10,220 eGenes vs the
  paper's 10,241; 13,540 sQTL introns, exact).
- Nominal cis tables are partitioned by chromosome and a 100-gene TSS bin, one row group per
  phenotype (gene for eQTL, intron for sQTL), so a gene page or a selected intron fetches a
  ~90 KB footer and its own row group. `gene_detail` holds the gene row, collapsed exons, and
  tested introns in one row per gene, partitioned the same way.
- eQTL nominal rows are kept for every tested gene; sQTL nominal rows only for the 13,540
  significant introns (`sqtl_nominal` in `pipeline/config.yaml`), which keeps the bucket near
  7 GB for the R2 free tier instead of 16. Every intron keeps its permutation row.
- The 21 eQTL and 4 sQTL colocalized genes are hard-coded from the authors' list until the coloc
  tables (PP.H4, sentinels) are shared; `coloc.parquet` is an empty stub for them.
- DCM GWAS: Jurgens 2024 biobanks-only meta-analysis (5,022 cases / 932,941 controls), the set
  the preprint's figures were drawn from although its Methods cite the full meta-analysis
  (9,365 / 946,368). Position-sorted parquet for LocusCompare plus a 5 Mb-binned table for
  the landing track; `pipeline/config.yaml` `dcm_gwas` picks the file, and the full-meta build
  is parked at `data/derived/_full/gwas_meta/`.
