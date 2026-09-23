# qtl-browser

Static browser for the TOPCHeF cis/trans eQTL and sQTL summary statistics (Murray et al. 2026,
medRxiv 10.64898/2026.01.12.26343934), with the Jurgens et al. 2024 DCM GWAS alongside for
colocalization views. There is no server, and the browser reads no parquet: every page fetches
binary pack files with plain HTTP range requests (`PACKS.md` explains them) and decodes them in
JavaScript. DuckDB-WASM holds only the tables built from those, for the plots and for paging.
Every file read at a byte offset is content-addressed: it lives under `immutable/` as
`<stem>.<sha16>.<ext>`, is cached for a year, and is never overwritten. `manifest.json` is the one
name whose contents change.

## Layout

| Path | What |
|---|---|
| `data/raw/` | `sources.yaml` lists every input (Zenodo QTL archives, GENCODE v34, dbSNP b157, DCM GWAS) with URLs, versions, and checksums; `download.py` fetches and verifies them. Everything else here is gitignored. |
| `pipeline/` | Python build that turns `data/raw/` into the browser's pack files, `search_index.arrow.zst` and `manifest.json`. Every pack and index is published to `data/derived/immutable/` under its content-addressed name, and `manifest.json` is the only thing that names those files. The tables the build makes on the way stay in `data/derived/_tables/` and are never uploaded. `config.yaml` holds paths, the significance rule, window sizes, the `r2:` bucket settings, and the `replaces:` map of old bucket paths. `figures.py` makes quick-look PNGs. |
| `ui/` | Vite + React + TypeScript app: Tailwind 4 and DaisyUI 5, DuckDB-WASM, Mosaic/vgplot plots. |
| `PACKS.md`* | Plain-language overview of the pack files: why they exist, which file answers which question, how a page reads them, the design choices, a glossary, and how to inspect or build packs with `pipeline/packtool.py`. Start here. |
| `SPEC.md` | The qtlb v1 store format (the qtlstore), byte by byte: store layout and object names, the 64-byte v1 header, the variant catalog and its identity digest, the annotation object, allele orientation, the experiment's results, search index, paged hits, trans results and GWAS, what `Store.validate` checks, and store maintenance. `pipeline/qtlstore.py`, `catalog.py`, `annotation.py`, `results.py` and `gwas.py` write it, with the codec `packfmt_v1.py`; `pipeline/CONTRACT.md` is the input tables they read. |
| `ui/bench/` | Playwright harness that loads fixed gene pages on the live site, a local preview, or a deploy rehearsal and records requests, bytes, and time until the locus plot is drawn; `results/` holds the committed baselines. The three browser smoke suites (`smoke_variant_page.mjs`, `smoke_gene_page.mjs`, `smoke_trans_tab.mjs`) live here too. |
| `plans/` | Dated plans with decision ledgers and implementation logs. |

\* `PACKS.md` and the v0 pack format are not in this repo. They live with the format's tests and
benchmarks in the `analysis` repo at `qtlb-format/docs/`: the v0 `SPEC.md` is in its git history
(commit `fe5a606`) and its measurements are `EVIDENCE.md`. The v0 build below (`manifest.json`,
`packfmt_v0.py`, `packcheck`) follows v0; references to "v0 `SPEC.md`" below mean that copy.

## Setup

```bash
uv sync                              # python deps
data/raw/download.py                 # ~45 GB of inputs; resumable, md5-verified
uv run python -m pipeline build      # -> data/derived/ (3.18 GB is uploaded; _tables/ adds 6.8 GB that is not)
uv run python -m pipeline validate   # 941 checks: counts vs the preprint, rsID text, every pack decoded, every published name and digest
uv run python -m pipeline.packtool header data/derived/immutable/eqtl.chr21.*.qbe   # look inside a pack (PACKS.md)
uv run python -m pipeline.packtool_example                                          # packs from tiny tables and back, no data needed

cd ui && npm install
npm run dev                          # serves ../data/derived at /data with Range support
npm run build && npm run preview     # production bundle, same data plugin
npm run bench -- --target live --label baseline    # gene page cost on the live site (ui/bench/README.md)
npm run bench:compare -- bench/results/baseline-live.json bench/results/<label>-live.json
```

`ui/.env.production` points production builds at the R2 bucket through `VITE_DATA_BASE`;
`VITE_DATA_BASE= npm run build` (empty) makes a bundle that reads `/data` on its own origin.

`uv sync` installs `gtars`, which the `refget_store` step uses to read a refgetstore. The
`reference:` block in `pipeline/config.yaml` names the reference FASTA, the store, and
`reference.collection`: the seqcol digest that pins which GRCh38 this release sits on. On a machine
with no store built yet, the first `uv run python -m pipeline build --step refget_store` builds one
from the FASTA and prints the collection digest to paste into `reference.collection`.

## Deploy

**Every command below that writes to the bucket needs Sam's explicit go-ahead.** He owns the bucket,
the R2 token, and the Workers project. `cors --print`, `inventory`, `budget`, `check`, any
`--dryrun`, and a bare `prune` are read-only and safe to run at any time.

A deploy is staged, rehearsed, then switched. New files go up before old ones come down, so the
site keeps serving the old build until one small file, `manifest.json`, is replaced.

```bash
# .env at the repo root (gitignored): R2_ACCESS_KEY_ID and R2_SECRET_ACCESS_KEY,
# an R2 API token with Object Read & Write on the bucket. Run from the repo root.
uv run python -m pipeline.upload cors                            # CORS rule from pipeline/config.yaml
uv run python -m pipeline.upload inventory -o bucket.json        # read-only listing; reuse it with --listing
uv run python -m pipeline.upload budget --listing bucket.json    # peak, final, GB-month against the free tier
uv run python -m pipeline.upload stage --dryrun                  # then without --dryrun: new immutable/ keys only
uv run python -m pipeline.upload check --staged                  # every staged file on the public URL

# rehearsal: a production bundle on a local preview, reading the real bucket through the staged
# manifest copy, so the new build is exercised end to end before anything live changes
(cd ui && npm run build && npm run preview)                      # leave it running
(cd ui && npm run bench -- --target rehearsal --label rehearsal --manifest immutable/manifest.<sha16>.json)

uv run python -m pipeline.upload release --dryrun                # then without --dryrun: the live switch
uv run python -m pipeline.upload check                           # live manifest equals local, plus the staged checks
# deploy the UI, then check the live site by hand
uv run python -m pipeline.upload prune                           # dry run: exactly which old keys would go
uv run python -m pipeline.upload prune --yes                     # the delete
uv run python -m pipeline.upload check --retired                 # the replaced prefixes are empty
```

- **Uploads follow `manifest.json`, never a directory walk.** The bucket is a copy of what the
  manifest names, so the underscore folders in `data/derived/` (`_tmp`, `_tables`, `_full`, `_old`,
  `_checks`, `_deploy`, `_retired`) cannot reach it by accident, and there is no exclude list to
  keep in step.
- **`stage` changes nothing the live site reads.** It uploads only new `immutable/` keys plus a
  content-addressed copy of the manifest, which is what the rehearsal points at. `release` uploads
  the changed plain-path JSON and then `manifest.json` last: that is the only moment visitors see a
  change, and `upload.py pointer <old manifest>` puts it back.
- **Headers are set at upload, never by a later server-side copy:**
  `public, max-age=31536000, immutable` under `immutable/`, `no-cache` for `manifest.json` and the
  plain-path JSON, `application/octet-stream` and no `Content-Encoding` for packs, whose byte
  offsets address the stored bytes. `upload.py` also sets the aws checksum variables, because
  aws-cli 2.23+ otherwise sends CRC headers that R2 rejects.
- **Sizes.** The rollout leaves 147 keys and 3,175,337,386 B in the bucket, 31.8% of R2's 10 GB free
  tier, replacing 7.26 GB of parquet. Because new files go up first, the switchover peaks at
  10.44 GB for as long as the old keys sit there; R2 bills the month's average, not the peak, so a
  few days cost nothing. `qtl_browser_for_sam.md` is the handover document for Sam, including the list of keys to remove.
- Bucket, endpoint, public URL, allowed origins, cache-control values, and the free-tier numbers are
  the `r2:` block in `pipeline/config.yaml`; the old bucket paths each new file takes over are the
  top-level `replaces:` map, which is what `prune` deletes and nothing else. Setting the CORS policy
  needs an Admin token or the dashboard (`upload.py cors --print` gives the JSON to paste).

The UI is a Cloudflare Workers static-assets project built from the repo: root directory
`ui`, build `npm run build`, with `ui/wrangler.jsonc` naming `dist` and the single-page-application
fallback. It deploys through Workers Builds on a push to `main`, or by hand with `npx wrangler
deploy` from `ui/`; a bad UI build rolls back from Workers Deployments, and a bad data build rolls
back with `upload.py pointer`. Workers caps assets at 25 MiB, so the DuckDB wasm modules are not
bundled; the app loads DuckDB-WASM's jsDelivr bundles, as drumbeat-viewer does.

## Data notes

- `manifest.json`, written last by the pipeline and published at the bucket root, is the
  provenance record: build time, pipeline commit, significance rule, source versions, and counts.
  `packs` names every pack file, the search index, the GWAS, rsID and variant indexes, and the dof
  the pack reader uses; `immutable` gives each published file its logical key, bytes, sha256, md5, and
  the old bucket prefixes it `replaces`; `precision` records the worst rounding error a reader can
  see. There is no `tables` block: every parquet table is a build intermediate. The About page shows
  the counts, source versions, build date, and rounding from it.
- Coordinates GRCh38; genes GENCODE v34; rsIDs from dbSNP by position and alleles.
- `manifest.json` also carries a `reference` block: the seqcol digest of the reference collection,
  each chromosome's refget sequence digest and length, and the result of checking every variant's
  alleles against those bases. It says which GRCh38 the positions are on, in a form another dataset
  or a refgetstore can be matched against. It is data about the release, not part of the pack
  format (v0 `SPEC.md` section 3; evidence in `EVIDENCE.md` A.10, analysis repo).
- eGene / sQTL intron: permutation p < 0.05, the preprint's wording (10,220 eGenes vs the
  paper's 10,241; 13,540 sQTL introns, exact).
- The gene page reads binary pack files (`immutable/`, format in v0 `SPEC.md`) with plain Range
  requests whose byte offsets come from `search_index`. A cold gene page sends 3 in parallel: the
  gene's eQTL block (gene row, exons, tested introns, eQTL rows, credible sets), its variants
  range, and its DCM GWAS window, plus one for its trans rows from the trans pack (none when the
  gene has no trans rows). Each sQTL intron the page shows adds one request. chrX genes skip the
  GWAS request, since the GWAS has no chrX rows.
- The sQTL pack holds nominal rows for all 80,750 tested introns, so the sQTL tab can show any of
  them. Every intron also keeps its permutation row.
- Each variants file holds the chromosome's cis-tested variants and then, as a second section,
  the variants seen only in trans, so every variant the study reports has an index.
- The variant page reads three files and no more: one range into the rsID index turns an rsID into
  a (chromosome, index), then the variant's page and its hits frame arrive together. The hits frame
  holds that variant's trans rows, the phenotypes it leads, and its credible-set memberships. The
  "Scan cis windows" button adds two range requests, one span of eQTL blocks and one of sQTL blocks.
- `search_index` is an Arrow IPC stream in one zstd frame, fetched whole at startup and handed
  straight to DuckDB, so the bundle needs no parquet reader and no DuckDB extension. Its schema
  metadata records the SHA-256 of every pack its offsets reach, and the `manifest` step refuses to
  write a manifest that disagrees with it, so index and packs can never go live out of step.
- **QTL values are rounded** to keep the packs small, and the manifest's `precision` block bounds
  how far: -log10 p is within 0.0016 (p within 0.37%), standard errors within 0.0047%, slopes within
  0.0045 of a standard error, and allele frequencies within 0.0000077. The About page and a note
  beside the cis table say so, and the cis download is named `*.rounded.csv`. Exact values are in
  the Zenodo release. GWAS values are not rounded: the GWAS pack is lossless at the source's own
  precision.
- The 21 eQTL and 4 sQTL colocalized genes are hard-coded from the authors' list until the coloc
  tables (PP.H4, sentinels) are shared.
- DCM GWAS: Jurgens 2024 biobanks-only meta-analysis (5,022 cases / 932,941 controls), the set
  the preprint's figures were drawn from although its Methods cite the full meta-analysis
  (9,365 / 946,368). The GWAS pack (position-sorted, one file per chromosome, with a small index
  the app loads at startup) feeds LocusCompare and takes over the old `gwas_dcm/` parquet, and a
  5 Mb-binned JSON feeds the landing track; `pipeline/config.yaml` `dcm_gwas` picks the file, and
  the full-meta build is parked at `data/derived/_full/gwas_meta/`.
