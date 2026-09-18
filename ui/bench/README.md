# bench

Measures what a gene page costs: requests, bytes, and time until the locus plot is drawn. It runs
the real app in headless Chromium (Playwright), so React, DuckDB-WASM in its worker, and the
Mosaic plots all run as they do for a visitor. Plans that change how the gene page reads data run
the same pages with a new `--label` and compare against the committed baseline.

```bash
npm run bench -- --target live --label baseline            # 3 runs, all pages, all scenarios
npm run bench -- --target preview --label baseline
npm run bench -- --target live --label x --runs 1 --pages FLNC-eqtl --scenario cold,nav
npm run bench:compare -- bench/results/baseline-live.json bench/results/<label>-live.json
npm run bench -- --target preview --label x --latency 90 --scenario cold,nav   # 90 ms added per request
npm run bench -- --target preview --label pair --scenario pair --pages MYOZ1-eqtl,SYNPO2L-sqtl
npm run bench -- --target rehearsal --label r2 --manifest immutable/manifest.<sha16>.json
```

- `--query k=v&...` adds parameters to every page path, and to `/genes` in `nav`. The app reads no
  flags today (Plan 2's `?pack` flag is gone); Plan 2's results used it to pick the data path.
- `--manifest <key>` names the manifest the target serves. It applies to `--target rehearsal` only,
  where the bucket holds a staged copy (`immutable/manifest.<sha16>.json`) that `manifest.json`
  itself has not been switched to yet.
- `--latency <ms>` adds that round-trip latency to every request of the tab with CDP
  `Network.emulateNetworkConditions`. Local serving hides the cost of serial requests, which
  dominates on R2. With latency on, the self-check also reports how many DuckDB worker XHRs
  waited at least 90% of it for headers, and stops if any did not.

Playwright is pinned in `ui/package.json`, which pins the Chromium build. After `npm install`, run
`npx playwright install chromium` once (the browser goes to `~/.cache/ms-playwright`).

## Pages

`pages.json` is fixed: FLNC and MYOZ1 eQTL, SYNPO2L and CAMK2D sQTL (the sQTL tab draws the first
significant intron's locus without a click). Do not edit existing entries; add new ones under new
ids so old results stay comparable. The chr7 pack prototype added AC011287.1 (the largest chr7 eQTL
window, 10,419 variants), PILRB (ENSG00000121716, the smallest chr7 p), and AC069288.1 on its sQTL
tab (tested for sQTL only). The gene page on packs (Plan 5) added SKI (chr1), TBC1D5 (chr3, one of
the genes whose variants range covers intron runs beyond its eQTL run), AC244197.2 (chrX, which has
no GWAS rows, so no GWAS request), AL031282.2 on its sQTL tab (tested for sQTL only), and HSPG2 on
its sQTL tab. The trans pack (Plan 8) added FLNC-trans and HHATL-trans-sqtl, which wait for the trans
table instead of the locus plot (`ready` below): FLNC's 87 trans eQTL rows, and HHATL's 65,065 trans
sQTL rows, the largest trans frame.

An entry may set `ready`, a CSS selector that replaces the locus scatter selector for that page
(see Readiness contract). The trans entries use `[data-trans-total="<n>"]`, which `TransTable`
sets on its root once the table has loaded, with `n` the unfiltered row count of that tab's QTL
type.

An entry may also set `click` (a CSS selector) with `pre_ready` (another selector). The harness
then waits for `pre_ready`, clicks `click`, and only then waits for `ready`, so one entry measures
a page plus the action taken on it. `variant-rs10824026-scan` is the only user: it waits for the
variant page's lists, presses "Scan cis windows", and is ready when the scan renders. Its ready
time therefore includes the page it acts on, so compare it with `variant-rs10824026`, not on its
own.

The variant entries use `[data-variant-lists="1"]` (the lead and credible-set sections have their
data, or show the outside-cis message) together with `[data-trans-total]`, and the scan entry uses
`[data-scan-ready="1"]`. `Variant.tsx` sets all three.

Data kinds (`hits`, `rsid_index`, `variant_index`) name the variant page's packs. `parquet_other`
is the catch-all for any remaining `.parquet` request and is matched last. No parquet is published
any more, so it is a tripwire: any value above zero means a page found one.

## Targets

- `live`: `https://qtl-browser.topchef.workers.dev`, which reads the R2 bucket named in
  `ui/.env.production`.
- `preview`: `http://localhost:4173`, serving `../data/derived` at `/data` with Range support
  (`vite.config.ts`). Build the bundle for local data first; a plain `npm run build` reads R2:

  ```bash
  cd ui && VITE_DATA_BASE= npm run build && npm run preview
  ```

  The harness does not start the server. It stops if `/data/manifest.json` does not answer 200,
  if any data request goes to a non-local host, or if any data response is 404.

  The local `/data` server now sends the same cache headers as the bucket (SPEC section 3):
  `public, max-age=31536000, immutable` for anything under `immutable/`, and `no-cache` for the
  rest. It also sends an `ETag` and a `Last-Modified`, as R2 does: Chromium stores a 206 only when
  the response carries a strong validator, so without them `warm` re-downloaded every pack range
  while reading whole files from cache. Every file a page reads at a byte offset has an `immutable/`
  name, so `warm` on preview reads its packs from the browser cache, the way `live` does. Cold runs
  stay cold because each one gets a fresh browser context. Results recorded before Plan 6 were taken
  with `no-store`, where `warm` on preview downloaded everything again, so do not compare `warm`
  across that line. The local server does not answer `If-None-Match`, so `warm`'s `manifest.json`
  fetch is a 200 here and a 304 against the bucket.

- `rehearsal`: the deploy rehearsal of Plan 6. A local preview at `http://localhost:4173` serving a
  production bundle (`npm run build` with `VITE_DATA_BASE` left at its `.env.production` value) that
  reads the real bucket, with `VITE_MANIFEST` pointing at the staged manifest copy named by
  `--manifest`. It is the reverse of the preview guard: the run stops if any data request goes to
  `localhost` or `127.0.0.1`, or if the app fetched a manifest other than the named key.

## Scenarios

- `cold`: new browser context (empty cache and storage), `page.goto(url)`. App bundle, DuckDB
  wasm from jsDelivr, the parquet extension, and data are all cold. Zero is navigation start.
- `warm`: the same context right after `cold`, `page.reload()`. HTTP caches are warm; the engine
  restarts. Zero is navigation start of the reloaded document.
- `nav`: new context, load `/genes` and wait for table rows (the engine and `search_index` are
  loaded), wait until the network has been quiet for 2 s, then navigate inside the app with
  `history.pushState` plus a synthetic `popstate`. Zero is `performance.now()` just before
  `pushState`. This is one gene's own data cost with the engine running. If React Router ignores
  the popstate, the harness falls back to typing the symbol into the genes filter and clicking the
  row (`--nav-mode click` forces it).
- `pair`: exactly two `--pages` ids, for example `--scenario pair --pages MYOZ1-eqtl,SYNPO2L-sqtl`.
  A new context loads the first page cold and waits for ready and idle, all unrecorded, then the
  navigation to the second page is measured with the same popstate (or click) code as `nav`. It is
  recorded under the page id `MYOZ1-eqtl>SYNPO2L-sqtl`. Both of those genes are on chr10, so the
  second page sends new ranges to pack URLs the browser has just cached, which is the pattern behind
  the 20 s Chrome cache-lock stall logged on 2026-09-04. Watch `longest data req ms`. `pair` runs on
  its own, not alongside `cold`, `warm` or `nav`.

Order per run: for each page, `cold` then `warm` in one context, then `nav` in a fresh context.
Headless, 1280x900, no throttling.

## Readiness contract

The page is ready when the locus scatter's host has lost `invisible` and holds drawn dots:

```js
document.querySelector('.plot-host:not(.invisible) svg g[aria-label="dot"] > *')
```

polled every animation frame (timeout 120 s). A `pages.json` entry with `ready` uses that selector instead, and the run
metadata records each page's selector (`ready_selectors`). Ready time is `performance.now()` read inside the
poll at the frame it matched. LocusCompare also uses `.plot-host` but mounts only after the scatter
is ready, so the selector cannot fire early. Any change to the gene page must keep this selector
meaning "first paint of the scatter with data", or update it in the same commit.

After ready, recording continues until no request is in flight and none started or ended for 2 s,
capped at 30 s (`to_idle`).

## What is recorded

Per request (raw log): start and end relative to the scenario's zero, URL, method, resource type,
`Range` header, status, `content-length`, `content-range`, `request.sizes()`, `request.timing()`,
`fromServiceWorker`, and a category: `data` (the data host or `/data/`), `app` (target origin),
`duckdb-cdn` (`cdn.jsdelivr.net`), `duckdb-ext` (`extensions.duckdb.org`), `other`. `blob:` and
`data:` URLs are logged as `local` and never counted.

- Bytes are transferred bytes: `responseBodySize + responseHeadersSize` from `request.sizes()`.
- From cache: status 200 or 206 with zero transferred bytes (Chromium reports a cache hit as a
  body size of minus the header size), or served by a service worker. The raw fields are kept, so
  the rule can be refined without re-running.
- Max concurrency: most `data` requests in flight at once, from the request intervals. The app
  reads no parquet, so the old DuckDB probe-plus-`HEAD` pairs are gone; what overlaps now is the
  app's own parallel reads, such as the startup fetches (`manifest.json`, the search index, the
  GWAS index and the variant index) and a gene page's block, variants and GWAS requests, which go
  out in one tick.
- Longest data request: `data_max_ms` is the largest `end_ms - start_ms` over the `data` requests
  of the window, and `data_max_url` is that request's URL. A stall shows up here and nowhere else in
  the summary, so it is the number the `pair` scenario exists to watch. The markdown prints the ms
  as `longest data req ms` in the first table and the file name as `slowest data request` in the
  second; `compare.mjs` prints both windows' ms and lists the two file names.
- Kind, for `data` requests, from the file name. Every file the browser reads at a byte offset lives
  flat under `immutable/` as `<stem>.<sha16>.<ext>` (SPEC section 3), so the kind is the stem:
  `eqtl_pack` (`immutable/eqtl.`), `variants` (`immutable/variants.`), `sqtl_pack`
  (`immutable/sqtl.`), `gwas_pack` (`immutable/gwas.`), `trans_pack` (`immutable/trans.`), `hits`,
  `gwas_index`, `rsid_index`, `variant_index`, `search_index`, `manifest` (`manifest.json` or a
  staged `immutable/manifest.<sha16>.json`), and `other` (anything else). The parquet kinds
  `gene_detail`, `eqtl_nominal`, `sqtl_nominal`, `gwas` (`gwas_dcm/`) and `trans` (`trans_pairs/`)
  are kept under their old names so older results stay comparable; every one of them should now read
  0. Results recorded before Plan 5 counted the manifest and search index under `other`, and results
  before Plan 6 matched the packs under `packs/<kind>/` names that no longer exist.

Per scenario, from the page after idle: the gene page marks `gene:hit`, `gene:detail`, and
`locus:drawn`, and the `pack:*` measures summed over the scenario: the `pack:block`,
`pack:variants`, `pack:gwas`, and `pack:intron` fetches, `pack:decode`, `pack:worker-wait` (the
locus insert waiting for the DuckDB worker to finish other work, such as the trans read), and
`pack:insert` (the GWAS window and the locus rows), and for the trans tables (Plan 8)
`pack:trans` (the frame fetch), `pack:trans-decode` (frame decode and Arrow build), and the trans
table's worker wait and insert (`pack:worker-wait` and `pack:insert` with part `trans`, reported
as `pack_trans_wait_ms` and `pack_trans_insert_ms` and left out of the locus totals). Results before Plan 5 also hold `data_pack`,
from an attribute the page no longer has.

Summary per page and scenario, as median, min, and max over runs, for `to_ready` (requests started
before ready) and `to_idle` (everything): request counts (`data`, all; for `data`: GET, HEAD, 206,
200, from cache), bytes (`data`, all, per category), max `data` concurrency, the longest single
`data` request and its URL (`data_max_ms`, `data_max_url`), requests and bytes per kind, whether the eQTL block and variants requests overlapped in time (`pack_overlap`), whether the
block, variants, and GWAS pack requests were all in flight at one moment (`pack_overlap_gwas`), the
ms from the end of the eQTL block response to the start of the first intron block request
(`intron_after_block_ms`); plus ms to ready, ms to
idle (last response end), ms from the first `data` request to the end of the last one started
before ready (clipped at ready), and the marks as ms from the scenario's zero (`gene:hit` to
`gene:detail` and to `locus:drawn` as well). The markdown adds a second table per scenario with the
kinds and marks.

Self-check, first on every invocation: a cold FLNC load must record at least one `data` request
with a `Range` header and status 206 (today these are the DuckDB worker's synchronous XHRs). If it
does not, the harness stops with "worker requests not captured"; the fallback is `puppeteer-core`
driving the Playwright-installed Chromium, enabling CDP `Network` on each worker target from
`page.on('workercreated')`.

Metadata: time, target, label, runs, Playwright and Chromium versions, the app bundle name from the
target's `index.html`, the data `manifest.json` `built` and `pipeline_commit` as served, repo HEAD
(and whether the tree was dirty), hostname, CPU count, network interface and `--note`, and the
1-minute load average at the start and end of every scenario. Timings on a busy machine are
indicative; counts and bytes are the stable part.

## Outputs

- `out/<label>-<target>-<iso>.json`: raw per-request logs (gitignored).
- `results/<label>-<target>.json` and `.md`: the summary (committed). The markdown has one table
  per scenario, one row per page and window.
- `compare.mjs <before.json> <after.json>` prints a markdown table per scenario (`cold`, `warm`,
  `nav`, `pair`) with each metric
  as `before -> after (change %)` for page ids present in both, with `–` for a metric a run did
  not record (rows neither run recorded are left out). Paste it into the plan log.
