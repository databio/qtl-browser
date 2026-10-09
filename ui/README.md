# ui — the general multi-study browser

A fork of the TOPCHeF app, now `../ui-topchef/` (2026-09-25, `plans/2026-09-25-ui-fork-topchef.md`),
with the same stack, store and decoders. `../ui-topchef/` is the site at https://topchef.databio.org.
This one deploys to its own Worker, `qtl-browser` (`wrangler.jsonc`), on every push to main that
touches `ui/` (`.github/workflows/deploy-ui.yml`). It has no custom domain yet.

Today `src/` is the TOPCHeF app's code, version 2 reader included, so it opens one collection like
the TOPCHeF app does. It stops being a copy with the redesign: navigation by collection and
experiment with no default entry point, gene and variant pages that start from where a gene or
variant is present and expand per experiment, coloc and LocusCompare against a GWAS the visitor
picks, and none of what only applies to TOPCHeF. The design is
`plans/2026-10-09-general-ui-brainstorm.md`. Until then a change belonging in both apps is applied
to `../ui-topchef/` and copied; `diff -rq src ../ui-topchef/src` shows what differs.

## Running it

```bash
npm install
QTL_DATA_DIR=../data/derived/store-release-v2 VITE_DATA_BASE= npm run dev   # a local store at /data, with Range support
VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser npm run dev           # the live store
```

`QTL_DATA_DIR` is read when the dev server starts, so pointing it at another store needs a
restart. `VITE_COLLECTION_OF` names the experiment whose collection opens (default
`topchef_ge_cis`) and `VITE_GWAS` the GWAS beside it (default `dcm_jurgens2024_biobanks`, empty for
none). The reader refuses a GWAS on another variant catalog, so GTEx opens with
`VITE_COLLECTION_OF=gtex_v8_heart_lv_ge_cis VITE_GWAS=`.

`npx tsc -b` type-checks. `npm run store-check` validates a local store against the Python
decoders (`bench/README.md`). The bench and smoke suites describe the TOPCHeF pages.
