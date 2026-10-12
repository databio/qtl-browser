# ui — the general multi-study browser

A fork of the TOPCHeF app, now `../ui-topchef/` (2026-09-25, `plans/2026-09-25-ui-fork-topchef.md`),
with the same stack, store and decoders. `../ui-topchef/` is the site at https://topchef.databio.org.
This one deploys to its own Worker, `qtl-browser` (`wrangler.jsonc`), on every push to main that
touches `ui/` (`.github/workflows/deploy-ui.yml`). It has no custom domain yet.

It browses every collection in the store (`plans/2026-10-10-general-ui-first-pass.md`):

| Path | Page |
|---|---|
| `/` | the collections as rows, studies and GWAS apart; nothing opens until one is chosen |
| `/about` | the browser and the store |
| `/c/<digest>` | a collection: its publication, facets and experiments, each a link |
| `/c/<digest>/e/<id>` | an experiment: search and Top loci (a cis scan), search (a trans scan), Top loci (a GWAS) |
| `/c/<digest>/e/<id>/genes`, `/gene/<g>`, `/variant/<v>`, `/region/<loc>` | that QTL experiment's pages |
| `/c/<digest>/e/<id>/gwas/<gwas id>/...` | the same, with that GWAS for coloc and LocusCompare, chosen on a gene page |

`<digest>` is the collection object's digest. Every page is one experiment's: the URL's prefix up
to the GWAS is the router's basename (`src/lib/scope.ts`), so a page load holds one experiment and
one GWAS, and choosing another is a full navigation. No GWAS is chosen until the visitor picks
one. A GWAS on another variant catalog is matched by position and alleles. The code is not shared
with `../ui-topchef/`: a fix to both is made in both.

## Running it

```bash
npm install
QTL_DATA_DIR=../data/derived/store-release-v2 VITE_DATA_BASE= npm run dev   # a local store at /data, with Range support
VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser npm run dev           # the live store
```

`QTL_DATA_DIR` is read when the dev server starts, so pointing it at another store needs a
restart. Open the dev fault page inside an experiment, `/c/<digest>/e/<id>/dev/faults`.

`npx tsc -b` type-checks. `npm run store-check` validates a local store against the Python
decoders (`bench/README.md`). The smoke suites open each of TOPCHeF's four experiments, and GTEx's eQTL experiment against the DCM GWAS.
