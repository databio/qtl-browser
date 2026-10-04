# ui — the TOPCHeF browser

This is the app deployed at https://topchef.databio.org. Vite, React, TypeScript, Tailwind 4 and
DaisyUI 5, with DuckDB-WASM and Mosaic for the plots.

It reads a qtlstore over HTTP range requests and decodes it in the browser. `lib/store.ts` resolves
the pointer documents and builds every data URL; `lib/store-decode.ts` decodes the binary objects;
`../SPEC.md` is the byte layout they both implement. DuckDB holds only what a page materializes (a
locus window, a GWAS window, trans rows), so the gene and variant pages start it and the other
pages do not.

```bash
npm install
QTL_DATA_DIR=/path/to/store VITE_DATA_BASE= npm run dev   # a local store at /data, with Range support
npm run dev                                               # the live store, from .env.production
```

| | |
|---|---|
| `npx tsc -b` | type-check. `tsc -p .` is a no-op here and always passes |
| `npm run store-check` | validate a local store against the Python decoders |
| `npm run build` | `tsc -b && vite build`; `.env.production` points at the B2 store |
| `bench/` | browser smoke suites and a gene-page cost harness (`bench/README.md`) |

`VITE_EXPERIMENT` selects the experiment, default `topchef`. Deploying is `npx wrangler deploy` with
root directory `ui` (`wrangler.jsonc`); the repo README has the full steps.

`../ui-general/` is the multi-study fork of this app, and will take over the `ui/` name and the
deployment once it has caught up. Until then the two can drift apart without anything noticing;
`diff -rq src ../ui-general/src` shows what has changed.
