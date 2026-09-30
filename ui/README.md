# ui — the general multi-study browser

Forked from `ui-topchef/` on 2026-09-25 (`plans/2026-09-25-ui-fork-topchef.md`). Same stack, same
store, same decoders. What changes here is what the app assumes about *which* study it is showing.

`ui-topchef/` is the deployed TOPCHeF site and stays working throughout. Both copies stay live —
features may land in either — so a fix that belongs in both is a decision each time. This directory is not
deployed and carries no `wrangler.jsonc`, so it cannot ship over that site by accident.

## Where it is going

The store already holds N experiments over M variant catalogs (SPEC.md sections 2 and 8). On B2
today: `topchef` (two phenotype types, trans, the DCM GWAS, GENCODE v34) and `gtex_v8_heart_lv`
(no trans, no GWAS, a fitted dof, GENCODE v39). Both catalogs anchor to the same seqcol collection,
so a site key `(seq_digest, pos, ref, alt)` means the same thing in both.

The browser is the only layer that assumes one study. The work, roughly in dependency order:

1. **Experiment in the route** (`/e/<id>/gene/<gene>`), not `VITE_EXPERIMENT` at build time.
   `store.json` already lists the experiments; nothing reads that list yet.
2. **Per-experiment caches.** `getStore`, `variantIndex`, `lookupDir`, `wholeGenes`, `wholeIndex`
   and the `memoBy(chr)` caches in `lib/store.ts` are singletons keyed to one store.
3. **Phenotype types from the data.** `EQTL_TYPE = 'ge'` and `SQTL_TYPE = 'leafcutter'` become
   whatever `results[].phenotype_type` holds.
4. **Study-specific content behind a config**: the coloc tables and DCM chips (`lib/coloc.ts`,
   `routes/Gene.tsx`, `components/ColocLoci.tsx`), the GWAS comparison labels
   (`components/LocusCompare.tsx`), the About and Home prose, the nav wordmark, and the
   `topchef`/`topchef-dark` themes.
5. **Cross-catalog lookup in TypeScript** (SPEC.md section 11). Python has it as
   `qtlstore crosscat`; every decoder the browser needs already exists, so this is assembly rather
   than new format work. It is what makes "this locus in the other study" reachable.

## Running it

```bash
npm install
QTL_DATA_DIR=/path/to/store VITE_DATA_BASE= npm run dev   # a local store at /data, with Range support
VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser npm run dev   # the live store
```

`VITE_EXPERIMENT` still selects the experiment (default `topchef`) until step 1 lands;
`VITE_EXPERIMENT=gtex_v8_heart_lv` is the second study to test against.

`npx tsc -b` type-checks. `npm run store-check` validates a local store against the Python
decoders. The bench suites under `bench/` came with the fork and still describe the TOPCHeF pages.
