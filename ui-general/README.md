# ui-general — the general multi-study browser

A fork of `ui/` (2026-09-25, `plans/2026-09-25-ui-fork-topchef.md`) with the same stack, store and
decoders. The difference is what it assumes about which study it is showing.

`ui/` remains the deployed TOPCHeF site. This directory is not deployed and has no
`wrangler.jsonc`, so it cannot be shipped over that site by accident. It will take over the `ui/`
name and the deployment once it has caught up.

Features may land in either copy, so a change belonging in both has to be applied twice, and
nothing checks that it was. `diff -rq src ../ui/src` shows what differs.

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
