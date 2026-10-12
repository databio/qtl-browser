/**
 * Fault injection for looking at every error surface, dev builds only. `?fault=<name>` on any URL
 * is read once at load and holds until the next full reload; `&only=<ext>` narrows an object fault
 * to objects of that extension (`qbe` results, `qbv` variants, `qbh` hits, `qbt` trans, `qbg`
 * GWAS, `qgi` GWAS index, `qbr` rsID index, `qbx` variant index, `arrow.zst` index parts and
 * tables, `json` collection documents). routes/DevFaults.tsx (`/dev/faults`) lists a trigger per
 * surface. A production build has `FAULT === 'none'` as a constant, so every hook folds away.
 */
export const FAULTS = {
  none: 'no fault',
  offline: 'every fetch rejects the way a lost connection does',
  slow: 'every fetch waits 2.5 s before answering (skeletons and spinners)',
  'store-404': 'store.json answers 404',
  'store-version': 'store.json claims format version 1',
  'store-anchor': 'store.json lists no experiment of the configured collection',
  'store-names': 'an experiment document names an object that is not content-addressed',
  'object-404': 'immutable objects answer 404 (narrow with only=<ext>)',
  'object-corrupt': 'immutable objects answer bytes of the right length that decode to nothing (only=<ext>)',
  'range-200': 'range requests answer 200 with the whole object, as a host without Range support does',
  db: 'the DuckDB-WASM engine fails to start',
  render: 'a component throws while rendering (no error boundary catches it)',
} as const
export type Fault = keyof typeof FAULTS

function read(): { fault: Fault; only: string | null } {
  const q = new URLSearchParams(window.location.search)
  const f = q.get('fault')
  return { fault: f && f in FAULTS ? f as Fault : 'none', only: q.get('only') }
}
const active = import.meta.env.DEV ? read() : { fault: 'none' as Fault, only: null }

export const FAULT: Fault = active.fault
/** The extension an object fault is narrowed to, or null for every object. */
export const FAULT_ONLY: string | null = active.only

/** Whether an object fault applies to this object name. */
export const faultHits = (name: string) => FAULT_ONLY === null || name.endsWith(`.${FAULT_ONLY}`)

export const SLOW_MS = 2500
export const delay = (ms: number) => new Promise<void>(r => setTimeout(r, ms))

/** The fetch-level faults, before a request leaves: a rejected promise or a wait. */
export async function beforeFetch(): Promise<void> {
  if (FAULT === 'offline') throw new TypeError('Failed to fetch')
  if (FAULT === 'slow') await delay(SLOW_MS)
}
