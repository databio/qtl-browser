import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'
import { EQTL_TYPE, getStoreInfo, type StoreInfo } from '@/lib/store'
import { significanceOf, type Significance } from '@/lib/significance'

const Ctx = createContext<StoreInfo | null>(null)
const ErrorCtx = createContext<unknown>(null)

/** Opens the store once at the app root (store.json, the experiment, catalog and annotation
 *  pointers, and the search index for the counts; the fetches are shared with the readers, see
 *  lib/store.ts); every page reads it through `useStoreInfo`. Without it no page can fetch
 *  anything, so a failure is the app's page-level error (App.tsx reads `useStoreError`). */
export function StoreProvider({ children }: { children: ReactNode }) {
  const [m, setM] = useState<StoreInfo | null>(null)
  const [error, setError] = useState<unknown>(null)
  useEffect(() => {
    let alive = true
    getStoreInfo()
      .then(x => { if (alive) setM(x) })
      .catch((e: unknown) => { console.error(e); if (alive) setError(e) })
    return () => { alive = false }
  }, [])
  return <Ctx.Provider value={m}><ErrorCtx.Provider value={error}>{children}</ErrorCtx.Provider></Ctx.Provider>
}

/** Null until the store has opened. */
export const useStoreInfo = () => useContext(Ctx)

/** The error the store failed to open with, or null. */
export const useStoreError = () => useContext(ErrorCtx)

/** How the study decides significance, for every label, tooltip and badge about it
 *  (lib/significance.ts): the eQTL cis scan's rule, or the first cis scan's when there is none.
 *  `assessed` is false until the store opens, so nothing renders a rule it has not read yet. */
export const useSignificance = (): Significance => {
  const m = useStoreInfo()
  return significanceOf(m?.cis.get(EQTL_TYPE) ?? m?.cis.values().next().value)
}
