import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'
import { getStoreInfo, type StoreInfo } from '@/lib/store'

const Ctx = createContext<StoreInfo | null>(null)

/** Opens the store once at the app root (store.json, the experiment, catalog and annotation
 *  pointers, and the search index for the counts; the fetches are shared with the readers, see
 *  lib/store.ts); every page reads it through `useStoreInfo`. Without it no page can fetch
 *  anything, so a failure is shown above the app rather than only logged. */
export function StoreProvider({ children }: { children: ReactNode }) {
  const [m, setM] = useState<StoreInfo | null>(null)
  const [error, setError] = useState<string | null>(null)
  useEffect(() => {
    let alive = true
    getStoreInfo()
      .then(x => { if (alive) setM(x) })
      .catch((e: unknown) => {
        console.error(e)
        if (alive) setError(e instanceof Error ? e.message : String(e))
      })
    return () => { alive = false }
  }, [])
  return (
    <Ctx.Provider value={m}>
      {error && <div role="alert" className="alert alert-error">{error}</div>}
      {children}
    </Ctx.Provider>
  )
}

/** Null until the store has opened. */
export const useStoreInfo = () => useContext(Ctx)
