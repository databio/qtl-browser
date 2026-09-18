import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'
import { getManifest, type Manifest } from '@/lib/manifest'

const Ctx = createContext<Manifest | null>(null)

/** Loads the manifest once at the app root (the fetch is shared with the pack reader, see
 *  lib/manifest.ts); every page reads it through `useManifest`. Without it no page can fetch
 *  anything, so a failure is shown above the app rather than only logged. */
export function ManifestProvider({ children }: { children: ReactNode }) {
  const [m, setM] = useState<Manifest | null>(null)
  const [error, setError] = useState<string | null>(null)
  useEffect(() => {
    let alive = true
    getManifest()
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

/** Null until the manifest has loaded. */
export const useManifest = () => useContext(Ctx)
