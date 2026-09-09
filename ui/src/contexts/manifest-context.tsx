import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'

/** Shape of data/derived/manifest.json as the pipeline's manifest step writes it. */
export interface Manifest {
  built: string
  pipeline_commit: string | null
  significance_rule: string
  counts: Record<string, number>
  sources: Record<string, { version: string; description: string }>
  tables: Record<string, { path: string; rows: number; bytes: number; files: number; columns: string[]; load: string }>
  assets?: Record<string, { path: string; bytes: number }>
  gwas_dcm: { file: string; n_cases: number; n_controls: number; variants: number } | null
}

const URL = `${(import.meta.env.VITE_DATA_BASE as string | undefined) ?? '/data'}/manifest.json`
const Ctx = createContext<Manifest | null>(null)

/** Fetches manifest.json once at the app root; every page reads it through `useManifest`. */
export function ManifestProvider({ children }: { children: ReactNode }) {
  const [m, setM] = useState<Manifest | null>(null)
  useEffect(() => {
    let alive = true
    fetch(URL)
      .then(r => { if (!r.ok) throw new Error(`manifest.json: HTTP ${r.status}`); return r.json() as Promise<Manifest> })
      .then(x => { if (alive) setM(x) })
      .catch(e => console.error(e))
    return () => { alive = false }
  }, [])
  return <Ctx.Provider value={m}>{children}</Ctx.Provider>
}

/** Null until the manifest has loaded. */
export const useManifest = () => useContext(Ctx)
