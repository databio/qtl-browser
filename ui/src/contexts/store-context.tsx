import { createContext, useContext, useEffect, useState, type ReactNode } from 'react'
import { EQTL_TYPE, getIndex, getStoreInfo, type ExperimentDoc, type IndexCollection, type StoreIndex, type StoreInfo } from '@/lib/store'
import { experimentLabel, shortName } from '@/lib/describe'
import { collectionPath, SCOPE } from '@/lib/scope'
import { significanceOf, type Significance } from '@/lib/significance'

const Ctx = createContext<StoreInfo | null>(null)
const IndexCtx = createContext<StoreIndex | null>(null)
const ErrorCtx = createContext<unknown>(null)

/** Opens the store at the app root. Every page gets the index (`store.json` and the collection,
 *  experiment and coloc documents); a page inside a QTL experiment also gets its study, which
 *  every reader needs (lib/store.ts). A failure is the app's page-level error (App.tsx reads
 *  `useStoreError`). */
export function StoreProvider({ children }: { children: ReactNode }) {
  const [m, setM] = useState<StoreInfo | null>(null)
  const [index, setIndex] = useState<StoreIndex | null>(null)
  const [error, setError] = useState<unknown>(null)
  useEffect(() => {
    let alive = true
    const fail = (e: unknown) => { console.error(e); if (alive) setError(e) }
    // inside an experiment the study starts opening at once, beside the index, so a gene page's
    // reads do not wait for every collection document; whether it is a QTL experiment of this
    // collection is known when the index lands, and an open that failed because it was not is dropped
    const study = SCOPE?.experiment ? getStoreInfo() : null
    study?.catch(() => {})
    getIndex().then(ix => {
      if (!alive) return
      setIndex(ix)
      const e = experimentIn(ix)
      if (study && e && e.kind !== 'gwas') study.then(x => { if (alive) setM(x) }, fail)
    }, fail)
    return () => { alive = false }
  }, [])
  return (
    <IndexCtx.Provider value={index}>
      <Ctx.Provider value={m}><ErrorCtx.Provider value={error}>{children}</ErrorCtx.Provider></Ctx.Provider>
    </IndexCtx.Provider>
  )
}

/** Null until the study has opened, and always null outside a QTL collection. */
export const useStoreInfo = () => useContext(Ctx)

/** Null until the index has loaded. */
export const useIndex = () => useContext(IndexCtx)

/** The page's collection from the index, or undefined (no collection, or not loaded yet). */
export function useCollection(): IndexCollection | undefined {
  const ix = useIndex()
  return SCOPE ? ix?.collections.get(SCOPE.collection) : undefined
}

/** The URL's experiment, when the index holds it in the URL's collection. */
const experimentIn = (ix: StoreIndex | null): ExperimentDoc | undefined =>
  SCOPE?.experiment ? ix?.collections.get(SCOPE.collection)?.experiments.find(e => e.id === SCOPE!.experiment) : undefined

/** The page's experiment from the index, or undefined (none in the URL, not in this collection, or
 *  not loaded yet). */
export const useExperiment = (): ExperimentDoc | undefined => experimentIn(useIndex())

/** The crumbs every page inside an experiment starts with: the landing page and the collection's
 *  page (both outside the experiment's basename, so full navigations), then the experiment's own. */
export function useCollectionCrumbs(): { to?: string; href?: string; label: string }[] {
  const c = useCollection(), e = useExperiment()
  return [{ href: '/', label: 'Collections' },
    { href: SCOPE ? collectionPath(SCOPE.digest) : '/', label: c ? shortName(c.doc) : '…' },
    { to: '/', label: e ? experimentLabel(e) : '…' }]
}

/** The error the store failed to open with, or null. */
export const useStoreError = () => useContext(ErrorCtx)

/** How the study decides significance, for every label, tooltip and badge about it
 *  (lib/significance.ts): the eQTL cis scan's rule, or the first cis scan's when there is none.
 *  `assessed` is false until the store opens, so nothing renders a rule it has not read yet. */
export const useSignificance = (): Significance => {
  const m = useStoreInfo()
  return significanceOf(m?.cis.get(EQTL_TYPE) ?? m?.cis.values().next().value)
}
