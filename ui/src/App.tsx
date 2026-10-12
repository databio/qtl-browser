import { lazy, Suspense, useEffect, useState, type ReactNode } from 'react'
import { Route, Routes, useNavigate } from 'react-router'
import { Toaster } from 'sonner'
import { TopNav } from '@/components/top-nav'
import { Page } from '@/components/page'
import { Empty } from '@/components/states'
import { ErrorBoundary, ErrorPage } from '@/components/error-page'
import Landing from '@/routes/Landing'
import Collection, { NotInStore } from '@/routes/Collection'
import Experiment from '@/routes/Experiment'
import Genes from '@/routes/Genes'
import Gene from '@/routes/Gene'
import Variant from '@/routes/Variant'
import Region from '@/routes/Region'
import About from '@/routes/About'
import { onDBFailure } from '@/lib/db'
import { useThemedFavicon } from '@/lib/favicon'
import { FAULT } from '@/lib/faults'
import { useCollection, useExperiment, useIndex, useStoreError } from '@/contexts/store-context'
import { SCOPE } from '@/lib/scope'
import { useTheme } from '@/contexts/theme-context'

// dev builds only: the fault-trigger page (lib/faults.ts); dead code in a production bundle
const DevFaults = import.meta.env.DEV ? lazy(() => import('@/routes/DevFaults')) : null
/** The `render` fault: a component that throws, to see the error boundary. */
function RenderFault(): never { throw new Error('fault: a component threw while rendering') }

function NotFound() {
  const navigate = useNavigate()
  return <ErrorPage title="There is no such page." action={{ label: 'Home', onClick: () => navigate('/') }} />
}

/** A page that needs the experiment's study: only a QTL experiment has one. */
function Study({ children }: { children: ReactNode }) {
  const ix = useIndex(), c = useCollection(), e = useExperiment()
  if (!ix) return children
  if (!c) return <NotInStore />
  if (!e) return <NotInStore what="experiment" />
  if (e.kind === 'gwas') return <Page><Empty label="A GWAS has no gene, variant or region pages of its own; choose it on a QTL experiment's gene page." /></Page>
  return children
}

/** The shell renders at once. DuckDB starts only on a page that needs it (db.ts); every query
 *  awaits the engine internally and every page shows its own skeleton meanwhile, so there is no
 *  global "starting" screen. Two failures take the page over (components/error-page.tsx): the
 *  store not opening, and the engine not starting on a page that needed it. */
export default function App() {
  const [dbError, setDbError] = useState<Error | null>(null)
  const storeError = useStoreError()
  const { mode } = useTheme()
  useThemedFavicon()

  useEffect(() => onDBFailure(setDbError), [])

  const dark = mode === 'dark' || (mode === 'system' && window.matchMedia('(prefers-color-scheme: dark)').matches)
  return (
    <div className="min-h-screen bg-base-100 text-base-content">
      <TopNav />
      {/* section-level failures (lib/toast.ts): persistent until closed, one per title */}
      <Toaster position="top-center" theme={dark ? 'dark' : 'light'} richColors closeButton icons={{ error: null }} />
      <main className="min-w-0">
        <ErrorBoundary>
          {import.meta.env.DEV && FAULT === 'render' && <RenderFault />}
          {storeError ? <Page><ErrorPage title="The data could not be opened." error={storeError} /></Page>
            : dbError ? <Page><ErrorPage title="The query engine could not start." error={dbError} /></Page>
            : (
              <Routes>
                {DevFaults && <Route path="/dev/faults" element={<Suspense fallback={null}><DevFaults /></Suspense>} />}
                {/* inside an experiment or on a collection's page (lib/scope.ts: the router's basename),
                    or on the landing page */}
                {SCOPE?.experiment ? <>
                  <Route path="/" element={<Experiment />} />
                  <Route path="/genes" element={<Study><Genes /></Study>} />
                  <Route path="/gene/:id" element={<Study><Gene /></Study>} />
                  <Route path="/variant/:id" element={<Study><Variant /></Study>} />
                  <Route path="/region/:loc" element={<Study><Region /></Study>} />
                </> : SCOPE ? <Route path="/" element={<Collection />} /> : <>
                  <Route path="/" element={<Landing />} />
                  <Route path="/about" element={<About />} />
                </>}
                <Route path="*" element={<Page><NotFound /></Page>} />
              </Routes>
            )}
        </ErrorBoundary>
      </main>
    </div>
  )
}
