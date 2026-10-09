import { lazy, Suspense, useEffect, useState } from 'react'
import { Route, Routes, useNavigate } from 'react-router'
import { Toaster } from 'sonner'
import { TopNav } from '@/components/top-nav'
import { Page } from '@/components/page'
import { ErrorBoundary, ErrorPage } from '@/components/error-page'
import Home from '@/routes/Home'
import Genes from '@/routes/Genes'
import Gene from '@/routes/Gene'
import Variant from '@/routes/Variant'
import Region from '@/routes/Region'
import About from '@/routes/About'
import { onDBFailure } from '@/lib/db'
import { useThemedFavicon } from '@/lib/favicon'
import { FAULT } from '@/lib/faults'
import { useStoreError } from '@/contexts/store-context'
import { useTheme } from '@/contexts/theme-context'

// dev builds only: the fault-trigger page (lib/faults.ts); dead code in a production bundle
const DevFaults = import.meta.env.DEV ? lazy(() => import('@/routes/DevFaults')) : null
/** The `render` fault: a component that throws, to see the error boundary. */
function RenderFault(): never { throw new Error('fault: a component threw while rendering') }

function NotFound() {
  const navigate = useNavigate()
  return <ErrorPage title="There is no such page." action={{ label: 'Home', onClick: () => navigate('/') }} />
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
                <Route path="/" element={<Home />} />
                <Route path="/genes" element={<Genes />} />
                <Route path="/gene/:id" element={<Gene />} />
                <Route path="/variant/:id" element={<Variant />} />
                <Route path="/region/:loc" element={<Region />} />
                <Route path="/about" element={<About />} />
                <Route path="*" element={<Page><NotFound /></Page>} />
              </Routes>
            )}
        </ErrorBoundary>
      </main>
    </div>
  )
}
