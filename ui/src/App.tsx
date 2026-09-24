import { useEffect, useState } from 'react'
import { Route, Routes } from 'react-router'
import { TopNav } from '@/components/top-nav'
import { Page } from '@/components/page'
import Home from '@/routes/Home'
import Genes from '@/routes/Genes'
import Gene from '@/routes/Gene'
import Variant from '@/routes/Variant'
import Region from '@/routes/Region'
import About from '@/routes/About'
import { onDBFailure } from '@/lib/db'
import { useThemedFavicon } from '@/lib/favicon'

/** The shell renders at once. DuckDB starts only on a page that needs it (db.ts); every query
 *  awaits the engine internally and every page shows its own skeleton meanwhile, so there is no
 *  global "starting" screen. Only a start failure needs a message. */
export default function App() {
  const [dbError, setDbError] = useState<string | null>(null)
  useThemedFavicon()

  useEffect(() => onDBFailure(e => setDbError(e.message)), [])

  return (
    <div className="min-h-screen bg-base-100 text-base-content">
      <TopNav />
      <main className="min-w-0">
        {dbError && <Page><div className="alert alert-error text-sm">The query engine could not start: {dbError}</div></Page>}
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/genes" element={<Genes />} />
          <Route path="/gene/:id" element={<Gene />} />
          <Route path="/variant/:id" element={<Variant />} />
          <Route path="/region/:loc" element={<Region />} />
          <Route path="/about" element={<About />} />
        </Routes>
      </main>
    </div>
  )
}
