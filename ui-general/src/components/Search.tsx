import { useEffect, useRef, useState } from 'react'
import { useNavigate } from 'react-router'
import { ArrowRight } from 'lucide-react'
import { searchGenes, type SearchHit } from '@/lib/queries'

const RS = /^rs\d+$/i
const POS = /^(chr)?([0-9]+|x|y)[:\s]([0-9,]+)$/i
const REGION = /^(chr)?([0-9]+|x|y)[:\s]([0-9,]+)[-–]([0-9,]+)$/i

/** Route a free-text query: rsID, chr:pos, chr:start-end, else null (gene typeahead). */
export function routeFor(q: string): string | null {
  const s = q.trim()
  if (!s) return null
  if (RS.test(s)) return `/variant/${s.toLowerCase()}`
  let m = REGION.exec(s)
  if (m) return `/region/chr${m[2].toUpperCase()}:${m[3].replace(/,/g, '')}-${m[4].replace(/,/g, '')}`
  m = POS.exec(s)
  if (m) return `/variant/chr${m[2].toUpperCase()}:${m[3].replace(/,/g, '')}`
  return null
}

/** Search bar (atlas SearchBox): bordered bar, no leading icon, circular submit; `hero`
 *  uses a primary-tinted border. Gene typeahead drops below the bar. */
export default function Search({ hero = false, autoFocus = false }: { hero?: boolean; autoFocus?: boolean }) {
  const [q, setQ] = useState('')
  const [hits, setHits] = useState<SearchHit[]>([])
  const [open, setOpen] = useState(false)
  const [active, setActive] = useState(0)
  const nav = useNavigate()
  const box = useRef<HTMLDivElement>(null)

  useEffect(() => {
    const s = q.trim()
    if (s.length < 2 || routeFor(s)) { setHits([]); return }
    let alive = true
    const t = setTimeout(() => searchGenes(s).then(h => { if (alive) { setHits(h); setActive(0) } }), 80)
    return () => { alive = false; clearTimeout(t) }
  }, [q])

  useEffect(() => {
    const onDoc = (e: MouseEvent) => { if (!box.current?.contains(e.target as Node)) setOpen(false) }
    document.addEventListener('mousedown', onDoc)
    return () => document.removeEventListener('mousedown', onDoc)
  }, [])

  function go(hit?: SearchHit) {
    const target = hit ? `/gene/${hit.gene_id}` : routeFor(q) ?? (hits[0] ? `/gene/${hits[0].gene_id}` : null)
    if (!target) return
    setOpen(false); setQ('')
    nav(target)
  }

  const empty = !q.trim()
  return (
    <div ref={box} className="relative w-full">
      <div className={`flex items-center gap-3 rounded-lg border-base-300 transition-colors focus-within:border-base-content/50 ${
        hero ? 'border-[1.5px] px-3.5 py-2.5' : 'border bg-base-100 px-3 py-1.5'}`}>
        <input
          className="flex-1 bg-transparent text-sm outline-none placeholder:text-base-content/50"
          placeholder="Gene symbol, Ensembl ID, rsID, chr:pos, or chr:start-end"
          value={q}
          autoFocus={autoFocus}
          onChange={e => { setQ(e.target.value); setOpen(true) }}
          onFocus={() => setOpen(true)}
          onKeyDown={e => {
            if (e.key === 'ArrowDown') { e.preventDefault(); setActive(a => Math.min(a + 1, hits.length - 1)) }
            else if (e.key === 'ArrowUp') { e.preventDefault(); setActive(a => Math.max(a - 1, 0)) }
            else if (e.key === 'Enter') { e.preventDefault(); go(hits[active]) }
            else if (e.key === 'Escape') setOpen(false)
          }}
        />
        <button type="button" aria-label="Search" onClick={() => go(hits[active])} disabled={empty}
          className={`flex shrink-0 cursor-pointer items-center justify-center rounded-full bg-primary text-primary-content transition-colors hover:bg-primary/90 disabled:cursor-not-allowed disabled:opacity-40 ${hero ? 'size-7' : 'size-6'}`}>
          <ArrowRight className="size-3.5" />
        </button>
      </div>
      {open && hits.length > 0 && (
        <ul className="absolute z-30 mt-2 max-h-80 w-full overflow-auto rounded-lg border border-base-300 bg-base-100 p-1 shadow-lg">
          {hits.map((h, i) => (
            <li key={h.gene_id}>
              <button type="button" onMouseEnter={() => setActive(i)} onClick={() => go(h)}
                className={`flex w-full items-center gap-2 rounded-md px-2 py-1.5 text-left text-sm transition-colors ${i === active ? 'bg-base-200' : 'hover:bg-base-200'}`}>
                <span className="font-medium">{h.symbol ?? h.gene_id}</span>
                <span className="min-w-0 truncate text-xs text-base-content/55">{h.gene_id} · {h.chr}:{h.tss.toLocaleString()}</span>
                <span className="ml-auto flex shrink-0 gap-1">
                  {h.is_egene && <span className="badge badge-primary badge-xs">eGene</span>}
                  {h.n_sqtl_sig > 0 && <span className="badge badge-secondary badge-xs">sQTL</span>}
                  {!h.tested && <span className="badge badge-ghost badge-xs">not tested</span>}
                </span>
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}
