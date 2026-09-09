import { useEffect, useState } from 'react'
import { Download, Search as SearchIcon } from 'lucide-react'
import { SortableTh, type SortState } from '@/components/sortable-th'
import { Pager } from '@/components/pager'
import { Empty, TableSkeleton } from '@/components/states'
import { fmtInt, fmtNum, fmtP, fmtPhenotype, fmtSlopeSE } from '@/lib/format'
import { transAll, transCount, transRows, type TransQuery, type TransRow } from '@/lib/queries'
import { downloadCSV } from '@/lib/csv'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

/** Sortable, filterable page through one gene's trans rows of one QTL type, off the trans
 *  table the gene page materialized. Each change is one local query plus a count, like
 *  CisTable. `table` null means the gene's rows are still loading. */
export default function TransTable({ table, qtlType, fileStem }: { table: string | null; qtlType: 'e' | 's'; fileStem: string }) {
  const rowLink = useRowLink()
  const variantPath = (r: TransRow) => `/variant/${r.rsid ?? `${r.variant_chr}:${r.position}`}`
  const [sort, setSort] = useState<SortState>({ by: 'pval', order: 'asc' })
  const [maxP, setMaxP] = useState('')
  const [search, setSearch] = useState('')
  const [offset, setOffset] = useState(0)
  const [pageSize, setPageSize] = useState(10)
  const [data, setData] = useState<TransRow[] | null>(null)
  const [total, setTotal] = useState(0)
  const [all, setAll] = useState<number | null>(null)     // unfiltered count, to tell "none" from "none match"
  const [busy, setBusy] = useState(false)

  useEffect(() => { setOffset(0); setData(null); setAll(null) }, [table])
  useEffect(() => setOffset(0), [sort, maxP, search, pageSize])

  const query = (): TransQuery => {
    const p = maxP === '' ? undefined : Number(maxP)
    return { table: table!, qtlType, maxP: Number.isFinite(p) ? p : undefined, search: search || undefined,
      orderBy: sort.by || 'pval', desc: sort.by ? sort.order === 'desc' : false, limit: pageSize, offset }
  }

  useEffect(() => {
    if (!table) return
    let alive = true
    setBusy(true)
    const q = query()
    const t = setTimeout(() => {
      Promise.all([transRows(q), transCount(q), all ?? transCount({ table, qtlType })])
        .then(([r, n, a]) => { if (alive) { setData(r); setTotal(n); setAll(a) } })
        .catch(() => {})   // the table was dropped under us (gene changed): the next effect run replaces it
        .finally(() => { if (alive) setBusy(false) })
    }, search ? 150 : 0)
    return () => { alive = false; clearTimeout(t) }
  }, [table, sort, maxP, search, offset, pageSize]) // eslint-disable-line react-hooks/exhaustive-deps

  async function exportCSV() {
    const cols = ['variant_chr', 'position', 'rsid', 'af', 'pval', 'beta', 'beta_se', 'r2']
    downloadCSV(`${fileStem}.csv`, await transAll(query()), qtlType === 's' ? ['phenotype_id', ...cols] : cols)
  }

  const skel = [...(qtlType === 's' ? [{ w: 'w-40' }] : []), { w: 'w-24' }, { w: 'w-20' }, { w: 'w-10', align: 'right' as const },
    { w: 'w-14', align: 'right' as const }, { w: 'w-20', align: 'right' as const }, { w: 'w-10', align: 'right' as const }]

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-center gap-2">
        <label className="input input-bordered input-sm flex h-8 w-56 items-center gap-2 rounded-lg">
          <SearchIcon className="size-4 shrink-0 opacity-50" />
          <input type="search" className="grow bg-transparent outline-none" placeholder="rsID or position" value={search} onChange={e => setSearch(e.target.value)} />
        </label>
        {busy && data && <span className="loading loading-spinner loading-xs text-base-content/40" />}
        <div className="flex-1" />
        <select className="select select-bordered select-sm h-8 rounded-lg" value={maxP} onChange={e => setMaxP(e.target.value)} title="p-value threshold">
          <option value="">All p</option>
          <option value="1e-6">p &lt; 1e-6</option>
          <option value="1e-8">p &lt; 1e-8</option>
          <option value="1e-10">p &lt; 1e-10</option>
        </select>
        <button className="btn btn-sm h-8 gap-1.5 rounded-lg border-base-300 font-medium" onClick={exportCSV} disabled={!table || total === 0}><Download className="size-3.5" /> CSV</button>
      </div>
      {data === null ? <TableSkeleton columns={skel} rows={3} /> : all === 0 ? <Empty label="No trans associations." /> : total === 0 ? <Empty label="No variants match." /> : (
        <>
          <div className="overflow-x-auto rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead>
                <tr>
                  {qtlType === 's' && <th>Intron</th>}
                  <SortableTh sortKey="position" label="Variant" sort={sort} onSort={setSort} defaultOrder="asc" />
                  <th>rsID</th>
                  <SortableTh sortKey="af" label="AF" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <SortableTh sortKey="pval" label="p" sort={sort} onSort={setSort} defaultOrder="asc" className="text-right" align="right" />
                  <SortableTh sortKey="beta" label="Beta ± SE" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <SortableTh sortKey="r2" label="r²" sort={sort} onSort={setSort} className="text-right" align="right" />
                </tr>
              </thead>
              <tbody>
                {data.map((r, i) => (
                  <tr key={i} className={`${ROW_LINK} hover:bg-base-200/60`} {...rowLink(variantPath(r))}>
                    {qtlType === 's' && <td className="tabular-nums text-base-content/60">{fmtPhenotype(r.phenotype_id)}</td>}
                    <td className="tabular-nums">{r.variant_chr}:{fmtInt(r.position)}</td>
                    <td><span className={`${ROW_LINK_TEXT} ${r.rsid ? '' : 'text-base-content/40'}`}>{r.rsid ?? `${r.variant_chr}:${r.position}`}</span></td>
                    <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.af)}</td>
                    <td className="text-right tabular-nums">{fmtP(r.pval)}</td>
                    <td className="text-right tabular-nums">{fmtSlopeSE(r.beta, r.beta_se)}</td>
                    <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.r2)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <Pager total={total} offset={offset} pageSize={pageSize} onPage={setOffset} onPageSize={setPageSize} />
        </>
      )}
    </div>
  )
}
