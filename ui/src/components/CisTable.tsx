import { useEffect, useState, type CSSProperties } from 'react'
import { Download, Search as SearchIcon } from 'lucide-react'
import { SortableTh, type SortState } from '@/components/sortable-th'
import { Pager } from '@/components/pager'
import { Empty, TableSkeleton } from '@/components/states'
import { fmtBp, fmtNum, fmtP, rsFromNumber } from '@/lib/format'
import { cisAll, cisCount, cisRows, type CisQuery, type CisRow } from '@/lib/queries'
import { csTint, useIsDark } from '@/lib/plot-theme'
import { downloadCSV } from '@/lib/csv'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

const SKEL = [{ w: 'w-24' }, { w: 'w-20' }, { w: 'w-10' }, { w: 'w-14', align: 'right' as const }, { w: 'w-10', align: 'right' as const },
  { w: 'w-10', align: 'right' as const }, { w: 'w-14', align: 'right' as const }, { w: 'w-12', align: 'right' as const }, { w: 'w-12', align: 'right' as const }, { w: 'w-10', align: 'right' as const }]

/** Sortable, filterable page through one gene's (or intron's) cis window, off the table the
 *  locus plot materialized. Each change is one local query plus a count. `table` is null
 *  while the window loads; `failed` when the locus could not be materialized. */
export default function CisTable({ table, failed, chr, qtlType, phenotypeId, fileStem }: {
  table: string | null; failed?: boolean; chr: string; qtlType: 'e' | 's'; phenotypeId?: string; fileStem: string
}) {
  const [sort, setSort] = useState<SortState>({ by: 'pval_nominal', order: 'asc' })
  const [maxP, setMaxP] = useState('')
  const [search, setSearch] = useState('')
  const [offset, setOffset] = useState(0)
  const [pageSize, setPageSize] = useState(10)
  const [data, setData] = useState<CisRow[] | null>(null)
  const [total, setTotal] = useState(0)
  const [busy, setBusy] = useState(false)

  useEffect(() => { setOffset(0); setData(null) }, [table])
  useEffect(() => setOffset(0), [sort, maxP, search, pageSize])

  const query = (): CisQuery => {
    const p = maxP === '' ? undefined : Number(maxP)
    return { table: table!, chr, qtlType, phenotypeId, maxP: Number.isFinite(p) ? p : undefined, search: search || undefined,
      orderBy: sort.by || 'pval_nominal', desc: sort.by ? sort.order === 'desc' : false, limit: pageSize, offset }
  }

  useEffect(() => {
    if (!table) return
    let alive = true
    setBusy(true)
    const q = query()
    const t = setTimeout(() => {
      Promise.all([cisRows(q), cisCount(q)]).then(([r, n]) => { if (alive) { setData(r); setTotal(n) } })
        .catch(() => {})   // the table was dropped under us (locus changed): the next effect run replaces it
        .finally(() => { if (alive) setBusy(false) })
    }, search ? 150 : 0)
    return () => { alive = false; clearTimeout(t) }
  }, [table, sort, maxP, search, offset, pageSize]) // eslint-disable-line react-hooks/exhaustive-deps

  const dark = useIsDark()
  const rowLink = useRowLink()
  const variantPath = (r: CisRow) => `/variant/${r.rs_number != null ? rsFromNumber(r.rs_number) : `${chr}:${r.position}`}`

  async function exportCSV() {
    const all = await cisAll(query())
    const cols = ['position', 'rsid', 'A1', 'A2', 'tss_distance', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se', 'pip', 'cs_id']
    downloadCSV(`${fileStem}.csv`, all.map(r => ({ ...r, rsid: rsFromNumber(r.rs_number), phenotype_id: phenotypeId })), qtlType === 's' ? ['phenotype_id', ...cols] : cols)
  }

  if (failed) return <Empty label="Could not load the cis window." />
  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-center gap-2">
        <label className="input input-bordered input-sm flex h-8 w-56 items-center gap-2 rounded-lg">
          <SearchIcon className="size-4 shrink-0 opacity-50" />
          <input type="search" className="grow bg-transparent outline-none" placeholder="rsID or position" value={search} onChange={e => setSearch(e.target.value)} />
        </label>
        {busy && data && <span className="loading loading-spinner loading-xs text-base-content/40" />}
        <div className="flex-1" />
        <select className="select select-bordered select-sm h-8 rounded-lg" value={maxP} onChange={e => setMaxP(e.target.value)} title="Nominal p-value threshold">
          <option value="">p ≤ 1</option>
          <option value="1e-3">p ≤ 1e-3</option>
          <option value="1e-5">p ≤ 1e-5</option>
          <option value="5e-8">p ≤ 5e-8</option>
        </select>
        <button className="btn btn-sm h-8 gap-1.5 rounded-lg border-base-300 font-medium" onClick={exportCSV} disabled={!table}><Download className="size-3.5" /> CSV</button>
      </div>
      {data === null ? <TableSkeleton columns={SKEL} rows={10} /> : total === 0 ? <Empty label="No variants match." /> : (
        <>
          <div className="overflow-x-auto rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead>
                <tr>
                  <SortableTh sortKey="position" label="Position" sort={sort} onSort={setSort} defaultOrder="asc" />
                  <th>rsID</th>
                  <th>A1/A2</th>
                  <SortableTh sortKey="tss_distance" label="TSS dist" sort={sort} onSort={setSort} defaultOrder="asc" className="text-right" align="right" />
                  <SortableTh sortKey="af" label="AF" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <SortableTh sortKey="ma_count" label="MA count" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <SortableTh sortKey="pval_nominal" label="p" sort={sort} onSort={setSort} defaultOrder="asc" className="text-right" align="right" />
                  <SortableTh sortKey="slope" label="Slope" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <th className="text-right">SE</th>
                  <SortableTh sortKey="pip" label="PIP" sort={sort} onSort={setSort} className="text-right" align="right" />
                </tr>
              </thead>
              <tbody>
                {data.map((r, i) => (
                  // inline style: credible-set rows take the set's plot color from the chart palette;
                  // the resting and hover tints are passed as variables so the hover state is a
                  // stronger tint rather than a filter (which does not show on a table row)
                  <tr key={i} className={`${ROW_LINK} ${r.cs_id != null ? 'bg-(--tint) hover:bg-(--tint-hover)' : 'hover:bg-base-200/60'}`}
                    style={{ '--tint': csTint(r.cs_id, dark), '--tint-hover': csTint(r.cs_id, dark, true) } as CSSProperties}
                    {...rowLink(variantPath(r))}>
                    <td className="tabular-nums">{r.position.toLocaleString()}</td>
                    <td>{r.rs_number != null ? <span className={ROW_LINK_TEXT}>{rsFromNumber(r.rs_number)}</span> : <span className={`text-base-content/40 ${ROW_LINK_TEXT}`}>{chr}:{r.position}</span>}</td>
                    <td className="font-mono text-xs text-base-content/60">{r.A1}/{r.A2}</td>
                    <td className="text-right tabular-nums text-base-content/60">{fmtBp(r.tss_distance)}</td>
                    <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.af)}</td>
                    <td className="text-right tabular-nums text-base-content/60">{r.ma_count}</td>
                    <td className="text-right tabular-nums">{fmtP(r.pval_nominal)}</td>
                    <td className="text-right tabular-nums">{fmtNum(r.slope)}</td>
                    <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.slope_se)}</td>
                    <td className="text-right tabular-nums">{r.pip != null ? fmtNum(r.pip) : ''}</td>
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
