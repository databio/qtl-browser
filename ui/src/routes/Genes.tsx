import { useEffect, useMemo, useState } from 'react'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { Segmented } from '@/components/segmented'
import { SortableTh, type SortState } from '@/components/sortable-th'
import { Pager } from '@/components/pager'
import { Empty, TableSkeleton } from '@/components/states'
import { fmtInt } from '@/lib/format'
import { rows } from '@/lib/db'
import type { SearchHit } from '@/lib/queries'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

type Filter = 'egenes' | 'sqtl' | 'tested' | 'all'
const CHR_ORDER = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY', 'chrM']

/** Browse every gene in the search index (in memory): filter, sort, page. */
export default function Genes() {
  const [all, setAll] = useState<SearchHit[] | null>(null)
  const [filter, setFilter] = useState<Filter>('egenes')
  const [needle, setNeedle] = useState('')
  const [sort, setSort] = useState<SortState>({ by: 'symbol', order: 'asc' })
  const [offset, setOffset] = useState(0)
  const [pageSize, setPageSize] = useState(25)
  const rowLink = useRowLink()

  useEffect(() => { rows<SearchHit>('SELECT * FROM search_index').then(setAll) }, [])
  useEffect(() => setOffset(0), [filter, needle, sort])

  const shown = useMemo(() => {
    if (!all) return []
    const n = needle.trim().toUpperCase()
    let out = all.filter(g =>
      (filter === 'all' || (filter === 'tested' && g.tested) || (filter === 'egenes' && g.is_egene) || (filter === 'sqtl' && g.n_sqtl_sig > 0)) &&
      (!n || (g.symbol ?? '').toUpperCase().startsWith(n) || g.gene_id.startsWith(n)))
    const dir = sort.order === 'asc' ? 1 : -1
    if (sort.by === 'symbol') out = out.sort((a, b) => dir * (a.symbol ?? a.gene_id).localeCompare(b.symbol ?? b.gene_id))
    else if (sort.by === 'position') out = out.sort((a, b) => dir * ((CHR_ORDER.indexOf(a.chr) - CHR_ORDER.indexOf(b.chr)) || (a.tss - b.tss)))
    else if (sort.by === 'n_sqtl_sig') out = out.sort((a, b) => dir * (a.n_sqtl_sig - b.n_sqtl_sig))
    return out
  }, [all, filter, needle, sort])

  const page = shown.slice(offset, offset + pageSize)
  return (
    <Page>
      <PageHeader title="Genes" meta={all ? fmtInt(shown.length) : undefined}
        description="All GENCODE v34 genes; tested genes carry eQTL and sQTL results." />
      <div className="space-y-2">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <Segmented value={filter} onChange={setFilter} options={[
            { value: 'egenes', label: 'eGenes' }, { value: 'sqtl', label: 'sQTL genes' },
            { value: 'tested', label: 'Tested' }, { value: 'all', label: 'All' }]} />
          <label className="input input-bordered input-sm flex h-8 w-64 items-center gap-2 rounded-lg">
            <input type="search" className="grow bg-transparent outline-none" placeholder="Filter by symbol or ID…" value={needle} onChange={e => setNeedle(e.target.value)} />
          </label>
        </div>
        {all === null ? <TableSkeleton columns={[{ w: 'w-20' }, { w: 'w-32' }, { w: 'w-28' }, { w: 'w-8', align: 'right' }, { w: 'w-14' }]} rows={12} /> :
          shown.length === 0 ? <Empty label="No genes match." /> : (
            <>
              <div className="overflow-x-auto rounded-lg border border-base-300">
                <table className="table table-sm">
                  <thead>
                    <tr>
                      <SortableTh sortKey="symbol" label="Gene" sort={sort} onSort={setSort} defaultOrder="asc" />
                      <th>Ensembl ID</th>
                      <SortableTh sortKey="position" label="TSS" sort={sort} onSort={setSort} defaultOrder="asc" />
                      <SortableTh sortKey="n_sqtl_sig" label="sQTL" sort={sort} onSort={setSort} className="text-right" align="right" />
                      <th>Status</th>
                    </tr>
                  </thead>
                  <tbody>
                    {page.map(g => (
                      <tr key={g.gene_id} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${g.gene_id}`)}>
                        <td className="font-medium"><span className={ROW_LINK_TEXT}>{g.symbol ?? g.gene_id}</span></td>
                        <td className="text-base-content/60">{g.gene_id}</td>
                        <td className="tabular-nums text-base-content/60">{g.chr}:{fmtInt(g.tss)}</td>
                        <td className="text-right tabular-nums">{g.n_sqtl_sig || ''}</td>
                        <td>
                          {g.is_egene && <span className="badge badge-primary badge-xs">eGene</span>}
                          {!g.tested && <span className="badge badge-ghost badge-xs">not tested</span>}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <Pager total={shown.length} offset={offset} pageSize={pageSize} onPage={setOffset} onPageSize={n => { setPageSize(n); setOffset(0) }} />
            </>
          )}
      </div>
    </Page>
  )
}
