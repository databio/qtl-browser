import { useEffect, useMemo, useState } from 'react'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { Segmented } from '@/components/segmented'
import { SortableTh, type SortState } from '@/components/sortable-th'
import { Pager } from '@/components/pager'
import { Empty, TableSkeleton } from '@/components/states'
import { ErrorPage } from '@/components/error-page'
import { fmtInt } from '@/lib/format'
import { useCollectionCrumbs, useSignificance, useStoreInfo } from '@/contexts/store-context'
import { EQTL_TYPE, SQTL_TYPE } from '@/lib/store'
import { allHits } from '@/lib/gene-index'
import type { SearchHit } from '@/lib/queries'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

type Filter = 'egenes' | 'sqtl' | 'tested' | 'all'
const CHR_ORDER = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY', 'chrM']

/** Browse every annotated gene (every chromosome's genes and search index part, read once): filter, sort, page. */
export default function Genes() {
  const sig = useSignificance()
  const home = useCollectionCrumbs()
  const m = useStoreInfo()
  const [all, setAll] = useState<SearchHit[] | null>(null)
  // the significance filters only exist for an experiment that assessed it, so the default moves to
  // `tested` rather than opening on a filter that can match nothing (lib/significance.ts)
  // and to the experiment's own scan: an eQTL scan opens on its eGenes, an sQTL scan on its sQTL
  // genes, a trans scan (no cis results) on every gene
  const hasE = m?.cis.has(EQTL_TYPE) ?? false, hasS = m?.cis.has(SQTL_TYPE) ?? false
  const [filter, setFilter] = useState<Filter>('egenes')
  useEffect(() => {
    if (!m) return
    setFilter(!hasE && !hasS ? 'all' : !sig.assessed ? 'tested' : hasE ? 'egenes' : 'sqtl')
  }, [m, hasE, hasS, sig.assessed])
  const [needle, setNeedle] = useState('')
  const [sort, setSort] = useState<SortState>({ by: 'symbol', order: 'asc' })
  const [offset, setOffset] = useState(0)
  const [pageSize, setPageSize] = useState(25)
  const rowLink = useRowLink()

  const [error, setError] = useState<unknown>(null)
  useEffect(() => { allHits().then(setAll).catch(e => { console.error(e); setError(e) }) }, [])
  useEffect(() => setOffset(0), [filter, needle, sort])

  const shown = useMemo(() => {
    if (!all) return []
    const n = needle.trim().toUpperCase()
    let out = all.filter(g =>
      (filter === 'all' || (filter === 'tested' && (hasE ? g.tested : g.has_results)) || (filter === 'egenes' && g.is_egene) || (filter === 'sqtl' && g.n_sqtl_sig > 0)) &&
      (!n || (g.symbol ?? '').toUpperCase().startsWith(n) || g.gene_id.startsWith(n)))
    const dir = sort.order === 'asc' ? 1 : -1
    if (sort.by === 'symbol') out = out.sort((a, b) => dir * (a.symbol ?? a.gene_id).localeCompare(b.symbol ?? b.gene_id))
    else if (sort.by === 'position') out = out.sort((a, b) => dir * ((CHR_ORDER.indexOf(a.chr) - CHR_ORDER.indexOf(b.chr)) || (a.tss - b.tss)))
    else if (sort.by === 'n_sqtl_sig') out = out.sort((a, b) => dir * (a.n_sqtl_sig - b.n_sqtl_sig))
    return out
  }, [all, filter, needle, sort, hasE])

  const page = shown.slice(offset, offset + pageSize)
  if (error) return <Page><ErrorPage title="The gene list could not be loaded." error={error} /></Page>
  return (
    <Page>
      <PageHeader crumbs={[...home, { label: 'Genes' }]} title="Genes" meta={all ? fmtInt(shown.length) : undefined}
        description={`All ${m?.annotation.id ?? ''} genes; tested genes carry this experiment's QTL results.`} />
      <div className="space-y-2">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <Segmented value={filter} onChange={setFilter} options={[
            ...(sig.assessed && hasE ? [{ value: 'egenes' as Filter, label: 'eGenes' }] : []),
            ...(sig.assessed && hasS ? [{ value: 'sqtl' as Filter, label: 'sQTL genes' }] : []),
            ...(hasE || hasS ? [{ value: 'tested' as Filter, label: 'Tested' }] : []), { value: 'all', label: 'All' }]} />
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
