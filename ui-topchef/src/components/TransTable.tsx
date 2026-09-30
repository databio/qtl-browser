import { useEffect, useState } from 'react'
import { Download, Search as SearchIcon } from 'lucide-react'
import { SortableTh, type SortState } from '@/components/sortable-th'
import { Pager } from '@/components/pager'
import { Empty, TableSkeleton } from '@/components/states'
import { fmtInt, fmtNum, fmtP, fmtPhenotype, fmtSlopeSE } from '@/lib/format'
import { transAll, transCount, transRows, type TransQuery, type TransRow } from '@/lib/queries'
import { downloadCSV, roundedCsvName } from '@/lib/csv'
import { transRoundingDetail, useTransRoundingFacts } from '@/lib/rounding'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

/** Sortable, filterable page through a materialized trans table. Each change is one local
 *  query plus a count, like CisTable. `table` null means the rows are still loading.
 *
 *  Keyed by gene (the gene page): rows are the variants of one QTL type, linking to the
 *  variant. Keyed by variant (the variant page): rows are the genes and introns of both
 *  types, linking to the gene, with a type badge and the gene's location instead of the
 *  variant's. */
export default function TransTable({ table, qtlType, keyedBy = 'gene', fileStem }: {
  table: string | null; qtlType?: 'e' | 's'; keyedBy?: 'gene' | 'variant'; fileStem: string
}) {
  const byVariant = keyedBy === 'variant'
  const rowLink = useRowLink()
  const rowPath = (r: TransRow) => byVariant
    ? `/gene/${r.gene_id}${r.qtl_type === 's' ? '?tab=sqtl' : ''}`
    : `/variant/${r.rsid ?? `${r.variant_chr}:${r.position}`}`
  const [sort, setSort] = useState<SortState>({ by: 'pval', order: 'asc' })
  const rounding = useTransRoundingFacts()
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
    return { table: table!, qtlType, keyedBy, maxP: Number.isFinite(p) ? p : undefined, search: search || undefined,
      orderBy: sort.by || 'pval', desc: sort.by ? sort.order === 'desc' : false, limit: pageSize, offset }
  }

  useEffect(() => {
    if (!table) return
    let alive = true
    setBusy(true)
    const q = query()
    const t = setTimeout(() => {
      Promise.all([transRows(q), transCount(q), all ?? transCount({ table, qtlType, keyedBy })])
        .then(([r, n, a]) => { if (alive) { setData(r); setTotal(n); setAll(a) } })
        .catch(() => {})   // the table was dropped under us (gene changed): the next effect run replaces it
        .finally(() => { if (alive) setBusy(false) })
    }, search ? 150 : 0)
    return () => { alive = false; clearTimeout(t) }
  }, [table, sort, maxP, search, offset, pageSize]) // eslint-disable-line react-hooks/exhaustive-deps

  async function exportCSV() {
    const stats = ['pval', 'beta', 'beta_se', 'r2']
    const cols = byVariant
      ? ['qtl_type', 'gene_id', 'symbol', 'phenotype_id', 'gene_chr', 'gene_tss', ...stats]
      : [...(qtlType === 's' ? ['phenotype_id'] : []), 'variant_chr', 'position', 'rsid', 'af', ...stats]
    downloadCSV(roundedCsvName(fileStem), await transAll(query()), cols)
  }

  const R = { align: 'right' as const }
  const skel = byVariant
    ? [{ w: 'w-10' }, { w: 'w-16' }, { w: 'w-40' }, { w: 'w-24' }, { w: 'w-14', ...R }, { w: 'w-20', ...R }, { w: 'w-10', ...R }]
    : [...(qtlType === 's' ? [{ w: 'w-40' }] : []), { w: 'w-24' }, { w: 'w-20' }, { w: 'w-10', ...R }, { w: 'w-14', ...R }, { w: 'w-20', ...R }, { w: 'w-10', ...R }]

  return (
    <div className="space-y-2" data-trans-total={data !== null ? (all ?? undefined) : undefined}>
      <div className="flex flex-wrap items-center gap-2">
        <label className="input input-bordered input-sm flex h-8 w-56 items-center gap-2 rounded-lg">
          <SearchIcon className="size-4 shrink-0 opacity-50" />
          <input type="search" className="grow bg-transparent outline-none" placeholder={byVariant ? 'Gene symbol or Ensembl ID' : 'rsID or position'} value={search} onChange={e => setSearch(e.target.value)} />
        </label>
        {busy && data && <span className="loading loading-spinner loading-xs text-base-content/40" />}
        <div className="flex-1" />
        {/* TensorQTL only wrote trans pairs with p < 1e-5, so the unfiltered option is labeled with that floor */}
        <select className="select select-bordered select-sm h-8 rounded-lg" value={maxP} onChange={e => setMaxP(e.target.value)} title="p-value threshold">
          <option value="">p ≤ 1e-5</option>
          <option value="1e-6">p ≤ 1e-6</option>
          <option value="1e-8">p ≤ 1e-8</option>
          <option value="1e-10">p ≤ 1e-10</option>
        </select>
        {/* as on the cis table: the note is in the section description, the title rides with the file */}
        <button className="btn btn-sm h-8 gap-1.5 rounded-lg border-base-300 font-medium" title={rounding ? transRoundingDetail(rounding) : undefined}
          onClick={exportCSV} disabled={!table || total === 0}><Download className="size-3.5" /> CSV</button>
      </div>
      {data === null ? <TableSkeleton columns={skel} rows={3} /> : all === 0 ? <Empty label="No trans associations." /> : total === 0 ? <Empty label={byVariant ? 'No genes match.' : 'No variants match.'} /> : (
        <>
          <div className="overflow-x-auto rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead>
                <tr>
                  {byVariant ? <>
                    <SortableTh sortKey="type" label="Type" sort={sort} onSort={setSort} defaultOrder="asc" />
                    <th>Gene</th>
                    <th>Phenotype</th>
                    <SortableTh sortKey="gene" label="Gene location" sort={sort} onSort={setSort} defaultOrder="asc" />
                  </> : <>
                    {qtlType === 's' && <th>Intron</th>}
                    <SortableTh sortKey="position" label="Variant" sort={sort} onSort={setSort} defaultOrder="asc" />
                    <th>rsID</th>
                    <SortableTh sortKey="af" label="AF" sort={sort} onSort={setSort} className="text-right" align="right" />
                  </>}
                  <SortableTh sortKey="pval" label="p" sort={sort} onSort={setSort} defaultOrder="asc" className="text-right" align="right" />
                  <SortableTh sortKey="beta" label="Beta ± SE" sort={sort} onSort={setSort} className="text-right" align="right" />
                  <SortableTh sortKey="r2" label="r²" sort={sort} onSort={setSort} className="text-right" align="right" />
                </tr>
              </thead>
              <tbody>
                {data.map((r, i) => (
                  <tr key={i} className={`${ROW_LINK} hover:bg-base-200/60`} {...rowLink(rowPath(r))}>
                    {byVariant ? <>
                      <td><span className={`badge badge-xs ${r.qtl_type === 'e' ? 'badge-primary' : 'badge-secondary'}`}>{r.qtl_type === 'e' ? 'eQTL' : 'sQTL'}</span></td>
                      <td className="font-medium"><span className={ROW_LINK_TEXT}>{r.symbol ?? r.gene_id}</span></td>
                      <td className="tabular-nums text-base-content/60">{r.qtl_type === 'e' ? r.gene_id : fmtPhenotype(r.phenotype_id)}</td>
                      <td className="tabular-nums text-base-content/60">{r.gene_chr}:{fmtInt(r.gene_tss)}</td>
                    </> : <>
                      {qtlType === 's' && <td className="tabular-nums text-base-content/60">{fmtPhenotype(r.phenotype_id)}</td>}
                      <td className="tabular-nums">{r.variant_chr}:{fmtInt(r.position)}</td>
                      <td><span className={`${ROW_LINK_TEXT} ${r.rsid ? '' : 'text-base-content/40'}`}>{r.rsid ?? `${r.variant_chr}:${r.position}`}</span></td>
                      <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.af)}</td>
                    </>}
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
