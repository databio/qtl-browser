import { useEffect, useState } from 'react'
import { useParams } from 'react-router'
import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { Empty, TableSkeleton } from '@/components/states'
import { ucsc } from '@/lib/links'
import { fmtInt } from '@/lib/format'
import { genesInRegion, type SearchHit } from '@/lib/queries'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

export default function Region() {
  const { loc = '' } = useParams()
  const m = /^(chr[0-9XY]+):(\d+)-(\d+)$/i.exec(loc)
  const [genes, setGenes] = useState<SearchHit[] | null>(null)
  const rowLink = useRowLink()
  useEffect(() => {
    if (!m) return
    setGenes(null)
    genesInRegion(m[1], Number(m[2]), Number(m[3])).then(setGenes)
  }, [loc]) // eslint-disable-line react-hooks/exhaustive-deps
  if (!m) return <Page><Empty label="Region must look like chr10:73000000-74500000." /></Page>
  const [, chr, s, e] = m
  return (
    <Page>
      <PageHeader crumbs={[{ label: 'Regions' }, { label: loc }]} title={<span className="tabular-nums">{chr}:{fmtInt(s)}-{fmtInt(e)}</span>} meta={`${fmtInt(Number(e) - Number(s))} bp`}
        description={<ExternalLink icon href={ucsc(chr, Number(s), Number(e))}>Open in UCSC</ExternalLink>} />
      {genes === null ? <TableSkeleton columns={[{ w: 'w-20' }, { w: 'w-32' }, { w: 'w-24', align: 'right' }, { w: 'w-20' }]} /> :
        genes.length === 0 ? <Empty label="No genes with a TSS in this region." /> : (
          <div className="overflow-x-auto rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead><tr><th>Gene</th><th>Ensembl ID</th><th className="text-right">TSS</th><th>Status</th></tr></thead>
              <tbody>
                {genes.map(g => (
                  <tr key={g.gene_id} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${g.gene_id}`)}>
                    <td className="font-medium"><span className={ROW_LINK_TEXT}>{g.symbol ?? g.gene_id}</span></td>
                    <td className="text-base-content/60">{g.gene_id}</td>
                    <td className="text-right tabular-nums">{fmtInt(g.tss)}</td>
                    <td className="space-x-1">
                      {g.is_egene && <span className="badge badge-primary badge-xs">eGene</span>}
                      {g.n_sqtl_sig > 0 && <span className="badge badge-secondary badge-xs">{g.n_sqtl_sig} sQTL</span>}
                      {!g.tested && <span className="badge badge-ghost badge-xs">not tested</span>}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
    </Page>
  )
}
