import { useEffect, useState } from 'react'
import { useParams } from 'react-router'
import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { KvTable } from '@/components/kv-table'
import { SectionPanel } from '@/components/section-panel'
import { DetailSkeleton, Empty, TableSkeleton } from '@/components/states'
import TransTable from '@/components/TransTable'
import { dbsnp, ucsc } from '@/lib/links'
import { fmtBp, fmtInt, fmtNum, fmtP, fmtPhenotype, fmtSlopeSE } from '@/lib/format'
import { dropTable, materialize } from '@/lib/db'
import { cisHitsAt, credibleSetsAt, leadGenesAt, leadPhenotypesAt, transAtSQL, variantByPosition, variantByRsid,
  type CisHit, type CredibleSetHit, type Gene, type SplicePhenotype, type VariantRow } from '@/lib/queries'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

const gnomad = (v: VariantRow) => `https://gnomad.broadinstitute.org/variant/${v.chr.replace('chr', '')}-${v.position}-${v.A2}-${v.A1}?dataset=gnomad_r4`
const ensemblVar = (rsid: string) => `https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${rsid}`
const openTargets = (v: VariantRow) => `https://platform.opentargets.org/variant/${v.chr.replace('chr', '')}_${v.position}_${v.A2}_${v.A1}`
/** Shown in place of the three cis sections for a variant seen only in the genome-wide trans scan. */
const OUTSIDE_CIS = 'This variant is more than 1 Mb from every tested gene, so it is outside every cis window.'

export default function Variant() {
  const { id = '' } = useParams()
  const [vars, setVars] = useState<VariantRow[] | null | undefined>(undefined)

  useEffect(() => {
    setVars(undefined)
    const rs = /^rs(\d+)$/i.exec(id)
    const pos = /^(chr[0-9XY]+):(\d+)$/i.exec(id)
    const q = rs ? variantByRsid(Number(rs[1])) : pos ? variantByPosition(pos[1], Number(pos[2])) : Promise.resolve([])
    q.then(v => setVars(v.length ? v : null))
  }, [id])

  if (vars === undefined) return <Page><DetailSkeleton kvRows={4} /></Page>
  const v = vars?.[0]
  return (
    <Page>
      <PageHeader crumbs={[{ label: 'Variants' }, { label: id }]} title={v?.rsid ?? id}
        meta={v ? <span className="tabular-nums">{v.chr}:{fmtInt(v.position)}</span> : undefined} />
      {!v ? (
        <div className="space-y-2">
          <Empty label={`${id} is not among the variants tested in TOPCHeF (MAF ≥ 0.01; cis windows within 1 Mb of a tested gene, or the genome-wide trans scan).`} />
          {/^rs\d+$/i.test(id) && <p className="px-4 text-sm"><ExternalLink icon href={dbsnp(id.toLowerCase())}>Look it up in dbSNP</ExternalLink></p>}
        </div>
      ) : <VariantBody vars={vars!} />}
    </Page>
  )
}

function VariantBody({ vars }: { vars: VariantRow[] }) {
  const v = vars[0]
  const rowLink = useRowLink()
  const [leads, setLeads] = useState<{ genes: Gene[]; phens: SplicePhenotype[]; cs: CredibleSetHit[] } | null>(null)
  // the variant's trans rows as an in-memory table, materialized once per variant and dropped
  // when it changes; the trans table pages off it like the gene page's does
  const [transTable, setTransTable] = useState<string | null>(null)
  const [scan, setScan] = useState<{ e: CisHit[]; s: CisHit[] } | null | 'running'>(null)

  useEffect(() => {
    let alive = true
    let table: string | null = null
    setLeads(null); setTransTable(null); setScan(null)
    Promise.all([leadGenesAt(v.chr, v.position), leadPhenotypesAt(v.chr, v.position), credibleSetsAt(v.chr, v.position)])
      .then(([genes, phens, cs]) => { if (alive) setLeads({ genes, phens, cs }) })
    materialize(transAtSQL(v.chr, v.position), 'trans')
      .then(t => { if (!alive) { dropTable(t); return } table = t; setTransTable(t) })
      .catch(e => console.error(e))
    return () => { alive = false; if (table) dropTable(table) }
  }, [v.chr, v.position])

  async function runScan() {
    setScan('running')
    const [e, s] = await Promise.all([cisHitsAt(v.chr, v.position, 'e'), cisHitsAt(v.chr, v.position, 's')])
    setScan({ e, s })
  }

  return (
    <div className="space-y-8">
      <div className="grid items-start gap-4 md:grid-cols-2">
        <KvTable rows={[
          { label: 'rsID', value: v.rsid ? <ExternalLink icon href={dbsnp(v.rsid)}>{v.rsid}</ExternalLink> : '—' },
          { label: 'Position', value: <span className="tabular-nums">{v.chr}:{fmtInt(v.position)} (GRCh38)</span> },
          { label: vars.length > 1 ? 'Alleles (A1 / A2)' : 'A1 / A2',
            value: v.A1 ? vars.map(x => `${x.A1} / ${x.A2}`).join(', ') : 'not reported (the trans eQTL file has no allele columns)' },
          { label: 'rsID match', value: v.match === 'exact' ? 'alleles match dbSNP'
            : v.match === 'position' ? (v.A1 ? 'position only (alleles differ from dbSNP record)' : 'by position (alleles not reported)')
            : 'no dbSNP record' },
        ]} />
        <KvTable rows={[
          { label: 'Links', value: <span className="flex flex-wrap gap-x-4">
            <ExternalLink icon href={ucsc(v.chr, v.position - 50, v.position + 50)}>UCSC</ExternalLink>
            {/* both URLs embed ref and alt, so they need the alleles */}
            {v.A1 && <ExternalLink icon href={gnomad(v)}>gnomAD</ExternalLink>}
            {v.A1 && <ExternalLink icon href={openTargets(v)}>Open Targets</ExternalLink>}
            {v.rsid && <ExternalLink icon href={ensemblVar(v.rsid)}>Ensembl</ExternalLink>}
          </span> },
          { label: 'A1', value: 'effect allele (minor allele in TOPCHeF)' },
          { label: 'A2', value: 'reference allele' },
        ]} />
      </div>

      <SectionPanel title="Lead variant for" description="Genes and splice phenotypes where this is the top cis association.">
        {!v.in_cis ? <Empty label={OUTSIDE_CIS} /> : leads === null ? <TableSkeleton columns={[{ w: 'w-10' }, { w: 'w-16' }, { w: 'w-40' }, { w: 'w-20', align: 'right' }, { w: 'w-14', align: 'right' }, { w: 'w-12', align: 'right' }]} rows={2} /> :
          leads.genes.length + leads.phens.length === 0 ? <Empty label="Not the lead variant for any gene or splice phenotype." /> : (
            <div className="overflow-x-auto rounded-lg border border-base-300">
              <table className="table table-sm">
                <thead><tr><th>Type</th><th>Gene</th><th>Phenotype</th><th className="text-right">Slope ± SE</th><th className="text-right">Perm p</th><th className="text-right">Status</th></tr></thead>
                <tbody>
                  {leads.genes.map(g => (
                    <tr key={g.gene_id} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${g.gene_id}`)}>
                      <td><span className="badge badge-primary badge-xs">eQTL</span></td>
                      <td className="font-medium"><span className={ROW_LINK_TEXT}>{g.symbol ?? g.gene_id}</span></td>
                      <td className="tabular-nums text-base-content/60">{g.gene_id}</td>
                      <td className="text-right tabular-nums">{fmtSlopeSE(g.slope, g.slope_se)}</td>
                      <td className="text-right tabular-nums">{fmtP(g.pval_perm)}</td>
                      <td className="text-right">{g.is_egene ? <span className="badge badge-primary badge-xs">eGene</span> : ''}</td>
                    </tr>
                  ))}
                  {leads.phens.map(p => (
                    <tr key={p.phenotype_id} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${p.gene_id}?tab=sqtl`)}>
                      <td><span className="badge badge-secondary badge-xs">sQTL</span></td>
                      <td className="font-medium"><span className={ROW_LINK_TEXT}>{p.symbol ?? p.gene_id}</span></td>
                      <td className="tabular-nums text-base-content/60">{fmtPhenotype(p.phenotype_id)}</td>
                      <td className="text-right tabular-nums">{fmtSlopeSE(p.slope, p.slope_se)}</td>
                      <td className="text-right tabular-nums">{fmtP(p.pval_perm)}</td>
                      <td className="text-right">{p.is_sqtl ? <span className="badge badge-secondary badge-xs">sQTL</span> : ''}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
      </SectionPanel>

      <SectionPanel title="Credible-set membership" description="SuSiE 95% credible sets containing this variant.">
        {!v.in_cis ? <Empty label={OUTSIDE_CIS} /> : leads === null ? <TableSkeleton columns={[{ w: 'w-10' }, { w: 'w-16' }, { w: 'w-40' }, { w: 'w-6' }, { w: 'w-10', align: 'right' }]} rows={2} /> :
          leads.cs.length === 0 ? <Empty label="Not in any credible set." /> : (
            <div className="overflow-x-auto rounded-lg border border-base-300">
              <table className="table table-sm">
                <thead><tr><th>Type</th><th>Gene</th><th>Phenotype</th><th>Set</th><th className="text-right">PIP</th></tr></thead>
                <tbody>
                  {leads.cs.map((c, i) => (
                    <tr key={i} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${c.gene_id}${c.qtl_type === 's' ? '?tab=sqtl' : ''}`)}>
                      <td><span className={`badge badge-xs ${c.qtl_type === 'e' ? 'badge-primary' : 'badge-secondary'}`}>{c.qtl_type === 'e' ? 'eQTL' : 'sQTL'}</span></td>
                      <td className="font-medium"><span className={ROW_LINK_TEXT}>{c.symbol ?? c.gene_id}</span></td>
                      <td className="tabular-nums text-base-content/60">{c.qtl_type === 'e' ? c.gene_id : fmtPhenotype(c.phenotype_id)}</td>
                      <td><span className="badge badge-ghost badge-sm">{c.cs_id}</span></td>
                      <td className="text-right tabular-nums font-medium">{fmtNum(c.pip)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
      </SectionPanel>

      <SectionPanel title="trans associations" description="Genes and splice phenotypes anywhere in the genome whose expression or splicing this variant associates with, outside their cis windows.">
        <TransTable table={transTable} keyedBy="variant" fileStem={`${v.rsid ?? `${v.chr}_${v.position}`}_trans`} />
      </SectionPanel>

      <SectionPanel title="All cis associations" description="Nominal statistics for every gene and splice phenotype whose window covers this variant. Reads every gene window overlapping the position, so it runs on request."
        action={v.in_cis && scan === null && <button className="btn btn-sm h-8 rounded-lg border-base-300 font-medium" onClick={runScan}>Scan cis windows</button>}>
        {!v.in_cis ? <Empty label={OUTSIDE_CIS} /> : scan === null ? <Empty label="Not scanned yet." /> : scan === 'running' ? <TableSkeleton columns={[{ w: 'w-16' }, { w: 'w-14', align: 'right' }, { w: 'w-10', align: 'right' }, { w: 'w-14', align: 'right' }, { w: 'w-20', align: 'right' }, { w: 'w-16', align: 'right' }]} rows={6} /> : (
          <div className="space-y-4">
            <HitTable title={`Expression (${scan.e.length})`} hits={scan.e} qtlType="e" />
            <HitTable title={`Splicing (${scan.s.length})`} hits={scan.s} qtlType="s" />
          </div>
        )}
      </SectionPanel>
    </div>
  )
}

function HitTable({ title, hits, qtlType }: { title: string; hits: CisHit[]; qtlType: 'e' | 's' }) {
  const rowLink = useRowLink()
  return (
    <div className="space-y-2">
      <h3 className="text-sm font-medium">{title}</h3>
      {hits.length === 0 ? <Empty label="No windows cover this variant." /> : (
        <div className="overflow-x-auto rounded-lg border border-base-300">
          <table className="table table-sm">
            <thead><tr><th>Gene</th>{qtlType === 's' && <th>Phenotype</th>}<th className="text-right">TSS dist</th><th className="text-right">AF</th><th className="text-right">p</th><th className="text-right">Slope ± SE</th><th className="text-right">PIP</th></tr></thead>
            <tbody>
              {hits.map((h, i) => (
                <tr key={i} className={`${ROW_LINK} hover:bg-base-200 ${h.pip != null ? 'bg-base-200/70' : ''}`} {...rowLink(`/gene/${h.gene_id}${qtlType === 's' ? '?tab=sqtl' : ''}`)}>
                  <td className="font-medium"><span className={ROW_LINK_TEXT}>{h.symbol ?? h.gene_id}</span></td>
                  {qtlType === 's' && <td className="tabular-nums text-base-content/60">{fmtPhenotype(h.phenotype_id ?? '')}</td>}
                  <td className="text-right tabular-nums text-base-content/60">{fmtBp(h.tss_distance)}</td>
                  <td className="text-right tabular-nums text-base-content/60">{fmtNum(h.af)}</td>
                  <td className="text-right tabular-nums">{fmtP(h.pval_nominal)}</td>
                  <td className="text-right tabular-nums">{fmtSlopeSE(h.slope, h.slope_se)}</td>
                  <td className="text-right tabular-nums">{h.pip != null ? `${fmtNum(h.pip)} (set ${h.cs_id})` : ''}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}
