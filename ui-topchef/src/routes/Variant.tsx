import { useEffect, useState } from 'react'
import { useParams } from 'react-router'
import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { KvTable } from '@/components/kv-table'
import { SectionPanel } from '@/components/section-panel'
import RoundingNote from '@/components/rounding-note'
import { DetailSkeleton, Empty, TableSkeleton, Unavailable } from '@/components/states'
import { ErrorPage } from '@/components/error-page'
import { toastError } from '@/lib/toast'
import TransTable from '@/components/TransTable'
import { CopyButton } from '@/components/copy-button'
import { dbsnp, ucsc } from '@/lib/links'
import { fmtBp, fmtBytes, fmtInt, fmtNum, fmtP, fmtPhenotype, fmtBetaSE } from '@/lib/format'
import { planScan, runScan, type CisHit, type ScanPlan } from '@/lib/cis-scan'
import type { VariantRecord } from '@/lib/store-decode'
import { csValues, hitPhenotypes, leadValues, loadHits, lookupRsid, nominalsAt, variantAt, variantAtPosition,
  type HitPhenotype, type Hits, type HitSet } from '@/lib/variant'
import { SQTL_TYPE } from '@/lib/store'
import { variantTransTable } from '@/lib/trans'
import { dropTable } from '@/lib/db'
import { useSignificance, useStoreInfo } from '@/contexts/store-context'
import { ROW_LINK, ROW_LINK_TEXT, useRowLink } from '@/lib/row-link'

const gnomad = (v: VariantRecord) => `https://gnomad.broadinstitute.org/variant/${v.chr.replace('chr', '')}-${v.position}-${v.A2}-${v.A1}?dataset=gnomad_r4`
const ensemblVar = (rsid: string) => `https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${rsid}`
const openTargets = (v: VariantRecord) => `https://platform.opentargets.org/variant/${v.chr.replace('chr', '')}_${v.position}_${v.A2}_${v.A1}`
/** Shown in place of the three cis sections for a variant seen only in the genome-wide trans scan. */
const OUTSIDE_CIS = 'This variant is more than 1 Mb from every tested gene, so it is outside every cis window.'

const rsidOf = (v: VariantRecord) => (v.rsNumber ? `rs${v.rsNumber}` : null)

/** The variant and its hits records: a page request and the chromosome's hits file once the
 *  variant is located (lib/variant.ts). */
async function resolve(id: string): Promise<{ v: VariantRecord; hits: Hits } | null> {
  const rs = /^rs(\d+)$/i.exec(id)
  const pos = /^(chr[0-9XYM]+):(\d+)$/i.exec(id)
  if (rs) {
    const ref = await lookupRsid(Number(rs[1]))
    if (!ref) return null
    // both need only (chr, vidx), so the page and the frame go out together
    const hitsP = loadHits(ref.chr, ref.vidx)
    const v = await variantAt(ref.chr, ref.vidx)
    return { v, hits: await hitsP }
  }
  if (pos) {
    const v = await variantAtPosition(pos[1], Number(pos[2]))
    if (!v) return null
    return { v, hits: await loadHits(v.chr, v.vidx) }
  }
  return null
}

export default function Variant() {
  const { id = '' } = useParams()
  const [found, setFound] = useState<{ v: VariantRecord; hits: Hits } | null | undefined>(undefined)
  // the lookup itself failed (not "no such variant"): the page cannot say either way
  const [error, setError] = useState<unknown>(null)

  useEffect(() => {
    let alive = true
    setFound(undefined); setError(null)
    resolve(id).then(r => { if (alive) setFound(r) }, e => { console.error(e); if (alive) setError(e) })
    return () => { alive = false }
  }, [id])

  if (error) return <Page><ErrorPage title={`${id} could not be looked up.`} error={error} /></Page>
  if (found === undefined) return <Page><DetailSkeleton kvRows={4} /></Page>
  return (
    <Page>
      <PageHeader crumbs={[{ label: 'Variants' }, { label: id }]} title={found ? (rsidOf(found.v) ?? id) : id}
        meta={found ? <span className="tabular-nums">{found.v.chr}:{fmtInt(found.v.position)}</span> : undefined} />
      {!found ? (
        <div className="space-y-2">
          <Empty label={`${id} is not among the variants tested in TOPCHeF (MAF ≥ 0.01; cis windows within 1 Mb of a tested gene, or the genome-wide trans scan).`} />
          {/^rs\d+$/i.test(id) && <p className="px-4 text-sm"><ExternalLink icon href={dbsnp(id.toLowerCase())}>Look it up in dbSNP</ExternalLink></p>}
        </div>
      ) : <VariantBody v={found.v} hits={found.hits} />}
    </Page>
  )
}

/** A phenotype named by a hits record's `ord`, with its gene. */
interface GeneRow { gene_id: string; symbol: string | null }

/** A record of one cis scan's hits: `set` is the scan it came from and `row` its index in that
 *  scan's frame, since ordinals and records are per experiment (SPEC section 8). */
interface CsRow { set: HitSet; row: number; gene: GeneRow; qtlType: 'e' | 's'; phenotypeId: string | null }
/** A lead row also carries the variant's stored beta and SE in that phenotype's block. */
interface LeadRow extends CsRow { beta: number | null; betaSe: number | null }

/** The two list sections, built from each cis scan's hits records, the phenotypes they name (that
 *  scan's search index part for their chromosome, and the genes), and the lead rows' blocks for
 *  their beta and SE, nearby blocks read together. */
async function buildLists(v: VariantRecord, hits: Hits): Promise<{ leads: LeadRow[]; cs: CsRow[] }> {
  const leadRows: (CsRow & { p: HitPhenotype })[] = [], cs: CsRow[] = []
  for (const set of hits.sets) {
    if (set.experiment.kind !== 'qtl_cis') continue
    const wanted = [...set.leads, ...set.cs]
    if (!wanted.length) continue
    const byOrd = await hitPhenotypes(set, wanted.map(r => set.frame.ord[r]))
    const of = (r: number): CsRow & { p: HitPhenotype } => {
      const p = byOrd.get(set.frame.ord[r])
      if (!p) throw new Error(`experiment ${set.experiment.id}: the search index has no phenotype with ord ${set.frame.ord[r]}`)
      const qtlType: 'e' | 's' = p.phenotype_type === SQTL_TYPE ? 's' : 'e'
      return { set, row: r, p, gene: { gene_id: p.gene_id ?? p.phenotype_id, symbol: p.symbol }, qtlType, phenotypeId: qtlType === 's' ? p.phenotype_id : null }
    }
    leadRows.push(...set.leads.map(of))
    cs.push(...set.cs.map(of))
  }
  if (!leadRows.length && !cs.length) return { leads: [], cs: [] }
  const nominal = await nominalsAt(leadRows.map(x => x.p), v.vidx)
  const leads = leadRows.map(x => {
    const n = nominal.get(`${x.p.experiment}:${x.p.ord}`)
    return { ...x, beta: n?.beta ?? null, betaSe: n?.se ?? null }
  })
  return { leads, cs }
}

function VariantBody({ v, hits }: { v: VariantRecord; hits: Hits }) {
  // `signif`, not `sig`: `sig` below is this variant's significant introns
  const signif = useSignificance()
  const rowLink = useRowLink()
  const [lists, setLists] = useState<{ leads: LeadRow[]; cs: CsRow[] } | null>(null)
  // the lead and credible-set lists are the page: their failure is page-level
  const [listsError, setListsError] = useState<unknown>(null)
  // the variant's trans rows as an in-memory table, built from its hits records and dropped when
  // the variant changes; the trans table pages off it like the gene page's does
  const [transTable, setTransTable] = useState<string | null>(null)
  const hasTrans = useStoreInfo()?.hasTrans ?? false
  const [plan, setPlan] = useState<ScanPlan | null>(null)
  const [scan, setScan] = useState<{ e: CisHit[]; s: CisHit[] } | null | 'running'>(null)
  const [allIntrons, setAllIntrons] = useState(false)

  useEffect(() => {
    let alive = true
    let table: string | null = null
    setLists(null); setListsError(null); setTransTable(null); setScan(null); setPlan(null); setAllIntrons(false)
    // the trans table needs the query engine (a 7 MB download); it starts once the lists have their
    // data, so the lists do not share the connection with it
    // a failed list is the page's error; the trans table is not attempted after it
    buildLists(v, hits).then(l => { if (!alive) return null; setLists(l); return variantTransTable(v, hits) }, e => { if (alive) setListsError(e); return null })
      .then(t => { if (t == null) return; if (!alive) { dropTable(t); return } table = t; setTransTable(t) })
      .catch(e => { if (alive) toastError('The trans associations could not be loaded.', e) })
    if (v.inCis) planScan(v).then(p => { if (alive) setPlan(p) }, e => { if (alive) toastError('The cis scan could not be planned.', e) })
    return () => { alive = false; if (table) dropTable(table) }
  }, [v, hits])

  async function doScan() {
    if (!plan) return
    setScan('running')
    try {
      setScan(await runScan(plan, v))
    } catch (e) {
      toastError('The cis scan failed.', e)
      setScan(null)
    }
  }

  const rsid = rsidOf(v)
  if (listsError) return <ErrorPage title={`${rsid ?? `${v.chr}:${v.position}`}'s results could not be loaded.`} error={listsError} />
  // harness hooks (bench/README.md): the lists are "ready" once they have their data, or once the
  // outside-cis message stands in for them
  const listsReady = !v.inCis || lists !== null
  const sig = scan && scan !== 'running' ? scan.s.filter(h => h.is_sqtl) : []
  const shown = scan && scan !== 'running' ? (allIntrons ? scan.s : sig) : []
  return (
    <div className="space-y-8" {...(listsReady ? { 'data-variant-lists': '1' } : {})}>
      <div className="grid items-start gap-4 md:grid-cols-2">
        <KvTable rows={[
          { label: 'rsID', value: rsid ? <span className="flex items-center justify-between gap-2">
            <ExternalLink icon href={dbsnp(rsid)} title="Open in dbSNP">{rsid}</ExternalLink>
            <CopyButton text={rsid} className="-my-1 -mr-1" />
          </span> : '—' },
          { label: 'Position', value: <span className="tabular-nums">{v.chr}:{fmtInt(v.position)} (GRCh38)</span> },
          { label: 'Links', value: <span className="flex flex-wrap gap-x-4">
            <ExternalLink icon href={ucsc(v.chr, v.position - 50, v.position + 50)}>UCSC</ExternalLink>
            <ExternalLink icon href={gnomad(v)}>gnomAD</ExternalLink>
            <ExternalLink icon href={openTargets(v)}>Open Targets</ExternalLink>
            {rsid && <ExternalLink icon href={ensemblVar(rsid)}>Ensembl</ExternalLink>}
          </span> },
          { label: 'rsID match', value: v.match === 'exact' ? 'alleles match dbSNP'
            : v.match === 'position' ? 'position only (alleles differ from dbSNP record)'
            : 'no dbSNP record' },
        ]} />
        {/* every site in the store has both alleles (SPEC section 5), so neither is ever missing */}
        <KvTable rows={[
          { label: 'A1 / A2', value: `${v.A1} / ${v.A2}` },
          { label: 'A1', value: 'alternate allele (carries the effect)' },
          { label: 'A2', value: 'reference allele' },
        ]} />
      </div>

      <SectionPanel title="Lead variant for" description="Genes and splice phenotypes where this is the top cis association.">
        {!v.inCis ? <Empty label={OUTSIDE_CIS} /> : lists === null ? <TableSkeleton columns={[{ w: 'w-10' }, { w: 'w-16' }, { w: 'w-40' }, { w: 'w-20', align: 'right' }, { w: 'w-14', align: 'right' }, { w: 'w-12', align: 'right' }]} rows={2} /> :
          lists.leads.length === 0 ? <Empty label="Not the lead variant for any gene or splice phenotype." /> : (
            <div className="overflow-x-auto rounded-lg border border-base-300">
              <table className="table table-sm">
                <thead><tr><th>Type</th><th>Gene</th><th>Phenotype</th><th className="text-right">beta ± SE</th><th className="text-right">Perm p</th><th className="text-right">Status</th></tr></thead>
                <tbody>
                  {lists.leads.map(l => {
                    const val = leadValues(l.set.frame, l.row)
                    return (
                      <tr key={`${l.set.experiment.id}:${l.row}`} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${l.gene.gene_id}${l.qtlType === 's' ? '?tab=sqtl' : ''}`)}>
                        <td><span className={`badge badge-xs ${l.qtlType === 'e' ? 'badge-primary' : 'badge-secondary'}`}>{l.qtlType === 'e' ? 'eQTL' : 'sQTL'}</span></td>
                        <td className="font-medium"><span className={ROW_LINK_TEXT}>{l.gene.symbol ?? l.gene.gene_id}</span></td>
                        <td className="tabular-nums text-base-content/60">{l.phenotypeId ? fmtPhenotype(l.phenotypeId) : l.gene.gene_id}</td>
                        <td className="text-right tabular-nums">{fmtBetaSE(l.beta, l.betaSe)}</td>
                        <td className="text-right tabular-nums">{fmtP(val.sigValue)}</td>
                        <td className="text-right">{signif.assessed && val.significant ? <span className={`badge badge-xs ${l.qtlType === 'e' ? 'badge-primary' : 'badge-secondary'}`}>{l.qtlType === 'e' ? 'eGene' : 'sQTL'}</span> : ''}</td>
                      </tr>
                    )
                  })}
                </tbody>
              </table>
            </div>
          )}
      </SectionPanel>

      <SectionPanel title="Credible-set membership" description="SuSiE 95% credible sets containing this variant.">
        {!v.inCis ? <Empty label={OUTSIDE_CIS} /> : lists === null ? <TableSkeleton columns={[{ w: 'w-10' }, { w: 'w-16' }, { w: 'w-40' }, { w: 'w-6' }, { w: 'w-10', align: 'right' }]} rows={2} /> :
          lists.cs.length === 0 ? <Empty label="Not in any credible set." /> : (
            <div className="overflow-x-auto rounded-lg border border-base-300">
              <table className="table table-sm">
                <thead><tr><th>Type</th><th>Gene</th><th>Phenotype</th><th>Set</th><th className="text-right">PIP</th></tr></thead>
                <tbody>
                  {lists.cs.map(c => {
                    const val = csValues(c.set.frame, c.row)
                    return (
                      <tr key={`${c.set.experiment.id}:${c.row}`} className={`${ROW_LINK} hover:bg-base-200`} {...rowLink(`/gene/${c.gene.gene_id}${c.qtlType === 's' ? '?tab=sqtl' : ''}`)}>
                        <td><span className={`badge badge-xs ${c.qtlType === 'e' ? 'badge-primary' : 'badge-secondary'}`}>{c.qtlType === 'e' ? 'eQTL' : 'sQTL'}</span></td>
                        <td className="font-medium"><span className={ROW_LINK_TEXT}>{c.gene.symbol ?? c.gene.gene_id}</span></td>
                        <td className="tabular-nums text-base-content/60">{c.phenotypeId ? fmtPhenotype(c.phenotypeId) : c.gene.gene_id}</td>
                        <td><span className="badge badge-ghost badge-sm">{val.csId}</span></td>
                        <td className="text-right tabular-nums font-medium">{fmtNum(val.pip)}</td>
                      </tr>
                    )
                  })}
                </tbody>
              </table>
            </div>
          )}
      </SectionPanel>

      <SectionPanel title="trans associations" description={<>Genes and splice phenotypes anywhere in the genome whose expression or splicing this variant associates with, outside their cis windows. <RoundingNote kind="trans" /></>}>
        {hasTrans ? <TransTable table={transTable} keyedBy="variant" fileStem={`${rsid ?? `${v.chr}_${v.position}`}_trans`} />
          : <Unavailable what="trans eQTL and sQTL results" />}
      </SectionPanel>

      <SectionPanel title="All cis associations" description="Nominal statistics for every gene and splice phenotype whose window covers this variant. It runs on request."
        action={v.inCis && scan === null && plan && <button data-scan-button className="btn btn-sm h-8 rounded-lg border-base-300 font-medium" onClick={doScan}>
          Scan cis windows ({fmtBytes(plan.bytes)}{plan.sPending ? ' + splicing' : ''})
        </button>}>
        {!v.inCis ? <Empty label={OUTSIDE_CIS} /> : scan === null ? <Empty label="Not scanned yet." /> : scan === 'running' ? <TableSkeleton columns={[{ w: 'w-16' }, { w: 'w-14', align: 'right' }, { w: 'w-10', align: 'right' }, { w: 'w-14', align: 'right' }, { w: 'w-20', align: 'right' }, { w: 'w-16', align: 'right' }]} rows={6} /> : (
          <div className="space-y-4" data-scan-ready="1">
            <HitTable title={`Expression (${scan.e.length})`} hits={scan.e} qtlType="e" />
            <HitTable title={`Splicing (${shown.length})`} hits={shown} qtlType="s"
              toggle={scan.s.length > sig.length ? { allIntrons, total: scan.s.length, onToggle: () => setAllIntrons(x => !x) } : undefined} />
          </div>
        )}
      </SectionPanel>
    </div>
  )
}

function HitTable({ title, hits, qtlType, toggle }: {
  title: string; hits: CisHit[]; qtlType: 'e' | 's'
  toggle?: { allIntrons: boolean; total: number; onToggle: () => void }
}) {
  const rowLink = useRowLink()
  return (
    <div className="space-y-2">
      <div className="flex items-baseline justify-between gap-4">
        <h3 className="text-sm font-medium">{title}</h3>
        {toggle && <button className="link text-sm" onClick={toggle.onToggle}>
          {toggle.allIntrons ? 'Show significant introns only' : `Show all ${toggle.total} tested introns`}
        </button>}
      </div>
      {hits.length === 0 ? <Empty label="No windows cover this variant." /> : (
        <div className="overflow-x-auto rounded-lg border border-base-300">
          <table className="table table-sm">
            <thead><tr><th>Gene</th>{qtlType === 's' && <th>Phenotype</th>}<th className="text-right">TSS dist</th><th className="text-right">AF</th><th className="text-right">p</th><th className="text-right">beta ± SE</th><th className="text-right">PIP</th></tr></thead>
            <tbody>
              {hits.map((h, i) => (
                <tr key={i} className={`${ROW_LINK} hover:bg-base-200 ${h.pip != null ? 'bg-base-200/70' : ''}`} {...rowLink(`/gene/${h.gene_id}${qtlType === 's' ? '?tab=sqtl' : ''}`)}>
                  <td className="font-medium"><span className={ROW_LINK_TEXT}>{h.symbol ?? h.gene_id}</span></td>
                  {qtlType === 's' && <td className="tabular-nums text-base-content/60">{fmtPhenotype(h.phenotype_id ?? '')}</td>}
                  <td className="text-right tabular-nums text-base-content/60">{fmtBp(h.tss_distance)}</td>
                  <td className="text-right tabular-nums text-base-content/60">{fmtNum(h.af)}</td>
                  <td className="text-right tabular-nums">{fmtP(h.pval_nominal)}</td>
                  <td className="text-right tabular-nums">{fmtBetaSE(h.beta, h.beta_se)}</td>
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
