import { useEffect, useState, type ReactNode } from 'react'
import { Link, useParams, useSearchParams } from 'react-router'
import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import RoundingNote from '@/components/rounding-note'
import { SectionPanel } from '@/components/section-panel'
import { KvTable } from '@/components/kv-table'
import { Tooltip } from '@/components/tooltip'
import { Segmented } from '@/components/segmented'
import { DetailSkeleton, Empty, TableSkeleton, TabSkeleton } from '@/components/states'
import { ErrorPage } from '@/components/error-page'
import { toastError } from '@/lib/toast'
import CredibleSetTable from '@/components/CredibleSetTable'
import CisTable from '@/components/CisTable'
import { type LocusColoc } from '@/lib/coloc-abf'
import TransTable from '@/components/TransTable'
import LocusPlot, { LocusLegend } from '@/components/LocusPlot'
import { colocListing, colocSource, useColocLists } from '@/lib/coloc'
import { GwasPicker } from '@/components/GwasTrack'
import { ensemblGene, geneCards, gtexGene, openTargetsGene, ucsc } from '@/lib/links'
import { useCollectionCrumbs, useIndex, useSignificance, useStoreInfo } from '@/contexts/store-context'
import { CopyButton } from '@/components/copy-button'
import { fmtBp, fmtInt, fmtNum, fmtP, fmtPhenotype, fmtBetaSE } from '@/lib/format'
import { resolveGene, type GeneDetail, type SplicePhenotype,
  type CredibleSetRow, type Gene as GeneRow, type SearchHit } from '@/lib/queries'
import { loadGene, type GenePack } from '@/lib/gene'
import { geneTransTable } from '@/lib/trans'
import { EQTL_TYPE, getStore, SQTL_TYPE } from '@/lib/store'
import { gwasName } from '@/lib/describe'

const getStoreHasTrans = () => getStore().then(s => s.hasTrans)
import { dropTable, getDB } from '@/lib/db'

type Tab = 'eqtl' | 'sqtl'

/** The intron a gene's sQTL tab opens on: its first significant intron, else the one with the
 *  smallest tested significance value, else the first. An experiment that assessed no significance
 *  has neither, so it opens on the first intron. */
function defaultIntron(phens: SplicePhenotype[]): SplicePhenotype | null {
  return phens.find(p => p.is_sqtl)
    ?? phens.reduce<SplicePhenotype | null>((best, p) => (p.sig_value != null && (best === null || p.sig_value < (best.sig_value ?? Infinity)) ? p : best), null)
    ?? phens[0] ?? null
}

export default function Gene() {
  const { id = '' } = useParams()
  const [params, setParams] = useSearchParams()
  const sig = useSignificance()
  const colocs = useColocLists()
  const m = useStoreInfo()
  const ix = useIndex()
  const home = useCollectionCrumbs()
  const [hit, setHit] = useState<SearchHit | null | undefined>(undefined)
  // a gene tested for sQTL but not eQTL opens on its sQTL tab; any other tab value (old
  // `?tab=trans` links) falls back to eQTL
  // a tab per cis scan the study has (both until it has opened)
  const scans = { eqtl: m ? m.cis.has(EQTL_TYPE) : true, sqtl: m ? m.cis.has(SQTL_TYPE) : true }
  const asked: Tab = params.get('tab') === 'sqtl' ? 'sqtl' : params.get('tab') === 'eqtl' ? 'eqtl' : (hit && !hit.tested && hit.has_results ? 'sqtl' : 'eqtl')
  const tab: Tab = scans[asked] ? asked : asked === 'eqtl' ? 'sqtl' : 'eqtl'
  const [gp, setGp] = useState<GenePack | null>(null)
  const [detail, setDetail] = useState<GeneDetail | null>(null)
  // the gene's trans rows as an in-memory table, read once per gene alongside the other requests
  // and dropped when the gene changes; both tabs' trans tables page off it
  const [transTable, setTransTable] = useState<string | null>(null)
  // the gene could not be resolved or its block could not be read: the page cannot show anything
  const [error, setError] = useState<unknown>(null)
  // the query engine (the locus plot's) downloads while the gene's first requests are in flight
  useEffect(() => { getDB().catch(() => {}) }, [])
  useEffect(() => {
    let alive = true
    let table: string | null = null
    setHit(undefined); setGp(null); setDetail(null); setTransTable(null); setError(null)
    const fail = (e: unknown) => { console.error(e); if (alive) setError(e) }
    resolveGene(id).then(async h => {
      if (!alive) return
      setHit(h)
      if (!h) return
      // a trans experiment has no cis block to read: its page is the gene's trans associations
      if ((await getStore()).cis.size === 0) {
        try {
          const t = await geneTransTable(h)
          if (!alive) { dropTable(t); return }
          table = t
          setTransTable(t)
        } catch (e) { if (alive) toastError('The trans associations could not be loaded.', e) }
        return
      }
      if (!h.has_results) return
      performance.mark('gene:hit', { detail: h.gene_id })
      // the eQTL block and variants range requests leave together (or come from the cache)
      const pack = loadGene(h)
      setGp(pack)
      const tabParam = new URLSearchParams(window.location.search).get('tab')
      // a page opening on its sQTL tab sends the introns' span at once
      if (tabParam === 'sqtl' || (tabParam !== 'eqtl' && !h.tested)) pack.splice().catch(() => {})
      pack.detail
        .then(d => {
          if (!alive) return
          setDetail(d)
          performance.mark('gene:detail', { detail: h.gene_id })
        })
        .catch(fail)
      if (!(await getStoreHasTrans())) return
      try {
        const t = await geneTransTable(h)
        if (!alive) { dropTable(t); return }
        table = t
        setTransTable(t)
      } catch (e) { if (alive) toastError('The trans associations could not be loaded.', e) }
    }, fail)
    return () => { alive = false; if (table) dropTable(table) }
  }, [id])

  if (error) return <Page><ErrorPage title={`${id} could not be loaded.`} error={error} /></Page>
  // the page-level skeleton follows the tab in the URL: only the eQTL tab opens with a locus plot
  if (hit === undefined) return <Page><DetailSkeleton plot={tab === 'eqtl'} kvRows={tab === 'eqtl' ? 9 : 4} /></Page>
  if (hit === null) return <Page><Empty label={`No gene matches “${id}”.`} /></Page>

  const sym = hit.symbol ?? hit.gene_id
  const coloc = (trait: string) => {
    const l = colocListing(colocs, hit.gene_id, trait)
    return l && m?.gwas ? `${trait} colocalizes with ${gwasName(m.gwas, ix)} (${l.doc.threshold}), as reported by ${colocSource(l)}` : null
  }
  const colocE = coloc('eQTL'), colocS = coloc('sQTL')
  // a trans experiment (no cis scan): the gene's trans associations of that scan, and nothing else
  const transOnly = m != null && m.cis.size === 0
  const transTypes = transOnly ? [...m!.trans.keys()].map(t => (t === SQTL_TYPE ? 's' : 'e') as 'e' | 's') : []
  const tabs = [
    { value: 'eqtl' as Tab, label: 'eQTL' },
    { value: 'sqtl' as Tab, label: `sQTL${sig.assessed && hit.n_sqtl_sig ? ` (${hit.n_sqtl_sig})` : ''}` },
  ].filter(t => scans[t.value])
  return (
    <Page>
      <PageHeader
        crumbs={[...home, { to: '/genes', label: 'Genes' }, { label: sym }]}
        title={sym}
        meta={hit.gene_id}
        description={<span className="mt-1 flex flex-wrap items-center gap-1.5">
          {sig.assessed && hit.is_egene && <Chip cls="badge-primary" tip={`Significant cis-eQTL: ${sig.ruleText}`}>eGene</Chip>}
          {sig.assessed && hit.n_sqtl_sig > 0 && <Chip cls="badge-secondary" tip={`Introns with a significant cis-sQTL (${sig.ruleText})`}>{hit.n_sqtl_sig} sQTL intron{hit.n_sqtl_sig > 1 ? 's' : ''}</Chip>}
          {colocE && <Chip cls="badge-accent" tip={colocE}>coloc · eQTL</Chip>}
          {colocS && <Chip cls="badge-accent badge-outline" tip={colocS}>coloc · sQTL</Chip>}
          {!transOnly && !hit.has_results && <Chip cls="badge-ghost" tip="Filtered out before QTL mapping (expression or mappability)">not tested</Chip>}
          {/* only an eQTL scan can say a gene was left out of it */}
          {scans.eqtl && m != null && hit.has_results && !hit.tested && <Chip cls="badge-ghost" tip="Filtered out of the expression analysis">no eQTL test</Chip>}
        </span>}
        actions={hit.has_results && tabs.length > 1 ? <Segmented nav value={tab} onChange={t => setParams({ tab: t })} options={tabs} /> : undefined}
      />
      {transOnly ? (
        <div className="space-y-8">
          {transTypes.map(t => <TransSection key={t} table={transTable} qtlType={t} fileStem={`${sym}_trans_${t}qtl`} />)}
        </div>
      ) : !hit.has_results ? (
        <Empty label={`${sym} is annotated in ${m?.annotation.id ?? 'the gene annotation'} but was not tested for QTL in this experiment.`} />
      ) : (
        <>
          {tab === 'eqtl' && (detail && gp ? <EqtlTab hit={hit} gp={gp} d={detail} transTable={transTable} /> : <TabSkeleton plot chr={hit.chr} />)}
          {tab === 'sqtl' && (detail && gp ? <><GeneTable g={detail.gene} /><SqtlTab hit={hit} gp={gp} d={detail} transTable={transTable} /></> : <TabSkeleton kvRows={4} />)}
        </>
      )}
    </Page>
  )
}

/** Status chip with a hover explanation, since the label alone is terse. */
function Chip({ cls, tip, children }: { cls: string; tip: string; children: ReactNode }) {
  return <Tooltip tip={tip}><span className={`badge badge-sm cursor-help ${cls}`}>{children}</span></Tooltip>
}

function geneRows(g: GeneRow, annotation: string | undefined) {
  return [
    { label: 'Location', value: <span className="tabular-nums">{g.chr}:{fmtInt(g.start)}-{fmtInt(g.end)} ({g.strand})</span> },
    { label: 'TSS', value: <span className="tabular-nums">{g.chr}:{fmtInt(g.tss)}</span> },
    { label: 'Biotype', value: g.biotype.replace(/_/g, ' ') },
    // the version is the study annotation's; the link resolves the unversioned ID to Ensembl's current model
    { label: 'Ensembl ID', value: <span className="flex items-center justify-between gap-2">
      <ExternalLink icon href={ensemblGene(g.gene_id)}
        title={`Open in Ensembl. The version shown is from ${annotation ? `GENCODE ${annotation}` : 'the study annotation'}.`}>
        {g.gene_id_version}</ExternalLink>
      <CopyButton text={g.gene_id} className="-my-1 -mr-1" />
    </span> },
    { label: 'Links', value: <span className="flex flex-wrap gap-x-4">
      <ExternalLink icon href={ucsc(g.chr, g.start, g.end)}>UCSC</ExternalLink>
      {g.symbol && <ExternalLink icon href={gtexGene(g.symbol)}>GTEx</ExternalLink>}
      <ExternalLink icon href={openTargetsGene(g.gene_id)}>Open Targets</ExternalLink>
      {g.symbol && <ExternalLink icon href={geneCards(g.symbol)}>GeneCards</ExternalLink>}
    </span> },
  ]
}

function GeneTable({ g }: { g: GeneRow }) {
  const annotation = useStoreInfo()?.annotationVersion ?? undefined
  return <div className="mb-8 grid items-start gap-4 md:grid-cols-2"><KvTable rows={geneRows(g, annotation)} /></div>
}

/**
 * coloc.abf against the chosen GWAS, computed here from the two windows the page has already
 * downloaded (lib/coloc-abf.ts) -- no extra request, and about a millisecond.
 *
 * How far the shared-variant posterior is spread is the locus plot's business, not this table's:
 * the dots are sized by it. That matters for reading PP.H4, because a posterior smeared over tens
 * of variants is the regime where the p12 prior does much of the work -- FLNC falls from 0.94 to
 * 0.60 at p12 1e-6 while SKI holds at 0.999.
 */
function ColocSection({ qtlType, gp, phenotypeId }: {
  qtlType: 'e' | 's'; gp: GenePack; phenotypeId?: string
}) {
  const gwas = useStoreInfo()?.gwas ?? null
  const ix = useIndex()
  const [res, setRes] = useState<LocusColoc | null | 'running' | 'failed'>('running')
  useEffect(() => {
    let alive = true
    setRes('running')
    gp.coloc(qtlType, phenotypeId)
      .then(r => { if (alive) setRes(r) })
      .catch((e: Error) => { if (alive) { setRes('failed'); toastError('The GWAS colocalization could not be computed.', e) } })
    return () => { alive = false }
  }, [gp, qtlType, phenotypeId])

  return (
    <SectionPanel title="GWAS colocalization"
      description={gwas ? `coloc.abf against ${gwasName(gwas, ix)}. Values are computed in-browser.` : 'Choose a GWAS to compute coloc.abf in-browser.'}
      action={<GwasPicker />}>
      {!gwas ? null
        : res === 'running' ? <TableSkeleton columns={[{ w: 'w-16' }, { w: 'w-12', align: 'right' }]} rows={5} />
        : res === 'failed' ? <Empty label="The GWAS window could not be loaded." />
        : res === null ? <Empty label="No variants in this window are shared with the GWAS." />
        : (
          /* labelWidth: these labels are longer than the page's other KvTables and the values are
             short numbers, so the label column gets the room the value column does not need */
          <div className="grid items-start gap-4 md:grid-cols-2">
            <KvTable align="right" labelWidth="w-52" rows={[
              { label: 'H0 · no causal variant', value: fmtNum(res.pp[0], 3) },
              { label: 'H1 · QTL only', value: fmtNum(res.pp[1], 3) },
              { label: 'H2 · GWAS only', value: fmtNum(res.pp[2], 3) },
              { label: 'H3 · distinct causal variants', value: fmtNum(res.pp[3], 3) },
              { label: 'H4 · shared causal variant',
                value: <span className="font-medium text-base-content">{fmtNum(res.pp[4], 3)}</span> },
            ]} />
            <KvTable align="right" labelWidth="w-52" rows={[
              { label: 'Top shared variant', value: (() => {
                const id = res.top.rsNumber ? `rs${res.top.rsNumber}` : `${gp.hit.chr}:${res.top.position}`
                return <Link className="link-quiet tabular-nums" to={`/variant/${id}`}>{id}</Link>
              })() },
              { label: 'Top variant posterior', value: fmtNum(res.top.snpPP4, 3) },
              { label: 'Variants shared with GWAS', value: `${fmtInt(res.nShared)} of ${fmtInt(res.nQtlRows)}` },
              // `colocLocus` still counts the rows it had to drop (nNotTested, nUnderflow,
              // nNoStats) and they stay on LocusColoc for diagnostics, but the table does not
              // report counts of failures: the shared-of-total row above already says what the
              // run covered.
              { label: 'Priors', value: <span className="tabular-nums">p1 {res.priors.p1} · p2 {res.priors.p2} · p12 {res.priors.p12}</span> },
            ]} />
          </div>
        )}
    </SectionPanel>
  )
}

/** The locus plot's materialized window, handed up so the cis table can page off it. */
type LocusTable = { name: string | null; failed: boolean }
const NO_TABLE: LocusTable = { name: null, failed: false }

function EqtlTab({ hit, gp, d, transTable }: { hit: SearchHit; gp: GenePack; d: GeneDetail; transTable: string | null }) {
  const sig = useSignificance()
  const g = d.gene
  const [cs, setCs] = useState<CredibleSetRow[] | null>(null)
  const [nVar, setNVar] = useState<number | null>(null)
  const [legend, setLegend] = useState<{ sets: string[]; coloc: boolean } | null>(null)
  const [actions, setActions] = useState<ReactNode>(null)
  const [locus, setLocus] = useState<LocusTable>(NO_TABLE)
  const annotation = useStoreInfo()?.annotationVersion ?? undefined
  useEffect(() => { setCs(null); setNVar(null) }, [hit])
  const sym = hit.symbol ?? hit.gene_id
  if (!g.tested) return <><GeneTable g={g} /><Empty label={`${sym} was not tested for cis-eQTL (filtered out by expression or mappability); see the sQTL tab.`} /></>
  return (
    <div className="space-y-8">
      <div className="grid items-start gap-4 md:grid-cols-2">
        <KvTable rows={[
          ...geneRows(g, annotation),
          { label: 'Lead variant', value: <Link className="link-quiet" to={`/variant/${g.lead_rsid ?? `${g.chr}:${g.lead_position}`}`}>{g.lead_rsid ?? `${g.chr}:${fmtInt(g.lead_position)}`}</Link> },
          { label: 'A1 / A2', value: `${g.lead_A1} / ${g.lead_A2}` },
          { label: 'A1 frequency', value: fmtNum(g.lead_af) },
        ]} />
        <KvTable align="right" rows={[
          { label: 'Lead position', value: <span className="tabular-nums">{g.chr}:{fmtInt(g.lead_position)}</span> },
          { label: 'Lead distance to TSS', value: fmtBp(g.lead_tss_distance) },
          { label: 'Variants tested', value: fmtInt(g.num_var) },
          { label: 'beta ± SE', value: fmtBetaSE(g.beta, g.beta_se) },
          { label: 'Nominal p', value: fmtP(g.pval_nominal) },
          // the metric is the experiment's, so its name is too; an experiment that assessed none
          // shows neither row rather than an empty one
          ...(sig.assessed ? [{ label: sig.label, value: fmtP(g.sig_value) }] : []),
          ...(g.pval_beta != null ? [{ label: 'Beta-approximated p', value: fmtP(g.pval_beta) }] : []),
          { label: 'Credible sets', value: String(g.n_credible_sets) },
        ]} />
      </div>
      <SectionPanel title="Locus"
        description={<span className="inline-flex items-center gap-3 tabular-nums"><span>{g.chr}:{fmtInt(g.tss - 1_000_000)}–{fmtInt(g.tss + 1_000_000)}{nVar != null && ` · ${fmtInt(nVar)} variants`}</span>{actions}</span>}
        action={legend && <LocusLegend {...legend} />}>
        <LocusPlot spec={{ hit, pack: gp, qtlType: 'e', tss: g.tss, exons: d.exons }} onCount={setNVar} onLegend={setLegend} onActions={setActions} onCredibleSets={setCs}
          onTable={(name, failed) => setLocus({ name, failed: !!failed })} />
      </SectionPanel>
      <SectionPanel title="SuSiE 95% credible sets">
        {cs === null ? <TableSkeleton columns={[{ w: 'w-4' }, { w: 'w-12' }, { w: 'w-8', align: 'right' }, { w: 'w-24' }, { w: 'w-10', align: 'right' }, { w: 'w-10', align: 'right' }, { w: 'w-16', align: 'right' }]} rows={2} /> : <CredibleSetTable rows={cs} />}
      </SectionPanel>
      <ColocSection qtlType="e" gp={gp} />
      <SectionPanel title="cis associations" description={<>Every variant within ±1 Mb of the TSS; rows tinted when the variant is in a credible set. Click a row to open the variant. <RoundingNote /></>}>
        <CisTable table={locus.name} failed={locus.failed} chr={hit.chr} qtlType="e" fileStem={`${sym}_cis_eqtl`} />
      </SectionPanel>
      <TransSection table={transTable} qtlType="e" fileStem={`${sym}_trans_eqtl`} />
    </div>
  )
}

/** The sQTL tab: its intron list waits for the span of the gene's intron blocks (one request). */
function SqtlTab(props: { hit: SearchHit; gp: GenePack; d: GeneDetail; transTable: string | null }) {
  const [phens, setPhens] = useState<SplicePhenotype[] | null | 'failed'>(null)
  useEffect(() => {
    let alive = true
    setPhens(null)
    props.gp.splice().then(p => { if (alive) setPhens(p) }, e => { if (alive) { setPhens('failed'); toastError('The splicing phenotypes could not be loaded.', e) } })
    return () => { alive = false }
  }, [props.gp])
  if (phens === null) return <TabSkeleton kvRows={4} />
  if (phens === 'failed') return <Empty label="The splicing phenotypes could not be loaded." />
  return <SqtlIntrons {...props} phens={phens} />
}

function SqtlIntrons({ hit, gp, d, transTable, phens }: { hit: SearchHit; gp: GenePack; d: GeneDetail; transTable: string | null; phens: SplicePhenotype[] }) {
  const sig = useSignificance()
  const [cs, setCs] = useState<CredibleSetRow[] | null>(null)
  const [selected, setSelected] = useState<string | null>(() => defaultIntron(phens)?.phenotype_id ?? null)
  // every tested intron has its sQTL rows in the pack; the list starts with the significant ones
  // unless there are none
  const [showAll, setShowAll] = useState(() => !phens.some(p => p.is_sqtl))
  const [nVar, setNVar] = useState<number | null>(null)
  const [legend, setLegend] = useState<{ sets: string[]; coloc: boolean } | null>(null)
  const [actions, setActions] = useState<ReactNode>(null)
  const [locus, setLocus] = useState<LocusTable>(NO_TABLE)
  useEffect(() => {
    setSelected(defaultIntron(phens)?.phenotype_id ?? null); setShowAll(!phens.some(p => p.is_sqtl)); setNVar(null); setCs(null)
  }, [hit, phens])
  if (!phens.length) return <Empty label="No splicing phenotypes were tested for this gene." />
  const significant = phens.filter(p => p.is_sqtl)
  const listed = showAll ? phens : significant
  const sel = phens.find(p => p.phenotype_id === selected) ?? null
  const sym = hit.symbol ?? hit.gene_id
  return (
    <div className="space-y-8">
      <SectionPanel title="Splice phenotypes"
        description={`${significant.length ? `${fmtInt(significant.length)} of ${fmtInt(phens.length)} tested introns with a significant sQTL` : `None of the ${fmtInt(phens.length)} tested introns has a significant sQTL`}. Click a row to load its locus.`}
        action={
          <label className="inline-flex shrink-0 cursor-pointer items-center gap-1.5 text-xs text-base-content/70">
            <input type="checkbox" className="toggle toggle-xs" checked={showAll} onChange={e => setShowAll(e.target.checked)} />
            Show all {fmtInt(phens.length)} tested intron{phens.length === 1 ? '' : 's'}
          </label>
        }>
        {listed.length === 0 ? <Empty label={`None of the ${fmtInt(phens.length)} tested introns has a significant sQTL.`} /> : (
          <div className="overflow-x-auto rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead><tr><th>Cluster</th><th>Intron</th><th>Lead variant</th><th className="text-right">beta ± SE</th><th className="text-right">{sig.short || '—'}</th><th className="text-right">Sets</th></tr></thead>
              <tbody>
                {listed.map(p => (
                  <tr key={p.phenotype_id} onClick={() => setSelected(p.phenotype_id)}
                    className={`cursor-pointer transition-colors ${p.phenotype_id === selected ? 'bg-base-200' : 'hover:bg-base-200/60'}`}>
                    <td className="font-mono text-xs text-base-content/60">{p.cluster_id}</td>
                    <td className="tabular-nums">{fmtInt(p.intron_start)}–{fmtInt(p.intron_end)} <span className="text-base-content/50">({fmtBp(p.intron_end - p.intron_start)})</span></td>
                    <td><Link className="link-quiet" to={`/variant/${p.lead_rsid ?? `${p.chr}:${p.lead_position}`}`}>{p.lead_rsid ?? `${p.chr}:${fmtInt(p.lead_position)}`}</Link></td>
                    <td className="text-right tabular-nums">{fmtBetaSE(p.beta, p.beta_se)}</td>
                    <td className="text-right tabular-nums">{fmtP(p.sig_value)}</td>
                    <td className="text-right tabular-nums">{p.n_credible_sets}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
      </SectionPanel>
      {sel && (
        <>
          <SectionPanel title="Locus"
            description={<span className="inline-flex items-center gap-3 tabular-nums"><span>intron {fmtPhenotype(sel.phenotype_id)}{nVar != null && ` · ${fmtInt(nVar)} variants`}</span>{actions}</span>}
            action={legend && <LocusLegend {...legend} />}>
            <LocusPlot spec={{ hit, pack: gp, qtlType: 's', phenotypeId: sel.phenotype_id, tss: sel.tss, exons: d.exons, intron: { start: sel.intron_start, end: sel.intron_end } }} onCount={setNVar} onLegend={setLegend} onActions={setActions} onCredibleSets={setCs}
              onTable={(name, failed) => setLocus({ name, failed: !!failed })} />
          </SectionPanel>
          <SectionPanel title="SuSiE 95% credible sets">
            {cs === null ? <TableSkeleton columns={[{ w: 'w-4' }, { w: 'w-12' }, { w: 'w-8', align: 'right' }, { w: 'w-24' }, { w: 'w-10', align: 'right' }, { w: 'w-10', align: 'right' }, { w: 'w-16', align: 'right' }]} rows={2} /> : <CredibleSetTable rows={cs} />}
          </SectionPanel>
          <ColocSection qtlType="s" gp={gp} phenotypeId={sel.phenotype_id} />
          <SectionPanel title="cis associations" description={<>Every variant within ±1 Mb of the TSS for <b className="font-medium text-base-content/80">this intron</b>; rows tinted when the variant is in a credible set. Click a row to open the variant. <RoundingNote /></>}>
            <CisTable table={locus.name} failed={locus.failed} chr={hit.chr} qtlType="s" phenotypeId={sel.phenotype_id} fileStem={`${sym}_${sel.cluster_id}_${sel.intron_start}_${sel.intron_end}_cis_sqtl`} />
          </SectionPanel>
        </>
      )}
      <TransSection table={transTable} qtlType="s" fileStem={`${sym}_trans_sqtl`} />
    </div>
  )
}

/** Variants outside the cis window associated with this gene, for one QTL type. The eQTL
 *  table needs no phenotype column (the phenotype is the gene); the sQTL table names the
 *  intron, since a gene's trans sQTL rows can hit different introns. */
function TransSection({ table, qtlType, fileStem }: { table: string | null; qtlType: 'e' | 's'; fileStem: string }) {
  const what = qtlType === 'e' ? 'expression' : 'splicing'
  // a study with no trans scan of this type has no section for it
  const scanned = useStoreInfo()?.trans.has(qtlType === 'e' ? EQTL_TYPE : SQTL_TYPE) ?? false
  if (!scanned) return null
  return (
    <SectionPanel title="trans associations"
      description={<>Every variant outside the cis window associated with this gene's {what}{qtlType === 's' && <> <b className="font-medium text-base-content/80">(any intron)</b></>}. Click a row to open the variant. <RoundingNote kind="trans" /></>}>
      <TransTable table={table} qtlType={qtlType} fileStem={fileStem} />
    </SectionPanel>
  )
}
