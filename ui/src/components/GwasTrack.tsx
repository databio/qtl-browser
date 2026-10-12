import { useEffect, useMemo, useState, type ReactNode } from 'react'
import { useNavigate } from 'react-router'
import { GenomeTrack } from '@/components/genome-track/genome-track'
import type { TrackBin, TrackLocus } from '@/components/genome-track/types'
import { SectionPanel } from '@/components/section-panel'
import { fetchChromSizes, type ChromSizes } from '@/lib/chrom-sizes'
import { fmtInt, fmtP } from '@/lib/format'
import { ChevronDown } from 'lucide-react'
import { colocLists, getStore, gwasBinsOf, lookupGenes, type GwasExperiment } from '@/lib/store'
import { gwasName, gwasSummary } from '@/lib/describe'
import type { GwasBin } from '@/lib/store-decode'
import { colocGenes, colocTrait } from '@/lib/coloc'
import { allHits } from '@/lib/gene-index'
import { chooseGwas, SCOPE } from '@/lib/scope'
import { toastError } from '@/lib/toast'
import { useIndex, useStoreInfo } from '@/contexts/store-context'

/** A marker on the track: a gene at its TSS. */
interface Marker { gene_id: string; symbol: string; chr: string; tss: number; trait: string }

const BIN_CAP = 20   // -log10 p; the strongest bins exceed it and are drawn clipped
const GWS_P = 5e-8

/** Series colors are theme tokens so a palette change flows through. */
const TRAIT_COLORS: Record<string, string> = {
  eQTL: 'var(--color-primary)',
  sQTL: 'var(--color-secondary)',
  both: 'var(--color-accent)',
}

/** The GWAS experiments a visitor can choose, as a dropdown of short names with each one's size
 *  dimmed; the ones a coloc list pairs with this experiment are tagged. A filter box appears once
 *  the list outgrows a glance. Choosing reopens the page (lib/scope.ts). */
export function GwasPicker() {
  const ix = useIndex()
  const [q, setQ] = useState('')
  if (!ix || !SCOPE?.experiment) return null
  const listed = new Set(ix.colocs.filter(l => l.qtl === SCOPE!.experiment).map(l => l.gwas))
  const chosen = ix.gwas.find(g => g.id === SCOPE!.gwas)
  const needle = q.trim().toLowerCase()
  const shown = ix.gwas.filter(g => !needle || [gwasName(g, ix), g.title ?? '', g.id].some(t => t.toLowerCase().includes(needle)))
  const pick = (id: string | null) => { (document.activeElement as HTMLElement | null)?.blur(); chooseGwas(id) }
  return (
    <div className="dropdown dropdown-end">
      <button tabIndex={0} className="btn btn-xs gap-1 rounded-lg border-base-300 font-medium">
        {chosen ? <>GWAS: {gwasName(chosen, ix)}</> : 'Choose a GWAS'} <ChevronDown className="size-3" />
      </button>
      <div tabIndex={0} className="dropdown-content z-20 mt-1 w-80 rounded-lg border border-base-300 bg-base-100 p-1 shadow-lg">
        {ix.gwas.length > FILTER_AT && (
          <input className="input input-xs mb-1 w-full" placeholder="Filter" value={q} onChange={e => setQ(e.target.value)} />
        )}
        <ul className="max-h-72 overflow-y-auto">
          {chosen && (
            <li><button className="w-full rounded px-2 py-1 text-left text-xs text-base-content/60 hover:bg-base-200" onClick={() => pick(null)}>None</button></li>
          )}
          {shown.map(g => (
            <li key={g.id}>
              <button onClick={() => pick(g.id)}
                className={`w-full rounded px-2 py-1 text-left hover:bg-base-200 ${g.id === SCOPE!.gwas ? 'bg-primary/10' : ''}`}>
                <span className="flex items-baseline gap-2 text-xs">
                  <span className="truncate font-medium">{gwasName(g, ix)}</span>
                  {listed.has(g.id) && <span className="badge badge-xs badge-accent badge-outline shrink-0">coloc list</span>}
                  <span className="ml-auto shrink-0 tabular-nums text-base-content/45">{gwasSummary(g)}</span>
                </span>
                {g.title && <span className="block truncate text-[11px] text-base-content/50">{g.title}</span>}
              </button>
            </li>
          ))}
          {shown.length === 0 && <li className="px-2 py-1 text-xs text-base-content/45">No GWAS matches.</li>}
        </ul>
      </div>
    </div>
  )
}

/** Past this many GWAS, the picker gets a filter box. */
const FILTER_AT = 6

/** The chosen GWAS's genome track for a QTL collection. Markers are the genes a coloc list for the
 *  pair names; without a list, each genome-wide significant bin's nearest significant gene of the
 *  collection, from the same whole-genome read the search box makes. */
export function StudyGwasTrack() {
  const s = useStoreInfo()
  const ix = useIndex()
  const [markers, setMarkers] = useState<{ list: Marker[]; listed: boolean } | null>(null)
  const [bins, setBins] = useState<GwasBin[] | null>(null)
  useEffect(() => {
    if (!s?.gwas) return
    let alive = true
    const g = s.gwas
    gwasBinsOf(g).then(b => { if (alive) setBins(b ?? []) }).catch(e => { toastError('The GWAS track could not be loaded.', e); if (alive) setBins([]) })
    studyMarkers(g).then(m => { if (alive) setMarkers(m) }).catch(e => { toastError('The track markers could not be loaded.', e); if (alive) setMarkers({ list: [], listed: false }) })
    return () => { alive = false }
  }, [s])
  if (!s) return null
  if (!s.gwas) return (
    <SectionPanel title="GWAS" description="Choose a GWAS to see its genome track, and coloc and LocusCompare on the gene pages." action={<GwasPicker />}>
      <div />
    </SectionPanel>
  )
  const what = markers?.listed
    ? 'Markers: genes a coloc list reports for this pair.'
    : 'Markers: the nearest significant gene to each genome-wide significant bin.'
  return (
    <GwasTrack title="GWAS" gwas={s.gwas} bins={bins} markers={markers?.list ?? null}
      description={<>{gwasName(s.gwas, ix)}. {what} Click a marker to open the gene.</>} action={<GwasPicker />} />
  )
}

async function studyMarkers(g: GwasExperiment): Promise<{ list: Marker[]; listed: boolean }> {
  const s = await getStore()
  const lists = await colocLists()
  if (lists.length) {
    const byGene = colocGenes(lists)
    const ids = [...byGene.keys()]
    const found = await Promise.all(ids.map(async id => (await lookupGenes(id)).filter(r => r.gene_id === id)))
    const list: Marker[] = []
    ids.forEach((id, i) => {
      const traits = [...new Set(byGene.get(id)!.map(colocTrait))]
      for (const r of found[i]) if (s.chroms.has(r.chr))
        list.push({ gene_id: id, symbol: r.name || id, chr: r.chr, tss: r.tss, trait: traits.length > 1 ? 'both' : traits[0]! })
    })
    return { list, listed: true }
  }
  const [bins, hits] = await Promise.all([gwasBinsOf(g), allHits()])
  const sig = hits.filter(h => h.is_egene || h.n_sqtl_sig > 0)
  const list: Marker[] = []
  const seen = new Set<string>()
  for (const b of bins ?? []) {
    if (b.n_gws === 0 || b.min_p > GWS_P) continue
    let best: (typeof sig)[number] | null = null
    for (const h of sig) if (h.chr === b.chr && Math.abs(h.tss - b.lead_position) <= 1_000_000
      && (!best || Math.abs(h.tss - b.lead_position) < Math.abs(best.tss - b.lead_position))) best = h
    if (!best || seen.has(best.gene_id)) continue
    seen.add(best.gene_id)
    list.push({ gene_id: best.gene_id, symbol: best.symbol ?? best.gene_id, chr: best.chr, tss: best.tss,
      trait: best.is_egene && best.n_sqtl_sig > 0 ? 'both' : best.is_egene ? 'eQTL' : 'sQTL' })
  }
  return { list, listed: false }
}

/** A GWAS collection's own track: its bins, no markers. */
export function GwasOnlyTrack({ gwas }: { gwas: GwasExperiment }) {
  const ix = useIndex()
  const [bins, setBins] = useState<GwasBin[] | null>(null)
  useEffect(() => {
    let alive = true
    gwasBinsOf(gwas).then(b => { if (alive) setBins(b ?? []) }).catch(e => { toastError('The GWAS track could not be loaded.', e); if (alive) setBins([]) })
    return () => { alive = false }
  }, [gwas])
  return <GwasTrack title={gwasName(gwas, ix)} gwas={gwas} bins={bins} markers={[]} />
}

function GwasTrack({ title, description, action, gwas, bins, markers }: {
  title: string; description?: ReactNode; action?: ReactNode
  gwas: GwasExperiment; bins: GwasBin[] | null; markers: Marker[] | null
}) {
  const ix = useIndex()
  const digest = ix?.store.collection_digest ?? null
  const [chrom, setChrom] = useState<ChromSizes | null>(null)
  const [chromError, setChromError] = useState<string | null>(null)
  const [hovered, setHovered] = useState<string | null>(null)
  const [skipped, setSkipped] = useState(0)
  const navigate = useNavigate()
  useEffect(() => {
    if (!digest) return
    fetchChromSizes(digest)
      .then(c => { const keep = c.names.map((n, i) => [n, c.lengths[i]!] as const).filter(([n]) => n !== 'chrX' && n !== 'chrY')
        setChrom({ names: keep.map(k => k[0]), lengths: keep.map(k => k[1]) }) })
      .catch((e: Error) => { setChromError(e.message); toastError('Chromosome sizes could not be fetched; the genome track is hidden.', e) })
  }, [digest])

  const binBp = gwas.bins?.bin_bp
  const trackBins: TrackBin[] = useMemo(() => (bins ?? []).map(b => ({
    chr: b.chr, start: b.bin_start, end: b.bin_end, value: -Math.log10(Math.max(b.min_p, 1e-300)), sig: b.n_gws > 0,
    title: `${b.chr}:${(b.bin_start / 1e6).toFixed(0)}–${(b.bin_end / 1e6).toFixed(0)} Mb · strongest p ${fmtP(b.min_p)} at ${b.lead_rsid ?? fmtInt(b.lead_position)}` +
      (b.n_gws > 0 ? ` · ${fmtInt(b.n_gws)} genome-wide significant variants` : '') + ` · ${fmtInt(b.n_variants)} tested`,
  })), [bins])
  const loci: TrackLocus[] = useMemo(() =>
    (markers ?? []).map(h => ({ id: h.gene_id, chr: h.chr, start: h.tss, end: h.tss, label: h.symbol, trait: h.trait })), [markers])

  const legend = (
    <span className="flex flex-wrap items-center justify-end gap-x-4 gap-y-1 text-xs text-base-content/60">
      {loci.length > 0 && (['eQTL', 'sQTL', 'both'] as const).map(t => (
        <span key={t} className="inline-flex items-center gap-1.5">
          {/* inline style: the swatch is the marker's series color, a theme token resolved by the browser */}
          <span className="inline-block size-2.5 rounded-full" style={{ backgroundColor: TRAIT_COLORS[t] }} />
          {t === 'both' ? 'eQTL and sQTL' : t}
        </span>
      ))}
      {trackBins.length > 0 && (
        <span className="inline-flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-1.5 bg-error" />
          p &lt; 5×10⁻⁸{binBp ? ` · ${binBp / 1e6} Mb bins` : ''}
        </span>
      )}
    </span>
  )

  return (
    <SectionPanel title={title} description={description}
      action={<span className="flex flex-col items-end gap-1.5">{action}{bins !== null && markers !== null && legend}</span>}>
      {!chromError && !chrom && <div className="h-[268px]" aria-busy="true" />}
      {chrom && (
        <div className="min-w-0">
          <GenomeTrack loci={loci} hoveredLocusId={hovered ?? undefined}
            onLocusSelect={id => { setHovered(null); if (id) navigate(`/gene/${id}`) }}
            onSkipped={setSkipped} chromNames={chrom.names} chromLengths={chrom.lengths} traitColors={TRAIT_COLORS} static
            loading={bins === null || markers === null}
            bins={trackBins} binCap={BIN_CAP} binHeight={144} binUnit="−log10 p" />
          {skipped > 0 && <div className="text-xs text-warning">{skipped} locus{skipped > 1 ? 'i' : ''} on chromosomes not in the reference could not be placed.</div>}
        </div>
      )}
    </SectionPanel>
  )
}
