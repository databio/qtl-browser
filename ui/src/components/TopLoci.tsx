import { useEffect, useMemo, useState } from 'react'
import { useNavigate } from 'react-router'
import { GenomeTrack } from '@/components/genome-track/genome-track'
import type { TrackLocus } from '@/components/genome-track/types'
import { SectionPanel } from '@/components/section-panel'
import { useIndex } from '@/contexts/store-context'
import { fetchChromSizes, type ChromSizes } from '@/lib/chrom-sizes'
import { topLoci } from '@/lib/gene-index'
import { ucsc } from '@/lib/links'
import { gwasBinsOf, phenotypeLabel, type GwasExperiment } from '@/lib/store'
import { toastError } from '@/lib/toast'

/** At most this many markers. */
const N = 20

/** Series colors are theme tokens so a palette change flows through. */
const COLORS: Record<string, string> = {
  eQTL: 'var(--color-primary)',
  sQTL: 'var(--color-secondary)',
  GWAS: 'var(--color-accent)',
}

/** A study's strongest loci on the genome: each gene's lead cis association, the strongest per
 *  1 Mb, eQTL and sQTL together (lib/gene-index.ts `topLoci`). A marker opens its gene. */
export function StudyTopLoci() {
  const navigate = useNavigate()
  const [loci, setLoci] = useState<(TrackLocus & { tab: string })[] | null>(null)
  useEffect(() => {
    let alive = true
    topLoci(N).then(rows => {
      if (!alive) return
      setLoci(rows.filter(r => r.gene_id).map(r => {
        const qtl = phenotypeLabel(r.phenotype_type)
        return { id: r.gene_id!, chr: r.chr, start: r.lead_pos, end: r.lead_pos, trait: qtl, tab: qtl === 'sQTL' ? 'sqtl' : 'eqtl',
          label: r.symbol ?? r.gene_id! }
      }))
    }, e => { toastError('The top loci could not be loaded.', e); if (alive) setLoci([]) })
    return () => { alive = false }
  }, [])
  const tabOf = useMemo(() => new Map((loci ?? []).map(l => [l.id, l.tab])), [loci])
  if (loci && loci.length === 0) return null
  return (
    <LociTrack loci={loci} series={['eQTL', 'sQTL']}
      description={`The ${N} strongest genes by lead cis association, one per 1 Mb. Click a marker to open the gene.`}
      onSelect={id => navigate(`/gene/${id}?tab=${tabOf.get(id) ?? 'eqtl'}`)} />
  )
}

/** A GWAS's strongest windows on the genome: the lead variant of each of its strongest bins (the
 *  bin summary, one small object). A marker opens the lead in UCSC. */
export function GwasTopLoci({ gwas }: { gwas: GwasExperiment }) {
  const [loci, setLoci] = useState<TrackLocus[] | null>(null)
  useEffect(() => {
    let alive = true
    gwasBinsOf(gwas).then(bins => {
      if (!alive) return
      setLoci([...(bins ?? [])].sort((a, b) => a.min_p - b.min_p).slice(0, N).map(b => ({
        id: `${b.chr}:${b.lead_position}`, chr: b.chr, start: b.lead_position, end: b.lead_position, trait: 'GWAS',
        label: b.lead_rsid ?? `${b.chr}:${b.lead_position}`,
      })))
    }, e => { toastError('The top loci could not be loaded.', e); if (alive) setLoci([]) })
    return () => { alive = false }
  }, [gwas])
  if (loci && loci.length === 0) return null
  const mb = gwas.bins ? `${gwas.bins.bin_bp / 1e6} Mb ` : ''
  return (
    <LociTrack loci={loci} series={[]}
      description={`The lead variants of the ${N} strongest ${mb}windows; nearby signals in one window show as one. Click a marker to open it in UCSC.`}
      onSelect={id => {
        const [chr, pos] = id.split(':')
        window.open(ucsc(chr!, Math.max(1, Number(pos) - 250_000), Number(pos) + 250_000), '_blank', 'noopener')
      }} />
  )
}

/** The whole-genome track with markers only: chr1-22, X and Y of the store's own reference. */
function LociTrack({ loci, series, description, onSelect }: {
  loci: TrackLocus[] | null; series: string[]; description: string; onSelect: (id: string) => void
}) {
  const digest = useIndex()?.store.collection_digest ?? null
  const [chrom, setChrom] = useState<ChromSizes | null>(null)
  const [chromError, setChromError] = useState(false)
  useEffect(() => {
    if (!digest) return
    fetchChromSizes(digest)
      .then(setChrom)
      .catch((e: Error) => { setChromError(true); toastError('Chromosome sizes could not be fetched; the genome track is hidden.', e) })
  }, [digest])
  const legend = series.length > 0 && (
    <span className="flex items-center gap-4 text-xs text-base-content/60">
      {series.map(t => (
        <span key={t} className="inline-flex items-center gap-1.5">
          {/* inline style: the swatch is the marker's series color, a theme token resolved by the browser */}
          <span className="inline-block size-2.5 rounded-full" style={{ backgroundColor: COLORS[t] }} />{t}
        </span>
      ))}
    </span>
  )
  if (chromError) return null
  return (
    <SectionPanel title="Top loci" description={description} action={legend || undefined}>
      {!chrom ? <div className="h-[120px]" aria-busy="true" /> : (
        <div className="min-w-0">
          <GenomeTrack loci={loci ?? []} onLocusSelect={id => { if (id) onSelect(id) }}
            chromNames={chrom.names} chromLengths={chrom.lengths} traitColors={COLORS} static loading={loci === null} />
        </div>
      )}
    </SectionPanel>
  )
}
