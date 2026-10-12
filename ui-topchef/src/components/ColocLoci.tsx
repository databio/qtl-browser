import { useEffect, useMemo, useState } from 'react'
import { useNavigate } from 'react-router'
import { GenomeTrack } from '@/components/genome-track/genome-track'
import type { TrackBin, TrackLocus } from '@/components/genome-track/types'
import { SectionPanel } from '@/components/section-panel'
import { fetchChromSizes, type ChromSizes } from '@/lib/chrom-sizes'
import { fmtInt, fmtP } from '@/lib/format'
import { colocLists, getStore, gwasBins, lookupGenes } from '@/lib/store'
import type { GwasBin } from '@/lib/store-decode'
import { colocGenes, colocTrait } from '@/lib/coloc'
import { Unavailable } from '@/components/states'
import { toastError } from '@/lib/toast'
import { useStoreInfo } from '@/contexts/store-context'

/** A coloc locus. */
interface ColocLocus { gene_id: string; symbol: string; chr: string; tss: number; trait: string }

const BIN_CAP = 20   // -log10 p; BAG3 and a couple of others exceed it and are drawn clipped

/** Series colors are theme tokens so a palette change flows through. */
const TRAIT_COLORS: Record<string, string> = {
  eQTL: 'var(--color-primary)',
  sQTL: 'var(--color-secondary)',
  both: 'var(--color-accent)',
}

/**
 * The paper's DCM-colocalized loci on a static whole-genome track. Every marker is labeled;
 * clicking one opens the gene page. The loci are the store's coloc lists for the study and its GWAS
 * (SPEC section 20), each gene placed at its TSS through the annotation's gene lookup. The bars are
 * the GWAS's strongest p per 5 Mb bin, from the bin summary the GWAS experiment's `bins` names: one
 * small whole-object read.
 */
async function colocLoci(): Promise<ColocLocus[]> {
  const s = await getStore()
  if (!s.gwas) return []   // the panel shows "not available" instead; do not spend lookups on it
  const lists = await colocLists()
  const byGene = colocGenes(lists)
  const ids = [...byGene.keys()]
  // one small read per gene, in parallel, after the lookup header (lookupGenes shares it)
  const found = await Promise.all(ids.map(async g => (await lookupGenes(g)).filter(r => r.gene_id === g)))
  const loci: ColocLocus[] = []
  ids.forEach((g, i) => {
    const traits = [...new Set(byGene.get(g)!.map(colocTrait))]
    for (const r of found[i]) if (s.chroms.has(r.chr))
      loci.push({ gene_id: g, symbol: r.name || g, chr: r.chr, tss: r.tss, trait: traits.length > 1 ? 'both' : traits[0]! })
  })
  return loci
}

export default function ColocLoci() {
  const s = useStoreInfo()
  const [chrom, setChrom] = useState<ChromSizes | null>(null)
  const [chromError, setChromError] = useState<string | null>(null)
  const [hits, setHits] = useState<ColocLocus[] | null>(null)
  const [gwas, setGwas] = useState<GwasBin[] | null>(null)
  const [hovered, setHovered] = useState<string | null>(null)
  const [skipped, setSkipped] = useState(0)
  const navigate = useNavigate()

  // all three start at mount, not after the store has opened: chromosome sizes come from seqcol and
  // need no store at all, and the other two wait on `getStore()` themselves. Gating the whole
  // component on the store instead put two pointer round trips in front of the seqcol request and
  // held the page blank for both of them.
  useEffect(() => {
    // autosomes only: the coloc loci and the DCM GWAS are both autosomal
    fetchChromSizes('GRCh38')
      .then(c => { const keep = c.names.map((n, i) => [n, c.lengths[i]!] as const).filter(([n]) => n !== 'chrX' && n !== 'chrY')
        setChrom({ names: keep.map(k => k[0]), lengths: keep.map(k => k[1]) }) })
      .catch((e: Error) => { setChromError(e.message); toastError('Chromosome sizes could not be fetched; the genome track is hidden.', e) })
    colocLoci().then(setHits).catch(e => { toastError('The colocalized loci could not be loaded.', e); setHits([]) })
    gwasBins().then(b => setGwas(b ?? [])).catch(e => { toastError('The GWAS track could not be loaded.', e); setGwas([]) })
  }, [])

  const bins: TrackBin[] = useMemo(() => (gwas ?? []).map(b => ({
    chr: b.chr, start: b.bin_start, end: b.bin_end, value: -Math.log10(Math.max(b.min_p, 1e-300)), sig: b.n_gws > 0,
    title: `${b.chr}:${(b.bin_start / 1e6).toFixed(0)}–${(b.bin_end / 1e6).toFixed(0)} Mb · strongest p ${fmtP(b.min_p)} at ${b.lead_rsid ?? fmtInt(b.lead_position)}` +
      (b.n_gws > 0 ? ` · ${fmtInt(b.n_gws)} genome-wide significant variants` : '') + ` · ${fmtInt(b.n_variants)} tested`,
  })), [gwas])

  const loci: TrackLocus[] = useMemo(() =>
    (hits ?? []).map(h => ({ id: h.gene_id, chr: h.chr, start: h.tss, end: h.tss, label: h.symbol, trait: h.trait })), [hits])

  const legend = (
    <span className="flex flex-wrap items-center justify-end gap-x-4 gap-y-1 text-xs text-base-content/60">
      {(['eQTL', 'sQTL', 'both'] as const).map(t => (
        <span key={t} className="inline-flex items-center gap-1.5">
          {/* inline style: the swatch is the marker's series color, a theme token resolved by the browser */}
          <span className="inline-block size-2.5 rounded-full" style={{ backgroundColor: TRAIT_COLORS[t] }} />
          {t === 'both' ? 'eQTL and sQTL' : t}
        </span>
      ))}
      {bins.length > 0 && (
        <span className="inline-flex items-center gap-1.5">
          <span className="inline-block h-2.5 w-1.5 bg-error" />
          DCM GWAS p &lt; 5×10⁻⁸
        </span>
      )}
    </span>
  )

  // an experiment without the DCM GWAS has no bars and no coloc to show (after the hooks, so the
  // fetches above still start at mount whatever the store turns out to hold)
  if (s && !s.hasGwas) return (
    <SectionPanel title="Loci colocalized with dilated cardiomyopathy risk" description="Single-locus coloc with the Jurgens et al. 2024 DCM GWAS, PP.H4 ≥ 0.8.">
      <Unavailable what="Colocalization and DCM GWAS data" />
    </SectionPanel>
  )

  return (
    <SectionPanel
      title="Loci colocalized with dilated cardiomyopathy risk"
      description={`Single-locus coloc with the Jurgens et al. 2024 DCM GWAS, PP.H4 ≥ 0.8${loci.length ? `, ${loci.length} loci` : ''}. Click a marker to open the gene.`}
      action={hits !== null && gwas !== null ? legend : undefined}>
      {/* the track draws as soon as chromosome sizes are known (cached after the first visit);
          markers and GWAS bars fade in on top when their queries return, without moving it.
          Until then the space is held empty at the track's height: 80 pad + 40 bar/labels + 4 + 144 */}
      {!chromError && !chrom && <div className="h-[268px]" aria-busy="true" />}
      {chrom && (
        <div className="min-w-0">
          <GenomeTrack loci={loci} hoveredLocusId={hovered ?? undefined}
            onLocusSelect={id => { setHovered(null); if (id) navigate(`/gene/${id}`) }}
            onSkipped={setSkipped} chromNames={chrom.names} chromLengths={chrom.lengths} traitColors={TRAIT_COLORS} static
            loading={hits === null || gwas === null}
            bins={bins} binCap={BIN_CAP} binHeight={144} binUnit="−log10 p" />
          {skipped > 0 && <div className="text-xs text-warning">{skipped} locus{skipped > 1 ? 'i' : ''} on chromosomes not in the reference could not be placed.</div>}
        </div>
      )}
    </SectionPanel>
  )
}
