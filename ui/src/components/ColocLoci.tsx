import { useEffect, useMemo, useState } from 'react'
import { useNavigate } from 'react-router'
import { GenomeTrack } from '@/components/genome-track/genome-track'
import type { TrackBin, TrackLocus } from '@/components/genome-track/types'
import { SectionPanel } from '@/components/section-panel'
import { fetchChromSizes, type ChromSizes } from '@/lib/chrom-sizes'
import { fmtInt, fmtP } from '@/lib/format'
import { DATA_BASE } from '@/lib/manifest'

/** The pipeline's coloc_loci.json and gwas_dcm_bins.json: plain fetches, no query engine. */
interface ColocLocus { gene_id: string; symbol: string; chr: string; tss: number; trait: 'eQTL' | 'sQTL' | 'both' }
interface GwasBin {
  chr: string; bin_start: number; bin_end: number; min_p: number; lead_position: number; lead_rsid: string | null
  lead_beta: number; lead_ea: string; n_gws: number; n_variants: number
}
const getJSON = <T,>(name: string) => fetch(`${DATA_BASE}/${name}`).then(r => { if (!r.ok) throw new Error(`${name}: ${r.status}`); return r.json() as Promise<T> })
/** gwas_dcm_bins.json is columnar (keys once, one array per column); turn it back into rows */
const fromColumns = <T,>(d: { n: number; columns: Record<string, unknown[]> }): T[] =>
  Array.from({ length: d.n }, (_, i) => Object.fromEntries(Object.entries(d.columns).map(([k, v]) => [k, v[i]])) as T)

const BIN_CAP = 20   // -log10 p; BAG3 and a couple of others exceed it and are drawn clipped

/** Series colors are theme tokens so a palette change flows through. */
const TRAIT_COLORS: Record<string, string> = {
  eQTL: 'var(--color-primary)',
  sQTL: 'var(--color-secondary)',
  both: 'var(--color-accent)',
}

/**
 * The paper's DCM-colocalized loci on a static whole-genome track. Every marker is labeled;
 * clicking one opens the gene page.
 */
export default function ColocLoci() {
  const [chrom, setChrom] = useState<ChromSizes | null>(null)
  const [chromError, setChromError] = useState<string | null>(null)
  const [hits, setHits] = useState<ColocLocus[] | null>(null)
  const [gwas, setGwas] = useState<GwasBin[] | null>(null)
  const [hovered, setHovered] = useState<string | null>(null)
  const [skipped, setSkipped] = useState(0)
  const navigate = useNavigate()

  useEffect(() => {
    // autosomes only: the coloc loci and the DCM GWAS are both autosomal
    fetchChromSizes('GRCh38')
      .then(c => { const keep = c.names.map((n, i) => [n, c.lengths[i]!] as const).filter(([n]) => n !== 'chrX' && n !== 'chrY')
        setChrom({ names: keep.map(k => k[0]), lengths: keep.map(k => k[1]) }) })
      .catch((e: Error) => setChromError(e.message))
    getJSON<ColocLocus[]>('coloc_loci.json').then(setHits).catch(() => setHits([]))
    getJSON<{ n: number; columns: Record<string, unknown[]> }>('gwas_dcm_bins.json').then(d => setGwas(fromColumns<GwasBin>(d))).catch(() => setGwas([]))
  }, [])

  const bins: TrackBin[] = useMemo(() => (gwas ?? []).map(b => ({
    chr: b.chr, start: b.bin_start, end: b.bin_end, value: -Math.log10(Math.max(b.min_p, 1e-300)), sig: b.n_gws > 0,
    title: `${b.chr}:${(b.bin_start / 1e6).toFixed(0)}–${(b.bin_end / 1e6).toFixed(0)} Mb · strongest p ${fmtP(b.min_p)} at ${b.lead_rsid ?? fmtInt(b.lead_position)}` +
      (b.n_gws > 0 ? ` · ${fmtInt(b.n_gws)} genome-wide significant variants` : '') + ` · ${fmtInt(b.n_variants)} tested`,
  })), [gwas])

  const loci: TrackLocus[] = useMemo(() =>
    (hits ?? []).map(h => ({ id: h.gene_id, chr: h.chr, start: h.tss, end: h.tss, label: h.symbol ?? h.gene_id, trait: h.trait })), [hits])

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

  return (
    <SectionPanel
      title="Loci colocalized with dilated cardiomyopathy risk"
      description={`Single-locus coloc with the Jurgens et al. 2024 DCM GWAS, PP.H4 > 0.8${loci.length ? `, ${loci.length} loci` : ''}. Click a marker to open the gene.`}
      action={hits !== null && gwas !== null ? legend : undefined}>
      {chromError && <div className="text-xs text-error">Chromosome sizes unavailable ({chromError}); the track is hidden.</div>}
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
