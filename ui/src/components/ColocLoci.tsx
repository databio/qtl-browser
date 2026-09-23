import { useEffect, useMemo, useState } from 'react'
import { useNavigate } from 'react-router'
import { GenomeTrack } from '@/components/genome-track/genome-track'
import type { TrackBin, TrackLocus } from '@/components/genome-track/types'
import { SectionPanel } from '@/components/section-panel'
import { fetchChromSizes, type ChromSizes } from '@/lib/chrom-sizes'
import { fmtInt, fmtP } from '@/lib/format'
import { genesIPC, getStore, gwasBins } from '@/lib/store'
import type { GwasBin } from '@/lib/store-decode'
import { tableFromIPC } from 'apache-arrow'
import { COLOC_EQTL_GENES, COLOC_SQTL_GENES } from '@/lib/coloc'
import { Unavailable } from '@/components/states'
import { useStoreInfo } from '@/contexts/store-context'

/** A coloc locus. */
interface ColocLocus { gene_id: string; symbol: string; chr: string; tss: number; trait: 'eQTL' | 'sQTL' | 'both' }

const BIN_CAP = 20   // -log10 p; BAG3 and a couple of others exceed it and are drawn clipped

/** Series colors are theme tokens so a palette change flows through. */
const TRAIT_COLORS: Record<string, string> = {
  eQTL: 'var(--color-primary)',
  sQTL: 'var(--color-secondary)',
  both: 'var(--color-accent)',
}

/**
 * The paper's DCM-colocalized loci on a static whole-genome track. Every marker is labeled;
 * clicking one opens the gene page. The loci are the authors' gene list (lib/coloc.ts) placed at
 * each gene's TSS from the store's annotation, shown when the experiment carries the DCM GWAS. The
 * bars are the GWAS's strongest p per 5 Mb bin, from the bin summary the experiment's `gwas.bins`
 * names (v0's `gwas_dcm_bins.json`, same bins and values): one small whole-object read.
 */
export default function ColocLoci() {
  const s = useStoreInfo()
  if (s && !s.hasGwas) return (
    <SectionPanel title="Loci colocalized with dilated cardiomyopathy risk" description="Single-locus coloc with the Jurgens et al. 2024 DCM GWAS, PP.H4 > 0.8.">
      <Unavailable what="Colocalization and DCM GWAS data" />
    </SectionPanel>
  )
  if (!s) return <div className="h-[268px]" aria-busy="true" />
  return <ColocTrack />
}

/** The coloc genes at their annotated TSS, on the store's chromosomes (plain JS over the gene table
 *  the app already fetched at boot; no query engine needed). */
async function colocLoci(): Promise<ColocLocus[]> {
  const [s, ipc] = await Promise.all([getStore(), genesIPC()])
  const t = tableFromIPC(ipc)
  const name = t.getChild('name')!, id = t.getChild('gene_id')!, chr = t.getChild('chr')!, tss = t.getChild('tss')!
  const e = new Set(COLOC_EQTL_GENES), q = new Set(COLOC_SQTL_GENES)
  const out: ColocLocus[] = []
  for (let i = 0; i < t.numRows; i++) {
    const sym = name.get(i) as string | null
    if (!sym || !(e.has(sym) || q.has(sym)) || !s.chroms.has(chr.get(i))) continue
    out.push({ gene_id: id.get(i), symbol: sym, chr: chr.get(i), tss: Number(tss.get(i)),
      trait: e.has(sym) && q.has(sym) ? 'both' : q.has(sym) ? 'sQTL' : 'eQTL' })
  }
  return out
}

/** The track itself. */
function ColocTrack() {
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
    colocLoci().then(setHits).catch(e => { console.error(e); setHits([]) })
    gwasBins().then(b => setGwas(b ?? [])).catch(e => { console.error(e); setGwas([]) })
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
