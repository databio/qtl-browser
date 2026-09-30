import { useEffect, useRef, useState } from 'react'
import * as vg from '@uwdata/vgplot'
import { genesInWindow, type Exon, type WindowGene } from '@/lib/queries'

export interface TrackSpec {
  chr: string
  geneId: string
  domain: [number, number]
  exons: Exon[]                              // collapsed model of the gene of interest
  intron?: { start: number; end: number }   // the selected sQTL intron
}

type Lane = WindowGene & { lane: number; mid: number; label: string }

/** Pack genes into lanes so bars and labels do not overlap. Label width is estimated in bp
 *  from the pixel scale so a long symbol reserves room. */
function packLanes(genes: WindowGene[], bpPerPx: number): Lane[] {
  const laneEnds: number[] = []
  const out: Lane[] = []
  for (const g of genes) {
    const label = `${g.symbol ?? g.gene_id}${g.strand === '+' ? ' →' : ' ←'}`
    const labelBp = (label.length * 6.5 + 12) * bpPerPx
    const mid = (g.start + g.end) / 2
    const left = Math.min(g.start, mid - labelBp / 2)
    const right = Math.max(g.end, mid + labelBp / 2)
    let lane = laneEnds.findIndex(e => e < left)
    if (lane === -1) { lane = laneEnds.length; laneEnds.push(right) } else laneEnds[lane] = right
    out.push({ ...g, lane, mid, label })
  }
  return out
}

/** Genes in the plotted window on a shared x axis: thin bars for neighbors, the collapsed
 *  exon model for the gene of interest, and the selected intron shaded. Carries the x axis
 *  for the locus plot above it. */
export default function GeneTrack({ spec, width, marginLeft, dark, shade }: {
  spec: TrackSpec; width: number; marginLeft: number; dark: boolean
  /** the locus plot's brushed interval, shaded across every lane so the reader sees which slice the detail shows */
  shade?: [number, number] | null
}) {
  const host = useRef<HTMLDivElement>(null)
  const [data, setData] = useState<{ genes: WindowGene[]; exons: Exon[] } | null>(null)
  const [error, setError] = useState<string | null>(null)

  const key = `${spec.chr}|${spec.geneId}|${spec.domain.join('-')}`
  useEffect(() => {
    let alive = true
    setData(null); setError(null)
    // neighbours come from the chromosome's genes (read when the gene was resolved); the gene's own exons come with the gene details
    genesInWindow(spec.chr, spec.domain[0], spec.domain[1])
      .then(genes => { if (alive) setData({ genes, exons: spec.exons }) })
      .catch((e: Error) => { console.error(e); if (alive) setError(`gene track query failed: ${e.message}`) })
    return () => { alive = false }
  }, [key]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    const el = host.current!
    el.replaceChildren()
    if (!data) return
    const plotW = width - marginLeft - 20
    const bpPerPx = (spec.domain[1] - spec.domain[0]) / plotW
    const lanes = packLanes(data.genes.filter(g => g.gene_id !== spec.geneId), bpPerPx)
    const me = data.genes.find(g => g.gene_id === spec.geneId)
    const nLanes = Math.max(1, lanes.length ? Math.max(...lanes.map(l => l.lane)) + 1 : 0)
    // lane 0 (top) is the gene of interest; neighbors below
    const neighbors = lanes.map(l => ({ ...l, y: -(l.lane + 1) }))
    const ink = dark ? '#c3c2b7' : '#52514e'
    const muted = dark ? '#5a4d4a' : '#b8b3b1'   // warm grays matched to the theme surfaces
    // the gene of interest is drawn in ink, not a hue: the hues above it belong to credible sets
    const accent = ink
    // one lane unit is laneH px: the plot area is sized from the y domain so the spacing
    // does not shrink as lanes are added. The gene of interest's exon boxes are taller than
    // a neighbor bar, so the first neighbor (whose label rises 9 px above its bar) sits an
    // extra GENE_GAP below it instead of running into the exons.
    const laneH = 18
    const GENE_GAP = 0.35
    const yMin = -nLanes - GENE_GAP - 0.9
    const height = 16 + Math.round(laneH * (0.9 - yMin)) + 8
    // plain columns only: vgplot treats string channels as fields of the array data
    const nb = neighbors.map(l => ({ start: l.start, end: l.end, mid: l.mid, label: l.label, y: l.y - GENE_GAP, y1: l.y - GENE_GAP - 0.12, y2: l.y - GENE_GAP + 0.12 }))
    const meRow = me ? [{ start: me.start, end: me.end, mid: (me.start + me.end) / 2, y1: -0.05, y2: 0.05,
      label: `${me.symbol ?? me.gene_id}${me.strand === '+' ? ' →' : ' ←'}` }] : []
    const exons = data.exons.map(x => ({ start: x.start, end: x.end, y1: -0.3, y2: 0.3 }))
    const marks: unknown[] = []
    if (shade) {
      marks.push(vg.rect([{ start: shade[0], end: shade[1], y1: yMin, y2: 0.9 }],
        { x1: 'start', x2: 'end', y1: 'y1', y2: 'y2', fill: ink, fillOpacity: 0.08 }))
    }
    if (spec.intron) {
      marks.push(vg.rect([{ start: spec.intron.start, end: spec.intron.end, y1: yMin, y2: 0.9 }],
        { x1: 'start', x2: 'end', y1: 'y1', y2: 'y2', fill: accent, fillOpacity: 0.12 }))
    }
    if (nb.length) {
      marks.push(vg.rect(nb, { x1: 'start', x2: 'end', y1: 'y1', y2: 'y2', fill: muted }))
      // label sits closer to its own bar than the bar is to the next row's label
      marks.push(vg.text(nb, { x: 'mid', y: 'y', text: 'label', dy: -8, fontSize: 10, fill: ink, fillOpacity: 0.8 }))
    }
    if (meRow.length) {
      marks.push(vg.rect(meRow, { x1: 'start', x2: 'end', y1: 'y1', y2: 'y2', fill: accent }))
      if (exons.length) marks.push(vg.rect(exons, { x1: 'start', x2: 'end', y1: 'y1', y2: 'y2', fill: accent }))
      marks.push(vg.text(meRow.map(r => ({ ...r, y: 0 })), { x: 'mid', y: 'y', text: 'label', dy: -13, fontSize: 11, fontWeight: 600, fill: ink }))
    }
    try {
      const plot = vg.plot(
        ...marks,
        vg.xDomain(spec.domain), vg.xAxis(null), vg.yDomain([yMin, 0.9]), vg.yAxis(null),
        vg.xInset(8), vg.width(width), vg.height(height), vg.marginLeft(marginLeft), vg.marginRight(20), vg.marginTop(16), vg.marginBottom(8),
        vg.style({ fontFamily: 'inherit', fontSize: '11px', color: ink, background: 'transparent' }),
      ) as HTMLElement
      el.replaceChildren(plot)
    } catch (e) {
      console.error(e)
      setError(`gene track render failed: ${(e as Error).message}`)
    }
  }, [data, width, marginLeft, dark, spec.intron?.start, spec.intron?.end, shade?.[0], shade?.[1]]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <div>
      {error && <div className="px-2 py-1 text-xs text-error">{error}</div>}
      <div ref={host} className={data ? '' : 'min-h-[60px]'} />
    </div>
  )
}
