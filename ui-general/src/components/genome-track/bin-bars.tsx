import type { TrackLayout } from './use-track-layout'
import type { TrackBin } from './types'
import { CHROM_GAP } from './chromosome-track'

/**
 * Downward bars, one per genomic window, hanging from `topY`. Values above `cap` are
 * drawn at full height with a notch so the reader knows they are clipped; the exact value
 * stays in the hover readout. `threshold` draws a dashed reference line in data units.
 */
export function BinBars({ bins, layout, topY, height, cap, threshold, unit }: {
  bins: TrackBin[]; layout: TrackLayout; topY: number; height: number; cap: number; threshold?: number; unit?: string
}) {
  const scale = (v: number) => (Math.min(v, cap) / cap) * height
  const last = layout.chroms[layout.chroms.length - 1]
  const xEnd = last ? layout.bpToPixel(last.name, last.length) : 0
  return (
    <g>
      {threshold != null && threshold < cap && (
        <line x1={0} x2={xEnd} y1={topY + scale(threshold)} y2={topY + scale(threshold)}
          className="stroke-error/50" strokeWidth={0.75} strokeDasharray="3,3" />
      )}
      {bins.map((b, i) => {
        const seg = layout.visibleChroms.find(c => c.chr === b.chr)
        if (!seg) return null
        // bars live inside their chromosome's drawn segment (same 2 px gap as the chromosome
        // bar), all at one pixel width: the window size in pixels, rounded, minus a 1 px gap.
        // Each bar is placed by its window index from the segment start, so two windows can
        // never round onto the same pixel. A window past the segment's end is not drawn.
        const segStart = Math.ceil(seg.startPx + CHROM_GAP / 2), segEnd = Math.floor(seg.startPx + seg.widthPx - CHROM_GAP / 2)
        const binPx = (b.end - b.start) / layout.bpPerPixel
        const nominal = Math.max(1.25, Math.round(binPx) - 1)   // never thinner than 1.25 px
        const x0 = segStart + Math.round((b.start / (b.end - b.start)) * binPx)
        const w = Math.min(nominal, segEnd - x0)
        if (w < 1) return null
        const h = scale(b.value)
        const clipped = b.value > cap
        return (
          // grows down from the baseline on mount, sweeping left to right across the genome
          // (inline style: the per-bar delay is computed from its position)
          <g key={i} pointerEvents="all" className="animate-grow-down origin-top [transform-box:fill-box]"
            style={{ animationDelay: `${xEnd ? (x0 / xEnd) * 250 : 0}ms` }}>
            <rect x={x0} y={topY} width={w} height={Math.max(h, 0.5)} rx={0.5}
              className={b.sig ? 'fill-error' : 'fill-base-content/20'} />
            {clipped && <rect x={x0} y={topY + height - 2} width={w} height={2} className="fill-base-100" />}
            {b.title && <title>{b.title}{clipped ? ` (bar clipped at ${cap}${unit ? ` ${unit}` : ''})` : ''}</title>}
          </g>
        )
      })}
    </g>
  )
}
