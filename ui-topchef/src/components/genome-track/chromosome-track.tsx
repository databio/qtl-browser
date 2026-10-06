import type { TrackLayout } from './use-track-layout'

export const CHROM_GAP = 2   // pixels between chromosome segments; the bin bars respect it too

/** Rounded segment per chromosome, alternating tints from the theme, label beneath. */
export function ChromosomeTrack({ layout, barY, barHeight }: { layout: TrackLayout; barY: number; barHeight: number }) {
  return (
    <g>
      {layout.visibleChroms.map(({ chr, startPx, widthPx }) => {
        const idx = layout.chroms.findIndex(c => c.name === chr)
        const segW = Math.max(widthPx - CHROM_GAP, 1)
        return (
          <g key={chr}>
            <rect x={startPx + CHROM_GAP / 2} y={barY} width={segW} height={barHeight} rx={barHeight / 2}
              className={idx % 2 === 0 ? 'fill-base-content/30' : 'fill-base-content/15'} />
            {widthPx > 6 && (
              <text x={startPx + widthPx / 2} y={barY + barHeight + 10} textAnchor="middle"
                className="fill-base-content/50" fontSize={8}>
                {chr.replace('chr', '')}
              </text>
            )}
          </g>
        )
      })}
    </g>
  )
}
