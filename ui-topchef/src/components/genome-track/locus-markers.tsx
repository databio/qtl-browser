import { useState } from 'react'
import type { TrackLayout } from './use-track-layout'
import { isCluster, type TrackItem, type TrackLocus } from './types'

const TRI_WIDTH = 7
const TRI_HEIGHT = 6
// always-on labels: the hovered one is magnified, everything else stays at base size
const LABEL_BASE = 7.5
const LABEL_HOVER = 14

/**
 * Downward triangles sitting on the chromosome bar. Colors are theme tokens: the series
 * color comes from `traitColors` (CSS color strings, e.g. `var(--color-primary)`); the
 * selected marker uses the text color. A cluster takes the union of its members' traits, so
 * mixed eQTL and sQTL members read as 'both'; only a cluster with no resolvable trait is muted.
 */
export function LocusMarkers({ items, layout, barY, selectedLocusId, hoveredLocusId, traitColors, showLabels = false, onLabelClick }: {
  items: TrackItem[]; layout: TrackLayout; barY: number
  selectedLocusId?: string; hoveredLocusId?: string; traitColors?: Record<string, string>
  /** Label every locus regardless of zoom; cluster members are labeled side by side, and
   *  hovering a label magnifies just that label. */
  showLabels?: boolean
  onLabelClick?: (id: string) => void
}) {
  const seriesColor = (trait?: string) => traitColors?.[trait ?? ''] ?? 'var(--color-primary)'
  const LABEL_STEP = 8 // horizontal footprint of one rotated label
  const [labelHover, setLabelHover] = useState<string | null>(null)
  const hov = labelHover ?? hoveredLocusId

  const hoverHandlers = (l: TrackLocus) => ({
    onMouseEnter: () => setLabelHover(l.id),
    onMouseLeave: () => setLabelHover((h: string | null) => (h === l.id ? null : h)),
    onClick: (e: React.MouseEvent) => { e.stopPropagation(); onLabelClick?.(l.id) },
  })

  /** An always-on, hoverable label. `size` is decided by the caller; `dim` fades a label
   *  whose cluster sibling is hovered so the magnified one reads cleanly. */
  const label = (l: TrackLocus, x: number, topY: number, size: number, front = false, dim = false) => (
    <text key={`${l.id}${front ? '-front' : ''}`} textAnchor="start" dominantBaseline="central" pointerEvents="all"
      className={`cursor-pointer ${l.id === selectedLocusId || l.id === hov ? 'fill-base-content' : 'fill-base-content/50'} ${dim ? 'opacity-25' : ''}`}
      fontSize={size} fontWeight={600}
      stroke="transparent" strokeWidth={4} /* invisible stroke widens the hit area around the glyphs */
      transform={`translate(${x}, ${topY - 3}) rotate(-90)`}
      {...hoverHandlers(l)}>
      {l.label}
    </text>
  )

  /** Invisible column behind a label reaching down over its triangle to the bar, so hovering
   *  the marker itself also magnifies the label and clicking it opens the gene. */
  const hitColumn = (l: TrackLocus, x: number, tipY: number) => {
    const labelLen = l.label.length * LABEL_BASE * 0.62 + 6
    const top = tipY - TRI_HEIGHT - 3 - labelLen
    return (
      <rect key={`${l.id}-hit`} x={x - LABEL_STEP / 2} y={top} width={LABEL_STEP} height={tipY - top}
        fill="transparent" pointerEvents="all" className="cursor-pointer" {...hoverHandlers(l)} />
    )
  }

  return (
    <g>
      {items.map((item, i) => {
        if (isCluster(item)) {
          const cx = item.centerPixel
          // a cluster's trait is the union of its members': one shared trait keeps it, and members
          // that cover eQTL and sQTL between them (or include 'both') make the cluster 'both'
          const traits = new Set(item.loci.map(l => l.trait).filter(Boolean))
          const trait = traits.size === 1 ? traits.values().next().value
            : traits.has('both') || (traits.has('eQTL') && traits.has('sQTL')) ? 'both' : undefined
          const color = trait ? seriesColor(trait) : 'var(--color-base-content)'
          const hovered = hov ? item.loci.find(l => l.id === hov) : undefined
          const grow = !!hovered && !showLabels   // static track: markers keep their size, only labels react
          const w = grow ? 9 : TRI_WIDTH, h = grow ? 8 : TRI_HEIGHT
          const tipY = barY - 2, topY = tipY - h
          return (
            <g key={`cluster-${i}`}>
              <path d={`M ${cx - w / 2} ${topY} L ${cx + w / 2} ${topY} L ${cx} ${tipY} Z`} fill={color} opacity={hovered ? 0.9 : trait ? 0.5 : 0.3} />
              {showLabels ? (() => {
                const hi = item.loci.findIndex(l => l.id === hov)
                const xOf = (j: number) => cx + (j - (item.loci.length - 1) / 2) * LABEL_STEP
                const sizeOf = (j: number) => (j === hi ? LABEL_HOVER : LABEL_BASE)
                return (
                  <>
                    {item.loci.map((l, j) => hitColumn(l, xOf(j), tipY))}
                    {item.loci.map((l, j) => label(l, xOf(j), topY, sizeOf(j), false, hi >= 0 && j !== hi))}
                    {/* re-draw the hovered label last so it sits above its neighbors */}
                    {hi >= 0 && label(item.loci[hi]!, xOf(hi), topY, LABEL_HOVER, true)}
                  </>
                )
              })() : hovered ? (
                <text textAnchor="start" dominantBaseline="central" className="fill-base-content/80" fontSize={8} fontWeight={600}
                  transform={`translate(${cx}, ${topY - 3}) rotate(-90)`}>{hovered.label}</text>
              ) : (
                <text x={cx} y={topY - 3} textAnchor="middle" className="fill-base-content/50" fontSize={7} fontWeight={600}>{item.count}</text>
              )}
            </g>
          )
        }
        const locus = item as TrackLocus
        const cx = layout.bpToPixel(locus.chr, (locus.start + locus.end) / 2)
        const isSelected = locus.id === selectedLocusId
        const isHovered = !isSelected && locus.id === hov
        const color = isSelected ? 'var(--color-base-content)' : seriesColor(locus.trait)
        const w = showLabels ? TRI_WIDTH : isSelected ? 10 : isHovered ? 9 : TRI_WIDTH
        const h = showLabels ? TRI_HEIGHT : isSelected ? 9 : isHovered ? 8 : TRI_HEIGHT
        const tipY = barY - 2, topY = tipY - h
        return (
          <g key={locus.id}>
            <path d={`M ${cx - w / 2} ${topY} L ${cx + w / 2} ${topY} L ${cx} ${tipY} Z`} fill={color} opacity={isSelected || isHovered ? 1 : 0.7} />
            {showLabels ? <>{hitColumn(locus, cx, tipY)}{label(locus, cx, topY, isHovered ? LABEL_HOVER : LABEL_BASE)}</> : (isSelected || isHovered || layout.bpPerPixel < 500_000) && (
              <text textAnchor="start" dominantBaseline="central"
                className={isSelected ? 'fill-base-content' : isHovered ? 'fill-base-content/80' : 'fill-base-content/50'}
                fontSize={isSelected ? 9 : isHovered ? 8 : 7} fontWeight={600}
                transform={`translate(${cx}, ${topY - 3}) rotate(-90)`}>
                {locus.label}
              </text>
            )}
          </g>
        )
      })}
    </g>
  )
}
