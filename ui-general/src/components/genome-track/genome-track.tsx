import { forwardRef, useCallback, useEffect, useImperativeHandle, useLayoutEffect, useMemo, useRef, useState } from 'react'
import { isCluster, type TrackBin, type TrackItem, type TrackLocus, type ViewState } from './types'
import { useTrackLayout } from './use-track-layout'
import { useGenomeZoom } from './use-genome-zoom'
import { ChromosomeTrack } from './chromosome-track'
import { LocusMarkers } from './locus-markers'
import { BinBars } from './bin-bars'
import { buildChromList, chromOffsets, clusterLoci, sortLociByPosition, toAbsolute } from '@/lib/genome-coords'

const LABEL_AREA = 16
const MARKER_AREA = 6
const GAP = 2
const BAR_Y = LABEL_AREA + MARKER_AREA + GAP
const BAR_HEIGHT = 4
const CHR_LABEL_AREA = 12
const ZOOM_PAD_TOP = 80          // extra interactive area above the content for easier pan/zoom
const CONTENT_HEIGHT = BAR_Y + BAR_HEIGHT + CHR_LABEL_AREA
const BIN_GAP = 4                // between chr labels and the downward bars
const BIN_AREA = 72              // height of the downward bar area when bins are given
const MIN_PIXEL_GAP = 8
const INSET_X = 8                // room for the half-width of a marker at either genome end

export type GenomeTrackHandle = {
  zoomToChrom: (chr: string) => void
  zoomToRegion: (chr: string, start: number, end: number) => void
  zoomIn: () => void
  zoomOut: () => void
  fullReset: () => void
  navigateLocus: (dir: 1 | -1) => void
  hasLoci: boolean
}

export type GenomeTrackProps = {
  loci: TrackLocus[]
  selectedLocusId?: string
  hoveredLocusId?: string
  onLocusSelect: (id: string | null) => void
  onViewChange?: (view: ViewState) => void
  /** Called with the number of loci dropped because their chromosome is unknown to the layout. */
  onSkipped?: (n: number) => void
  chromNames: string[]
  chromLengths: number[]
  traitColors?: Record<string, string>
  highlightRegion?: { chr: string; start: number; end: number } | null
  /** Whole-genome view only: no pan/zoom, no keyboard shortcuts, no auto-zoom on select. */
  static?: boolean
  /** Downward bars per genomic window beneath the chromosome bar (Miami-style lower half). */
  bins?: TrackBin[]
  binCap?: number
  binHeight?: number
  binThreshold?: number
  binUnit?: string
  className?: string
  /** data still on its way: draw dimmed with a highlight sweeping across; full opacity when false */
  loading?: boolean
}

/** Whole-genome bar with locus markers, d3 pan/zoom, keyboard navigation (ported from pegasus-v2f-ui). */
export const GenomeTrack = forwardRef<GenomeTrackHandle, GenomeTrackProps>(function GenomeTrack(
  { loci, selectedLocusId, hoveredLocusId, onLocusSelect, onViewChange: onViewChangeProp, onSkipped, chromNames, chromLengths, traitColors, highlightRegion,
    static: isStatic = false, bins, binCap = 20, binHeight = BIN_AREA, binThreshold, binUnit, className, loading = false }, ref,
) {
  // an empty `bins` array still reserves the bar area, so the track does not move when bins arrive
  const TOTAL_HEIGHT = ZOOM_PAD_TOP + CONTENT_HEIGHT + (bins ? BIN_GAP + binHeight : 0)
  const containerRef = useRef<HTMLDivElement>(null)
  const svgRef = useRef<SVGSVGElement>(null)
  // 0 until measured; the SVG is not drawn before then, so the first painted frame is already
  // at the real width (a default width here showed for one frame and then jumped)
  const [containerWidth, setContainerWidth] = useState(0)
  const [view, setView] = useState<ViewState>({ startBp: 0, endBp: 1 })

  useLayoutEffect(() => {
    const el = containerRef.current
    if (!el) return
    const w = el.getBoundingClientRect().width
    if (w > 0) setContainerWidth(w)
    const obs = new ResizeObserver(entries => {
      const w = entries[0]?.contentRect.width
      if (w && w > 0) setContainerWidth(w)
    })
    obs.observe(el)
    return () => obs.disconnect()
  }, [])

  // the genome maps onto [INSET_X, containerWidth - INSET_X] so end markers are not clipped
  const plotWidth = Math.max(1, containerWidth - 2 * INSET_X)
  const layout = useTrackLayout(plotWidth, view, chromNames, chromLengths)
  const onViewChange = useCallback((v: ViewState) => { setView(v); onViewChangeProp?.(v) }, [onViewChangeProp])
  const { zoomTo, resetZoom, zoomIn, zoomOut } = useGenomeZoom(svgRef, { totalLength: layout.totalLength, containerWidth: plotWidth, onViewChange, enabled: !isStatic })

  // Loci on chromosomes the layout does not know cannot be placed (toAbsolute throws), so
  // they are dropped here and the count is reported to the parent instead of only logged.
  const plottableLoci = useMemo(() => loci.filter(l => layout.offsets.has(l.chr)), [loci, layout.offsets])
  useEffect(() => { onSkipped?.(loci.length - plottableLoci.length) }, [loci.length, plottableLoci.length, onSkipped])

  const sortedLoci = useMemo(() => sortLociByPosition(plottableLoci, layout.offsets), [plottableLoci, layout.offsets])
  const trackItems = useMemo(() => clusterLoci(sortedLoci, layout.bpToPixel, MIN_PIXEL_GAP), [sortedLoci, layout.bpToPixel])

  // offsets from props, not from the view-dependent layout, so zoom animations cannot feed back
  const stableOffsets = useMemo(() => chromOffsets(buildChromList(chromNames, chromLengths)).offsets, [chromNames, chromLengths])

  const zoomToLocus = useCallback((locus: TrackLocus) => {
    if (!stableOffsets.has(locus.chr)) return
    const start = toAbsolute(locus.chr, locus.start, stableOffsets)
    const end = toAbsolute(locus.chr, locus.end, stableOffsets)
    const span = Math.max(end - start, 5_000_000)
    const mid = (start + end) / 2
    zoomTo(mid - span * 3, mid + span * 3)
  }, [stableOffsets, zoomTo])

  const navigateLocus = useCallback((direction: 1 | -1) => {
    if (sortedLoci.length === 0) return
    const cur = selectedLocusId ? sortedLoci.findIndex(l => l.id === selectedLocusId) : -1
    const next = direction === 1
      ? sortedLoci[cur < sortedLoci.length - 1 ? cur + 1 : 0]!
      : sortedLoci[cur > 0 ? cur - 1 : sortedLoci.length - 1]!
    onLocusSelect(next.id)
    zoomToLocus(next)
  }, [sortedLoci, selectedLocusId, onLocusSelect, zoomToLocus])

  const zoomToChrom = useCallback((chr: string) => {
    if (!chr) { resetZoom(); return }
    const offset = stableOffsets.get(chr)
    const info = buildChromList(chromNames, chromLengths).find(c => c.name === chr)
    if (offset === undefined || !info) return
    zoomTo(offset, offset + info.length)
  }, [stableOffsets, chromNames, chromLengths, zoomTo, resetZoom])

  const zoomToRegion = useCallback((chr: string, start: number, end: number) => {
    if (!stableOffsets.has(chr)) return
    zoomTo(toAbsolute(chr, start, stableOffsets), toAbsolute(chr, end, stableOffsets))
  }, [stableOffsets, zoomTo])

  const hitTest = useCallback((px: number, py: number): TrackItem | null => {
    // static: labels (and the hit column over each triangle) carry hover and click themselves
    if (isStatic) return null
    const contentY = py - ZOOM_PAD_TOP
    if (contentY < 0 || contentY > BAR_Y) return null
    let best: { item: TrackItem; dist: number } | null = null
    for (const item of trackItems) {
      const cx = isCluster(item) ? item.centerPixel : layout.bpToPixel(item.chr, (item.start + item.end) / 2)
      const dist = Math.abs(px - cx)
      if (dist < 7 && (!best || dist < best.dist)) best = { item, dist }
    }
    return best?.item ?? null
  }, [trackItems, layout, isStatic])

  const handleClick = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    const svg = svgRef.current
    if (!svg) return
    const rect = svg.getBoundingClientRect()
    const hit = hitTest(e.clientX - rect.left - INSET_X, e.clientY - rect.top)
    if (!hit) return
    if (isCluster(hit)) {
      if (isStatic) return   // nothing to zoom into; the linked table reveals cluster members on hover
      const first = hit.loci[0]!, last = hit.loci[hit.loci.length - 1]!
      if (!stableOffsets.has(first.chr)) return
      const absStart = toAbsolute(first.chr, first.start, stableOffsets)
      const absEnd = toAbsolute(last.chr, last.end, stableOffsets)
      const span = Math.max(absEnd - absStart, 5_000_000)
      const mid = (absStart + absEnd) / 2
      zoomTo(mid - span * 2, mid + span * 2)
    } else {
      onLocusSelect(hit.id)
    }
  }, [hitTest, onLocusSelect, stableOffsets, zoomTo, isStatic])

  const handleMouseMove = useCallback((e: React.MouseEvent<SVGSVGElement>) => {
    const svg = svgRef.current
    if (!svg) return
    const rect = svg.getBoundingClientRect()
    if (isStatic) return
    svg.style.cursor = hitTest(e.clientX - rect.left - INSET_X, e.clientY - rect.top) ? 'pointer' : 'grab'
  }, [hitTest, isStatic])

  // Auto-zoom to the selected locus; deselecting returns to the full genome. The zoom
  // function is read through a ref so the effect's dependency list is complete without
  // re-running on every render of zoomToLocus.
  const zoomToLocusRef = useRef(zoomToLocus)
  zoomToLocusRef.current = zoomToLocus
  const zoomedTo = useRef<{ id: string; width: number } | null>(null)
  useEffect(() => {
    if (isStatic) return
    if (!selectedLocusId) {
      if (zoomedTo.current !== null) resetZoom()
      zoomedTo.current = null
      return
    }
    if (zoomedTo.current?.id === selectedLocusId && zoomedTo.current.width === containerWidth) return
    const locus = sortedLoci.find(l => l.id === selectedLocusId)
    if (locus) {
      zoomToLocusRef.current(locus)
      zoomedTo.current = { id: selectedLocusId, width: containerWidth }
    }
  }, [selectedLocusId, sortedLoci, containerWidth, resetZoom, isStatic])

  const fullReset = useCallback(() => {
    zoomedTo.current = null
    onLocusSelect(null)
    resetZoom()
  }, [onLocusSelect, resetZoom])

  useEffect(() => {
    if (isStatic) return
    const handler = (e: KeyboardEvent) => {
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement || e.target instanceof HTMLSelectElement) return
      if (e.key === 'ArrowLeft') { e.preventDefault(); navigateLocus(-1) }
      else if (e.key === 'ArrowRight') { e.preventDefault(); navigateLocus(1) }
      else if (e.key === '=' || e.key === '+') { e.preventDefault(); zoomIn() }
      else if (e.key === '-') { e.preventDefault(); zoomOut() }
      else if (e.key === 'Escape') fullReset()
    }
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [navigateLocus, zoomIn, zoomOut, fullReset, isStatic])

  useImperativeHandle(ref, () => ({ zoomToChrom, zoomToRegion, zoomIn, zoomOut, fullReset, navigateLocus, hasLoci: sortedLoci.length > 0 }),
    [zoomToChrom, zoomToRegion, zoomIn, zoomOut, fullReset, navigateLocus, sortedLoci.length])

  const highlight = highlightRegion && stableOffsets.has(highlightRegion.chr) ? (() => {
    const a = layout.bpToPixel(highlightRegion.chr, highlightRegion.start)
    const b = layout.bpToPixel(highlightRegion.chr, highlightRegion.end)
    const lo = Math.max(0, Math.min(a, b)), hi = Math.min(plotWidth, Math.max(a, b))
    return hi > lo ? { lo, w: Math.max(2, hi - lo) } : null
  })() : null

  return (
    // inline style: the height holds the space on the unmeasured first frame, before the SVG exists
    // dimmed while its data loads, then the whole track goes to full opacity at once
    <div ref={containerRef} data-genome-track className={`relative transition-opacity duration-300 ${loading ? 'opacity-50' : 'opacity-100'} ${className ?? ''}`} style={{ height: TOTAL_HEIGHT }}>
      {containerWidth > 0 && <svg ref={svgRef} width={containerWidth} height={TOTAL_HEIGHT} className="select-none overflow-visible" onClick={handleClick} onMouseMove={handleMouseMove}>
        {/* hit surface first so labels (which opt back into pointer events) stay on top */}
        <rect width={containerWidth} height={TOTAL_HEIGHT} fill="transparent" />
        <g pointerEvents="none" transform={`translate(${INSET_X},${ZOOM_PAD_TOP})`}>
          {highlight && (
            <>
              {/* inline style: SVG blur filter has no utility class */}
              <rect x={highlight.lo} width={highlight.w} y={BAR_Y - 4} height={BAR_HEIGHT + 8} rx={4} fill="var(--color-primary)" fillOpacity={0.22} style={{ filter: 'blur(3px)' }} />
              <rect x={highlight.lo} width={highlight.w} y={BAR_Y - 1} height={BAR_HEIGHT + 2} rx={2} fill="var(--color-primary)" fillOpacity={0.35} />
            </>
          )}
          {bins && bins.length > 0 && (
            <BinBars bins={bins} layout={layout} topY={CONTENT_HEIGHT + BIN_GAP} height={binHeight} cap={binCap} threshold={binThreshold} unit={binUnit} />
          )}
          <ChromosomeTrack layout={layout} barY={BAR_Y} barHeight={BAR_HEIGHT} />
          {/* mounts when loci arrive, so the rise-in runs once as the data lands */}
          {trackItems.length > 0 && (
            <g className="animate-rise-in">
              <LocusMarkers items={trackItems} layout={layout} barY={BAR_Y} selectedLocusId={selectedLocusId} hoveredLocusId={hoveredLocusId} traitColors={traitColors}
                showLabels={isStatic} onLabelClick={id => onLocusSelect(id)} />
            </g>
          )}
        </g>
      </svg>}
    </div>
  )
})
