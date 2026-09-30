import { useCallback, useEffect, useRef, type RefObject } from 'react'
import { zoom, zoomIdentity, type ZoomBehavior } from 'd3-zoom'
import { select } from 'd3-selection'
import 'd3-transition'
import type { ViewState } from './types'

type Options = { totalLength: number; containerWidth: number; onViewChange: (view: ViewState) => void; enabled?: boolean }
export type ZoomControls = {
  zoomTo: (startBp: number, endBp: number) => void
  resetZoom: () => void
  zoomIn: () => void
  zoomOut: () => void
}

/**
 * Bridges d3-zoom with React: d3 only tracks the transform; every zoom event becomes a
 * ViewState so React re-renders at the right positions. k=1, tx=0 is the whole genome.
 */
export function useGenomeZoom(svgRef: RefObject<SVGSVGElement | null>, { totalLength, containerWidth, onViewChange, enabled = true }: Options): ZoomControls {
  const zoomRef = useRef<ZoomBehavior<SVGSVGElement, unknown> | null>(null)
  const totalRef = useRef(totalLength)
  const widthRef = useRef(containerWidth)
  const viewChangeRef = useRef(onViewChange)
  const rafRef = useRef<number | null>(null)
  const pendingView = useRef<ViewState | null>(null)
  totalRef.current = totalLength
  widthRef.current = containerWidth
  viewChangeRef.current = onViewChange

  const transformToView = useCallback((t: { k: number; x: number }): ViewState => {
    const total = totalRef.current, w = widthRef.current
    const bpPerPx = total / (w * t.k)
    const startBp = -t.x * bpPerPx
    return { startBp: Math.max(0, startBp), endBp: Math.min(total, startBp + total / t.k) }
  }, [])

  useEffect(() => {
    const svg = svgRef.current
    if (!svg || totalLength === 0 || containerWidth === 0) return
    if (!enabled) {
      // static track: whole genome, no d3 behavior attached
      zoomRef.current = null
      viewChangeRef.current({ startBp: 0, endBp: totalLength })
      return
    }
    const maxZoom = Math.max(totalLength / 500_000, 10)
    const zb = zoom<SVGSVGElement, unknown>()
      .scaleExtent([1, maxZoom])
      .translateExtent([[0, 0], [containerWidth, 0]])
      .extent([[0, 0], [containerWidth, 0]])
      .filter(event => {
        if (event.type === 'dblclick') return false
        // horizontal-dominant wheel events are handled below as pan
        if (event.type === 'wheel' && Math.abs(event.deltaX) > Math.abs(event.deltaY)) return false
        return true
      })
      .on('zoom', event => {
        pendingView.current = transformToView(event.transform)
        if (rafRef.current === null) {
          rafRef.current = requestAnimationFrame(() => {
            rafRef.current = null
            viewChangeRef.current(pendingView.current!)
          })
        }
      })
    const sel = select(svg)
    sel.call(zb)
    sel.on('dblclick.zoom', null)
    zoomRef.current = zb
    const onWheel = (e: WheelEvent) => {
      if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
        e.preventDefault()
        sel.call(zb.translateBy, -e.deltaX / 2, 0)
      }
    }
    svg.addEventListener('wheel', onWheel, { passive: false })
    viewChangeRef.current({ startBp: 0, endBp: totalLength })
    return () => {
      sel.on('.zoom', null)
      svg.removeEventListener('wheel', onWheel)
      if (rafRef.current !== null) { cancelAnimationFrame(rafRef.current); rafRef.current = null }
    }
  }, [svgRef, totalLength, containerWidth, transformToView, enabled])

  const zoomTo = useCallback((startBp: number, endBp: number) => {
    const svg = svgRef.current, zb = zoomRef.current
    if (!svg || !zb) return
    const total = totalRef.current, w = widthRef.current
    const k = total / (endBp - startBp)
    const tx = -(startBp * w * k) / total
    select(svg).transition().duration(650).call(zb.transform, zoomIdentity.translate(tx, 0).scale(k))
  }, [svgRef])

  const resetZoom = useCallback(() => {
    const svg = svgRef.current, zb = zoomRef.current
    if (!svg || !zb) return
    select(svg).transition().duration(650).call(zb.transform, zoomIdentity)
  }, [svgRef])

  const zoomIn = useCallback(() => {
    const svg = svgRef.current, zb = zoomRef.current
    if (svg && zb) select(svg).transition().duration(300).call(zb.scaleBy, 2)
  }, [svgRef])

  const zoomOut = useCallback(() => {
    const svg = svgRef.current, zb = zoomRef.current
    if (svg && zb) select(svg).transition().duration(300).call(zb.scaleBy, 0.5)
  }, [svgRef])

  return { zoomTo, resetZoom, zoomIn, zoomOut }
}
