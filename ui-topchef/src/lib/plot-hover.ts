import { useEffect, type RefObject } from 'react'
import type { Selection } from '@uwdata/mosaic-core'

/**
 * Which linked plot is under the pointer. Decided geometrically on every pointer move rather
 * than with enter/leave events: the plots' overflow-visible labels and rings can hang over a
 * neighbor and keep enter/leave from firing where the eye expects.
 */
const CLASS = 'plot-hovered'

export function onPlotPointerMove(e: PointerEvent | React.PointerEvent) {
  const x = e.clientX, y = e.clientY
  document.querySelectorAll<HTMLElement>('.plot-host').forEach(el => {
    const r = el.getBoundingClientRect()
    const inside = x >= r.left && x <= r.right && y >= r.top && y <= r.bottom
    el.classList.toggle(CLASS, inside)
  })
}

export function clearPlotHover() {
  document.querySelectorAll<HTMLElement>(`.plot-host.${CLASS}`).forEach(el => el.classList.remove(CLASS))
}

/** What the hover overlay needs for one variant of the window, keyed by position. */
export interface HoverRow { rs_number: number | null; nlp: number; gwas_nlp: number | null; label: string }
export type HoverLookup = (position: number) => HoverRow | undefined

/** The rendered Plot SVG exposes its scales; Mosaic uses the same accessor. */
type PlotSVG = SVGSVGElement & { scale?: (name: string) => { apply(v: number): number } | undefined }
const SVG_NS = 'http://www.w3.org/2000/svg'
const MARK = 'hover-mark'
const FONT = 10   // px; the label's line height is 1em, as Plot's text mark draws it

/**
 * Linked hover ring and label as a DOM overlay. The plots' `nearest` interactor publishes the
 * hovered variant's position into `link`; nothing in Mosaic consumes it, so a hover issues no
 * query and rebuilds no plot. This hook listens to the selection instead and appends one `<g>`
 * to the rendered SVG, placed with the SVG's own scales, so a hover costs a few DOM updates.
 * The `<g>` is not under the pointer (`pointer-events: none`) and overflow is visible, so a
 * label near an edge can hang past the plot as before.
 */
export function useHoverOverlay(host: RefObject<HTMLElement | null>, link: Selection | null, lookup: HoverLookup,
  xField: 'position' | 'gwas_nlp', ink: string, surface: string) {
  useEffect(() => {
    if (!link) return
    const draw = () => {
      const svg = host.current?.querySelector('svg') as PlotSVG | null
      const old = svg?.querySelector(`:scope > g.${MARK}`) ?? null
      const v = link.value
      const row = typeof v === 'number' ? lookup(v) : undefined
      const xv = row == null ? null : xField === 'position' ? (v as number) : row.gwas_nlp
      const sx = svg?.scale?.('x'), sy = svg?.scale?.('y')
      if (!svg || !sx || !sy || row == null || xv == null) { old?.remove(); return }
      const px = sx.apply(xv), py = sy.apply(row.nlp)
      if (!Number.isFinite(px) || !Number.isFinite(py)) { old?.remove(); return }

      const g = document.createElementNS(SVG_NS, 'g')
      g.setAttribute('class', MARK)
      g.setAttribute('pointer-events', 'none')
      const ring = document.createElementNS(SVG_NS, 'circle')
      ring.setAttribute('cx', String(px)); ring.setAttribute('cy', String(py)); ring.setAttribute('r', '5.5')
      ring.setAttribute('fill', 'none'); ring.setAttribute('stroke', ink); ring.setAttribute('stroke-width', '3')
      g.append(ring)
      // the label, laid out like Plot's text mark with textAnchor middle, lineAnchor bottom, dy -12
      const text = document.createElementNS(SVG_NS, 'text')
      text.setAttribute('class', 'hover-label')
      text.setAttribute('transform', `translate(${px},${py - 12})`)
      text.setAttribute('text-anchor', 'middle'); text.setAttribute('font-size', String(FONT))
      text.setAttribute('fill', ink); text.setAttribute('stroke', surface); text.setAttribute('stroke-width', '5')
      text.setAttribute('stroke-linejoin', 'round'); text.setAttribute('paint-order', 'stroke')
      const lines = row.label.split('\n')
      lines.forEach((line, i) => {
        const tspan = document.createElementNS(SVG_NS, 'tspan')
        tspan.setAttribute('x', '0')
        if (i === 0) tspan.setAttribute('y', `${1 - lines.length}em`)
        else tspan.setAttribute('dy', '1em')
        tspan.textContent = line
        text.append(tspan)
      })
      g.append(text)
      if (old) old.replaceWith(g)
      else svg.append(g)
    }
    link.addEventListener('value', draw)
    draw()
    return () => {
      link.removeEventListener('value', draw)
      host.current?.querySelector(`svg > g.${MARK}`)?.remove()
    }
  }, [host, link, lookup, xField, ink, surface])
}
