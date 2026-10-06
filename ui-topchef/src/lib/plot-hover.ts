import { useEffect, type RefObject } from 'react'
import type { Selection } from '@uwdata/mosaic-core'
import { DOT_R } from '@/lib/plot-theme'

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

/** What the hover overlay needs for one variant of the window, keyed by position. `cs` and
 *  `colocPp` are the dot's symbol and size channels, so the outline can take the point's own
 *  shape and size instead of a fixed circle that a large point swallows. */
export interface HoverRow {
  rs_number: number | null; nlp: number; gwas_nlp: number | null; label: string
  cs: string; colocPp: number | null
}
export type HoverLookup = (position: number) => HoverRow | undefined

/** The rendered Plot SVG exposes its scales; Mosaic uses the same accessor. Narrowed per scale at
 *  the call site: `r` is absent when the radius is a constant, `symbol` when nothing is shaped. */
type PlotSVG = SVGSVGElement & { scale?: (name: string) => unknown }
type NumScale = { apply(v: number): number }
/** A d3/Plot symbol: `draw` writes one symbol of the given *area* into a path context. */
type SymbolScale = { apply(v: string): { draw(context: PathSink, area: number): void } | undefined }

const SVG_NS = 'http://www.w3.org/2000/svg'
const MARK = 'hover-mark'
const FONT = 10   // px; the label's line height is 1em, as Plot's text mark draws it
/** How far outside the point's own edge the outline sits, in px of radius. */
const GAP = 3

/**
 * The handful of calls d3's symbol `draw` makes, collected as an SVG path string. It stands in for
 * d3-path, which is in the tree only as Plot's own dependency.
 *
 * `arc` is only ever `symbolCircle`'s full sweep (`moveTo(r, 0)` then 0 to tau), so it is emitted
 * as two half arcs and the angle arguments are ignored; `rect` is only ever `symbolSquare`. Writing
 * the path rather than cloning the hovered dot's own `<path>` out of the SVG keeps the overlay from
 * depending on how Plot orders and labels its marks in the DOM.
 */
class PathSink {
  private d = ''
  moveTo(x: number, y: number) { this.d += `M${x},${y}` }
  lineTo(x: number, y: number) { this.d += `L${x},${y}` }
  rect(x: number, y: number, w: number, h: number) { this.d += `M${x},${y}h${w}v${h}h${-w}Z` }
  arc(x: number, y: number, r: number) { this.d += `A${r},${r},0,1,1,${x - r},${y}A${r},${r},0,1,1,${x + r},${y}` }
  closePath() { this.d += 'Z' }
  toString() { return this.d }
}

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
      const sx = svg?.scale?.('x') as NumScale | undefined, sy = svg?.scale?.('y') as NumScale | undefined
      if (!svg || !sx || !sy || row == null || xv == null) { old?.remove(); return }
      const px = sx.apply(xv), py = sy.apply(row.nlp)
      if (!Number.isFinite(px) || !Number.isFinite(py)) { old?.remove(); return }

      // The outline takes the hovered dot's own geometry, so a large coloc point is circled rather
      // than hidden under a fixed ring. Both come from the plot's own scales: `r` is a channel only
      // when coloc sizes the points (`dotOptions`), and Plot draws a symbol at the area of a circle
      // of that radius (@observablehq/plot marks/dot.js: `S[i].draw(p, R[i] * R[i] * Math.PI)`).
      const sr = svg.scale?.('r') as NumScale | undefined
      const ssym = svg.scale?.('symbol') as SymbolScale | undefined
      const ringR = (sr ? sr.apply(row.colocPp ?? 0) : DOT_R) + GAP
      const sym = ssym?.apply(row.cs)
      const p = new PathSink()
      if (sym) sym.draw(p, ringR * ringR * Math.PI)
      else { p.moveTo(ringR, 0); p.arc(0, 0, ringR) }   // no symbol scale: the plain circle, as before

      const g = document.createElementNS(SVG_NS, 'g')
      g.setAttribute('class', MARK)
      g.setAttribute('pointer-events', 'none')
      // placed by transform and drawn around the origin, which is how Plot places the dots themselves
      const ring = document.createElementNS(SVG_NS, 'path')
      ring.setAttribute('transform', `translate(${px},${py})`); ring.setAttribute('d', String(p))
      ring.setAttribute('fill', 'none'); ring.setAttribute('stroke', ink); ring.setAttribute('stroke-width', '3')
      ring.setAttribute('stroke-linejoin', 'round')
      g.append(ring)
      // the label, laid out like Plot's text mark with textAnchor middle, lineAnchor bottom, and
      // lifted clear of the outline rather than the old fixed 12px
      const text = document.createElementNS(SVG_NS, 'text')
      text.setAttribute('class', 'hover-label')
      text.setAttribute('transform', `translate(${px},${py - ringR - 5})`)
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
