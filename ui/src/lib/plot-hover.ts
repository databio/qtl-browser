import { useEffect, type RefObject } from 'react'
import type { Selection } from '@uwdata/mosaic-core'
import { DOT_R } from '@/lib/plot-theme'

/**
 * Which linked plot is under the pointer. Decided geometrically on every pointer move rather
 * than with enter/leave events: the plots' overflow-visible labels and rings can hang over a
 * neighbor and keep enter/leave from firing where the eye expects.
 */
const CLASS = 'plot-hovered'

const within = (r: DOMRect, x: number, y: number) => x >= r.left && x <= r.right && y >= r.top && y <= r.bottom

/**
 * Forget which variant Mosaic's `nearest` interactor last published for a plot.
 *
 * It stashes that as `valueIndex` on the plot's `<svg>` (d3 binds `this` to the node) and skips the
 * update when the pointer lands on the same index again, but its own `pointerleave` handler clears
 * the *selection* and leaves the *index* alone (`@uwdata/mosaic-plot` `interactors/Nearest.js`).
 * So after the pointer leaves the plot -- to another window, or just to another part of the page --
 * and comes back onto the same variant, the interactor publishes nothing, no value event fires, and
 * the overlay never redraws: that one variant cannot be re-hovered until a different one has been.
 * Clearing the index wherever the pointer is not makes the next entry publish.
 */
function forgetNearest(host: Element, x?: number, y?: number) {
  const svg = host.querySelector('svg') as (SVGSVGElement & { valueIndex?: number }) | null
  if (svg && !(x !== undefined && y !== undefined && within(svg.getBoundingClientRect(), x, y))) svg.valueIndex = -1
}

export function onPlotPointerMove(e: PointerEvent | React.PointerEvent) {
  const x = e.clientX, y = e.clientY
  document.querySelectorAll<HTMLElement>('.plot-host').forEach(el => {
    el.classList.toggle(CLASS, within(el.getBoundingClientRect(), x, y))
    forgetNearest(el, x, y)
  })
}

export function clearPlotHover() {
  document.querySelectorAll<HTMLElement>('.plot-host').forEach(el => {
    el.classList.remove(CLASS)
    forgetNearest(el)
  })
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
const MARK = 'hover-mark'     // the ring and the label together
const LABEL = 'hover-label'   // the label alone; app.css hides it in the plot not under the pointer
const FONT = 10               // px; the label's line height is 1em, as Plot's text mark draws it
const HALO = 6                // px of stroke behind the label's glyphs, so it reaches 3px out
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
 * query and rebuilds no plot. This hook listens to the selection instead, so a hover costs a few
 * DOM updates. Neither piece is under the pointer (`pointer-events: none`) and overflow is
 * visible, so a label near an edge can hang past the plot.
 *
 * Both are SVG, appended to the rendered plot as one `<g>` and placed with the SVG's own scales.
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

      // The label, laid out like Plot's text mark with textAnchor middle and lineAnchor bottom,
      // lifted clear of the outline. Two <text> elements over the same lines: the first draws only
      // the halo, a round-jointed stroke of the surface colour, the second only the glyphs.
      //
      // One element with `paint-order: stroke` is the obvious way and is wrong: a browser paints a
      // line box's stroke and then its glyphs, line by line, so the halo of one line lands on top
      // of the line above it. Splitting the layers is the only way the halo is always underneath,
      // and SVG is the only place a text stroke can round its joins at all -- CSS has no
      // stroke-linejoin for `-webkit-text-stroke`, which miters and spikes at every corner.
      const lab = document.createElementNS(SVG_NS, 'g')
      lab.setAttribute('class', LABEL)
      lab.setAttribute('transform', `translate(${px},${py - ringR - 5})`)
      const lines = row.label.split('\n')
      const layer = (halo: boolean) => {
        const t = document.createElementNS(SVG_NS, 'text')
        t.setAttribute('text-anchor', 'middle'); t.setAttribute('font-size', String(FONT))
        if (halo) {
          t.setAttribute('fill', 'none'); t.setAttribute('stroke', surface)
          t.setAttribute('stroke-width', String(HALO)); t.setAttribute('stroke-linejoin', 'round')
        } else t.setAttribute('fill', ink)
        lines.forEach((line, i) => {
          const tspan = document.createElementNS(SVG_NS, 'tspan')
          tspan.setAttribute('x', '0')
          if (i === 0) tspan.setAttribute('y', `${1 - lines.length}em`)
          else tspan.setAttribute('dy', '1em')
          tspan.textContent = line
          t.append(tspan)
        })
        return t
      }
      lab.append(layer(true), layer(false))
      g.append(lab)
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
