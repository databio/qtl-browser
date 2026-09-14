import { useCallback, useEffect, useRef, useState, type ReactNode } from 'react'
import { createPortal } from 'react-dom'
import * as vg from '@uwdata/vgplot'
import { Selection } from '@uwdata/mosaic-core'
import { dropTable, getCoordinator, getDB, lit, materialize, parquet } from '@/lib/db'
import { CS_COLORS, CS_DOMAIN, CS_SWATCH_CLIP, CS_SYMBOLS, isDark } from '@/lib/plot-theme'
import type { SearchHit } from '@/lib/queries'
import { CompareSkeleton, LocusSkeleton } from '@/components/plot-skeleton'
import { nominalFile, nominalRows, type CredibleSetRow, type Exon } from '@/lib/queries'
import GeneTrack from '@/components/GeneTrack'
import LocusCompare from '@/components/LocusCompare'
import { clearPlotHover, onPlotPointerMove, useHoverOverlay, type HoverLookup, type HoverRow } from '@/lib/plot-hover'
import { useOpenPath } from '@/lib/row-link'
import { rsFromNumber } from '@/lib/format'
import ExportMenu from '@/components/ExportMenu'

export interface LocusSpec {
  hit: SearchHit
  qtlType: 'e' | 's'
  phenotypeId?: string
  tss: number
  exons: Exon[]                              // collapsed model of the gene, from gene_detail
  intron?: { start: number; end: number }
}
const MARGIN_LEFT = 48
const SCATTER_H = 290    // height of the locus scatter; the LocusCompare square matches it
// localStorage: 'shown' when the gene track is toggled on; off by default. The key was renamed
// when the default flipped, so browsers that had 'shown' written under the old key start off.
const TRACK_KEY = 'topchef-gene-track-v2'
export const PLOT_MARGIN_TOP = 20
export const PLOT_MARGIN_BOTTOM = 36
// Brush magnifier (plans/2026-09-08-locus-brush-zoom.md), parked: with this off no brush
// interactor is added, so the Selection never gets a value and the popup, the gene-track
// shade, and the LocusCompare filter stay inert. Flip to true to bring it back.
const BRUSH_MAGNIFIER = false

/**
 * Linked hover for the locus scatter and the LocusCompare panel. A `nearest` interactor on
 * each plot publishes the hovered variant's position into one shared selection. No Mosaic mark
 * reads that selection: the ring and label are drawn by `useHoverOverlay` from the selection's
 * value event, so a hover issues no query and rebuilds no plot, and hovering a point in either
 * panel highlights and labels the same variant in both. Replaces Plot's built-in tip.
 */
export function hoverInteractor(link: Selection) {
  return vg.nearest({ as: link, channels: ['position'], fields: ['position'], maxRadius: 24 })
}
export const INK = { light: '#52514e', dark: '#c3c2b7' }
export const SURFACE = { light: '#ffffff', dark: '#1b1a1a' }

/** The position the linked hover selection currently holds (the `nearest` interactor's single
 *  field), mirrored into React state: null when no variant is under the pointer. */
export function useHoveredVariant(link: Selection | null): number | null {
  const [pos, setPos] = useState<number | null>(null)
  useEffect(() => {
    setPos(null)
    if (!link) return
    const read = () => { const v = link.value; setPos(typeof v === 'number' && Number.isFinite(v) ? v : null) }
    link.addEventListener('value', read)
    return () => link.removeEventListener('value', read)
  }, [link])
  return pos
}

/** Handlers and cursor class that make a plot host open the hovered variant's page on click:
 *  plain click navigates, cmd/ctrl/shift or middle click opens a tab, like the table rows.
 *
 *  The click is assembled from pointerdown and pointerup on the host rather than from the
 *  browser's click event: Mosaic swaps the whole SVG whenever it redraws (resize, theme, and
 *  the nearest interactor re-publishes on the first pointer event of each new SVG, pointerdown
 *  included), and a mousedown whose element is detached by mouseup never becomes a click. */
export interface VariantClick {
  className: string
  onPointerDown: (e: React.PointerEvent) => void
  onPointerUp: (e: React.PointerEvent) => void
}
const CLICK_SLOP = 4   // px of pointer travel between down and up beyond which it is a drag, not a click
export function useVariantClick(link: Selection | null, href: (position: number) => string | null): VariantClick {
  const open = useOpenPath()
  const hovered = useHoveredVariant(link)
  const press = useRef<{ x: number; y: number; pos: number; button: number } | null>(null)
  return {
    className: hovered !== null ? 'cursor-pointer' : '',
    onPointerDown: e => {
      const pos = link?.value
      press.current = typeof pos === 'number' ? { x: e.clientX, y: e.clientY, pos, button: e.button } : null
    },
    onPointerUp: e => {
      const p = press.current
      press.current = null
      if (!p || p.button !== e.button || Math.hypot(e.clientX - p.x, e.clientY - p.y) > CLICK_SLOP) return
      const to = href(p.pos)
      if (!to) return
      e.preventDefault()
      open(to, e)
    },
  }
}

/** One cis window as a table for the plots: -log10 p, credible-set class, a tooltip label, and
 *  the DCM GWAS statistics for variants present there (matched on position and alleles in
 *  either orientation, GWAS beta re-signed to the QTL effect allele A1). Ordered so
 *  credible-set variants are drawn last (on top). */
function locusSQL(spec: LocusSpec): string {
  const where = [`q.gene_id = ${lit(spec.hit.gene_id)}`, spec.phenotypeId ? `q.phenotype_id = ${lit(spec.phenotypeId)}` : null]
    .filter(Boolean).join(' AND ')
  const lo = spec.tss - 1_000_000, hi = spec.tss + 1_000_000
  return `
    SELECT q.position,
           -- p underflows to 0 for a few extreme variants: place them just above the largest finite value
           CASE WHEN q.pval_nominal = 0 THEN max(-log10(nullif(q.pval_nominal, 0))) OVER () * 1.05 ELSE -log10(q.pval_nominal) END AS nlp,
           q.pval_nominal = 0 AS clipped,
           q.pval_nominal, q.slope, q.slope_se, q.af, q.pip, q.cs_id, q.rs_number, q.A1, q.A2,
           q.tss_distance, q.ma_samples, q.ma_count,
           coalesce(q.cs_id::VARCHAR, 'none') AS cs,
           g.p AS gwas_p, -log10(g.p) AS gwas_nlp,
           CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END AS gwas_beta,
           coalesce('rs' || q.rs_number, q.position::VARCHAR) || '  ' || q.A1 || '/' || q.A2
             || chr(10) || CASE WHEN q.pval_nominal = 0 THEN 'p = 0 (underflow; drawn above the maximum)' ELSE 'p = ' || format('{:.2e}', q.pval_nominal) END
             || chr(10) || 'slope ' || format('{:.3f}', q.slope) || ' ± ' || format('{:.3f}', q.slope_se)
             || chr(10) || 'AF ' || format('{:.3f}', q.af)
             || CASE WHEN q.pip IS NULL THEN '' ELSE chr(10) || 'PIP ' || format('{:.3f}', q.pip) || ' (set ' || q.cs_id || ')' END
             || CASE WHEN g.p IS NULL THEN '' ELSE chr(10) || 'DCM GWAS p = ' || format('{:.2e}', g.p) || ', beta ' || format('{:+.3f}', CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END) || ' (A1 as effect allele)' END AS label
    FROM ${nominalRows([nominalFile(spec.hit, spec.qtlType)])} q
    LEFT JOIN (SELECT * FROM ${parquet(`gwas_dcm/chr=${spec.hit.chr}/data.parquet`)} WHERE position BETWEEN ${lo} AND ${hi}) g
      ON g.position = q.position AND ((g.ea = q.A1 AND g.nea = q.A2) OR (g.ea = q.A2 AND g.nea = q.A1))
    WHERE ${where}
    -- the GWAS lists some indels in both allele orientations as separate records: keep one per
    -- QTL variant, preferring the orientation that matches the QTL alleles as written
    QUALIFY row_number() OVER (PARTITION BY q.position, q.A1, q.A2 ORDER BY (g.ea = q.A1) DESC NULLS LAST, g.p) = 1
    ORDER BY q.cs_id IS NOT NULL, q.position`
}

/** -log10 p against position for one cis window. Dots colored by credible-set membership,
 *  PIP as opacity, a TSS rule, and Observable Plot's nearest-point tip. */
export default function LocusPlot({ spec, onCount, onLegend, onActions, onCredibleSets, onTable }: {
  spec: LocusSpec
  onCount?: (n: number) => void
  onLegend?: (sets: string[] | null) => void
  /** Receives the header controls (export menu, gene-track toggle) so the parent can place them. */
  onActions?: (actions: ReactNode | null) => void
  /** credible-set members of this locus, read from the materialized window (no extra fetch) */
  onCredibleSets?: (rows: CredibleSetRow[] | null) => void
  /** the materialized window's table name once it exists (null while loading or after it is
   *  dropped; `failed` when the materialize threw), so the cis table can page off it */
  onTable?: (name: string | null, failed?: boolean) => void
}) {
  const host = useRef<HTMLDivElement>(null)
  const [state, setState] = useState<'loading' | 'ready' | 'error'>('loading')
  // which locus the ready state belongs to: a new spec renders once before its effect flips
  // state to loading, and the gene track must not redraw for the new intron in that frame
  const [readyFor, setReadyFor] = useState<string | null>(null)
  const [width, setWidth] = useState(0)
  const [tableName, setTableName] = useState<string | null>(null)
  const [link, setLink] = useState<Selection | null>(null)
  // Transient magnifier. Dragging on the scatter draws a brush whose interval lands in
  // `brushSel`; a floating detail plot (x domain bound to the Selection), the LocusCompare
  // panel, and the gene-track shade all follow it through Mosaic while the mouse is down. On
  // release the brush and the Selection are reset and everything snaps back. React mirrors the
  // value into `brushed` only to mount the popup, shade the track, and mute the hover tooltip.
  const [brushSel, setBrushSel] = useState<Selection | null>(null)
  const [brushed, setBrushed] = useState<[number, number] | null>(null)
  // the interval interactor instance of the current plot, so release can clear the brush graphic
  const brushInteractor = useRef<{ reset(): void } | null>(null)
  // the scatter's screen rect at pointer-down: the popup is positioned from it
  const anchor = useRef<DOMRect | null>(null)
  const [yMax, setYMax] = useState(1)
  const [dark, setDark] = useState(isDark)
  // the window by position: what the hover overlay draws and what a click on a hovered dot
  // needs to resolve its page synchronously (a query at click time would fall outside the
  // user gesture for new tabs). Loaded once with the locus table.
  const hoverIndex = useRef<Map<number, HoverRow>>(new Map())
  const lookup = useCallback<HoverLookup>(pos => hoverIndex.current.get(pos), [])
  const variantHref = (pos: number) => {
    const row = hoverIndex.current.get(pos)
    return row ? `/variant/${row.rs_number != null ? rsFromNumber(row.rs_number) : `${spec.hit.chr}:${pos}`}` : null
  }
  const click = useVariantClick(link, variantHref)
  useHoverOverlay(host, link, lookup, 'position', dark ? INK.dark : INK.light, dark ? SURFACE.dark : SURFACE.light)
  // the gene track under the scatter is off by default; the choice is kept across pages
  const [showTrack, setShowTrack] = useState(() => localStorage.getItem(TRACK_KEY) === 'shown')
  useEffect(() => { localStorage.setItem(TRACK_KEY, showTrack ? 'shown' : 'hidden') }, [showTrack])

  // colors are baked into the SVG, so redraw when the theme flips
  useEffect(() => {
    const obs = new MutationObserver(() => setDark(isDark()))
    obs.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] })
    return () => obs.disconnect()
  }, [])

  // The left column's width drives the scatter width; both plots redraw on resize without
  // re-materializing the locus table (data and drawing are separate effects).
  const column = useRef<HTMLDivElement>(null)
  useEffect(() => {
    const el = column.current
    if (!el) return
    let t: number | null = null
    const obs = new ResizeObserver(entries => {
      const w = entries[0]?.contentRect.width
      if (!w) return
      if (t !== null) window.clearTimeout(t)
      t = window.setTimeout(() => setWidth(Math.max(320, Math.round(w))), 80)
    })
    obs.observe(el)
    return () => { obs.disconnect(); if (t !== null) window.clearTimeout(t) }
  }, [])

  // 1. data: materialize the window once per locus
  const key = `${spec.hit.gene_id}|${spec.qtlType}|${spec.phenotypeId ?? ''}`
  useEffect(() => {
    let alive = true
    let table: string | null = null
    setState('loading')
    setTableName(null)
    setLink(null)
    setBrushSel(null)
    hoverIndex.current = new Map()
    onLegend?.(null)
    onCredibleSets?.(null)
    onTable?.(null)
    ;(async () => {
      try {
        await getCoordinator()
        table = await materialize(locusSQL(spec))
        if (!alive) return
        onTable?.(table)
        const { con } = await getDB()
        const agg = (await con.query(`SELECT count(*) AS n, max(nlp) AS ymax FROM ${table}`)).toArray()[0]
        const sets = (await con.query(`SELECT DISTINCT cs FROM ${table} WHERE cs <> 'none' ORDER BY cs`)).toArray().map(r => String(r.cs))
        if (!alive) return
        onCount?.(Number(agg.n))
        onLegend?.(sets)
        hoverIndex.current = new Map((await con.query(`SELECT position, rs_number, nlp, gwas_nlp, label FROM ${table}`)).toArray()
          .map(r => [Number(r.position), {
            rs_number: r.rs_number == null ? null : Number(r.rs_number),
            nlp: Number(r.nlp), gwas_nlp: r.gwas_nlp == null ? null : Number(r.gwas_nlp), label: String(r.label),
          }]))
        if (!alive) return
        if (onCredibleSets) {
          // the window already holds every variant's set and PIP: the credible-set table comes
          // from it instead of a second range read of credible_sets.parquet
          const cs = (await con.query(`
            SELECT position, A1, A2, CASE WHEN rs_number IS NULL THEN NULL ELSE 'rs' || rs_number END AS rsid, af, cs_id, pip
            FROM ${table} WHERE cs_id IS NOT NULL ORDER BY cs_id, pip DESC`)).toArray()
          if (!alive) return
          onCredibleSets(cs.map(r => {
            const o: Record<string, unknown> = {}
            for (const [k, v] of Object.entries(r.toJSON())) o[k] = typeof v === 'bigint' ? Number(v) : v
            return { ...o, qtl_type: spec.qtlType, phenotype_id: spec.phenotypeId ?? spec.hit.gene_id, chr: spec.hit.chr } as CredibleSetRow
          }))
        }
        // one explicit y domain shared with the LocusCompare panel so the two y axes coincide
        setYMax(Math.max(1, Number(agg.ymax)) * 1.04)
        // empty: true → no hovered variant means the highlight layers draw nothing (an empty
        // selection otherwise means "no filter", which rings every point)
        setLink(Selection.single({ empty: true }))
        setBrushSel(Selection.intersect())
        setTableName(table)
      } catch (e) {
        console.error(e)
        if (alive) { setState('error'); onTable?.(null, true) }
      }
    })()
    return () => {
      alive = false
      onTable?.(null)
      if (table) dropTable(table)
    }
  }, [key]) // eslint-disable-line react-hooks/exhaustive-deps

  // mirror the brush Selection into React state, one update per frame at most
  useEffect(() => {
    setBrushed(null)
    if (!brushSel) return
    let raf: number | null = null
    const read = () => {
      raf = null
      const v = brushSel.value as [number, number] | undefined
      setBrushed(v && Number.isFinite(v[0]) && Number.isFinite(v[1]) ? [v[0], v[1]] : null)
    }
    const cb = () => { if (raf === null) raf = requestAnimationFrame(read) }
    brushSel.addEventListener('value', cb)
    return () => { brushSel.removeEventListener('value', cb); if (raf !== null) cancelAnimationFrame(raf) }
  }, [brushSel])

  // while brushed: drop the hover tooltip, and end the brush on release. The resets are
  // deferred a tick so d3-brush's own mouseup handling (which follows pointerup) runs first
  // and cannot re-publish the final extent after we cleared it.
  useEffect(() => {
    if (!brushed) return
    link?.reset()
    const end = () => setTimeout(() => { brushInteractor.current?.reset(); brushSel?.reset() }, 0)
    window.addEventListener('pointerup', end)
    window.addEventListener('pointercancel', end)
    return () => { window.removeEventListener('pointerup', end); window.removeEventListener('pointercancel', end) }
  }, [brushed !== null]) // eslint-disable-line react-hooks/exhaustive-deps

  // 2. drawing: redraw whenever the table, the width, or the theme changes
  useEffect(() => {
    const el = host.current
    if (!el || !tableName || !link || !brushSel || width === 0) return
    const colors = dark ? CS_COLORS.dark : CS_COLORS.light
    const ink = dark ? INK.dark : INK.light
    try {
      const plot = vg.plot(
        // drawn first so it sits under the dots and the hover label
        vg.ruleX([spec.tss], { stroke: ink, strokeOpacity: 0.6, strokeDasharray: '2,3' }),
        vg.dot(vg.from(tableName), {
          x: 'position', y: 'nlp', fill: 'cs', symbol: 'cs', r: 3.5,
          fillOpacity: vg.sql`CASE WHEN cs = 'none' THEN 0.35 ELSE 0.45 + 0.4 * pip END`,
          channels: { position: 'position' },
        }),
        // interactors bind to the mark added just before them: both the brush and the nearest
        // interactor below take their x field from the data dots
        ...(BRUSH_MAGNIFIER ? [vg.intervalX({ as: brushSel, brush: { fill: ink, fillOpacity: 0.08, stroke: ink, strokeOpacity: 0.5 } })] : []),
        hoverInteractor(link),
        vg.xDomain([spec.tss - 1_000_000, spec.tss + 1_000_000]), vg.yLabel('QTL −log₁₀ p'),
        vg.xLabel(`${spec.hit.chr} position (Mb)`), vg.xTickFormat((d: number) => (d / 1e6).toFixed(2)),
        vg.colorDomain([...CS_DOMAIN]), vg.colorRange(colors),
        vg.symbolDomain([...CS_DOMAIN]), vg.symbolRange(CS_SYMBOLS),
        // fillOpacity is a channel, so Plot scales it; without a fixed domain it stretches to
        // the data maximum and a gene with no credible sets draws every point fully opaque
        vg.opacityDomain([0, 1]),
        vg.xInset(8), vg.yDomain([0, yMax]), vg.yGrid(true),
        vg.width(width), vg.height(SCATTER_H), vg.marginLeft(MARGIN_LEFT), vg.marginRight(20), vg.marginTop(PLOT_MARGIN_TOP), vg.marginBottom(PLOT_MARGIN_BOTTOM),
        vg.style({ fontFamily: 'inherit', fontSize: '11px', color: ink, background: 'transparent' }),
      ) as HTMLElement
      el.replaceChildren(plot)
      // vgplot exposes the Plot instance on the element; keep its brush interactor for release
      const interactors = (plot as unknown as { value?: { interactors?: { selection: unknown; reset(): void }[] } }).value?.interactors ?? []
      brushInteractor.current = interactors.find(i => i.selection === brushSel) ?? null
      setState('ready')
      setReadyFor(key)
    } catch (e) {
      console.error(e)
      setState('error')
    }
    return () => { el.replaceChildren(); brushInteractor.current = null }
  }, [tableName, link, brushSel, width, dark, yMax, spec.tss, spec.hit.chr])

  const compareCol = useRef<HTMLDivElement>(null)
  const stem = `${spec.hit.symbol ?? spec.hit.gene_id}${spec.phenotypeId ? '_' + spec.phenotypeId.split(':').slice(0, 3).join('_') : ''}`
  useEffect(() => {
    // the buttons stay in place while a locus loads, disabled, so the header does not reflow
    onActions?.(
      <>
        {BRUSH_MAGNIFIER && <span className={`text-xs text-base-content/45 ${state === 'ready' ? '' : 'opacity-50'}`}>drag on the plot to magnify</span>}
        <label className={`inline-flex items-center gap-1.5 ${state === 'ready' ? 'cursor-pointer' : 'opacity-50'}`} title={showTrack ? 'Hide the gene track' : 'Show the gene track'}>
          <input type="checkbox" className="toggle toggle-xs" checked={showTrack} disabled={state !== 'ready'} onChange={e => setShowTrack(e.target.checked)} />
          Gene track
        </label>
        <ExportMenu disabled={state !== 'ready'} background={dark ? SURFACE.dark : SURFACE.light} targets={[
          { label: showTrack ? 'Locus plot with gene track' : 'Locus plot', name: `${stem}_locus`, el: () => column.current },
          { label: 'QTL versus GWAS', name: `${stem}_locuscompare`, el: () => compareCol.current },
        ]} />
      </>
    )
  }, [state, dark, stem, showTrack]) // eslint-disable-line react-hooks/exhaustive-deps

  // the popup sits above the scatter when there is room, else below it; it never takes the pointer
  const rect = anchor.current
  const popupStyle = rect ? {
    left: rect.left, width: rect.width,
    top: rect.top - DETAIL_H - 24 >= 8 ? rect.top - DETAIL_H - 20 : rect.bottom + 12,
  } : undefined

  return (
    // while brushing, every plot ignores the pointer so the nearest-variant tooltip stops
    // chasing the cursor; d3-brush tracks the drag through window listeners, so it is unaffected
    <div className={`flex flex-col gap-6 md:flex-row md:gap-2 ${brushed ? '[&_.plot-host_svg]:pointer-events-none' : ''}`}>
      {/* the minimum height holds room for the skeleton only; once drawn the column is as tall
          as its content, so a hidden gene track leaves no blank strip under the scatter */}
      <div ref={column} className={`relative min-w-0 flex-1 ${state === 'loading' ? 'min-h-[340px]' : ''}`}>
        {state === 'loading' && <LocusSkeleton chr={spec.hit.chr} />}
        {state === 'error' && <div className="p-4 text-sm text-error">Could not draw the locus.</div>}
        <div ref={host} className={`plot-host ${state === 'ready' ? '' : 'invisible'} ${click.className}`}
          onPointerMove={onPlotPointerMove} onPointerLeave={clearPlotHover}
          onPointerUp={click.onPointerUp}
          onPointerDown={e => { anchor.current = host.current?.getBoundingClientRect() ?? null; click.onPointerDown(e) }} />
        {BRUSH_MAGNIFIER && brushed && brushSel && tableName && state === 'ready' && readyFor === key && popupStyle && createPortal(
          <div className="pointer-events-none fixed z-50 rounded-lg border border-base-300 bg-base-100 p-1 shadow-lg" style={popupStyle}>
            <LocusDetail table={tableName} brush={brushSel} dark={dark} width={popupStyle.width - 8} yMax={yMax} tss={spec.tss} chr={spec.hit.chr} />
          </div>,
          document.body,
        )}
        {showTrack && state === 'ready' && readyFor === key && width > 0 && (
          <GeneTrack spec={{ chr: spec.hit.chr, geneId: spec.hit.gene_id, domain: [spec.tss - 1_000_000, spec.tss + 1_000_000], exons: spec.exons, intron: spec.intron }}
            width={width} marginLeft={MARGIN_LEFT} dark={dark} shade={brushed} />
        )}
      </div>
      {/* the right column is reserved from the start so the scatter measures its final width;
          the panel is a square the height of the scatter so the two plots share a top and bottom */}
      <div ref={compareCol} className="shrink-0" style={{ width: SCATTER_H }}>
        {state === 'ready' && tableName && link && brushSel
          ? <LocusCompare table={tableName} dark={dark} size={Math.min(SCATTER_H, Math.max(width, 200))} yDomain={[0, yMax]} link={link} brush={brushSel} click={click} lookup={lookup} />
          : <CompareSkeleton />}
      </div>
    </div>
  )
}

const DETAIL_H = 200

/** The brushed slice of the locus: the same encoding as the overview, reading the same table
 *  filtered by the brush Selection, with its x domain bound to that Selection so Mosaic
 *  re-queries and re-scales it as the brush moves. No hover: it exists only during the drag. */
function LocusDetail({ table, brush, dark, width, yMax, tss, chr }: {
  table: string; brush: Selection; dark: boolean; width: number; yMax: number; tss: number; chr: string
}) {
  const host = useRef<HTMLDivElement>(null)
  useEffect(() => {
    const el = host.current
    if (!el || width === 0) return
    const colors = dark ? CS_COLORS.dark : CS_COLORS.light
    const ink = dark ? '#c3c2b7' : '#52514e'
    try {
      const plot = vg.plot(
        vg.ruleX([tss], { stroke: ink, strokeOpacity: 0.6, strokeDasharray: '2,3' }),
        vg.dot(vg.from(table, { filterBy: brush }), {
          x: 'position', y: 'nlp', fill: 'cs', symbol: 'cs', r: 3.5,
          fillOpacity: vg.sql`CASE WHEN cs = 'none' THEN 0.35 ELSE 0.45 + 0.4 * pip END`,
          channels: { position: 'position' },
        }),
        vg.xDomain(brush), vg.yLabel('QTL −log₁₀ p'),
        vg.xLabel(`${chr} position (Mb), brushed region`), vg.xTickFormat((d: number) => (d / 1e6).toFixed(3)),
        vg.colorDomain([...CS_DOMAIN]), vg.colorRange(colors),
        vg.symbolDomain([...CS_DOMAIN]), vg.symbolRange(CS_SYMBOLS),
        vg.opacityDomain([0, 1]),
        vg.xInset(8), vg.yDomain([0, yMax]), vg.yGrid(true),
        vg.width(width), vg.height(DETAIL_H), vg.marginLeft(MARGIN_LEFT), vg.marginRight(20), vg.marginTop(PLOT_MARGIN_TOP), vg.marginBottom(PLOT_MARGIN_BOTTOM),
        vg.style({ fontFamily: 'inherit', fontSize: '11px', color: ink, background: 'transparent' }),
      ) as HTMLElement
      el.replaceChildren(plot)
    } catch (e) {
      console.error(e)
    }
    return () => { el.replaceChildren() }
  }, [table, brush, dark, width, yMax, tss, chr])
  return <div ref={host} />
}

/** Legend for the credible-set encoding; rendered by the parent so it can sit in the section header. */
export function LocusLegend({ sets }: { sets: string[] }) {
  const [dark, setDark] = useState(isDark)
  useEffect(() => {
    const obs = new MutationObserver(() => setDark(isDark()))
    obs.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] })
    return () => obs.disconnect()
  }, [])
  const colors = dark ? CS_COLORS.dark : CS_COLORS.light
  return (
    <div className="flex flex-wrap items-center justify-end gap-x-4 gap-y-1 text-xs text-base-content/60">
      {['none', ...sets].map(d => {
        const i = CS_DOMAIN.indexOf(d as (typeof CS_DOMAIN)[number])
        const shape = CS_SYMBOLS[i]
        return (
          <span key={d} className="inline-flex items-center gap-1.5">
            {/* inline style: swatch color and shape are data values from the chart palette, not theme tokens */}
            <span className={`inline-block size-2.5 ${shape === 'circle' ? 'rounded-full' : ''}`}
              style={{ backgroundColor: colors[i], clipPath: CS_SWATCH_CLIP[shape] }} />
            {d === 'none' ? 'not in a credible set' : `credible set ${d}`}
          </span>
        )
      })}
      <span className="inline-flex items-center gap-1.5"><span className="inline-block h-3 border-l border-dashed border-base-content/60" /> TSS</span>
    </div>
  )
}
