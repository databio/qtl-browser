import { useEffect, useRef, useState } from 'react'
import * as vg from '@uwdata/vgplot'
import { getDB } from '@/lib/db'
import { CS_COLORS, CS_DOMAIN, CS_SYMBOLS } from '@/lib/plot-theme'
import { linkedHoverMarks, PLOT_MARGIN_BOTTOM, PLOT_MARGIN_TOP, SURFACE } from '@/components/LocusPlot'
import type { Selection } from '@uwdata/mosaic-core'
import { clearPlotHover, onPlotPointerMove } from '@/lib/plot-hover'
import { Empty } from '@/components/states'

/**
 * LocusCompare: QTL −log10 p against DCM GWAS −log10 p for every variant of the window
 * present in both, from the already-materialized locus table. Same credible-set color and
 * shape encoding as the locus scatter, hover tooltips only. A colocalized locus streaks along
 * the diagonal; independent signals form an L.
 */
export default function LocusCompare({ table, dark, size = 320, yDomain, link, brush }: {
  table: string; dark: boolean; size?: number; yDomain: [number, number]; link: Selection
  /** the locus plot's brush: the panel shows only the brushed slice while one is drawn */
  brush: Selection
}) {
  const host = useRef<HTMLDivElement>(null)
  const [n, setN] = useState<number | null>(null)

  useEffect(() => {
    let alive = true
    const el = host.current!
    el.replaceChildren()
    setN(null)
    ;(async () => {
      const { con } = await getDB()
      const count = Number((await con.query(`SELECT count(*) AS n FROM ${table} WHERE gwas_p IS NOT NULL`)).toArray()[0].n)
      if (!alive) return
      setN(count)
      if (count === 0) return
      const colors = dark ? CS_COLORS.dark : CS_COLORS.light
      const ink = dark ? '#c3c2b7' : '#52514e'
      const plot = vg.plot(
        vg.dot(vg.from(table, { filterBy: brush }), {
          x: 'gwas_nlp', y: 'nlp', fill: 'cs', symbol: 'cs', r: 3.5,
          fillOpacity: vg.sql`CASE WHEN cs = 'none' THEN 0.35 ELSE 0.45 + 0.4 * pip END`,
          channels: { position: 'position' },
        }),
        ...linkedHoverMarks(table, link, 'gwas_nlp', 'nlp', ink, dark ? SURFACE.dark : SURFACE.light),
        vg.xLabel('DCM GWAS −log₁₀ p'), vg.yLabel('QTL −log₁₀ p'),
        vg.colorDomain([...CS_DOMAIN]), vg.colorRange(colors),
        vg.symbolDomain([...CS_DOMAIN]), vg.symbolRange(CS_SYMBOLS),
        vg.opacityDomain([0, 1]),   // see LocusPlot: fillOpacity is scaled, pin the domain
        // y axis identical to the locus scatter: same domain, height, margins, insets
        vg.xZero(true), vg.yDomain(yDomain), vg.xGrid(true), vg.yGrid(true), vg.xInset(8),
        vg.width(size), vg.height(size), vg.marginLeft(48), vg.marginBottom(PLOT_MARGIN_BOTTOM), vg.marginRight(20), vg.marginTop(PLOT_MARGIN_TOP),
        vg.style({ fontFamily: 'inherit', fontSize: '11px', color: ink, background: 'transparent' }),
      ) as HTMLElement
      if (!alive) return
      el.replaceChildren(plot)
    })().catch(e => console.error(e))
    return () => { alive = false; el.replaceChildren() }
  }, [table, dark, size, yDomain[0], yDomain[1], link, brush]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <div>
      <div ref={host} className="plot-host" onPointerMove={onPlotPointerMove} onPointerLeave={clearPlotHover} />
      {n === 0 && <Empty label="No variants in this window are present in the DCM GWAS." />}
    </div>
  )
}
