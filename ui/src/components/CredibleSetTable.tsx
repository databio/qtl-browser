import { Fragment, useEffect, useState, type CSSProperties } from 'react'
import { Link } from 'react-router'
import { ChevronRight } from 'lucide-react'
import { Empty } from '@/components/states'
import { fmtInt, fmtNum } from '@/lib/format'
import type { CredibleSetRow } from '@/lib/queries'
import { csTint, useIsDark } from '@/lib/plot-theme'

/** One row per credible set, closed by default; expanding lists the member variants. */
export default function CredibleSetTable({ rows }: { rows: CredibleSetRow[] }) {
  const [open, setOpen] = useState<Set<number>>(new Set())
  const dark = useIsDark()
  useEffect(() => setOpen(new Set()), [rows])
  if (!rows.length) return <Empty label="No SuSiE credible set." />

  const sets = new Map<number, CredibleSetRow[]>()
  for (const r of rows) sets.set(r.cs_id, [...(sets.get(r.cs_id) ?? []), r])
  const toggle = (id: number) => setOpen(o => { const n = new Set(o); if (n.has(id)) n.delete(id); else n.add(id); return n })
  const variantLink = (r: CredibleSetRow) => (
    <Link className="link-quiet" to={`/variant/${r.rsid ?? `${r.chr}:${r.position}`}`}>
      {r.rsid ?? <span className="tabular-nums">{r.chr}:{fmtInt(r.position)}</span>}
    </Link>
  )

  return (
    <div className="overflow-x-auto rounded-lg border border-base-300">
      <table className="table table-sm">
        <thead>
          <tr><th className="w-8" /><th>Set</th><th className="text-right">Variants</th><th>Top variant</th><th className="text-right">Top PIP</th><th className="text-right">Σ PIP</th><th className="text-right">Span</th></tr>
        </thead>
        <tbody>
          {[...sets.entries()].sort((a, b) => a[0] - b[0]).map(([id, members]) => {
            const sorted = [...members].sort((a, b) => b.pip - a.pip)
            const top = sorted[0]
            const sum = members.reduce((s, r) => s + r.pip, 0)
            const lo = Math.min(...members.map(r => r.position)), hi = Math.max(...members.map(r => r.position))
            const isOpen = open.has(id)
            return (
              <Fragment key={id}>
                {/* inline style: the row tint is the set's plot color from the chart palette, not a theme token */}
                <tr className="cursor-pointer transition-colors bg-(--tint) hover:bg-(--tint-hover)" onClick={() => toggle(id)} aria-expanded={isOpen}
                  style={{ '--tint': csTint(id, dark), '--tint-hover': csTint(id, dark, true) } as CSSProperties}>
                  <td><ChevronRight className={`size-4 text-base-content/40 transition-transform ${isOpen ? 'rotate-90' : ''}`} /></td>
                  <td className="font-medium">Set {id}</td>
                  <td className="text-right tabular-nums">{members.length}</td>
                  <td onClick={e => e.stopPropagation()}>{variantLink(top)}</td>
                  <td className="text-right tabular-nums">{fmtNum(top.pip)}</td>
                  <td className="text-right tabular-nums text-base-content/60">{fmtNum(sum, 2)}</td>
                  <td className="text-right tabular-nums text-base-content/60">{fmtInt(hi - lo)} bp</td>
                </tr>
                {isOpen && (
                  <tr className="bg-base-200/40">
                    <td />
                    <td colSpan={6} className="p-0">
                      <table className="table table-xs">
                        <thead><tr><th>Position</th><th>rsID</th><th>A1/A2</th><th className="text-right">AF</th><th className="text-right">PIP</th></tr></thead>
                        <tbody>
                          {sorted.map((r, i) => (
                            <tr key={i}>
                              <td className="tabular-nums">{fmtInt(r.position)}</td>
                              <td>{variantLink(r)}</td>
                              <td className="font-mono text-xs text-base-content/60">{r.A1}/{r.A2}</td>
                              <td className="text-right tabular-nums text-base-content/60">{fmtNum(r.af)}</td>
                              <td className="text-right tabular-nums">{fmtNum(r.pip)}</td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </td>
                  </tr>
                )}
              </Fragment>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}
