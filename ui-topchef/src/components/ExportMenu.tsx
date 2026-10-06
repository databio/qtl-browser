import { Download } from 'lucide-react'
import { exportPNG, exportSVG } from '@/lib/export-plot'

/** Small dropdown offering PNG and SVG exports of one or more plot containers. */
export default function ExportMenu({ targets, background, disabled = false }: {
  targets: { label: string; name: string; el: () => HTMLElement | null }[]
  background: string
  disabled?: boolean
}) {
  // -my-1: the 24 px button sits inside a 16 px description line without making that line taller
  if (disabled) {
    return (
      <button className="btn btn-xs -my-1 gap-1 rounded-lg border-base-300 font-medium" disabled title="Export plots">
        <Download className="size-3" /> Export
      </button>
    )
  }
  const run = (t: (typeof targets)[number], fmt: 'png' | 'svg') => {
    const el = t.el()
    if (!el) return
    if (fmt === 'svg') exportSVG(el, t.name, background)
    else exportPNG(el, t.name, background).catch(e => console.error(e))
    ;(document.activeElement as HTMLElement | null)?.blur()   // close the focus-driven dropdown
  }
  return (
    <div className="dropdown dropdown-end -my-1">
      <button tabIndex={0} className="btn btn-xs gap-1 rounded-lg border-base-300 font-medium" title="Export plots">
        <Download className="size-3" /> Export
      </button>
      <ul tabIndex={0} className="dropdown-content z-20 mt-1 w-max min-w-64 rounded-lg border border-base-300 bg-base-100 p-1 shadow-lg">
        {targets.map(t => (
          <li key={t.name} className="flex items-center justify-between gap-4 px-2 py-1 text-xs">
            <span className="whitespace-nowrap">{t.label}</span>
            <span className="flex shrink-0 gap-3">
              <a href="#" className="link-quiet" onClick={e => { e.preventDefault(); run(t, 'png') }}>PNG</a>
              <a href="#" className="link-quiet" onClick={e => { e.preventDefault(); run(t, 'svg') }}>SVG</a>
            </span>
          </li>
        ))}
      </ul>
    </div>
  )
}
