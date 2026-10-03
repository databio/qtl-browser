import { useState } from 'react'
import { ChevronLeft, ChevronRight, RotateCcw, ZoomIn, ZoomOut } from 'lucide-react'

const iconBtn = 'rounded p-1 transition-colors hover:bg-base-200 hover:text-base-content'

/** Compact navigation for the genome track: chromosome select, region input, zoom, locus stepper. */
export function TrackControls({ chromNames, onChromSelect, onRegionInput, onZoomIn, onZoomOut, onReset, onPrevLocus, onNextLocus, hasLoci }: {
  chromNames: string[]
  onChromSelect: (chr: string) => void
  onRegionInput: (chr: string, start: number, end: number) => void
  onZoomIn: () => void; onZoomOut: () => void; onReset: () => void
  onPrevLocus: () => void; onNextLocus: () => void; hasLoci: boolean
}) {
  const [regionText, setRegionText] = useState('')
  const submit = () => {
    const parsed = parseRegion(regionText)
    if (parsed) { onRegionInput(parsed.chr, parsed.start, parsed.end); setRegionText('') }
  }
  return (
    <div className="flex items-center gap-2 text-base-content/50">
      <div className="join">
        <select className="select select-bordered select-xs join-item w-24" onChange={e => onChromSelect(e.target.value)} defaultValue="">
          <option value="">All chr</option>
          {chromNames.map(n => <option key={n} value={n}>{n.replace('chr', 'Chr ')}</option>)}
        </select>
        <input type="text" className="input input-bordered input-xs join-item w-32 placeholder:text-base-content/30" placeholder="chr2:150M-160M"
          value={regionText} onChange={e => setRegionText(e.target.value)} onKeyDown={e => { if (e.key === 'Enter') submit() }} />
      </div>
      <div className="flex items-center">
        <button className={iconBtn} onClick={onZoomOut} title="Zoom out (-)"><ZoomOut className="size-3.5" /></button>
        <button className={iconBtn} onClick={onZoomIn} title="Zoom in (+)"><ZoomIn className="size-3.5" /></button>
        <button className={iconBtn} onClick={onReset} title="Reset (Esc)"><RotateCcw className="size-3.5" /></button>
      </div>
      {hasLoci && (
        <div className="flex items-center">
          <button className={iconBtn} onClick={onPrevLocus} title="Previous locus (←)"><ChevronLeft className="size-3.5" /></button>
          <button className={iconBtn} onClick={onNextLocus} title="Next locus (→)"><ChevronRight className="size-3.5" /></button>
        </div>
      )}
    </div>
  )
}

/** "chr2:150000000-160000000", "chr2:150M-160M", "2:150m-160m", "chrX:1k-5k". */
export function parseRegion(text: string): { chr: string; start: number; end: number } | null {
  const m = /^(?:chr)?([0-9]{1,2}|x|y)\s*:\s*([0-9.,]+\s*[mk]?)\s*[-–]\s*([0-9.,]+\s*[mk]?)$/i.exec(text.trim())
  if (!m) return null
  const chr = `chr${m[1]!.toUpperCase()}`
  const start = parseBp(m[2]!), end = parseBp(m[3]!)
  if (!Number.isFinite(start) || !Number.isFinite(end) || start >= end) return null
  return { chr, start, end }
}

function parseBp(s: string): number {
  const t = s.replace(/[,\s]/g, '').toLowerCase()
  if (t.endsWith('m')) return Math.round(parseFloat(t) * 1_000_000)
  if (t.endsWith('k')) return Math.round(parseFloat(t) * 1_000)
  return parseInt(t, 10)
}
