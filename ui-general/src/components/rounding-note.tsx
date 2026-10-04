import { Tooltip } from '@/components/tooltip'
import { roundingDetail, transRoundingDetail, useRoundingFacts, useTransRoundingFacts } from '@/lib/rounding'

/**
 * The rounding disclosure for a table of per-variant values, written to sit inside a section's
 * description sentence. Only the phrase shows; the error bounds and where the exact values are kept
 * come on hover, so the sentence it sits in stays short.
 *
 * `kind` picks which bounds the tooltip quotes: the cis results sets' (`precision`) or the trans
 * objects' (`trans.precision`), which are quantized on their own scales (lib/rounding.ts).
 * Renders nothing until the store has opened, so the sentence around it must read without it.
 */
export default function RoundingNote({ kind = 'cis' }: { kind?: 'cis' | 'trans' }) {
  const cis = useRoundingFacts()
  const trans = useTransRoundingFacts()
  const detail = kind === 'trans' ? (trans && transRoundingDetail(trans)) : (cis && roundingDetail(cis))
  if (!detail) return null
  return (
    <Tooltip tip={detail}>
      {/* no rule under it: the app's other tooltip trigger (routes/Gene.tsx) is cursor-help and
          nothing else */}
      <span className="cursor-help text-base-content/75">Values are rounded.</span>
    </Tooltip>
  )
}
