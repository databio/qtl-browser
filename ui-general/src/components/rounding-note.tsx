import ExternalLink from '@/components/ExternalLink'
import { Tooltip } from '@/components/tooltip'
import { ZENODO } from '@/lib/links'
import { roundingDetail, transRoundingDetail, useRoundingFacts, useTransRoundingFacts } from '@/lib/rounding'

/**
 * The rounding disclosure for a table of per-variant values, written to sit inside a section's
 * description sentence: a short phrase that carries the error bounds on hover, with the Zenodo link
 * left outside it and clickable, since a tooltip is plain text and cannot hold one.
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
    <>
      <Tooltip tip={detail}>
        {/* no rule under it: the app's other tooltip trigger (routes/Gene.tsx) is cursor-help and
            nothing else, and a line here would sit next to a real link and read as one */}
        <span className="cursor-help text-base-content/75">Values are rounded</span>
      </Tooltip>
      {' '}(exact values on <ExternalLink href={ZENODO}>Zenodo</ExternalLink>).
    </>
  )
}
