import ExternalLink from '@/components/ExternalLink'
import { ZENODO } from '@/lib/links'
import { useRoundingFacts } from '@/lib/rounding'

/** One line beside a table of rounded values, written from the experiment's `precision` blocks
 *  (lib/rounding.ts). Renders nothing until the store has opened. */
export default function RoundingNote({ className = '' }: { className?: string }) {
  const f = useRoundingFacts()
  if (!f) return null
  return (
    <span className={`text-xs text-base-content/60 ${className}`}>
      p-values, slopes, standard errors, and allele frequencies are rounded (p within {f.pPct}%). Exact values:{' '}
      <ExternalLink className="underline" href={ZENODO}>Zenodo</ExternalLink>.
    </span>
  )
}
