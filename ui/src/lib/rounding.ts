/**
 * How much the stored per-variant values were rounded, as text. Every number comes from the
 * experiment's per-results-set `precision` blocks (SPEC.md sections 8 and 11), the worst-case
 * bounds the builder computed over every block, so the site never claims more accuracy than that.
 *
 * Each error is the worse of the eQTL and sQTL results sets, rounded away from zero to two
 * significant figures. About.tsx spells them all out; components/rounding-note.tsx is the one-line
 * version beside the cis table. v1 records no bound on the rebuilt slope, so `slopeSe` is null.
 */
import { useStoreInfo } from '@/contexts/store-context'
import type { PrecisionKind, StoreInfo } from '@/lib/store'

const SIG = 2
/** The exponent of x's leading digit: 0.0016 -> -3. */
const exp10 = (x: number) => Math.floor(Math.log10(x))
/** x rounded away from zero to SIG significant figures. */
const roundUp = (x: number) => { const f = 10 ** (SIG - 1 - exp10(x)); return Math.ceil(x * f) / f }
/** Plain decimals, never exponent notation: 7.7e-6 prints as 0.0000077. */
const dec = (x: number) => x.toFixed(Math.max(0, SIG - 1 - exp10(x)))
const fig = (x: number) => dec(roundUp(x))

/** Every rounding error as text, ready to drop into a sentence. */
export interface RoundingFacts {
  /** absolute error on -log10 p */
  nlp: string
  /** the p-value error that follows from it, in percent */
  pPct: string
  /** relative error on the standard error, in percent */
  sePct: string
  /** slope error, in standard errors (null: not recorded) */
  slopeSe: string | null
  /** absolute error on allele frequency */
  af: string
}

/** Null until the store has opened, or if the experiment lacks an eQTL or sQTL `precision` block. */
export function roundingFacts(m: StoreInfo | null | undefined): RoundingFacts | null {
  const p = m?.precision
  if (!p?.eqtl || !p.sqtl) return null
  const worst = (k: keyof PrecisionKind) => Math.max(p.eqtl[k], p.sqtl[k])
  // the p-value error follows from the rounded -log10 p, so the two numbers agree with each other
  const nlp = roundUp(worst('neglog10p_max_error'))
  return {
    nlp: dec(nlp),
    pPct: fig((10 ** nlp - 1) * 100),
    sePct: fig(worst('slope_se_max_rel_error') * 100),
    slopeSe: null,
    af: fig(p.af_max_error),
  }
}

export const useRoundingFacts = () => roundingFacts(useStoreInfo())

/** The note's own sentence without the link, for a `title` attribute. */
export const roundingText = (f: RoundingFacts) =>
  `p-values, slopes, standard errors, and allele frequencies are rounded (p within ${f.pPct}%). Exact values: Zenodo.`
