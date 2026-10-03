/**
 * How much the stored per-variant values were rounded, as text. Every number comes from the
 * experiment's per-results-set `precision` blocks (SPEC.md sections 8 and 13), so the site never
 * claims more accuracy than the builder recorded. The -log10 p, standard error and allele frequency
 * errors are worst-case bounds; the slope error is measured, by decoding every block again after
 * encoding and comparing each rebuilt slope to the source beta.
 *
 * Each error is the worse of the eQTL and sQTL results sets, rounded away from zero to two
 * significant figures. About.tsx spells them all out; components/rounding-note.tsx is the one-line
 * version beside the cis table. `slopeSe` is null when either set has no dof, which is also when no
 * slope is shown at all (SPEC section 8).
 */
import { useStoreInfo } from '@/contexts/store-context'
import type { PrecisionKind, StoreInfo } from '@/lib/store'

/** the `precision` keys that are plain numbers in both results sets, so the worse of the two is
 *  just a max (the slope error may be null, and is taken separately) */
type BoundedKey = Exclude<keyof PrecisionKind, 'slope_max_error_over_se'>

const SIG = 2
/** The exponent of x's leading digit: 0.0016 -> -3. */
const exp10 = (x: number) => Math.floor(Math.log10(x))
/** x rounded away from zero to SIG significant figures. */
const roundUp = (x: number) => { const f = 10 ** (SIG - 1 - exp10(x)); return Math.ceil(x * f) / f }
/** Plain decimals, never exponent notation: 7.7e-6 prints as 0.0000077. */
const dec = (x: number) => x.toFixed(Math.max(0, SIG - 1 - exp10(x)))
const fig = (x: number) => x > 0 ? dec(roundUp(x)) : '0'   // exp10(0) is -Infinity; a zero bound prints as 0

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
  const worst = (k: BoundedKey) => Math.max(p.eqtl[k], p.sqtl[k])
  // the p-value error follows from the rounded -log10 p, so the two numbers agree with each other
  const nlp = roundUp(worst('neglog10p_max_error'))
  const se = p.eqtl.slope_max_error_over_se, sq = p.sqtl.slope_max_error_over_se
  return {
    nlp: dec(nlp),
    pPct: fig((10 ** nlp - 1) * 100),
    sePct: fig(worst('slope_se_max_rel_error') * 100),
    slopeSe: se != null && sq != null ? fig(Math.max(se, sq)) : null,
    af: fig(p.af_max_error),
  }
}

export const useRoundingFacts = () => roundingFacts(useStoreInfo())

/** The note's own sentence without the link, for a `title` attribute. */
export const roundingText = (f: RoundingFacts) =>
  `p-values, slopes, standard errors, and allele frequencies are rounded (p within ${f.pPct}%). Exact values: Zenodo.`

/** What the cis note's tooltip spells out. The phrase around it carries the Zenodo link, so this
 *  says only how far each value can be off. */
export const roundingDetail = (f: RoundingFacts) =>
  `p-values are within ${f.pPct}% of the source, standard errors within ${f.sePct}%` +
  `${f.slopeSe ? `, slopes within ${f.slopeSe} standard errors` : ''}, and allele frequencies within ${f.af}.`

/** The same for a trans table, whose values are quantized on their own scales (SPEC section 9): the
 *  effect size has its own code rather than being rebuilt from a standard error, and the standard
 *  error and r² are derived from the two stored values, so they carry both errors. */
export interface TransRoundingFacts { pPct: string; beta: string; af: string }

export function transRoundingFacts(m: StoreInfo | null | undefined): TransRoundingFacts | null {
  const sets = [...(m?.results.values() ?? [])].map(r => r.trans?.precision).filter(p => p != null)
  if (!sets.length) return null
  const worst = (k: string) => Math.max(...sets.map(p => p[k] ?? 0))
  const nlp = roundUp(worst('neglog10p_max_error'))
  return { pPct: fig((10 ** nlp - 1) * 100), beta: fig(worst('beta_max_error')), af: fig(worst('af_max_error')) }
}

export const useTransRoundingFacts = () => transRoundingFacts(useStoreInfo())

export const transRoundingDetail = (f: TransRoundingFacts) =>
  `p-values are within ${f.pPct}% of the source, effect sizes within ${f.beta}, and allele frequencies ` +
  `within ${f.af}. Standard errors and r² are derived from the stored p-value and effect size.`
