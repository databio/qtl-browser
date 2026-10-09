/**
 * How much the stored per-variant values were quantized, as text. Every number comes from the
 * experiments' `precision` blocks (SPEC.md section 13), so the site never claims more accuracy than
 * the builder recorded. Each is a worst-case bound, `theoretical_limit`, which any reader could
 * recompute from a block header; the builder also measured the actual worst error over every row,
 * and it sits below the bound in every experiment.
 *
 * The cis facts are the worse of the cis scans, the trans facts the worse of the trans scans, and
 * the GWAS its own, each rounded away from zero to two significant figures. About.tsx spells them
 * out; components/rounding-note.tsx is the one-line version beside a table.
 */
import { useStoreInfo } from '@/contexts/store-context'
import type { Bounds, StoreInfo } from '@/lib/store'

const SIG = 2
/** The exponent of x's leading digit: 0.0016 -> -3. */
const exp10 = (x: number) => Math.floor(Math.log10(x))
/** x rounded away from zero to SIG significant figures. */
const roundUp = (x: number) => { const f = 10 ** (SIG - 1 - exp10(x)); return Math.ceil(x * f) / f }
/** Plain decimals, never exponent notation: 7.7e-6 prints as 0.0000077. */
const dec = (x: number) => x.toFixed(Math.max(0, SIG - 1 - exp10(x)))
const fig = (x: number) => x > 0 ? dec(roundUp(x)) : '0'   // exp10(0) is -Infinity; a zero bound prints as 0

/** Every quantization error as text, ready to drop into a sentence. */
export interface RoundingFacts {
  /** absolute error on -log10 p */
  nlp: string
  /** the p-value error that follows from it, in percent */
  pPct: string
  /** error on beta, in units of the row's own standard error */
  betaSe: string
  /** relative error on the standard error, in percent */
  sePct: string
  /** absolute error on allele frequency */
  af: string
}

function facts(b: Bounds | null | undefined): RoundingFacts | null {
  if (!b) return null
  // the p-value error follows from the rounded -log10 p, so the two numbers agree with each other
  const nlp = roundUp(b.neglog10p)
  return { nlp: dec(nlp), pPct: fig((10 ** nlp - 1) * 100), betaSe: fig(b.beta_over_se), sePct: fig(b.se_rel * 100), af: fig(b.af) }
}

/** Null until the store has opened, or when the study has no scan of that kind. */
export const roundingFacts = (m: StoreInfo | null | undefined) => facts(m?.precision.cis)
export const transRoundingFacts = (m: StoreInfo | null | undefined) => facts(m?.precision.trans)
export const gwasRoundingFacts = (m: StoreInfo | null | undefined) => facts(m?.precision.gwas)

export const useRoundingFacts = () => roundingFacts(useStoreInfo())
export const useTransRoundingFacts = () => transRoundingFacts(useStoreInfo())

/** The note's own sentence, for a `title` attribute. */
export const roundingText = (f: RoundingFacts) =>
  `p-values, betas, standard errors, and allele frequencies are rounded (p within ${f.pPct}%). Exact values: Zenodo.`

/** What a note's tooltip spells out: how far each stored value can be from the source. The phrase
 *  around it carries the Zenodo link. */
export const roundingDetail = (f: RoundingFacts) =>
  `p-values are within ${f.pPct}% of the source, betas within ${f.betaSe} standard errors, standard errors within ` +
  `${f.sePct}%, and allele frequencies within ${f.af}. Exact values are on Zenodo.`

/** The same for a trans table; the trans scans carry their own bounds, and r² is formed from the
 *  stored beta and standard error. */
export const transRoundingDetail = (f: RoundingFacts) =>
  `p-values are within ${f.pPct}% of the source, betas within ${f.betaSe} standard errors, standard errors within ` +
  `${f.sePct}%, and allele frequencies within ${f.af}. r² is computed from the stored beta and standard error. Exact values are on Zenodo.`
