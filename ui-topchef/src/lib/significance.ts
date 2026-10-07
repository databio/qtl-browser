import type { ExperimentDoc } from '@/lib/store'

/**
 * How the experiment decided what counts as a signal, resolved once for the whole UI.
 *
 * The rule lives in the experiment pointer (SPEC section 8), so it differs per study: TOPCHeF and
 * GTEx both test tensorQTL's permutation p, while a study whose source ran PLINK2 `--glm` or
 * MatrixEQTL has no permutation p at all and tests something else, or nothing. Every user-facing
 * string about significance comes from here, because a literal "permutation p < 0.05" is wrong the
 * moment a second kind of study is in the store -- and wrong in a way that reads as a result.
 *
 * `assessed` false means the source assessed no significance. That is not "nothing was
 * significant": the badges, counts and filters are hidden rather than shown as zero, and the index
 * carries null rather than false for every phenotype.
 */
export interface SignificanceRule { column: string; op: string; threshold: number; label?: string }

export interface Significance {
  /** the experiment declared a rule, so there is something to show */
  assessed: boolean
  /** the metric's display name, e.g. 'Permutation p' */
  label: string
  /** a column-header-sized version, e.g. 'Perm p' */
  short: string
  /** the rule as a phrase for a tooltip, e.g. 'permutation p < 0.05' */
  ruleText: string
}

/** Fallback names for the metrics we ship adapters for, used when a rule carries no `label`.
 *  A rule written by a newer adapter brings its own, so this list does not have to grow. */
const LABELS: Record<string, string> = {
  p_perm: 'Permutation p',
  p_beta: 'Beta-approximated p',
  p_bonferroni: 'Bonferroni p',
  fdr: 'FDR',
}

/** 'Permutation p' -> 'Perm p': the first word abbreviated, since table headers are narrow. */
const SHORT: Record<string, string> = {
  p_perm: 'Perm p',
  p_beta: 'Beta p',
  p_bonferroni: 'Bonf p',
  fdr: 'FDR',
}

const NOT_ASSESSED: Significance = { assessed: false, label: '', short: '', ruleText: '' }

/** The experiment's rule, resolved. `assessed` is false for a null rule and also before the store
 *  has opened, so a component renders no badge rather than a wrong one on the first frame.
 *  `useSignificance` in contexts/store-context.tsx is how a component gets it. */
export function significanceOf(exp: Pick<ExperimentDoc, 'significance'> | null | undefined): Significance {
  const rule = exp?.significance
  if (!rule || !rule.column) return NOT_ASSESSED
  const label = rule.label || LABELS[rule.column] || rule.column
  return {
    assessed: true,
    label,
    short: SHORT[rule.column] || label,
    // lower-cased so it reads inside a sentence: 'Significant cis-eQTL: permutation p < 0.05'
    ruleText: `${label.charAt(0).toLowerCase()}${label.slice(1)} ${rule.op} ${rule.threshold}`,
  }
}
