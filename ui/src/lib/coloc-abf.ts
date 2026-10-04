/**
 * coloc.abf (Giambartolomei et al. 2014) over a locus, from values the qtlstore already carries.
 *
 * Nothing here needs a new pipeline output. Per variant the approximate Bayes factor needs only a
 * z statistic and a standard error, and both traits store those:
 *
 *   QTL   `z = sign * t(p, dof)`, `se` -- the SE code decodes directly and `tFromNlp` already
 *         exists for the slope, so the block's two u16 codes are exactly the inputs. The slope
 *         itself, the one value in the format that is reconstructed rather than decoded, is not on
 *         this path at all.
 *   GWAS  `z = beta / se`, both lossless at the source's printed precision (SPEC section 10).
 *
 * No LD matrix: ABF assumes a single causal variant in the region and never conditions, which is
 * the whole reason this is feasible in a browser. The multi-causal-variant version (coloc.susie)
 * needs both traits' SuSiE alpha/lbf matrices, and the store keeps only credible-set membership and
 * PIP, so it is not reconstructible from a qtlstore.
 *
 * Everything is computed in log space with a logsumexp, because the ABFs at a real locus span
 * hundreds of orders of magnitude.
 */

/** Prior variance on the effect size, and the three configuration priors. coloc's defaults:
 *  W 0.15^2 for a quantitative trait, 0.2^2 on the log-odds scale for a case-control one. */
export interface AbfPriors {
  /** prior variance for trait 1 (the QTL) */
  W1: number
  /** prior variance for trait 2 (the GWAS) */
  W2: number
  /** prior that a variant is causal for trait 1 only */
  p1: number
  /** prior for trait 2 only */
  p2: number
  /** prior that one variant is causal for both */
  p12: number
}

export const DEFAULT_PRIORS: AbfPriors = { W1: 0.15 ** 2, W2: 0.2 ** 2, p1: 1e-4, p2: 1e-4, p12: 1e-5 }

export interface ColocResult {
  /** variants that contributed (both traits present, both finite) */
  n: number
  /** posteriors for H0 (neither), H1 (trait 1 only), H2 (trait 2 only), H3 (both, distinct causal
   *  variants), H4 (both, shared causal variant) */
  pp: [number, number, number, number, number]
  /** per-variant posterior that it is the shared causal variant, given H4; sums to 1 */
  snpPP4: Float64Array
  /** index into the contributing variants of the largest `snpPP4` */
  best: number
}

/**
 * Wakefield's log approximate Bayes factor for one variant.
 *
 * `r` is the shrinkage, the share of the total variance that is prior rather than sampling. The
 * log is of `1 - r`, not `r`, so a variant with a tiny standard error (r near 1) gets a large
 * negative first term and the evidence has to come from z^2.
 */
export function logAbf(z: number, se: number, W: number): number {
  const r = W / (W + se * se)
  return 0.5 * (Math.log1p(-r) + r * z * z)
}

/** log(sum(exp(x))) without overflowing: the ABFs at a strong locus reach e^300. */
export function logSumExp(x: ArrayLike<number>): number {
  let m = -Infinity
  for (let i = 0; i < x.length; i++) if (x[i]! > m) m = x[i]!
  if (!Number.isFinite(m)) return m           // all -Infinity (no variants), or an Infinity to propagate
  let s = 0
  for (let i = 0; i < x.length; i++) s += Math.exp(x[i]! - m)
  return m + Math.log(s)
}

/**
 * The five posteriors over a locus. `z1/se1` is one trait and `z2/se2` the other, aligned variant
 * for variant: element i of all four arrays must be the same site.
 *
 * H3 looks like it needs every pair of variants and does not:
 * `sum_{i != j} ABF1_i ABF2_j = (sum ABF1)(sum ABF2) - sum_i ABF1_i ABF2_i`, so the whole thing is
 * linear in the number of variants.
 */
export function colocAbf(z1: ArrayLike<number>, se1: ArrayLike<number>, z2: ArrayLike<number>,
  se2: ArrayLike<number>, priors: AbfPriors = DEFAULT_PRIORS): ColocResult {
  const n = z1.length
  if (se1.length !== n || z2.length !== n || se2.length !== n)
    throw new Error(`coloc.abf: arrays of ${n}, ${se1.length}, ${z2.length}, ${se2.length}; they must be aligned`)
  if (!n) throw new Error('coloc.abf: no variants in the locus')
  const l1 = new Float64Array(n), l2 = new Float64Array(n), l12 = new Float64Array(n)
  for (let i = 0; i < n; i++) {
    l1[i] = logAbf(z1[i]!, se1[i]!, priors.W1)
    l2[i] = logAbf(z2[i]!, se2[i]!, priors.W2)
    l12[i] = l1[i]! + l2[i]!
  }
  const S1 = logSumExp(l1), S2 = logSumExp(l2), S12 = logSumExp(l12)

  // log(exp(S1 + S2) - exp(S12)), the subtraction done against the larger term. S1 + S2 >= S12
  // always, by Cauchy-Schwarz on the ABF vectors, so the difference is never negative; a zero here
  // means a single variant carries the whole locus for both traits, and H3 is then genuinely 0.
  const hi = Math.max(S1 + S2, S12)
  const diff = Math.exp(S1 + S2 - hi) - Math.exp(S12 - hi)
  const lH3 = diff > 0 ? hi + Math.log(diff) : -Infinity

  const lp = [
    0,                                                   // H0: no causal variant, ABF 1
    Math.log(priors.p1) + S1,
    Math.log(priors.p2) + S2,
    Math.log(priors.p1) + Math.log(priors.p2) + lH3,
    Math.log(priors.p12) + S12,
  ]
  const norm = logSumExp(lp)
  const pp = lp.map(x => Math.exp(x - norm)) as ColocResult['pp']

  const snpPP4 = new Float64Array(n)
  let best = 0
  for (let i = 0; i < n; i++) {
    snpPP4[i] = Math.exp(l12[i]! - S12)
    if (snpPP4[i]! > snpPP4[best]!) best = i
  }
  return { n, pp, snpPP4, best }
}

// ---- assembling a locus from store reads ------------------------------------------------------

/** One trait's per-variant statistics at a site, keyed so two traits can be aligned. `rsNumber`
 *  rides along from the catalog side (0 = no dbSNP record) so the shared variant can be named. */
export interface TraitRow { position: number; ref: string; alt: string; z: number; se: number; rsNumber?: number }

/** What a gene page shows: the posteriors, how well the shared variant is pinned down, and what the
 *  run cost in variants. */
export interface LocusColoc {
  pp: ColocResult['pp']
  /** variants both studies carry, and the QTL block's row count they came from */
  nShared: number
  nQtlRows: number
  /** QTL rows the phenotype did not test (nlp code 65535): a legitimate exclusion, and a run can
   *  hold them even though it is a contiguous vidx range (SPEC section 8). */
  nNotTested: number
  /** QTL rows whose p underflowed to 0 (code 65534), so no z exists. Not benign: an underflowed row
   *  is among the *strongest* signals in the window, and dropping one silently would bias coloc. */
  nUnderflow: number
  /** rows dropped for anything else -- a null SE code, or a results set with no dof at all */
  nNoStats: number
  /** how many variants the shared-variant posterior needs to reach 0.95 — 1 means the variant is
   *  identified, 35 means PP.H4 is a statement about the region and not about any variant */
  credible95: number
  /** the largest `snpPP4` and where it sits. `rsNumber` is 0 where the catalog has no dbSNP record,
   *  and the caller falls back to the position. */
  top: { position: number; ref: string; alt: string; rsNumber: number; snpPP4: number }
  /** `snpPP4` per QTL block row, NaN where the GWAS has no row for that site. NaN is not zero: the
   *  variant was never in the comparison, so a plot must keep it visually out of play rather than
   *  draw it as a rejected candidate. */
  perRow: Float64Array
  /** block row index of the top shared variant */
  topRow: number
  priors: AbfPriors
}

/** The QTL side of a locus: `z` is the stored t statistic, which is what ABF wants, so the slope
 *  (the one reconstructed value in the format) is never formed. NaN where the row has no usable
 *  statistics, which `alignTraits` then drops. */
export function qtlRows(block: { nRows: number; nlp: Float64Array; se: Float64Array; slope: Float64Array; negative: Uint8Array },
  run: { position: Int32Array | Uint32Array; ref: string[]; alt: string[]; rsNumber: Uint32Array }): TraitRow[] {
  const out: TraitRow[] = new Array(block.nRows)
  for (let i = 0; i < block.nRows; i++) {
    const se = block.se[i]!
    // slope / se == sign * t(p, dof); NaN propagates from a null or underflowed row
    const z = (block.negative[i] ? -1 : 1) * Math.abs(block.slope[i]! / se)
    out[i] = { position: run.position[i]!, ref: run.ref[i]!, alt: run.alt[i]!, rsNumber: run.rsNumber[i]!, z, se }
  }
  return out
}

/** The GWAS side: `z = beta / se`, both lossless at the source's printed precision. Ordered by
 *  descending N so that where a site carries two rows (SPEC section 10) the better-powered one wins
 *  the dedupe in `alignTraits`. */
export function gwasRows(g: { rows: number; position: Int32Array; ref: string[]; alt: string[];
  beta: Float64Array; se: Float64Array; n: Int32Array }): TraitRow[] {
  const order = Array.from({ length: g.rows }, (_, i) => i).sort((a, b) => g.n[b]! - g.n[a]!)
  return order.map(i => ({ position: g.position[i]!, ref: g.ref[i]!, alt: g.alt[i]!,
    z: g.beta[i]! / g.se[i]!, se: g.se[i]! }))
}

/** coloc.abf over one phenotype's window against the GWAS window the page already holds. Null when
 *  there is no GWAS, no overlap, or nothing to run on. */
export function colocLocus(
  block: Parameters<typeof qtlRows>[0] | null,
  run: Parameters<typeof qtlRows>[1] | null,
  gwas: Parameters<typeof gwasRows>[0] | null,
  priors: AbfPriors = DEFAULT_PRIORS,
): LocusColoc | null {
  if (!block || !run || !gwas?.rows) return null
  const q = qtlRows(block, run)
  // the three reasons a row has no z, kept apart: NaN -log10 p is "not tested", Infinity is p = 0,
  // and anything else left is a null SE or a results set with no dof
  let nNotTested = 0, nUnderflow = 0, nNoStats = 0
  for (let i = 0; i < block.nRows; i++) {
    if (Number.isFinite(q[i]!.z) && Number.isFinite(q[i]!.se)) continue
    if (Number.isNaN(block.nlp[i]!)) nNotTested++
    else if (!Number.isFinite(block.nlp[i]!)) nUnderflow++
    else nNoStats++
  }
  const { z1, se1, z2, se2, sites, aIndex } = alignTraits(q, gwasRows(gwas))
  if (!sites.length) return null
  const r = colocAbf(z1, se1, z2, se2, priors)
  const order = Array.from(r.snpPP4.keys()).sort((a, b) => r.snpPP4[b]! - r.snpPP4[a]!)
  let cum = 0, credible95 = 0
  for (const i of order) { cum += r.snpPP4[i]!; credible95++; if (cum >= 0.95) break }
  const perRow = new Float64Array(block.nRows).fill(NaN)
  for (let i = 0; i < aIndex.length; i++) perRow[aIndex[i]!] = r.snpPP4[i]!
  const t = sites[r.best]!
  return { pp: r.pp, nShared: sites.length, nQtlRows: block.nRows, nNotTested, nUnderflow, nNoStats, credible95,
    top: { position: t.position, ref: t.ref, alt: t.alt, rsNumber: t.rsNumber ?? 0, snpPP4: r.snpPP4[r.best]! },
    perRow, topRow: aIndex[r.best]!, priors }
}

/**
 * Align two traits on `(position, ref, alt)`.
 *
 * Both sides of a v1 store are oriented to the reference with effects ALT-relative (SPEC section
 * 7), so a site key matches without any strand or effect-allele reasoning and no sign flip is ever
 * needed -- the thing that makes cross-study comparison unsafe elsewhere.
 *
 * The DCM GWAS lists 796,531 indels twice, once per allele order, with different N and statistics
 * that the meta-analysis did not merge (SPEC section 10), so after orientation one site can carry
 * two rows. coloc must see each site once: this keeps the row the caller ordered first, and
 * `gwasRows` below orders by descending N, so the better-powered measurement wins. The locus plot
 * resolves the same duplication by smallest p, so a lead variant shown on the plot and the row used
 * here can differ for those indels.
 */
export function alignTraits(a: TraitRow[], b: TraitRow[]): {
  z1: Float64Array; se1: Float64Array; z2: Float64Array; se2: Float64Array; sites: TraitRow[]
  /** index into `a` of each kept site, so a caller can put results back on its own rows */
  aIndex: Int32Array
} {
  const key = (r: TraitRow) => `${r.position}:${r.ref}:${r.alt}`
  const seen = new Map<string, TraitRow>()
  for (const r of b) if (!seen.has(key(r))) seen.set(key(r), r)
  const z1: number[] = [], se1: number[] = [], z2: number[] = [], se2: number[] = [], sites: TraitRow[] = []
  const idx: number[] = []
  const used = new Set<string>()
  for (let i = 0; i < a.length; i++) {
    const r = a[i]!
    const k = key(r)
    if (used.has(k)) continue
    const m = seen.get(k)
    if (!m) continue
    if (!Number.isFinite(r.z) || !Number.isFinite(r.se) || !Number.isFinite(m.z) || !Number.isFinite(m.se)) continue
    if (r.se <= 0 || m.se <= 0) continue
    used.add(k)
    z1.push(r.z); se1.push(r.se); z2.push(m.z); se2.push(m.se); sites.push(r); idx.push(i)
  }
  return { z1: new Float64Array(z1), se1: new Float64Array(se1), z2: new Float64Array(z2),
    se2: new Float64Array(se2), sites, aIndex: Int32Array.from(idx) }
}
