import { phenotypeLabel, type CollectionDoc, type ExperimentDoc, type Facets, type GwasExperiment, type IndexCollection, type StoreIndex } from './store'
import { fmtInt } from './format'

/** What the page's collection is, once the index has said: `qtl` opens a study, `gwas` and
 *  `parent` have pages of their own, `missing` is a digest this store does not hold. */
export type CollectionKind = 'qtl' | 'gwas' | 'parent' | 'missing'

export function collectionKind(c: IndexCollection | undefined): CollectionKind {
  if (!c) return 'missing'
  if (c.experiments.some(e => e.kind !== 'gwas')) return 'qtl'
  return c.experiments.length ? 'gwas' : 'parent'
}

// a consortium is its own name; a person is their surname
const authorName = (c: CollectionDoc): string | null => {
  const a = c.citation?.first_author
  if (typeof a !== 'string') return null
  return /consortium/i.test(a) ? a.replace(/^the\s+/i, '') : a.split(' ').pop()!
}

/** What lists, crumbs and pickers call a collection: its title, which the store keeps short
 *  ("TOPCHeF LV"; SPEC section 8). */
export const shortName = (c: CollectionDoc): string => c.title

/** What lists and pickers call a GWAS: its collection's short name, with the experiment id when
 *  the collection holds more than one GWAS. */
export function gwasName(g: GwasExperiment, ix: StoreIndex | null): string {
  const c = ix?.collections.get(g.collection)
  if (!c) return g.title || g.id
  const n = c.experiments.filter(e => e.kind === 'gwas').length
  return n > 1 ? `${shortName(c.doc)} · ${g.id}` : shortName(c.doc)
}

/** "Murray et al. 2026, medRxiv": the collection's citation, short. */
export function citationLine(c: CollectionDoc): string {
  const ct = c.citation ?? {}
  const a = authorName(c)
  const who = a ? (/consortium/i.test(a) ? a : `${a} et al.`) : null
  return [who && ct.year ? `${who} ${ct.year}` : who ?? (ct.year ? String(ct.year) : null), typeof ct.venue === 'string' ? ct.venue : null]
    .filter(Boolean).join(', ')
}

export const doiUrl = (c: CollectionDoc): string | null =>
  typeof c.citation?.doi === 'string' ? `https://doi.org/${c.citation.doi}` : null

/** "cis eQTL", "trans sQTL", "GWAS". */
export function experimentLabel(e: ExperimentDoc): string {
  if (e.kind === 'gwas') return 'GWAS'
  return `${e.kind === 'qtl_cis' ? 'cis' : 'trans'} ${phenotypeLabel(e.phenotype_type)}`
}

/** What an experiment holds, in a few words: phenotypes (and significant ones) or GWAS variants. */
export function experimentSize(e: ExperimentDoc): string {
  if (e.kind === 'gwas') return `${fmtInt(e.n_rows)} variants`
  if (e.kind === 'qtl_trans') return `${fmtInt(e.n_phenotypes)} phenotypes · ${fmtInt(e.n_rows)} pairs`
  const sig = e.counts?.significant
  return `${fmtInt(e.n_phenotypes)} phenotypes${sig != null ? ` · ${fmtInt(sig)} significant` : ''}`
}

/** "eQTL, sQTL · cis, trans": a study's experiments in one line. */
export function studySummary(c: IndexCollection): string {
  const qtl = c.experiments.filter(e => e.kind !== 'gwas')
  const types = [...new Set(qtl.map(e => phenotypeLabel(e.phenotype_type)))]
  const scans = [...new Set(qtl.map(e => (e.kind === 'qtl_cis' ? 'cis' : 'trans')))]
  return [types.join(', '), scans.join(', ')].filter(Boolean).join(' · ')
}

/** "12,504,079 variants · 5,022 cases": a GWAS in one line. */
export function gwasSummary(g: GwasExperiment): string {
  const cases = g.source?.n_cases
  return `${fmtInt(g.n_rows)} variants${typeof cases === 'number' ? ` · ${fmtInt(cases)} cases` : ''}`
}

/** The facets some experiments state (a collection's, or one), each value once, in this order: a
 *  study's scans usually agree, and where they do not every value is listed. */
const FACET_ORDER: (keyof Facets)[] = ['trait', 'organ', 'tissue', 'cell_type', 'assay', 'ancestry', 'n_samples', 'sex', 'condition']

export function facetValues(exps: ExperimentDoc[]): { key: keyof Facets; values: string[] }[] {
  return FACET_ORDER.map(key => {
    const vs = [...new Set(exps.map(e => e.facets?.[key]).filter(v => v != null)
      .map(v => (key === 'n_samples' ? fmtInt(v as number) : String(v))))]
    return { key, values: vs }
  }).filter(f => f.values.length > 0)
}

/** "heart · left ventricle · bulk RNA-seq · multi-ancestry": the facets on one line. `sex: both`
 *  says nothing a reader would miss, so it is left out here. */
export const facetLine = (exps: ExperimentDoc[]): string =>
  facetValues(exps).filter(f => !(f.key === 'sex' && f.values.join() === 'both'))
    .flatMap(f => (f.key === 'n_samples' ? f.values.map(v => `n = ${v}`) : f.values)).filter((v, i, a) => a.indexOf(v) === i).join(' · ')

const FACET_LABEL: Record<keyof Facets, string> = {
  trait: 'Trait', organ: 'Organ', tissue: 'Tissue', cell_type: 'Cell type', assay: 'Assay', ancestry: 'Ancestry',
  n_samples: 'Samples', sex: 'Sex', condition: 'Condition',
}

/** The facets as key/value rows, one per facet stated, its values joined. */
export const facetRows = (exps: ExperimentDoc[]): { label: string; value: string }[] =>
  facetValues(exps).map(f => ({ label: FACET_LABEL[f.key], value: f.values.join(', ') }))
