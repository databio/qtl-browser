/**
 * Genes and phenotypes, per chromosome, in plain JS (SPEC sections 6 and 8). A page reads only the
 * chromosomes it shows: that chromosome's genes and exon models (annotation) and, per experiment it
 * shows, that experiment's search index part, each one whole-object read per session that the
 * browser cache keeps after. A gene named in a URL is found through the gene lookup (two small range
 * reads), so a first gene page never reads a whole-genome table. The gene list and the search box
 * need every chromosome; they read the whole-genome genes table and the cis experiments' whole
 * search indexes instead, the first time they are used, and chromosomes not yet read come from
 * those after.
 *
 * One experiment is one scan and numbers its own rows (SPEC section 8), so a phenotype row here
 * carries the experiment it came from, and an `ord` is only meaningful with it. `SearchHit` joins a
 * gene's annotation row to its cis phenotypes on the same chromosome: its eQTL phenotype (block,
 * run, significance) and its sQTL counts. A gene's phenotypes are read from the parts of its
 * annotation chromosome; a phenotype the experiment places on another chromosome than its gene is
 * not shown on the gene page (TOPCHeF has none).
 *
 * No query engine here: pages that only list or look up genes never start DuckDB.
 */
import type { Table } from 'apache-arrow'
import { chromExonModels, chromGenes, EQTL_TYPE, experimentPart, getStore, lookupGenes, partKey, partOfOrd, SQTL_TYPE,
  unplacedIndex, wholeExperimentIndex, wholeGenes, type QtlExperiment, type Store } from './store'
import type { Exon, SearchHit, WindowGene } from './queries'

/** A search index row (SPEC section 8), with the experiment it belongs to. A cis scan's rows carry
 *  the block and run columns; a trans scan's the frame columns; the others are null. */
export interface PhenotypeRow {
  experiment: string
  ord: number; phenotype_type: string; phenotype_id: string; gene_id: string | null; chr: string | null
  /** null when the scan assessed no significance, which is not false (lib/significance.ts) */
  significant: boolean | null
  /** the value the scan's rule tested, whichever `permuted` column it names */
  sig_value: number | null
  blk_off: number | null; blk_len: number | null; var_start: number | null; n_var: number | null
  var_off: number | null; var_len: number | null; w_lo: number | null; w_hi: number | null
  /** a cis run's strongest row: its decoded -log10 p (Infinity for a p of 0) and position; null
   *  when no row has a p, or the index predates the columns */
  lead_nlp: number | null; lead_pos: number | null
  trans_off: number | null; trans_len: number | null; n_trans: number | null
}

/** An annotation genes row (SPEC section 6). */
export interface GeneRow { gene_id: string; version: number; name: string; biotype: string; chr: string; tss: number; strand: string; start: number; end: number }

const SHARED_COLS = ['ord', 'phenotype_type', 'phenotype_id', 'gene_id', 'chr', 'significant', 'sig_value']
const CIS_COLS = ['blk_off', 'blk_len', 'var_start', 'n_var', 'var_off', 'var_len', 'w_lo', 'w_hi']
const LEAD_COLS = ['lead_nlp', 'lead_pos']
const TRANS_COLS = ['trans_off', 'trans_len', 'n_trans']
const GENE_COLS = ['gene_id', 'version', 'name', 'biotype', 'chr', 'tss', 'strand', 'start', 'end']

/** Plain objects from an Arrow table's columns (nulls as null). `cols` must exist; `optional` come
 *  out null; `ifPresent` are read when the table has them and null when it does not (columns a
 *  store built before them lacks). */
function rowsOf<T>(t: Table | null, cols: string[], optional: string[] = [], extra: Record<string, unknown> = {},
  ifPresent: string[] = []): T[] {
  if (!t) return []
  const vecs = cols.map(c => {
    const v = t.getChild(c)
    if (!v) throw new Error(`table has no column ${c}`)
    return v
  })
  const maybe = ifPresent.map(c => [c, t.getChild(c)] as const)
  const out: T[] = new Array(t.numRows)
  for (let i = 0; i < t.numRows; i++) {
    const o: Record<string, unknown> = { ...extra }
    for (let j = 0; j < cols.length; j++) o[cols[j]] = vecs[j].get(i) ?? null
    for (const [c, v] of maybe) o[c] = v ? v.get(i) ?? null : null
    for (const c of optional) o[c] = null
    out[i] = o as T
  }
  return out
}

/** A QTL experiment's index rows: its own columns, the other scope's as null, tagged with its id. */
function phenotypeRows(e: QtlExperiment, t: Table | null): PhenotypeRow[] {
  return e.kind === 'qtl_cis'
    ? rowsOf<PhenotypeRow>(t, [...SHARED_COLS, ...CIS_COLS], TRANS_COLS, { experiment: e.id }, LEAD_COLS)
    : rowsOf<PhenotypeRow>(t, [...SHARED_COLS, ...TRANS_COLS], [...CIS_COLS, ...LEAD_COLS], { experiment: e.id })
}

function memoBy<T>(make: (key: string) => Promise<T>): (key: string) => Promise<T> {
  const m = new Map<string, Promise<T>>()
  return key => {
    let p = m.get(key)
    if (!p) {
      const q = make(key)
      m.set(key, q)
      q.catch(() => { if (m.get(key) === q) m.delete(key) })
      p = q
    }
    return p
  }
}

const group = <T>(xs: T[], key: (x: T) => string | null) => {
  const m = new Map<string, T[]>()
  for (const x of xs) { const k = key(x); if (k != null) { const a = m.get(k); if (a) a.push(x); else m.set(k, [x]) } }
  return m
}

/** The whole-genome tables by chromosome, once `allHits` has asked for them (null before): the
 *  genes table and each cis experiment's whole index. Rows keep table order, which is the
 *  per-chromosome objects' order too (SPEC sections 6 and 8). */
let whole: Promise<{ genes: Map<string, GeneRow[]>; phen: Map<string, PhenotypeRow[]> }> | null = null
function readWhole() {
  if (!whole) {
    const p = getStore().then(async s => {
      const cis = [...s.cis.values()]
      const [g, ...idx] = await Promise.all([wholeGenes(), ...cis.map(e => wholeExperimentIndex(e.id))])
      const rows = cis.flatMap((e, i) => phenotypeRows(e, idx[i]))
      return { genes: group(rowsOf<GeneRow>(g, GENE_COLS), r => r.chr), phen: group(rows, r => r.chr) }
    })
    whole = p
    p.catch(() => { if (whole === p) whole = null })
  }
  return whole
}

/** One chromosome's genes, by gene id: its genes object, or the whole-genome table once read. */
const genesOf = memoBy(async chr => {
  const genes = whole ? (await whole).genes.get(chr) ?? [] : rowsOf<GeneRow>(await chromGenes(chr), GENE_COLS)
  return { genes, byId: new Map(genes.map(g => [g.gene_id, g])) }
})

/** One chromosome's **cis** search index rows across the study's cis experiments, by gene id: their
 *  parts for that chromosome, or the whole indexes once read. `byOrd` is per experiment. */
const phenotypesOf = memoBy(async chr => {
  const s = await getStore()
  const rows = whole ? (await whole).phen.get(chr) ?? [] : await cisPartRows(s, chr)
  return { rows, byGene: group(rows, r => r.gene_id), byOrd: ordMap(rows) }
})

async function cisPartRows(s: Store, chr: string): Promise<PhenotypeRow[]> {
  const cis = [...s.cis.values()]
  const parts = await Promise.all(cis.map(e => experimentPart(partKey(e.id, chr))))
  return cis.flatMap((e, i) => phenotypeRows(e, parts[i]))
}

const ordMap = (rows: PhenotypeRow[]) => {
  const m = new Map<string, Map<number, PhenotypeRow>>()
  for (const r of rows) { let x = m.get(r.experiment); if (!x) { x = new Map(); m.set(r.experiment, x) } x.set(r.ord, r) }
  return m
}

/** One trans experiment's rows on one chromosome part, or its unplaced residue (`chr` = 'unplaced'). */
const transPartRows = memoBy(async key => {
  const [expId, chr] = key.split('\t')
  const s = await getStore()
  const e = s.experiments.get(expId)
  if (!e || e.kind !== 'qtl_trans') throw new Error(`no trans experiment ${expId}`)
  const t = chr === 'unplaced' ? await unplacedIndex(expId) : await experimentPart(partKey(expId, chr))
  const rows = phenotypeRows(e, t)
  return { rows, byGene: group(rows, r => r.gene_id), byOrd: new Map(rows.map(r => [r.ord, r])) }
})

/** The gene's search hit from its annotation row and its cis phenotypes on that chromosome.
 *  `tested`: its eQTL phenotype has rows; `is_egene`: that phenotype is significant -- null both
 *  without an eQTL phenotype and when the scan assessed no significance, which `tested` and
 *  `significanceOf().assessed` tell apart; `has_results`: any phenotype has a block. */
function hitOf(g: GeneRow, ps: PhenotypeRow[]): SearchHit {
  const e = ps.find(p => p.phenotype_type === EQTL_TYPE && p.blk_off != null) ?? null
  const sq = ps.filter(p => p.phenotype_type === SQTL_TYPE && p.blk_off != null)
  return {
    gene_id: g.gene_id, symbol: g.name, chr: g.chr, tss: g.tss, start: g.start, end: g.end, strand: g.strand, biotype: g.biotype,
    gene_version: g.version,
    ord: e?.ord ?? null, blk_off: e?.blk_off ?? null, blk_len: e?.blk_len ?? null, var_start: e?.var_start ?? null, n_var: e?.n_var ?? null,
    var_off: e?.var_off ?? null, var_len: e?.var_len ?? null, w_lo: e?.w_lo ?? null, w_hi: e?.w_hi ?? null,
    tested: e != null && e.n_var != null, is_egene: e ? e.significant : null, sig_value: e?.sig_value ?? null,
    n_sqtl_sig: sq.filter(p => p.significant).length, n_sqtl: sq.length, has_results: e != null || sq.length > 0,
  }
}

/** Every annotated gene on one catalog chromosome as a search hit, in annotation order (empty off the catalog). */
export const chromHits = memoBy(async (chr): Promise<SearchHit[]> => {
  const s = await getStore()
  if (!s.chroms.has(chr)) return []
  const [g, p] = await Promise.all([genesOf(chr), phenotypesOf(chr)])
  return g.genes.map(x => hitOf(x, p.byGene.get(x.gene_id) ?? []))
})

/** Every annotated gene on every catalog chromosome (the gene list and the search box): the
 *  whole-genome genes table and the cis indexes, read once (chromosomes already read keep theirs). */
export async function allHits(): Promise<SearchHit[]> {
  const s = await getStore()
  await readWhole()
  return (await Promise.all([...s.chroms.keys()].map(chromHits))).flat()
}

/** A study locus for the collection page: a phenotype's strongest cis row, with its gene. */
export interface TopLocus {
  gene_id: string | null; symbol: string | null; phenotype_type: string; phenotype_id: string
  chr: string; lead_pos: number; lead_nlp: number
}

/** A study's strongest loci across its cis scans, by each phenotype's lead -log10 p (the index's
 *  `lead_nlp`; `sig_value` cannot rank them, since a permutation p floors). One per gene, and none
 *  within `spacing` bp of a stronger one on its chromosome, so one region does not fill the list.
 *  Reads what the search box reads: the whole-genome genes table and each cis scan's whole index.
 *  Empty when the indexes predate the lead columns. */
export async function topLoci(n = 20, spacing = 1_000_000): Promise<TopLocus[]> {
  const w = await readWhole()
  const name = new Map<string, string>()
  for (const gs of w.genes.values()) for (const g of gs) name.set(g.gene_id, g.name)
  const rows = [...w.phen.values()].flat().filter(r => r.lead_nlp != null && r.lead_pos != null && r.chr != null)
  // Infinity (a p of 0) sorts first; the comparison avoids Infinity - Infinity
  rows.sort((a, b) => (b.lead_nlp! > a.lead_nlp! ? 1 : b.lead_nlp! < a.lead_nlp! ? -1 : 0))
  const out: TopLocus[] = []
  const genes = new Set<string>()
  for (const r of rows) {
    if (out.length >= n) break
    if (r.gene_id && genes.has(r.gene_id)) continue
    if (out.some(x => x.chr === r.chr && Math.abs(x.lead_pos - r.lead_pos!) < spacing)) continue
    if (r.gene_id) genes.add(r.gene_id)
    out.push({ gene_id: r.gene_id, symbol: r.gene_id ? name.get(r.gene_id) ?? null : null, phenotype_type: r.phenotype_type,
      phenotype_id: r.phenotype_id, chr: r.chr!, lead_pos: r.lead_pos!, lead_nlp: r.lead_nlp! })
  }
  return out
}

/** The gene a URL names by Ensembl id or symbol (ASCII case ignored), preferring a gene with results,
 *  then an eQTL-tested one; null when none is on a catalog chromosome. */
export async function resolveGene(id: string): Promise<SearchHit | null> {
  const s = await getStore()
  const found = (await lookupGenes(id)).filter(r => s.chroms.has(r.chr))
  const hits = (await Promise.all(found.map(async r => (await chromHits(r.chr)).find(h => h.gene_id === r.gene_id) ?? null)))
    .filter((h): h is SearchHit => h != null)
  hits.sort((a, b) => Number(b.has_results) - Number(a.has_results) || Number(b.tested) - Number(a.tested))
  return hits[0] ?? null
}

/** Genes whose symbol or id starts with `q` (case ignored), best first: with results, eQTL-tested,
 *  eGenes, then shorter and alphabetical symbols. */
export async function searchGenes(q: string, limit = 12): Promise<SearchHit[]> {
  const Q = q.toUpperCase()
  const rank = (x: boolean | null) => (x === true ? 0 : x === false ? 1 : 2)
  return (await allHits())
    .filter(h => (h.symbol != null && h.symbol.toUpperCase().startsWith(Q)) || h.gene_id.toUpperCase().startsWith(Q))
    .sort((a, b) => rank(a.has_results) - rank(b.has_results) || rank(a.tested) - rank(b.tested) || rank(a.is_egene) - rank(b.is_egene)
      || (a.symbol ?? '').length - (b.symbol ?? '').length || ((a.symbol ?? '') < (b.symbol ?? '') ? -1 : (a.symbol ?? '') > (b.symbol ?? '') ? 1 : 0))
    .slice(0, limit)
}

/** Genes with a TSS in [start, end] on `chr`, by TSS. */
export async function genesInRegion(chr: string, start: number, end: number): Promise<SearchHit[]> {
  return (await chromHits(chr)).filter(h => h.tss >= start && h.tss <= end).sort((a, b) => a.tss - b.tss)
}

/** Genes overlapping [lo, hi] on `chr`, by start: the gene track's lanes. */
export async function genesInWindow(chr: string, lo: number, hi: number): Promise<WindowGene[]> {
  return (await chromHits(chr)).filter(h => h.end >= lo && h.start <= hi).sort((a, b) => a.start - b.start)
    .map(h => ({ gene_id: h.gene_id, symbol: h.symbol, start: h.start, end: h.end, strand: h.strand, tss: h.tss, biotype: h.biotype }))
}

/** The gene's collapsed exon model (the union of its transcripts' exons, merged by the builder). */
export async function exonModel(geneId: string, chr: string): Promise<Exon[]> {
  const t = await chromExonModels(chr)
  if (!t) return []
  const ids = t.getChild('gene_id')!
  for (let i = 0; i < t.numRows; i++) {
    if (ids.get(i) !== geneId) continue
    const s = t.getChild('exon_starts')!.get(i), e = t.getChild('exon_ends')!.get(i)
    const out: Exon[] = []
    for (let k = 0; k < (s?.length ?? 0); k++) out.push({ start: Number(s.get(k)), end: Number(e.get(k)) })
    return out
  }
  return []
}

/** The gene's cis phenotypes on `chr` (every cis experiment), by (experiment, ord). */
export async function phenotypesOfGene(geneId: string, chr: string): Promise<PhenotypeRow[]> {
  return [...((await phenotypesOf(chr)).byGene.get(geneId) ?? [])]
    .sort((a, b) => (a.experiment < b.experiment ? -1 : a.experiment > b.experiment ? 1 : a.ord - b.ord))
}

/** The gene's phenotypes with trans rows, from each trans experiment: its part for the gene's
 *  chromosome (where the builder files a phenotype by its gene's chromosome), then the unplaced
 *  residue for any phenotype of the gene the annotation could not place. */
export async function transPhenotypesOfGene(geneId: string, chr: string): Promise<PhenotypeRow[]> {
  const s = await getStore()
  const out: PhenotypeRow[] = []
  await Promise.all([...s.trans.values()].map(async e => {
    const [part, unplaced] = await Promise.all([transPartRows(partKey(e.id, chr)),
      e.search_index_unplaced ? transPartRows(partKey(e.id, 'unplaced')) : null])
    out.push(...(part.byGene.get(geneId) ?? []), ...(unplaced?.byGene.get(geneId) ?? []))
  }))
  return out.filter(r => r.trans_off != null)
}

/** Phenotypes by ord within one experiment (its hits files name them so), read from the parts that
 *  hold them. */
export async function phenotypesByOrd(e: QtlExperiment, ords: number[]): Promise<Map<number, PhenotypeRow>> {
  const parts = new Map<string, number[]>()
  for (const o of new Set(ords)) {
    const p = partOfOrd(e, o)
    if (p === undefined) throw new Error(`experiment ${e.id}: the search index has no phenotype with ord ${o}`)
    const a = parts.get(p)
    if (a) a.push(o); else parts.set(p, [o])
  }
  const out = new Map<number, PhenotypeRow>()
  await Promise.all([...parts].map(async ([p, os]) => {
    const byOrd = e.kind === 'qtl_cis'
      ? (await phenotypesOf(p)).byOrd.get(e.id) ?? new Map<number, PhenotypeRow>()
      : (await transPartRows(partKey(e.id, p))).byOrd
    for (const o of os) { const r = byOrd.get(o); if (r) out.set(o, r) }
  }))
  return out
}

/** Phenotypes of one type on `chr` whose run covers `vidx`, by block offset (the cis scan). */
export async function phenotypesCovering(phenotypeType: string, chr: string, vidx: number): Promise<PhenotypeRow[]> {
  return (await phenotypesOf(chr)).rows
    .filter(p => p.phenotype_type === phenotypeType && p.var_start != null && p.n_var != null && p.var_start <= vidx && vidx < p.var_start + p.n_var)
    .sort((a, b) => a.blk_off! - b.blk_off!)
}

/** A gene's annotation name, chromosome and TSS: from its chromosome's genes when `chr` is known
 *  (one read serves every gene there), else through the gene lookup; null when unannotated. */
export async function geneInfo(geneId: string, chr: string | null): Promise<{ name: string; chr: string; tss: number } | null> {
  if (chr) {
    const g = (await genesOf(chr)).byId.get(geneId)
    if (g) return g
  }
  return (await lookupGenes(geneId)).find(r => r.gene_id === geneId) ?? null
}
