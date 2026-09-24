/**
 * Genes and phenotypes, per chromosome, in plain JS (SPEC sections 6 and 8). A page reads only the
 * chromosomes it shows: that chromosome's genes and exon models (annotation) and its search index
 * part (experiment), each one whole-object read per session that the browser cache keeps after. A
 * gene named in a URL is found through the gene lookup (two small range reads), so a first gene page
 * never reads a whole-genome table. The gene list and the search box need every chromosome; they
 * read the whole-genome genes table and search index instead (two requests rather than one per
 * chromosome), the first time they are used, and chromosomes not yet read come from those after.
 *
 * `SearchHit` joins a gene's annotation row to its phenotypes on the same chromosome: its eQTL
 * phenotype (block, run, significance) and its sQTL counts. A gene's phenotypes are read from the
 * part of its annotation chromosome; a phenotype the experiment places on another chromosome than
 * its gene is not shown on the gene page (TOPCHeF and GTEx have none).
 *
 * No query engine here: pages that only list or look up genes never start DuckDB.
 */
import type { Table } from 'apache-arrow'
import { chromExonModels, chromGenes, chromIndex, EQTL_TYPE, getStore, lookupGenes, partOfOrd, SQTL_TYPE, transOnlyIndex, wholeGenes,
  wholeIndex } from './store'
import type { Exon, SearchHit, WindowGene } from './queries'

/** A search index row (SPEC section 8, `results.INDEX_SCHEMA`). */
export interface PhenotypeRow {
  ord: number; phenotype_type: string; phenotype_id: string; gene_id: string | null; chr: string | null
  significant: boolean; p_perm: number | null
  blk_off: number | null; blk_len: number | null; var_start: number | null; n_var: number | null
  var_off: number | null; var_len: number | null; w_lo: number | null; w_hi: number | null
  trans_off: number | null; trans_len: number | null; n_trans: number | null
}

/** An annotation genes row (SPEC section 6). */
export interface GeneRow { gene_id: string; version: number; name: string; biotype: string; chr: string; tss: number; strand: string; start: number; end: number }

const PHENOTYPE_COLS = ['ord', 'phenotype_type', 'phenotype_id', 'gene_id', 'chr', 'significant', 'p_perm', 'blk_off', 'blk_len',
  'var_start', 'n_var', 'var_off', 'var_len', 'w_lo', 'w_hi', 'trans_off', 'trans_len', 'n_trans']
const GENE_COLS = ['gene_id', 'version', 'name', 'biotype', 'chr', 'tss', 'strand', 'start', 'end']

/** Plain objects from an Arrow table's columns (nulls as null). */
function rowsOf<T>(t: Table | null, cols: string[]): T[] {
  if (!t) return []
  const vecs = cols.map(c => {
    const v = t.getChild(c)
    if (!v) throw new Error(`table has no column ${c}`)
    return v
  })
  const out: T[] = new Array(t.numRows)
  for (let i = 0; i < t.numRows; i++) {
    const o: Record<string, unknown> = {}
    for (let j = 0; j < cols.length; j++) o[cols[j]] = vecs[j].get(i) ?? null
    out[i] = o as T
  }
  return out
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

/** The whole-genome tables by chromosome, once `allHits` has asked for them (null before). Rows keep
 *  table order, which is the per-chromosome objects' order too (SPEC sections 6 and 8). */
let whole: Promise<{ genes: Map<string, GeneRow[]>; phen: Map<string, PhenotypeRow[]> }> | null = null
function readWhole() {
  if (!whole) {
    const p = Promise.all([wholeGenes(), wholeIndex()]).then(([g, i]) => ({
      genes: group(rowsOf<GeneRow>(g, GENE_COLS), r => r.chr), phen: group(rowsOf<PhenotypeRow>(i, PHENOTYPE_COLS), r => r.chr) }))
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

/** One chromosome's search index rows, by gene id: its part, or the whole index once read. */
const phenotypesOf = memoBy(async chr => {
  const rows = whole ? (await whole).phen.get(chr) ?? [] : rowsOf<PhenotypeRow>(await chromIndex(chr), PHENOTYPE_COLS)
  return { rows, byGene: group(rows, r => r.gene_id), byOrd: new Map(rows.map(r => [r.ord, r])) }
})

const transOnly = (() => {
  let p: Promise<{ rows: PhenotypeRow[]; byOrd: Map<number, PhenotypeRow> }> | null = null
  return () => {
    if (!p) {
      const q = transOnlyIndex().then(t => { const rows = rowsOf<PhenotypeRow>(t, PHENOTYPE_COLS); return { rows, byOrd: new Map(rows.map(r => [r.ord, r])) } })
      p = q
      q.catch(() => { if (p === q) p = null })
    }
    return p
  }
})()

/** The gene's search hit from its annotation row and its phenotypes on that chromosome: the rules of
 *  v1's former `search_index` table. `tested`: its eQTL phenotype has rows; `is_egene`: that
 *  phenotype is significant (null without one); `has_results`: any phenotype has a block. */
function hitOf(g: GeneRow, ps: PhenotypeRow[]): SearchHit {
  const e = ps.find(p => p.phenotype_type === EQTL_TYPE && p.blk_off != null) ?? null
  const sq = ps.filter(p => p.phenotype_type === SQTL_TYPE && p.blk_off != null)
  return {
    gene_id: g.gene_id, symbol: g.name, chr: g.chr, tss: g.tss, start: g.start, end: g.end, strand: g.strand, biotype: g.biotype,
    gene_version: g.version,
    ord: e?.ord ?? null, blk_off: e?.blk_off ?? null, blk_len: e?.blk_len ?? null, var_start: e?.var_start ?? null, n_var: e?.n_var ?? null,
    var_off: e?.var_off ?? null, var_len: e?.var_len ?? null, w_lo: e?.w_lo ?? null, w_hi: e?.w_hi ?? null,
    tested: e != null && e.n_var != null, is_egene: e ? e.significant : null,
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
 *  whole-genome genes table and search index, read once (chromosomes already read keep theirs). */
export async function allHits(): Promise<SearchHit[]> {
  const s = await getStore()
  await readWhole()
  return (await Promise.all([...s.chroms.keys()].map(chromHits))).flat()
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

/** The gene's phenotypes on `chr`, by ord. */
export async function phenotypesOfGene(geneId: string, chr: string): Promise<PhenotypeRow[]> {
  return [...((await phenotypesOf(chr)).byGene.get(geneId) ?? [])].sort((a, b) => a.ord - b.ord)
}

/** The gene's phenotypes with trans rows: those on `chr` and its trans-only phenotypes. */
export async function transPhenotypesOfGene(geneId: string, chr: string): Promise<PhenotypeRow[]> {
  const [p, t] = await Promise.all([phenotypesOf(chr), transOnly()])
  return [...(p.byGene.get(geneId) ?? []), ...t.rows.filter(r => r.gene_id === geneId)].filter(r => r.trans_off != null)
}

/** Phenotypes by ord (the hits files name them so), read from the parts that hold them. */
export async function phenotypesByOrd(ords: number[]): Promise<Map<number, PhenotypeRow>> {
  const s = await getStore()
  const parts = new Map<string | null, number[]>()
  for (const o of new Set(ords)) {
    const p = partOfOrd(s, o)
    if (p === undefined) throw new Error(`the search index has no phenotype with ord ${o}`)
    const a = parts.get(p)
    if (a) a.push(o); else parts.set(p, [o])
  }
  const out = new Map<number, PhenotypeRow>()
  await Promise.all([...parts].map(async ([p, os]) => {
    const byOrd = p === null ? (await transOnly()).byOrd : (await phenotypesOf(p)).byOrd
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
