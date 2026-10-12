/**
 * The trans tables, from the qtlstore (SPEC.md sections 8 and 9), as in-memory DuckDB tables with
 * the `TransRow` columns, so `transRows`, `transCount` and `transAll` page off them.
 *
 * - Gene page: the frames of the gene's phenotypes (its eQTL phenotype and every intron, the
 *   trans-only ones included) in each trans scan's object. The builder writes a gene's frames back
 *   to back (SPEC section 9), so they are one request per object: the smallest `trans_off` to the
 *   largest `trans_off + trans_len`. The SE is stored in the frame; r² is formed from beta, SE and
 *   the scan's dof.
 * - Variant page: each trans scan's hits records for the variant, one per trans association. A
 *   hits record carries -log10 p and beta but no SE (SPEC section 8), so the SE and r² here are
 *   rebuilt through the inverse t with the scan's dof -- the one derived statistic left in the
 *   reader -- and can differ in the last digits from the gene page's stored SE for the same row.
 */
import { Int64, makeVector, Table, tableToIPC, Utf8, vectorFromArray } from 'apache-arrow'
import { getDB, insertArrow, tableName } from './db'
import { geneInfo, phenotypesByOrd, transPhenotypesOfGene } from './gene-index'
import { EQTL_TYPE, fetchDecoded, getStore, SQTL_TYPE, transFile, type Store } from './store'
import { decodeTransFrame, HIT_TRANS, transR2, transSeFromHit, type VariantRecord } from './store-decode'
import type { Hits } from './variant'
import type { SearchHit } from './queries'

/** The trans table's columns and types, for a table with no rows. */
const TRANS_DDL = 'qtl_type VARCHAR, phenotype_id VARCHAR, gene_id VARCHAR, symbol VARCHAR, gene_chr VARCHAR, gene_tss BIGINT, ' +
  'variant_chr VARCHAR, position INTEGER, rsid VARCHAR, af FLOAT, pval DOUBLE, beta FLOAT, beta_se FLOAT, r2 FLOAT'

interface TransRows {
  qtl: string[]; phenotype: string[]; gene: (string | null)[]; symbol: (string | null)[]; geneChr: (string | null)[]; tss: (bigint | null)[]
  variantChr: string[]; position: number[]; rsid: (string | null)[]; af: number[]; pval: number[]; beta: number[]; se: number[]; r2: number[]
}
const emptyRows = (): TransRows => ({ qtl: [], phenotype: [], gene: [], symbol: [], geneChr: [], tss: [], variantChr: [], position: [],
  rsid: [], af: [], pval: [], beta: [], se: [], r2: [] })

async function insertRows(r: TransRows, detail: Record<string, unknown>): Promise<string> {
  const name = tableName('trans')
  const n = r.qtl.length
  if (!n) { await insertArrow(name, TRANS_DDL, detail); return name }
  for (const p of r.position) if (p > 0x7fffffff) throw new Error(`trans row position ${p} does not fit INTEGER`)
  const utf = (a: (string | null)[]) => vectorFromArray(a, new Utf8())
  const f32 = (a: number[]) => makeVector(Float32Array.from(a))
  const ipc = tableToIPC(new Table({
    qtl_type: utf(r.qtl), phenotype_id: utf(r.phenotype), gene_id: utf(r.gene), symbol: utf(r.symbol), gene_chr: utf(r.geneChr),
    gene_tss: vectorFromArray(r.tss, new Int64()), variant_chr: utf(r.variantChr), position: makeVector(Int32Array.from(r.position)),
    rsid: utf(r.rsid), af: f32(r.af), pval: makeVector(Float64Array.from(r.pval)), beta: f32(r.beta), beta_se: f32(r.se), r2: f32(r.r2),
  }), 'stream')
  await insertArrow(name, ipc, detail)
  // NaN stands for "none" in the typed columns; the trans queries expect SQL NULLs
  const { con } = await getDB()
  await con.query(`UPDATE ${name} SET af = CASE WHEN isnan(af) THEN NULL ELSE af END,
    beta_se = CASE WHEN isnan(beta_se) THEN NULL ELSE beta_se END, r2 = CASE WHEN isnan(r2) THEN NULL ELSE r2 END`)
  return name
}

const qtlOf = (phenotypeType: string) => (phenotypeType === SQTL_TYPE ? 's' : 'e')

interface FrameRow { phenotype_type: string; phenotype_id: string; trans_off: number; trans_len: number; n_trans: number }

/** Every trans row of one gene's phenotypes; the caller drops the table. */
export async function geneTransTable(hit: SearchHit): Promise<string> {
  const s = await getStore()
  const frames = (await transPhenotypesOfGene(hit.gene_id, hit.chr) as FrameRow[])
    .sort((a, b) => (a.phenotype_type < b.phenotype_type ? -1 : a.phenotype_type > b.phenotype_type ? 1 : a.trans_off - b.trans_off))
  const out = emptyRows()
  const detail = { gene_id: hit.gene_id }
  // one request per trans scan, both in flight at once; rows go out eQTL first
  const perType = await Promise.all([EQTL_TYPE, SQTL_TYPE].map(async type => {
    const mine = frames.filter(f => f.phenotype_type === type)
    const tf = transFile(s, type)
    if (!mine.length) return null
    if (!tf) throw new Error(`phenotypes of ${hit.gene_id} have ${type} trans frames but the study has no ${type} trans scan`)
    // the gene's frames are contiguous (sorted by trans_off above), so one range covers them all
    for (let k = 1; k < mine.length; k++)
      if (mine[k].trans_off !== mine[k - 1].trans_off + mine[k - 1].trans_len)
        throw new Error(`${tf.name}: the trans frames of ${hit.gene_id} are not contiguous (${mine[k].phenotype_id})`)
    const off = mine[0].trans_off, len = mine[mine.length - 1].trans_off + mine[mine.length - 1].trans_len - off
    const decoded = await fetchDecoded('store:trans', detail, tf.name, off, len, (bytes, what) =>
      mine.map(f => decodeTransFrame(bytes.subarray(f.trans_off - off, f.trans_off - off + f.trans_len), `${what} ${f.phenotype_id}`)),
      'store:trans-decode')
    return { type, tf, mine, decoded }
  }))
  for (const t of perType) {
    if (!t) continue
    const { type, tf, mine, decoded } = t
    mine.forEach((f, k) => {
      const fr = decoded[k]
      if (fr.n !== f.n_trans) throw new Error(`${tf.name}: ${f.phenotype_id} frame has ${fr.n} rows, the search index says ${f.n_trans}`)
      for (let i = 0; i < fr.n; i++) {
        out.qtl.push(qtlOf(type)); out.phenotype.push(f.phenotype_id); out.gene.push(hit.gene_id); out.symbol.push(hit.symbol)
        out.geneChr.push(hit.chr); out.tss.push(BigInt(hit.tss))
        out.variantChr.push(chromName(s, fr.ordinal[i])); out.position.push(fr.position[i])
        out.rsid.push(fr.rsNumber[i] ? `rs${fr.rsNumber[i]}` : null); out.af.push(fr.af[i])
        out.pval.push(fr.pval[i]); out.beta.push(fr.beta[i]); out.se.push(fr.se[i]); out.r2.push(transR2(fr.beta[i], fr.se[i], tf.dof))
      }
    })
  }
  return insertRows(out, { ...detail, part: 'trans' })
}

function chromName(s: Store, ordinal: number): string {
  const c = s.catalog.chromosomes[ordinal - 1]
  if (!c) throw new Error(`trans row names chromosome ordinal ${ordinal}, not in catalog ${s.catalog.id}`)
  return c.name
}

interface HitGene { ord: number; phenotype_type: string; phenotype_id: string; gene_id: string | null; symbol: string | null; chr: string | null; tss: number | null }

/** One variant's trans rows, from each trans scan's hits records (kind 2); the caller drops the table. */
export async function variantTransTable(v: VariantRecord, hits: Hits): Promise<string> {
  const out = emptyRows()
  for (const set of hits.sets) {
    if (set.experiment.kind !== 'qtl_trans' || !set.trans.length) continue
    const f = set.frame, e = set.experiment
    const ords = [...new Set(set.trans.map(r => f.ord[r]))]
    // each phenotype from its part of this scan's index (the unplaced part for those with no
    // chromosome), and its gene's symbol, chromosome and TSS from the annotation
    const genes = new Map<number, HitGene>()
    await Promise.all([...(await phenotypesByOrd(e, ords)).values()].map(async p => {
      const g = p.gene_id ? await geneInfo(p.gene_id, p.chr) : null
      genes.set(p.ord, { ord: p.ord, phenotype_type: p.phenotype_type, phenotype_id: p.phenotype_id, gene_id: p.gene_id,
        symbol: g?.name ?? null, chr: g?.chr ?? null, tss: g?.tss ?? null })
    }))
    const tCache = new Map<number, number>()
    for (const r of set.trans) {
      if (f.kind[r] !== HIT_TRANS) continue
      const g = genes.get(f.ord[r])
      if (!g) throw new Error(`experiment ${e.id}: the search index has no phenotype with ord ${f.ord[r]}`)
      const nlp = f.value[r], beta = f.beta[r]
      const { se, r2 } = transSeFromHit(nlp, beta, e.dof, tCache)
      out.qtl.push(qtlOf(e.phenotype_type)); out.phenotype.push(g.phenotype_id); out.gene.push(g.gene_id ?? g.phenotype_id)
      out.symbol.push(g.symbol); out.geneChr.push(g.chr); out.tss.push(g.tss == null ? null : BigInt(g.tss))
      out.variantChr.push(v.chr); out.position.push(v.position); out.rsid.push(v.rsNumber ? `rs${v.rsNumber}` : null); out.af.push(v.af)
      out.pval.push(nlp === Infinity ? 0 : Math.pow(10, -nlp)); out.beta.push(beta); out.se.push(se); out.r2.push(r2)
    }
  }
  return insertRows(out, { variant: v.vidx, part: 'trans' })
}
