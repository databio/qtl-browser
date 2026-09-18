/**
 * The gene page's trans tables from the trans pack (SPEC.md section 12): one Range request for the
 * gene's frame at `trans_off`/`trans_len` in `manifest.packs.files.trans[chr]`, decoded on the main
 * thread, and inserted into DuckDB with the columns and types of the trans parquet (`TransRow`), so
 * `transRows`, `transCount`, and `transAll` page off it unchanged. A gene without trans rows gets the
 * same table, empty, and sends no request.
 */
import { makeVector, Table, tableToIPC, Utf8, vectorFromArray } from 'apache-arrow'
import { dropTable, getDB, lit, tableName } from './db'
import { fetchDecoded, insertArrow, packFile, packManifest } from './pack'
import { decodeTransFrame, deriveTrans, PackError, TRANS_CHROMS, transPhenotypeIds, type TransFrame } from './pack-decode'
import type { Hits } from './variant-pack'
import type { VariantRecord } from './pack-decode'
import type { SearchHit } from './queries'

/** The trans parquet's columns and types (DESCRIBE), for a gene without trans rows. */
const TRANS_DDL = 'qtl_type VARCHAR, phenotype_id VARCHAR, gene_id VARCHAR, symbol VARCHAR, gene_chr VARCHAR, gene_tss BIGINT, ' +
  'variant_chr VARCHAR, position INTEGER, rsid VARCHAR, af FLOAT, pval DOUBLE, beta FLOAT, beta_se FLOAT, r2 FLOAT'

/** A decoded frame as an Arrow IPC stream in the TRANS_DDL column order; the gene columns come from `hit`. */
function transIPC(hit: SearchHit, f: TransFrame): Uint8Array {
  const n = f.nE + f.nS
  const position = new Int32Array(n)
  for (let i = 0; i < n; i++) {
    if (f.position[i] > 0x7fffffff) throw new PackError(`${hit.gene_id} trans row ${i}: position ${f.position[i]} does not fit INTEGER`)
    position[i] = f.position[i]
  }
  const same = (v: string | null) => vectorFromArray(new Array<string | null>(n).fill(v), new Utf8())
  return tableToIPC(new Table({
    qtl_type: vectorFromArray(new Array<string>(n).fill('e', 0, f.nE).fill('s', f.nE), new Utf8()),
    phenotype_id: vectorFromArray(transPhenotypeIds(f, hit), new Utf8()),
    gene_id: same(hit.gene_id),
    symbol: same(hit.symbol),
    gene_chr: same(hit.chr),
    gene_tss: makeVector(new BigInt64Array(n).fill(BigInt(hit.tss))),
    variant_chr: vectorFromArray(Array.from(f.variantChr, c => TRANS_CHROMS[c]), new Utf8()),
    position: makeVector(position),
    rsid: vectorFromArray(Array.from(f.rsNumber, x => (x ? `rs${x}` : null)), new Utf8()),
    af: makeVector(Float32Array.from(f.af)),
    pval: makeVector(f.pval),
    beta: makeVector(Float32Array.from(f.beta)),
    beta_se: makeVector(Float32Array.from(f.betaSe)),
    r2: makeVector(Float32Array.from(f.r2)),
  }), 'stream')
}

/** The gene's trans rows as a new in-memory table (request timed as pack:trans, decode and Arrow
 *  build as pack:trans-decode); the caller drops it. */
export async function geneTransTable(hit: SearchHit): Promise<string> {
  let ipc: Uint8Array | null = null
  if (hit.trans_len != null) {
    if (hit.trans_off == null) throw new Error(`${hit.gene_id}: search_index has trans_len but no trans_off`)
    const off = hit.trans_off, len = hit.trans_len
    const m = await packManifest()
    ipc = await fetchDecoded('pack:trans', { gene_id: hit.gene_id }, packFile(m.files.trans, hit.chr, 'trans'), off, len,
      (bytes, what) => transIPC(hit, decodeTransFrame(bytes, { e: m.dof.eqtl, s: m.dof.sqtl }, what)), 'pack:trans-decode')
  }
  const name = tableName('trans')
  await insertArrow(name, ipc ?? TRANS_DDL, { gene_id: hit.gene_id, part: 'trans' })
  return name
}


// ---- the variant page's trans rows, from its hits frame (SPEC.md section 13) -----------------

/** The kind 0-1 rows of one variant, as an Arrow IPC stream: the row's gene `ord` and its intron
 *  fields, plus the values derived from p and beta with the frame's scales, exactly as section 12
 *  derives them for the gene-keyed pack. */
function hitsTransIPC(hits: Hits, dof: { e: number; s: number }): Uint8Array {
  const { frame: f, trans } = hits
  const n = trans.length
  const ord = new Uint16Array(n), qtl = new Array<string>(n), strand = new Array<string>(n)
  const intronStart = new Uint32Array(n), intronEnd = new Uint32Array(n), cluster = new Uint32Array(n)
  const pval = new Float64Array(n), beta = new Float32Array(n), betaSe = new Float32Array(n), r2 = new Float32Array(n)
  const tByCode = [new Map<number, number>(), new Map<number, number>()]
  for (let i = 0; i < n; i++) {
    const r = trans[i], k = f.kind[r], d = k === 0 ? dof.e : dof.s
    ord[i] = f.ord[r]
    qtl[i] = k === 0 ? 'e' : 's'
    strand[i] = (f.flags[r] & 1) ? '-' : '+'
    intronStart[i] = f.intronStart[r]; intronEnd[i] = f.intronEnd[r]; cluster[i] = f.cluster[r]
    const v = deriveTrans(f.v1[r], f.v3[r], f.trans.nlpMax, f.trans.betaMax, d, tByCode[k])
    pval[i] = v.pval; beta[i] = v.beta; betaSe[i] = v.betaSe; r2[i] = v.r2
  }
  return tableToIPC(new Table({
    ord: makeVector(ord),
    qtl_type: vectorFromArray(qtl, new Utf8()),
    strand: vectorFromArray(strand, new Utf8()),
    intron_start: makeVector(intronStart),
    intron_end: makeVector(intronEnd),
    cluster: makeVector(cluster),
    pval: makeVector(pval),
    beta: makeVector(beta),
    beta_se: makeVector(betaSe),
    r2: makeVector(r2),
  }), 'stream')
}

/** One variant's trans rows as a new in-memory table with the trans parquet's columns and types
 *  (`TransRow`), so `TransTable keyedBy="variant"` and the trans queries stay as they are. The
 *  gene columns come from `search_index` by `ord`; the variant columns come from `v`. The caller
 *  drops the table. */
export async function variantTransTable(v: VariantRecord, hits: Hits): Promise<string> {
  const name = tableName('trans')
  if (!hits.trans.length) {
    await insertArrow(name, TRANS_DDL, { variant: v.vidx, part: 'trans' })
    return name
  }
  const m = await packManifest()
  const detail = { variant: v.vidx, part: 'trans' }
  const t0 = performance.now()
  const ipc = hitsTransIPC(hits, { e: m.dof.eqtl, s: m.dof.sqtl })
  performance.measure('pack:hits-trans-build', { start: t0, end: performance.now(), detail })
  const tmp = tableName('tmp_hits')
  await insertArrow(tmp, ipc, { ...detail, part: 'trans-tmp' })
  const { con } = await getDB()
  try {
    await con.query(`CREATE TABLE ${name} AS
      SELECT t.qtl_type,
             CASE WHEN t.qtl_type = 'e' THEN s.gene_id
                  ELSE s.chr || ':' || t.intron_start || ':' || t.intron_end || ':clu_' || t.cluster || '_' || t.strand
                       || ':' || s.gene_id || '.' || s.gene_version END AS phenotype_id,
             s.gene_id, s.symbol, s.chr AS gene_chr, CAST(s.tss AS BIGINT) AS gene_tss,
             ${lit(v.chr)} AS variant_chr, CAST(${v.position} AS INTEGER) AS position,
             ${v.rsNumber ? lit(`rs${v.rsNumber}`) : 'CAST(NULL AS VARCHAR)'} AS rsid,
             CAST(${Number.isNaN(v.af) ? 'NULL' : v.af} AS FLOAT) AS af,
             CAST(t.pval AS DOUBLE) AS pval, CAST(t.beta AS FLOAT) AS beta,
             CAST(t.beta_se AS FLOAT) AS beta_se, CAST(t.r2 AS FLOAT) AS r2
      FROM ${tmp} t JOIN search_index s ON s.ord = t.ord`)
  } finally {
    await dropTable(tmp)
  }
  return name
}
