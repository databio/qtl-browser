/**
 * DuckDB-WASM bootstrap. One database, one connection, shared by plain queries and by the
 * Mosaic coordinator.
 *
 * The browser reads no parquet at all. Every page's data comes from qtlstore objects fetched with
 * plain Range requests (gene.ts, variant.ts, cis-scan.ts) and decoded in JS; DuckDB holds only
 * in-memory tables built from those, plus three from Arrow objects inserted directly at boot:
 *
 * - `phenotypes`: the experiment's search index, one row per phenotype (SPEC section 8);
 * - `genes`: the annotation's gene table (SPEC section 6);
 * - `search_index`: one row per annotated gene on a catalog chromosome, the join the gene list,
 *   search, region and gene pages read (its eQTL phenotype's block and run, eGene and sQTL counts).
 */
import * as duckdb from '@duckdb/duckdb-wasm'
import { EQTL_TYPE, genesIPC, getStore, searchIndexIPC, SQTL_TYPE } from './store'
import { startVariantIndex } from './variant'
import { loadGwasIndex } from './gene'

export type Row = Record<string, unknown>

let dbPromise: Promise<{ db: duckdb.AsyncDuckDB; con: duckdb.AsyncDuckDBConnection }> | null = null

async function boot() {
  // the pointers, the search index and the gene table are plain fetches started now, alongside
  // the wasm download; the variant index is too, and is ready by the time a variant page needs it
  const tables = Promise.all([getStore(), searchIndexIPC(), genesIPC()])
  startVariantIndex()
  loadGwasIndex().catch(() => {})
  // wasm and worker from jsDelivr (the 36 MB module is over the Workers asset limit). The
  // worker script is cross-origin, so it is loaded through a same-origin blob shim.
  const bundle = await duckdb.selectBundle(duckdb.getJsDelivrBundles())
  const workerUrl = URL.createObjectURL(new Blob([`importScripts("${bundle.mainWorker}");`], { type: 'text/javascript' }))
  const worker = new Worker(workerUrl)
  const logger = new duckdb.VoidLogger()
  const db = new duckdb.AsyncDuckDB(logger, worker)
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker)
  URL.revokeObjectURL(workerUrl)
  const con = await db.connect()
  // nothing here needs an extension; fail loudly rather than reach out to extensions.duckdb.org
  await con.query(`SET autoinstall_known_extensions = false`).catch(() => {})
  await con.query(`SET autoload_known_extensions = false`).catch(() => {})
  const [s, index, genes] = await tables
  await con.insertArrowFromIPCStream(index, { name: 'phenotypes', create: true })
  await con.insertArrowFromIPCStream(genes, { name: 'genes', create: true })
  const chroms = [...s.chroms.keys()].map(lit).join(', ')
  // Trans-only phenotypes (chr and blk_* null) have no cis block, so they are left out here: a
  // gene with only trans results is "not tested" on its page, as in v0. `tested`: the gene has an
  // eQTL phenotype with rows; `is_egene`: that phenotype passes the
  // experiment's significance rule (null when not eQTL-tested); `has_results`: any phenotype of
  // the gene has a block, so the gene page has something to show
  await con.query(`CREATE TABLE search_index AS
    WITH e AS (SELECT * FROM phenotypes WHERE phenotype_type = ${lit(EQTL_TYPE)} AND blk_off IS NOT NULL),
         s AS (SELECT gene_id, count(*)::INTEGER AS n_sqtl, (count(*) FILTER (WHERE significant))::INTEGER AS n_sqtl_sig
               FROM phenotypes WHERE phenotype_type = ${lit(SQTL_TYPE)} AND blk_off IS NOT NULL GROUP BY gene_id)
    SELECT g.gene_id, g.name AS symbol, g.chr, g.tss, g.start, g."end", g.strand, g.biotype, g.version AS gene_version,
           e.ord, e.blk_off, e.blk_len, e.var_start, e.n_var, e.var_off, e.var_len, e.w_lo, e.w_hi,
           e.n_var IS NOT NULL AS tested,
           CASE WHEN e.ord IS NULL THEN NULL ELSE e.significant END AS is_egene,
           coalesce(s.n_sqtl_sig, 0) AS n_sqtl_sig, coalesce(s.n_sqtl, 0) AS n_sqtl,
           e.ord IS NOT NULL OR s.n_sqtl IS NOT NULL AS has_results
    FROM genes g LEFT JOIN e USING (gene_id) LEFT JOIN s USING (gene_id)
    WHERE g.chr IN (${chroms})`)
  return { db, con }
}

export function getDB() {
  if (!dbPromise) dbPromise = boot()
  return dbPromise
}

/** Run SQL and return plain JS objects (BigInt -> number). */
export async function rows<T extends Row = Row>(sql: string): Promise<T[]> {
  const { con } = await getDB()
  const table = await con.query(sql)
  const out: T[] = []
  for (const r of table) {
    const o: Row = {}
    for (const [k, v] of Object.entries(r.toJSON())) o[k] = typeof v === 'bigint' ? Number(v) : v
    out.push(o as T)
  }
  return out
}

export async function one<T extends Row = Row>(sql: string): Promise<T | null> {
  const r = await rows<T>(sql)
  return r[0] ?? null
}

/** SQL string literal escaping for the few user-controlled strings we interpolate. */
export function lit(s: string): string {
  return `'${s.replace(/'/g, "''")}'`
}

// ---- Mosaic ---------------------------------------------------------------------------------
import { coordinator, wasmConnector } from '@uwdata/mosaic-core'

let mosaicReady: Promise<void> | null = null

/** Point the global Mosaic coordinator at our DuckDB instance (once). */
export function getCoordinator() {
  if (!mosaicReady) {
    mosaicReady = getDB().then(({ db, con }) => {
      coordinator().databaseConnector(wasmConnector({ duckdb: db, connection: con }))
    })
  }
  return mosaicReady.then(() => coordinator())
}

/** Creates the DuckDB table `name` from an Arrow IPC stream, or empty from a column list
 *  (`position INTEGER, ...`). The DuckDB worker runs one request at a time, so a no-op query first
 *  separates waiting for it (another table's insert, say) from the insert itself: measured as
 *  store:worker-wait and store:insert, with `detail.part` naming the table. */
export async function insertArrow(name: string, source: Uint8Array | string, detail: Record<string, unknown>): Promise<void> {
  const { con } = await getDB()
  const tw = performance.now()
  await con.query('SELECT 1')
  const t0 = performance.now()
  performance.measure('store:worker-wait', { start: tw, end: t0, detail })
  // an IPC stream, not insertArrowTable: duckdb-wasm serializes a Table with its own apache-arrow copy (^17, the app has ^21)
  if (typeof source === 'string') await con.query(`CREATE TABLE ${name} (${source})`)
  else await con.insertArrowFromIPCStream(source, { name, create: true })
  performance.measure('store:insert', { start: t0, end: performance.now(), detail })
}

let tableSeq = 0
/** A fresh in-memory table name; every table the app creates takes its name from here. */
export const tableName = (prefix: string) => `${prefix}_${++tableSeq}`

export async function dropTable(name: string) {
  const { con } = await getDB()
  await con.query(`DROP TABLE IF EXISTS ${name}`).catch(() => {})
}
