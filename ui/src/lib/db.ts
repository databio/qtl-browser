/**
 * DuckDB-WASM bootstrap. One database, one connection, shared by plain queries and by the
 * Mosaic coordinator.
 *
 * DuckDB holds only the tables a page materializes from qtlstore reads: a gene's locus window and
 * GWAS window (gene.ts), trans rows (trans.ts), for the plots and the paged tables. Gene lookups,
 * lists and search are plain JS (gene-index.ts), so the engine is started only by a page that
 * draws or pages such a table: the gene and variant pages start it at mount, alongside their
 * first data requests, and nothing starts it on Home, the gene list or a region.
 */
import * as duckdb from '@duckdb/duckdb-wasm'

export type Row = Record<string, unknown>

let dbPromise: Promise<{ db: duckdb.AsyncDuckDB; con: duckdb.AsyncDuckDBConnection }> | null = null
const failureListeners = new Set<(e: Error) => void>()

/** Called with the error when the engine fails to start (App shows it above the page). */
export function onDBFailure(cb: (e: Error) => void): () => void {
  failureListeners.add(cb)
  return () => { failureListeners.delete(cb) }
}

async function boot() {
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
  return { db, con }
}

/** The engine, started on first use; a failed start is reported to `onDBFailure` listeners and
 *  forgotten, so the next call retries. */
export function getDB() {
  if (!dbPromise) {
    const p = boot()
    dbPromise = p
    p.catch((e: Error) => {
      if (dbPromise === p) dbPromise = null
      for (const cb of failureListeners) cb(e)
    })
  }
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
