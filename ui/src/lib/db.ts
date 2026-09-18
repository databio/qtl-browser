/**
 * DuckDB-WASM bootstrap. One database, one connection, shared by plain queries and by the
 * Mosaic coordinator.
 *
 * The browser reads no parquet at all. Every page's data comes from pack files fetched with plain
 * Range requests (pack.ts, trans-pack.ts, variant-pack.ts, cis-scan.ts) and decoded in JS;
 * DuckDB holds only in-memory tables built from those, plus `search_index`, which arrives as an
 * Arrow IPC stream in one zstd frame (SPEC section 6) and is inserted directly. So there is no
 * file registration, no HTTP filesystem, and no parquet extension to load.
 */
import * as duckdb from '@duckdb/duckdb-wasm'
import { dataUrl, getManifest } from './manifest'
import { decodeSearchIndex } from './pack-decode'
import { loadGwasIndex } from './pack'
import { startVariantIndex } from './variant-pack'

export type Row = Record<string, unknown>

let dbPromise: Promise<{ db: duckdb.AsyncDuckDB; con: duckdb.AsyncDuckDBConnection }> | null = null

async function boot() {
  // the search index is one plain fetch started now, alongside the wasm download: one request,
  // and no DuckDB involvement until the engine exists
  const indexBytes = getManifest().then(async m => {
    const path = m.packs.search_index
    const r = await fetch(dataUrl(path))
    if (!r.ok) throw new Error(`${path}: ${r.status}`)
    return decodeSearchIndex(new Uint8Array(await r.arrayBuffer()), path)
  })
  // that manifest fetch is the app's only one; the GWAS index and the variant index are plain
  // fetches too, ready by the time a page needs them. They never wait on DuckDB, and a failure
  // surfaces where they are used.
  loadGwasIndex().catch(() => {})
  startVariantIndex()
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
  await con.insertArrowFromIPCStream(await indexBytes, { name: 'search_index', create: true })
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

let tableSeq = 0
/** A fresh in-memory table name; every table the app creates takes its name from here. */
export const tableName = (prefix: string) => `${prefix}_${++tableSeq}`

export async function dropTable(name: string) {
  const { con } = await getDB()
  await con.query(`DROP TABLE IF EXISTS ${name}`).catch(() => {})
}
