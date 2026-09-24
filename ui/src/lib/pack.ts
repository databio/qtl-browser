/**
 * Every gene-page read, from the pack format (SPEC.md). Opening a gene starts three plain Range
 * requests in the same tick, with byte offsets from `search_index` and file names and dof from
 * `manifest.packs`: the gene's eQTL block (details, eQTL rows, credible sets), its variants range
 * (covering its eQTL run and every intron run), and its GWAS window (through the GWAS index loaded at
 * startup). Each sQTL intron the page shows adds one request for its block, and the gene page
 * starts one more for the gene's trans frame (trans-pack.ts). Decoding lives in pack-decode.ts;
 * rows reach DuckDB as Arrow IPC streams (insertArrow).
 */
import { makeVector, Table, tableToIPC, Utf8, vectorFromArray } from 'apache-arrow'
import { dropTable, getDB, tableName } from './db'
import { DATA_BASE, getManifest, type PackManifest } from './manifest'
import { csMembers, decodeGwasRange, decodeResultBlock, decodeVariantRange, gwasRange, parseGwasIndex, readerColumns, sliceRun,
  type GwasColumns, type GwasIndex, type ResultBlock, type VariantRange } from './pack-decode'
import type { CredibleSetRow, Gene, GeneDetail, SearchHit, SplicePhenotype } from './queries'

export interface GenePack {
  hit: SearchHit
  /** request 1: the gene's eQTL block, `blk_off`/`blk_len` in packs.files.eqtl[chr] */
  block: Promise<ResultBlock>
  /** request 2: the gene's variants range, `var_off`/`var_len` in packs.files.variants[chr] */
  variants: Promise<VariantRange>
  /** request 3: the GWAS rows in [w_lo, w_hi] as a DuckDB table (empty on chrX); dropped when the gene leaves the cache */
  gwas: Promise<string>
  /** the block's details: no extra request */
  detail: Promise<GeneDetail>
  /** +1 request per intron, memoized */
  intron(phenotypeId: string): Promise<ResultBlock>
}

/** Readers accept exactly the format the manifest names. */
export async function packManifest(): Promise<PackManifest> {
  const p = (await getManifest()).packs
  if (p?.format !== 'qtlb' || p.version !== 0) throw new Error(`manifest.json: packs are ${p?.format} version ${p?.version}; this reader reads qtlb version 0`)
  return p
}

export function packFile(files: Record<string, string>, chr: string, kind: string): string {
  const f = files[chr]
  if (!f) throw new Error(`manifest.json: no ${kind} pack for ${chr}`)
  return f
}

async function rangeFetch(path: string, offset: number, length: number): Promise<Uint8Array> {
  const range = `bytes=${offset}-${offset + length - 1}`
  const r = await fetch(`${DATA_BASE}/${path}`, { headers: { Range: range } })
  // a 200 means the server ignored Range and is sending the whole pack: fail instead of downloading it
  if (r.status !== 206) throw new Error(`${path} ${range}: expected 206, got ${r.status}`)
  const buf = new Uint8Array(await r.arrayBuffer())
  if (buf.length !== length) throw new Error(`${path} ${range}: got ${buf.length} bytes, expected ${length}`)
  return buf
}

/** One range request timed as `mark` (pack:block, pack:variants, pack:gwas, pack:intron, pack:trans), then
 *  its decode as `decodeMark`. */
export async function fetchDecoded<T>(mark: string, detail: Record<string, unknown>, path: string, offset: number, length: number,
  decode: (bytes: Uint8Array, what: string) => T, decodeMark = 'pack:decode'): Promise<T> {
  const t0 = performance.now()
  const bytes = await rangeFetch(path, offset, length)
  const t1 = performance.now()
  performance.measure(mark, { start: t0, end: t1, detail })
  const out = decode(bytes, `${path} bytes ${offset}+${length}`)
  performance.measure(decodeMark, { start: t1, end: performance.now(), detail: { ...detail, part: mark.slice('pack:'.length) } })
  return out
}

/** Creates the DuckDB table `name` from an Arrow IPC stream, or empty from a column list
 *  (`position INTEGER, ...`). The DuckDB worker runs one request at a time, so a no-op query first
 *  separates waiting for it (another table's insert, say) from the insert itself: measured as
 *  pack:worker-wait and pack:insert, with `detail.part` naming the table. */
export async function insertArrow(name: string, source: Uint8Array | string, detail: Record<string, unknown>): Promise<void> {
  const { con } = await getDB()
  const tw = performance.now()
  await con.query('SELECT 1')
  const t0 = performance.now()
  performance.measure('pack:worker-wait', { start: tw, end: t0, detail })
  // an IPC stream, not insertArrowTable: duckdb-wasm serializes a Table with its own apache-arrow copy (^17, the app has ^21)
  if (typeof source === 'string') await con.query(`CREATE TABLE ${name} (${source})`)
  else await con.insertArrowFromIPCStream(source, { name, create: true })
  performance.measure('pack:insert', { start: t0, end: performance.now(), detail })
}

let gwasIndex: Promise<GwasIndex> | null = null
/** The GWAS index (SPEC section 11): one fetch per session, started by db.ts boot(). */
export function loadGwasIndex(): Promise<GwasIndex> {
  if (!gwasIndex) {
    const p = packManifest().then(async m => {
      const r = await fetch(`${DATA_BASE}/${m.gwas_index}`)
      if (!r.ok) throw new Error(`${m.gwas_index}: HTTP ${r.status}`)
      const index = parseGwasIndex(new Uint8Array(await r.arrayBuffer()), m.gwas_index)
      if (index.blockRows !== m.gwas_block_rows) throw new Error(`${m.gwas_index}: ${index.blockRows} rows per block, manifest says ${m.gwas_block_rows}`)
      return index
    })
    gwasIndex = p
    p.catch(() => { if (gwasIndex === p) gwasIndex = null })
  }
  return gwasIndex
}

/** What the gene page shows besides the locus, from the block's details JSON (SPEC section 5). */
function detailFromBlock(hit: SearchHit, block: ResultBlock): GeneDetail {
  const d = block.details as { gene?: Record<string, unknown>; exons?: [number, number][]; splice?: Record<string, unknown>[] }
  if (!d.gene || !Array.isArray(d.exons) || !Array.isArray(d.splice)) throw new Error(`${hit.gene_id}: details lack gene, exons, or splice`)
  return {
    gene: { ...d.gene, chr: hit.chr } as unknown as Gene,
    exons: d.exons.map(([start, end]) => ({ start, end })),
    splice: d.splice.map(p => ({ ...p, gene_id: hit.gene_id, symbol: hit.symbol, chr: hit.chr, tss: hit.tss }) as unknown as SplicePhenotype),
  }
}

/** The column names and types the locus SQL reads from a GWAS window. */
const GWAS_DDL = 'position INTEGER, ea VARCHAR, nea VARCHAR, beta FLOAT, se FLOAT, p DOUBLE, eaf FLOAT, rsid VARCHAR, n INTEGER'

/** Request 3: decode the gene's GWAS window once and insert it; a chromosome without GWAS rows gets the same table, empty. */
async function gwasTable(hit: SearchHit, m: PackManifest, index: GwasIndex): Promise<string> {
  const detail = { gene_id: hit.gene_id }
  const file = m.files.gwas[hit.chr]
  if (!file !== !index.chroms.has(hit.chr)) throw new Error(`${hit.chr}: the manifest's GWAS files and the GWAS index disagree`)
  const lo = hit.w_lo, hi = hit.w_hi
  let cols: GwasColumns | null = null
  if (file && lo != null && hi != null) {
    const range = gwasRange(index, hit.chr, lo, hi)
    if (range) cols = await fetchDecoded('pack:gwas', detail, file, range.off, range.len, (b, what) => decodeGwasRange(b, index, hit.chr, range, lo, hi, what))
  }
  const ipc = cols?.rows ? tableToIPC(new Table({
    position: makeVector(cols.position), ea: vectorFromArray(cols.ea, new Utf8()), nea: vectorFromArray(cols.nea, new Utf8()),
    beta: makeVector(Float32Array.from(cols.beta)), se: makeVector(Float32Array.from(cols.se)), p: makeVector(cols.p),
    eaf: makeVector(Float32Array.from(cols.eaf)),
    rsid: vectorFromArray(Array.from(cols.rsNumber, x => (x ? `rs${x}` : null)), new Utf8()),
    n: makeVector(cols.n),
  }), 'stream') : null
  const name = tableName('gwas')
  await insertArrow(name, ipc ?? GWAS_DDL, { ...detail, part: 'gwas' })
  return name
}

function openGene(hit: SearchHit): GenePack {
  const m = packManifest()
  const detail0 = { gene_id: hit.gene_id }
  const block = m.then(p => {
    if (hit.blk_off == null || hit.blk_len == null) throw new Error(`${hit.gene_id}: search_index has no block for this gene`)
    return fetchDecoded('pack:block', detail0, packFile(p.files.eqtl, hit.chr, 'eQTL'), hit.blk_off, hit.blk_len,
      (b, what) => decodeResultBlock(b, 2, p.dof.eqtl, hit, what))
  })
  const variants = m.then(p => {
    if (hit.var_off == null || hit.var_len == null) throw new Error(`${hit.gene_id}: search_index has no variants range for this gene`)
    return fetchDecoded('pack:variants', detail0, packFile(p.files.variants, hit.chr, 'variants'), hit.var_off, hit.var_len,
      (b, what) => decodeVariantRange(b, hit, what))
  })
  // w_lo and w_hi come from search_index, so the GWAS request waits on nothing the other two fetch
  const gwas = Promise.all([m, loadGwasIndex()]).then(([p, index]) => gwasTable(hit, p, index))
  const detail = block.then(b => detailFromBlock(hit, b))
  const introns = new Map<string, Promise<ResultBlock>>()
  const intron = (id: string) => {
    const have = introns.get(id)
    if (have) return have
    const p = Promise.all([m, detail]).then(([pk, d]) => {
      const s = d.splice.find(x => x.phenotype_id === id)
      if (!s) throw new Error(`${hit.gene_id}: ${id} is not one of the gene's tested introns`)
      return fetchDecoded('pack:intron', { ...detail0, phenotype_id: id }, packFile(pk.files.sqtl, hit.chr, 'sQTL'), s.blk_off, s.blk_len,
        (b, what) => decodeResultBlock(b, 3, pk.dof.sqtl, { blk_len: s.blk_len }, what))
    })
    introns.set(id, p)
    p.catch(() => { if (introns.get(id) === p) introns.delete(id) })
    return p
  }
  // consumers see each rejection when they await; these handlers only keep an unused one quiet
  for (const p of [block, variants, gwas, detail]) p.catch(() => {})
  return { hit, block, variants, gwas, detail, intron }
}

// the last few genes stay open, so tab switches and back navigation send no request
const CACHE_SIZE = 4
const genes = new Map<string, GenePack>()
const release = (gp: GenePack) => { gp.gwas.then(dropTable, () => {}) }

/** Opens a gene (starting requests 1 to 3) or returns it from the cache. */
export function loadGene(hit: SearchHit): GenePack {
  const key = hit.gene_id
  const have = genes.get(key)
  if (have) { genes.delete(key); genes.set(key, have); return have }
  const gp = openGene(hit)
  genes.set(key, gp)
  const forget = () => { if (genes.get(key) === gp) { genes.delete(key); release(gp) } }
  for (const p of [gp.block, gp.variants, gp.gwas]) p.catch(forget)
  while (genes.size > CACHE_SIZE) {
    const [k, old] = genes.entries().next().value as [string, GenePack]
    genes.delete(k)
    release(old)
  }
  return gp
}

const blockFor = (gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<ResultBlock> =>
  qtlType === 'e' ? gp.block
    : phenotypeId ? gp.intron(phenotypeId) : Promise.reject(new Error(`${gp.hit.gene_id}: an sQTL locus needs an intron`))

/** The in-band null codes as SQL NULLs, so `locusSQL` reads the pack rows with the parquet column types. */
const packSource = (raw: string) => `(SELECT position, A1, A2,
  CASE WHEN rs_number = 0 THEN NULL ELSE rs_number::BIGINT END AS rs_number,
  tss_distance,
  CASE WHEN isnan(af) THEN NULL ELSE af END AS af,
  CASE WHEN ma_samples < 0 THEN NULL ELSE ma_samples END AS ma_samples,
  CASE WHEN ma_count < 0 THEN NULL ELSE ma_count END AS ma_count,
  CASE WHEN isnan(pval_nominal) THEN NULL ELSE pval_nominal END AS pval_nominal,
  CASE WHEN isnan(slope) THEN NULL ELSE slope END AS slope,
  CASE WHEN isnan(slope_se) THEN NULL ELSE slope_se END AS slope_se,
  CASE WHEN cs_id < 0 THEN NULL ELSE pip END AS pip,
  CASE WHEN cs_id < 0 THEN NULL ELSE cs_id END AS cs_id
  FROM ${raw})`

/** One cis window as a table for the plots: -log10 p, credible-set class, a tooltip label, and
 *  the DCM GWAS statistics for variants present there (matched on position and alleles in
 *  either orientation, GWAS beta re-signed to the QTL effect allele A1). Both inputs already hold
 *  exactly the phenotype's rows and the window's GWAS rows. Ordered so credible-set variants are
 *  drawn last (on top). */
const locusSQL = (qtl: string, gwas: string) => `
    SELECT q.position,
           -- p underflows to 0 for a few extreme variants: place them just above the largest finite value
           coalesce(-log10(nullif(q.pval_nominal, 0)), max(-log10(nullif(q.pval_nominal, 0))) OVER () * 1.05) AS nlp,
           q.pval_nominal = 0 AS clipped,
           q.pval_nominal, q.slope, q.slope_se, q.af, q.pip, q.cs_id, q.rs_number, q.A1, q.A2,
           q.tss_distance, q.ma_samples, q.ma_count,
           coalesce(q.cs_id::VARCHAR, 'none') AS cs,
           g.p AS gwas_p, -log10(g.p) AS gwas_nlp,
           CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END AS gwas_beta,
           coalesce('rs' || q.rs_number, q.position::VARCHAR) || '  ' || q.A1 || '/' || q.A2
             || chr(10) || CASE WHEN q.pval_nominal = 0 THEN 'p = 0 (underflow; drawn above the maximum)' ELSE 'p = ' || format('{:.2e}', q.pval_nominal) END
             -- the slope is rebuilt from p and SE, so a p = 0 row has none; a NULL here would NULL the whole label
             || chr(10) || CASE WHEN q.slope IS NULL THEN 'slope not recoverable (p underflow)' ELSE 'slope ' || format('{:.3f}', q.slope) || ' ± ' || format('{:.3f}', q.slope_se) END
             || chr(10) || 'AF ' || format('{:.3f}', q.af)
             || CASE WHEN q.pip IS NULL THEN '' ELSE chr(10) || 'PIP ' || format('{:.3f}', q.pip) || ' (set ' || q.cs_id || ')' END
             || CASE WHEN g.p IS NULL THEN '' ELSE chr(10) || 'DCM GWAS p = ' || format('{:.2e}', g.p) || ', beta ' || format('{:+.3f}', CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END) || ' (A1 as effect allele)' END AS label
    FROM ${qtl} q
    LEFT JOIN ${gwas} g
      ON g.position = q.position AND ((g.ea = q.A1 AND g.nea = q.A2) OR (g.ea = q.A2 AND g.nea = q.A1))
    -- the GWAS lists some indels in both allele orientations as separate records: keep one per
    -- QTL variant, preferring the orientation that matches the QTL alleles as written
    QUALIFY row_number() OVER (PARTITION BY q.position, q.A1, q.A2 ORDER BY (g.ea = q.A1) DESC NULLS LAST, g.p) = 1
    ORDER BY q.cs_id IS NOT NULL, q.position`

/** The locus table of the gene's eQTL rows or one intron's sQTL rows, joined to the gene's GWAS
 *  window, as an in-memory DuckDB table; the caller drops it. */
export async function locusTable(gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<string> {
  const [block, range, gwas] = await Promise.all([blockFor(gp, qtlType, phenotypeId), gp.variants, gp.gwas])
  const detail = { gene_id: gp.hit.gene_id, phenotype_id: phenotypeId }
  const t0 = performance.now()
  if (block.varStart == null) throw new Error(`${gp.hit.gene_id}: the eQTL block has no rows`)
  const c = readerColumns(block, sliceRun(range, block.varStart, block.nRows, block, `${gp.hit.gene_id} variants range`))
  const ipc = tableToIPC(new Table({
    position: makeVector(c.position), A1: vectorFromArray(c.a1, new Utf8()), A2: vectorFromArray(c.a2, new Utf8()),
    rs_number: makeVector(c.rsNumber), tss_distance: makeVector(c.tssDistance), af: makeVector(c.af),
    ma_samples: makeVector(c.maSamples), ma_count: makeVector(c.maCount), pval_nominal: makeVector(c.pval),
    slope: makeVector(c.slope), slope_se: makeVector(c.se), pip: makeVector(c.pip), cs_id: makeVector(c.csId),
  }), 'stream')
  const t1 = performance.now()
  performance.measure('pack:decode', { start: t0, end: t1, detail: { ...detail, part: 'columns' } })
  const { con } = await getDB()
  const q = tableName('q')
  try {
    await insertArrow(q, ipc, { ...detail, part: 'locus' })
    const name = tableName('locus')
    await con.query(`CREATE TABLE ${name} AS ${locusSQL(packSource(q), gwas)}`)
    return name
  } finally {
    await con.query(`DROP TABLE IF EXISTS ${q}`).catch(() => {})
  }
}

/** Every credible-set membership of the gene's eQTL rows or one intron, in the credible-set table's
 *  row shape; a variant in two sets is listed under both. */
export async function credibleSets(gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<CredibleSetRow[]> {
  const [block, range] = await Promise.all([blockFor(gp, qtlType, phenotypeId), gp.variants])
  if (block.varStart == null) return []
  const run = sliceRun(range, block.varStart, block.nRows, block, `${gp.hit.gene_id} variants range`)
  return csMembers(block)
    .map(({ row, pip, csId }) => ({
      qtl_type: qtlType, phenotype_id: phenotypeId ?? gp.hit.gene_id, chr: gp.hit.chr,
      position: run.position[row], A1: run.a1[row], A2: run.a2[row], rsid: run.rsNumber[row] ? `rs${run.rsNumber[row]}` : null,
      af: run.af[row], cs_id: csId, pip,
    }))
    .sort((a, b) => a.cs_id - b.cs_id || b.pip - a.pip)
}
