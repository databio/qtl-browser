/**
 * Every variant-page read, from the qtlstore (SPEC.md sections 5 and 8). The catalog's variant
 * index (`.qbx`) is one whole-object fetch begun by db.ts boot() and never awaited by boot; it
 * holds every variants-file page offset and each rsID block's first number, so a variant costs a
 * small fixed number of requests:
 *
 *   by rsID:    rsID block(s) (1) -> variants page + hits frame (2, in parallel)
 *   by chr:pos: variants page (1, 2 on a cis miss) -> hits frame (1)
 *
 * The hits file is paged by vidx (SPEC section 8): its header and frame table are read once per
 * chromosome, then the one frame holding the variant (no request when that frame is empty). Each lead row also reads its phenotype's block, for the slope and SE
 * the lead list prints. The cis scan is cis-scan.ts, behind its button.
 */
import { rows, type Row } from './db'
import { fetchDecoded, fetchObject, getStore, rangeObject, resultsFile, variantsFile } from './store'
import { ALL, countBelow, decodeHitsFrame, decodeHitsTable, decodeRsidRecords, decodeVariantIndex, decodeVariantPage, HEADER_LEN, hitsOf,
  hitsTableLen, HIT_CS, HIT_LEAD, KIND, pageOfVidx, RSID_RECORD_LEN, rsidFind, scanBlockRow,
  type HitsFrame, type HitsTable, type VariantIndex, type VariantRecord } from './store-decode'

/** Variants per hits frame the builder writes (SPEC section 8); the header is checked against it. */
const HITS_FRAME_VARIANTS = 1024

let index: Promise<VariantIndex> | null = null

/** Starts the whole-object fetch once; safe to call before anything needs it. */
export function startVariantIndex(): void {
  variantIndex().catch(() => {})
}

/** The catalog's variant index: one fetch per session, checked against the catalog pointer. A
 *  failure is forgotten so the next call retries. */
export function variantIndex(): Promise<VariantIndex> {
  if (!index) {
    const p = getStore().then(async s => {
      const name = s.catalog.vidx
      const bytes = await fetchObject(name)
      const t0 = performance.now()
      const idx = decodeVariantIndex(bytes, s.catalog.chromosomes.map(c => c.name), s.catalog.collection_digest, name)
      performance.measure('store:variant-index', { start: t0, end: performance.now(), detail: { part: 'variant-index' } })
      if (idx.pageSize !== s.catalog.page_size) throw new Error(`${name}: page size ${idx.pageSize}, catalog says ${s.catalog.page_size}`)
      for (const c of s.catalog.chromosomes) {
        const x = idx.chroms.get(c.name)!
        if (x.nCis !== c.n_cis || x.nCis + x.nTrans !== c.count) throw new Error(`${name}: ${c.name} counts differ from the catalog pointer`)
      }
      return idx
    })
    index = p
    p.catch(() => { if (index === p) index = null })
  }
  return index
}

function chromIndex(idx: VariantIndex, chr: string) {
  const c = idx.chroms.get(chr)
  if (!c) throw new Error(`the variant index holds no ${chr}`)
  return c
}

// ---- caches: the last 8 pages and hits frames; every hits frame table read this session ----------------------------

const CACHE_SIZE = 8

function cache<T>(size: number): { get(key: string, make: () => Promise<T>): Promise<T> } {
  const m = new Map<string, Promise<T>>()
  return {
    get(key, make) {
      const hit = m.get(key)
      if (hit) { m.delete(key); m.set(key, hit); return hit }        // most recent last
      const p = make()
      p.catch(() => { if (m.get(key) === p) m.delete(key) })         // a rejected promise is not kept
      m.set(key, p)
      for (const k of m.keys()) { if (m.size <= size) break; m.delete(k) }
      return p
    },
  }
}

const pages = cache<VariantRecord[]>(CACHE_SIZE)
const hitsTables = cache<{ name: string; table: HitsTable }>(64)
const hitsFrames = cache<HitsFrame>(CACHE_SIZE)

/** Page `k` of a chromosome's variants file (both sections). */
function page(chr: string, k: number): Promise<VariantRecord[]> {
  return pages.get(`${chr}:${k}`, async () => {
    const [s, idx] = [await getStore(), await variantIndex()]
    const c = chromIndex(idx, chr)
    if (!(k >= 0 && k < c.nPagesCis + c.nPagesTrans)) throw new Error(`${chr}: no variants page ${k}`)
    const off = c.pageOff[k], len = c.pageOff[k + 1] - off
    const f = variantsFile(s, chr)
    return fetchDecoded('store:variant-page', { chr, page: k }, f.name, f.expect, off, len,
      (bytes, what) => decodeVariantPage(bytes, chr, c.nCis, what), 'store:variant-page-decode')
  })
}

/** An rs number to a chromosome and vidx, usually with one range request (`rsidFind`). rs_number
 *  may repeat in v1 (dbSNP puts several allele pairs under one rsID); the page shows the first site
 *  of the run, in (ordinal, vidx) order. */
export async function lookupRsid(rs: number): Promise<{ chr: string; vidx: number } | null> {
  if (!Number.isInteger(rs) || rs < 1 || rs > 0xffffffff) return null
  const [s, idx] = [await getStore(), await variantIndex()]
  const B = idx.rsidBlockRecords
  const name = s.catalog.rsid
  const hit = await rsidFind(idx.rsidFirst, rs, b => fetchDecoded('store:rsid', { rs, block: b }, name,
    { kind: KIND.rsid, chrom: ALL, seqDigest: s.catalog.collection_digest, count: idx.rsidN }, HEADER_LEN + b * B * RSID_RECORD_LEN,
    Math.min(B, idx.rsidN - b * B) * RSID_RECORD_LEN, (bytes, what) => decodeRsidRecords(bytes, what), 'store:rsid-decode'))
  if (!hit) return null
  const chr = idx.names[hit.ordinal - 1]
  if (!chr) throw new Error(`${name}: rs${rs} names chromosome ordinal ${hit.ordinal}`)
  return { chr, vidx: hit.vidx }
}

/** The variant at `vidx`: one page request (cached). */
export async function variantAt(chr: string, vidx: number): Promise<VariantRecord> {
  const idx = await variantIndex()
  const c = chromIndex(idx, chr)
  const k = pageOfVidx(c, idx.pageSize, vidx)
  const recs = await page(chr, k)
  const v = recs[vidx - recs[0].vidx]
  if (!v || v.vidx !== vidx) throw new Error(`${chr} page ${k} does not hold vidx ${vidx}`)
  return v
}

/** The first variant at a position: the cis section first, then the trans-only section on a miss.
 *  A position can hold several sites (different alleles); the first in vidx order is shown. */
export async function variantAtPosition(chr: string, pos: number): Promise<VariantRecord | null> {
  const idx = await variantIndex()
  const c = idx.chroms.get(chr)
  if (!c) return null
  for (const [from, to] of [[0, c.nPagesCis], [c.nPagesCis, c.nPagesCis + c.nPagesTrans]] as [number, number][]) {
    if (from >= to) continue
    // the last page starting at or below pos; sites at pos may begin on the page before it when
    // that page ends with the same position
    const k = countBelow(c.pageFirstPos, pos, true, from, to) - 1
    if (k < from) continue
    if (k > from && c.pageFirstPos[k] === pos) {
      const hit = (await page(chr, k - 1)).find(r => r.position === pos)
      if (hit) return hit
    }
    const hit = (await page(chr, k)).find(r => r.position === pos)
    if (hit) return hit
  }
  return null
}

/** One chromosome's hits file header and frame table: one range request, once per session. */
function hitsTable(chr: string): Promise<{ name: string; table: HitsTable }> {
  return hitsTables.get(chr, async () => {
    const s = await getStore()
    const name = s.experiment.hits[chr]
    const c = s.chroms.get(chr)
    if (!name || !c) throw new Error(`experiment ${s.experiment.id}: no hits file for ${chr}`)
    // sized from the catalog's variant count and the builder's frame size; the header confirms both
    const t0 = performance.now()
    const bytes = await rangeObject(name, 0, hitsTableLen(c.count, HITS_FRAME_VARIANTS))
    performance.measure('store:hits-table', { start: t0, end: performance.now(), detail: { chr } })
    return { name, table: decodeHitsTable(bytes, { chrom: chr, seqDigest: c.seq_digest, nVariants: c.count }, name) }
  })
}

/** One variant's hits records, as indices into `frame`: group leads, credible-set members, and
 *  trans associations. */
export interface Hits { frame: HitsFrame; leads: number[]; cs: number[]; trans: number[] }

export async function loadHits(chr: string, vidx: number): Promise<Hits> {
  const { name, table } = await hitsTable(chr)
  const F = table.frameVariants
  const g = Math.floor(vidx / F)
  if (!(g >= 0 && g < table.frameOff.length - 1)) throw new Error(`${chr}: no hits frame ${g}`)
  const off = table.frameOff[g], len = table.frameOff[g + 1] - off
  const frame = await hitsFrames.get(`${chr}:${g}`, async () => {
    if (!len) return decodeHitsFrame(new Uint8Array(0), g * F, F)          // an empty frame: no request
    // the header was checked with the table, so the frame is a plain range read
    const t0 = performance.now()
    const bytes = await rangeObject(name, off, len)
    const t1 = performance.now()
    performance.measure('store:hits', { start: t0, end: t1, detail: { chr, frame: g } })
    const out = decodeHitsFrame(bytes, g * F, F, `${name} frame ${g}`)
    performance.measure('store:hits-decode', { start: t1, end: performance.now(), detail: { chr, part: 'hits' } })
    return out
  })
  const [lo, hi] = hitsOf(frame, vidx)
  const leads: number[] = [], cs: number[] = [], trans: number[] = []
  for (let r = lo; r < hi; r++) (frame.kind[r] === HIT_LEAD ? leads : frame.kind[r] === HIT_CS ? cs : trans).push(r)
  return { frame, leads, cs, trans }
}

/** A `phenotypes` row with its gene's symbol, for the lead and credible-set lists. */
export interface HitPhenotype extends Row {
  ord: number; phenotype_type: string; phenotype_id: string; gene_id: string | null; symbol: string | null; chr: string
  blk_off: number; blk_len: number; var_start: number | null; n_var: number | null
}

/** The phenotypes the records name, by ord (one local query). */
export async function hitPhenotypes(ords: number[]): Promise<Map<number, HitPhenotype>> {
  if (!ords.length) return new Map()
  const got = await rows<HitPhenotype>(`SELECT p.ord, p.phenotype_type, p.phenotype_id, p.gene_id, g.name AS symbol, p.chr,
      p.blk_off, p.blk_len, p.var_start, p.n_var
    FROM phenotypes p LEFT JOIN genes g USING (gene_id) WHERE p.ord IN (${[...new Set(ords)].map(Number).join(',')})`)
  return new Map(got.map(g => [g.ord, g]))
}

/** A lead record: its permutation p and significance, from the hits file. */
export function leadValues(h: Hits, r: number) {
  const v = h.frame.value[r]
  return { pvalPerm: Number.isNaN(v) ? null : v, significant: (h.frame.flags[r] & 1) !== 0 }
}

/** A credible-set record: PIP and set id. */
export function csValues(h: Hits, r: number) {
  return { pip: h.frame.value[r], csId: h.frame.csId[r] }
}

/** The variant's nominal slope and SE in one phenotype's block: one range request for the block.
 *  Null when the variant is not one of the block's rows; the slope is null when the results set
 *  has no dof. */
export async function nominalAt(p: HitPhenotype, vidx: number): Promise<{ slope: number | null; se: number | null } | null> {
  if (p.var_start == null || p.n_var == null || vidx < p.var_start || vidx >= p.var_start + p.n_var) return null
  const s = await getStore()
  const f = resultsFile(s, p.phenotype_type, p.chr)
  const r = await fetchDecoded('store:lead-block', { ord: p.ord }, f.name, f.expect, p.blk_off, p.blk_len,
    (b, what) => scanBlockRow(b, vidx - p.var_start!, f.dof, what))
  const nan = (x: number) => (Number.isNaN(x) ? null : x)
  return { slope: nan(r.slope), se: nan(r.se) }
}
