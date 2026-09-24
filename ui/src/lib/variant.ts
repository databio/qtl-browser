/**
 * Every variant-page read, from the qtlstore (SPEC.md sections 5 and 8). The catalog's variant
 * index (`.qbx`) is one whole-object fetch, the first time a page needs it; it holds every
 * variants-file page offset and each rsID block's first number, so a variant costs a small fixed
 * number of requests:
 *
 *   by rsID:    rsID block(s) (1) -> variants page + hits frame (2, in parallel)
 *   by chr:pos: variants page (1, 2 on a cis miss) -> hits frame (1)
 *
 * The hits file is paged by vidx (SPEC section 8): its header and frame table are read once per
 * chromosome, then the one frame holding the variant (no request when that frame is empty). Each lead row also reads its phenotype's block, for the slope and SE
 * the lead list prints. The cis scan is cis-scan.ts, behind its button.
 */
import type { Row } from './db'
import { geneInfo, phenotypesByOrd } from './gene-index'
import { fetchDecoded, fetchObject, getStore, rangeObject, resultsFile, variantsFile } from './store'
import { countBelow, decodeHitsFrame, decodeHitsTable, decodeRsidRecords, decodeVariantIndex, decodeVariantPage, HEADER_LEN, hitsOf,
  hitsTableLen, HIT_CS, HIT_LEAD, pageOfVidx, RSID_RECORD_LEN, rsidFind, scanBlockRow,
  type HitsFrame, type HitsTable, type VariantIndex, type VariantRecord } from './store-decode'

/** Variants per hits frame the builder writes (SPEC section 8); the header is checked against it. */
const HITS_FRAME_VARIANTS = 1024

let index: Promise<VariantIndex> | null = null

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
    return fetchDecoded('store:variant-page', { chr, page: k }, f.name, off, len,
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
  const hit = await rsidFind(idx.rsidFirst, rs, b => fetchDecoded('store:rsid', { rs, block: b }, name, HEADER_LEN + b * B * RSID_RECORD_LEN,
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

/** The phenotypes the records name, by ord, with their genes' symbols: the search index part of
 *  each phenotype's chromosome and that chromosome's genes (gene-index.ts), cached for the session. */
export async function hitPhenotypes(ords: number[]): Promise<Map<number, HitPhenotype>> {
  if (!ords.length) return new Map()
  const byOrd = await phenotypesByOrd(ords)
  const out = new Map<number, HitPhenotype>()
  await Promise.all([...byOrd.values()].map(async p => {
    if (p.chr == null || p.blk_off == null || p.blk_len == null) return      // trans-only: no block, never a lead or set member
    const g = p.gene_id ? await geneInfo(p.gene_id, p.chr) : null
    out.set(p.ord, { ord: p.ord, phenotype_type: p.phenotype_type, phenotype_id: p.phenotype_id, gene_id: p.gene_id,
      symbol: g?.name ?? null, chr: p.chr, blk_off: p.blk_off, blk_len: p.blk_len, var_start: p.var_start, n_var: p.n_var })
  }))
  return out
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

/** Blocks closer than this in one results file are read as one range (the bytes between are
 *  cheaper than another request). */
const MERGE_GAP = 64 * 1024

/** The variant's nominal slope and SE in each phenotype's block, by ord: the blocks of one results
 *  file are read in runs, a run joining blocks less than MERGE_GAP apart, one request per run. Null
 *  for a phenotype the variant is not a row of; the slope is null when the results set has no dof. */
export async function nominalsAt(ps: HitPhenotype[], vidx: number): Promise<Map<number, { slope: number | null; se: number | null } | null>> {
  const s = await getStore()
  const out = new Map<number, { slope: number | null; se: number | null } | null>()
  const byFile = new Map<string, { f: ReturnType<typeof resultsFile>; ps: HitPhenotype[] }>()
  for (const p of ps) {
    if (p.var_start == null || p.n_var == null || vidx < p.var_start || vidx >= p.var_start + p.n_var) { out.set(p.ord, null); continue }
    const f = resultsFile(s, p.phenotype_type, p.chr)
    const g = byFile.get(f.name)
    if (g) g.ps.push(p); else byFile.set(f.name, { f, ps: [p] })
  }
  const nan = (x: number) => (Number.isNaN(x) ? null : x)
  await Promise.all([...byFile.values()].flatMap(({ f, ps: group }) => {
    group.sort((a, b) => a.blk_off - b.blk_off)
    const runs: HitPhenotype[][] = []
    for (const p of group) {
      const run = runs[runs.length - 1]
      const end = run ? Math.max(...run.map(x => x.blk_off + x.blk_len)) : 0
      if (run && p.blk_off - end <= MERGE_GAP) run.push(p); else runs.push([p])
    }
    return runs.map(async run => {
      const off = run[0].blk_off, len = Math.max(...run.map(x => x.blk_off + x.blk_len)) - off
      await fetchDecoded('store:lead-block', { ords: run.map(p => p.ord) }, f.name, off, len, (b, what) => {
        for (const p of run) {
          const r = scanBlockRow(b.subarray(p.blk_off - off, p.blk_off - off + p.blk_len), vidx - p.var_start!, f.dof, `${what} ord ${p.ord}`)
          out.set(p.ord, { slope: nan(r.slope), se: nan(r.se) })
        }
      })
    })
  }))
  return out
}
