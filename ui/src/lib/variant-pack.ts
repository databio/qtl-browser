/**
 * Every variant-page read, from the pack format (SPEC.md sections 13 to 15). The startup file
 * `variant_index.qbx` is one whole-file fetch begun by db.ts boot() and never awaited by boot; it
 * holds every variants-file page offset, every hits frame offset, and each rsID block's first
 * number, so a variant costs a small fixed number of range requests:
 *
 *   by rsID:   rsid_index (1) -> variants page + hits frame (2, in parallel)  = 3
 *   by chr:pos: variants page (1) -> hits frame (1)                           = 2
 *
 * Decoding lives in pack-decode.ts. The cis scan is cis-scan.ts, behind its button.
 */
import { DATA_BASE, getManifest } from './manifest'
import { fetchDecoded, packFile, packManifest } from './pack'
import { decodeHitsFrame, decodeVariantIndex, decodeVariantPage, NLP_MAXQ, BETA_MAXQ, rsidInBlock, VIDX_CHROMS,
  type HitsFrame, type VariantIndex, type VariantRecord } from './pack-decode'

/** One variant's rows in its frame, split into the three lists the page shows. */
export interface Hits {
  frame: HitsFrame
  /** row numbers in `frame`, by kind group: trans (0-1), lead (2-3), credible set (4-5) */
  trans: number[]
  leads: number[]
  cs: number[]
}

let index: Promise<VariantIndex> | null = null

/** Starts the whole-file fetch once; safe to call before anything needs it. */
export function startVariantIndex(): void {
  variantIndex().catch(() => {})
}

/** The startup file (SPEC section 15): one fetch per session. A failure is forgotten so the next
 *  call retries, as the GWAS index does. */
export function variantIndex(): Promise<VariantIndex> {
  if (!index) {
    const p = packManifest().then(async m => {
      const r = await fetch(`${DATA_BASE}/${m.variant_index}`)
      if (!r.ok) throw new Error(`${m.variant_index}: HTTP ${r.status}`)
      const t0 = performance.now()
      const idx = decodeVariantIndex(new Uint8Array(await r.arrayBuffer()), m.variant_index)
      performance.measure('pack:variant-index', { start: t0, end: performance.now(), detail: { part: 'variant-index' } })
      if (idx.pageSize !== m.variant_page_size || idx.frameVariants !== m.hits_frame_variants || idx.rsidBlockRecords !== m.rsid_block_records)
        throw new Error(`${m.variant_index}: page size, frame variants, or rsID block records differ from the manifest`)
      if (idx.rsidNRecords !== (await getManifest()).packs.counts.rsid_index.records)
        throw new Error(`${m.variant_index}: rsid_n_records differs from the manifest`)
      return idx
    })
    index = p
    p.catch(() => { if (index === p) index = null })
  }
  return index
}

function chromIndex(idx: VariantIndex, chr: string) {
  const c = idx.chroms.get(chr)
  if (!c) throw new Error(`variant_index holds no ${chr}`)
  return c
}

/** Index of the last entry of the sorted array at or below `x`, or -1. */
function lastAtOrBelow(a: Uint32Array, x: number, from = 0, to = a.length): number {
  let lo = from, hi = to - 1, best = -1
  while (lo <= hi) {
    const mid = (lo + hi) >> 1
    if (a[mid] <= x) { best = mid; lo = mid + 1 } else hi = mid - 1
  }
  return best
}

// ---- caches: the last 8 pages and the last 8 frames, as pack.ts caches genes -----------------

const CACHE_SIZE = 8

function cache<T>(): { get(key: string, make: () => Promise<T>): Promise<T> } {
  const m = new Map<string, Promise<T>>()
  return {
    get(key, make) {
      const hit = m.get(key)
      if (hit) { m.delete(key); m.set(key, hit); return hit }        // most recent last
      const p = make()
      p.catch(() => { if (m.get(key) === p) m.delete(key) })         // a rejected promise is not kept
      m.set(key, p)
      for (const k of m.keys()) { if (m.size <= CACHE_SIZE) break; m.delete(k) }
      return p
    },
  }
}

const pages = cache<VariantRecord[]>()
const frames = cache<HitsFrame>()

/** Page `k` of a chromosome's variants file (both sections; see SPEC section 15). */
function page(chr: string, k: number): Promise<VariantRecord[]> {
  return pages.get(`${chr}:${k}`, async () => {
    const [m, idx] = [await packManifest(), await variantIndex()]
    const c = chromIndex(idx, chr)
    if (!(k >= 0 && k < c.nPagesCis + c.nPagesTrans)) throw new Error(`${chr}: no variants page ${k}`)
    const off = c.pageOff[k], len = c.pageOff[k + 1] - off
    return fetchDecoded('pack:variant-page', { chr, page: k }, packFile(m.files.variants, chr, 'variants'), off, len,
      (bytes, what) => decodeVariantPage(bytes, chr, c.nCis, what), 'pack:variant-page-decode')
  })
}

/** SPEC section 14: turn an rs number into a chromosome and a vidx with one range request. */
export async function lookupRsid(rs: number): Promise<{ chr: string; vidx: number } | null> {
  if (!Number.isInteger(rs) || rs < 1 || rs > 0xffffffff) return null
  const [m, idx] = [await packManifest(), await variantIndex()]
  const b = lastAtOrBelow(idx.rsidFirst, rs)
  if (b < 0) return null
  const B = idx.rsidBlockRecords
  const n = Math.min(B, idx.rsidNRecords - b * B)
  const got = await fetchDecoded('pack:rsid', { rs, block: b }, m.rsid_index, 32 + b * B * 8, n * 8,
    (bytes, what) => rsidInBlock(bytes, rs, what), 'pack:rsid-decode')
  if (!got) return null
  const chr = VIDX_CHROMS[got.chrOrdinal]
  if (!chr) throw new Error(`${m.rsid_index}: rs${rs} names chromosome ordinal ${got.chrOrdinal}`)
  return { chr, vidx: got.vidx }
}

/** The variant at `vidx`: one page request (cached). */
export async function variantAt(chr: string, vidx: number): Promise<VariantRecord> {
  const idx = await variantIndex()
  const c = chromIndex(idx, chr)
  const P = idx.pageSize
  const k = vidx < c.nCis ? Math.floor(vidx / P) : c.nPagesCis + Math.floor((vidx - c.nCis) / P)
  const recs = await page(chr, k)
  const v = recs[vidx - recs[0].vidx]
  if (!v || v.vidx !== vidx) throw new Error(`${chr} page ${k} does not hold vidx ${vidx}`)
  return v
}

/** The variant at a position: the cis section first, then the trans-only section on a miss.
 *  `(chr, position)` is unique over every variant, so one position is one variant. */
export async function variantAtPosition(chr: string, pos: number): Promise<VariantRecord | null> {
  const idx = await variantIndex()
  const c = chromIndex(idx, chr)
  // the two sections are searched one after the other: cis is the common case, so its miss is the
  // only one that costs a second request
  for (const [from, to] of [[0, c.nPagesCis], [c.nPagesCis, c.nPagesCis + c.nPagesTrans]] as [number, number][]) {
    if (from >= to) continue
    const k = lastAtOrBelow(c.pageFirstPos, pos, from, to)
    if (k < from) continue                                   // before the section's first position
    const hit = (await page(chr, k)).find(r => r.position === pos)
    if (hit) return hit
  }
  return null
}

/** SPEC section 13: the variant's hits frame (cached), split into its three lists. */
export async function loadHits(chr: string, vidx: number): Promise<Hits> {
  const [m, idx] = [await packManifest(), await variantIndex()]
  const c = chromIndex(idx, chr)
  const F = idx.frameVariants
  const g = Math.floor(vidx / F)
  if (!(g >= 0 && g < c.hitsOff.length - 1)) throw new Error(`${chr}: no hits frame ${g}`)
  const frame = await frames.get(`${chr}:${g}`, () => {
    const off = c.hitsOff[g], len = c.hitsOff[g + 1] - off
    return fetchDecoded('pack:hits', { chr, frame: g }, packFile(m.files.hits, chr, 'hits'), off, len,
      (bytes, what) => decodeHitsFrame(bytes, g * F, what), 'pack:hits-decode')
  })
  const slot = vidx - frame.firstVidx
  if (!(slot >= 0 && slot < frame.nVariants)) throw new Error(`${chr}: vidx ${vidx} is not in frame ${g}`)
  const trans: number[] = [], leads: number[] = [], cs: number[] = []
  for (let r = frame.rowStart[slot]; r < frame.rowStart[slot + 1]; r++) {
    const k = frame.kind[r]
    ;(k <= 1 ? trans : k <= 3 ? leads : cs).push(r)
  }
  return { frame, trans, leads, cs }
}

// ---- the values the page prints, decoded from one row's codes (SPEC section 13) --------------

/** A lead row (kinds 2 and 3): perm p, slope and its SE, and whether the phenotype is significant. */
export function leadValues(f: HitsFrame, r: number) {
  return {
    pvalPerm: Math.pow(10, -((f.v1[r] * f.lead.nlpMax) / NLP_MAXQ)),
    slopeSe: (f.v2[r] * f.lead.seMax) / 65535,
    slope: (f.v3[r] * f.lead.slopeMax) / BETA_MAXQ,
    significant: (f.flags[r] & 2) !== 0,
  }
}

/** A credible-set row (kinds 4 and 5). */
export function csValues(f: HitsFrame, r: number) {
  return { pip: f.v1[r] / 65535, csId: f.v2[r] }
}

/** An sQTL row's phenotype id, rebuilt from its intron fields and its gene (SPEC section 13). */
export function hitsPhenotypeId(f: HitsFrame, r: number, gene: { chr: string; gene_id: string; gene_version: number | null }): string {
  return `${gene.chr}:${f.intronStart[r]}:${f.intronEnd[r]}:clu_${f.cluster[r]}_${(f.flags[r] & 1) ? '-' : '+'}:${gene.gene_id}.${gene.gene_version}`
}
