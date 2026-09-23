/**
 * Decoder for the qtlb format, version 1 (SPEC.md at the repo root): the 64-byte file header,
 * variant pages, the variant index, rsID index records, result blocks with their v1 details, paged
 * hits files, trans frames, GWAS blocks and the GWAS index, and the inverse Student t that rebuilds
 * a cis slope and a trans SE. The Python reference decoders are `pipeline/catalog.py`
 * (`decode_file`, `decode_vidx`, `decode_rsid`, `rsid_lookup`), `pipeline/results.py`
 * (`read_block`, `read_trans`, `hits_frame`, `decode_hits`), `pipeline/gwas.py` (`read_window`, `read_bins`)
 * and the codec `pipeline/packfmt_v1.py`.
 *
 * Pure functions only: no DuckDB, no import.meta.env, no relative imports, and only erasable
 * TypeScript, so `npm run store-check` (Node) runs this exact file against the Python decoders.
 * Bytes that break a SPEC rule throw a PackError naming the object and the rule.
 */
import { decompress } from 'fzstd'
import { tableFromIPC } from 'apache-arrow'

export class PackError extends Error {}

/** The one format version this file reads (`store.json` `format_version`, header byte 5). */
export const FORMAT_VERSION = 1
/** Header kinds v1 writes (SPEC section 3). */
export const KIND = { variants: 1, results: 2, gwas: 4, gwasIndex: 5, trans: 6, hits: 7, rsid: 8, variantIndex: 9 } as const
export const HEADER_LEN = 64
/** Header chromosome of the catalog-wide objects (variant index, rsID index). */
export const ALL = 'all'
/** `^[A-Za-z0-9_-]{32}$`: a sha512t24u digest. */
export const DIGEST = /^[A-Za-z0-9_-]{32}$/
/** `<digest>.<ext>`, the only names a pointer may give an object (SPEC section 3). */
export const OBJECT_NAME = /^[A-Za-z0-9_-]{32}\.[a-z0-9]+(\.[a-z0-9]+)*$/

if (new Uint8Array(new Uint16Array([1]).buffer)[0] !== 1) throw new PackError('qtlb files are little-endian; this platform is not')

const NLP_NULL = 65535, NLP_ZERO = 65534
export const NLP_MAXQ = 65533
const SE_NULL = 0xffff, SE_SIGN = 0x8000, SE_QMASK = 0x7fff, SE_MAXQ = 32766
const AF_NULL = 65535
export const AF_MAXQ = 65534, COUNT_NULL = 65535
const SNP_REF = ['', 'A', 'A', 'A', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T', 'T']
const SNP_ALT = ['', 'C', 'G', 'T', 'A', 'G', 'T', 'A', 'C', 'T', 'A', 'C', 'G']
const SNP_PAIRS = new Set(SNP_REF.slice(1).map((a, i) => `${a}\t${SNP_ALT[i + 1]}`))
const MAGIC_FILE = 0x424c5451      // "QTLB" read as a little-endian u32
const MAGIC_BLOCK = 0x30424751     // "QGB0"
const MAGIC_ZSTD = 0xfd2fb528
const MAGIC_VIDX = 0x31585651      // "QVX1"
const NO_VAR_START = 0xffffffff
// variant record flags (SPEC section 5): bit 0 alt_is_minor, bits 1-2 match code, bits 3-7 zero
const FLAG_ALT_IS_MINOR = 0x01, FLAG_MATCH_SHIFT = 1, FLAG_MATCH_MASK = 0x06, FLAG_RESERVED = 0xf8
const MATCH_NAMES = ['none', 'exact', 'position'] as const
export type MatchName = typeof MATCH_NAMES[number]
const UTF8 = new TextDecoder('utf-8', { fatal: true, ignoreBOM: true })

/** Typed-array views need a byteOffset that is a multiple of 4 (a Node Buffer from a pool may not have one). */
const aligned = (b: Uint8Array) => (b.byteOffset % 4 === 0 ? b : b.slice())
const view = (b: Uint8Array) => new DataView(b.buffer, b.byteOffset, b.byteLength)

// ---- zstd framing (SPEC section 1) --------------------------------------------------------------

/** One zstd frame: content size present, checksum flag on, no dictionary, and no bytes after the
 *  frame (found by walking the block headers). The XXH64 value itself is not verified. */
export function unframe(b: Uint8Array, expect: number | null, what: string): Uint8Array {
  const bad = (rule: string) => new PackError(`${what}: zstd frame ${rule}`)
  if (b.length < 9 || view(b).getUint32(0, true) !== MAGIC_ZSTD) throw bad('magic is missing')
  const fhd = b[4]
  const fcsFlag = fhd >> 6, single = (fhd >> 5) & 1
  if ((fhd >> 3) & 1) throw bad('has the reserved header bit set')
  if (fhd & 3) throw bad('names a dictionary')
  if (!((fhd >> 2) & 1)) throw bad('has no content checksum')
  if (fcsFlag === 0 && !single) throw bad('has no content size')
  let p = 5 + (single ? 0 : 1)
  const fcsLen = [single ? 1 : 0, 2, 4, 8][fcsFlag]
  if (p + fcsLen > b.length) throw bad('header is truncated')
  const dv = view(b)
  const size = fcsLen === 1 ? b[p] : fcsLen === 2 ? dv.getUint16(p, true) + 256 : fcsLen === 4 ? dv.getUint32(p, true)
    : dv.getUint32(p + 4, true) * 2 ** 32 + dv.getUint32(p, true)
  p += fcsLen
  for (;;) {
    if (p + 3 > b.length) throw bad('block header is truncated')
    const h = b[p] | (b[p + 1] << 8) | (b[p + 2] << 16)
    const type = (h >> 1) & 3
    if (type === 3) throw bad('uses the reserved block type')
    p += 3 + (type === 1 ? 1 : h >>> 3)
    if (h & 1) break
  }
  p += 4
  if (p !== b.length) throw bad(`ends at byte ${p} of ${b.length} (bytes after the frame, or truncated)`)
  if (expect !== null && size !== expect) throw bad(`content size ${size} != expected ${expect}`)
  const out = decompress(b, new Uint8Array(size))
  if (out.length !== size) throw bad(`decoded ${out.length} bytes, the frame says ${size}`)
  return out
}

/** An `.arrow.zst` object (annotation tables, search index): the Arrow IPC stream inside its one
 *  zstd frame, ready for DuckDB's `insertArrowFromIPCStream` or apache-arrow's `tableFromIPC`. */
export function decodeArrowObject(bytes: Uint8Array, what = 'arrow object'): Uint8Array {
  return unframe(aligned(bytes), null, what)
}

// ---- the 64-byte v1 file header (SPEC section 4) ------------------------------------------------

export interface FileHeader { kind: number; chrom: string; count: number; pageSize: number; nCis: number; seqDigest: string }

/** `qtlstore.parse_file_header`: magic, version 1, header length 64, reserved bytes 60-63 zero,
 *  chromosome non-empty ASCII without an inner NUL, seq_digest a 32-character digest. */
export function parseFileHeader(bytes: Uint8Array, what = 'file header'): FileHeader {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (bytes.length < HEADER_LEN) throw bad(`${bytes.length} bytes, need ${HEADER_LEN}`)
  const b = aligned(bytes)
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_FILE) throw bad('magic is not QTLB')
  if (b[5] !== FORMAT_VERSION) throw bad(`version ${b[5]}, this reader reads ${FORMAT_VERSION}`)
  if (dv.getUint16(6, true) !== HEADER_LEN) throw bad(`header length ${dv.getUint16(6, true)} != ${HEADER_LEN}`)
  if (b[60] || b[61] || b[62] || b[63]) throw bad('reserved bytes 60-63 are not zero')
  let end = 16
  while (end > 8 && b[end - 1] === 0) end--
  const name = b.subarray(8, end)
  if (!name.length || name.some(x => x === 0 || x >= 0x80)) throw bad('chromosome field is empty, not ASCII, or holds a NUL')
  const seqDigest = String.fromCharCode(...b.subarray(28, 60))
  if (!DIGEST.test(seqDigest)) throw bad(`seq_digest ${JSON.stringify(seqDigest)} is not a 32-character digest`)
  return { kind: b[4], chrom: String.fromCharCode(...name), count: dv.getUint32(16, true), pageSize: dv.getUint32(20, true),
    nCis: dv.getUint32(24, true), seqDigest }
}

/** A header must name the kind, chromosome and sequence (or collection) digest the pointer promised. */
export function checkFileHeader(h: FileHeader, expect: { kind: number; chrom: string; seqDigest: string }, what = 'file header'): void {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (h.kind !== expect.kind) throw bad(`kind ${h.kind}, expected ${expect.kind}`)
  if (h.chrom !== expect.chrom) throw bad(`chromosome "${h.chrom}", expected "${expect.chrom}"`)
  if (h.seqDigest !== expect.seqDigest) throw bad(`seq_digest ${h.seqDigest}, expected ${expect.seqDigest}`)
}

// ---- variant pages (SPEC section 5, kind 1) -----------------------------------------------------

interface VariantColumns {
  n: number
  position: Uint32Array
  rsNumber: Uint32Array        // 0 = no rsID
  afCode: Uint16Array          // 65535 = null
  af: Float64Array             // ALT frequency, code / 65534; NaN null
  maSamples: Uint16Array; maCount: Uint16Array   // 65535 = null
  match: Uint8Array            // 0 none, 1 exact, 2 position
  altIsMinor: Uint8Array       // flags bit 0
}
/** Every record of a run of whole pages; record i has vidx `first + i`. */
export interface VariantRange extends VariantColumns { first: number; ref: string[]; alt: string[] }
/** Exactly one phenotype's rows [varStart, varStart + n) of a variants range. */
export interface VariantRun extends VariantColumns { varStart: number; ref: string[]; alt: string[] }

/** Allele codes and heap: codes 1-12 are SNPs, code 0 takes the next `ref\talt\n` heap record. */
function decodeAlleles(codes: Uint8Array, heap: Uint8Array, bad: (rule: string) => PackError): { ref: string[]; alt: string[] } {
  const n = codes.length
  const ref: string[] = new Array(n), alt: string[] = new Array(n)
  let h = 0
  for (let i = 0; i < n; i++) {
    const c = codes[i]
    if (c > 12) throw bad(`reserved allele code ${c} at record ${i}`)
    if (c) { ref[i] = SNP_REF[c]; alt[i] = SNP_ALT[c]; continue }
    let tab = -1, j = h
    for (; j < heap.length && heap[j] !== 10; j++) {
      if (heap[j] >= 0x80) throw bad('heap is not ASCII')
      if (heap[j] === 9) { if (tab >= 0) throw bad('heap record has two tabs'); tab = j }
    }
    if (j >= heap.length || tab < 0) throw bad(`heap record for record ${i} is missing its tab or newline`)
    if (tab === h || tab + 1 === j) throw bad(`heap record for record ${i} has an empty allele`)
    ref[i] = String.fromCharCode(...heap.subarray(h, tab))
    alt[i] = String.fromCharCode(...heap.subarray(tab + 1, j))
    if (SNP_PAIRS.has(`${ref[i]}\t${alt[i]}`)) throw bad(`SNP ${ref[i]}>${alt[i]} is in the heap instead of a code`)
    h = j + 1
  }
  if (h !== heap.length) throw bad('heap holds bytes beyond its code-0 records')
  return { ref, alt }
}

/** Walk the pages of a byte range of a variants file and decode every record in them. Pages must
 *  follow each other in vidx, and positions must not decrease within the range (a range the reader
 *  asks for never crosses the cis/trans-only boundary, which restarts positions). */
export function decodeVariantRange(bytes: Uint8Array, expectLen: number | null, what = 'variants range'): VariantRange {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (expectLen !== null && bytes.length !== expectLen) throw bad(`range is ${bytes.length} bytes, expected ${expectLen}`)
  const b = aligned(bytes)
  const dv = view(b)
  const pages: { first: number; n: number; payload: Uint8Array; ref: string[]; alt: string[] }[] = []
  let off = 0, nextVidx = -1, lastPos = 0, total = 0
  while (off < b.length) {
    const where = `page at byte ${off}`
    const pbad = (rule: string) => bad(`${where}: ${rule}`)
    if (b.length - off < 12) throw pbad('header is truncated')
    const stored = dv.getUint32(off, true), first = dv.getUint32(off + 4, true), n = dv.getUint16(off + 8, true)
    const codec = b[off + 10]
    if (b[off + 11]) throw pbad('reserved byte is not zero')
    if (codec > 1) throw pbad(`codec ${codec}`)
    if (n === 0) throw pbad('no records')
    if (nextVidx >= 0 && first !== nextVidx) throw pbad(`first_vidx ${first} does not follow ${nextVidx}`)
    const end = off + 12 + stored, stop = off + Math.ceil((12 + stored) / 4) * 4
    if (stop > b.length) throw pbad(`stored length ${stored} runs past the range`)
    for (let i = end; i < stop; i++) if (b[i]) throw pbad('padding is not zero')
    const payload = aligned(codec === 0 ? b.subarray(off + 12, end) : unframe(b.subarray(off + 12, end), null, `${what} ${where}`))
    if (payload.length < 4 + 16 * n) throw pbad('payload shorter than 4 + 16n')
    const heapLen = view(payload).getUint32(0, true)
    if (payload.length !== 4 + 16 * n + heapLen) throw pbad(`payload ${payload.length} != 4 + 16n + heap_len`)
    const deltas = new Uint32Array(payload.buffer, payload.byteOffset + 4, n)
    let acc = deltas[0]
    if (acc < 1 || acc < lastPos) throw pbad('positions below 1 or decreasing')
    for (let i = 1; i < n; i++) acc += deltas[i]
    if (acc > 0xffffffff) throw pbad('positions above 2^32 - 1')
    lastPos = acc
    const flags = payload.subarray(4 + 15 * n, 4 + 16 * n)
    for (let i = 0; i < n; i++)
      if (flags[i] & FLAG_RESERVED || (flags[i] & FLAG_MATCH_MASK) >> FLAG_MATCH_SHIFT > 2) throw pbad(`reserved flag bits or match code 3 at record ${i}`)
    const { ref, alt } = decodeAlleles(payload.subarray(4 + 14 * n, 4 + 15 * n), payload.subarray(4 + 16 * n), pbad)
    pages.push({ first, n, payload, ref, alt })
    nextVidx = first + n
    total += n
    off = stop
  }
  if (!pages.length) throw bad('no pages')

  const r: VariantRange = { first: pages[0].first, n: total, position: new Uint32Array(total), ref: new Array(total), alt: new Array(total),
    rsNumber: new Uint32Array(total), afCode: new Uint16Array(total), af: new Float64Array(total),
    maSamples: new Uint16Array(total), maCount: new Uint16Array(total), match: new Uint8Array(total), altIsMinor: new Uint8Array(total) }
  let k = 0
  for (const { n, payload: p, ref, alt } of pages) {
    const at = (o: number) => p.byteOffset + o
    const deltas = new Uint32Array(p.buffer, at(4), n)
    const af = new Uint16Array(p.buffer, at(4 + 8 * n), n)
    r.rsNumber.set(new Uint32Array(p.buffer, at(4 + 4 * n), n), k)
    r.afCode.set(af, k)
    r.maSamples.set(new Uint16Array(p.buffer, at(4 + 10 * n), n), k)
    r.maCount.set(new Uint16Array(p.buffer, at(4 + 12 * n), n), k)
    let acc = 0
    for (let i = 0; i < n; i++, k++) {
      acc = i === 0 ? deltas[0] : acc + deltas[i]
      r.position[k] = acc
      r.af[k] = af[i] === AF_NULL ? NaN : af[i] / AF_MAXQ
      const f = p[4 + 15 * n + i]
      r.match[k] = (f & FLAG_MATCH_MASK) >> FLAG_MATCH_SHIFT
      r.altIsMinor[k] = f & FLAG_ALT_IS_MINOR
      r.ref[k] = ref[i]; r.alt[k] = alt[i]
    }
  }
  return r
}

/** One phenotype's run inside a decoded range, whose end positions must equal the block's
 *  pos_first and pos_last. The typed arrays are views on the range. */
export function sliceRun(range: VariantRange, varStart: number, n: number, block: { posFirst: number; posLast: number },
  what = 'variants range'): VariantRun {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const i0 = varStart - range.first
  if (!(n >= 1) || i0 < 0 || i0 + n > range.n)
    throw bad(`run ${varStart}..${varStart + n - 1} is not inside the decoded records ${range.first}..${range.first + range.n - 1}`)
  const position = range.position.subarray(i0, i0 + n)
  if (position[0] !== block.posFirst || position[n - 1] !== block.posLast)
    throw bad(`positions ${position[0]}..${position[n - 1]} != block pos_first ${block.posFirst}, pos_last ${block.posLast}`)
  return { varStart, n, position, ref: range.ref.slice(i0, i0 + n), alt: range.alt.slice(i0, i0 + n),
    rsNumber: range.rsNumber.subarray(i0, i0 + n), afCode: range.afCode.subarray(i0, i0 + n), af: range.af.subarray(i0, i0 + n),
    maSamples: range.maSamples.subarray(i0, i0 + n), maCount: range.maCount.subarray(i0, i0 + n),
    match: range.match.subarray(i0, i0 + n), altIsMinor: range.altIsMinor.subarray(i0, i0 + n) }
}

/** One record of a variants-file page, for the variant page. `A1` is ALT (the effect allele) and
 *  `A2` is REF, the labels the pages print. */
export interface VariantRecord {
  chr: string; vidx: number; position: number; A1: string; A2: string
  rsNumber: number; af: number; match: MatchName; inCis: boolean
}

/** One variants-file page on its own; `nCis` decides `inCis`. */
export function decodeVariantPage(bytes: Uint8Array, chr: string, nCis: number, what = 'variants page'): VariantRecord[] {
  const r = decodeVariantRange(bytes, null, what)
  const out: VariantRecord[] = new Array(r.n)
  for (let i = 0; i < r.n; i++) {
    const vidx = r.first + i
    out[i] = { chr, vidx, position: r.position[i], A1: r.alt[i], A2: r.ref[i], rsNumber: r.rsNumber[i],
      af: r.af[i], match: MATCH_NAMES[r.match[i]], inCis: vidx < nCis }
  }
  return out
}

// ---- variant index (kind 9) and rsID index (kind 8), SPEC section 5 -----------------------------

export interface ChromVariantIndex {
  name: string; ordinal: number
  nCis: number; nTrans: number; nPagesCis: number; nPagesTrans: number
  /** one extra entry, the file size, so page k is bytes pageOff[k] .. pageOff[k + 1] - 1 */
  pageOff: Uint32Array; pageFirstPos: Uint32Array
}
export interface VariantIndex {
  pageSize: number; rsidBlockRecords: number; rsidN: number; rsidFirst: Uint32Array
  chroms: Map<string, ChromVariantIndex>
  /** chromosome names by ordinal - 1 (the catalog's chromosome table order) */
  names: string[]
}

/** The whole `.qbx` object. `names` is the catalog pointer's chromosome table, in order: the file
 *  does not store names. `collection` is the catalog's collection digest. */
export function decodeVariantIndex(bytes: Uint8Array, names: string[], collection: string, what = 'variant index'): VariantIndex {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  const h = parseFileHeader(b, what)
  checkFileHeader(h, { kind: KIND.variantIndex, chrom: ALL, seqDigest: collection }, what)
  const p = aligned(unframe(b.subarray(HEADER_LEN), null, what))
  if (p.length < 24) throw bad('payload shorter than its 24-byte header')
  const dv = view(p)
  if (dv.getUint32(0, true) !== MAGIC_VIDX) throw bad('magic is not QVX1')
  const nChrom = dv.getUint32(4, true), pageSize = dv.getUint32(8, true)
  const rsidBlockRecords = dv.getUint32(12, true), rsidN = dv.getUint32(16, true), rsidBlocks = dv.getUint32(20, true)
  if (nChrom !== h.count || nChrom !== names.length) throw bad(`n_chrom ${nChrom}, header count ${h.count}, catalog table ${names.length}`)
  if (pageSize !== h.pageSize || pageSize < 1) throw bad(`page size ${pageSize}, header ${h.pageSize}`)
  if (rsidBlockRecords < 1) throw bad('rsid_block_records is 0')
  if (rsidBlocks !== Math.ceil(rsidN / rsidBlockRecords)) throw bad(`rsid_blocks ${rsidBlocks} does not follow from ${rsidN} records of ${rsidBlockRecords}`)
  let o = 24
  const u32s = (n: number, name: string): Uint32Array => {
    if (o + 4 * n > p.length) throw bad(`${name} runs past the payload`)
    const a = new Uint32Array(p.buffer, p.byteOffset + o, n)
    o += 4 * n
    return a
  }
  const chroms = new Map<string, ChromVariantIndex>()
  names.forEach((name, i) => {
    const [nCis, nTrans, nPagesCis, nPagesTrans] = u32s(4, `${name} counts`)
    if (nPagesCis !== Math.ceil(nCis / pageSize) || nPagesTrans !== Math.ceil(nTrans / pageSize))
      throw bad(`${name}: page counts do not follow from ${nCis} cis and ${nTrans} trans-only sites of ${pageSize}`)
    const nPages = nPagesCis + nPagesTrans
    const pageOff = u32s(nPages + 1, `${name} page_off`), pageFirstPos = u32s(nPages, `${name} page_first_position`)
    if (pageOff[0] !== HEADER_LEN) throw bad(`${name}: the first page does not start at byte ${HEADER_LEN}`)
    for (let k = 1; k <= nPages; k++) if (pageOff[k] <= pageOff[k - 1] || pageOff[k] % 4) throw bad(`${name}: page offsets do not increase by whole 4-byte units`)
    for (let k = 0; k < nPages; k++) {
      if (pageFirstPos[k] < 1) throw bad(`${name}: a page's first position is 0`)
      if (k !== 0 && k !== nPagesCis && pageFirstPos[k] < pageFirstPos[k - 1]) throw bad(`${name}: first positions decrease inside a section`)
    }
    if (chroms.has(name)) throw bad(`chromosome ${name} appears twice in the catalog table`)
    chroms.set(name, { name, ordinal: i + 1, nCis, nTrans, nPagesCis, nPagesTrans, pageOff, pageFirstPos })
  })
  const rsidFirst = u32s(rsidBlocks, 'rsid_first')
  for (let i = 1; i < rsidBlocks; i++) if (rsidFirst[i] < rsidFirst[i - 1]) throw bad('rsid_first decreases')
  if (o !== p.length) throw bad(`${p.length - o} bytes after the payload`)
  return { pageSize, rsidBlockRecords, rsidN, rsidFirst, chroms, names }
}

/** Page holding `vidx` (SPEC section 5, "Page of a vidx"). */
export function pageOfVidx(c: ChromVariantIndex, pageSize: number, vidx: number): number {
  return vidx < c.nCis ? Math.floor(vidx / pageSize) : c.nPagesCis + Math.floor((vidx - c.nCis) / pageSize)
}

/** Number of entries of the sorted array `a` in [from, to) below `x` (strict) or at most `x`. */
export function countBelow(a: Uint32Array, x: number, orEqual: boolean, from = 0, to = a.length): number {
  let lo = from, hi = to
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (a[mid] < x || (orEqual && a[mid] === x)) lo = mid + 1
    else hi = mid
  }
  return lo
}

export interface RsidRecord { rsNumber: number; vidx: number; ordinal: number }

/** The first rsID block that can hold `rs`, or -1 when no block can. SPEC section 5 says to take
 *  the last block whose first number is at or below `rs`; with repeated rs numbers the run can
 *  start one block earlier (when that block's first number equals `rs`), so this takes the last
 *  block whose first number is strictly below `rs`, else block 0. A reader scans forward from it
 *  while the next block's first number is still at or below `rs`. */
export function rsidStartBlock(rsidFirst: Uint32Array, rs: number): number {
  if (!rsidFirst.length || rsidFirst[0] > rs) return -1
  return Math.max(0, countBelow(rsidFirst, rs, false) - 1)
}

/** Scan consecutive rsID blocks from `rsidStartBlock` for the first record of `rs`'s run.
 *  `read(b)` returns block b's records; it is called for as few blocks as the run needs. */
export async function rsidFind(rsidFirst: Uint32Array, rs: number,
  read: (b: number) => RsidRecord[] | Promise<RsidRecord[]>): Promise<RsidRecord | null> {
  for (let b = rsidStartBlock(rsidFirst, rs); b >= 0 && b < rsidFirst.length && rsidFirst[b] <= rs; b++) {
    const recs = await read(b)
    const hit = recs.find(r => r.rsNumber === rs)
    if (hit) return hit
    if (recs.length && recs[recs.length - 1].rsNumber > rs) return null
  }
  return null
}
export const RSID_RECORD_LEN = 12

/** Records of one or more consecutive rsID blocks (12 bytes each, no header). Checks the zero
 *  field and the (rs_number, ordinal, vidx) order. */
export function decodeRsidRecords(bytes: Uint8Array, what = 'rsID block'): RsidRecord[] {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (bytes.length % RSID_RECORD_LEN) throw bad(`${bytes.length} bytes is not a whole number of 12-byte records`)
  const dv = view(bytes)
  const out: RsidRecord[] = []
  for (let o = 0; o < bytes.length; o += RSID_RECORD_LEN) {
    const r = { rsNumber: dv.getUint32(o, true), vidx: dv.getUint32(o + 4, true), ordinal: dv.getUint16(o + 8, true) }
    if (dv.getUint16(o + 10, true)) throw bad(`record ${o / RSID_RECORD_LEN}: the zero field is not zero`)
    if (r.rsNumber === 0 || r.ordinal === 0) throw bad(`record ${o / RSID_RECORD_LEN}: rs_number or ordinal is 0`)
    const q = out[out.length - 1]
    if (q && (r.rsNumber - q.rsNumber || r.ordinal - q.ordinal || r.vidx - q.vidx) <= 0)
      throw bad(`record ${o / RSID_RECORD_LEN}: not strictly ascending by (rs_number, ordinal, vidx)`)
    out.push(r)
  }
  return out
}

// ---- result blocks (SPEC section 8, kind 2) -------------------------------------------------------

/** What the search index says about a block: its length, and the run (null when it has no rows). */
export interface BlockExpect { blk_len: number | null; n_var?: number | null; var_start?: number | null }

/** Details JSON `v: 1` (SPEC section 8): study fields only, never annotation. */
export interface Details {
  v: 1; phenotype_type: string; phenotype_id: string; phenotype_object_id: string; gene_id: string | null
  has_nominal: boolean; n_nominal: number; n_credible_sets: number
  extra: Record<string, unknown>
  group: null | {
    lead_phenotype_id: string; n_variants: number | null; p_perm: number | null; p_beta: number | null; significant: boolean
    lead: { chr: string; pos: number; ref: string; alt: string }
  }
}

export interface ResultBlock {
  nRows: number; varStart: number | null; posFirst: number; posLast: number
  nlpMax: number; lseMin: number; lseMax: number
  details: Details
  nlpCode: Uint16Array; seCode: Uint16Array
  nlp: Float64Array            // -log10 p; NaN null, Infinity for p = 0
  pval: Float64Array           // NaN null, 0 for p = 0
  se: Float64Array             // exp(lse_min + q * step); NaN null
  negative: Uint8Array         // the SE code's sign bit (0 on null rows)
  /** sign * se * t(p, dof); NaN where p is null or 0, the SE code is null, or the results set has
   *  no dof (SPEC section 8: `dof: null` means no slope, and no stand-in dof is ever used) */
  slope: Float64Array
  /** credible-set records, one per membership, sorted by (row, cs_id) */
  csRow: Uint32Array; csPip: Float32Array; csId: Uint8Array
}

/** A block's details frame, one zstd frame of JSON with `"v": 1`. */
export function decodeDetails(frame: Uint8Array, expectLen: number, what = 'block'): Details {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const raw = unframe(frame, expectLen, `${what} details`)
  if (raw[0] === 0xef && raw[1] === 0xbb && raw[2] === 0xbf) throw bad('details JSON starts with a byte-order mark')
  let parsed: unknown
  try { parsed = JSON.parse(UTF8.decode(raw)) } catch (e) { throw bad(`details JSON: ${(e as Error).message}`) }   // JSON.parse rejects NaN and Infinity
  if (parsed === null || typeof parsed !== 'object' || Array.isArray(parsed) || (parsed as { v?: unknown }).v !== 1)
    throw bad('details are not a JSON object with "v": 1')
  return parsed as Details
}

/** `packfmt_v1.decode_gene_block` with `details_version` 1, then the SPEC section 11 values with the
 *  results set's `dof`; `dof` null leaves every slope NaN. */
export function decodeResultBlock(bytes: Uint8Array, dof: number | null, expect: BlockExpect, what = 'result block'): ResultBlock {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  if (b.length < 64) throw bad('shorter than the 64-byte block header')
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_BLOCK) throw bad('magic is not QGB0')
  const blkLen = dv.getUint32(4, true), n = dv.getUint32(8, true), varStart = dv.getUint32(12, true)
  const anchor = dv.getInt32(16, true), posFirst = dv.getUint32(20, true), posLast = dv.getUint32(24, true)
  const nCs = dv.getUint32(28, true)
  const nlpMax = dv.getFloat64(32, true), lseMin = dv.getFloat64(40, true), lseMax = dv.getFloat64(48, true)
  const dz = dv.getUint32(56, true), dl = dv.getUint32(60, true)
  if (blkLen !== b.length || (expect.blk_len !== null && blkLen !== expect.blk_len))
    throw bad(`length field ${blkLen} != range ${b.length} / expected blk_len ${expect.blk_len}`)
  const body = 64 + 4 * n + 12 * nCs + dz
  if (Math.ceil(body / 4) * 4 !== blkLen) throw bad(`64 + 4 n_rows + 12 n_cs + details_zlen = ${body}, padded to 4, != block length ${blkLen}`)
  for (let i = body; i < blkLen; i++) if (b[i]) throw bad('padding is not zero')
  if (!dz) throw bad('a v1 block needs a details frame (details_zlen above 0)')
  if (anchor !== 0) throw bad(`anchor ${anchor}: v1 blocks store 0 (the TSS comes from the annotation)`)
  if (expect.n_var !== undefined && n !== (expect.n_var ?? 0)) throw bad(`n_rows ${n} != search index n_var ${expect.n_var}`)
  if (n === 0) {
    if (varStart !== NO_VAR_START || posFirst || posLast || nCs || nlpMax !== 0 || lseMin !== 0 || lseMax !== 0)
      throw bad('a block with no rows needs var_start 0xFFFFFFFF and zero positions, credible sets, and scales')
    if (expect.var_start != null) throw bad(`n_rows 0 but search index var_start ${expect.var_start}`)
  } else {
    if (expect.var_start !== undefined && varStart !== expect.var_start) throw bad(`var_start ${varStart} != search index var_start ${expect.var_start}`)
    if (varStart + n - 1 >= NO_VAR_START) throw bad('var_start + n_rows - 1 reaches 0xFFFFFFFF')
    if (!(posFirst >= 1 && posFirst <= posLast)) throw bad(`pos_first ${posFirst}, pos_last ${posLast} out of order`)
  }
  if (!(Number.isFinite(nlpMax) && nlpMax >= 0 && Number.isFinite(lseMin) && Number.isFinite(lseMax) && lseMin <= lseMax))
    throw bad(`scales nlp_max ${nlpMax}, lse_min ${lseMin}, lse_max ${lseMax} are not valid`)

  const codes = new Uint16Array(b.buffer, b.byteOffset + 64, 2 * n)
  const nlpCode = new Uint16Array(n), seCode = new Uint16Array(n)
  let maxFinite = -1, lqMin = Infinity, lqMax = -1
  for (let i = 0; i < n; i++) {
    const q = (nlpCode[i] = codes[2 * i]), s = (seCode[i] = codes[2 * i + 1])
    if (q <= NLP_MAXQ && q > maxFinite) maxFinite = q
    if (s !== SE_NULL) {
      const lq = s & SE_QMASK
      if (lq === SE_QMASK) throw bad('SE code 0x7FFF is invalid')
      if (lq < lqMin) lqMin = lq
      if (lq > lqMax) lqMax = lq
    }
  }
  if (nlpMax > 0 ? maxFinite !== NLP_MAXQ : maxFinite > 0) throw bad('-log10 p codes break the scale rule')
  if (lqMax < 0 ? lseMin !== 0 || lseMax !== 0 : lseMin === lseMax ? lqMax !== 0 : lqMin !== 0 || lqMax !== SE_MAXQ)
    throw bad('log(SE) codes break the scale rule')

  const csRow = new Uint32Array(nCs), csPip = new Float32Array(nCs), csId = new Uint8Array(nCs)
  for (let k = 0, o = 64 + 4 * n; k < nCs; k++, o += 12) {
    csRow[k] = dv.getUint32(o, true); csPip[k] = dv.getFloat32(o + 4, true); csId[k] = b[o + 8]
    if (b[o + 9] || b[o + 10] || b[o + 11]) throw bad('credible-set padding is not zero')
    if (csRow[k] >= n) throw bad(`credible-set row ${csRow[k]} is not below n_rows ${n}`)
    if (k > 0 && (csRow[k] < csRow[k - 1] || (csRow[k] === csRow[k - 1] && csId[k] <= csId[k - 1])))
      throw bad('credible-set (row, cs_id) pairs are not strictly ascending')
    if (!(csPip[k] >= 0 && csPip[k] <= 1) || csId[k] > 127) throw bad('credible-set pip outside [0, 1] or cs_id above 127')
  }

  const details = decodeDetails(b.subarray(64 + 4 * n + 12 * nCs, body), dl, what)

  const nlp = new Float64Array(n), pval = new Float64Array(n), se = new Float64Array(n), slope = new Float64Array(n)
  const negative = new Uint8Array(n)
  const nlpStep = nlpMax / NLP_MAXQ, seStep = (lseMax - lseMin) / SE_MAXQ
  const tByCode = new Map<number, number>()
  for (let i = 0; i < n; i++) {
    const q = nlpCode[i], s = seCode[i]
    if (q === NLP_NULL) { nlp[i] = NaN; pval[i] = NaN }
    else if (q === NLP_ZERO) { nlp[i] = Infinity; pval[i] = 0 }
    else { nlp[i] = q * nlpStep; pval[i] = Math.pow(10, -nlp[i]) }
    if (s === SE_NULL) { se[i] = NaN; slope[i] = NaN; continue }
    se[i] = Math.exp(lseMin + (s & SE_QMASK) * seStep)
    negative[i] = (s & SE_SIGN) ? 1 : 0
    if (dof === null || q >= NLP_ZERO) { slope[i] = NaN; continue }
    let t = tByCode.get(q)
    if (t === undefined) { t = tFromNlp(nlp[i], dof); tByCode.set(q, t) }
    slope[i] = (negative[i] ? -1 : 1) * se[i] * t
  }
  return { nRows: n, varStart: n ? varStart : null, posFirst, posLast, nlpMax, lseMin, lseMax,
    details, nlpCode, seCode, nlp, pval, se, negative, slope, csRow, csPip, csId }
}

/** Every credible-set membership of a block, for the credible-set table (a row in two sets appears twice). */
export function csMembers(block: ResultBlock): { row: number; pip: number; csId: number }[] {
  return Array.from(block.csRow, (row, k) => ({ row, pip: block.csPip[k], csId: block.csId[k] }))
}

/** A block's rows joined to its variant run, as the locus table's typed columns with in-band
 *  nulls (NaN for floats, -1 for counts and cs_id, 0 for rs_number). `A1` is ALT, `A2` REF, and
 *  tss_distance is position minus the annotation's TSS. A row in two credible sets takes its
 *  higher-PIP record, and the lower cs_id on a tie. */
export interface ReaderColumns {
  n: number
  position: Int32Array; a1: string[]; a2: string[]; rsNumber: Uint32Array; tssDistance: Int32Array
  af: Float32Array; maSamples: Int16Array; maCount: Int16Array
  pval: Float64Array; slope: Float32Array; se: Float32Array; pip: Float32Array; csId: Int8Array
}

export function readerColumns(block: ResultBlock, run: VariantRun, tss: number): ReaderColumns {
  const n = block.nRows
  if (run.n !== n || run.varStart !== block.varStart) throw new PackError(`reader: variant run ${run.varStart}+${run.n} does not match block ${block.varStart}+${n}`)
  const c: ReaderColumns = { n, position: new Int32Array(n), a1: run.alt, a2: run.ref, rsNumber: run.rsNumber, tssDistance: new Int32Array(n),
    af: new Float32Array(n), maSamples: new Int16Array(n), maCount: new Int16Array(n), pval: block.pval,
    slope: Float32Array.from(block.slope), se: Float32Array.from(block.se), pip: new Float32Array(n).fill(NaN), csId: new Int8Array(n).fill(-1) }
  for (let i = 0; i < n; i++) {
    const pos = run.position[i], d = pos - tss
    if (pos > 0x7fffffff || d > 0x7fffffff || d < -0x80000000) throw new PackError('reader: position or tss_distance does not fit INTEGER')
    c.position[i] = pos; c.tssDistance[i] = d
    c.af[i] = run.af[i]
    const ms = run.maSamples[i], mc = run.maCount[i]
    if ((ms !== COUNT_NULL && ms > 32767) || (mc !== COUNT_NULL && mc > 32767)) throw new PackError('reader: a count does not fit SMALLINT')
    c.maSamples[i] = ms === COUNT_NULL ? -1 : ms
    c.maCount[i] = mc === COUNT_NULL ? -1 : mc
  }
  for (let k = 0; k < block.csRow.length; k++) {
    const row = block.csRow[k]
    if (c.csId[row] < 0 || block.csPip[k] > c.pip[row]) { c.pip[row] = block.csPip[k]; c.csId[row] = block.csId[k] }
  }
  return c
}

/** The block header fields a span walk needs (the cis scan reads block headers only). */
export function blockHeader(span: Uint8Array, at: number, what = 'span'): { blkLen: number; nRows: number; varStart: number } {
  if (at + 64 > span.length) throw new PackError(`${what}: a block header at byte ${at} runs past the span`)
  const dv = new DataView(span.buffer, span.byteOffset + at, 64)
  if (dv.getUint32(0, true) !== MAGIC_BLOCK) throw new PackError(`${what}: the block at byte ${at} does not start with QGB0`)
  return { blkLen: dv.getUint32(4, true), nRows: dv.getUint32(8, true), varStart: dv.getUint32(12, true) }
}

/** One block read through its header, the one raw pair of `row`, and that row's best credible-set
 *  record; nothing is decompressed. `dof` null gives a NaN slope. */
export function scanBlockRow(block: Uint8Array, row: number, dof: number | null, what = 'span block'):
  { pval: number; se: number; slope: number; pip: number | null; csId: number | null } {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(block)
  if (b.length < 64) throw bad('shorter than the 64-byte block header')
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_BLOCK) throw bad('magic is not QGB0')
  const blkLen = dv.getUint32(4, true), n = dv.getUint32(8, true), nCs = dv.getUint32(28, true)
  const nlpMax = dv.getFloat64(32, true), lseMin = dv.getFloat64(40, true), lseMax = dv.getFloat64(48, true)
  if (blkLen !== b.length) throw bad(`length field ${blkLen} != slice ${b.length}`)
  if (!(row >= 0 && row < n)) throw bad(`row ${row} is not below n_rows ${n}`)
  if (!(Number.isFinite(nlpMax) && nlpMax >= 0 && Number.isFinite(lseMin) && Number.isFinite(lseMax) && lseMin <= lseMax))
    throw bad(`scales nlp_max ${nlpMax}, lse_min ${lseMin}, lse_max ${lseMax} are not valid`)
  if (64 + 4 * n + 12 * nCs > blkLen) throw bad('rows and credible sets run past the block')

  const q = dv.getUint16(64 + 4 * row, true), s = dv.getUint16(64 + 4 * row + 2, true)
  const nlp = q === NLP_NULL ? NaN : q === NLP_ZERO ? Infinity : q * (nlpMax / NLP_MAXQ)
  const pval = q === NLP_NULL ? NaN : q === NLP_ZERO ? 0 : Math.pow(10, -nlp)
  let se = NaN, slope = NaN
  if (s !== SE_NULL) {
    const lq = s & SE_QMASK
    if (lq === SE_QMASK) throw bad('SE code 0x7FFF is invalid')
    se = Math.exp(lseMin + lq * ((lseMax - lseMin) / SE_MAXQ))
    if (dof !== null && q < NLP_ZERO) slope = ((s & SE_SIGN) ? -1 : 1) * se * tFromNlp(nlp, dof)
  }
  let pip: number | null = null, csId: number | null = null
  for (let k = 0, o = 64 + 4 * n; k < nCs; k++, o += 12) {
    const r = dv.getUint32(o, true)
    if (r < row) continue
    if (r > row) break
    const v = dv.getFloat32(o + 4, true), id = b[o + 8]
    if (!(v >= 0 && v <= 1) || id > 127) throw bad('credible-set pip outside [0, 1] or cs_id above 127')
    if (pip === null || v > pip) { pip = v; csId = id }
  }
  return { pval, se, slope, pip, csId }
}

// ---- hits (SPEC section 8, kind 7) -----------------------------------------------------------------

export const HIT_LEAD = 0, HIT_CS = 1, HIT_TRANS = 2
export const HIT_RECORD_LEN = 20

/** A hits file's header and frame table: frame g is bytes frameOff[g] .. frameOff[g + 1] - 1 and
 *  holds vidx g * frameVariants .. (g + 1) * frameVariants - 1; equal offsets mean no records. */
export interface HitsTable { count: number; frameVariants: number; nVariants: number; frameOff: Uint32Array }

/** Bytes needed to read a hits file's header and table: the header, then 4 (n_frames + 1). */
export const hitsTableLen = (nVariants: number, frameVariants: number) => HEADER_LEN + 4 * (Math.ceil(nVariants / frameVariants) + 1)

/** `packfmt_v1.hits_frame_table` with the header checks: `bytes` holds at least the header and the
 *  table (`hitsTableLen`); `expect.nVariants` is the chromosome's variant count in the catalog. */
export function decodeHitsTable(bytes: Uint8Array, expect: { chrom: string; seqDigest: string; nVariants: number }, what = 'hits file'): HitsTable {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  const h = parseFileHeader(b, what)
  checkFileHeader(h, { kind: KIND.hits, chrom: expect.chrom, seqDigest: expect.seqDigest }, what)
  if (h.pageSize < 1) throw bad('frame size (page size) is 0')
  if (h.nCis !== expect.nVariants) throw bad(`variant count ${h.nCis}, the catalog says ${expect.nVariants}`)
  const nFrames = Math.ceil(h.nCis / h.pageSize)
  if (b.length < hitsTableLen(h.nCis, h.pageSize)) throw bad('the frame table is truncated')
  const frameOff = new Uint32Array(b.buffer, b.byteOffset + HEADER_LEN, nFrames + 1).slice()
  if (frameOff[0] !== HEADER_LEN + 4 * (nFrames + 1)) throw bad('the frame table does not start right after itself')
  for (let g = 1; g <= nFrames; g++) if (frameOff[g] < frameOff[g - 1]) throw bad('frame offsets decrease')
  return { count: h.count, frameVariants: h.pageSize, nVariants: h.nCis, frameOff }
}

/** One frame's records as parallel columns, sorted by (vidx, kind, ord, cs_id). Kind 0 is a group's
 *  lead (`value` p_perm, NaN when null; flags bit 0 significant), kind 1 a credible-set member
 *  (`value` PIP), kind 2 a trans association (`value` -log10 p, +Infinity for p = 0; `beta` ALT). */
export interface HitsFrame {
  count: number
  vidx: Uint32Array; ord: Uint32Array; value: Float32Array; beta: Float32Array; kind: Uint8Array; csId: Uint8Array; flags: Uint8Array
}

/** `packfmt_v1.decode_hits_frame`: a zero-length frame has no records. */
export function decodeHitsFrame(bytes: Uint8Array, firstVidx: number, frameVariants: number, what = 'hits frame'): HitsFrame {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const empty = (n: number): HitsFrame => ({ count: n, vidx: new Uint32Array(n), ord: new Uint32Array(n), value: new Float32Array(n),
    beta: new Float32Array(n), kind: new Uint8Array(n), csId: new Uint8Array(n), flags: new Uint8Array(n) })
  if (!bytes.length) return empty(0)
  const p = aligned(unframe(aligned(bytes), null, what))
  if (!p.length || p.length % HIT_RECORD_LEN) throw bad(`${p.length} bytes is not a positive multiple of 20`)
  const n = p.length / HIT_RECORD_LEN, dv = view(p), out = empty(n)
  for (let i = 0, o = 0; i < n; i++, o += HIT_RECORD_LEN) {
    const v = dv.getUint32(o, true), ord = dv.getUint32(o + 4, true), beta = dv.getFloat32(o + 12, true)
    const k = p[o + 16], cs = p[o + 17], f = p[o + 18]
    if (v < firstVidx || v >= firstVidx + frameVariants) throw bad(`record ${i}: vidx ${v} outside the frame`)
    if (k > HIT_TRANS || p[o + 19] || (f && k !== HIT_LEAD) || f & 0xfe) throw bad(`record ${i}: bad kind, flags or pad`)
    if (Number.isFinite(beta) && k !== HIT_TRANS) throw bad(`record ${i}: beta set on a lead or credible-set record`)
    if (i > 0 && (v - out.vidx[i - 1] || k - out.kind[i - 1] || ord - out.ord[i - 1] || cs - out.csId[i - 1]) < 0)
      throw bad(`record ${i}: not sorted by (vidx, kind, ord, cs_id)`)
    out.vidx[i] = v; out.ord[i] = ord; out.value[i] = dv.getFloat32(o + 8, true); out.beta[i] = beta
    out.kind[i] = k; out.csId[i] = cs; out.flags[i] = f
  }
  return out
}

/** Record indices [lo, hi) of one vidx in a frame. */
export function hitsOf(f: HitsFrame, vidx: number): [number, number] {
  const lo = countBelow(f.vidx, vidx, false)
  return [lo, countBelow(f.vidx, vidx, true, lo)]
}

// ---- trans frames (SPEC section 9, kind 6) ------------------------------------------------------

const MAGIC_TRANS = 0x32545451     // "QTT2"
export const BETA_MAXQ = 32767

/** One phenotype's trans rows (`packfmt_v1.decode_trans_frame`), sorted by (ordinal, pos, ref, alt).
 *  No vidx: links use the rsID or chr:pos. */
export interface TransFrame {
  n: number; nlpMax: number; betaMax: number
  ordinal: Uint8Array; position: Uint32Array; rsNumber: Uint32Array
  afCode: Uint16Array; af: Float64Array            // ALT, NaN null
  nlpCode: Uint16Array; nlp: Float64Array; pval: Float64Array
  betaCode: Int16Array; beta: Float64Array
  ref: string[]; alt: string[]
}

export function decodeTransFrame(bytes: Uint8Array, what = 'trans frame'): TransFrame {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const p = aligned(unframe(aligned(bytes), null, what))
  if (p.length < 32) throw bad('shorter than its 32-byte header')
  const dv = view(p)
  if (dv.getUint32(0, true) !== MAGIC_TRANS) throw bad('magic is not QTT2')
  const n = dv.getUint32(4, true), nlpMax = dv.getFloat64(8, true), betaMax = dv.getFloat64(16, true)
  const heapLen = dv.getUint32(24, true)
  if (dv.getUint32(28, true) || n < 1) throw bad('reserved field is not zero, or no rows')
  if (p.length !== 32 + 16 * n + heapLen) throw bad(`payload ${p.length} != 32 + 16n + heap_len`)
  if (!(Number.isFinite(nlpMax) && nlpMax >= 0 && Number.isFinite(betaMax) && betaMax >= 0)) throw bad('bad scales')
  // every column starts at a multiple of 2 after the 32-byte header, and the u32 ones at a multiple of 4
  const at = (o: number) => p.byteOffset + o
  const delta = new Uint32Array(p.buffer, at(32), n), rsNumber = new Uint32Array(p.buffer, at(32 + 4 * n), n)
  const afCode = new Uint16Array(p.buffer, at(32 + 8 * n), n), nlpCode = new Uint16Array(p.buffer, at(32 + 10 * n), n)
  const betaCode = new Int16Array(p.buffer, at(32 + 12 * n), n)
  const ordinal = p.subarray(32 + 14 * n, 32 + 15 * n), codes = p.subarray(32 + 15 * n, 32 + 16 * n)
  const position = new Uint32Array(n), af = new Float64Array(n), nlp = new Float64Array(n), pval = new Float64Array(n), beta = new Float64Array(n)
  const nlpStep = nlpMax / NLP_MAXQ, betaStep = betaMax / BETA_MAXQ
  let pos = 0
  for (let i = 0; i < n; i++) {
    if (ordinal[i] < 1 || (i > 0 && ordinal[i] < ordinal[i - 1])) throw bad(`row ${i}: chromosome ordinals not >= 1 and non-decreasing`)
    if (codes[i] > 12 || betaCode[i] === -32768 || nlpCode[i] === NLP_NULL) throw bad(`row ${i}: reserved allele, beta or -log10 p code`)
    pos = i === 0 || ordinal[i] !== ordinal[i - 1] ? delta[i] : pos + delta[i]
    if (pos > 0xffffffff) throw bad(`row ${i}: position above 2^32 - 1`)
    position[i] = pos
    af[i] = afCode[i] === AF_NULL ? NaN : afCode[i] / AF_MAXQ
    if (nlpCode[i] === NLP_ZERO) { nlp[i] = Infinity; pval[i] = 0 }
    else { nlp[i] = nlpCode[i] * nlpStep; pval[i] = Math.pow(10, -nlp[i]) }
    beta[i] = betaCode[i] * betaStep
  }
  const { ref, alt } = decodeAlleles(codes, p.subarray(32 + 16 * n), bad)
  return { n, nlpMax, betaMax, ordinal, position, rsNumber, afCode, af, nlpCode, nlp, pval, betaCode, beta, ref, alt }
}

/** A trans row's SE and r2 from its -log10 p and beta: se = |beta| / t(p, dof), r2 = t^2 / (t^2 + dof);
 *  both NaN without a dof. */
export function transSe(nlp: number, beta: number, dof: number | null, tCache?: Map<number, number>): { se: number; r2: number } {
  if (dof === null || Number.isNaN(nlp)) return { se: NaN, r2: NaN }
  let t = tCache?.get(nlp)
  if (t === undefined) { t = tFromNlp(nlp, dof); tCache?.set(nlp, t) }
  return { se: Math.abs(beta) / t, r2: t === Infinity ? 1 : (t * t) / (t * t + dof) }
}

// ---- GWAS (SPEC section 10, kinds 4 and 5) ---------------------------------------------------------

const GWAS_SCALE = 10000
// the double nearest each power of ten, as the reader divides by it (p_exp is -128 to -3)
const POW10 = Array.from({ length: 129 }, (_, k) => Number(`1e${k}`))

/** The GWAS index (`gwas.decode_index`): the n table, and per chromosome each block's first
 *  position and end offset in that chromosome's `.qbg`. */
export interface GwasIndex {
  blockRows: number; nValues: Uint32Array
  chroms: Map<string, { firstPosition: Uint32Array; endOffset: Uint32Array }>
}

export function decodeGwasIndex(bytes: Uint8Array, collection: string, what = 'GWAS index'): GwasIndex {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  const h = parseFileHeader(b, what)
  checkFileHeader(h, { kind: KIND.gwasIndex, chrom: ALL, seqDigest: collection }, what)
  const p = aligned(unframe(b.subarray(HEADER_LEN), null, what))
  let off = 0
  const u32s = (k: number) => {
    if (off + 4 * k > p.length) throw bad('payload ends early')
    const a = new Uint32Array(p.slice(off, off + 4 * k).buffer)
    off += 4 * k
    return a
  }
  const nValues = u32s(u32s(1)[0])
  if (nValues.length < 1 || nValues.length > 255 || nValues.some((x, i) => i > 0 && x <= nValues[i - 1]))
    throw bad('n values must be 1 to 255 ascending distinct values')
  const nChroms = u32s(1)[0]
  if (nChroms !== h.count) throw bad(`${nChroms} chromosomes, header count ${h.count}`)
  const chroms = new Map<string, { firstPosition: Uint32Array; endOffset: Uint32Array }>()
  for (let c = 0; c < nChroms; c++) {
    if (off + 8 > p.length) throw bad('payload ends early')
    let end = 8
    while (end > 0 && p[off + end - 1] === 0) end--
    const raw = p.subarray(off, off + end)
    off += 8
    if (!raw.length || raw.some(x => x === 0 || x >= 0x80)) throw bad('bad chromosome name')
    const name = String.fromCharCode(...raw)
    if (chroms.has(name)) throw bad(`chromosome ${name} appears twice`)
    const nb = u32s(1)[0]
    const firstPosition = u32s(nb), endOffset = u32s(nb)
    if (nb < 1 || firstPosition[0] < 1 || endOffset[0] <= HEADER_LEN ||
        firstPosition.some((x, i) => i > 0 && x < firstPosition[i - 1]) || endOffset.some((x, i) => i > 0 && x <= endOffset[i - 1]))
      throw bad(`${name}: blocks must be non-empty, first positions non-decreasing, end offsets increasing past byte 64`)
    chroms.set(name, { firstPosition, endOffset })
  }
  if (off !== p.length) throw bad(`${p.length - off} bytes after the last chromosome`)
  return { blockRows: h.pageSize, nValues, chroms }
}

/** One range request of a GWAS file: bytes [off, off + len), holding blocks firstBlock..lastBlock. */
export interface GwasRange { off: number; len: number; firstBlock: number; lastBlock: number }

/** The window rule (`packfmt_v1.gwas_window`), or null when no block starts at or before `hi`. */
export function gwasRange(index: GwasIndex, chr: string, lo: number, hi: number): GwasRange | null {
  if (lo > hi) throw new PackError(`GWAS window ${chr}:${lo}-${hi}: lo is above hi`)
  const c = index.chroms.get(chr)
  if (!c) return null
  const end = countBelow(c.firstPosition, hi, true) - 1
  if (end < 0) return null
  const start = Math.max(countBelow(c.firstPosition, lo, false) - 1, 0)
  const off = start === 0 ? HEADER_LEN : c.endOffset[start - 1]
  return { off, len: c.endOffset[end] - off, firstBlock: start, lastBlock: end }
}

/** GWAS rows with lo <= position <= hi, in file order. `ref`/`alt` oriented to the reference;
 *  `beta` and `af` describe ALT. */
export interface GwasColumns {
  rows: number
  position: Int32Array; ref: string[]; alt: string[]
  beta: Float64Array; se: Float64Array; af: Float64Array; p: Float64Array
  rsNumber: Uint32Array; n: Int32Array
}

/** Decode the blocks of one GWAS range (`packfmt_v1.decode_gwas_block` per block) and keep the rows
 *  in the window. */
export function decodeGwasRange(bytes: Uint8Array, index: GwasIndex, chr: string, range: GwasRange, lo: number, hi: number,
  what = `GWAS ${chr}`): GwasColumns {
  const bad = (rule: string) => new PackError(`${what} bytes ${range.off}+${range.len}: ${rule}`)
  const c = index.chroms.get(chr)
  if (!c) throw bad('the index has no such chromosome')
  const { firstBlock: s, lastBlock: e } = range
  if (!(s >= 0 && s <= e && e < c.firstPosition.length)) throw bad(`blocks ${s}..${e} outside the index`)
  const start = s === 0 ? HEADER_LEN : c.endOffset[s - 1]
  if (range.off !== start || range.len !== c.endOffset[e] - start || bytes.length !== range.len) throw bad('the range does not match its blocks')
  const nv = index.nValues
  const blocks: { n: number; p: Uint8Array; pos: Uint32Array; a: { ref: string[]; alt: string[] } }[] = []
  let keep = 0
  for (let k = s; k <= e; k++) {
    const kbad = (rule: string) => bad(`block ${k}: ${rule}`)
    const from = (k === 0 ? HEADER_LEN : c.endOffset[k - 1]) - range.off, to = c.endOffset[k] - range.off
    const p = aligned(unframe(bytes.subarray(from, to), null, `${what} block ${k}`))
    if (p.length < 8) throw kbad('payload is shorter than its 8-byte header')
    const dv = view(p)
    const n = dv.getUint32(0, true), heapLen = dv.getUint32(4, true)
    if (n < 1) throw kbad('no rows')
    if (p.length !== 8 + 21 * n + heapLen) throw kbad(`payload ${p.length} != 8 + 21n + heap_len`)
    const at = (o: number) => p.byteOffset + o
    const deltas = new Uint32Array(p.buffer, at(8), n)
    const pos = new Uint32Array(n)
    let acc = 0
    for (let i = 0; i < n; i++) {
      acc = i === 0 ? deltas[0] : acc + deltas[i]
      pos[i] = acc
      if (acc >= lo && acc <= hi) keep++
    }
    if (pos[0] < 1 || acc > 0x7fffffff) throw kbad('positions below 1 or beyond INTEGER')
    if (pos[0] !== c.firstPosition[k]) throw kbad(`first position ${pos[0]} != index ${c.firstPosition[k]}`)
    const af = new Uint16Array(p.buffer, at(8 + 14 * n), n), mant = new Uint16Array(p.buffer, at(8 + 16 * n), n)
    const exp = new Int8Array(p.buffer, at(8 + 18 * n), n), nCode = p.subarray(8 + 19 * n, 8 + 20 * n)
    for (let i = 0; i < n; i++) {
      if (af[i] > GWAS_SCALE) throw kbad(`af code ${af[i]} above 10000 at row ${i}`)
      if (mant[i] < 1000 || mant[i] > 9999) throw kbad(`p mantissa ${mant[i]} outside 1000..9999 at row ${i}`)
      if (exp[i] > -3 || (exp[i] === -3 && mant[i] !== 1000)) throw kbad(`p above 1 at row ${i}`)
      if (nCode[i] >= nv.length) throw kbad(`n code ${nCode[i]} beyond the n table at row ${i}`)
    }
    const a = decodeAlleles(p.subarray(8 + 20 * n, 8 + 21 * n), p.subarray(8 + 21 * n), kbad)
    blocks.push({ n, p, pos, a })
  }
  const out: GwasColumns = { rows: keep, position: new Int32Array(keep), ref: new Array(keep), alt: new Array(keep),
    beta: new Float64Array(keep), se: new Float64Array(keep), af: new Float64Array(keep), p: new Float64Array(keep),
    rsNumber: new Uint32Array(keep), n: new Int32Array(keep) }
  let j = 0
  for (const { n, p, pos, a } of blocks) {
    const at = (o: number) => p.byteOffset + o
    const beta = new Int32Array(p.buffer, at(8 + 4 * n), n), rs = new Uint32Array(p.buffer, at(8 + 8 * n), n)
    const se = new Uint16Array(p.buffer, at(8 + 12 * n), n), af = new Uint16Array(p.buffer, at(8 + 14 * n), n)
    const mant = new Uint16Array(p.buffer, at(8 + 16 * n), n), exp = new Int8Array(p.buffer, at(8 + 18 * n), n)
    for (let i = 0; i < n; i++) {
      if (pos[i] < lo || pos[i] > hi) continue
      out.position[j] = pos[i]; out.ref[j] = a.ref[i]; out.alt[j] = a.alt[i]
      out.beta[j] = beta[i] / GWAS_SCALE; out.se[j] = se[i] / GWAS_SCALE; out.af[j] = af[i] / GWAS_SCALE
      out.p[j] = mant[i] / POW10[-exp[i]]
      out.rsNumber[j] = rs[i]
      out.n[j] = nv[p[8 + 19 * n + i]]
      j++
    }
  }
  return out
}

// ---- GWAS bin summary (SPEC section 10, "Bin summary") -------------------------------------------

/** One bin of the GWAS bin summary: v0's gwas_dcm_bins.json row. */
export interface GwasBin {
  chr: string; bin_start: number; bin_end: number; min_p: number; lead_position: number; lead_rsid: string | null
  lead_beta: number; lead_ea: string; n_gws: number; n_variants: number
}
/** `gwas.BINS_SCHEMA`: column names and Arrow types, in order. */
export const GWAS_BINS_SCHEMA: [keyof GwasBin, string][] = [['chr', 'Utf8'], ['bin_start', 'Uint32'], ['bin_end', 'Uint32'],
  ['min_p', 'Float64'], ['lead_position', 'Uint32'], ['lead_rsid', 'Utf8'], ['lead_beta', 'Float64'], ['lead_ea', 'Utf8'],
  ['n_gws', 'Uint32'], ['n_variants', 'Uint32']]

/** The bin summary object (one zstd frame around an Arrow IPC stream) as rows, after checking its
 *  schema, the pointer's bin count and width, and that every bin holds its lead position. */
export function decodeGwasBins(bytes: Uint8Array, expect: { n_bins: number; bin_bp: number }, what = 'GWAS bins'): GwasBin[] {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const t = tableFromIPC(unframe(aligned(bytes), null, what))
  const got = t.schema.fields.map(f => `${f.name}:${String(f.type)}`).join(',')
  if (got !== GWAS_BINS_SCHEMA.map(([n, ty]) => `${n}:${ty}`).join(',')) throw bad(`schema ${got}`)
  if (t.numRows !== expect.n_bins) throw bad(`${t.numRows} bins, the pointer says ${expect.n_bins}`)
  const cols = GWAS_BINS_SCHEMA.map(([n]) => [n, t.getChild(n)!] as const)
  const out: GwasBin[] = []
  for (let i = 0; i < t.numRows; i++) {
    const b = Object.fromEntries(cols.map(([n, c]) => [n, c.get(i)])) as unknown as GwasBin
    if (b.bin_start % expect.bin_bp || b.bin_end !== b.bin_start + expect.bin_bp || b.lead_position < b.bin_start || b.lead_position >= b.bin_end)
      throw bad(`bin ${b.chr}:${b.bin_start}-${b.bin_end} is not one ${expect.bin_bp} bp bin holding its lead ${b.lead_position}`)
    out.push(b)
  }
  return out
}

// ---- inverse Student t ----------------------------------------------------------------------------

const HALF_LN_2PI = 0.5 * Math.log(2 * Math.PI)

/** ln Gamma(x) for x > 0: shift up to x >= 15, then Stirling's series through x^-13 (truncation below 1e-19). */
function lgamma(x: number): number {
  let shift = 0
  while (x < 15) { shift += Math.log(x); x += 1 }
  const r = 1 / (x * x)
  const series = (1 / 12 + r * (-1 / 360 + r * (1 / 1260 + r * (-1 / 1680 + r * (1 / 1188 + r * (-691 / 360360 + r / 156)))))) / x
  return (x - 0.5) * Math.log(x) - x + HALF_LN_2PI + series - shift
}

/** Continued fraction of the regularized incomplete beta (modified Lentz). */
function betacf(a: number, b: number, x: number): number {
  const TINY = 1e-300
  const qab = a + b, qap = a + 1, qam = a - 1
  let c = 1, d = 1 - (qab * x) / qap
  if (Math.abs(d) < TINY) d = TINY
  d = 1 / d
  let h = d
  for (let m = 1; m <= 10000; m++) {
    const m2 = 2 * m
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2))
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY
    c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY
    d = 1 / d; h *= d * c
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2))
    d = 1 + aa * d; if (Math.abs(d) < TINY) d = TINY
    c = 1 + aa / c; if (Math.abs(c) < TINY) c = TINY
    d = 1 / d
    const del = d * c
    h *= del
    if (Math.abs(del - 1) <= 1e-15) break
  }
  return h
}

interface TConsts { nu: number; a: number; lnB: number; lnPdf0: number; split: number }
const T_CONSTS = new Map<number, TConsts>()
function tConsts(dof: number): TConsts {
  let k = T_CONSTS.get(dof)
  if (!k) {
    const a = dof / 2
    k = { nu: dof, a, lnB: lgamma(a) + lgamma(0.5) - lgamma(a + 0.5),
      lnPdf0: lgamma((dof + 1) / 2) - lgamma(dof / 2) - 0.5 * Math.log(dof * Math.PI), split: (a + 1) / (a + 2.5) }
    T_CONSTS.set(dof, k)
  }
  return k
}

/** ln of the two-sided p-value of |t| = t: I_x(nu/2, 1/2) with x = nu / (nu + t^2), in logs throughout. */
function logP(t: number, k: TConsts): number {
  if (t === 0) return 0
  const r = (t * t) / k.nu
  const lnx = -Math.log1p(r)                          // ln(nu / (nu + t^2))
  const lny = 2 * Math.log(t) - Math.log(k.nu + t * t) // ln(t^2 / (nu + t^2))
  const x = Math.exp(lnx)
  if (x < k.split) return k.a * lnx + 0.5 * lny - Math.log(k.a) - k.lnB + Math.log(betacf(k.a, 0.5, x))
  // I_x(a, b) = 1 - I_y(b, a)
  const iy = Math.exp(0.5 * lny + k.a * lnx - Math.log(0.5) - k.lnB + Math.log(betacf(0.5, k.a, -Math.expm1(lnx))))
  return Math.log1p(-iy)
}

const logPdf = (t: number, k: TConsts) => k.lnPdf0 - ((k.nu + 1) / 2) * Math.log1p((t * t) / k.nu)

// normal quantile (Acklam) for a starting point
const NA = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239]
const NB = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1]
const NC = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783]
const ND = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416]

/** |z| whose two-sided normal p is exp(L). */
function zFromLogP(L: number): number {
  const lq = L - Math.LN2                             // ln of the one-sided tail
  if (lq < Math.log(0.02425)) {
    const s = Math.sqrt(-2 * lq)
    return -(((((NC[0] * s + NC[1]) * s + NC[2]) * s + NC[3]) * s + NC[4]) * s + NC[5]) / ((((ND[0] * s + ND[1]) * s + ND[2]) * s + ND[3]) * s + 1)
  }
  const q = Math.expm1(L) / 2                         // one-sided tail - 0.5, without cancellation near p = 1
  const r = q * q
  return -((((((NA[0] * r + NA[1]) * r + NA[2]) * r + NA[3]) * r + NA[4]) * r + NA[5]) * q) / (((((NB[0] * r + NB[1]) * r + NB[2]) * r + NB[3]) * r + NB[4]) * r + 1)
}

/** |t| whose two-sided Student t p-value with `dof` degrees of freedom is 10^-nlp: the reader's
 *  max(0, -stdtrit(dof, p / 2)), computed in log space so it holds from nlp = 0 to beyond 320. */
export function tFromNlp(nlp: number, dof: number): number {
  if (Number.isNaN(nlp)) return NaN
  if (nlp <= 0) return 0
  if (nlp === Infinity) return Infinity
  const k = tConsts(dof)
  const L = -nlp * Math.LN10
  const z = zFromLogP(L)
  let t = Math.max(z + (z * z * z + z) / (4 * dof), Number.MIN_VALUE)
  let lo = 0, hi = Infinity
  for (let it = 0; it < 200; it++) {
    const lp = logP(t, k)
    const g = lp - L
    if (g === 0) return t
    if (g > 0) lo = t
    else hi = t
    const deriv = -2 * Math.exp(logPdf(t, k) - lp)   // d ln p / dt
    let next = t - g / deriv
    if (!(next > lo && next < hi)) next = hi === Infinity ? 2 * t + 1 : (lo + hi) / 2
    if (Math.abs(next - t) <= 1e-12 * next) return next
    if (hi < Infinity && hi - lo <= 1e-15 * hi) return next
    t = next
  }
  return t
}
