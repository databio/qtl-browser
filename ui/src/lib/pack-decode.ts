/**
 * Decoder for the qtlb pack format, version 0 (SPEC.md): eQTL and sQTL result blocks, variant
 * pages, GWAS blocks and the GWAS index, trans frames, the per-phenotype reader columns (SPEC
 * section 8), and the inverse Student t that derives the cis slope and the trans SE and r2.
 *
 * Pure functions only: no DuckDB, no import.meta.env, no relative imports, and only erasable
 * TypeScript, so `npm run pack-check` (Node) runs this exact file against the pipeline's reference
 * rows. A range that breaks a rule of SPEC section 7, 11, or 12 throws a PackError naming the range
 * and the rule.
 */
import { decompress } from 'fzstd'

export class PackError extends Error {}

/** The format and version this file reads. `manifest.json` must name exactly these (SPEC section 3). */
export const PACK_FORMAT = 'qtlb'
export const PACK_VERSION = 0

if (new Uint8Array(new Uint16Array([1]).buffer)[0] !== 1) throw new PackError('qtlb packs are little-endian; this platform is not')

/** What a result block must match (SPEC section 7 steps 3 and 5): the range length for both kinds,
 *  and for a gene's kind 2 block the `search_index` run fields. */
export interface BlockExpect { blk_len: number | null; n_var?: number | null; var_start?: number | null }

export interface ResultBlock {
  kind: 2 | 3
  nRows: number; varStart: number | null; anchor: number; posFirst: number; posLast: number
  nlpMax: number; lseMin: number; lseMax: number
  /** kind 2: the gene details JSON (SPEC section 5), { v, gene, exons, splice }; kind 3: null */
  details: Record<string, unknown> | null
  nlpCode: Uint16Array; seCode: Uint16Array
  nlp: Float64Array            // -log10 p; NaN null, Infinity for p = 0
  pval: Float64Array           // NaN null, 0 for p = 0
  se: Float64Array             // exp(lse_min + q * step); NaN null
  slope: Float64Array          // sign * se * t(p, dof); NaN where p is null or 0, or the SE code is null
  /** credible-set records, one per membership, sorted by (row, cs_id) */
  csRow: Uint32Array; csPip: Float32Array; csId: Uint8Array
}

interface VariantColumns {
  n: number
  position: Uint32Array
  rsNumber: Uint32Array        // 0 = no rsID
  afCode: Uint16Array          // 65535 = null
  af: Float64Array             // code / 65534; NaN null
  maSamples: Uint16Array; maCount: Uint16Array   // 65535 = null
  match: Uint8Array            // 0 none, 1 exact, 2 position
}
/** Every record of a fetched variants range; record i has vidx `first + i`. A1 and A2 are null for a
 *  record with flags bit 2 (alleles not reported), which only a trans-only variant carries. */
export interface VariantRange extends VariantColumns { first: number; a1: (string | null)[]; a2: (string | null)[] }
/** Exactly one phenotype's rows [varStart, varStart + n) of a variants range: cis records, all with alleles. */
export interface VariantRun extends VariantColumns { varStart: number; a1: string[]; a2: string[] }

/** SPEC section 8 as typed arrays, with in-band nulls: NaN for floats, -1 for counts and cs_id, 0 for rs_number. */
export interface ReaderColumns {
  n: number
  position: Int32Array; a1: string[]; a2: string[]; rsNumber: Uint32Array; tssDistance: Int32Array
  af: Float32Array; maSamples: Int16Array; maCount: Int16Array
  pval: Float64Array; slope: Float32Array; se: Float32Array; pip: Float32Array; csId: Int8Array
}

/** The GWAS index (SPEC section 11, kind 5). */
export interface GwasIndex {
  blockRows: number
  nValues: Uint32Array
  chroms: Map<string, { firstPosition: Uint32Array; endOffset: Uint32Array }>
}
/** One range request of a GWAS pack: bytes [off, off + len), holding blocks firstBlock..lastBlock. */
export interface GwasRange { off: number; len: number; firstBlock: number; lastBlock: number }
/** GWAS rows with lo <= position <= hi, in file order. */
export interface GwasColumns {
  rows: number
  position: Int32Array; ea: string[]; nea: string[]
  beta: Float64Array; se: Float64Array; eaf: Float64Array; p: Float64Array
  rsNumber: Uint32Array        // 0 = no rsID
  n: Int32Array
}

const NLP_NULL = 65535, NLP_ZERO = 65534
export const NLP_MAXQ = 65533
const SE_NULL = 0xffff, SE_SIGN = 0x8000, SE_QMASK = 0x7fff, SE_MAXQ = 32766
const AF_NULL = 65535
export const AF_MAXQ = 65534, COUNT_NULL = 65535
const SNP_A1 = ['', 'A', 'A', 'A', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T', 'T']
const SNP_A2 = ['', 'C', 'G', 'T', 'A', 'G', 'T', 'A', 'C', 'T', 'A', 'C', 'G']
const SNP_PAIRS = new Set(SNP_A1.slice(1).map((a, i) => `${a}\t${SNP_A2[i + 1]}`))
const MAGIC_FILE = 0x424c5451      // "QTLB" read as a little-endian u32
const MAGIC_BLOCK = 0x30424751     // "QGB0"
const MAGIC_ZSTD = 0xfd2fb528
const MAGIC_TRANS = 0x30545451     // "QTT0"
export const BETA_MAXQ = 32767
const FILE_HEADER_LEN = 32
const GWAS_SCALE = 10000
const GWAS_ROW_BYTES = 21
// the double nearest each power of ten, as SPEC section 11 divides by it (p_exp is -128 to -3)
const POW10 = Array.from({ length: 129 }, (_, k) => Number(`1e${k}`))
const UTF8 = new TextDecoder('utf-8', { fatal: true, ignoreBOM: true })

/** Typed-array views need a byteOffset that is a multiple of 4 (a Node Buffer from a pool may not have one). */
const aligned = (b: Uint8Array) => (b.byteOffset % 4 === 0 ? b : b.slice())
const view = (b: Uint8Array) => new DataView(b.buffer, b.byteOffset, b.byteLength)
const pad4 = (x: number) => Math.ceil(x / 4) * 4

/** One zstd frame (SPEC section 2): content size present, checksum flag on, no dictionary, and no
 *  bytes after the frame (found by walking the block headers). The XXH64 value itself is not
 *  verified (SPEC section 7). Returns the decompressed bytes. */
function unframe(b: Uint8Array, expect: number | null, what: string): Uint8Array {
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

/** Allele codes and heap (SPEC section 4), shared by variant pages and GWAS blocks. A variant page
 *  passes its flags: a record with bit 2 (alleles not reported) has code 0, no heap record, and null alleles. */
function decodeAlleles(codes: Uint8Array, heap: Uint8Array, bad: (rule: string) => PackError,
  flags?: Uint8Array): { a1: (string | null)[]; a2: (string | null)[] } {
  const n = codes.length
  const a1: (string | null)[] = new Array(n), a2: (string | null)[] = new Array(n)
  let h = 0
  for (let i = 0; i < n; i++) {
    const c = codes[i]
    if (c > 12) throw bad(`reserved allele code ${c} at record ${i}`)
    if (flags && flags[i] & 4) {
      if (c) throw bad(`record ${i} has flags bit 2 (alleles not reported) with allele code ${c}`)
      a1[i] = null; a2[i] = null
      continue
    }
    if (c) { a1[i] = SNP_A1[c]; a2[i] = SNP_A2[c]; continue }
    let tab = -1, j = h
    for (; j < heap.length && heap[j] !== 10; j++) {
      if (heap[j] >= 0x80) throw bad('heap is not ASCII')
      if (heap[j] === 9) { if (tab >= 0) throw bad('heap record has two tabs'); tab = j }
    }
    if (j >= heap.length || tab < 0) throw bad(`heap record for record ${i} is missing its tab or newline`)
    a1[i] = String.fromCharCode(...heap.subarray(h, tab))
    a2[i] = String.fromCharCode(...heap.subarray(tab + 1, j))
    if (SNP_PAIRS.has(`${a1[i]}\t${a2[i]}`)) throw bad(`SNP ${a1[i]}/${a2[i]} is in the heap instead of a code`)
    h = j + 1
  }
  if (h !== heap.length) throw bad('heap holds bytes beyond its code-0 records')
  return { a1, a2 }
}

/** SPEC section 5: a kind 2 block's details frame, one zstd frame of JSON. The cis scan reads this
 *  straight out of a span without decoding the block's rows, so it lives on its own. */
export function decodeDetails(frame: Uint8Array, expectLen: number, what = 'block'): Record<string, unknown> {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const raw = unframe(frame, expectLen, `${what} details`)
  if (raw[0] === 0xef && raw[1] === 0xbb && raw[2] === 0xbf) throw bad('details JSON starts with a byte-order mark')
  let parsed: unknown
  try { parsed = JSON.parse(UTF8.decode(raw)) } catch (e) { throw bad(`details JSON: ${(e as Error).message}`) }   // JSON.parse rejects NaN and Infinity
  if (parsed === null || typeof parsed !== 'object' || Array.isArray(parsed) || (parsed as { v?: unknown }).v !== 0)
    throw bad('details are not a JSON object with "v": 0')
  return parsed as Record<string, unknown>
}

/** SPEC section 7 step 3 for a kind 2 (eQTL, with details) or kind 3 (sQTL) block, then the section 5
 *  values with `dof` (manifest packs.dof.eqtl or .sqtl). */
export function decodeResultBlock(bytes: Uint8Array, kind: 2 | 3, dof: number, expect: BlockExpect, what = 'result block'): ResultBlock {
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
  if (blkLen !== b.length || blkLen !== expect.blk_len) throw bad(`length field ${blkLen} != range ${b.length} / expected blk_len ${expect.blk_len}`)
  const body = 64 + 4 * n + 12 * nCs + dz
  if (Math.ceil(body / 4) * 4 !== blkLen) throw bad(`64 + 4 n_rows + 12 n_cs + details_zlen = ${body}, padded to 4, != block length ${blkLen}`)
  for (let i = body; i < blkLen; i++) if (b[i]) throw bad('padding is not zero')
  if (kind === 2) {
    if (!dz) throw bad('a kind 2 block needs a details frame (details_zlen above 0)')
    if (n !== (expect.n_var ?? 0)) throw bad(`n_rows ${n} != search_index n_var ${expect.n_var}`)
    if (n === 0) {
      if (varStart !== 0xffffffff || anchor || posFirst || posLast || nCs || nlpMax !== 0 || lseMin !== 0 || lseMax !== 0 || expect.var_start != null)
        throw bad('a block with no rows needs var_start 0xFFFFFFFF and zero anchor, positions, credible sets, and scales')
    } else if (varStart !== expect.var_start) throw bad(`var_start ${varStart} != search_index var_start ${expect.var_start}`)
  } else {
    if (dz || dl) throw bad('a kind 3 block has no details (details_zlen and details_len 0)')
    if (n === 0) throw bad('a kind 3 block needs at least one row')
  }
  if (n > 0) {
    if (varStart + n - 1 >= 0xffffffff) throw bad('var_start + n_rows - 1 reaches 0xFFFFFFFF')
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

  const details = kind === 2 ? decodeDetails(b.subarray(64 + 4 * n + 12 * nCs, body), dl, what) : null

  const nlp = new Float64Array(n), pval = new Float64Array(n), se = new Float64Array(n), slope = new Float64Array(n)
  const nlpStep = nlpMax / NLP_MAXQ, seStep = (lseMax - lseMin) / SE_MAXQ
  const tByCode = new Map<number, number>()
  for (let i = 0; i < n; i++) {
    const q = nlpCode[i], s = seCode[i]
    if (q === NLP_NULL) { nlp[i] = NaN; pval[i] = NaN }
    else if (q === NLP_ZERO) { nlp[i] = Infinity; pval[i] = 0 }
    else { nlp[i] = q * nlpStep; pval[i] = Math.pow(10, -nlp[i]) }
    if (s === SE_NULL) { se[i] = NaN; slope[i] = NaN; continue }
    se[i] = Math.exp(lseMin + (s & SE_QMASK) * seStep)
    if (q >= NLP_ZERO) { slope[i] = NaN; continue }
    let t = tByCode.get(q)
    if (t === undefined) { t = tFromNlp(nlp[i], dof); tByCode.set(q, t) }
    slope[i] = ((s & SE_SIGN) ? -1 : 1) * se[i] * t
  }
  return { kind, nRows: n, varStart: n ? varStart : null, anchor, posFirst, posLast, nlpMax, lseMin, lseMax,
    details, nlpCode, seCode, nlp, pval, se, slope, csRow, csPip, csId }
}

/** Every credible-set membership of a block, for the credible-set table (a row in two sets appears twice). */
export function csMembers(block: ResultBlock): { row: number; pip: number; csId: number }[] {
  return Array.from(block.csRow, (row, k) => ({ row, pip: block.csPip[k], csId: block.csId[k] }))
}

/** SPEC section 7 step 4: walk the pages of a gene's variants range and decode every record in them. */
export function decodeVariantRange(bytes: Uint8Array, expect: { var_len: number | null }, what = 'variants range'): VariantRange {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (bytes.length !== expect.var_len) throw bad(`range is ${bytes.length} bytes, expected var_len ${expect.var_len}`)
  const b = aligned(bytes)
  const dv = view(b)
  const pages: { first: number; n: number; payload: Uint8Array; a1: (string | null)[]; a2: (string | null)[] }[] = []
  let off = 0, codec0 = -1, nextVidx = -1, lastPos = 0, total = 0
  while (off < b.length) {
    const where = `page at byte ${off}`
    const pbad = (rule: string) => bad(`${where}: ${rule}`)
    if (b.length - off < 12) throw pbad('header is truncated')
    const stored = dv.getUint32(off, true), first = dv.getUint32(off + 4, true), n = dv.getUint16(off + 8, true)
    const codec = b[off + 10]
    if (b[off + 11]) throw pbad('reserved byte is not zero')
    if (codec > 1 || (codec0 >= 0 && codec !== codec0)) throw pbad(`codec ${codec} (first page ${codec0})`)
    codec0 = codec
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
    for (let i = 0; i < n; i++) if (flags[i] & 0xf8 || (flags[i] & 3) === 3) throw pbad(`reserved flags at record ${i}`)
    const { a1, a2 } = decodeAlleles(payload.subarray(4 + 14 * n, 4 + 15 * n), payload.subarray(4 + 16 * n), pbad, flags)
    pages.push({ first, n, payload, a1, a2 })
    nextVidx = first + n
    total += n
    off = stop
  }
  if (!pages.length) throw bad('no pages')

  const r: VariantRange = { first: pages[0].first, n: total, position: new Uint32Array(total), a1: new Array(total), a2: new Array(total),
    rsNumber: new Uint32Array(total), afCode: new Uint16Array(total), af: new Float64Array(total),
    maSamples: new Uint16Array(total), maCount: new Uint16Array(total), match: new Uint8Array(total) }
  let k = 0
  for (const { n, payload: p, a1, a2 } of pages) {
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
      r.match[k] = p[4 + 15 * n + i] & 3
      r.a1[k] = a1[i]; r.a2[k] = a2[i]
    }
  }
  return r
}

/** SPEC section 7 steps 4 and 5: one phenotype's run inside a decoded range, whose end positions
 *  must equal the block's pos_first and pos_last, and whose records all have alleles (a run never
 *  reaches n_cis). The typed arrays are views on the range. */
export function sliceRun(range: VariantRange, varStart: number, n: number, block: { posFirst: number; posLast: number },
  what = 'variants range'): VariantRun {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const i0 = varStart - range.first
  if (!(n >= 1) || i0 < 0 || i0 + n > range.n)
    throw bad(`run ${varStart}..${varStart + n - 1} is not inside the decoded records ${range.first}..${range.first + range.n - 1}`)
  const position = range.position.subarray(i0, i0 + n)
  if (position[0] !== block.posFirst || position[n - 1] !== block.posLast)
    throw bad(`positions ${position[0]}..${position[n - 1]} != block pos_first ${block.posFirst}, pos_last ${block.posLast}`)
  const a1 = range.a1.slice(i0, i0 + n), a2 = range.a2.slice(i0, i0 + n)
  const missing = a1.indexOf(null)
  if (missing >= 0) throw bad(`record ${varStart + missing} has no alleles (flags bit 2): a gene or intron run never reaches n_cis`)
  return { varStart, n, position, a1: a1 as string[], a2: a2 as string[],
    rsNumber: range.rsNumber.subarray(i0, i0 + n), afCode: range.afCode.subarray(i0, i0 + n), af: range.af.subarray(i0, i0 + n),
    maSamples: range.maSamples.subarray(i0, i0 + n), maCount: range.maCount.subarray(i0, i0 + n), match: range.match.subarray(i0, i0 + n) }
}

/** SPEC section 8: the block joined to its variant run by row, with in-band nulls. A row in two
 *  credible sets takes its higher-PIP record, and the lower cs_id on a tie. */
export function readerColumns(block: ResultBlock, run: VariantRun): ReaderColumns {
  const n = block.nRows
  if (run.n !== n || run.varStart !== block.varStart) throw new PackError(`reader: variant run ${run.varStart}+${run.n} does not match block ${block.varStart}+${n}`)
  const c: ReaderColumns = { n, position: new Int32Array(n), a1: run.a1, a2: run.a2, rsNumber: run.rsNumber, tssDistance: new Int32Array(n),
    af: new Float32Array(n), maSamples: new Int16Array(n), maCount: new Int16Array(n), pval: block.pval,
    slope: Float32Array.from(block.slope), se: Float32Array.from(block.se), pip: new Float32Array(n).fill(NaN), csId: new Int8Array(n).fill(-1) }
  for (let i = 0; i < n; i++) {
    const pos = run.position[i], tss = pos - block.anchor
    if (pos > 0x7fffffff || tss > 0x7fffffff || tss < -0x80000000) throw new PackError('reader: position or tss_distance does not fit INTEGER')
    c.position[i] = pos; c.tssDistance[i] = tss
    c.af[i] = run.af[i]
    const ms = run.maSamples[i], mc = run.maCount[i]
    if ((ms !== COUNT_NULL && ms > 32767) || (mc !== COUNT_NULL && mc > 32767)) throw new PackError('reader: a count does not fit SMALLINT')
    c.maSamples[i] = ms === COUNT_NULL ? -1 : ms
    c.maCount[i] = mc === COUNT_NULL ? -1 : mc
  }
  // records are sorted by (row, cs_id), so a strict > keeps the lower cs_id on a PIP tie
  for (let k = 0; k < block.csRow.length; k++) {
    const row = block.csRow[k]
    if (c.csId[row] < 0 || block.csPip[k] > c.pip[row]) { c.pip[row] = block.csPip[k]; c.csId[row] = block.csId[k] }
  }
  return c
}

// ---- GWAS (SPEC section 11) -----------------------------------------------------------------------

/** An 8-byte ASCII name, zero-padded. */
function asciiName(b: Uint8Array, bad: (rule: string) => PackError): string {
  let end = b.indexOf(0)
  if (end < 0) end = b.length
  for (let i = end; i < b.length; i++) if (b[i]) throw bad('name is not zero-padded')
  if (end === 0 || b.subarray(0, end).some(x => x >= 0x80)) throw bad('name is empty or not ASCII')
  return String.fromCharCode(...b.subarray(0, end))
}

/** The kind 5 file (header, then one zstd frame to the end of the file). */
export function parseGwasIndex(bytes: Uint8Array, what = 'GWAS index'): GwasIndex {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  if (b.length < FILE_HEADER_LEN) throw bad('shorter than the 32-byte file header')
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_FILE) throw bad('magic is not QTLB')
  if (b[4] !== 5 || b[5] !== 0 || dv.getUint16(6, true) !== FILE_HEADER_LEN)
    throw bad(`kind ${b[4]}, version ${b[5]}, header length ${dv.getUint16(6, true)}; expected 5, 0, 32`)
  if (asciiName(b.subarray(8, 16), bad) !== 'all') throw bad('header chromosome is not "all"')
  const count = dv.getUint32(16, true), blockRows = dv.getUint32(20, true)
  if (blockRows < 1 || blockRows > 65535) throw bad(`rows per block ${blockRows} outside 1..65535`)
  for (let i = 24; i < FILE_HEADER_LEN; i++) if (b[i]) throw bad('reserved header bytes are not zero')
  const p = unframe(b.subarray(FILE_HEADER_LEN), null, what)
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
  if (nChroms !== count) throw bad(`${nChroms} chromosomes, header count ${count}`)
  const chroms = new Map<string, { firstPosition: Uint32Array; endOffset: Uint32Array }>()
  for (let c = 0; c < nChroms; c++) {
    if (off + 8 > p.length) throw bad('payload ends early')
    const name = asciiName(p.subarray(off, off + 8), bad)
    off += 8
    if (chroms.has(name)) throw bad(`chromosome ${name} appears twice`)
    const nb = u32s(1)[0]
    const firstPosition = u32s(nb), endOffset = u32s(nb)
    if (nb < 1 || firstPosition[0] < 1 || endOffset[0] <= FILE_HEADER_LEN ||
        firstPosition.some((x, i) => i > 0 && x < firstPosition[i - 1]) || endOffset.some((x, i) => i > 0 && x <= endOffset[i - 1]))
      throw bad(`${name}: blocks must be non-empty, first positions non-decreasing from 1, end offsets increasing past byte 32`)
    chroms.set(name, { firstPosition, endOffset })
  }
  if (off !== p.length) throw bad(`${p.length - off} bytes after the last chromosome`)
  return { blockRows, nValues, chroms }
}

/** Number of entries of the sorted array `a` below `x` (strict) or at most `x`. */
function countBelow(a: Uint32Array, x: number, orEqual: boolean): number {
  let lo = 0, hi = a.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (a[mid] < x || (orEqual && a[mid] === x)) lo = mid + 1
    else hi = mid
  }
  return lo
}

/** SPEC section 11's window rule: the one byte range holding every row with lo <= position <= hi,
 *  or null when the chromosome has no GWAS rows (chrX) or no block starts at or before hi. */
export function gwasRange(index: GwasIndex, chr: string, lo: number, hi: number): GwasRange | null {
  if (lo > hi) throw new PackError(`GWAS window ${chr}:${lo}-${hi}: lo is above hi`)
  const c = index.chroms.get(chr)
  if (!c) return null
  const end = countBelow(c.firstPosition, hi, true) - 1
  if (end < 0) return null
  const start = Math.max(countBelow(c.firstPosition, lo, false) - 1, 0)
  const off = start === 0 ? FILE_HEADER_LEN : c.endOffset[start - 1]
  return { off, len: c.endOffset[end] - off, firstBlock: start, lastBlock: end }
}

/** Decode the blocks of one GWAS range (SPEC section 11 reader rules) and keep the rows with
 *  lo <= position <= hi. */
export function decodeGwasRange(bytes: Uint8Array, index: GwasIndex, chr: string, range: GwasRange, lo: number, hi: number,
  what = `GWAS ${chr}`): GwasColumns {
  const bad = (rule: string) => new PackError(`${what} bytes ${range.off}+${range.len}: ${rule}`)
  const c = index.chroms.get(chr)
  if (!c) throw bad('the index has no such chromosome')
  const { firstBlock: s, lastBlock: e } = range
  const nBlocks = c.firstPosition.length
  if (!(s >= 0 && s <= e && e < nBlocks)) throw bad(`blocks ${s}..${e} outside the index's ${nBlocks}`)
  const start = s === 0 ? FILE_HEADER_LEN : c.endOffset[s - 1]
  if (range.off !== start || range.len !== c.endOffset[e] - start) throw bad('the range does not match its blocks in the index')
  if (bytes.length !== range.len) throw bad(`got ${bytes.length} bytes`)
  const nv = index.nValues
  const blocks: { n: number; p: Uint8Array; pos: Uint32Array; a: { a1: string[]; a2: string[] } }[] = []
  let lastPos = 0, keep = 0
  for (let k = s; k <= e; k++) {
    const kbad = (rule: string) => bad(`block ${k}: ${rule}`)
    const from = (k === 0 ? FILE_HEADER_LEN : c.endOffset[k - 1]) - range.off, to = c.endOffset[k] - range.off
    const p = aligned(unframe(bytes.subarray(from, to), null, `${what} block ${k}`))
    if (p.length < 8) throw kbad('payload is shorter than its 8-byte header')
    const dv = view(p)
    const n = dv.getUint32(0, true), heapLen = dv.getUint32(4, true)
    if (n === 0) throw kbad('no rows')
    if (k === nBlocks - 1 ? n > index.blockRows : n !== index.blockRows)
      throw kbad(`${n} rows; every block but the chromosome's last holds ${index.blockRows}`)
    if (p.length !== 8 + GWAS_ROW_BYTES * n + heapLen) throw kbad(`payload ${p.length} != 8 + 21n + heap_len`)
    const at = (o: number) => p.byteOffset + o
    const deltas = new Uint32Array(p.buffer, at(8), n)
    const pos = new Uint32Array(n)
    let acc = 0
    for (let i = 0; i < n; i++) {
      acc = i === 0 ? deltas[0] : acc + deltas[i]
      pos[i] = acc
      if (acc >= lo && acc <= hi) keep++
    }
    if (pos[0] !== c.firstPosition[k]) throw kbad(`first position ${pos[0]} != index ${c.firstPosition[k]}`)
    if (pos[0] < 1 || pos[0] < lastPos || acc > 0x7fffffff) throw kbad('positions below 1, decreasing, or beyond INTEGER')
    lastPos = acc
    const eaf = new Uint16Array(p.buffer, at(8 + 14 * n), n), mant = new Uint16Array(p.buffer, at(8 + 16 * n), n)
    const exp = new Int8Array(p.buffer, at(8 + 18 * n), n), nCode = p.subarray(8 + 19 * n, 8 + 20 * n)
    for (let i = 0; i < n; i++) {
      if (eaf[i] > GWAS_SCALE) throw kbad(`eaf code ${eaf[i]} above 10000 at row ${i}`)
      if (mant[i] < 1000 || mant[i] > 9999) throw kbad(`p mantissa ${mant[i]} outside 1000..9999 at row ${i}`)
      if (exp[i] > -3 || (exp[i] === -3 && mant[i] !== 1000)) throw kbad(`p above 1 at row ${i}`)
      if (nCode[i] >= nv.length) throw kbad(`n code ${nCode[i]} beyond the ${nv.length}-value n table at row ${i}`)
    }
    // no flags, so every allele is set
    const a = decodeAlleles(p.subarray(8 + 20 * n, 8 + 21 * n), p.subarray(8 + 21 * n), kbad) as { a1: string[]; a2: string[] }
    blocks.push({ n, p, pos, a })
  }

  const out: GwasColumns = { rows: keep, position: new Int32Array(keep), ea: new Array(keep), nea: new Array(keep),
    beta: new Float64Array(keep), se: new Float64Array(keep), eaf: new Float64Array(keep), p: new Float64Array(keep),
    rsNumber: new Uint32Array(keep), n: new Int32Array(keep) }
  let j = 0
  for (const { n, p, pos, a } of blocks) {
    const at = (o: number) => p.byteOffset + o
    const beta = new Int32Array(p.buffer, at(8 + 4 * n), n), rs = new Uint32Array(p.buffer, at(8 + 8 * n), n)
    const se = new Uint16Array(p.buffer, at(8 + 12 * n), n), eaf = new Uint16Array(p.buffer, at(8 + 14 * n), n)
    const mant = new Uint16Array(p.buffer, at(8 + 16 * n), n), exp = new Int8Array(p.buffer, at(8 + 18 * n), n)
    for (let i = 0; i < n; i++) {
      if (pos[i] < lo || pos[i] > hi) continue
      out.position[j] = pos[i]; out.ea[j] = a.a1[i]; out.nea[j] = a.a2[i]
      out.beta[j] = beta[i] / GWAS_SCALE; out.se[j] = se[i] / GWAS_SCALE; out.eaf[j] = eaf[i] / GWAS_SCALE
      out.p[j] = mant[i] / POW10[-exp[i]]
      out.rsNumber[j] = rs[i]
      out.n[j] = nv[p[8 + 19 * n + i]]
      j++
    }
  }
  return out
}

// ---- trans pack (SPEC section 12) -----------------------------------------------------------------

/** One gene's trans frame. Rows are the eQTL rows (0 .. nE - 1), then the sQTL rows grouped by intron. */
export interface TransFrame {
  nE: number; nS: number; k: number; nlpMax: number; betaMax: number
  /** the intron table, strictly ascending by (start, end, cluster, strand); strand 0 '+', 1 '-' */
  intronStart: Uint32Array; intronEnd: Uint32Array; cluster: Uint32Array; strand: Uint8Array
  variantChr: Uint8Array       // 1..22, 23 = chrX (TRANS_CHROMS)
  position: Uint32Array        // absolute
  rsNumber: Uint32Array        // 0 = no rsID
  af: Float64Array             // code / 65534
  nlp: Float64Array; pval: Float64Array; beta: Float64Array
  betaSe: Float64Array         // |beta| / t, with t from -log10 p and the row type's dof
  r2: Float64Array             // t^2 / (t^2 + dof)
  intron: Uint8Array           // length nS: each sQTL row's index into the intron table
}

/** Variant chromosome names by trans `variant_chr` code (index 0 is unused). */
export const TRANS_CHROMS = ['', ...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX']

/** One gene's frame (its `trans_off`, `trans_len` bytes) under the reader rules of SPEC section 12,
 *  with `dof` from manifest packs.dof: `e` for eQTL rows, `s` for sQTL rows. */
export function decodeTransFrame(bytes: Uint8Array, dof: { e: number; s: number }, what = 'trans frame'): TransFrame {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const p = aligned(unframe(bytes, null, what))
  if (p.length < 32) throw bad(`decompresses to ${p.length} bytes, fewer than 32`)
  const dv = view(p)
  if (dv.getUint32(0, true) !== MAGIC_TRANS) throw bad('magic is not QTT0')
  const nE = dv.getUint32(4, true), nS = dv.getUint32(8, true), k = dv.getUint16(12, true)
  const nlpMax = dv.getFloat64(16, true), betaMax = dv.getFloat64(24, true)
  if (dv.getUint16(14, true)) throw bad('reserved field is not zero')
  const n = nE + nS
  if (n === 0) throw bad('n_e + n_s is 0')
  if (nS === 0 ? k !== 0 : k < 1 || k > 255) throw bad(`k ${k} with n_s ${nS}: k is 0 without sQTL rows, else 1 to 255`)
  if (!(Number.isFinite(nlpMax) && nlpMax >= 0 && Number.isFinite(betaMax) && betaMax >= 0))
    throw bad(`scales nlp_max ${nlpMax}, beta_max ${betaMax} are not finite and at least 0`)
  const R = pad4(32 + 13 * k), C = pad4(R + 14 * n), used = C + n + nS
  if (p.length !== pad4(used)) throw bad(`payload is ${p.length} bytes; n_e ${nE}, n_s ${nS}, k ${k} lay out ${pad4(used)}`)
  const zeros = (from: number, to: number) => { for (let i = from; i < to; i++) if (p[i]) throw bad(`padding byte ${i} is not zero`) }
  zeros(32 + 13 * k, R); zeros(R + 14 * n, C); zeros(used, p.length)

  const at = (o: number) => p.byteOffset + o
  const intronStart = new Uint32Array(p.buffer, at(32), k), intronEnd = new Uint32Array(p.buffer, at(32 + 4 * k), k)
  const cluster = new Uint32Array(p.buffer, at(32 + 8 * k), k), strand = p.subarray(32 + 12 * k, 32 + 13 * k)
  for (let j = 0; j < k; j++) {
    if (strand[j] > 1) throw bad(`intron ${j}: strand code ${strand[j]} above 1`)
    if (j > 0 && (intronStart[j] - intronStart[j - 1] || intronEnd[j] - intronEnd[j - 1] || cluster[j] - cluster[j - 1] || strand[j] - strand[j - 1]) <= 0)
      throw bad(`intron table is not strictly ascending by (start, end, cluster, strand) at intron ${j}`)
  }
  const intron = p.subarray(C + n, used)
  for (let j = 0; j < nS; j++) {
    const prev = j ? intron[j - 1] : -1
    if (intron[j] >= k) throw bad(`sQTL row ${j}: intron index ${intron[j]} is not below k ${k}`)
    if (intron[j] < prev) throw bad(`sQTL row ${j}: intron indices decrease`)
    if (intron[j] > prev + 1) throw bad(`intron ${prev + 1} of the table has no row`)
  }
  if (nS && intron[nS - 1] !== k - 1) throw bad(`intron ${intron[nS - 1] + 1} of the table has no row`)

  const posEntry = new Uint32Array(p.buffer, at(R), n), rsNumber = new Uint32Array(p.buffer, at(R + 4 * n), n)
  const afCode = new Uint16Array(p.buffer, at(R + 8 * n), n), nlpCode = new Uint16Array(p.buffer, at(R + 10 * n), n)
  const betaCode = new Int16Array(p.buffer, at(R + 12 * n), n), variantChr = p.subarray(C, C + n)
  const position = new Uint32Array(n), af = new Float64Array(n), nlp = new Float64Array(n), pval = new Float64Array(n)
  const beta = new Float64Array(n), betaSe = new Float64Array(n), r2 = new Float64Array(n)
  const nlpStep = nlpMax / NLP_MAXQ, betaStep = betaMax / BETA_MAXQ
  // t depends only on the nlp code and the row type, and a gene's rows repeat codes
  const tE = new Map<number, number>(), tS = new Map<number, number>()
  let pos = 0
  for (let i = 0; i < n; i++) {
    const c = variantChr[i]
    if (c < 1 || c > 23) throw bad(`row ${i}: variant_chr ${c} is outside 1 to 23`)
    const runStart = i === 0 || i === nE || (i > nE && intron[i - nE] !== intron[i - nE - 1])
    if (!runStart && c < variantChr[i - 1]) throw bad(`row ${i}: variant_chr decreases inside a run`)
    if (posEntry[i] === 0) throw bad(`row ${i}: position entry is 0`)
    pos = runStart || c !== variantChr[i - 1] ? posEntry[i] : pos + posEntry[i]
    if (pos > 0xffffffff) throw bad(`row ${i}: position exceeds 2^32 - 1`)
    position[i] = pos
    const q = nlpCode[i]
    if (q > NLP_MAXQ) throw bad(`row ${i}: nlp code ${q} above 65533`)
    if (afCode[i] > AF_MAXQ) throw bad(`row ${i}: af code ${afCode[i]} above 65534`)
    if (betaCode[i] === -32768) throw bad(`row ${i}: beta code -32768`)
    const d = i < nE ? dof.e : dof.s, cache = i < nE ? tE : tS
    af[i] = afCode[i] / AF_MAXQ
    nlp[i] = q * nlpStep
    pval[i] = Math.pow(10, -nlp[i])
    beta[i] = betaCode[i] * betaStep
    let t = cache.get(q)
    if (t === undefined) { t = tFromNlp(nlp[i], d); cache.set(q, t) }
    betaSe[i] = Math.abs(beta[i]) / t
    r2[i] = (t * t) / (t * t + d)
  }
  return { nE, nS, k, nlpMax, betaMax, intronStart, intronEnd, cluster, strand, variantChr, position, rsNumber, af, nlp, pval, beta, betaSe, r2, intron }
}

/** Each row's phenotype_id: the gene id for eQTL rows, and for sQTL rows
 *  `chr:start:end:clu_<cluster>_<strand>:gene_id.gene_version`, with the gene's (the file's) chromosome
 *  and `gene_version` from search_index. */
export function transPhenotypeIds(f: TransFrame, gene: { chr: string; gene_id: string; gene_version: number | null }): string[] {
  if (f.nS && gene.gene_version == null) throw new PackError(`${gene.gene_id}: search_index has no gene_version to rebuild its sQTL phenotype ids`)
  const introns = Array.from({ length: f.k }, (_, j) =>
    `${gene.chr}:${f.intronStart[j]}:${f.intronEnd[j]}:clu_${f.cluster[j]}_${f.strand[j] ? '-' : '+'}:${gene.gene_id}.${gene.gene_version}`)
  const out = new Array<string>(f.nE + f.nS).fill(gene.gene_id, 0, f.nE)
  for (let j = 0; j < f.nS; j++) out[f.nE + j] = introns[f.intron[j]]
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

/** Derive pval, beta, betaSe, r2 from quantized nlp and beta codes (shared by the gene-keyed trans
 *  decoder and the variant-page hits decoder). */
export function deriveTrans(
  nlpCode: number, betaCode: number, nlpMax: number, betaMax: number, dof: number,
  tCache: Map<number, number>,
): { pval: number; beta: number; betaSe: number; r2: number } {
  const nlp = nlpCode * (nlpMax / NLP_MAXQ)
  const pval = Math.pow(10, -nlp)
  const beta = betaCode * (betaMax / BETA_MAXQ)
  let t = tCache.get(nlpCode)
  if (t === undefined) { t = tFromNlp(nlp, dof); tCache.set(nlpCode, t) }
  const betaSe = Math.abs(beta) / t
  const r2 = (t * t) / (t * t + dof)
  return { pval, beta, betaSe, r2 }
}

// ---- the variant page: hits pack, rsID index, variant index (SPEC sections 13 to 15) ---------

const MAGIC_HITS = 0x30485651      // "QVH0"
const MAGIC_VARIANT_INDEX = 0x30585651   // "QVX0"
/** SPEC section 14: chr_ordinal 1..22 then 23 = chrX, the same order as the variant index. */
export const VIDX_CHROMS = TRANS_CHROMS
const MATCH_NAMES = ['none', 'exact', 'position'] as const

/** One chromosome's offsets in the startup file. `pageOff` and `hitsOff` each hold one extra
 *  entry, the file size, so a range is `off[k] .. off[k + 1] - 1`. */
export interface ChromVariantIndex {
  nCis: number; nTransOnly: number; nPagesCis: number; nPagesTrans: number
  pageOff: Uint32Array; pageFirstPos: Uint32Array; hitsOff: Uint32Array
}

/** SPEC section 15, kind 9: the whole file, fetched at startup. */
export interface VariantIndex {
  pageSize: number; frameVariants: number; rsidBlockRecords: number; rsidNRecords: number
  chroms: Map<string, ChromVariantIndex>; rsidFirst: Uint32Array
}

/** One record of a variants-file page (SPEC section 4). `A1`/`A2` are null on a trans-only record
 *  whose alleles the source never reported. */
export interface VariantRecord {
  chr: string; vidx: number; position: number; A1: string | null; A2: string | null
  rsNumber: number; af: number; match: 'exact' | 'position' | 'none'; inCis: boolean
}

/** SPEC section 13: one decoded hits frame. The row arrays are parallel; a variant's rows are
 *  `rowStart[i] .. rowStart[i + 1] - 1` for its slot `i = vidx - firstVidx`. */
export interface HitsFrame {
  firstVidx: number; nVariants: number; rowStart: Uint32Array
  trans: { nlpMax: number; betaMax: number }
  lead: { nlpMax: number; seMax: number; slopeMax: number }
  kind: Uint8Array; flags: Uint8Array; ord: Uint16Array
  v1: Uint16Array; v2: Uint16Array; v3: Int16Array
  intronStart: Uint32Array; intronEnd: Uint32Array; cluster: Uint32Array
}

/** The 32-byte file header (SPEC section 3), checked against the kind and chromosome expected. */
function fileHeader(b: Uint8Array, kind: number, chrom: string, what: string): { count: number; pageSize: number } {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (b.length < FILE_HEADER_LEN) throw bad('shorter than the 32-byte file header')
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_FILE) throw bad('magic is not QTLB')
  if (b[4] !== kind) throw bad(`kind ${b[4]} is not ${kind}`)
  const name = asciiName(b.subarray(8, 16), bad)
  if (name !== chrom) throw bad(`chromosome "${name}" is not "${chrom}"`)
  return { count: dv.getUint32(16, true), pageSize: dv.getUint32(20, true) }
}

/** SPEC section 6: the search index, an Arrow IPC stream in one zstd frame. Returns the stream
 *  bytes for DuckDB's `insertArrowFromIPCStream`; there is no parquet reader in the bundle. */
export function decodeSearchIndex(bytes: Uint8Array, what = 'search index'): Uint8Array {
  return unframe(aligned(bytes), null, what)
}

/** SPEC section 15: the startup file, header plus one zstd frame to the end. */
export function decodeVariantIndex(bytes: Uint8Array, what = 'variant index'): VariantIndex {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(bytes)
  const h = fileHeader(b, 9, 'all', what)
  const p = aligned(unframe(b.subarray(FILE_HEADER_LEN), null, what))
  if (p.length < 32) throw bad('payload shorter than its 32-byte header')
  const dv = view(p)
  if (dv.getUint32(0, true) !== MAGIC_VARIANT_INDEX) throw bad('magic is not QVX0')
  const nChrom = dv.getUint16(4, true)
  if (dv.getUint16(6, true)) throw bad('reserved field is not zero')
  const pageSize = dv.getUint32(8, true), frameVariants = dv.getUint32(12, true)
  const rsidBlockRecords = dv.getUint32(16, true), rsidNRecords = dv.getUint32(20, true)
  const rsidNBlocks = dv.getUint32(24, true)
  if (dv.getUint32(28, true)) throw bad('reserved field is not zero')
  if (nChrom !== h.count || nChrom !== 23) throw bad(`n_chrom ${nChrom} is not the header count ${h.count} or not 23`)
  if (pageSize !== h.pageSize) throw bad(`page size ${pageSize} is not the header's ${h.pageSize}`)
  if (!(pageSize >= 1 && frameVariants >= 1 && rsidBlockRecords >= 1)) throw bad('a size field is zero')
  if (rsidNBlocks !== Math.ceil(rsidNRecords / rsidBlockRecords))
    throw bad(`rsid_n_blocks ${rsidNBlocks} does not follow from ${rsidNRecords} records of ${rsidBlockRecords}`)

  const chroms = new Map<string, ChromVariantIndex>()
  let o = 32
  const u32s = (n: number, name: string): Uint32Array => {
    if (o + 4 * n > p.length) throw bad(`${name} runs past the payload`)
    const a = new Uint32Array(p.buffer, p.byteOffset + o, n)
    o += 4 * n
    return a
  }
  for (let c = 1; c <= nChrom; c++) {
    const chr = VIDX_CHROMS[c]
    if (o + 20 > p.length) throw bad(`${chr} header runs past the payload`)
    const nCis = dv.getUint32(o, true), nTransOnly = dv.getUint32(o + 4, true)
    const nPagesCis = dv.getUint32(o + 8, true), nPagesTrans = dv.getUint32(o + 12, true)
    const nFrames = dv.getUint32(o + 16, true)
    o += 20
    if (nPagesCis !== Math.ceil(nCis / pageSize) || nPagesTrans !== Math.ceil(nTransOnly / pageSize))
      throw bad(`${chr}: page counts do not follow from ${nCis} cis and ${nTransOnly} trans-only variants of ${pageSize}`)
    if (nFrames !== Math.ceil((nCis + nTransOnly) / frameVariants))
      throw bad(`${chr}: n_frames ${nFrames} does not follow from ${nCis + nTransOnly} variants of ${frameVariants}`)
    const nPages = nPagesCis + nPagesTrans
    const pageOff = u32s(nPages + 1, `${chr} page_off`)
    const pageFirstPos = u32s(nPages, `${chr} page_first_position`)
    const hitsOff = u32s(nFrames + 1, `${chr} hits_off`)
    if (pageOff[0] !== FILE_HEADER_LEN || hitsOff[0] !== FILE_HEADER_LEN)
      throw bad(`${chr}: the first page or frame does not start at byte ${FILE_HEADER_LEN}`)
    for (let i = 1; i <= nPages; i++) if (pageOff[i] <= pageOff[i - 1]) throw bad(`${chr}: page offsets do not increase`)
    for (let i = 1; i <= nFrames; i++) if (hitsOff[i] <= hitsOff[i - 1]) throw bad(`${chr}: hits frame offsets do not increase`)
    // positions never decrease within a section; the two sections are independent
    for (let i = 0; i < nPages; i++) {
      if (pageFirstPos[i] < 1) throw bad(`${chr}: a page's first position is 0`)
      if (i !== 0 && i !== nPagesCis && pageFirstPos[i] < pageFirstPos[i - 1]) throw bad(`${chr}: first positions decrease inside a section`)
    }
    chroms.set(chr, { nCis, nTransOnly, nPagesCis, nPagesTrans, pageOff, pageFirstPos, hitsOff })
  }
  const rsidFirst = u32s(rsidNBlocks, 'rsid_first')
  for (let i = 1; i < rsidNBlocks; i++) if (rsidFirst[i] <= rsidFirst[i - 1]) throw bad('rsid_first does not strictly increase')
  if (o !== p.length) throw bad(`${p.length - o} bytes left over after the payload`)
  return { pageSize, frameVariants, rsidBlockRecords, rsidNRecords, chroms, rsidFirst }
}

/** SPEC section 4: one variants-file page on its own, through the same walk the gene page uses.
 *  `nCis` decides `inCis`; the records are vidx `first .. first + n - 1`. */
export function decodeVariantPage(bytes: Uint8Array, chr: string, nCis: number, what = 'variants page'): VariantRecord[] {
  const r = decodeVariantRange(bytes, { var_len: bytes.length }, what)
  const out: VariantRecord[] = new Array(r.n)
  for (let i = 0; i < r.n; i++) {
    const vidx = r.first + i
    out[i] = { chr, vidx, position: r.position[i], A1: r.a1[i], A2: r.a2[i], rsNumber: r.rsNumber[i],
      af: r.af[i], match: MATCH_NAMES[r.match[i]], inCis: vidx < nCis }
  }
  return out
}

/** SPEC section 13: one hits frame, checked against every reader rule of that section. */
export function decodeHitsFrame(bytes: Uint8Array, firstVidx: number, what = 'hits frame'): HitsFrame {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const p = aligned(unframe(aligned(bytes), null, what))
  if (p.length < 48) throw bad(`payload is ${p.length} bytes, fewer than 48`)
  const dv = view(p)
  if (dv.getUint32(0, true) !== MAGIC_HITS) throw bad('magic is not QVH0')
  const V = dv.getUint16(4, true)
  if (dv.getUint16(6, true)) throw bad('reserved field is not zero')
  if (V < 1) throw bad('n_variants is 0')
  const transNlpMax = dv.getFloat64(8, true), transBetaMax = dv.getFloat64(16, true)
  const permNlpMax = dv.getFloat64(24, true), seMax = dv.getFloat64(32, true), slopeMax = dv.getFloat64(40, true)
  for (const [name, s] of [['trans_nlp_max', transNlpMax], ['trans_beta_max', transBetaMax], ['perm_nlp_max', permNlpMax],
    ['se_max', seMax], ['slope_max', slopeMax]] as [string, number][])
    if (!(Number.isFinite(s) && s >= 0)) throw bad(`${name} ${s} is not finite and at least 0`)

  if (48 + 2 * V > p.length) throw bad('the count column runs past the payload')
  const counts = new Uint16Array(p.buffer, p.byteOffset + 48, V)
  const rowStart = new Uint32Array(V + 1)
  for (let i = 0; i < V; i++) rowStart[i + 1] = rowStart[i] + counts[i]
  const R = rowStart[V]
  const A = pad4(48 + 2 * V), B = pad4(A + 2 * R)
  if (p.length !== B + 20 * R) throw bad(`payload is ${p.length} bytes, the layout needs ${B + 20 * R}`)
  for (let i = 48 + 2 * V; i < A; i++) if (p[i]) throw bad('padding after the count column is not zero')
  for (let i = A + 2 * R; i < B; i++) if (p[i]) throw bad('padding after the flags column is not zero')

  const kind = p.subarray(A, A + R), flags = p.subarray(A + R, A + 2 * R)
  const at = (o: number) => p.byteOffset + o
  const ord = new Uint16Array(p.buffer, at(B), R)
  const v1 = new Uint16Array(p.buffer, at(B + 2 * R), R)
  const v2 = new Uint16Array(p.buffer, at(B + 4 * R), R)
  const v3 = new Int16Array(p.buffer, at(B + 6 * R), R)
  const intronStart = new Uint32Array(p.buffer, at(B + 8 * R), R)
  const intronEnd = new Uint32Array(p.buffer, at(B + 12 * R), R)
  const cluster = new Uint32Array(p.buffer, at(B + 16 * R), R)

  let maxNlp = -1, maxBeta = -1, maxPermNlp = -1, maxSe = -1, maxSlope = -1
  for (let r = 0; r < R; r++) {
    const k = kind[r], f = flags[r]
    if (k > 5) throw bad(`row ${r}: kind ${k} is above 5`)
    if (f & 0xfc) throw bad(`row ${r}: reserved flag bits are set`)
    if ((f & 1) && !(k & 1)) throw bad(`row ${r}: flags bit 0 (strand) on even kind ${k}`)
    if ((f & 2) && k !== 2 && k !== 3) throw bad(`row ${r}: flags bit 1 (significant) on kind ${k}`)
    if (v3[r] === -32768) throw bad(`row ${r}: signed code -32768`)
    // kinds 0-3 carry an nlp code in v1; kinds 4-5 carry a pip code, which has no such limit
    if (k <= 3 && v1[r] > NLP_MAXQ) throw bad(`row ${r}: -log10 p code ${v1[r]} is above ${NLP_MAXQ}`)
    if (k <= 1) {
      if (v2[r]) throw bad(`row ${r}: v2 is nonzero on a trans row`)
      if (v1[r] > maxNlp) maxNlp = v1[r]
      const a = Math.abs(v3[r]); if (a > maxBeta) maxBeta = a
    } else if (k <= 3) {
      if (v1[r] > maxPermNlp) maxPermNlp = v1[r]
      if (v2[r] > maxSe) maxSe = v2[r]
      const a = Math.abs(v3[r]); if (a > maxSlope) maxSlope = a
    } else {
      if (v3[r]) throw bad(`row ${r}: v3 is nonzero on a credible-set row`)
    }
    if (!(k & 1) && (intronStart[r] || intronEnd[r] || cluster[r])) throw bad(`row ${r}: an intron field is nonzero on even kind ${k}`)
  }
  const scale = (name: string, max: number, top: number, s: number) => {
    if (max < 0) return                                  // no row of those kinds
    if (s > 0 ? max !== top : max !== 0) throw bad(`${name} ${s} breaks the scale rule (largest code ${max}, expected ${s > 0 ? top : 0})`)
  }
  scale('trans_nlp_max', maxNlp, NLP_MAXQ, transNlpMax)
  scale('trans_beta_max', maxBeta, BETA_MAXQ, transBetaMax)
  scale('perm_nlp_max', maxPermNlp, NLP_MAXQ, permNlpMax)
  scale('se_max', maxSe, 65535, seMax)
  scale('slope_max', maxSlope, BETA_MAXQ, slopeMax)

  for (let i = 0; i < V; i++) {
    for (let r = rowStart[i] + 1; r < rowStart[i + 1]; r++) {
      if (kind[r] < kind[r - 1]) throw bad(`variant slot ${i}: kind decreases at row ${r}`)
      if (kind[r] === kind[r - 1] && v1[r] > v1[r - 1]) throw bad(`variant slot ${i}: v1 increases inside kind ${kind[r]} at row ${r}`)
    }
  }
  return { firstVidx, nVariants: V, rowStart, trans: { nlpMax: transNlpMax, betaMax: transBetaMax },
    lead: { nlpMax: permNlpMax, seMax, slopeMax }, kind, flags, ord, v1, v2, v3, intronStart, intronEnd, cluster }
}

/** SPEC section 14 step 3: binary-search one rsID block. Null when the block does not hold `rs`. */
export function rsidInBlock(block: Uint8Array, rs: number, what = 'rsid block'): { chrOrdinal: number; vidx: number } | null {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  if (block.length === 0 || block.length % 8) throw bad(`${block.length} bytes is not a whole number of 8-byte records`)
  const b = aligned(block)
  const a = new Uint32Array(b.buffer, b.byteOffset, block.length / 4)
  let lo = 0, hi = block.length / 8 - 1
  while (lo <= hi) {
    const mid = (lo + hi) >> 1, got = a[2 * mid]
    if (got === rs) {
      const ref = a[2 * mid + 1]
      return { chrOrdinal: ref >>> 27, vidx: ref & 0x7ffffff }
    }
    if (got < rs) lo = mid + 1; else hi = mid - 1
  }
  return null
}

/** SPEC section 7's cis scan: one covering block of an eQTL or sQTL span, read through its header,
 *  the one raw pair of `row`, and that row's credible-set record. Nothing else is decompressed. */
export function scanBlockRow(block: Uint8Array, row: number, dof: number, what = 'span block'):
  { anchor: number; pval: number; se: number; slope: number; pip: number | null; csId: number | null } {
  const bad = (rule: string) => new PackError(`${what}: ${rule}`)
  const b = aligned(block)
  if (b.length < 64) throw bad('shorter than the 64-byte block header')
  const dv = view(b)
  if (dv.getUint32(0, true) !== MAGIC_BLOCK) throw bad('magic is not QGB0')
  const blkLen = dv.getUint32(4, true), n = dv.getUint32(8, true)
  const anchor = dv.getInt32(16, true), nCs = dv.getUint32(28, true)
  const nlpMax = dv.getFloat64(32, true), lseMin = dv.getFloat64(40, true), lseMax = dv.getFloat64(48, true)
  if (blkLen !== b.length) throw bad(`length field ${blkLen} != slice ${b.length}`)
  if (!(row >= 0 && row < n)) throw bad(`row ${row} is not below n_rows ${n}`)
  if (!(Number.isFinite(nlpMax) && nlpMax >= 0 && Number.isFinite(lseMin) && Number.isFinite(lseMax) && lseMin <= lseMax))
    throw bad(`scales nlp_max ${nlpMax}, lse_min ${lseMin}, lse_max ${lseMax} are not valid`)
  if (64 + 4 * n + 12 * nCs > blkLen) throw bad('rows and credible sets run past the block')

  const q = dv.getUint16(64 + 4 * row, true), s = dv.getUint16(64 + 4 * row + 2, true)
  if (q > NLP_MAXQ && q !== NLP_ZERO && q !== NLP_NULL) throw bad(`-log10 p code ${q} is reserved`)
  const nlp = q === NLP_NULL ? NaN : q === NLP_ZERO ? Infinity : (q * nlpMax) / NLP_MAXQ
  const pval = q === NLP_NULL ? NaN : q === NLP_ZERO ? 0 : Math.pow(10, -nlp)
  let se = NaN, slope = NaN
  if (s !== SE_NULL) {
    const lq = s & SE_QMASK
    if (lq === SE_QMASK) throw bad('SE code 0x7FFF is invalid')
    se = Math.exp(lseMin + lq * ((lseMax - lseMin) / SE_MAXQ))
    if (q < NLP_ZERO) slope = ((s & SE_SIGN) ? -1 : 1) * se * tFromNlp(nlp, dof)
  }
  // the block's credible-set records are ascending by (row, cs_id); take this row's best PIP
  let pip: number | null = null, csId: number | null = null
  for (let k = 0, o = 64 + 4 * n; k < nCs; k++, o += 12) {
    const r = dv.getUint32(o, true)
    if (r < row) continue
    if (r > row) break
    const v = dv.getFloat32(o + 4, true), id = b[o + 8]
    if (!(v >= 0 && v <= 1) || id > 127) throw bad('credible-set pip outside [0, 1] or cs_id above 127')
    if (pip === null || v > pip) { pip = v; csId = id }
  }
  return { anchor, pval, se, slope, pip, csId }
}
