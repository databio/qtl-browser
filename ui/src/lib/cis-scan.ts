/**
 * The variant page's "Scan cis windows" button: every nominal cis result at one variant, from the
 * packs (SPEC.md section 7, "The cis scan"). Two range requests, never on page load, because a span
 * can reach about 35 MB.
 *
 * The scan reads block headers only. A block that covers the variant gives up its one raw pair at a
 * fixed offset and that row's credible-set record; nothing else in a span is decompressed, except
 * each window gene's details, which name the intron blocks the sQTL span is built from. That makes
 * this case (ii) of the plan's step 1.2: the sQTL span cannot be planned until the eQTL span has
 * arrived, so the button reports the eQTL bytes "+ splicing".
 *
 * **The window rule is `w_lo <= position <= w_hi`, and the details of every window gene are read,
 * not only a covering block's.** A gene that is sQTL-tested but not eQTL-tested has an empty block
 * (`n_rows` 0) that can never cover a variant; skipping its details loses its introns. At
 * rs10824026 that is 4 of 46 genes and 33 of 239 introns. `pipeline/steps_pack_variant.py`'s
 * `cis_scan` is the same algorithm, checked against the raw Zenodo files.
 */
import { lit, rows, type Row } from './db'
import { fetchDecoded, packFile, packManifest } from './pack'
import { decodeDetails, PackError, scanBlockRow, type VariantRecord } from './pack-decode'

/** One nominal result at the variant. `phenotype_id` and `is_sqtl` are set on splicing rows only. */
export interface CisHit {
  gene_id: string; symbol: string | null; phenotype_id?: string; is_sqtl?: boolean
  tss_distance: number; pval_nominal: number; slope: number; slope_se: number; af: number
  pip: number | null; cs_id: number | null
}

/** A block inside a span, named by the phenotype it belongs to. `ref` carries whatever the caller
 *  needs to build the hit: a gene row for the eQTL span, a splice entry for the sQTL span. */
export interface SpanBlock { off: number; len: number; ref: GeneRef | IntronRef }
export interface Span { file: string; off: number; len: number; blocks: SpanBlock[] }

interface GeneRef { kind: 'gene'; gene_id: string; symbol: string | null }
interface IntronRef { kind: 'intron'; gene_id: string; symbol: string | null; phenotype_id: string; is_sqtl: boolean }

/** What the button knows before any request: the eQTL span, and that splicing follows it. */
export interface ScanPlan {
  chr: string; vidx: number; position: number
  e: Span | null
  s: Span | null
  bytes: number
  /** the sQTL span is not known until the eQTL span's details are decoded (case (ii)) */
  sPending: boolean
}

interface WindowGene extends Row { gene_id: string; symbol: string | null; blk_off: number; blk_len: number }

/** Yield to the event loop so a long span does not freeze the page. */
const YIELD_EVERY = 256
const yieldToPage = () => new Promise<void>(r => setTimeout(r, 0))

/**
 * The eQTL span for a variant, from `search_index` alone: no network. Every gene whose window
 * covers the position has a block, and the span runs from the first of those blocks to the end of
 * the last. Blocks the span crosses that belong to no window gene are read and skipped.
 */
export async function planScan(v: VariantRecord): Promise<ScanPlan> {
  const m = await packManifest()
  const genes = await rows<WindowGene>(`SELECT gene_id, symbol, blk_off, blk_len FROM search_index
    WHERE chr = ${lit(v.chr)} AND blk_off IS NOT NULL AND w_lo <= ${v.position} AND ${v.position} <= w_hi
    ORDER BY blk_off`)
  const base = { chr: v.chr, vidx: v.vidx, position: v.position, s: null, sPending: true }
  if (!genes.length) return { ...base, e: null, bytes: 0, sPending: false }
  const off = genes[0].blk_off
  const len = genes[genes.length - 1].blk_off + genes[genes.length - 1].blk_len - off
  const blocks = genes.map(g => ({ off: g.blk_off, len: g.blk_len,
    ref: { kind: 'gene', gene_id: g.gene_id, symbol: g.symbol } as GeneRef }))
  return { ...base, e: { file: packFile(m.files.eqtl, v.chr, 'eqtl'), off, len, blocks }, bytes: len }
}

/** The 64-byte block header fields the scan walks by (SPEC section 5). */
function header(span: Uint8Array, at: number, what: string) {
  if (at + 64 > span.length) throw new PackError(`${what}: a block header at byte ${at} runs past the span`)
  const dv = new DataView(span.buffer, span.byteOffset + at, 64)
  if (dv.getUint32(0, true) !== 0x30424751) throw new PackError(`${what}: the block at byte ${at} does not start with QGB0`)
  return { blkLen: dv.getUint32(4, true), nRows: dv.getUint32(8, true), varStart: dv.getUint32(12, true) }
}

/**
 * Walk one span, keeping only the rows of blocks that cover `vidx` and are named. `details` is
 * called for every named block whether it covers or not, which is how the eQTL span names the
 * intron blocks of the sQTL span.
 */
async function walkSpan(span: Uint8Array, s: Span, vidx: number, dof: number,
  make: (ref: SpanBlock['ref'], row: ReturnType<typeof scanBlockRow>) => CisHit,
  details?: (ref: SpanBlock['ref'], block: Uint8Array) => void): Promise<CisHit[]> {
  const named = new Map(s.blocks.map(b => [b.off, b.ref]))
  const out: CisHit[] = []
  let at = 0, seen = 0
  while (at < span.length) {
    const h = header(span, at, s.file)
    if (h.blkLen < 64 || at + h.blkLen > span.length) throw new PackError(`${s.file}: the block at byte ${at} claims ${h.blkLen} bytes`)
    const ref = named.get(s.off + at)
    if (ref) {
      const block = span.subarray(at, at + h.blkLen)
      if (details) details(ref, block)
      if (h.nRows > 0 && vidx >= h.varStart && vidx < h.varStart + h.nRows)
        out.push(make(ref, scanBlockRow(block, vidx - h.varStart, dof, `${s.file} block at ${s.off + at}`)))
    }
    at += h.blkLen
    if (++seen % YIELD_EVERY === 0) await yieldToPage()
  }
  return out
}

/** One span as a single range request, timed as `mark`. */
function fetchSpan(mark: string, s: Span, detail: Record<string, unknown>): Promise<Uint8Array> {
  return fetchDecoded(mark, detail, s.file, s.off, s.len, bytes => bytes, `${mark}-decode`)
}

/** Run the plan: the eQTL span, then the sQTL span its details name. */
export async function runScan(plan: ScanPlan, v: VariantRecord): Promise<{ e: CisHit[]; s: CisHit[] }> {
  if (!plan.e) return { e: [], s: [] }
  const m = await packManifest()
  const detail = { chr: plan.chr, vidx: plan.vidx }
  const eBytes = await fetchSpan('pack:scan-eqtl', plan.e, detail)

  // every window gene's details name its introns, with each intron's block in the sQTL pack
  const introns: SpanBlock[] = []
  const eHits = await walkSpan(eBytes, plan.e, plan.vidx, m.dof.eqtl,
    (ref, r) => hit(ref, r, v),
    (ref, block) => {
      const g = ref as GeneRef
      for (const sp of spliceEntries(block, g.gene_id)) {
        introns.push({ off: sp.blk_off, len: sp.blk_len,
          ref: { kind: 'intron', gene_id: g.gene_id, symbol: g.symbol, phenotype_id: sp.phenotype_id, is_sqtl: sp.is_sqtl } })
      }
    })

  let sHits: CisHit[] = []
  if (introns.length) {
    const off = Math.min(...introns.map(b => b.off))
    const end = Math.max(...introns.map(b => b.off + b.len))
    const sSpan: Span = { file: packFile(m.files.sqtl, plan.chr, 'sqtl'), off, len: end - off, blocks: introns }
    plan.s = sSpan
    plan.sPending = false
    const sBytes = await fetchSpan('pack:scan-sqtl', sSpan, detail)
    sHits = await walkSpan(sBytes, sSpan, plan.vidx, m.dof.sqtl, (ref, r) => hit(ref, r, v))
  }
  const byP = (a: CisHit, b: CisHit) => (a.pval_nominal ?? Infinity) - (b.pval_nominal ?? Infinity)
  return { e: eHits.sort(byP), s: sHits.sort(byP) }
}

function hit(ref: SpanBlock['ref'], r: ReturnType<typeof scanBlockRow>, v: VariantRecord): CisHit {
  const base = { gene_id: ref.gene_id, symbol: ref.symbol, tss_distance: v.position - r.anchor,
    pval_nominal: r.pval, slope: r.slope, slope_se: r.se, af: v.af, pip: r.pip, cs_id: r.csId }
  return ref.kind === 'intron' ? { ...base, phenotype_id: ref.phenotype_id, is_sqtl: ref.is_sqtl } : base
}

/** The `splice` entries of a gene's details, decompressed straight out of the span. */
function spliceEntries(block: Uint8Array, geneId: string): { phenotype_id: string; is_sqtl: boolean; blk_off: number; blk_len: number }[] {
  const dv = new DataView(block.buffer, block.byteOffset, 64)
  const n = dv.getUint32(8, true), nCs = dv.getUint32(28, true)
  const dz = dv.getUint32(56, true), dl = dv.getUint32(60, true)
  if (!dz) return []                      // a kind 3 block has no details; a kind 2 block always does
  const start = 64 + 4 * n + 12 * nCs
  const parsed = decodeDetails(block.subarray(start, start + dz), dl, geneId)
  const splice = (parsed as { splice?: unknown }).splice
  if (!Array.isArray(splice)) throw new PackError(`${geneId}: details have no splice list`)
  return splice as { phenotype_id: string; is_sqtl: boolean; blk_off: number; blk_len: number }[]
}
