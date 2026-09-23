/**
 * The variant page's "Scan cis windows" button: every nominal cis result at one variant, from the
 * qtlstore (SPEC.md section 9's rule: a phenotype covers a site when `var_start <= vidx <
 * var_start + n_var`). Never on page load, because a span can reach several MB.
 *
 * The covering phenotypes come from the search index alone, with no request. For each phenotype
 * type the scan reads one span of its results file, from the first covering block to the end of
 * the last; both spans go out together. Only block headers are read, plus each covering block's
 * one raw pair at the variant's row and that row's credible-set record; nothing is decompressed.
 */
import { lit, rows, type Row } from './db'
import { EQTL_TYPE, fetchDecoded, getStore, resultsFile, SQTL_TYPE } from './store'
import { blockHeader, PackError, scanBlockRow, type VariantRecord } from './store-decode'

/** One nominal result at the variant. `phenotype_id` and `is_sqtl` are set on splicing rows only. */
export interface CisHit {
  gene_id: string; symbol: string | null; phenotype_id?: string; is_sqtl?: boolean
  tss_distance: number; pval_nominal: number; slope: number; slope_se: number; af: number
  pip: number | null; cs_id: number | null
}

/** A block inside a span, named by the phenotype it belongs to. */
export interface SpanBlock { off: number; len: number; ref: PhenotypeRef }
export interface Span { file: string; off: number; len: number; blocks: SpanBlock[] }

interface PhenotypeRef { gene_id: string; symbol: string | null; tss: number | null; phenotype_id: string; significant: boolean }

/** What the button knows before any request: the span per QTL type and their total size. */
export interface ScanPlan {
  chr: string; vidx: number; position: number
  e: Span | null
  s: Span | null
  bytes: number
  /** kept for the button label; v1 knows both spans before any request, so it is always false */
  sPending: boolean
}

interface Covering extends Row {
  phenotype_id: string; gene_id: string | null; symbol: string | null; tss: number | null; significant: boolean
  blk_off: number; blk_len: number
}

/** Yield to the event loop so a long span does not freeze the page. */
const YIELD_EVERY = 256
const yieldToPage = () => new Promise<void>(r => setTimeout(r, 0))

async function planSpan(v: VariantRecord, phenotypeType: string): Promise<Span | null> {
  const s = await getStore()
  if (!s.results.get(phenotypeType)?.files[v.chr]) return null
  const cover = await rows<Covering>(`SELECT p.phenotype_id, p.gene_id, g.name AS symbol, g.tss, p.significant, p.blk_off, p.blk_len
    FROM phenotypes p LEFT JOIN genes g USING (gene_id)
    WHERE p.phenotype_type = ${lit(phenotypeType)} AND p.chr = ${lit(v.chr)}
      AND p.var_start <= ${v.vidx} AND ${v.vidx} < p.var_start + p.n_var
    ORDER BY p.blk_off`)
  if (!cover.length) return null
  const off = cover[0].blk_off
  const len = cover[cover.length - 1].blk_off + cover[cover.length - 1].blk_len - off
  return { file: resultsFile(s, phenotypeType, v.chr).name, off, len, blocks: cover.map(c => ({ off: c.blk_off, len: c.blk_len,
    ref: { gene_id: c.gene_id ?? c.phenotype_id, symbol: c.symbol, tss: c.tss, phenotype_id: c.phenotype_id, significant: c.significant } })) }
}

/** The spans for a variant, from the search index alone: no network. */
export async function planScan(v: VariantRecord): Promise<ScanPlan> {
  const [e, s] = await Promise.all([planSpan(v, EQTL_TYPE), planSpan(v, SQTL_TYPE)])
  return { chr: v.chr, vidx: v.vidx, position: v.position, e, s, bytes: (e?.len ?? 0) + (s?.len ?? 0), sPending: false }
}

/** Walk one span, keeping the row of every named block that covers `vidx`. */
async function walkSpan(span: Uint8Array, s: Span, vidx: number, dof: number | null,
  make: (ref: PhenotypeRef, row: ReturnType<typeof scanBlockRow>) => CisHit): Promise<CisHit[]> {
  const named = new Map(s.blocks.map(b => [b.off, b.ref]))
  const out: CisHit[] = []
  let at = 0, seen = 0
  while (at < span.length) {
    const h = blockHeader(span, at, s.file)
    if (h.blkLen < 64 || at + h.blkLen > span.length) throw new PackError(`${s.file}: the block at byte ${at} claims ${h.blkLen} bytes`)
    const ref = named.get(s.off + at)
    if (ref && h.nRows > 0 && vidx >= h.varStart && vidx < h.varStart + h.nRows)
      out.push(make(ref, scanBlockRow(span.subarray(at, at + h.blkLen), vidx - h.varStart, dof, `${s.file} block at ${s.off + at}`)))
    at += h.blkLen
    if (++seen % YIELD_EVERY === 0) await yieldToPage()
  }
  if (out.length !== s.blocks.length) throw new PackError(`${s.file}: ${out.length} of the ${s.blocks.length} covering blocks cover vidx ${vidx}`)
  return out
}

async function scanType(plan: ScanPlan, span: Span | null, phenotypeType: string, v: VariantRecord, splicing: boolean): Promise<CisHit[]> {
  if (!span) return []
  const st = await getStore()
  const f = resultsFile(st, phenotypeType, plan.chr)
  const bytes = await fetchDecoded(splicing ? 'store:scan-sqtl' : 'store:scan-eqtl', { chr: plan.chr, vidx: plan.vidx },
    f.name, f.expect, span.off, span.len, b => b, splicing ? 'store:scan-sqtl-decode' : 'store:scan-eqtl-decode')
  return walkSpan(bytes, span, plan.vidx, f.dof, (ref, r) => {
    const base = { gene_id: ref.gene_id, symbol: ref.symbol, tss_distance: ref.tss == null ? NaN : v.position - ref.tss,
      pval_nominal: r.pval, slope: r.slope, slope_se: r.se, af: v.af, pip: r.pip, cs_id: r.csId }
    return splicing ? { ...base, phenotype_id: ref.phenotype_id, is_sqtl: ref.significant } : base
  })
}

/** Run the plan: both spans in parallel. */
export async function runScan(plan: ScanPlan, v: VariantRecord): Promise<{ e: CisHit[]; s: CisHit[] }> {
  const [e, s] = await Promise.all([scanType(plan, plan.e, EQTL_TYPE, v, false), scanType(plan, plan.s, SQTL_TYPE, v, true)])
  const byP = (a: CisHit, b: CisHit) => (Number.isNaN(a.pval_nominal) ? Infinity : a.pval_nominal) - (Number.isNaN(b.pval_nominal) ? Infinity : b.pval_nominal)
  return { e: e.sort(byP), s: s.sort(byP) }
}
