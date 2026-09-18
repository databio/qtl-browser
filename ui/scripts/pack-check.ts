/**
 * Round trip of the browser decoder (src/lib/pack-decode.ts) against the pipeline's reference files.
 * Run `uv run python -m pipeline validate` first; it writes data/derived/_tmp/pack_check/.
 *
 *   npm run pack-check
 *
 * - eQTL: every gene of <chr>_index.json through its block, its variants range (decodeVariantRange,
 *   sliceRun), and its details, against <chr>_rows.arrow and <chr>_details.json.
 * - sQTL: every intron of sqtl/introns.json through its kind 3 block and its gene's variants range,
 *   against sqtl/rows.arrow; its gene's details must point at the block.
 * - GWAS: every window of gwas/windows.json through the window rule and the block decoder, against
 *   gwas/rows.arrow.
 * - trans: every frame of trans/frames.json through decodeTransFrame and transPhenotypeIds, against
 *   trans/rows.arrow with the tolerances in frames.json; the decoder's reader rules on broken copies of
 *   one frame; and the trans-only variant pages of trans/variants.json against trans/variant_rows.arrow.
 * - tFromNlp against scipy (t_grid.arrow).
 * The QTL rows use validate's limits (tolerances.json, the per-row bounds). Prints the worst errors
 * and the decode times of the largest intron and the densest GWAS window; exits 1 on any failure.
 */
import { closeSync, fstatSync, openSync, readdirSync, readFileSync, readSync } from 'node:fs'
import { dirname, join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { constants as zlibConstants, zstdCompressSync, zstdDecompressSync } from 'node:zlib'
import { tableFromIPC, type Table } from 'apache-arrow'
import { csMembers, decodeGwasRange, decodeHitsFrame, decodeResultBlock, decodeTransFrame, decodeVariantIndex, decodeVariantRange,
  gwasRange, PackError, parseGwasIndex, readerColumns, rsidInBlock, scanBlockRow,
  sliceRun, tFromNlp, TRANS_CHROMS, transPhenotypeIds, VIDX_CHROMS,
  type ReaderColumns, type ResultBlock, type VariantRange, type VariantRun } from '../src/lib/pack-decode.ts'

const HERE = dirname(fileURLToPath(import.meta.url))
const DERIVED = join(HERE, '..', '..', 'data', 'derived')
const CHECK = join(DERIVED, '_tmp', 'pack_check')
const tol = JSON.parse(readFileSync(join(CHECK, 'tolerances.json'), 'utf8'))
const packs = JSON.parse(readFileSync(join(DERIVED, 'manifest.json'), 'utf8')).packs
const failures: string[] = []
const fail = (msg: string) => { failures.push(msg); console.log(`FAIL ${msg}`) }
const pass = (ok: boolean, msg: string) => (ok ? console.log(`PASS ${msg}`) : fail(msg))
const g = (x: number) => x.toPrecision(3)
const errMsg = (e: unknown) => (e instanceof Error ? e.message : String(e))

function col(t: Table, name: string): Float64Array {
  const v = t.getChild(name)
  if (!v) throw new Error(`reference rows lack ${name}`)
  const out = new Float64Array(v.length)
  for (let i = 0; i < v.length; i++) out[i] = v.isValid(i) ? Number(v.get(i)) : NaN
  return out
}
const strCol = (t: Table, name: string) => Array.from({ length: t.numRows }, (_, i) => String(t.getChild(name)!.get(i)))

/** Byte ranges of the derived files, without reading whole packs into memory. */
const fds = new Map<string, number>()
function readRange(rel: string, off: number, len: number): Uint8Array {
  let fd = fds.get(rel)
  if (fd === undefined) { fd = openSync(join(DERIVED, rel), 'r'); fds.set(rel, fd) }
  const out = new Uint8Array(len)
  for (let got = 0; got < len;) {
    const r = readSync(fd, out, got, len - got, off + got)
    if (r === 0) throw new Error(`${rel}: file ends before byte ${off + len}`)
    got += r
  }
  return out
}

/** First run, then the median of five more. */
function timeIt(f: () => void): { first: number; median: number } {
  let t = performance.now()
  f()
  const first = performance.now() - t
  const runs: number[] = []
  for (let i = 0; i < 5; i++) { t = performance.now(); f(); runs.push(performance.now() - t) }
  runs.sort((a, b) => a - b)
  return { first, median: runs[2] }
}
const ms = (x: { first: number; median: number }) => `${x.first.toFixed(1)} ms (warm median ${x.median.toFixed(1)} ms)`

// ---- QTL rows against reference rows ---------------------------------------------------------------

interface Ref { id: string[]; S: Record<string, Float64Array>; A1: string[]; A2: string[]; seB: Float64Array; slB: Float64Array; rows: number }
function loadRef(rowsPath: string, boundsPath: string, idCol: string, label: string): Ref {
  const rows = tableFromIPC(readFileSync(rowsPath)), bounds = tableFromIPC(readFileSync(boundsPath))
  const S = Object.fromEntries(['position', 'rs_number', 'tss_distance', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se', 'pip', 'cs_id']
    .map(k => [k, col(rows, k)])) as Record<string, Float64Array>
  const id = strCol(rows, idCol), bId = strCol(bounds, idCol), bPos = col(bounds, 'position')
  pass(bounds.numRows === rows.numRows && bId.every((x, i) => x === id[i] && bPos[i] === S.position[i]),
    `${label}: bounds rows align with reference rows (${rows.numRows.toLocaleString()})`)
  return { id, S, A1: strCol(rows, 'A1'), A2: strCol(rows, 'A2'), seB: col(bounds, 'se_bound'), slB: col(bounds, 'slope_bound'), rows: rows.numRows }
}

function newStats() {
  return { exactBad: new Map<string, number>(), rows: 0, afMax: 0, nlpMax: 0, nlpRatio: 0, seMax: 0, slMax: 0, seRatio: 0, slRatio: 0, sanityBad: 0,
    bands: (tol.bands as number[]).map((lo, i, a) => ({ lo, hi: i + 1 < a.length ? a[i + 1] : Infinity, n: 0, se: 0, sl: 0, rse: 0, rsl: 0 })) }
}
type Stats = ReturnType<typeof newStats>

/** One phenotype's decoded rows against reference rows start.. (same order: vidx). */
function compareRows(st: Stats, block: ResultBlock, run: VariantRun, c: ReaderColumns, ref: Ref, start: number) {
  const { S, A1, A2 } = ref
  const bump = (k: string) => st.exactBad.set(k, (st.exactBad.get(k) ?? 0) + 1)
  for (let i = 0; i < block.nRows; i++) {
    const k = start + i
    st.rows++
    if (c.position[i] !== S.position[k] || c.a1[i] !== A1[k] || c.a2[i] !== A2[k]) bump('position, A1, A2')
    if (c.rsNumber[i] === 0 ? !Number.isNaN(S.rs_number[k]) : c.rsNumber[i] !== S.rs_number[k]) bump('rs_number')
    if (c.tssDistance[i] !== S.tss_distance[k]) bump('tss_distance')
    for (const [name, got] of [['ma_samples', c.maSamples[i]], ['ma_count', c.maCount[i]]] as const) {
      const s = S[name][k]
      if (Number.isNaN(s) ? got !== -1 : got !== s) bump(name)
    }
    const pipS = S.pip[k], csS = S.cs_id[k]
    if (Number.isNaN(pipS) ? !Number.isNaN(c.pip[i]) : c.pip[i] !== Math.fround(pipS)) bump('pip')
    if (Number.isNaN(csS) ? c.csId[i] !== -1 : c.csId[i] !== csS) bump('cs_id')
    const p = S.pval_nominal[k], sl = S.slope[k], se = S.slope_se[k]
    if (Number.isNaN(block.pval[i]) !== Number.isNaN(p)) bump('null p')
    if (Number.isNaN(block.se[i]) !== (Number.isNaN(se) || Number.isNaN(sl))) bump('null SE')
    if (Number.isNaN(block.slope[i]) !== (Number.isNaN(sl) || Number.isNaN(p) || p === 0 || Number.isNaN(se))) bump('null slope')
    if (Number.isNaN(c.slope[i]) !== Number.isNaN(block.slope[i]) || Number.isNaN(c.se[i]) !== Number.isNaN(block.se[i])) bump('reader column nulls')
    if (Number.isNaN(S.af[k]) ? !Number.isNaN(run.af[i]) : !(Math.abs(run.af[i] - S.af[k]) <= tol.af_tol)) bump('af')
    if (!Number.isNaN(S.af[k])) st.afMax = Math.max(st.afMax, Math.abs(run.af[i] - S.af[k]))
    if (!(p > 0) || Number.isNaN(se) || Number.isNaN(sl)) continue
    const nlpErr = Math.abs(block.nlp[i] + Math.log10(p)), nlpLim = block.nlpMax / tol.nlp_half_step_divisor + tol.nlp_slack
    st.nlpMax = Math.max(st.nlpMax, nlpErr); st.nlpRatio = Math.max(st.nlpRatio, nlpErr / nlpLim)
    const eSe = Math.abs(block.se[i] - se), eSl = Math.abs(block.slope[i] - sl)
    const rse = eSe / ref.seB[k], rsl = eSl / ref.slB[k]
    st.seMax = Math.max(st.seMax, eSe); st.slMax = Math.max(st.slMax, eSl); st.seRatio = Math.max(st.seRatio, rse); st.slRatio = Math.max(st.slRatio, rsl)
    if (!(Number.isFinite(block.se[i]) && block.se[i] > 0 && Number.isFinite(block.slope[i])) ||
        (sl !== 0 && (block.slope[i] < 0 || Object.is(block.slope[i], -0)) !== sl < 0)) st.sanityBad++
    const t = Math.abs(sl) / se
    const band = st.bands.find(b => t >= b.lo && t < b.hi)!
    band.n++; band.se = Math.max(band.se, eSe); band.sl = Math.max(band.sl, eSl); band.rse = Math.max(band.rse, rse); band.rsl = Math.max(band.rsl, rsl)
  }
}

function report(st: Stats, label: string) {
  const { exactBad } = st
  pass(exactBad.size === 0, `${label}: exact on ${st.rows.toLocaleString()} rows: position, alleles, rs_number, tss_distance, counts, pip, cs_id, null pattern, af within ${g(tol.af_tol)}${exactBad.size ? ` (${[...exactBad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
  pass(st.nlpRatio <= 1, `${label}: -log10 p within half the step: max error ${g(st.nlpMax)}, max error / limit ${st.nlpRatio.toFixed(4)}`)
  pass(st.seRatio <= tol.bound_factor && st.slRatio <= tol.bound_factor, `${label}: SE and slope within ${tol.bound_factor} x the per-row bound: max error / bound ${st.seRatio.toFixed(4)} (SE), ${st.slRatio.toFixed(4)} (slope); worst SE error ${g(st.seMax)}, slope ${g(st.slMax)}`)
  pass(st.sanityBad === 0, `${label}: every non-null row has a finite SE > 0 and a finite slope with the source's sign (${st.sanityBad} bad)`)
  console.log(`  worst errors: af ${g(st.afMax)}, -log10 p ${g(st.nlpMax)}, SE ${g(st.seMax)}, slope ${g(st.slMax)}`)
  for (const b of st.bands) if (b.n) console.log(`  |t| [${b.lo}, ${b.hi}): ${b.n.toLocaleString()} rows, SE max ${g(b.se)} (error / bound ${b.rse.toFixed(4)}), slope max ${g(b.sl)} (error / bound ${b.rsl.toFixed(4)})`)
}

/** Reference rows grouped by id: id -> [first, end). */
function groups(ids: string[], label: string): Map<string, [number, number]> {
  const at = new Map<string, [number, number]>()
  for (let i = 0; i < ids.length;) {
    let j = i
    while (j < ids.length && ids[j] === ids[i]) j++
    if (at.has(ids[i])) fail(`${label}: reference rows of ${ids[i]} are not contiguous`)
    at.set(ids[i], [i, j])
    i = j
  }
  return at
}

// ---- eQTL -------------------------------------------------------------------------------------------

const P = packs.variant_page_size as number
const wider: string[] = []
let eqtlGenes = 0, emptyBlocks = 0
for (const f of readdirSync(CHECK).filter(x => x.endsWith('_index.json')).sort()) {
  const chrom = f.replace('_index.json', '')
  const index = JSON.parse(readFileSync(join(CHECK, f), 'utf8'))
  const details = JSON.parse(readFileSync(join(CHECK, `${chrom}_details.json`), 'utf8'))
  const ref = loadRef(join(CHECK, `${chrom}_rows.arrow`), join(CHECK, `${chrom}_bounds.arrow`), 'gene_id', chrom)
  const at = groups(ref.id, chrom)
  const st = newStats()
  let detailBad = 0, splicePointers = 0, covered = 0
  for (const row of index.rows) {
    eqtlGenes++
    try {
      const block = decodeResultBlock(readRange(index.files.eqtl, row.blk_off, row.blk_len), 2, index.dof, row, `${index.files.eqtl} ${row.symbol}`)
      // <chr>_details.json holds each gene's details object ({ v, gene, exons, splice }) as validate rebuilt it
      if (JSON.stringify(block.details) !== JSON.stringify(details[row.gene_id])) detailBad++
      splicePointers += ((block.details!.splice as { blk_off?: number; blk_len?: number }[])).filter(s => Number.isInteger(s.blk_off) && s.blk_len! > 0).length
      const range = decodeVariantRange(readRange(index.files.variants, row.var_off, row.var_len), row, `${index.files.variants} ${row.symbol}`)
      const [a, b] = at.get(row.gene_id) ?? [0, 0]
      covered += b - a
      const n = block.nRows
      if (n !== b - a) { fail(`${row.symbol}: block has ${n} rows, reference ${b - a}`); continue }
      if (!n) { emptyBlocks++; continue }
      // a union range reaches past the pages of the eQTL run when an intron run does
      if (range.first < Math.floor(row.var_start / P) * P || range.first + range.n > (Math.floor((row.var_start + n - 1) / P) + 1) * P) wider.push(row.symbol)
      const run = sliceRun(range, block.varStart!, n, block, `${index.files.variants} ${row.symbol}`)
      compareRows(st, block, run, readerColumns(block, run), ref, a)
    } catch (e) { fail(`${chrom} ${row.symbol}: ${errMsg(e)}`) }
  }
  pass(covered === ref.rows, `${chrom}: every reference row belongs to a sample gene (${covered} of ${ref.rows})`)
  pass(detailBad === 0, `${chrom}: details JSON (with ${splicePointers} splice blk_off/blk_len pointers) equals the reference for ${index.rows.length - detailBad} of ${index.rows.length} genes`)
  report(st, chrom)
}
pass(wider.includes('TBC1D5'), `eQTL: ${eqtlGenes} sample genes decode, ${emptyBlocks} sQTL-only with an empty block; ${wider.length} with a variants range wider than their eQTL run's pages (${wider.join(', ')}), TBC1D5 among them`)

// ---- sQTL -------------------------------------------------------------------------------------------

{
  const SQ = join(CHECK, 'sqtl')
  const info = JSON.parse(readFileSync(join(SQ, 'introns.json'), 'utf8'))
  pass(info.dof === packs.dof.sqtl, `sQTL: reference dof ${info.dof} equals manifest packs.dof.sqtl ${packs.dof.sqtl}`)
  const ref = loadRef(join(SQ, 'rows.arrow'), join(SQ, 'bounds.arrow'), 'phenotype_id', 'sQTL')
  const at = groups(ref.id, 'sQTL')
  const genes = new Map<string, Record<string, number | string | null>>(info.genes.map((x: { gene_id: string }) => [x.gene_id, x]))
  const ranges = new Map<string, VariantRange>()
  const splice = new Map<string, { phenotype_id: string; blk_off: number; blk_len: number }[]>()
  const st = newStats()
  let covered = 0, pointerBad = 0, twoSetRows = 0
  const twoSetIntrons: string[] = []
  for (const it of info.introns) {
    try {
      const gene = genes.get(it.gene_id)!
      const files = info.files[it.chrom]
      const label = `${gene.symbol} ${it.phenotype_id}`
      if (!splice.has(it.gene_id)) {
        const eq = decodeResultBlock(readRange(packs.files.eqtl[it.chrom], gene.blk_off as number, gene.blk_len as number), 2, packs.dof.eqtl,
          gene as { blk_len: number }, `${packs.files.eqtl[it.chrom]} ${gene.symbol}`)
        splice.set(it.gene_id, eq.details!.splice as { phenotype_id: string; blk_off: number; blk_len: number }[])
        ranges.set(it.gene_id, decodeVariantRange(readRange(files.variants, gene.var_off as number, gene.var_len as number),
          gene as { var_len: number }, `${files.variants} ${gene.symbol}`))
      }
      const s = splice.get(it.gene_id)!.find(x => x.phenotype_id === it.phenotype_id)
      if (!s || s.blk_off !== it.blk_off || s.blk_len !== it.blk_len) pointerBad++
      const block = decodeResultBlock(readRange(files.sqtl, it.blk_off, it.blk_len), 3, info.dof, { blk_len: it.blk_len }, `${files.sqtl} ${label}`)
      const [a, b] = at.get(it.phenotype_id) ?? [0, 0]
      covered += b - a
      if (block.nRows !== b - a) { fail(`${label}: block has ${block.nRows} rows, reference ${b - a}`); continue }
      const run = sliceRun(ranges.get(it.gene_id)!, block.varStart!, block.nRows, block, `${files.variants} ${label}`)
      compareRows(st, block, run, readerColumns(block, run), ref, a)
      const perRow = new Map<number, number>()
      for (const m of csMembers(block)) perRow.set(m.row, (perRow.get(m.row) ?? 0) + 1)
      const doubles = [...perRow.values()].filter(x => x > 1).length
      if (doubles) { twoSetRows += doubles; twoSetIntrons.push(`${label} (${doubles} rows)`) }
    } catch (e) { fail(`sQTL ${it.phenotype_id}: ${errMsg(e)}`) }
  }
  pass(covered === ref.rows, `sQTL: ${info.introns.length} introns of ${info.genes.length} genes on ${Object.keys(info.files).join(', ')} cover every reference row (${covered} of ${ref.rows})`)
  pass(pointerBad === 0, `sQTL: each intron's block equals the blk_off/blk_len in its gene's details (${pointerBad} differ)`)
  report(st, 'sQTL')
  pass(twoSetRows === 12 && twoSetIntrons.length === 1 && twoSetIntrons[0].startsWith('LINC01954'),
    `sQTL: rows with two credible-set records keep both (SPEC: 12 on one LINC01954 intron): ${twoSetIntrons.join('; ') || 'none'}`)
}

// ---- GWAS -------------------------------------------------------------------------------------------

{
  const GW = join(CHECK, 'gwas')
  const info = JSON.parse(readFileSync(join(GW, 'windows.json'), 'utf8'))
  const rows = tableFromIPC(readFileSync(join(GW, 'rows.arrow')))
  const index = parseGwasIndex(readFileSync(join(DERIVED, info.index)), info.index)
  pass(index.blockRows === info.block_rows && index.blockRows === packs.gwas_block_rows && index.nValues.join() === info.n_values.join(),
    `GWAS index: ${index.chroms.size} chromosomes, ${index.blockRows} rows per block, n table ${index.nValues.join(', ')}`)
  const win = col(rows, 'window'), pos = col(rows, 'position'), rs = col(rows, 'rs_number'), nn = col(rows, 'n')
  const beta = col(rows, 'beta'), se = col(rows, 'se'), eaf = col(rows, 'eaf'), pv = col(rows, 'p')
  const ea = strCol(rows, 'ea'), nea = strCol(rows, 'nea')
  const at = groups(Array.from(win, String), 'GWAS')
  let rangeBad = 0, decoded = 0, pRel = 0, withRange = 0
  const bad = new Map<string, number>(), exact = { beta: 0, se: 0, eaf: 0, p: 0 }
  const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
  for (const w of info.windows) {
    const label = `window ${w.id} (${w.name}) ${w.chrom}:${w.lo}-${w.hi}`
    try {
      const r = gwasRange(index, w.chrom, w.lo, w.hi)
      if ((r === null) !== (w.byte_start === null) || (r && (r.off !== w.byte_start || r.off + r.len !== w.byte_end))) {
        rangeBad++
        fail(`${label}: window rule gives ${r ? `${r.off}..${r.off + r.len}` : 'no range'}, validate ${w.byte_start}..${w.byte_end}`)
        continue
      }
      const [a, b] = at.get(String(w.id)) ?? [0, 0]
      if (b - a !== w.rows) fail(`${label}: ${b - a} reference rows, windows.json says ${w.rows}`)
      if (!r) continue
      withRange++
      const c = decodeGwasRange(readRange(info.files[w.chrom], r.off, r.len), index, w.chrom, r, w.lo, w.hi, `${info.files[w.chrom]} ${label}`)
      if (c.rows !== b - a) { fail(`${label}: decoded ${c.rows} rows, reference ${b - a}`); continue }
      for (let i = 0; i < c.rows; i++) {
        const k = a + i
        decoded++
        if (c.position[i] !== pos[k] || c.ea[i] !== ea[k] || c.nea[i] !== nea[k]) bump('position, ea, nea')
        if (c.rsNumber[i] !== (Number.isNaN(rs[k]) ? 0 : rs[k])) bump('rs_number')
        if (c.n[i] !== nn[k]) bump('n')
        for (const [name, got, want] of [['beta', c.beta[i], beta[k]], ['se', c.se[i], se[k]], ['eaf', c.eaf[i], eaf[k]]] as const) {
          // SPEC section 9: the codes hold the printed 4-decimal values exactly
          if (Math.round(got * 1e4) !== Math.round(want * 1e4) || Math.abs(got - want) > 1e-12 * Math.max(1, Math.abs(want))) bump(name)
          if (got === want) exact[name]++
        }
        const rel = Math.abs(c.p[i] - pv[k]) / pv[k]
        pRel = Math.max(pRel, rel)
        if (!(rel <= 1e-12)) bump('p')
        if (c.p[i] === pv[k]) exact.p++
      }
    } catch (e) { fail(`GWAS ${label}: ${errMsg(e)}`) }
  }
  const kinds = [...new Set(info.windows.map((w: { name: string }) => w.name))].join(', ')
  pass(rangeBad === 0, `GWAS: the window rule gives validate's byte range for all ${info.windows.length} windows (${kinds}); ${withRange} have a range`)
  pass(bad.size === 0 && decoded === rows.numRows, `GWAS: ${decoded.toLocaleString()} of ${rows.numRows.toLocaleString()} rows equal the source: position, alleles, rs_number, n exact; beta, se, eaf equal after 4-decimal rounding (bit-exact ${exact.beta}, ${exact.se}, ${exact.eaf}); p within 1e-12 relative (max ${g(pRel)}, bit-exact ${exact.p})${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)

  // the densest window a gene can ask for: the 2 Mb span (w_hi - w_lo <= 2,000,000) with the largest byte range
  let best = { chrom: '', lo: 0, hi: 0, len: 0 }
  for (const [chrom, c] of index.chroms) {
    for (const lo of c.firstPosition) {
      const r = gwasRange(index, chrom, lo, lo + 2_000_000)
      if (r && r.len > best.len) best = { chrom, lo, hi: lo + 2_000_000, len: r.len }
    }
  }
  const r = gwasRange(index, best.chrom, best.lo, best.hi)!
  const bytes = readRange(packs.files.gwas[best.chrom], r.off, r.len)
  let n = 0
  const t = timeIt(() => { n = decodeGwasRange(bytes, index, best.chrom, r, best.lo, best.hi).rows })
  console.log(`  Node decode of the densest 2 Mb GWAS window ${best.chrom}:${best.lo}-${best.hi} (${r.lastBlock - r.firstBlock + 1} blocks, ${r.len.toLocaleString()} B, ${n.toLocaleString()} rows): ${ms(t)}`)
}

// ---- decode time of the largest intron --------------------------------------------------------------

{
  let best = { chrom: '', off: 0, len: 0, n: 0, varStart: 0 }
  const head = new Uint8Array(64), hv = new DataView(head.buffer)
  for (const [chrom, path] of Object.entries(packs.files.sqtl as Record<string, string>)) {
    const fd = openSync(join(DERIVED, path), 'r'), size = fstatSync(fd).size
    for (let off = 32; off < size;) {
      readSync(fd, head, 0, 64, off)
      const len = hv.getUint32(4, true), n = hv.getUint32(8, true)
      if (n > best.n) best = { chrom, off, len, n, varStart: hv.getUint32(12, true) }
      if (len < 64) throw new Error(`${path}: block at ${off} has length ${len}`)
      off += len
    }
    closeSync(fd)
  }
  // the variant pages covering the intron's run
  const vpath = packs.files.variants[best.chrom] as string
  const fd = openSync(join(DERIVED, vpath), 'r'), size = fstatSync(fd).size
  const ph = new Uint8Array(12), pv = new DataView(ph.buffer)
  let vOff = -1, vEnd = -1
  for (let off = 32; off < size && vEnd < 0;) {
    readSync(fd, ph, 0, 12, off)
    const stored = pv.getUint32(0, true), first = pv.getUint32(4, true), n = pv.getUint16(8, true)
    const stop = off + Math.ceil((12 + stored) / 4) * 4
    if (vOff < 0 && best.varStart < first + n) vOff = off
    if (best.varStart + best.n - 1 < first + n) vEnd = stop
    off = stop
  }
  closeSync(fd)
  const blockBytes = readRange(packs.files.sqtl[best.chrom], best.off, best.len), varBytes = readRange(vpath, vOff, vEnd - vOff)
  let block: ResultBlock | null = null, range: VariantRange | null = null
  const tb = timeIt(() => { block = decodeResultBlock(blockBytes, 3, packs.dof.sqtl, { blk_len: best.len }) })
  const tv = timeIt(() => { range = decodeVariantRange(varBytes, { var_len: vEnd - vOff }) })
  const tc = timeIt(() => { const b = block!; readerColumns(b, sliceRun(range!, b.varStart!, b.nRows, b)) })
  console.log(`  Node decode of the largest intron (${best.chrom} sQTL block at byte ${best.off}, ${best.n.toLocaleString()} rows, ${best.len.toLocaleString()} B; SPEC A.5: PPT2): ` +
    `block ${ms(tb)}; its ${(vEnd - vOff).toLocaleString()} B of variant pages ${ms(tv)}; run + reader columns ${ms(tc)}`)
}

// ---- trans pack -------------------------------------------------------------------------------------

{
  const TR = join(CHECK, 'trans')
  const info = JSON.parse(readFileSync(join(TR, 'frames.json'), 'utf8'))
  const tl = info.tolerances
  const dof = { e: packs.dof.eqtl as number, s: packs.dof.sqtl as number }
  pass(info.dof.eqtl === dof.e && info.dof.sqtl === dof.s && JSON.stringify(info.files) === JSON.stringify(packs.files.trans),
    `trans: reference dof (${info.dof.eqtl}, ${info.dof.sqtl}) and ${Object.keys(info.files).length} files equal manifest packs.dof and packs.files.trans`)
  const rows = tableFromIPC(readFileSync(join(TR, 'rows.arrow')))
  const gid = strCol(rows, 'gene_id'), qtl = strCol(rows, 'qtl_type'), pheno = strCol(rows, 'phenotype_id')
  const vchr = strCol(rows, 'variant_chr'), gchr = strCol(rows, 'gene_chr'), rsid = rows.getChild('rsid')!
  const pos = col(rows, 'position'), rs = col(rows, 'rs_number'), af = col(rows, 'af'), pv = col(rows, 'pval')
  const be = col(rows, 'beta'), se = col(rows, 'beta_se'), r2 = col(rows, 'r2')
  const at = groups(gid, 'trans')
  const bad = new Map<string, number>()
  const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
  const worst = { af: 0, nlp: 0, nlpRatio: 0, beta: 0, betaRatio: 0, se: 0, r2: 0 }
  let covered = 0, eRows = 0, sRows = 0
  for (const gene of info.genes) {
    const label = `${gene.symbol} (${gene.reason}, ${gene.chr})`
    try {
      const f = decodeTransFrame(readRange(gene.file, gene.trans_off, gene.trans_len), dof, `${gene.file} ${gene.symbol} bytes ${gene.trans_off}+${gene.trans_len}`)
      const [a, b] = at.get(gene.gene_id) ?? [0, 0]
      const n = f.nE + f.nS
      if (f.nE !== gene.n_e || f.nS !== gene.n_s || f.k !== gene.k || f.nlpMax !== gene.nlp_max || f.betaMax !== gene.beta_max || n !== b - a) {
        fail(`${label}: frame n_e ${f.nE}, n_s ${f.nS}, k ${f.k}, nlp_max ${f.nlpMax}, beta_max ${f.betaMax}; frames.json ${gene.n_e}, ${gene.n_s}, ${gene.k}, ${gene.nlp_max}, ${gene.beta_max}; ${b - a} reference rows`)
        continue
      }
      covered += n; eRows += f.nE; sRows += f.nS
      const ids = transPhenotypeIds(f, gene)
      for (let i = 0; i < n; i++) {
        const r = a + i
        if ((i < f.nE ? 'e' : 's') !== qtl[r]) bump('qtl_type')
        if (ids[i] !== pheno[r]) bump('phenotype_id')
        if (gchr[r] !== gene.chr) bump('gene_chr')
        if (TRANS_CHROMS[f.variantChr[i]] !== vchr[r] || f.position[i] !== pos[r]) bump('variant_chr, position')
        if (f.rsNumber[i] === 0 ? !Number.isNaN(rs[r]) : f.rsNumber[i] !== rs[r]) bump('rs_number')
        if ((f.rsNumber[i] ? `rs${f.rsNumber[i]}` : null) !== (rsid.get(r) ?? null)) bump('rsid')
        const eAf = Math.abs(f.af[i] - af[r])
        worst.af = Math.max(worst.af, eAf)
        if (!(eAf <= tl.af_tol)) bump('af')
        const eN = Math.abs(f.nlp[i] + Math.log10(pv[r])), limN = f.nlpMax / tl.nlp_half_step_divisor + tl.nlp_slack
        worst.nlp = Math.max(worst.nlp, eN); worst.nlpRatio = Math.max(worst.nlpRatio, eN / limN)
        if (!(eN <= limN)) bump('-log10 p')
        const eB = Math.abs(f.beta[i] - be[r]), limB = f.betaMax / tl.beta_half_step_divisor + tl.beta_rel_f32 * Math.abs(be[r])
        worst.beta = Math.max(worst.beta, eB); worst.betaRatio = Math.max(worst.betaRatio, eB / limB)
        if (!(eB <= limB)) bump('beta')
        const eS = Math.abs(f.betaSe[i] - se[r]) / se[r]
        worst.se = Math.max(worst.se, eS)
        if (!(eS <= tl.beta_se_rel)) bump('beta_se')
        const eR = Math.abs(f.r2[i] - r2[r])
        worst.r2 = Math.max(worst.r2, eR)
        if (!(eR <= tl.r2_abs)) bump('r2')
      }
    } catch (e) { fail(`trans ${label}: ${errMsg(e)}`) }
  }
  const reasons = [...new Set(info.genes.map((x: { reason: string }) => x.reason))].join(', ')
  pass(covered === rows.numRows, `trans: ${info.genes.length} frames (${reasons}) decode to all ${rows.numRows.toLocaleString()} reference rows (${eRows.toLocaleString()} eQTL, ${sRows.toLocaleString()} sQTL): ${covered.toLocaleString()} matched`)
  pass(bad.size === 0, `trans: qtl_type, phenotype_id, gene_chr, variant_chr, position, rs_number, rsid exact; af within ${g(tl.af_tol)}; -log10 p and beta within half a step (beta plus float32); beta_se within ${tl.beta_se_rel * 100}% relative; r2 within ${tl.r2_abs}${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
  console.log(`  worst errors: af ${g(worst.af)}, -log10 p ${g(worst.nlp)} (error / limit ${worst.nlpRatio.toFixed(4)}), beta ${g(worst.beta)} (error / limit ${worst.betaRatio.toFixed(4)}), beta_se relative ${g(worst.se)}, r2 ${g(worst.r2)}`)

  const big = info.genes.reduce((x: { trans_len: number }, y: { trans_len: number }) => (y.trans_len > x.trans_len ? y : x))
  const bigBytes = readRange(big.file, big.trans_off, big.trans_len)
  let bigRows = 0
  const tt = timeIt(() => { const f = decodeTransFrame(bigBytes, dof); bigRows = f.nE + f.nS })
  console.log(`  Node decode of the largest trans frame (${big.symbol}, ${big.trans_len.toLocaleString()} B, ${bigRows.toLocaleString()} rows): ${ms(tt)}`)

  // reader rules: FLNC's frame (eQTL and sQTL rows, 3 introns) broken one rule at a time and framed again
  const flnc = info.genes.find((x: { symbol: string }) => x.symbol === 'FLNC')
  const raw = new Uint8Array(zstdDecompressSync(readRange(flnc.file, flnc.trans_off, flnc.trans_len)))
  const reframe = (p: Uint8Array, checksum: number) =>
    new Uint8Array(zstdCompressSync(p, { params: { [zlibConstants.ZSTD_c_checksumFlag]: checksum, [zlibConstants.ZSTD_c_contentSizeFlag]: 1 } }))
  const hv = new DataView(raw.buffer)
  const nE = hv.getUint32(4, true), nS = hv.getUint32(8, true), k = hv.getUint16(12, true), n = nE + nS
  const R = Math.ceil((32 + 13 * k) / 4) * 4, C = Math.ceil((R + 14 * n) / 4) * 4
  type Case = { rule: string; want: RegExp; edit?: (p: Uint8Array, d: DataView) => Uint8Array | void; checksum?: number }
  const cases: Case[] = [
    { rule: 'content checksum', want: /no content checksum/, checksum: 0 },
    { rule: 'shorter than 32', want: /fewer than 32/, edit: p => p.slice(0, 28) },
    { rule: 'magic', want: /magic is not QTT0/, edit: p => { p[3] = 0x31 } },
    { rule: 'reserved field', want: /reserved field/, edit: (_, d) => d.setUint16(14, 1, true) },
    { rule: 'no rows', want: /n_e \+ n_s is 0/, edit: (_, d) => { d.setUint32(4, 0, true); d.setUint32(8, 0, true) } },
    { rule: 'k 0 with sQTL rows', want: /k 0 with n_s/, edit: (_, d) => d.setUint16(12, 0, true) },
    { rule: 'nlp_max not finite', want: /not finite/, edit: (_, d) => d.setFloat64(16, NaN, true) },
    { rule: 'payload length', want: /lay out/, edit: p => p.slice(0, p.length - 4) },
    { rule: 'nonzero padding', want: /padding byte/, edit: p => { p[R - 1] = 1 } },
    { rule: 'strand code', want: /strand code 2/, edit: p => { p[32 + 12 * k] = 2 } },
    { rule: 'intron table order', want: /strictly ascending/, edit: (p, d) => {
      for (const o of [0, 4, 8]) d.setUint32(32 + o * k + 4, d.getUint32(32 + o * k, true), true)
      p[32 + 12 * k + 1] = p[32 + 12 * k]
    } },
    { rule: 'nlp code 65534', want: /nlp code 65534/, edit: (_, d) => d.setUint16(R + 10 * n, 65534, true) },
    { rule: 'af code 65535', want: /af code 65535/, edit: (_, d) => d.setUint16(R + 8 * n, 65535, true) },
    { rule: 'beta code -32768', want: /beta code -32768/, edit: (_, d) => d.setInt16(R + 12 * n, -32768, true) },
    { rule: 'variant_chr 24', want: /variant_chr 24/, edit: p => { p[C] = 24 } },
    { rule: 'variant_chr decreasing', want: /decreases inside a run/, edit: p => { p[C] = 23 } },
    { rule: 'position entry 0', want: /position entry is 0/, edit: (_, d) => d.setUint32(R, 0, true) },
    { rule: 'intron index k', want: /is not below k/, edit: p => { p[C + n] = k } },
    { rule: 'intron indices decreasing', want: /intron indices decrease/, edit: p => { p[C + n + nS - 1] = 0 } },
    { rule: 'intron without rows', want: /has no row/, edit: p => { p.fill(0, C + n, C + n + nS) } },
  ]
  let rejected = 0
  try { decodeTransFrame(reframe(raw, 1), dof, 'FLNC framed again') } catch (e) { fail(`trans: FLNC's payload framed again does not decode: ${errMsg(e)}`) }
  for (const c of cases) {
    let p = raw.slice()
    p = c.edit?.(p, new DataView(p.buffer)) ?? p
    try {
      decodeTransFrame(reframe(p, c.checksum ?? 1), dof, `FLNC with ${c.rule}`)
      fail(`trans reader rule "${c.rule}": the frame was accepted`)
    } catch (e) {
      if (e instanceof PackError && c.want.test(e.message)) rejected++
      else fail(`trans reader rule "${c.rule}": threw ${errMsg(e)}`)
    }
  }
  pass(rejected === cases.length, `trans: the decoder rejects ${rejected} of ${cases.length} broken copies of FLNC's frame, naming the rule (${cases.map(c => c.rule).join('; ')})`)
}

// ---- trans-only variant pages -----------------------------------------------------------------------

{
  const info = JSON.parse(readFileSync(join(CHECK, 'trans', 'variants.json'), 'utf8'))
  const vr = tableFromIPC(readFileSync(join(CHECK, 'trans', 'variant_rows.arrow')))
  const rid = col(vr, 'range'), vidx = col(vr, 'vidx'), vpos = col(vr, 'position'), vrs = col(vr, 'rs_number'), vaf = col(vr, 'af'), vcode = col(vr, 'af_code')
  const A1 = vr.getChild('A1')!, A2 = vr.getChild('A2')!, noAlleles = vr.getChild('no_alleles')!, match = strCol(vr, 'match')
  const MATCH: Record<string, number> = { none: 0, exact: 1, position: 2 }
  for (const r of info.ranges) {
    const label = `${r.chrom} trans-only page at vidx ${r.first_vidx} (n_cis ${r.n_cis})`
    try {
      const head = new DataView(readRange(r.file, 0, 32).buffer)
      const len = r.byte_end - r.byte_start
      const range = decodeVariantRange(readRange(r.file, r.byte_start, len), { var_len: len }, `${r.file} bytes ${r.byte_start}+${len}`)
      const bad = new Map<string, number>()
      const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
      if (head.getUint32(16, true) !== r.count || head.getUint32(24, true) !== r.n_cis || r.file !== packs.files.variants[r.chrom]) bump('file header count or n_cis, or manifest file')
      if (range.first !== r.first_vidx || range.n !== r.n || range.first < r.n_cis) bump('first_vidx, n, or section')
      const ks = Array.from({ length: vr.numRows }, (_, k) => k).filter(k => rid[k] === r.id)
      if (ks.length !== range.n) bump('reference record count')
      let nulls = 0
      ks.forEach((k, i) => {
        const no = noAlleles.get(k) === true
        if (no) nulls++
        if (vidx[k] !== range.first + i || range.position[i] !== vpos[k]) bump('vidx, position')
        if (no ? range.a1[i] !== null || range.a2[i] !== null : range.a1[i] !== A1.get(k) || range.a2[i] !== A2.get(k)) bump('alleles')
        if (range.rsNumber[i] === 0 ? !Number.isNaN(vrs[k]) : range.rsNumber[i] !== vrs[k]) bump('rs_number')
        if (range.match[i] !== MATCH[match[k]]) bump('match')
        if (range.afCode[i] !== vcode[k] || !(Math.abs(range.af[i] - vaf[k]) <= tol.af_tol)) bump('af')
        if (range.maSamples[i] !== 65535 || range.maCount[i] !== 65535) bump('ma_samples, ma_count null')
      })
      // a gene or intron run never reaches n_cis, so a run holding a record without alleles is rejected
      const j = range.a1.indexOf(null)
      let runRejected = false
      try { sliceRun(range, range.first + j, 1, { posFirst: range.position[j], posLast: range.position[j] }) }
      catch (e) { runRejected = e instanceof PackError && /no alleles/.test(e.message) }
      pass(bad.size === 0 && nulls === r.no_alleles && j >= 0 && runRejected,
        `${label}: ${range.n} records equal the reference (position, alleles, rs_number, match, af code, null counts), ${nulls} with flags bit 2 read as null alleles (reference ${r.no_alleles}); a run over one is rejected ${runRejected}${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
    } catch (e) { fail(`${label}: ${errMsg(e)}`) }
  }
}

// ---- the variant page: hits frames, rsID index, variant index, cis scan (SPEC sections 13 to 15) -----

const VARIANT = join(CHECK, 'variant')

// the whole startup file: every count and offset against the pipeline's own reading of it
let vindex: ReturnType<typeof decodeVariantIndex> | null = null
{
  const ref = JSON.parse(readFileSync(join(VARIANT, 'variant_index.json'), 'utf8'))
  try {
    const bytes = readFileSync(join(DERIVED, ref.file))
    vindex = decodeVariantIndex(new Uint8Array(bytes), ref.file)
    const bad: string[] = []
    if (vindex.pageSize !== ref.page_size || vindex.frameVariants !== ref.frame_variants
      || vindex.rsidBlockRecords !== ref.rsid_block_records || vindex.rsidNRecords !== ref.rsid_n_records
      || vindex.rsidFirst.length !== ref.rsid_n_blocks) bad.push('header constants')
    if (vindex.chroms.size !== Object.keys(ref.chroms).length) bad.push('chromosome count')
    for (const [chr, r] of Object.entries(ref.chroms) as [string, Record<string, number>][]) {
      const c = vindex.chroms.get(chr)
      if (!c) { bad.push(`${chr} missing`); continue }
      if (c.nCis !== r.n_cis || c.nTransOnly !== r.n_trans_only || c.nPagesCis !== r.n_pages_cis
        || c.nPagesTrans !== r.n_pages_trans || c.hitsOff.length - 1 !== r.n_frames) bad.push(`${chr} counts`)
      if (c.pageOff[0] !== r.page_off_first || c.pageOff[c.pageOff.length - 1] !== r.page_off_last
        || c.hitsOff[0] !== r.hits_off_first || c.hitsOff[c.hitsOff.length - 1] !== r.hits_off_last) bad.push(`${chr} first/last offsets`)
    }
    pass(bad.length === 0, `variant_index.qbx decodes to ${vindex.chroms.size} chromosomes, ${vindex.rsidFirst.length.toLocaleString()} rsID blocks `
      + `and ${vindex.rsidNRecords.toLocaleString()} records, matching the reference${bad.length ? ` (${bad.slice(0, 3).join(', ')})` : ''}`)
  } catch (e) { fail(`variant_index.qbx: ${errMsg(e)}`) }
}

// the first and last rsID block: every record, and the binary search that finds it
{
  const ref = JSON.parse(readFileSync(join(VARIANT, 'rsid.json'), 'utf8'))
  const rr = tableFromIPC(readFileSync(join(VARIANT, 'rsid_rows.arrow')))
  const rblock = col(rr, 'block'), rrec = col(rr, 'record'), rnum = col(rr, 'rs_number'), rvidx = col(rr, 'vidx')
  const rchr = strCol(rr, 'chr')
  for (const b of ref.blocks) {
    const label = `rsid_index.qbr block ${b.block} (${b.records.toLocaleString()} records at byte ${b.byte_start})`
    try {
      const block = readRange(ref.file, b.byte_start, b.byte_end - b.byte_start)
      const ks = Array.from({ length: rr.numRows }, (_, k) => k).filter(k => rblock[k] === b.block)
      const bad = new Map<string, number>()
      const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
      if (ks.length !== b.records) bump('reference record count')
      for (const k of ks) {
        const got = rsidInBlock(block, rnum[k])
        if (!got) { bump('not found'); continue }
        if (VIDX_CHROMS[got.chrOrdinal] !== rchr[k] || got.vidx !== rvidx[k]) bump('chromosome or vidx')
        if (rrec[k] === 0 && vindex && vindex.rsidFirst[b.block] !== rnum[k]) bump('rsid_first')
      }
      // a number between two records, and one past the end, are misses
      const gap = ks.find(k => rrec[k] > 0 && rnum[k] > rnum[k - 1] + 1)
      const misses = (gap === undefined || rsidInBlock(block, rnum[gap] - 1) === null) && rsidInBlock(block, 0) === null
      pass(bad.size === 0 && misses, `${label}: every record is found by binary search with its own chromosome and vidx, `
        + `and a number the block does not hold returns null${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
    } catch (e) { fail(`${label}: ${errMsg(e)}`) }
  }
}

// three hits frames, row by row, and the rejection cases
{
  const ref = JSON.parse(readFileSync(join(VARIANT, 'hits.json'), 'utf8'))
  const hr = tableFromIPC(readFileSync(join(VARIANT, 'hits_rows.arrow')))
  const hRef = col(hr, 'ref'), hVidx = col(hr, 'vidx'), hKind = col(hr, 'kind'), hGene = col(hr, 'gene')
  const hPval = col(hr, 'pval'), hBeta = col(hr, 'beta'), hSlope = col(hr, 'slope'), hSe = col(hr, 'slope_se')
  const hPip = col(hr, 'pip'), hCs = col(hr, 'cs_id')
  const hStart = col(hr, 'intron_start'), hEnd = col(hr, 'intron_end'), hClu = col(hr, 'cluster')
  const hStrand = strCol(hr, 'strand'), hSig = hr.getChild('significant')!
  const T = ref.tolerances
  let firstFrame: Uint8Array | null = null
  for (const f of ref.frames) {
    const label = `${f.chrom} hits frame ${f.frame} (${f.why}, ${f.rows.toLocaleString()} rows in ${(f.byte_end - f.byte_start).toLocaleString()} B)`
    try {
      const bytes = readRange(f.file, f.byte_start, f.byte_end - f.byte_start)
      if (firstFrame === null) firstFrame = bytes
      const fr = decodeHitsFrame(bytes, f.first_vidx, label)
      const bad = new Map<string, number>()
      const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
      if (fr.nVariants !== f.n_variants || fr.rowStart[fr.nVariants] !== f.rows) bump('n_variants or row count')
      if (fr.trans.nlpMax !== f.scales.trans_nlp_max || fr.trans.betaMax !== f.scales.trans_beta_max
        || fr.lead.nlpMax !== f.scales.perm_nlp_max || fr.lead.seMax !== f.scales.se_max || fr.lead.slopeMax !== f.scales.slope_max) bump('scales')
      const ks = Array.from({ length: hr.numRows }, (_, k) => k).filter(k => hRef[k] === f.id)
      if (ks.length !== f.rows) bump('reference row count')
      ks.forEach((k, r) => {
        const kind = fr.kind[r]
        if (hVidx[k] !== f.first_vidx + slotOf(fr, r) || kind !== hKind[k] || fr.ord[r] !== hGene[k]) bump('vidx, kind, or gene ord')
        if (kind <= 1) {
          const nlp = (fr.v1[r] * fr.trans.nlpMax) / 65533
          if (!within(Math.pow(10, -nlp), hPval[k], nlp, fr.trans.nlpMax / T.nlp_half_step_divisor + T.nlp_slack)) bump('trans p')
          if (!(Math.abs((fr.v3[r] * fr.trans.betaMax) / 32767 - hBeta[k]) <= fr.trans.betaMax / T.beta_half_step_divisor + 1e-12)) bump('trans beta')
        } else if (kind <= 3) {
          const nlp = (fr.v1[r] * fr.lead.nlpMax) / 65533
          if (!within(Math.pow(10, -nlp), hPval[k], nlp, fr.lead.nlpMax / T.nlp_half_step_divisor + T.nlp_slack)) bump('lead perm p')
          if (!(Math.abs((fr.v2[r] * fr.lead.seMax) / 65535 - hSe[k]) <= fr.lead.seMax / T.se_half_step_divisor + 1e-12)) bump('lead slope_se')
          if (!(Math.abs((fr.v3[r] * fr.lead.slopeMax) / 32767 - hSlope[k]) <= fr.lead.slopeMax / T.slope_half_step_divisor + 1e-12)) bump('lead slope')
          if (((fr.flags[r] & 2) !== 0) !== (hSig.get(k) === true)) bump('significant flag')
        } else {
          if (!(Math.abs(fr.v1[r] / 65535 - hPip[k]) <= T.pip_abs)) bump('pip')
          if (fr.v2[r] !== hCs[k]) bump('cs_id')
        }
        if (kind & 1) {
          if (fr.intronStart[r] !== hStart[k] || fr.intronEnd[r] !== hEnd[k] || fr.cluster[r] !== hClu[k]
            || ((fr.flags[r] & 1) ? '-' : '+') !== hStrand[k]) bump('intron fields')
        } else if (fr.intronStart[r] || fr.intronEnd[r] || fr.cluster[r]) bump('intron fields set on an even kind')
      })
      pass(bad.size === 0, `${label}: every row decodes to the reference (vidx, kind, gene ord, intron fields, and values within `
        + `SPEC section 13)${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
    } catch (e) { fail(`${label}: ${errMsg(e)}`) }
  }
  if (firstFrame) pass(...rejects(firstFrame, ref.frames[0].first_vidx))
}

/** The variant slot a decoded row belongs to. */
function slotOf(fr: ReturnType<typeof decodeHitsFrame>, row: number): number {
  let lo = 0, hi = fr.nVariants - 1
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1
    if (fr.rowStart[mid] <= row) lo = mid; else hi = mid - 1
  }
  return lo
}

/** A decoded value is within its half step, compared on the -log10 p scale where p is tiny. */
function within(got: number, want: number, nlp: number, limit: number): boolean {
  if (Number.isNaN(want)) return Number.isNaN(got)
  if (want === 0 || nlp > 300) return Math.abs(Math.log10(got || Number.MIN_VALUE) + Math.log10(want || Number.MIN_VALUE)) <= limit
  return Math.abs(-Math.log10(got) - -Math.log10(want)) <= limit
}

/** Every rule of SPEC section 13 that a broken frame must trip. The frame is re-compressed after
 *  each edit, because the decoder decompresses before it checks anything. */
function rejects(frame: Uint8Array, firstVidx: number): [boolean, string] {
  const good = zstdDecompressSync(frame)
  const reframe = (p: Uint8Array) => new Uint8Array(zstdCompressSync(p, { params: { [zlibConstants.ZSTD_c_checksumFlag]: 1 } }))
  const edits: [string, (p: Uint8Array) => Uint8Array | null][] = [
    ['magic', p => { p[0] ^= 0xff; return p }],
    ['reserved field', p => { p[6] = 1; return p }],
    ['n_variants 0', p => { p[4] = 0; p[5] = 0; return p }],
    ['scale not finite', p => { new DataView(p.buffer, p.byteOffset).setFloat64(8, NaN, true); return p }],
    ['negative scale', p => { new DataView(p.buffer, p.byteOffset).setFloat64(16, -1, true); return p }],
    ['payload length', p => p.subarray(0, p.length - 4)],
    ['too short', p => p.subarray(0, 40)],
  ]
  const caught: string[] = []
  for (const [name, edit] of edits) {
    const p = edit(new Uint8Array(good))
    if (!p) continue
    let threw = false
    try { decodeHitsFrame(reframe(p), firstVidx, `broken (${name})`) } catch (e) { threw = e instanceof PackError }
    if (threw) caught.push(name)
  }
  // a frame whose bytes are not a zstd frame at all
  let rawThrew = false
  try { decodeHitsFrame(new Uint8Array(good), firstVidx, 'not a zstd frame') } catch (e) { rawThrew = e instanceof PackError }
  if (rawThrew) caught.push('not a zstd frame')
  const want = edits.length + 1
  return [caught.length === want, `hits: the decoder rejects ${caught.length} of ${want} broken copies of the reference frame, naming the rule (${caught.join('; ')})`]
}

// the cis scan: the same two spans the pipeline walked, block by block
{
  const ref = JSON.parse(readFileSync(join(VARIANT, 'scan.json'), 'utf8'))
  const sr = tableFromIPC(readFileSync(join(VARIANT, 'scan_rows.arrow')))
  const sScan = strCol(sr, 'scan'), sType = strCol(sr, 'qtl_type'), sPid = strCol(sr, 'phenotype_id')
  const sP = col(sr, 'pval_nominal'), sSlope = col(sr, 'slope'), sSe = col(sr, 'slope_se')
  const sNlpMax = col(sr, 'nlp_max'), sLseMin = col(sr, 'lse_min'), sLseMax = col(sr, 'lse_max')
  for (const scan of ref.scans) {
    for (const [type, span, dof] of [['eqtl', scan.eqtl, packs.dof.eqtl], ['sqtl', scan.sqtl, packs.dof.sqtl]] as [string, Record<string, number>, number][]) {
      const label = `cis scan ${scan.id} ${type} span (${span.blocks} blocks in ${span.bytes.toLocaleString()} B)`
      try {
        const bytes = readRange(span.file as unknown as string, span.byte_start, span.bytes)
        // walk the span by its block headers, exactly as cis-scan.ts does
        const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength)
        const got = new Map<number, ReturnType<typeof scanBlockRow>>()
        let at = 0, blocks = 0, covering = 0
        while (at < bytes.length) {
          const blkLen = dv.getUint32(at + 4, true), nRows = dv.getUint32(at + 8, true), varStart = dv.getUint32(at + 12, true)
          blocks++
          if (nRows > 0 && scan.vidx >= varStart && scan.vidx < varStart + nRows) {
            got.set(span.byte_start + at, scanBlockRow(bytes.subarray(at, at + blkLen), scan.vidx - varStart, dof, label))
            covering++
          }
          at += blkLen
        }
        const ks = Array.from({ length: sr.numRows }, (_, k) => k).filter(k => sScan[k] === scan.id && sType[k] === type)
        const bad = new Map<string, number>()
        const bump = (k: string) => bad.set(k, (bad.get(k) ?? 0) + 1)
        if (blocks !== span.blocks) bump('block count')
        if (covering !== ks.length) bump('covering block count')
        const rowsByP = [...got.values()].sort((a, b) => a.pval - b.pval)
        const wantByP = ks.slice().sort((a, b) => sP[a] - sP[b])
        wantByP.forEach((k, i) => {
          const r = rowsByP[i]
          if (!r) { bump('missing row'); return }
          const nlp = -Math.log10(sP[k])
          if (!within(r.pval, sP[k], nlp, sNlpMax[k] / ref.nlp_half_step_divisor + 1e-9)) bump('p')
          const hLse = (sLseMax[k] - sLseMin[k]) / 65532
          const seBound = sSe[k] * (Math.expm1(hLse) + 1e-12)
          if (!(Math.abs(r.se - sSe[k]) <= ref.bound_factor * seBound + 1e-12)) bump('slope_se')
          // slope follows from p and SE, so its bound grows with the p step; compare on the SE scale
          if (!(Math.abs(r.slope - sSlope[k]) <= ref.bound_factor * (seBound * (Math.abs(sSlope[k] / sSe[k]) + 1) + sSe[k]) + 1e-9)) bump('slope')
          if (sPid[k] === '') bump('reference phenotype id')
        })
        pass(bad.size === 0, `${label}: ${covering} covering blocks decode to the reference rows, values within SPEC section 9`
          + `${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
      } catch (e) { fail(`${label}: ${errMsg(e)}`) }
    }
  }
}

// ---- inverse t against scipy ------------------------------------------------------------------------

{
  const grid = tableFromIPC(readFileSync(join(CHECK, 't_grid.arrow')))
  const nlp = col(grid, 'nlp'), dof = col(grid, 'dof'), tRef = col(grid, 't')
  let worstRel = 0, worstAbs = 0, bad = 0, worstAt = ''
  const t0 = performance.now()
  for (let i = 0; i < nlp.length; i++) {
    const t = tFromNlp(nlp[i], dof[i])
    const err = Math.abs(t - tRef[i])
    if (tRef[i] >= tol.t_grid.small_t) {
      const rel = err / tRef[i]
      if (rel > worstRel) { worstRel = rel; worstAt = `nlp ${nlp[i]} dof ${dof[i]}` }
      if (!(rel <= tol.t_grid.rel)) bad++
    } else {
      worstAbs = Math.max(worstAbs, err)
      if (!(err <= tol.t_grid.abs)) bad++
    }
  }
  const perCall = (performance.now() - t0) / nlp.length
  pass(bad === 0, `tFromNlp matches scipy stdtrit on ${nlp.length.toLocaleString()} points: worst relative error ${g(worstRel)} (${worstAt}) where t >= ${tol.t_grid.small_t}, worst absolute ${g(worstAbs)} below; ${bad} over the limit; ${(perCall * 1000).toFixed(1)} us per call`)
  // below the smallest normal double p scipy is no reference: t must stay finite and increase with nlp
  let prev = 0, monotone = true
  for (let x = 300; x <= 340; x += 0.25) {
    const t = tFromNlp(x, 435)
    if (!(Number.isFinite(t) && t > prev)) monotone = false
    prev = t
  }
  pass(monotone, `tFromNlp is finite and increasing for nlp 300 to 340 (t at 340: ${g(tFromNlp(340, 435))})`)
}

for (const fd of fds.values()) closeSync(fd)
if (failures.length) { console.log(`pack-check: ${failures.length} check(s) failed`); process.exit(1) }
console.log('pack-check: all checks passed')
