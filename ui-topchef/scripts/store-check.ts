/**
 * Round trip of the browser decoder (src/lib/store-decode.ts) against the Python decoders
 * (pipeline/catalog.py, pipeline/results.py, pipeline/gwas.py, pipeline/packfmt_v2.py) on one
 * collection of a local version 2 store.
 *
 *   uv run python ui-topchef/scripts/store_reference.py <store> <ref.json> [--anchor topchef_ge_cis] [--gwas dcm_jurgens2024_biobanks]
 *   npm run store-check -- <store> <ref.json>
 *
 * The reference holds a seeded sample of every object family in full and exact totals over every
 * row (its docstring says which); this side decodes everything and recomputes both:
 * - variants files: every site of every chromosome (both sections, walked with the TS variant
 *   index) summed against `catalog.decode_file`, and the sampled pages column by column, one page
 *   at a time through decodeVariantPage too;
 * - the variant index against `catalog.decode_vidx`; every rsID record summed against `decode_rsid`,
 *   the sampled blocks exactly, and rsidFind for a sample of rs numbers (repeats and block edges
 *   included) against `catalog.rsid_lookup`;
 * - hits files of every experiment, paged: the frame table against `packfmt_v1.hits_frame_table`,
 *   every frame's records summed against `results.decode_hits`, the sampled frames exactly;
 * - the sampled trans frames against `results.read_trans` (codes, positions, alleles, credible sets
 *   exactly; p, beta and SE to rounding), and the layout from the search index: the frames tile each
 *   trans object and each gene's frames are one contiguous range (what the gene page reads in one request);
 * - GWAS windows (the sampled blocks' runs and the smallest chromosome whole) against `gwas.read_window`,
 *   every chromosome whole summed against it; every row of the bin summary against `gwas.read_bins`;
 * - the sampled result blocks against `results.read_block`: codes, credible sets, details and
 *   positions exactly; -log10 p, p, beta and SE as the largest relative difference; scanBlockRow
 *   row by row against decodeResultBlock;
 * - what a first gene page reads: the gene lookup's key normalization and bucket hash on pinned
 *   vectors (the same ones as pipeline/test_annotation.py), its directory and the rows for a sample
 *   of keys, misses included; every chromosome's genes and exon models; every experiment's search
 *   index parts (the trans experiments' unplaced part included) with their rows and ord runs; the
 *   collection document; and the coloc lists for the collection and its GWAS, with every listed gene
 *   placed by the gene lookup as the Home track places it.
 * Integers, strings and codes must match exactly. Exits 1 on any failure.
 */
import { closeSync, openSync, readFileSync, readSync, statSync } from 'node:fs'
import { join } from 'node:path'
import { tableFromIPC } from 'apache-arrow'
import { seqcolDigest } from '../src/lib/chrom-sizes.ts'
import { checkFileHeader, checkNames, decodeArrowObject, decodeLookupDir, lookupBucket, lookupDirLen, lookupKey, decodeGwasBins, decodeGwasIndex, decodeGwasRange, decodeHitsFrame, decodeHitsTable, decodeResultBlock, decodeRsidRecords,
  decodeTransFrame, decodeVariantIndex, decodeVariantPage, decodeVariantRange, FORMAT_VERSION, gwasRange, HEADER_LEN, hitsOf, hitsTableLen, KIND,
  parseFileHeader, readerColumns, RSID_RECORD_LEN, rsidFind, scanBlockRow, sliceRun,
  type RsidRecord } from '../src/lib/store-decode.ts'

const [STORE, REF] = process.argv.slice(2)
if (!STORE || !REF) { console.error('usage: npm run store-check -- <store dir> <reference.json>'); process.exit(2) }
const ref = JSON.parse(readFileSync(REF, 'utf8'))
const imm = (name: string) => join(STORE, 'immutable', name)
const failures: string[] = []
const check = (ok: boolean, msg: string) => { if (!ok) failures.push(msg); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }
const g = (x: number) => x.toPrecision(3)

function readRange(name: string, off: number, len: number): Uint8Array {
  const fd = openSync(imm(name), 'r')
  try {
    const out = new Uint8Array(len)
    for (let got = 0; got < len;) {
      const r = readSync(fd, out, got, len - got, off + got)
      if (r === 0) throw new Error(`${name}: ends before byte ${off + len}`)
      got += r
    }
    return out
  } finally { closeSync(fd) }
}
const whole = (name: string) => new Uint8Array(readFileSync(imm(name)))
const pointer = (level: string, id: string) => JSON.parse(readFileSync(join(STORE, level, `${id}.json`), 'utf8'))

/** First index where two sequences differ, or -1. */
function firstDiff(a: ArrayLike<unknown>, b: ArrayLike<unknown>): number {
  if (a.length !== b.length) return Math.min(a.length, b.length)
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return i
  return -1
}
/** Python floats come as numbers, null (NaN) or "inf"/"-inf". */
const py = (x: number | string | null) => (x === null ? NaN : x === 'inf' ? Infinity : x === '-inf' ? -Infinity : x)
/** The largest relative difference, and whether NaN/infinite entries sit in the same places. */
function relDiff(ts: ArrayLike<number>, pyv: (number | string | null)[]): { max: number; exact: number; placed: boolean } {
  let max = 0, exact = 0, placed = ts.length === pyv.length
  for (let i = 0; i < Math.min(ts.length, pyv.length); i++) {
    const a = ts[i], b = py(pyv[i])
    if (!Number.isFinite(a) || !Number.isFinite(b)) { if (!Object.is(a, b) && !(Number.isNaN(a) && Number.isNaN(b))) placed = false; else exact++; continue }
    if (a === b) { exact++; continue }
    max = Math.max(max, Math.abs(a - b) / Math.max(Math.abs(b), 1e-300))
  }
  return { max, exact, placed }
}
/** Exact integer sum (every total here sits well below 2^53). */
const isum = (a: ArrayLike<number>) => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s }
/** Sum of the finite entries and how many were not; compared to Python's fsum within `TOL`. */
function fsum(a: ArrayLike<number>, map: (x: number) => number = x => x): { sum: number; bad: number } {
  let sum = 0, bad = 0
  for (let i = 0; i < a.length; i++) { const x = map(a[i]); if (Number.isFinite(x)) sum += x; else bad++ }
  return { sum, bad }
}
const TOL = 1e-8
const near = (a: number, b: number) => Math.abs(a - b) <= TOL * Math.max(Math.abs(b), 1e-300)

// ---- store.json, the collection, the experiments ----------------------------------------------------
const storeDoc = JSON.parse(readFileSync(join(STORE, 'store.json'), 'utf8'))
const expIds = Object.keys(ref.experiments) as string[]
const docs = Object.fromEntries(expIds.map(e => [e, pointer('experiments', e)]))
const cisIds = expIds.filter(e => ref.experiments[e].kind === 'qtl_cis'), transIds = expIds.filter(e => ref.experiments[e].kind === 'qtl_trans')
{
  const kinds = expIds.filter(e => storeDoc.experiments?.[e]?.kind !== ref.experiments[e].kind || storeDoc.experiments[e].collection !== ref.collection.name)
  check(storeDoc.format_version === FORMAT_VERSION && kinds.length === 0,
    `store.json: format version ${storeDoc.format_version}; ${expIds.length} experiments of ${ref.collection.name} with the reference's kinds (${cisIds.length} cis, ${transIds.length} trans)${kinds.length ? `; off: ${kinds.join(', ')}` : ''}`)
  const cdoc = JSON.parse(readFileSync(imm(ref.collection.name), 'utf8'))
  check(JSON.stringify(cdoc) === JSON.stringify(ref.collection.doc) && typeof cdoc.title === 'string',
    `collection ${ref.collection.name}: "${cdoc.title}" (${cdoc.citation?.first_author ?? 'no citation'}, ${cdoc.citation?.year ?? '?'})`)
}

const catalog = pointer(ref.catalog_dir, ref.catalog)
const names: string[] = catalog.chromosomes.map((c: { name: string }) => c.name)
{
  // the app's own document check (openStore runs it on every pointer it reads), so a provenance
  // string that only looks like an object name fails here and not in the browser
  const annotation = pointer('annotations', docs[expIds[0]].annotation)
  const all: [string, unknown][] = [...expIds.map(e => [`experiments/${e}.json`, docs[e]] as [string, unknown]),
    ...(ref.gwas_experiment ? [[`experiments/${ref.gwas_experiment}.json`, pointer('experiments', ref.gwas_experiment)] as [string, unknown]] : []),
    [`${ref.catalog_dir}/${ref.catalog}.json`, catalog], [`annotations/${annotation.id}.json`, annotation]]
  const bad = all.flatMap(([where, d]) => { try { checkNames(d, where); return [] } catch (e) { return [(e as Error).message] } })
  check(bad.length === 0, `openStore's name check accepts all ${all.length} pointer documents${bad.length ? `: ${bad.join('; ')}` : ''}`)
}

// ---- variant index ----------------------------------------------------------------------------
// the genome track fetches its chromosome lengths from seqcol by a digest of its own; it has to be
// the collection the catalog's positions are anchored to, or the track is drawn against a reference
// the data does not use
check(seqcolDigest() === catalog.collection_digest,
  `genome track collection (lib/chrom-sizes.ts) is the catalog's ${catalog.collection_digest}${seqcolDigest() === catalog.collection_digest ? '' : `, but it cites ${seqcolDigest()}`}`)

const vx = decodeVariantIndex(whole(catalog.vidx), names, catalog.collection_digest, catalog.vidx)
{
  const p = ref.vidx
  let bad = firstDiff(vx.rsidFirst, p.rsid_first) >= 0 || vx.pageSize !== p.page_size || vx.rsidBlockRecords !== p.rsid_block_records || vx.rsidN !== p.rsid_n
  for (const n of names) {
    const a = vx.chroms.get(n)!, b = p.chroms[n]
    bad ||= a.nCis !== b.n_cis || a.nTrans !== b.n_trans || firstDiff(a.pageOff, b.page_off) >= 0 || firstDiff(a.pageFirstPos, b.page_first_position) >= 0
  }
  check(!bad, `variant index ${catalog.vidx}: ${names.length} chromosomes, page size ${vx.pageSize}, ${vx.rsidFirst.length} rsID blocks equal decode_vidx`)
}

// ---- variants files ---------------------------------------------------------------------------
let sites = 0
for (const c of catalog.chromosomes as { name: string; file: string; seq_digest: string; count: number; n_cis: number }[]) {
  const want = ref.chroms[c.name]
  const h = parseFileHeader(readRange(c.file, 0, HEADER_LEN), c.file)
  checkFileHeader(h, { kind: KIND.variants, chrom: c.name, seqDigest: c.seq_digest }, c.file)
  check(h.count === want.header.count && h.nCis === want.header.n_cis && h.pageSize === want.header.page_size && h.seqDigest === want.header.seq_digest,
    `${c.name} ${c.file}: header kind 1, count ${h.count}, n_cis ${h.nCis}, page size ${h.pageSize}, seq_digest ${h.seqDigest} equal parse_file_header`)
  const ci = vx.chroms.get(c.name)!
  const cols = { pos: [] as number[], ref: [] as string[], alt: [] as string[], af_code: [] as number[], af: [] as number[], rs_number: [] as number[],
    ma_samples: [] as number[], ma_count: [] as number[], match: [] as number[], alt_is_minor: [] as number[] }
  for (const [p0, p1] of [[0, ci.nPagesCis], [ci.nPagesCis, ci.nPagesCis + ci.nPagesTrans]]) {
    if (p0 === p1) continue
    const off = ci.pageOff[p0], len = ci.pageOff[p1] - off
    const r = decodeVariantRange(readRange(c.file, off, len), len, `${c.file} pages ${p0}..${p1 - 1}`)
    const add = <T,>(to: T[], from: ArrayLike<T>) => { for (let i = 0; i < from.length; i++) to.push(from[i]) }
    add(cols.pos, r.position); add(cols.ref, r.ref); add(cols.alt, r.alt); add(cols.af_code, r.afCode); add(cols.af, r.af)
    add(cols.rs_number, r.rsNumber); add(cols.ma_samples, r.maSamples); add(cols.ma_count, r.maCount)
    add(cols.match, r.match); add(cols.alt_is_minor, r.altIsMinor)
  }
  const t = want.totals
  const got = { pos: isum(cols.pos), af_code: isum(cols.af_code), rs_number: isum(cols.rs_number), ma_samples: isum(cols.ma_samples),
    ma_count: isum(cols.ma_count), match: isum(cols.match), alt_is_minor: isum(cols.alt_is_minor),
    allele_len: isum(cols.ref.map(x => x.length)) + isum(cols.alt.map(x => x.length)),
    allele_first: isum(cols.ref.map(x => x.charCodeAt(0))) + isum(cols.alt.map(x => x.charCodeAt(0))) }
  const offTotals = Object.keys(got).filter(k => (got as Record<string, number>)[k] !== t[k])
  check(offTotals.length === 0 && cols.pos.length === c.count && c.count === want.count,
    `${c.name}: all ${cols.pos.length} sites sum to decode_file's totals (pos, af code, rs_number, counts, match, alt_is_minor, alleles)${offTotals.length ? `: ${offTotals.join(', ')}` : ''}`)
  // the sampled pages: the same columns, and one page at a time as the variant page reads them
  const badPages: string[] = []
  const nPages = ci.nPagesCis + ci.nPagesTrans
  for (const [k, p] of Object.entries(want.pages as Record<string, { first_vidx: number; n: number; pos: number[]; ref: string[]; alt: string[]; af_code: number[]
    af: (number | null)[]; rs_number: number[]; ma_samples: number[]; ma_count: number[]; match: number[]; alt_is_minor: number[] }>)) {
    const lo = p.first_vidx, hi = lo + p.n
    const diffs = (Object.keys(cols) as (keyof typeof cols)[]).filter(key => {
      const v = cols[key].slice(lo, hi)
      if (key === 'af') return (v as number[]).some((x, j) => !(Object.is(x, py(p.af[j])) || (Number.isNaN(x) && p.af[j] === null)))
      return firstDiff(v, p[key] as unknown[]) >= 0
    })
    if (diffs.length) badPages.push(`page ${k}: ${diffs.join(', ')}`)
    const kk = Number(k)
    const recs = decodeVariantPage(readRange(c.file, ci.pageOff[kk], ci.pageOff[kk + 1] - ci.pageOff[kk]), c.name, ci.nCis)
    if (recs.length !== p.n || recs.some((r, j) => r.vidx !== lo + j || r.position !== p.pos[j] || r.A1 !== p.alt[j] || r.A2 !== p.ref[j] || r.inCis !== r.vidx < ci.nCis))
      badPages.push(`page ${k} single`)
  }
  check(badPages.length === 0 && Number(want.n_pages) === nPages, `${c.name}: ${Object.keys(want.pages).length} sampled pages of ${nPages} equal decode_file column by column, and decode alone (A1 = alt, A2 = ref, inCis by n_cis)${badPages.length ? `: ${badPages.slice(0, 3).join('; ')}` : ''}`)
  sites += cols.pos.length
}

// ---- hits: every experiment, every chromosome ---------------------------------------------------
{
  let files = 0, records = 0
  const bad: string[] = []
  const kinds = [0, 0, 0]
  for (const e of expIds) for (const c of catalog.chromosomes as { name: string; file: string; seq_digest: string; count: number }[]) {
    const hw = ref.hits[e][c.name]
    const what = `${e} ${c.name}`
    const hb = whole(hw.file)
    const tab = decodeHitsTable(hb.subarray(0, hitsTableLen(c.count, 1024)), { chrom: c.name, seqDigest: c.seq_digest, nVariants: c.count }, hw.file)
    if (firstDiff(tab.frameOff, hw.frame_off) >= 0 || tab.count !== hw.header.count || tab.frameVariants !== hw.header.page_size) bad.push(`${what}: frame table`)
    const sum = { vidx: 0, ord: 0, cs_id: 0, flags: 0, value: 0, value_bad: 0, abs_beta: 0, beta_bad: 0 }
    const k = [0, 0, 0]
    let n = 0, probed = false
    for (let f = 0; f + 1 < tab.frameOff.length; f++) {
      const fr = decodeHitsFrame(hb.subarray(tab.frameOff[f], tab.frameOff[f + 1]), f * tab.frameVariants, tab.frameVariants, `${hw.file} frame ${f}`)
      n += fr.count
      for (let i = 0; i < fr.count; i++) {
        sum.vidx += fr.vidx[i]; sum.ord += fr.ord[i]; sum.cs_id += fr.csId[i]; sum.flags += fr.flags[i]; k[fr.kind[i]]++
        if (Number.isFinite(fr.value[i])) sum.value += fr.value[i]; else sum.value_bad++
        if (Number.isFinite(fr.beta[i])) sum.abs_beta += Math.abs(fr.beta[i]); else sum.beta_bad++
      }
      const w = hw.frames[String(f)]
      if (w) {
        const d = (['vidx', 'ord', 'kind', 'cs_id', 'flags'] as const).filter(key => firstDiff(key === 'kind' ? fr.kind : key === 'cs_id' ? fr.csId : fr[key], w[key]) >= 0)
        const v = relDiff(fr.value, w.value), b = relDiff(fr.beta, w.beta)
        if (d.length || v.max !== 0 || !v.placed || b.max !== 0 || !b.placed) bad.push(`${what} frame ${f}: ${d.join(', ') || 'value or beta'}`)
        if (!probed && fr.count) {
          probed = true
          const v0 = fr.vidx[Math.floor(fr.count / 2)], [lo, hi] = hitsOf(fr, v0)
          if (!(hi > lo && Array.from(fr.vidx.subarray(lo, hi)).every(x => x === v0) && (lo === 0 || fr.vidx[lo - 1] < v0) && (hi === fr.count || fr.vidx[hi] > v0)))
            bad.push(`${what} frame ${f}: hitsOf(${v0})`)
        }
      }
    }
    const t = hw.totals
    if (n !== hw.n || n !== tab.count || firstDiff(k, hw.kinds) >= 0 || sum.vidx !== t.vidx || sum.ord !== t.ord || sum.cs_id !== t.cs_id || sum.flags !== t.flags
      || sum.value_bad !== t.value_not_finite || sum.beta_bad !== t.beta_not_finite || !near(sum.value, t.value) || !near(sum.abs_beta, t.abs_beta))
      bad.push(`${what}: totals`)
    files++; records += n; for (let i = 0; i < 3; i++) kinds[i] += k[i]
  }
  check(bad.length === 0, `hits: ${files} files, ${records} records (leads ${kinds[0]}, credible sets ${kinds[1]}, trans ${kinds[2]}): frame tables, sums over every record and the sampled frames equal decode_hits${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : ''}`)
}

// ---- rsID index -------------------------------------------------------------------------------
{
  const R = ref.rsid
  const all = decodeRsidRecords(readRange(R.file, HEADER_LEN, R.n * RSID_RECORD_LEN), R.file)
  const h = parseFileHeader(readRange(R.file, 0, HEADER_LEN), R.file)
  checkFileHeader(h, { kind: KIND.rsid, chrom: 'all', seqDigest: catalog.collection_digest }, R.file)
  const B = vx.rsidBlockRecords
  const block = (b: number): RsidRecord[] => all.slice(b * B, (b + 1) * B)
  const t = R.totals
  const sums = all.reduce((s, r) => { s.rs += r.rsNumber; s.v += r.vidx; s.o += r.ordinal; return s }, { rs: 0, v: 0, o: 0 })
  const badBlocks = Object.entries(R.blocks as Record<string, { rs_number: number[]; vidx: number[]; ordinal: number[] }>)
    .filter(([b, w]) => { const got = block(Number(b)); return firstDiff(got.map(r => r.rsNumber), w.rs_number) >= 0 || firstDiff(got.map(r => r.vidx), w.vidx) >= 0 || firstDiff(got.map(r => r.ordinal), w.ordinal) >= 0 })
  check(h.count === all.length && all.length === R.n && statSync(imm(R.file)).size === HEADER_LEN + R.n * RSID_RECORD_LEN
    && sums.rs === t.rs_number && sums.v === t.vidx && sums.o === t.ordinal && badBlocks.length === 0,
    `rsID index ${R.file}: all ${all.length} records sum to decode_rsid's totals, ${Object.keys(R.blocks).length} sampled blocks equal; header kind 8 on the collection digest${badBlocks.length ? `; blocks off: ${badBlocks.map(x => x[0]).join(', ')}` : ''}`)
  // rs numbers to try: a spread, repeated ones, both sides of each block edge, and misses. The
  // records are sorted by rs number, so the first of a run is a lower-bound search (a Map over
  // nine million records exceeds V8's limit).
  const first = (rs: number): RsidRecord | null => {
    let lo = 0, hi = all.length
    while (lo < hi) { const mid = (lo + hi) >>> 1; if (all[mid].rsNumber < rs) lo = mid + 1; else hi = mid }
    return lo < all.length && all[lo].rsNumber === rs ? all[lo] : null
  }
  const repeats = all.filter((r, i) => i > 0 && all[i - 1].rsNumber === r.rsNumber).map(r => r.rsNumber)
  const edges = Array.from(vx.rsidFirst).flatMap((x, b) => [x, all[Math.max(0, b * B - 1)].rsNumber])
  const tries = [...new Set([...all.filter((_, i) => i % 9973 === 0).map(r => r.rsNumber), ...repeats.filter((_, i) => i % 97 === 0), ...edges, 1, 2, all[all.length - 1].rsNumber + 1,
    all[0].rsNumber - 1, all[5].rsNumber + 1])]
  let bad = 0
  for (const rs of tries) {
    const got = await rsidFind(vx.rsidFirst, rs, block)
    const want = first(rs)
    if (!(got === want || (got && want && got.vidx === want.vidx && got.ordinal === want.ordinal))) bad++
  }
  check(bad === 0, `rsidFind: ${tries.length} rs numbers (${new Set(repeats).size} repeated, block edges, misses) find the first record of their run (${bad} wrong)`)

  // against catalog.rsid_lookup on the real index
  let badLookup = 0, n = 0
  for (const [rs, want] of Object.entries(ref.rsid_lookup) as [string, [number, number][]][]) {
    n++
    const got = await rsidFind(vx.rsidFirst, Number(rs), block)
    if (!(want.length ? got && got.ordinal === want[0][0] && got.vidx === want[0][1] : got === null)) badLookup++
  }
  check(badLookup === 0, `rsidFind equals catalog.rsid_lookup's first record on ${n} rs numbers (${badLookup} wrong)`)
}

// a repeated rs number whose run straddles block edges (no repeats occur in the smoke stores):
// blocks of 3 records [5 7 7 | 7 7 9 | 9 11 12]; the run of 7 starts in block 0, which SPEC's
// "last block whose first number is at or below" rule (block 1) would miss
{
  const recs = [5, 7, 7, 7, 7, 9, 9, 11, 12].map((rs, i) => ({ rsNumber: rs, vidx: i, ordinal: 1 }))
  const first = Uint32Array.from([5, 7, 9])
  const read = (b: number) => recs.slice(3 * b, 3 * b + 3)
  const got = await Promise.all([7, 9, 5, 12, 8, 4, 13].map(rs => rsidFind(first, rs, read)))
  check(got.map(r => r?.vidx ?? null).join(',') === '1,5,0,8,,,',
    `rsidFind on a synthetic straddling run: rs7 -> vidx 1, rs9 -> 5, rs5 -> 0, rs12 -> 8, misses null (got ${got.map(r => r?.vidx ?? 'null').join(',')})`)
}

// ---- trans frames -----------------------------------------------------------------------------
{
  const bad: string[] = []
  let rows = 0, cs = 0
  const worst = { p: 0, beta: 0, se: 0 }
  const headers = new Set<string>()
  for (const t of ref.trans) {
    if (!headers.has(t.file)) {
      headers.add(t.file)
      checkFileHeader(parseFileHeader(readRange(t.file, 0, HEADER_LEN), t.file), { kind: KIND.trans, chrom: 'all', seqDigest: catalog.collection_digest }, t.file)
    }
    const fr = decodeTransFrame(readRange(t.file, t.trans_off, t.trans_len), `${t.file} ${t.phenotype_id}`)
    const fail = (k: string) => bad.push(`${t.experiment} ${t.phenotype_id}: ${k}`)
    if (fr.n !== t.n_trans || fr.n !== t.pos.length) fail('row count')
    if (fr.nlpMax !== py(t.nlp_max) || fr.betaMax !== py(t.beta_max) || fr.lseMin !== py(t.lse_min) || fr.lseMax !== py(t.lse_max)) fail('scales')
    for (const [k, a, b] of [['pos', fr.position, t.pos], ['rs_number', fr.rsNumber, t.rs_number], ['af_code', fr.afCode, t.af_code],
      ['nlp_code', fr.nlpCode, t.nlp_code], ['beta_code', fr.betaCode, t.beta_code], ['se_code', fr.seCode, t.se_code],
      ['ref', fr.ref, t.ref], ['alt', fr.alt, t.alt], ['cs_row', fr.csRow, t.cs_row], ['cs_id', fr.csId, t.cs_id]] as const)
      if (firstDiff(a as ArrayLike<unknown>, b as unknown[]) >= 0) fail(k)
    if (fr.csRow.length !== t.n_cs || relDiff(fr.csPip, t.cs_pip).max !== 0) fail('credible sets')
    if (firstDiff(Array.from(fr.ordinal, o => names[o - 1]), t.chr) >= 0) fail('chromosome')
    if (relDiff(fr.af, t.af).max !== 0 || relDiff(fr.nlp, t.nlp).max !== 0) fail('af or -log10 p')
    const p = relDiff(fr.pval, t.p), be = relDiff(fr.beta, t.beta), se = relDiff(fr.se, t.se)
    if (!p.placed || !be.placed || !se.placed) fail('null or infinite placement')
    worst.p = Math.max(worst.p, p.max); worst.beta = Math.max(worst.beta, be.max); worst.se = Math.max(worst.se, se.max)
    rows += fr.n; cs += fr.csRow.length
  }
  check(bad.length === 0, `${ref.trans.length} sampled trans frames, ${rows} rows, ${cs} credible-set records: positions, chromosomes, rsID, af, codes, alleles and credible sets equal read_trans${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : ''}`)
  check(worst.p <= 1e-13 && worst.beta <= 1e-13 && worst.se <= 1e-13,
    `trans values: largest relative difference p ${g(worst.p)}, beta ${g(worst.beta)}, SE ${g(worst.se)} (all three stored; nothing derived)`)
  // layout, from the search index: per object the frames tile it from byte 64; per gene the frames are back to back
  const layoutBad: string[] = []
  let genes = 0, maxFrames = 0, frames = 0
  for (const e of transIds) {
    const file = docs[e].file as string
    const fs: { off: number; len: number; gene: string | null }[] = []
    for (const part of Object.values(ref.index_parts[e] as Record<string, { trans_off: (number | null)[]; trans_len: (number | null)[]; gene_id: (string | null)[] }>))
      part.trans_off.forEach((off, i) => { if (off !== null) fs.push({ off, len: part.trans_len[i]!, gene: part.gene_id[i] }) })
    fs.sort((a, b) => a.off - b.off)
    let end = HEADER_LEN
    for (const f of fs) { if (f.off !== end) { layoutBad.push(`${e}: gap or overlap at ${f.off}`); break } end = f.off + f.len }
    if (end !== statSync(imm(file)).size) layoutBad.push(`${e}: frames end at ${end}, not at the end of ${file}`)
    const byGene = new Map<string, { off: number; len: number }[]>()
    for (const f of fs) if (f.gene !== null) { if (!byGene.has(f.gene)) byGene.set(f.gene, []); byGene.get(f.gene)!.push(f) }
    for (const [gene, xs] of byGene) {
      genes++; maxFrames = Math.max(maxFrames, xs.length)
      if (xs.some((x, i) => i > 0 && x.off !== xs[i - 1].off + xs[i - 1].len)) layoutBad.push(`${e}: ${gene} frames not contiguous`)
    }
    frames += fs.length
  }
  check(layoutBad.length === 0, `trans layout: ${transIds.length} objects tiled by their ${frames} frames; ${genes} (object, gene) pairs each one contiguous range (up to ${maxFrames} frames)${layoutBad.length ? `: ${layoutBad.slice(0, 4).join('; ')}` : ''}`)
}

// ---- GWAS -------------------------------------------------------------------------------------
{
  const G = ref.gwas
  if (!G) check(!ref.gwas_experiment, 'no GWAS experiment named')
  else {
    const gw = pointer('experiments', G.experiment)
    const index = decodeGwasIndex(whole(gw.index), catalog.collection_digest, gw.index)
    const seqOf = (chr: string) => catalog.chromosomes.find((c: { name: string }) => c.name === chr).seq_digest
    const headers = new Set<string>()
    const read = (chr: string, lo: number, hi: number) => {
      const file = gw.files[chr]
      if (!headers.has(file)) { headers.add(file); checkFileHeader(parseFileHeader(readRange(file, 0, HEADER_LEN), file), { kind: KIND.gwas, chrom: chr, seqDigest: seqOf(chr) }, file) }
      const range = gwasRange(index, chr, lo, hi)
      return range ? decodeGwasRange(readRange(file, range.off, range.len), index, chr, range, lo, hi) : null
    }
    const bad: string[] = []
    let rows = 0, worst = 0, worstSe = 0
    for (const w of G.windows) {
      const got = read(w.chr, w.lo, w.hi)
      const fail = (k: string) => bad.push(`${w.chr}:${w.lo}-${w.hi}: ${k}`)
      if ((got?.rows ?? 0) !== w.pos.length) { fail(`rows ${got?.rows ?? 0} != ${w.pos.length}`); continue }
      if (!got) continue
      for (const [k, a, b] of [['pos', got.position, w.pos], ['ref', got.ref, w.ref], ['alt', got.alt, w.alt], ['n', got.n, w.n],
        ['rs_number', got.rsNumber, w.rs_number], ['cs_row', got.csRow, w.cs_row], ['cs_id', got.csId, w.cs_id]] as const)
        if (firstDiff(a as ArrayLike<unknown>, b as unknown[]) >= 0) fail(k)
      if (relDiff(got.csPip, w.cs_pip).max !== 0) fail('cs_pip')
      // the SE is exp() of its code, and V8's exp can sit one ulp from numpy's; the rest are exact
      for (const [k, a, b] of [['beta', got.beta, w.beta], ['se', got.se, w.se], ['af', got.af, w.af], ['p', got.p, w.p]] as const) {
        const d = relDiff(a, b)
        if (!d.placed) fail(`${k} placement`)
        if (k === 'se') worstSe = Math.max(worstSe, d.max); else worst = Math.max(worst, d.max)
      }
      rows += got.rows
    }
    check(bad.length === 0 && worst === 0 && worstSe <= 1e-15, `${G.windows.length} GWAS windows (${rows} rows, the smallest chromosome whole): positions, alleles, rsID, n, credible sets exact; beta, af, p bit-identical to read_window, SE within ${g(worstSe)}${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : worst ? ` (beta/af/p off by ${g(worst)})` : ''}`)
    // every chromosome whole, summed
    const badTotals: string[] = []
    let all = 0
    for (const [chr, t] of Object.entries(G.chrom_totals as Record<string, { rows: number; pos: number; n: number; rs_number: number; abs_beta: number; se: number; af: number; af_not_finite: number; nlp: number; nlp_not_finite: number; n_cs: number }>)) {
      const got = read(chr, 1, 0x7fffffff)
      if (!got) { badTotals.push(`${chr}: no range`); continue }
      const af = fsum(got.af), nlp = fsum(got.nlp), beta = fsum(got.beta, Math.abs), se = fsum(got.se)
      if (got.rows !== t.rows || isum(got.position) !== t.pos || isum(got.n) !== t.n || isum(got.rsNumber) !== t.rs_number || got.csRow.length !== t.n_cs
        || af.bad !== t.af_not_finite || nlp.bad !== t.nlp_not_finite || !near(af.sum, t.af) || !near(nlp.sum, t.nlp) || !near(beta.sum, t.abs_beta) || !near(se.sum, t.se))
        badTotals.push(chr)
      all += got.rows
    }
    check(badTotals.length === 0 && all === gw.n_rows, `GWAS ${G.experiment}: all ${all} rows of ${Object.keys(G.chrom_totals).length} chromosomes sum to read_window's totals (n_rows ${gw.n_rows})${badTotals.length ? `: ${badTotals.join(', ')}` : ''}`)
    // the bin summary
    const b = gw.bins ?? null
    if (!b || !ref.gwas_bins) check(!b && !ref.gwas_bins, 'no GWAS bin summary (pointer and reference agree)')
    else {
      const got = decodeGwasBins(whole(b.file), b, b.file)
      const want = ref.gwas_bins.rows as Record<string, unknown>[]
      const badBins = got.length !== want.length ? [`${got.length} rows != ${want.length}`]
        : got.flatMap((r, i) => Object.entries(want[i]).filter(([k, v]) => !Object.is((r as unknown as Record<string, unknown>)[k], v)).map(([k]) => `row ${i} ${k}`))
      const chroms = [...new Set(got.map(r => r.chr))]
      check(badBins.length === 0 && want.length > 0, `GWAS bins ${b.file}: ${got.length} bins of ${b.bin_bp / 1e6} Mb on ${chroms.length} chromosomes, every field equal to read_bins${badBins.length ? `: ${badBins.slice(0, 4).join('; ')}` : ''}`)
      let threw = false
      try { decodeGwasBins(whole(b.file), { ...b, n_bins: b.n_bins + 1 }, b.file) } catch { threw = true }
      check(threw, 'decodeGwasBins refuses a bin count the pointer does not match')
    }
  }
}

// ---- result blocks ----------------------------------------------------------------------------
{
  let exact = 0
  const bad: string[] = []
  const worst = { nlp: 0, pval: 0, beta: 0, se: 0 }
  let rows = 0, cs = 0, scanRows = 0, scanBad = 0, withCs = 0
  const files = new Set<string>()
  for (const b of ref.blocks) {
    const r = b.row
    const what = `${b.experiment} ${r.phenotype_id}`
    if (!files.has(b.file)) {
      files.add(b.file)
      const h = parseFileHeader(readRange(b.file, 0, HEADER_LEN), b.file)
      const c = catalog.chromosomes.find((x: { name: string }) => x.name === r.chr)
      checkFileHeader(h, { kind: KIND.results, chrom: r.chr, seqDigest: c.seq_digest }, b.file)
    }
    const bytes = readRange(b.file, r.blk_off, r.blk_len)
    const t = decodeResultBlock(bytes, { blk_len: r.blk_len, n_var: r.n_var, var_start: r.var_start }, what)
    const fail = (k: string) => bad.push(`${what}: ${k}`)
    if (t.nRows !== b.n_rows || t.varStart !== b.var_start || t.posFirst !== b.pos_first || t.posLast !== b.pos_last) fail('header')
    if (t.nlpMax !== py(b.nlp_max) || t.betaMax !== py(b.beta_max) || t.lseMin !== py(b.lse_min) || t.lseMax !== py(b.lse_max)) fail('scales')
    if (firstDiff(t.nlpCode, b.nlp_code) >= 0 || firstDiff(t.betaCode, b.beta_code) >= 0 || firstDiff(t.seCode, b.se_code) >= 0) fail('codes')
    if (t.csRow.length !== b.n_cs || firstDiff(t.csRow, b.cs_row) >= 0 || firstDiff(t.csId, b.cs_id) >= 0 || relDiff(t.csPip, b.cs_pip).max !== 0) fail('credible sets')
    if (JSON.stringify(t.details) !== JSON.stringify(b.details)) fail('details')
    for (const k of ['nlp', 'pval', 'beta', 'se'] as const) {
      const d = relDiff(t[k], b[k])
      if (!d.placed) fail(`${k} null/infinite placement`)
      worst[k] = Math.max(worst[k], d.max)
      exact += d.exact
    }
    rows += t.nRows; cs += t.csRow.length; if (t.csRow.length) withCs++
    // the variant run and the reader columns
    if (t.nRows && r.var_off != null) {
      const range = decodeVariantRange(readRange(ref.chroms[r.chr].file, r.var_off, r.var_len), r.var_len, `${what} variants`)
      const run = sliceRun(range, t.varStart!, t.nRows, t, what)
      const cols = readerColumns(t, run, 1000)
      if (cols.tssDistance[0] !== run.position[0] - 1000 || cols.a1[0] !== run.alt[0] || cols.beta[0] !== Math.fround(t.beta[0]) && !(Number.isNaN(t.beta[0]) && Number.isNaN(cols.beta[0]))) fail('reader columns')
      if (run.position[0] < t.posFirst || run.position[t.nRows - 1] > t.posLast) fail('run positions outside the block')
    }
    // the cis scan's one-row read, row by row
    for (let i = 0; i < t.nRows; i += Math.max(1, Math.floor(t.nRows / 50))) {
      const s = scanBlockRow(bytes, i, what)
      scanRows++
      const same = (x: number, y: number) => Object.is(x, y) || (Number.isNaN(x) && Number.isNaN(y))
      if (!same(s.pval, t.pval[i]) || !same(s.se, t.se[i]) || !same(s.beta, t.beta[i])) scanBad++
    }
  }
  check(bad.length === 0, `${ref.blocks.length} blocks of ${cisIds.length} cis experiments (${withCs} with credible sets, ${rows} rows, ${cs} credible-set records): header, scales, codes, credible sets, details, runs equal read_block${bad.length ? `: ${bad.slice(0, 5).join('; ')}` : ''}`)
  check(worst.nlp <= 1e-15 && worst.pval <= 1e-13 && worst.beta <= 1e-15 && worst.se <= 1e-13,
    `values: ${exact} of ${4 * rows} -log10 p, p, beta and SE values bit-identical; largest relative difference -log10 p ${g(worst.nlp)}, p ${g(worst.pval)}, beta ${g(worst.beta)}, SE ${g(worst.se)}`)
  check(scanBad === 0, `scanBlockRow equals decodeResultBlock on ${scanRows} sampled rows`)
}

// ---- the gene lookup, per-chromosome genes and exon models, search index parts -------------------
{
  // [key as given, normalized, bucket of 1024]: pipeline/test_annotation.py LOOKUP_VECTORS
  const vectors: [string, string, number][] = [['', '', 453], ['a', 'A', 716], ['FLNC', 'FLNC', 208], ['flnc', 'FLNC', 208],
    ['ENSG00000128591', 'ENSG00000128591', 414], ['HLA-DRB1', 'HLA-DRB1', 280], ['Y_RNA', 'Y_RNA', 828], ['Ä', 'Ä', 834]]
  const off = vectors.filter(([raw, key, b]) => lookupKey(raw) !== key || lookupBucket(key, 1024) !== b)
  check(off.length === 0 && lookupBucket('A', 2 ** 32) === 0xc40bf6cc, `lookupKey and lookupBucket give the pinned vectors (FNV-1a 32)${off.length ? `: ${JSON.stringify(off)}` : ''}`)

  const L = ref.lookup
  const dir = decodeLookupDir(readRange(L.file, 0, lookupDirLen(L.n_buckets)), L.identity, L.file)
  const rowsFor = (raw: string) => {
    const key = lookupKey(raw), b = lookupBucket(key, dir.nBuckets)
    const len = dir.off[b + 1] - dir.off[b]
    if (!len) return []
    return tableFromIPC(decodeArrowObject(readRange(L.file, dir.off[b], len), `${L.file} bucket ${b}`)).toArray()
      .map(r => r.toJSON()).filter(r => r.key === key).map(r => ({ key: r.key, gene_id: r.gene_id, name: r.name, chr: r.chr, tss: Number(r.tss) }))
  }
  const bad = Object.entries(L.keys as Record<string, unknown[]>).filter(([k, want]) => JSON.stringify(rowsFor(k)) !== JSON.stringify(want))
  const misses = Object.values(L.keys as Record<string, unknown[]>).filter(x => !x.length).length
  check(dir.nBuckets === L.n_buckets && bad.length === 0, `gene lookup: ${dir.nBuckets} buckets; ${Object.keys(L.keys).length} keys (${misses} misses) give the rows Python gives${bad.length ? `: ${bad.slice(0, 3).map(x => x[0]).join(', ')}` : ''}`)

  let nGenes = 0, nModels = 0
  const badChrom: string[] = []
  for (const [c, want] of Object.entries(ref.annotation_chroms as Record<string, { genes: string; exon_models: string; n: number; gene_id: string[]; tss: number[]; models: { gene_id: string; exon_starts: number[]; exon_ends: number[] }[] }>)) {
    const gt = tableFromIPC(decodeArrowObject(whole(want.genes), want.genes))
    const m = tableFromIPC(decodeArrowObject(whole(want.exon_models), want.exon_models))
    const ids = gt.getChild('gene_id')!.toArray() as string[], tss = Array.from(gt.getChild('tss')!.toArray() as Int32Array)
    const mid = m.getChild('gene_id')!.toArray() as string[]
    if (gt.numRows !== want.n || firstDiff(ids, want.gene_id) >= 0 || firstDiff(tss, want.tss) >= 0 || firstDiff(mid, want.gene_id) >= 0) badChrom.push(`${c} genes`)
    const byId = new Map(mid.map((x, i) => [x, i]))
    for (const w of want.models) {
      const i = byId.get(w.gene_id)!
      const st = Array.from(m.getChild('exon_starts')!.get(i)?.toArray() ?? []), en = Array.from(m.getChild('exon_ends')!.get(i)?.toArray() ?? [])
      if (firstDiff(st, w.exon_starts) >= 0 || firstDiff(en, w.exon_ends) >= 0) badChrom.push(`${c} ${w.gene_id} exons`)
      nModels++
    }
    nGenes += gt.numRows
  }
  check(badChrom.length === 0, `${Object.keys(ref.annotation_chroms).length} chromosomes: ${nGenes} genes rows and ${nModels} sampled exon models equal the Python decode${badChrom.length ? `: ${badChrom.slice(0, 3).join(', ')}` : ''}`)

  // every experiment's search index parts: rows, ords and runs, and the per-kind offsets; the parts
  // of one experiment together hold every phenotype once
  const badPart: string[] = []
  let nParts = 0, nRows = 0
  for (const e of expIds) {
    const doc = docs[e]
    const parts = ref.index_parts[e] as Record<string, { file: string; rows: number; ords: [number, number][]; ord: number[]; phenotype_id: string[]; chr: string[]
      blk_off?: (number | null)[]; trans_off?: (number | null)[]; trans_len?: (number | null)[]; n_trans?: (number | null)[]; gene_id?: (string | null)[] }>
    const seen = new Set<number>()
    for (const [c, want] of Object.entries(parts)) {
      const p = c === 'unplaced' ? doc.search_index_unplaced : doc.search_index_parts[c]
      if (!p || p.file !== want.file || p.rows !== want.rows) { badPart.push(`${e} ${c}: pointer`); continue }
      const t = tableFromIPC(decodeArrowObject(whole(want.file), want.file))
      const col = (k: string) => Array.from({ length: t.numRows }, (_, i) => t.getChild(k)!.get(i) ?? null)
      const ords = col('ord') as number[]
      const runs: [number, number][] = []
      for (const o of ords) { const r = runs[runs.length - 1]; if (r && o === r[1] + 1) r[1] = o; else runs.push([o, o]) }
      ords.forEach(o => seen.add(o))
      const keys = doc.kind === 'qtl_cis' ? ['blk_off'] : ['trans_off', 'trans_len', 'n_trans', 'gene_id']
      if (t.numRows !== want.rows || firstDiff(ords, want.ord) >= 0 || JSON.stringify(runs) !== JSON.stringify(want.ords)
        || firstDiff(col('phenotype_id'), want.phenotype_id) >= 0 || firstDiff(col('chr'), want.chr) >= 0
        || keys.some(k => firstDiff(col(k), (want as unknown as Record<string, unknown[]>)[k]) >= 0)) badPart.push(`${e} ${c}`)
      nRows += t.numRows; nParts++
    }
    if (seen.size !== doc.n_phenotypes || [...seen].some(o => o < 0 || o >= doc.n_phenotypes)) badPart.push(`${e}: parts hold ${seen.size} ords, n_phenotypes ${doc.n_phenotypes}`)
  }
  check(badPart.length === 0, `${nParts} search index parts of ${expIds.length} experiments: ${nRows} rows, ords, runs and offsets equal the Python decode; each experiment's parts hold its ordinals once${badPart.length ? `: ${badPart.slice(0, 4).join('; ')}` : ''}`)

  // the coloc lists (SPEC section 20): each pointer and object as Python reads them, and every listed
  // gene placed by the gene lookup the Home track uses, on a chromosome of the catalog
  const colocs = (ref.colocs ?? []) as { doc: Record<string, unknown>; rows: { gene_id: string; phenotype_id: string | null }[]; lookup: Record<string, unknown[]> }[]
  const badColoc: string[] = []
  const chromNames = new Set(Object.keys(ref.chroms))
  for (const c of colocs) {
    const id = c.doc.id as string
    const doc = pointer('colocs', id)
    const rows = JSON.parse(new TextDecoder().decode(whole(doc.object)))
    if (JSON.stringify(doc) !== JSON.stringify(c.doc) || JSON.stringify(rows) !== JSON.stringify(c.rows) || rows.length !== doc.n) badColoc.push(`${id}: pointer or rows`)
    for (const g of Object.keys(c.lookup)) {
      const got = rowsFor(g).filter(r => r.gene_id === g)
      if (JSON.stringify(rowsFor(g)) !== JSON.stringify(c.lookup[g]) || !got.length || !got.every(r => chromNames.has(r.chr))) badColoc.push(`${id}: ${g} placement`)
    }
  }
  check(badColoc.length === 0, `${colocs.length} coloc lists: pointers and rows equal Python's, ${colocs.reduce((a, c) => a + Object.keys(c.lookup).length, 0)} listed genes placed by the gene lookup${badColoc.length ? `: ${badColoc.slice(0, 4).join('; ')}` : ''}`)
}

console.log(failures.length ? `store-check: ${failures.length} failed (${sites} sites checked)` : `store-check: all checks passed (${sites} sites)`)
process.exit(failures.length ? 1 : 0)
