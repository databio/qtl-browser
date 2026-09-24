/**
 * Round trip of the browser decoder (src/lib/store-decode.ts) against the Python decoders
 * (pipeline/catalog.py, pipeline/results.py, pipeline/packfmt.py) on one experiment of a local store.
 *
 *   uv run python ui/scripts/store_reference.py <store> <ref.json> [--experiment topchef]
 *   npm run store-check -- <store> <ref.json>
 *
 * - variants files: every site of every chromosome (both sections, walked with the TS variant
 *   index), against `catalog.decode_file`; one page at a time through decodeVariantPage too;
 * - the variant index against `catalog.decode_vidx`, and every rsID record against `decode_rsid`,
 *   plus a first-of-run lookup (rsidFind) for a sample of rs numbers, repeats and block edges included;
 * - hits files, paged: the frame table against `packfmt_v1.hits_frame_table` and every frame's
 *   records against `results.decode_hits` (value and beta exact as f32), empty frames included;
 * - rsidFind against `catalog.rsid_lookup` (first record of the run) for a sample of rs numbers;
 * - every trans frame against `results.read_trans` (codes, positions, alleles exactly; p, beta and
 *   the derived SE to rounding), and the layout: the frames tile each trans object and each gene's
 *   frames are one contiguous range (what the gene page reads in one request);
 * - GWAS windows (each sampled block's run, and one whole chromosome) against `gwas.read_window`;
 * - every row of the GWAS bin summary against `gwas.read_bins`, every field exact;
 * - the sampled result blocks against `results.read_block`: codes, credible sets, details and
 *   positions exactly; -log10 p, p, SE and the rebuilt slope as the largest relative difference;
 *   scanBlockRow row by row against decodeResultBlock; and `dof` null leaving every slope NaN.
 * - what a first gene page reads (SPEC sections 6 and 8): the gene lookup's key normalization and
 *   bucket hash on pinned vectors (the same ones as pipeline/test_annotation.py), its directory and
 *   the rows for a sample of keys, misses included; every chromosome's genes and exon models; every
 *   search index part's rows and ord runs; and the Home page's coloc table (lib/coloc.ts) against
 *   the lookup when the store has the annotation the table was read from.
 * Integers, strings and codes must match exactly. Exits 1 on any failure.
 */
import { closeSync, openSync, readFileSync, readSync } from 'node:fs'
import { join } from 'node:path'
import { tableFromIPC } from 'apache-arrow'
import { COLOC_ANNOTATION, COLOC_LOCI } from '../src/lib/coloc.ts'
import { checkFileHeader, decodeArrowObject, decodeLookupDir, lookupBucket, lookupDirLen, lookupKey, decodeGwasBins, decodeGwasIndex, decodeGwasRange, decodeHitsFrame, decodeHitsTable, decodeResultBlock, decodeRsidRecords,
  decodeTransFrame, decodeVariantIndex, decodeVariantPage, decodeVariantRange, gwasRange, HEADER_LEN, hitsOf, hitsTableLen, KIND,
  parseFileHeader, readerColumns, RSID_RECORD_LEN, rsidFind, scanBlockRow, sliceRun, transSe,
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

const catalog = JSON.parse(readFileSync(join(STORE, ref.catalog_dir, `${ref.catalog}.json`), 'utf8'))
const names: string[] = catalog.chromosomes.map((c: { name: string }) => c.name)

// ---- variant index ----------------------------------------------------------------------------
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
  const diffs = Object.entries(cols).map(([k, v]) => {
    const w = k === 'af' ? (want.af as (number | null)[]).map(py) : want[k]
    const i = k === 'af' ? v.findIndex((x, j) => !(Object.is(x, w[j]) || (Number.isNaN(x) && Number.isNaN(w[j])))) : firstDiff(v, w)
    return i >= 0 || v.length !== w.length ? `${k} at ${i}` : null
  }).filter(Boolean)
  check(diffs.length === 0 && cols.pos.length === c.count,
    `${c.name}: all ${cols.pos.length} sites (pos, ref, alt, af code and value, rs_number, counts, match, alt_is_minor) equal decode_file${diffs.length ? `: ${diffs.join(', ')}` : ''}`)
  // pages one at a time, as the variant page reads them
  let pageBad = 0
  const nPages = ci.nPagesCis + ci.nPagesTrans
  for (let k = 0; k < nPages; k += Math.max(1, Math.floor(nPages / 25))) {
    const recs = decodeVariantPage(readRange(c.file, ci.pageOff[k], ci.pageOff[k + 1] - ci.pageOff[k]), c.name, ci.nCis)
    for (const r of recs)
      if (r.position !== want.pos[r.vidx] || r.A1 !== want.alt[r.vidx] || r.A2 !== want.ref[r.vidx] || r.inCis !== r.vidx < ci.nCis) pageBad++
  }
  check(pageBad === 0, `${c.name}: sampled single pages decode to the same sites (A1 = alt, A2 = ref, inCis by n_cis)`)

  // hits, paged: the table from the first bytes, then every frame on its own
  const hw = want.hits
  const hb = whole(hw.file)
  const tab = decodeHitsTable(hb.subarray(0, hitsTableLen(c.count, 1024)), { chrom: c.name, seqDigest: c.seq_digest, nVariants: c.count }, hw.file)
  const hc = { vidx: [] as number[], ord: [] as number[], value: [] as number[], beta: [] as number[], kind: [] as number[], cs_id: [] as number[], flags: [] as number[] }
  let empty = 0
  for (let g = 0; g + 1 < tab.frameOff.length; g++) {
    const fr = decodeHitsFrame(hb.subarray(tab.frameOff[g], tab.frameOff[g + 1]), g * tab.frameVariants, tab.frameVariants, `${hw.file} frame ${g}`)
    if (!fr.count) empty++
    for (let i = 0; i < fr.count; i++) {
      hc.vidx.push(fr.vidx[i]); hc.ord.push(fr.ord[i]); hc.value.push(fr.value[i]); hc.beta.push(fr.beta[i])
      hc.kind.push(fr.kind[i]); hc.cs_id.push(fr.csId[i]); hc.flags.push(fr.flags[i])
    }
    if (fr.count && g === Math.floor(tab.frameOff.length / 2)) {
      const v0 = fr.vidx[Math.floor(fr.count / 2)], [lo, hi] = hitsOf(fr, v0)
      check(hi > lo && Array.from(fr.vidx.subarray(lo, hi)).every(x => x === v0) && (lo === 0 || fr.vidx[lo - 1] < v0) && (hi === fr.count || fr.vidx[hi] > v0),
        `${c.name}: hitsOf(vidx ${v0}) finds its ${hi - lo} records in frame ${g}`)
    }
  }
  const hd = (['vidx', 'ord', 'kind', 'cs_id', 'flags'] as const).filter(k => firstDiff(hc[k], hw[k]) >= 0)
  const hv = relDiff(hc.value, hw.value), hbeta = relDiff(hc.beta, hw.beta)
  const kinds = [0, 1, 2].map(k => hc.kind.filter(x => x === k).length)
  check(firstDiff(tab.frameOff, hw.frame_off) < 0 && tab.count === hw.header.count && tab.frameVariants === hw.header.page_size && hc.vidx.length === tab.count,
    `${c.name} hits ${hw.file}: header (count ${tab.count}, ${tab.frameVariants} variants per frame, ${tab.nVariants} variants) and ${tab.frameOff.length - 1}-frame table equal hits_frame_table (${empty} empty frames)`)
  check(hd.length === 0 && hv.max === 0 && hv.placed && hbeta.max === 0 && hbeta.placed,
    `${c.name} hits: ${hc.vidx.length} records (leads ${kinds[0]}, credible sets ${kinds[1]}, trans ${kinds[2]}) equal decode_hits, value and beta exact${hd.length ? `; differ: ${hd.join(', ')}` : ''}`)
}

// ---- rsID index -------------------------------------------------------------------------------
{
  const R = ref.rsid
  const all = decodeRsidRecords(readRange(R.file, HEADER_LEN, R.rs_number.length * RSID_RECORD_LEN), R.file)
  const h = parseFileHeader(readRange(R.file, 0, HEADER_LEN), R.file)
  checkFileHeader(h, { kind: KIND.rsid, chrom: 'all', seqDigest: catalog.collection_digest }, R.file)
  check(firstDiff(all.map(r => r.rsNumber), R.rs_number) < 0 && firstDiff(all.map(r => r.vidx), R.vidx) < 0 && firstDiff(all.map(r => r.ordinal), R.ordinal) < 0 && h.count === all.length,
    `rsID index ${R.file}: all ${all.length} records equal decode_rsid, header kind 8 on the collection digest`)
  const B = vx.rsidBlockRecords
  const block = (b: number): RsidRecord[] => all.slice(b * B, (b + 1) * B)
  // rs numbers to try: a spread, every repeated one, both sides of each block edge, and misses
  const first = new Map<number, RsidRecord>()
  for (const r of all) if (!first.has(r.rsNumber)) first.set(r.rsNumber, r)
  const repeats = all.filter((r, i) => i > 0 && all[i - 1].rsNumber === r.rsNumber).map(r => r.rsNumber)
  const edges = Array.from(vx.rsidFirst).flatMap((x, b) => [x, all[Math.max(0, b * B - 1)].rsNumber])
  const tries = [...new Set([...all.filter((_, i) => i % 997 === 0).map(r => r.rsNumber), ...repeats, ...edges, 1, 2, all[all.length - 1].rsNumber + 1,
    all[0].rsNumber - 1, all[5].rsNumber + 1])]
  let bad = 0
  for (const rs of tries) {
    const got = await rsidFind(vx.rsidFirst, rs, block)
    const want = first.get(rs) ?? null
    if (!(got === want || (got && want && got.vidx === want.vidx && got.ordinal === want.ordinal))) bad++
  }
  check(bad === 0, `rsidFind: ${tries.length} rs numbers (${new Set(repeats).size} repeated, block edges, misses) find the first record of their run (${bad} wrong)`)
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

// rsidFind against catalog.rsid_lookup on the real index
{
  const R = ref.rsid
  const all = decodeRsidRecords(readRange(R.file, HEADER_LEN, R.rs_number.length * RSID_RECORD_LEN), R.file)
  const B = vx.rsidBlockRecords
  let bad = 0, n = 0
  for (const [rs, want] of Object.entries(ref.rsid_lookup) as [string, [number, number][]][]) {
    n++
    const got = await rsidFind(vx.rsidFirst, Number(rs), b => all.slice(b * B, (b + 1) * B))
    if (!(want.length ? got && got.ordinal === want[0][0] && got.vidx === want[0][1] : got === null)) bad++
  }
  check(bad === 0, `rsidFind equals catalog.rsid_lookup's first record on ${n} rs numbers (${bad} wrong)`)
}

// ---- trans frames -----------------------------------------------------------------------------
{
  let bad: string[] = [], rows = 0
  const worst = { p: 0, beta: 0, se: 0 }
  const headers = new Set<string>()
  for (const t of ref.trans) {
    if (!headers.has(t.file)) {
      headers.add(t.file)
      checkFileHeader(parseFileHeader(readRange(t.file, 0, HEADER_LEN), t.file), { kind: KIND.trans, chrom: 'all', seqDigest: catalog.collection_digest }, t.file)
    }
    const fr = decodeTransFrame(readRange(t.file, t.trans_off, t.trans_len), `${t.file} ${t.phenotype_id}`)
    const fail = (k: string) => bad.push(`${t.phenotype_id}: ${k}`)
    if (fr.n !== t.n_trans || fr.n !== t.pos.length) fail('row count')
    if (fr.nlpMax !== py(t.nlp_max) || fr.betaMax !== py(t.beta_max)) fail('scales')
    for (const [k, a, b] of [['pos', fr.position, t.pos], ['rs_number', fr.rsNumber, t.rs_number],
      ['af_code', fr.afCode, t.af_code], ['nlp_code', fr.nlpCode, t.nlp_code], ['ref', fr.ref, t.ref], ['alt', fr.alt, t.alt]] as const)
      if (firstDiff(a as ArrayLike<unknown>, b as unknown[]) >= 0) fail(k)
    if (firstDiff(Array.from(fr.ordinal, o => names[o - 1]), t.chr) >= 0) fail('chromosome')
    if (relDiff(fr.af, t.af).max !== 0 || relDiff(fr.nlp, t.nlp).max !== 0) fail('af or -log10 p')
    const p = relDiff(fr.pval, t.p), be = relDiff(fr.beta, t.beta)
    const se = relDiff(Array.from(fr.nlp, (x, i) => transSe(x, fr.beta[i], t.dof).se), t.se)
    if (!p.placed || !be.placed || !se.placed) fail('null or infinite placement')
    worst.p = Math.max(worst.p, p.max); worst.beta = Math.max(worst.beta, be.max); worst.se = Math.max(worst.se, se.max)
    rows += fr.n
  }
  check(bad.length === 0, `${ref.trans.length} trans frames, ${rows} rows: positions, chromosomes, rsID, af, codes and alleles equal read_trans${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : ''}`)
  check(worst.p <= 1e-13 && worst.beta <= 1e-13 && worst.se <= 1e-9,
    `trans values: largest relative difference p ${g(worst.p)}, beta ${g(worst.beta)}, SE (TS inverse t against scipy) ${g(worst.se)}`)
  // layout: per object the frames tile it from byte 64; per (object, gene) the frames are back to back
  const byFile = new Map<string, { off: number; len: number; gene: string | null }[]>()
  for (const t of ref.trans) {
    if (!byFile.has(t.file)) byFile.set(t.file, [])
    byFile.get(t.file)!.push({ off: t.trans_off, len: t.trans_len, gene: t.gene_id })
  }
  const layoutBad: string[] = []
  let genes = 0, maxFrames = 0
  for (const [file, fs] of byFile) {
    fs.sort((a, b) => a.off - b.off)
    let end = HEADER_LEN
    for (const f of fs) { if (f.off !== end) { layoutBad.push(`${file}: gap or overlap at ${f.off}`); break } end = f.off + f.len }
    if (end !== readFileSync(imm(file)).length) layoutBad.push(`${file}: frames end at ${end}, not at the end of the object`)
    const byGene = new Map<string, { off: number; len: number }[]>()
    for (const f of fs) if (f.gene !== null) { if (!byGene.has(f.gene)) byGene.set(f.gene, []); byGene.get(f.gene)!.push(f) }
    for (const [gene, xs] of byGene) {
      genes++; maxFrames = Math.max(maxFrames, xs.length)
      if (xs.some((x, i) => i > 0 && x.off !== xs[i - 1].off + xs[i - 1].len)) layoutBad.push(`${file}: ${gene} frames not contiguous`)
    }
  }
  check(layoutBad.length === 0, `trans layout: ${byFile.size} objects tiled by their frames; ${genes} (object, gene) pairs each one contiguous range (up to ${maxFrames} frames)${layoutBad.length ? `: ${layoutBad.slice(0, 4).join('; ')}` : ''}`)
}

// ---- GWAS -------------------------------------------------------------------------------------
{
  const exp = JSON.parse(readFileSync(join(STORE, 'experiments', `${ref.experiment}.json`), 'utf8'))
  const gw = exp.gwas
  if (!gw) check(ref.gwas.length === 0, 'no GWAS in this experiment')
  else {
    const index = decodeGwasIndex(whole(gw.index), catalog.collection_digest, gw.index)
    let bad: string[] = [], rows = 0, worst = 0
    for (const w of ref.gwas) {
      const file = gw.files[w.chr]
      const h = parseFileHeader(readRange(file, 0, HEADER_LEN), file)
      checkFileHeader(h, { kind: KIND.gwas, chrom: w.chr, seqDigest: catalog.chromosomes.find((c: { name: string }) => c.name === w.chr).seq_digest }, file)
      const range = gwasRange(index, w.chr, w.lo, w.hi)
      const got = range ? decodeGwasRange(readRange(file, range.off, range.len), index, w.chr, range, w.lo, w.hi) : null
      const fail = (k: string) => bad.push(`${w.chr}:${w.lo}-${w.hi}: ${k}`)
      if ((got?.rows ?? 0) !== w.pos.length) { fail(`rows ${got?.rows ?? 0} != ${w.pos.length}`); continue }
      if (!got) continue
      for (const [k, a, b] of [['pos', got.position, w.pos], ['ref', got.ref, w.ref], ['alt', got.alt, w.alt], ['n', got.n, w.n],
        ['rs_number', got.rsNumber, w.rs_number]] as const)
        if (firstDiff(a as ArrayLike<unknown>, b as unknown[]) >= 0) fail(k)
      for (const [k, a, b] of [['beta', got.beta, w.beta], ['se', got.se, w.se], ['af', got.af, w.af], ['p', got.p, w.p]] as const) {
        const d = relDiff(a, b)
        if (!d.placed) fail(`${k} placement`)
        worst = Math.max(worst, d.max)
      }
      rows += got.rows
    }
    check(bad.length === 0 && worst === 0, `${ref.gwas.length} GWAS windows (${rows} rows, one whole chromosome): positions, alleles, rsID, n exact; beta, se, af, p bit-identical to read_window${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : ` (worst ${g(worst)})`}`)
  }
}

// ---- GWAS bin summary -------------------------------------------------------------------------
{
  const exp = JSON.parse(readFileSync(join(STORE, 'experiments', `${ref.experiment}.json`), 'utf8'))
  const b = exp.gwas?.bins ?? null
  if (!b || !ref.gwas_bins) check(!b && !ref.gwas_bins, 'no GWAS bin summary in this experiment (pointer and reference agree)')
  else {
    const got = decodeGwasBins(whole(b.file), b, b.file)
    const want = ref.gwas_bins.rows as Record<string, unknown>[]
    const bad = got.length !== want.length ? [`${got.length} rows != ${want.length}`]
      : got.flatMap((r, i) => Object.entries(want[i]).filter(([k, v]) => !Object.is((r as unknown as Record<string, unknown>)[k], v)).map(([k]) => `row ${i} ${k}`))
    const chroms = [...new Set(got.map(r => r.chr))]
    check(bad.length === 0 && want.length > 0, `GWAS bins ${b.file}: ${got.length} bins of ${b.bin_bp / 1e6} Mb on ${chroms.join(', ')}, every field equal to read_bins${bad.length ? `: ${bad.slice(0, 4).join('; ')}` : ''}`)
    let threw = false
    try { decodeGwasBins(whole(b.file), { ...b, n_bins: b.n_bins + 1 }, b.file) } catch { threw = true }
    check(threw, 'decodeGwasBins refuses a bin count the pointer does not match')
  }
}

// ---- result blocks ----------------------------------------------------------------------------
{
  let exact = 0, bad: string[] = []
  const worst = { nlp: 0, pval: 0, se: 0, slope: 0 }
  let rows = 0, cs = 0, scanRows = 0, scanBad = 0, withCs = 0
  const files = new Set<string>()
  for (const b of ref.blocks) {
    const r = b.row
    const what = `${r.phenotype_type} ${r.phenotype_id}`
    if (!files.has(b.file)) {
      files.add(b.file)
      const h = parseFileHeader(readRange(b.file, 0, HEADER_LEN), b.file)
      const c = catalog.chromosomes.find((x: { name: string }) => x.name === r.chr)
      checkFileHeader(h, { kind: KIND.results, chrom: r.chr, seqDigest: c.seq_digest }, b.file)
    }
    const bytes = readRange(b.file, r.blk_off, r.blk_len)
    const t = decodeResultBlock(bytes, b.dof, { blk_len: r.blk_len, n_var: r.n_var, var_start: r.var_start }, what)
    const fail = (k: string) => bad.push(`${what}: ${k}`)
    if (t.nRows !== b.n_rows || t.varStart !== b.var_start || t.posFirst !== b.pos_first || t.posLast !== b.pos_last) fail('header')
    if (t.nlpMax !== py(b.nlp_max) || t.lseMin !== py(b.lse_min) || t.lseMax !== py(b.lse_max)) fail('scales')
    if (firstDiff(t.nlpCode, b.nlp_code) >= 0 || firstDiff(t.seCode, b.se_code) >= 0) fail('codes')
    if (firstDiff(t.csRow, b.cs_row) >= 0 || firstDiff(t.csId, b.cs_id) >= 0 || relDiff(t.csPip, b.cs_pip).max !== 0) fail('credible sets')
    if (JSON.stringify(t.details) !== JSON.stringify(b.details)) fail('details')
    for (const k of ['nlp', 'pval', 'se', 'slope'] as const) {
      const d = relDiff(t[k], b[k])
      if (!d.placed) fail(`${k} null/infinite placement`)
      worst[k] = Math.max(worst[k], d.max)
      if (k !== 'slope') exact += d.exact
    }
    rows += t.nRows; cs += t.csRow.length; if (t.csRow.length) withCs++
    // the variant run and the reader columns
    if (t.nRows && r.var_off != null) {
      const range = decodeVariantRange(readRange(ref.chroms[r.chr].file, r.var_off, r.var_len), r.var_len, `${what} variants`)
      const run = sliceRun(range, t.varStart!, t.nRows, t, what)
      const want = ref.chroms[r.chr]
      if (firstDiff(run.position, want.pos.slice(t.varStart!, t.varStart! + t.nRows)) >= 0) fail('run positions')
      const cols = readerColumns(t, run, 1000)
      if (cols.a1[0] !== want.alt[t.varStart!] || cols.tssDistance[0] !== run.position[0] - 1000) fail('reader columns')
    }
    // the cis scan's one-row read, row by row
    for (let i = 0; i < t.nRows; i += Math.max(1, Math.floor(t.nRows / 50))) {
      const s = scanBlockRow(bytes, i, b.dof, what)
      scanRows++
      const same = (x: number, y: number) => Object.is(x, y) || (Number.isNaN(x) && Number.isNaN(y))
      if (!same(s.pval, t.pval[i]) || !same(s.se, t.se[i]) || !same(s.slope, t.slope[i])) scanBad++
    }
    // dof null: the stored -log10 p and SE stand, and no slope is rebuilt
    const n = decodeResultBlock(bytes, null, { blk_len: r.blk_len }, what)
    if (!n.slope.every(Number.isNaN) || firstDiff(n.pval, t.pval) >= 0 || firstDiff(n.se, t.se) >= 0) fail('dof null')
  }
  check(bad.length === 0, `${ref.blocks.length} blocks (${withCs} with credible sets, ${rows} rows, ${cs} credible-set records): header, scales, codes, credible sets, details, run positions equal read_block${bad.length ? `: ${bad.slice(0, 5).join('; ')}` : ''}`)
  check(worst.nlp <= 1e-15 && worst.pval <= 1e-13 && worst.se <= 1e-13,
    `values: ${exact} of ${3 * rows} -log10 p, p and SE values bit-identical; largest relative difference -log10 p ${g(worst.nlp)}, p ${g(worst.pval)}, SE ${g(worst.se)}`)
  check(worst.slope <= 1e-9, `slope (TS inverse t against scipy stdtrit): largest relative difference ${g(worst.slope)}`)
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
    const g = tableFromIPC(decodeArrowObject(whole(want.genes), want.genes))
    const m = tableFromIPC(decodeArrowObject(whole(want.exon_models), want.exon_models))
    const ids = g.getChild('gene_id')!.toArray() as string[], tss = Array.from(g.getChild('tss')!.toArray() as Int32Array)
    const mid = m.getChild('gene_id')!.toArray() as string[]
    if (g.numRows !== want.n || firstDiff(ids, want.gene_id) >= 0 || firstDiff(tss, want.tss) >= 0 || firstDiff(mid, want.gene_id) >= 0) badChrom.push(`${c} genes`)
    const byId = new Map(mid.map((x, i) => [x, i]))
    for (const w of want.models) {
      const i = byId.get(w.gene_id)!
      const st = Array.from(m.getChild('exon_starts')!.get(i)?.toArray() ?? []), en = Array.from(m.getChild('exon_ends')!.get(i)?.toArray() ?? [])
      if (firstDiff(st, w.exon_starts) >= 0 || firstDiff(en, w.exon_ends) >= 0) badChrom.push(`${c} ${w.gene_id} exons`)
      nModels++
    }
    nGenes += g.numRows
  }
  check(badChrom.length === 0, `${Object.keys(ref.annotation_chroms).length} chromosomes: ${nGenes} genes rows and ${nModels} sampled exon models equal the Python decode${badChrom.length ? `: ${badChrom.slice(0, 3).join(', ')}` : ''}`)

  const badPart: string[] = []
  let nRows = 0
  for (const [c, want] of Object.entries(ref.index_parts as Record<string, { file: string; rows: number; ords: [number, number][]; ord: number[]; phenotype_id: string[]; blk_off: (number | null)[]; n_trans: (number | null)[] }>)) {
    const t = tableFromIPC(decodeArrowObject(whole(want.file), want.file))
    const col = (k: string) => Array.from({ length: t.numRows }, (_, i) => t.getChild(k)!.get(i) ?? null)
    const ords = col('ord') as number[]
    const runs: [number, number][] = []
    for (const o of ords) { const r = runs[runs.length - 1]; if (r && o === r[1] + 1) r[1] = o; else runs.push([o, o]) }
    if (t.numRows !== want.rows || firstDiff(ords, want.ord) >= 0 || JSON.stringify(runs) !== JSON.stringify(want.ords)
      || firstDiff(col('phenotype_id'), want.phenotype_id) >= 0 || firstDiff(col('blk_off'), want.blk_off) >= 0 || firstDiff(col('n_trans'), want.n_trans) >= 0) badPart.push(c)
    nRows += t.numRows
  }
  check(badPart.length === 0, `${Object.keys(ref.index_parts).length} search index parts: ${nRows} rows, ords and runs equal the Python decode${badPart.length ? `: ${badPart.join(', ')}` : ''}`)

  if (L.identity === COLOC_ANNOTATION) {
    const off = Object.entries(COLOC_LOCI).filter(([sym, loci]) =>
      JSON.stringify(rowsFor(sym).filter(r => r.name === sym).map(r => ({ gene_id: r.gene_id, chr: r.chr, tss: r.tss }))) !== JSON.stringify(loci))
    check(off.length === 0, `coloc table (lib/coloc.ts): ${Object.keys(COLOC_LOCI).length} symbols sit where the store's lookup puts them${off.length ? `: ${off.map(x => x[0]).join(', ')}` : ''}`)
  } else {
    check(false, `coloc table (lib/coloc.ts) is pinned to annotation ${COLOC_ANNOTATION}; this store's is ${L.identity}: re-read COLOC_LOCI from it`)
  }
}

console.log(failures.length ? `store-check: ${failures.length} failed` : 'store-check: all checks passed')
process.exit(failures.length ? 1 : 0)
