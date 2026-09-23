// Smoke check of the trans tables on the qtlb v1 store (SPEC.md sections 8 and 9), against the local
// preview serving the chr21/chr22 smoke store. Gene page: both tabs page off the frames of the gene's
// phenotypes in the trans objects, one range request per object (a gene's frames are contiguous); the rows exported as CSV must equal `results.read_trans` for the
// same phenotypes (written by `ui/scripts/store_reference.py`, path in SMOKE_REF). Variant page: the
// hits frame's kind 2 records.
//   SMOKE_REF=/tmp/ref.json node bench/smoke_trans_tab.mjs
import { chromium } from '../node_modules/playwright/index.mjs'
import { readFileSync } from 'node:fs'

const BASE = process.env.SMOKE_BASE ?? 'http://localhost:4173'
const REF = process.env.SMOKE_REF ? JSON.parse(readFileSync(process.env.SMOKE_REF, 'utf8')) : null
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }
const exp = await (await fetch(`${BASE}/data/experiments/topchef.json`)).json()
const TRANS = new Set(exp.results.map(r => r.trans?.file).filter(Boolean))
check(TRANS.size === 2, `experiment ${exp.id} has ${TRANS.size} trans objects (ge, leafcutter)`)

const browser = await chromium.launch({ headless: true })
async function open(path) {
  const ctx = await browser.newContext({ viewport: { width: 1280, height: 900 }, acceptDownloads: true })
  const page = await ctx.newPage()
  const errors = [], reqs = []
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()) })
  page.on('pageerror', e => errors.push(String(e)))
  page.on('request', r => { const u = new URL(r.url()); if (u.pathname.startsWith('/data/')) reqs.push(u.pathname) })
  await page.goto(BASE + path)
  return { ctx, page, errors, reqs }
}
const table = page => page.locator('[data-trans-total]')
async function csv(page) {
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 60_000 }), table(page).locator('button:has-text("CSV")').first().click()])
  const text = readFileSync(await dl.path(), 'utf8').trim().split('\n')
  const header = text[0].split(',')
  return { header, rows: text.slice(1).map(line => Object.fromEntries(line.split(',').map((v, i) => [header[i], v]))) }
}
/** Exported rows against read_trans rows of the gene's phenotypes of one type. */
function compare(got, ids, label) {
  if (!REF) return check(true, `${label}: ${got.length} rows (no SMOKE_REF: values not compared)`)
  const want = REF.trans.filter(t => ids.has(t.phenotype_id)).flatMap(t => t.pos.map((p, i) => ({ id: t.phenotype_id, chr: t.chr[i], pos: p,
    rs: t.rs_number[i] ? `rs${t.rs_number[i]}` : '', p: t.p[i], beta: t.beta[i], se: t.se[i] })))
  const key = r => `${r.id}|${r.chr}:${r.pos}`
  const byKey = new Map(want.map(r => [key(r), r]))
  let bad = 0, worst = { p: 0, beta: 0, se: 0 }
  for (const r of got) {
    const w = byKey.get(key({ id: r.phenotype_id ?? [...ids][0], chr: r.variant_chr, pos: Number(r.position) }))
    if (!w || (r.rsid || '') !== w.rs) { bad++; continue }
    const rel = (a, b) => Math.abs(Number(a) - b) / Math.max(Math.abs(b), 1e-300)
    worst.p = Math.max(worst.p, rel(r.pval, w.p)); worst.beta = Math.max(worst.beta, rel(r.beta, w.beta)); worst.se = Math.max(worst.se, rel(r.beta_se, w.se))
  }
  // the CSV holds floats as printed (beta and SE as f32), so compare to f32 precision
  check(bad === 0 && got.length === want.length && worst.p < 1e-12 && worst.beta < 1e-6 && worst.se < 1e-6,
    `${label}: ${got.length} exported rows equal read_trans (${want.length}); worst relative p ${worst.p.toPrecision(2)}, beta ${worst.beta.toPrecision(2)}, SE ${worst.se.toPrecision(2)}`)
}

// ---- SMARCB1: 7 trans eQTL rows, 34 trans sQTL rows ----
{
  const { ctx, page, errors, reqs } = await open('/gene/SMARCB1')
  await page.waitForSelector('[data-trans-total="7"]', { timeout: 60_000 })
  const n = reqs.filter(p => TRANS.has(p.split('/').pop())).length
  check(n === 4, `SMARCB1: ${n} trans object requests (per phenotype type: the header and one range for the gene's frames)`)
  const e = await csv(page)
  check(e.header.join(',') === 'variant_chr,position,rsid,af,pval,beta,beta_se,r2', `SMARCB1: trans eQTL CSV header ${e.header.join(',')}`)
  compare(e.rows, new Set(['ENSG00000099956']), 'SMARCB1 trans eQTL CSV')
  const first = await table(page).locator('tbody tr td').first().innerText()
  await table(page).locator('th button:text-is("p")').click()
  await page.waitForTimeout(400)
  check(first !== await table(page).locator('tbody tr td').first().innerText(), 'SMARCB1: sorting the trans table by p flips the first row')
  await table(page).locator('th button:text-is("p")').click()
  await page.getByText(/^sQTL/).first().click()
  await page.waitForSelector('[data-trans-total="34"]', { timeout: 60_000 })
  const s = await csv(page)
  const ids = new Set(s.rows.map(r => r.phenotype_id))
  check(s.header[0] === 'phenotype_id' && [...ids].every(x => x.endsWith('ENSG00000099956.20') || x.includes('ENSG00000099956')), `SMARCB1: trans sQTL CSV names ${ids.size} introns of the gene`)
  compare(s.rows, ids, 'SMARCB1 trans sQTL CSV')
  check(errors.length === 0, `SMARCB1: no console errors${errors.length ? ` (${errors.slice(0, 2).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- MICAL3: 140 trans sQTL rows across its introns, no trans eQTL rows ----
{
  const { ctx, page, errors, reqs } = await open('/gene/MICAL3?tab=sqtl')
  await page.waitForSelector('[data-trans-total="140"]', { timeout: 60_000 })
  const n = reqs.filter(p => TRANS.has(p.split('/').pop())).length
  check(n === 2, `MICAL3: ${n} trans object requests (sQTL only: the header and one range over all its introns' frames)`)
  const s = await csv(page)
  compare(s.rows, new Set(s.rows.map(r => r.phenotype_id)), 'MICAL3 trans sQTL CSV')
  await page.getByText(/^eQTL$/).first().click()
  await page.waitForSelector('[data-trans-total="0"]', { timeout: 60_000 })
  check(await page.getByText('No trans associations.').count() === 1, 'MICAL3: eQTL tab shows an empty trans table')
  check(errors.length === 0, `MICAL3: no console errors${errors.length ? ` (${errors.slice(0, 2).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- variant page: rs4819361's 10 trans associations from its hits frame ----
{
  const { ctx, page, errors, reqs } = await open('/variant/rs4819361')
  await page.waitForSelector('[data-trans-total="10"]', { timeout: 60_000 })
  check(!reqs.some(p => TRANS.has(p.split('/').pop())), 'rs4819361: no trans object request (the hits frame holds the rows)')
  const v = await csv(page)
  // every row names a phenotype whose trans frame holds this variant with the same p and beta
  if (REF) {
    let bad = 0
    for (const r of v.rows) {
      const t = REF.trans.find(x => x.phenotype_id === r.phenotype_id)
      // the variant-keyed CSV lists genes, not the variant: rs4819361 is chr21:44,000,956
      const i = t ? t.pos.findIndex((p, k) => p === 44000956 && t.chr[k] === 'chr21') : -1
      // hits carry the source's -log10 p and beta as f32; the frame holds them quantized, so they
      // agree within the frame's rounding bounds (SPEC section 13: nlp_max / 131066, beta_max / 65534)
      if (i < 0 || Math.abs(-Math.log10(Number(r.pval)) - t.nlp[i]) > t.nlp_max / 131066 + 1e-6 ||
          Math.abs(Number(r.beta) - t.beta[i]) > t.beta_max / 65534 + 1e-6 * Math.abs(t.beta[i])) bad++
    }
    check(bad === 0 && v.rows.length === 10, `rs4819361: its ${v.rows.length} trans rows match the phenotypes' trans frames (${bad} off)`)
  }
  check(errors.length === 0, `rs4819361: no console errors${errors.length ? ` (${errors.slice(0, 2).join(' | ')})` : ''}`)
  await ctx.close()
}

await browser.close()
const failed = results.filter(r => !r[0])
console.log(failed.length ? `smoke-trans: ${failed.length} of ${results.length} failed` : `smoke-trans: all ${results.length} checks passed`)
process.exit(failed.length ? 1 : 0)
