// Plan 8 step 8 smoke check against the local preview: the gene page's trans tables on the trans pack.
// Click-through for FLNC, HHATL (the largest frame), and a chrM gene, with the rows the page shows
// compared against the trans parquet rows in trans_smoke_rows.json (written next to this script).
//   node data/derived/_tmp/pack_check/smoke_trans_tab.mjs
import { chromium } from '/home/nsheff/Dropbox/workspaces/assistant/tasks/qtl-browser/ui/node_modules/playwright/index.mjs'
import { readFileSync } from 'node:fs'

const BASE = 'http://localhost:4173'
const HERE = '/home/nsheff/Dropbox/workspaces/assistant/tasks/qtl-browser/data/derived/_tmp/pack_check'
const REF = JSON.parse(readFileSync(`${HERE}/trans_smoke_rows.json`, 'utf8')).genes
const E_COLS = ['variant_chr', 'position', 'rsid', 'af', 'pval', 'beta', 'beta_se', 'r2']
const S_COLS = ['phenotype_id', ...E_COLS]
const AF_TOL = 7.63e-6, P_REL = 5e-3, SE_REL = 0.01, R2_ABS = 0.002
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }

const browser = await chromium.launch({ headless: true })

async function open(path) {
  const ctx = await browser.newContext({ viewport: { width: 1280, height: 900 }, acceptDownloads: true })
  const page = await ctx.newPage()
  const errors = [], reqs = []
  let inflight = 0, last = Date.now()
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()) })
  page.on('pageerror', e => errors.push(String(e)))
  const net = new Set()
  ctx.on('request', r => {
    const u = new URL(r.url())
    if (!u.protocol.startsWith('http')) return
    net.add(r); inflight++; last = Date.now()
    if (u.pathname.startsWith('/data/')) reqs.push({ path: u.pathname, range: r.headers()['range'] ?? null })
  })
  const done = r => { if (net.delete(r)) { inflight--; last = Date.now() } }
  ctx.on('requestfinished', done)
  ctx.on('requestfailed', done)
  await page.goto(BASE + path)
  const idle = async (quiet = 1500) => {
    const t0 = Date.now()
    while (Date.now() - t0 < 30_000) { if (inflight === 0 && Date.now() - last >= quiet) return true; await page.waitForTimeout(100) }
    return false
  }
  return { ctx, page, errors, reqs, idle }
}

// packs live flat under immutable/ as <stem>.<chr>.<sha16>.<ext> (SPEC section 3)
const TRANS_PACK = /\/immutable\/trans\./
const kinds = reqs => ({
  pack: reqs.filter(r => TRANS_PACK.test(r.path)),
  parquet: reqs.filter(r => r.path.includes('/trans_pairs/')),
})
const transTable = page => page.locator('[data-trans-total]')
const transTotal = async page => Number(await transTable(page).getAttribute('data-trans-total'))
const waitTotal = (page, n) => page.waitForSelector(`[data-trans-total="${n}"]`, { timeout: 60_000 })

async function csv(page) {
  const [dl] = await Promise.all([
    page.waitForEvent('download', { timeout: 60_000 }),
    transTable(page).locator('button:has-text("CSV")').first().click(),
  ])
  const text = readFileSync(await dl.path(), 'utf8').trim().split('\n')
  const header = text[0].split(',')
  return { header, rows: text.slice(1).map(line => Object.fromEntries(line.split(',').map((v, i) => [header[i], v]))) }
}

/** Every exported row against the parquet reference row for the same variant. */
function compare(got, ref, cols, betaMax, label) {
  // sQTL rows are unique only by intron and variant: a variant can hit more than one of the gene's introns
  const key = r => `${cols.includes('phenotype_id') ? `${r.phenotype_id}|` : ''}${r.variant_chr}:${r.position}`
  const want = new Map(ref.map(r => [key(r), r]))
  const bad = new Map()
  const bump = k => bad.set(k, (bad.get(k) ?? 0) + 1)
  const worst = { af: 0, p: 0, beta: 0, se: 0, r2: 0 }
  for (const row of got) {
    const r = want.get(key(row))
    if (!r) { bump('row not in the reference'); continue }
    if ((row.rsid || null) !== (r.rsid ?? null)) bump('rsid')
    if (cols.includes('phenotype_id') && row.phenotype_id !== r.phenotype_id) bump('phenotype_id')
    const e = (k, v) => Math.abs(Number(row[k]) - v)
    worst.af = Math.max(worst.af, e('af', r.af)); if (!(e('af', r.af) <= AF_TOL)) bump('af')
    const pRel = e('pval', r.pval) / r.pval
    worst.p = Math.max(worst.p, pRel); if (!(pRel <= P_REL)) bump('pval')
    const bLim = betaMax / 65534 + 1.2e-7 * Math.abs(r.beta) + 1.2e-7 * Math.abs(r.beta)
    worst.beta = Math.max(worst.beta, e('beta', r.beta)); if (!(e('beta', r.beta) <= bLim)) bump('beta')
    const seRel = e('beta_se', r.beta_se) / r.beta_se
    worst.se = Math.max(worst.se, seRel); if (!(seRel <= SE_REL)) bump('beta_se')
    worst.r2 = Math.max(worst.r2, e('r2', r.r2)); if (!(e('r2', r.r2) <= R2_ABS)) bump('r2')
  }
  const g = x => x.toPrecision(3)
  check(bad.size === 0 && got.length === ref.length,
    `${label}: ${got.length} exported rows (reference ${ref.length}) match the trans parquet: variant, rsID${cols.includes('phenotype_id') ? ', intron id' : ''} exact; ` +
    `worst af ${g(worst.af)}, p ${g(worst.p)} relative, beta ${g(worst.beta)}, beta_se ${g(worst.se)} relative, r2 ${g(worst.r2)}` +
    `${bad.size ? ` (${[...bad].map(([k, v]) => `${k}: ${v}`).join(', ')})` : ''}`)
}

// ---- FLNC: both tabs off one frame ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/FLNC')
  const ref = REF.FLNC
  await waitTotal(page, ref.n_e)
  await idle()
  const k = kinds(reqs)
  check(k.pack.length === 1 && k.parquet.length === 0 && k.pack[0].range === `bytes=6987563-${6987563 + ref.trans_len - 1}`,
    `FLNC: ${k.pack.length} trans pack request (${k.pack[0]?.range}, ${ref.trans_len} bytes), ${k.parquet.length} trans parquet requests`)
  compare((await csv(page)).rows, ref.e, E_COLS, ref.beta_max, 'FLNC trans eQTL CSV')
  const first = await transTable(page).locator('tbody tr td').first().innerText()
  await transTable(page).locator('th button:text-is("p")').click()
  await page.waitForTimeout(400)
  const flipped = await transTable(page).locator('tbody tr td').first().innerText()
  check(first !== flipped, `FLNC: sorting the trans table by p flips the first row (${first} -> ${flipped})`)
  await transTable(page).locator('th button:text-is("p")').click()
  const rsid = ref.e.find(r => r.rsid)?.rsid
  await transTable(page).locator('input[placeholder="rsID or position"]').fill(rsid)
  await page.waitForTimeout(600)
  const found = await transTable(page).locator('tbody tr').count()
  check(found === 1, `FLNC: searching the trans table for ${rsid} finds ${found} row`)
  await transTable(page).locator('input[placeholder="rsID or position"]').fill('')
  const before = reqs.length
  await page.getByText(/^sQTL/).first().click()
  await waitTotal(page, ref.n_s)
  await idle()
  check(reqs.length === before + 1 && kinds(reqs).pack.length === 1,
    `FLNC: the sQTL tab shows its ${ref.n_s} trans sQTL rows off the same frame (${reqs.length - before} new request, the intron block)`)
  const sqtl = await csv(page)
  check(sqtl.header.join(',') === S_COLS.join(','), `FLNC: trans sQTL CSV header ${sqtl.header.join(',')}`)
  compare(sqtl.rows, ref.s, S_COLS, ref.beta_max, 'FLNC trans sQTL CSV')
  const ids = [...new Set(sqtl.rows.map(r => r.phenotype_id))].sort()
  check(ids.length === ref.s_phenotypes.length && ids.every((x, i) => x === ref.s_phenotypes[i]),
    `FLNC: the ${ids.length} rebuilt intron ids equal the parquet's (${ids[0]})`)
  check(errors.length === 0, `FLNC: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- HHATL: the largest frame, on the sQTL tab ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/HHATL?tab=sqtl')
  const ref = REF.HHATL
  await waitTotal(page, ref.n_s)
  await idle()
  const k = kinds(reqs)
  check(k.pack.length === 1 && k.parquet.length === 0, `HHATL: ${k.pack.length} trans pack request (${ref.trans_len} bytes), ${k.parquet.length} trans parquet requests`)
  const marks = await page.evaluate(() => Object.fromEntries(['pack:trans', 'pack:trans-decode', 'pack:worker-wait', 'pack:insert']
    .map(n => [n, performance.getEntriesByType('measure').filter(m => m.name === n).map(m => Math.round(m.duration))])))
  check(marks['pack:trans'].length === 1 && marks['pack:trans-decode'].length === 1,
    `HHATL: measures pack:trans ${marks['pack:trans']} ms, pack:trans-decode ${marks['pack:trans-decode']} ms, worker wait ${marks['pack:worker-wait']} ms, insert ${marks['pack:insert']} ms`)
  const before = reqs.length
  await page.getByText(/^eQTL$/).first().click()
  await waitTotal(page, ref.n_e)
  check(reqs.length === before && kinds(reqs).pack.length === 1, `HHATL: the eQTL tab shows its ${ref.n_e} trans eQTL rows off the same frame (${reqs.length - before} new requests)`)
  check(errors.length === 0, `HHATL: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- MT-ATP6 (chrM): a gene with trans rows but no gene block stays hidden (T4) ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/MT-ATP6')
  await page.getByText(/was not tested for QTL/).first().waitFor({ timeout: 30_000 })
  await idle()
  const k = kinds(reqs)
  // the app-level startup files (the GWAS index, the variant index and the search index) load on
  // every page; only per-gene pack reads matter here
  const packs = reqs.filter(r => r.path.includes('/immutable/') && !/\/(gwas_index|variant_index|rsid_index|search_index)\./.test(r.path))
  check(k.pack.length === 0 && k.parquet.length === 0 && packs.length === 0 && (await transTable(page).count()) === 0,
    `MT-ATP6 (chrM, trans rows but no gene block): hidden, ${packs.length} per-gene pack requests`)
  check(errors.length === 0, `MT-ATP6: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

await browser.close()
const failed = results.filter(r => !r[0])
console.log(failed.length ? `smoke-trans: ${failed.length} of ${results.length} failed` : `smoke-trans: all ${results.length} checks passed`)
process.exit(failed.length ? 1 : 0)
