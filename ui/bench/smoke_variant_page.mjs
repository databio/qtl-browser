// Smoke check of the variant page on the qtlb v1 store, against the local preview serving the
// chr21/chr22 smoke store (bench/README.md, "Local v1 store"). Checks the lead and credible-set
// lists, the trans table (from the paged hits file), the cis scan, and the request counts per lookup path.
//   node bench/smoke_variant_page.mjs
import { chromium } from '../node_modules/playwright/index.mjs'

const BASE = process.env.SMOKE_BASE ?? 'http://localhost:4173'
const LISTS = '[data-variant-lists="1"]'
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }

const exp = await (await fetch(`${BASE}/data/experiments/topchef.json`)).json()
const cat = await (await fetch(`${BASE}/data/variant_catalogs/${exp.catalog}.json`)).json()
const names = new Map([
  ...Object.values(exp.results.find(r => r.phenotype_type === 'ge').files).map(f => [f, 'ge']),
  ...Object.values(exp.results.find(r => r.phenotype_type === 'leafcutter').files).map(f => [f, 'leafcutter']),
  ...cat.chromosomes.map(c => [c.file, 'variants']), ...Object.values(exp.hits).map(f => [f, 'hits']), [cat.rsid, 'rsid'],
])
/** Per-variant requests by kind; boot's whole-object reads (pointers, search index, genes, variant index) are left out. */
const kinds = reqs => reqs.map(p => names.get(p.split('/').pop())).filter(Boolean).sort().join(',')

const browser = await chromium.launch({ headless: true })

async function open(path) {
  const ctx = await browser.newContext({ viewport: { width: 1280, height: 900 } })
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
    if (u.pathname.startsWith('/data/')) reqs.push(u.pathname)
  })
  const done = r => { if (net.delete(r)) { inflight--; last = Date.now() } }
  ctx.on('requestfinished', done)
  ctx.on('requestfailed', done)
  await page.goto(BASE + path)
  const idle = async (quiet = 1500) => {
    const t0 = Date.now()
    while (Date.now() - t0 < 60_000) { if (inflight === 0 && Date.now() - last >= quiet) return true; await page.waitForTimeout(100) }
    return false
  }
  return { ctx, page, errors, reqs, idle }
}
const table = (page, head) => page.locator(`table:has(th:has-text("${head}"))`)
const rowsOf = loc => loc.locator('tbody tr').evaluateAll(trs => trs.map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim())))

// ---- rs34599497 (chr22): lead of 21 phenotypes, in 50 credible sets ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs34599497')
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const leads = await rowsOf(table(page, 'Perm p'))
  const cs = await rowsOf(table(page, 'PIP'))
  const ld = leads.filter(r => r[0] === 'eQTL' && r[1] === 'AP000346.2')[0]
  check(leads.length === 21 && cs.length === 50 && ld && /^-0\.412 ± 0\.032$/.test(ld[3]) && ld[4] === '1.0e-4' && ld[5] === 'eGene',
    `rs34599497: ${leads.length} lead rows, ${cs.length} credible-set rows; AP000346.2 eQTL lead slope ${ld?.[3]}, perm p ${ld?.[4]}`)
  // rsID block and page (no header reads), the hits file's table and frame, then the lead blocks for
  // their slopes: its 1 eQTL lead block, and its 20 intron lead blocks in 4 runs of blocks < 64 KB apart
  const k = kinds(reqs)
  const want = ['ge', ...Array(4).fill('leafcutter'), 'hits', 'hits', 'rsid', 'variants'].sort().join(',')
  check(k === want, `rs34599497: requests ${k.split(',').length}: rsID block, variants page, hits table + frame, lead blocks in runs (got ${k})`)
  check((await page.locator('h1').first().innerText()).trim() === 'rs34599497', 'rs34599497: header title')
  check(await page.getByText('chr22:23,680,950 (GRCh38)').count() === 1 && await page.getByText('T / C', { exact: true }).count() === 1,
    'rs34599497: position chr22:23,680,950, A1 / A2 = T / C (ALT / REF)')
  check(await page.locator('[data-trans-total="0"]').count() === 1, 'rs34599497: trans table present with no rows (its hits have no kind 2 record)')

  const before = reqs.length
  const label = (await page.locator('[data-scan-button]').innerText()).trim()
  check(/^Scan cis windows \([\d.]+ MB\)$/.test(label), `rs34599497: scan button reads "${label}" (both spans known before any request)`)
  await page.locator('[data-scan-button]').click()
  await page.waitForSelector('[data-scan-ready="1"]', { timeout: 120_000 })
  await idle()
  check(kinds(reqs.slice(before)) === 'ge,leafcutter', `rs34599497 scan: one eQTL span and one sQTL span (got ${kinds(reqs.slice(before))})`)
  const eTitle = (await page.getByRole('heading', { name: /^Expression \(/ }).innerText()).trim()
  const sTitle = (await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()).trim()
  const toggle = (await page.locator('button:has-text("Show all")').innerText()).trim()
  // search index: 46 ge and 267 leafcutter phenotypes cover chr22 vidx 22469, 96 of the introns significant
  check(eTitle === 'Expression (46)' && sTitle === 'Splicing (96)' && toggle === 'Show all 267 tested introns',
    `rs34599497 scan: ${eTitle}, ${sTitle}, "${toggle}" (search index: 46, 96 of 267)`)
  await page.locator('button:has-text("Show all")').click()
  check((await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()).trim() === 'Splicing (267)', 'rs34599497 scan: the toggle shows all 267')
  check(errors.length === 0, `rs34599497: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- the same variant by chr:pos: no rsID request ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/chr22:23680950')
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const k = kinds(reqs)
  check(!k.includes('rsid') && k.split(',').filter(x => x === 'hits').length === 2 && k.split(',').filter(x => x === 'variants').length === 1, `chr22:23680950: no rsID request; variants page, hits table + frame (got ${k})`)
  check((await page.locator('h1').first().innerText()).trim() === 'rs34599497', 'chr22:23680950: resolves to rs34599497')
  check(errors.length === 0, `chr22:23680950: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs4819361 (chr21): a cis variant with 10 trans associations ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs4819361')
  await page.waitForSelector(`${LISTS} [data-trans-total="10"]`, { timeout: 60_000 })
  await idle()
  const rows = await rowsOf(page.locator('[data-trans-total] table'))
  check(rows.length === 10 && rows.every(r => r.length > 3), `rs4819361: trans table lists ${rows.length} genes and introns from its hits frame (kind 2)`)
  check(!kinds(reqs).includes('trans'), 'rs4819361: the variant page reads no trans object (the hits frame carries the rows)')
  check(errors.length === 0, `rs4819361: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs457868 (chr21): a trans-only site, outside every cis window, with 1 trans association ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs457868')
  await page.waitForSelector(`${LISTS} [data-trans-total="1"]`, { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === 'hits,hits,rsid,variants', `rs457868: rsID block, variants page, hits table + frame only (got ${kinds(reqs)})`)
  check((await page.locator('[data-scan-button]').count()) === 0, 'rs457868: no scan button (outside every cis window)')
  const outside = await page.getByText('outside every cis window', { exact: false }).count()
  check(outside === 3, `rs457868: the outside-cis message stands in for all three cis sections (${outside} of 3); trans table has its 1 row`)
  check(errors.length === 0, `rs457868: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs1: below the first rsID block, answered from the variant index alone ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs1')
  await page.waitForSelector('text=is not among the variants tested', { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === '', `rs1: no per-variant request (got ${kinds(reqs) || 'none'})`)
  check((await page.getByText('Look it up in dbSNP').count()) === 1, 'rs1: offers the dbSNP link')
  check(errors.length === 0, `rs1: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- a malformed id, and a chromosome the store does not hold ----
for (const id of ['not-a-variant', 'chr1:1000000']) {
  const { ctx, page, errors, reqs, idle } = await open(`/variant/${id}`)
  await page.waitForSelector('text=is not among the variants tested', { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === '', `${id}: no per-variant request (got ${kinds(reqs) || 'none'})`)
  check(errors.length === 0, `${id}: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

await browser.close()
const failed = results.filter(([ok]) => !ok)
console.log(`\n${results.length - failed.length}/${results.length} checks passed`)
if (failed.length) process.exit(1)
