// Smoke check of the variant page on a version 2 store, against the local preview serving the
// genome-wide local store (bench/README.md, "Local store"). Checks the lead and credible-set lists,
// the trans table (from the paged hits files), the cis scan, and the request counts per lookup path.
// Each page opens inside one experiment (src/lib/scope.ts), so a variant's eQTL, sQTL and trans
// results are on four pages; where a number spans them the checks add the pages up. The numbers
// are the genome-wide TOPCHeF store's.
//   node bench/smoke_variant_page.mjs
import { chromium } from '../node_modules/playwright/index.mjs'

const BASE = process.env.SMOKE_BASE ?? 'http://localhost:4173'
const LISTS = '[data-variant-lists="1"]'
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }

const store = await (await fetch(`${BASE}/data/store.json`)).json()
// the study the app opens (its default VITE_COLLECTION_OF and VITE_GWAS), not every cohort in the store
const studyIds = Object.keys(store.experiments).filter(id =>
  store.experiments[id].collection === store.experiments.topchef_ge_cis.collection || id === 'dcm_jurgens2024_biobanks')
const at = id => `/c/${store.experiments[id].collection.replace(/\.json$/, '')}/e/${id}`
const E = at('topchef_ge_cis'), S = at('topchef_leafcutter_cis'), TE = at('topchef_ge_trans'), TS = at('topchef_leafcutter_trans')
const docs = Object.values(await Promise.all(studyIds.map(async id => (await fetch(`${BASE}/data/experiments/${id}.json`)).json())))
const cisOf = type => docs.find(d => d.kind === 'qtl_cis' && d.phenotype_type === type)
const cat = await (await fetch(`${BASE}/data/variant_catalogs/${cisOf('ge').catalog}.json`)).json()
const qtl = docs.filter(d => d.kind === 'qtl_cis' || d.kind === 'qtl_trans')
const names = new Map([
  ...Object.values(cisOf('ge').files).map(f => [f, 'ge']),
  ...Object.values(cisOf('leafcutter').files).map(f => [f, 'leafcutter']),
  ...cat.chromosomes.map(c => [c.file, 'variants']), ...qtl.flatMap(d => Object.values(d.hits).map(f => [f, 'hits'])), [cat.rsid, 'rsid'],
])
/** Per-variant requests by kind; boot's whole-object reads (pointers, search index, genes, variant index) are left out. */
const kinds = reqs => reqs.map(p => names.get(p.split('/').pop())).filter(Boolean).sort().join(',')
// a hits file is read as its frame table then the variant's frame, for the page's one experiment;
// an empty frame sends no request
const HITS = n => Array(n).fill('hits')
check(qtl.length === 4, `${qtl.length} QTL experiments (hits files per chromosome: ${qtl.length})`)

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

// ---- rs34599497 (chr22): lead of 21 phenotypes (1 eQTL, 20 introns), in 50 credible sets, 1 trans association ----
{
  // the eQTL experiment: its one lead, its credible sets, its scan
  const { ctx, page, errors, reqs, idle } = await open(`${E}/variant/rs34599497`)
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const leadsE = await rowsOf(table(page, 'Perm p'))
  const csE = await rowsOf(table(page, 'PIP'))
  const ld = leadsE.filter(r => r[0] === 'eQTL' && r[1] === 'AP000346.2')[0]
  check(leadsE.length === 1 && ld && /^-0\.412 ± 0\.032$/.test(ld[3]) && ld[4] === '1.0e-4' && ld[5] === 'eGene',
    `rs34599497 (eQTL): ${leadsE.length} lead row; AP000346.2 lead beta ${ld?.[3]}, perm p ${ld?.[4]}`)
  // rsID block and page (no header reads), the experiment's hits table and frame, its lead block
  const k = kinds(reqs)
  check(k === ['ge', ...HITS(2), 'rsid', 'variants'].sort().join(','), `rs34599497 (eQTL): rsID block, variants page, hits table + frame, lead block (got ${k})`)
  check((await page.locator('h1').first().innerText()).trim() === 'rs34599497', 'rs34599497: header title')
  check(await page.getByText('chr22:23,680,950 (GRCh38)').count() === 1 && await page.getByText('T / C', { exact: true }).count() === 1,
    'rs34599497: position chr22:23,680,950, A1 / A2 = T / C (ALT / REF)')
  check(await page.locator('[data-trans-total]').count() === 0, 'rs34599497 (eQTL): no trans section in a cis experiment')
  let before = reqs.length
  await page.locator('[data-scan-button]').click()
  await page.waitForSelector('[data-scan-ready="1"]', { timeout: 120_000 })
  await idle()
  check(kinds(reqs.slice(before)) === 'ge', `rs34599497 (eQTL) scan: one eQTL span (got ${kinds(reqs.slice(before))})`)
  const eTitle = (await page.getByRole('heading', { name: /^Expression \(/ }).innerText()).trim()
  // search index: 46 ge and 267 leafcutter phenotypes cover chr22 vidx 22469, 96 of the introns significant
  check(eTitle === 'Expression (46)', `rs34599497 (eQTL) scan: ${eTitle} (search index: 46)`)

  // the sQTL experiment: its 20 intron leads, in 7 runs of blocks < 64 KB apart (v2 block offsets:
  // four of the gaps are 75-77 KB, just over the limit)
  const n0 = reqs.length
  await page.goto(`${BASE}${S}/variant/rs34599497`)
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const leadsS = await rowsOf(table(page, 'Perm p'))
  const csS = await rowsOf(table(page, 'PIP'))
  const kS = kinds(reqs.slice(n0))
  check(leadsS.length === 20 && csE.length + csS.length === 50,
    `rs34599497 (sQTL): ${leadsS.length} lead rows; credible-set rows ${csE.length} + ${csS.length} over both cis experiments`)
  check(kS === [...Array(7).fill('leafcutter'), ...HITS(2), 'rsid', 'variants'].sort().join(','), `rs34599497 (sQTL): lead blocks in 7 runs (got ${kS})`)
  before = reqs.length
  await page.locator('[data-scan-button]').click()
  await page.waitForSelector('[data-scan-ready="1"]', { timeout: 120_000 })
  await idle()
  check(kinds(reqs.slice(before)) === 'leafcutter', `rs34599497 (sQTL) scan: one sQTL span (got ${kinds(reqs.slice(before))})`)
  const sTitle = (await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()).trim()
  const toggle = (await page.locator('button:has-text("Show all")').innerText()).trim()
  check(sTitle === 'Splicing (96)' && toggle === 'Show all 267 tested introns', `rs34599497 (sQTL) scan: ${sTitle}, "${toggle}" (96 of 267)`)
  await page.locator('button:has-text("Show all")').click()
  check((await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()).trim() === 'Splicing (267)', 'rs34599497 scan: the toggle shows all 267')

  // its one trans record, in the eQTL trans scan
  await page.goto(`${BASE}${TE}/variant/rs34599497`)
  const t1 = await page.locator('[data-trans-total="1"]').waitFor({ timeout: 60_000 }).then(() => true, () => false)
  check(t1, 'rs34599497 (eQTL trans): the trans table has its 1 row')
  check(errors.length === 0, `rs34599497: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- the same variant by chr:pos: no rsID request ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${E}/variant/chr22:23680950`)
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const k = kinds(reqs)
  check(!k.includes('rsid') && k.split(',').filter(x => x === 'hits').length === 2 && k.split(',').filter(x => x === 'variants').length === 1, `chr22:23680950: no rsID request; variants page, hits table + frame (got ${k})`)
  check((await page.locator('h1').first().innerText()).trim() === 'rs34599497', 'chr22:23680950: resolves to rs34599497')
  check(errors.length === 0, `chr22:23680950: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs4819361 (chr21): a cis variant with 10 trans associations (2 eQTL, 8 sQTL) ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${TE}/variant/rs4819361`)
  await page.waitForSelector(`${LISTS} [data-trans-total="2"]`, { timeout: 60_000 })
  const rowsE = await rowsOf(page.locator('[data-trans-total] table'))
  await page.goto(`${BASE}${TS}/variant/rs4819361`)
  await page.waitForSelector(`${LISTS} [data-trans-total="8"]`, { timeout: 60_000 })
  await idle()
  const rows = [...rowsE, ...await rowsOf(page.locator('[data-trans-total] table'))]
  check(rows.length === 10 && rows.every(r => r.length > 3), `rs4819361: the two trans experiments list ${rows.length} genes and introns from its hits frames (2 + 8)`)
  check(!kinds(reqs).includes('trans'), 'rs4819361: the variant page reads no trans object (the hits frames carry the rows)')
  check(errors.length === 0, `rs4819361: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs457868 (chr21): a trans-only site, outside every cis window, with 1 trans association ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${E}/variant/rs457868`)
  await page.waitForSelector(LISTS, { timeout: 60_000 })
  await idle()
  const k = kinds(reqs)
  check(!k.includes('ge') && !k.includes('leafcutter') && k.includes('rsid') && k.includes('variants'),
    `rs457868 (eQTL): rsID block, variants page, hits only (got ${k})`)
  check((await page.locator('[data-scan-button]').count()) === 0, 'rs457868: no scan button (outside every cis window)')
  const outside = await page.getByText('outside every cis window', { exact: false }).count()
  check(outside === 3, `rs457868: the outside-cis message stands in for all three cis sections (${outside} of 3)`)
  let trans = 0
  for (const where of [TE, TS]) {
    await page.goto(`${BASE}${where}/variant/rs457868`)
    await page.waitForSelector(`${LISTS} [data-trans-total]`, { timeout: 60_000 })
    trans += Number(await page.locator('[data-trans-total]').getAttribute('data-trans-total'))
  }
  check(trans === 1, `rs457868: its 1 trans association across the two trans experiments (${trans})`)
  check(errors.length === 0, `rs457868: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs1: below the first rsID block, answered from the variant index alone ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${E}/variant/rs1`)
  await page.waitForSelector('text=is not among the variants tested', { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === '', `rs1: no per-variant request (got ${kinds(reqs) || 'none'})`)
  check((await page.getByText('Look it up in dbSNP').count()) === 1, 'rs1: offers the dbSNP link')
  check(errors.length === 0, `rs1: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- a malformed id, and a chromosome the store does not hold ----
for (const id of ['not-a-variant', 'chrY:1000000']) {
  const { ctx, page, errors, reqs, idle } = await open(`${E}/variant/${id}`)
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
