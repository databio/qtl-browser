// Smoke check of the gene page on the qtlb v1 store, against the local preview serving the
// chr21/chr22 two-experiment smoke store (bench/README.md, "Local v1 store"). Adapted from the v0
// suite: PDXK's panels and cis table, request counts per tab and intron, the "show all tested
// introns" toggle, an sQTL-only gene, an untested gene, the sQTL CSV, the trans tables and the
// QTL-versus-GWAS panel, the Home colocalization track, the About counts, and dead code.
//   node bench/smoke_gene_page.mjs
import { chromium } from '../node_modules/playwright/index.mjs'
import { readFileSync } from 'node:fs'
import { spawnSync } from 'node:child_process'
import { fileURLToPath } from 'node:url'

const BASE = process.env.SMOKE_BASE ?? 'http://localhost:4173'
const UI = fileURLToPath(new URL('..', import.meta.url))
const READY = '.plot-host:not(.invisible) svg g[aria-label="dot"] > *'
const ROUNDING = /are rounded \(p within [\d.]+%\)\. Exact values: Zenodo\./
const CIS_COLS = ['position', 'rsid', 'A1', 'A2', 'tss_distance', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se', 'pip', 'cs_id']
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }
const badPaths = []
let dataTotal = 0

// the object each per-gene request reads, by name from the store's pointers
const store = await (await fetch(`${BASE}/data/store.json`)).json()
const exp = await (await fetch(`${BASE}/data/experiments/topchef.json`)).json()
const res = Object.fromEntries(exp.results.map(r => [r.phenotype_type, r.files]))
const cat = await (await fetch(`${BASE}/data/variant_catalogs/${exp.catalog}.json`)).json()
const VARIANTS = new Set(cat.chromosomes.map(c => c.file))
const TRANS = new Set(exp.results.map(r => r.trans?.file).filter(Boolean))
const GWAS = new Set(Object.values(exp.gwas?.files ?? {}))
const kindOf = p => {
  const name = p.split('/').pop()
  if (Object.values(res.ge).includes(name)) return 'ge'
  if (Object.values(res.leafcutter).includes(name)) return 'leafcutter'
  if (VARIANTS.has(name)) return 'variants'
  if (TRANS.has(name)) return 'trans'
  if (GWAS.has(name)) return 'gwas'
  return 'other'
}
check(store.format_version === 1, `store.json is format version 1 (${store.name})`)

const browser = await chromium.launch({ headless: true })

async function open(path, { ready = true } = {}) {
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
    if (u.pathname.startsWith('/data/')) {
      reqs.push(u.pathname); dataTotal++
      // only pointers and content-addressed objects are ever read
      if (!/^\/data\/(store\.json|(experiments|variant_catalogs|annotations)\/[\w.-]+\.json|immutable\/[A-Za-z0-9_-]{32}\.[a-z0-9.]+)$/.test(u.pathname)) badPaths.push(`${path}: ${u.pathname}`)
    }
  })
  const done = r => { if (net.delete(r)) { inflight--; last = Date.now() } }
  ctx.on('requestfinished', done)
  ctx.on('requestfailed', done)
  await page.goto(BASE + path)
  if (ready) await page.waitForSelector(READY, { timeout: 60_000 })
  const idle = async (quiet = 2000) => {
    const t0 = Date.now()
    while (Date.now() - t0 < 30_000) { if (inflight === 0 && Date.now() - last >= quiet) return true; await page.waitForTimeout(100) }
    return false
  }
  return { ctx, page, errors, reqs, idle }
}
const count = (reqs, kind, from = 0) => reqs.slice(from).filter(p => kindOf(p) === kind).length

const cisTable = page => page.locator('table:has(th:has-text("TSS dist"))')
const pagerText = page => page.locator('span.tabular-nums:has-text(" of ")').first().innerText()
async function cisRows(page) {
  await page.waitForTimeout(400)
  return await cisTable(page).locator('tbody tr').evaluateAll(trs => trs.map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim())))
}
async function waitPager(page, pred) {
  for (let i = 0; i < 50; i++) { const t = await pagerText(page); if (pred(t)) return t; await page.waitForTimeout(100) }
  return await pagerText(page)
}
const intronTable = page => page.locator('table:has(th:has-text("Perm p"))')
const intronRows = page => intronTable(page).locator('tbody tr').evaluateAll(trs => trs.map(tr => ({
  text: [...tr.querySelectorAll('td')].map(td => td.textContent.trim()), selected: tr.classList.contains('bg-base-200') })))
const showAll = page => page.locator('label:has-text("Show all") input[type=checkbox]').first()
const drawnCount = page => page.evaluate(() => performance.getEntriesByName('locus:drawn', 'mark').length)
async function waitDrawn(page, before) {
  await page.waitForFunction(n => performance.getEntriesByName('locus:drawn', 'mark').length > n, before, { timeout: 60_000 })
  await page.waitForSelector(READY, { timeout: 60_000 })
}
const locusHeader = page => page.evaluate(() => [...document.querySelectorAll('span')].map(s => s.textContent.trim())
  .find(t => /^intron .* · [\d,]+ variants$/.test(t)) ?? null)
const nVariants = h => (h ? Number(/· ([\d,]+) variants/.exec(h)[1].replace(/,/g, '')) : null)
async function downloadCsv(page) {
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 30_000 }), page.locator('button:has-text("CSV")').first().click()])
  return readFileSync(await dl.path(), 'utf8').trim().split('\n')
}

// ---- PDXK (chr21): eQTL tab ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/PDXK')
  const lead = await page.locator('tr:has-text("Lead variant") td').nth(1).innerText()
  check(lead.trim() === 'rs117208074', `PDXK: lead variant from the block details and variants range: ${lead.trim()}`)
  await page.waitForFunction(() => document.querySelectorAll('.plot-host svg').length >= 2, null, { timeout: 30_000 })
  check((await page.locator('.plot-host svg').count()) >= 2, 'PDXK: locus plot and the QTL-versus-GWAS panel both drawn')
  const csRows = await page.locator('table:has(th:has-text("Top PIP")) tbody tr').count()
  check(csRows === 3, `PDXK: credible-set table has ${csRows} set rows (block details say 3)`)
  await cisTable(page).locator('tbody tr').first().waitFor({ timeout: 30_000 })
  const total = await pagerText(page)
  check(/of 5,827/.test(total), `PDXK: cis table pager reads "${total}" (search index n_var 5,827)`)
  const asc = await cisRows(page)
  await cisTable(page).locator('th button:text-is("p")').click()
  const desc = await cisRows(page)
  check(asc[0][0] !== desc[0][0], `PDXK: sorting by p flips the first row (${asc[0][1]} -> ${desc[0][1]})`)
  await cisTable(page).locator('th button:text-is("p")').click()
  await page.selectOption('select[title="Nominal p-value threshold"]', '1e-5')
  const filtered = await waitPager(page, t => !/of 5,827/.test(t))
  check(!/of 5,827/.test(filtered), `PDXK: p <= 1e-5 filter narrows the table: "${filtered}"`)
  await page.selectOption('select[title="Nominal p-value threshold"]', '')
  await waitPager(page, t => /of 5,827/.test(t))
  await page.fill('input[placeholder="rsID or position"]', 'rs117208074')
  const searched = await waitPager(page, t => /of 1$/.test(t))
  check(/of 1$/.test(searched), `PDXK: rsID search rs117208074 finds one row: "${searched}"`)
  await page.fill('input[placeholder="rsID or position"]', '')
  await waitPager(page, t => /of 5,827/.test(t))
  await page.locator('button:text-is("Next")').first().click()
  const paged = await waitPager(page, t => t.startsWith('11–'))
  check(paged.startsWith('11–20'), `PDXK: paging Next shows "${paged}"`)
  const csv = await downloadCsv(page)
  check(csv.length === 5828 && csv[0] === CIS_COLS.join(','), `PDXK: CSV export has ${csv.length - 1} rows, header ${csv[0]}`)
  check(await page.getByText(ROUNDING).first().isVisible(), 'PDXK: rounding note by the CSV button')
  // PDXK's eQTL phenotype has no trans rows; its introns have 2 (search index n_trans)
  check(await page.locator('[data-trans-total="0"]').count() === 1 && await page.getByText('No trans associations.').count() === 1,
    'PDXK: the eQTL trans table is present and empty')
  const toggle = page.locator('label:has-text("Gene track") input[type=checkbox]')
  const before = await page.locator('svg').count()
  await toggle.click()
  await page.waitForTimeout(1500)
  const after = await page.locator('svg').count()
  check(await toggle.isChecked() && after > before, `PDXK: gene track toggle adds the track (${before} -> ${after} svg elements)`)
  await toggle.click()
  await idle()
  // one header read plus one range for each ranged object
  const k = Object.fromEntries(['ge', 'variants', 'leafcutter', 'gwas', 'trans'].map(x => [x, count(reqs, x)]))
  check(k.ge === 2 && k.variants === 2 && k.leafcutter === 0 && k.gwas === 2 && k.trans === 2,
    `PDXK: requests eQTL results ${k.ge} (header + block), variants ${k.variants} (header + range), GWAS ${k.gwas} (header + window), trans ${k.trans} (header + frames), sQTL results ${k.leafcutter}`)
  const n0 = reqs.length
  const d = await drawnCount(page)
  await page.getByText(/^sQTL/).first().click()
  await waitDrawn(page, d); await idle()
  const sqtlText = await page.locator('text=/tested introns/').first().innerText()
  check(/tested introns/.test(sqtlText) && count(reqs, 'leafcutter', n0) === 2 && reqs.length - n0 === 2,
    `PDXK: sQTL tab adds ${reqs.length - n0} requests (sQTL header + introns span): "${sqtlText.slice(0, 70)}"`)
  check(await page.locator('[data-trans-total="2"]').count() === 1, 'PDXK: the sQTL trans table has its 2 rows from the frames read when the gene opened')
  check(errors.length === 0, `PDXK: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- MICAL3 (chr22): requests per intron choice and tab switch, show all, sQTL CSV ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/MICAL3?tab=sqtl')
  await idle()
  const label = (await page.locator('label:has-text("Show all")').first().innerText()).trim()
  const N = Number(/Show all ([\d,]+) tested/.exec(label)[1].replace(/,/g, ''))
  const onAtStart = await showAll(page).isChecked()
  const sig = await intronRows(page)
  check(!onAtStart && sig.length === 9 && N === 33 && sig[0].selected, `MICAL3: toggle off at start, ${sig.length} significant of ${N} introns, first selected`)
  await showAll(page).check()
  await page.waitForTimeout(300)
  const all = await intronRows(page)
  const key = r => `${r.text[0]} ${r.text[1]}`
  const sigKeys = new Set(sig.map(key))
  const idx = all.findIndex(r => !sigKeys.has(key(r)))
  const n1 = reqs.length
  let d = await drawnCount(page)
  await intronTable(page).locator('tbody tr').nth(idx).locator('td').first().click()
  await waitDrawn(page, d); await idle()
  const header = await locusHeader(page)
  check(all.length === N && (await intronRows(page))[idx].selected && nVariants(header) > 0 && reqs.length === n1,
    `MICAL3: non-significant intron ${key(all[idx])} (perm p ${all[idx].text[4]}) draws "${header}" with ${reqs.length - n1} requests (its block came with the span)`)
  await cisTable(page).locator('tbody tr').first().waitFor({ timeout: 30_000 })
  const csv = await downloadCsv(page)
  check(csv[0] === ['phenotype_id', ...CIS_COLS].join(',') && csv.length - 1 === nVariants(header), `MICAL3: sQTL CSV has ${csv.length - 1} rows, header ${csv[0]}`)
  d = await drawnCount(page)
  const n2 = reqs.length
  await page.getByText(/^eQTL$/).first().click()
  await waitDrawn(page, d); await idle()
  check(reqs.length === n2 && count(reqs, 'ge') === 2, `MICAL3: sQTL -> eQTL tab sends ${reqs.length - n2} requests (the eQTL block left when the gene opened)`)
  d = await drawnCount(page)
  const n3 = reqs.length
  await page.getByText(/^sQTL/).first().click()
  await waitDrawn(page, d); await idle()
  check(reqs.length === n3, `MICAL3: back to the sQTL tab sends ${reqs.length - n3}`)
  check(await page.locator('[data-trans-total="140"]').count() === 1, 'MICAL3: sQTL trans table has 140 rows (search index n_trans over its introns)')
  check(errors.length === 0, `MICAL3: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- SMARCB1: trans rows on both tabs, and the GWAS panel with shared variants ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/SMARCB1')
  await page.waitForSelector('[data-trans-total="7"]', { timeout: 60_000 })
  await idle()
  const dots = await page.locator('.plot-host svg').nth(1).locator('g[aria-label="dot"] > *').count()
  check(dots > 0, `SMARCB1: the QTL-versus-GWAS panel draws ${dots} shared variants`)
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 30_000 }), page.locator('[data-trans-total] button:has-text("CSV")').first().click()])
  const csv = readFileSync(await dl.path(), 'utf8').trim().split('\n')
  check(csv.length === 8, `SMARCB1: trans eQTL CSV has ${csv.length - 1} rows, header ${csv[0]}`)
  const n0 = reqs.length
  const d = await drawnCount(page)
  await page.getByText(/^sQTL/).first().click()
  await waitDrawn(page, d); await idle()
  check(await page.locator('[data-trans-total="34"]').count() === 1 && count(reqs, 'trans', n0) === 0,
    `SMARCB1: sQTL trans table has 34 rows with no new trans request`)
  check(errors.length === 0, `SMARCB1: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- GUSBP11: sQTL only, opens on its first significant intron ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/GUSBP11')
  await idle()
  const rows = await intronRows(page)
  const sel = rows.findIndex(r => r.selected)
  const header = await locusHeader(page)
  const chip = await page.getByText('no eQTL test').count()
  check(sel === 0 && nVariants(header) > 0 && chip === 1, `GUSBP11: opens on its sQTL tab with the first significant intron selected (${rows.length} listed), "no eQTL test" chip: "${header}"`)
  check(count(reqs, 'ge') === 0 && count(reqs, 'leafcutter') === 2 && count(reqs, 'variants') === 2,
    `GUSBP11: requests eQTL ${count(reqs, 'ge')}, sQTL ${count(reqs, 'leafcutter')}, variants ${count(reqs, 'variants')}`)
  check(errors.length === 0, `GUSBP11: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- FP565260.4: annotated on chr21, no phenotype ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/FP565260.4', { ready: false })
  await page.getByText(/was not tested for QTL/).first().waitFor({ timeout: 30_000 })
  await idle()
  const per = reqs.filter(p => ['ge', 'leafcutter', 'variants', 'gwas', 'trans'].includes(kindOf(p))).length
  check(per === 0, `FP565260.4: "not tested" with ${per} per-gene requests`)
  check(errors.length === 0, `FP565260.4: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- Home, genes list, region, About ----
{
  const { ctx, page, errors, idle } = await open('/', { ready: false })
  const sentence = page.locator('[data-sentence]', { hasText: /eGenes and [\d,]+ significant sQTL introns/ }).first()
  await sentence.waitFor({ state: 'attached', timeout: 30_000 })
  await idle()
  const line = await sentence.textContent()
  check(/^372 eGenes and 512 significant sQTL introns across 686 tested genes/.test(line), `Home: counts from the search index: "${line.slice(0, 70)}"`)
  const loci = await page.getByText(/PP\.H4 > 0\.8, 4 loci/).count()
  check(loci === 1, 'Home: the colocalization track places the 4 coloc genes on chr21/chr22 (VPREB3, MAP3K7CL, MMP11, SMARCB1)')
  await page.goto(`${BASE}/genes`)
  await page.locator('table tbody tr').first().waitFor({ timeout: 30_000 })
  const meta = await page.locator('h1 + *').first().innerText().catch(() => '')
  check((await page.locator('table tbody tr').count()) > 0, `Genes: the eGene list renders (${meta.trim()})`)
  await page.goto(`${BASE}/region/chr21:43000000-44000000`)
  const pdxk = await page.getByText('PDXK', { exact: true }).first().waitFor({ timeout: 30_000 }).then(() => true, () => false)
  check(pdxk, 'Region chr21:43-44 Mb lists PDXK')
  await page.goto(`${BASE}/about`)
  await page.getByText('Genes tested').first().waitFor({ timeout: 30_000 })
  const text = await page.locator('body').innerText()
  check(/Genes tested\s+686/.test(text) && /DCM GWAS variants\s+349,950/.test(text) && /trans pairs\s+15,048/.test(text) && /5,022 cases/.test(text),
    'About: genes tested 686, DCM GWAS variants 349,950, trans pairs 15,048, 5,022 cases')
  check(errors.length === 0, `Home/Genes/Region/About: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

await browser.close()

check(badPaths.length === 0, `every /data/ request is a pointer or an immutable object (${dataTotal} seen)${badPaths.length ? `: ${badPaths.slice(0, 3).join(', ')}` : ''}`)
const g = spawnSync('grep', ['-rn', 'manifest\\|pack-decode\\|trans-pack\\|variant-pack\\|lib/pack\'', 'src'], { cwd: UI, encoding: 'utf8' })
check(g.status === 1 && !g.stdout.trim(), `dead v0 code grep finds nothing in ui/src${g.stdout.trim() ? `: ${g.stdout.trim().split('\n').slice(0, 3).join(' | ')}` : ''}`)

const failed = results.filter(r => !r[0])
console.log(failed.length ? `smoke: ${failed.length} of ${results.length} failed` : `smoke: all ${results.length} checks passed`)
process.exit(failed.length ? 1 : 0)
