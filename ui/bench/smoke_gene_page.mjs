// Smoke check of the gene page on a version 2 store, against the local preview serving the
// genome-wide local store (bench/README.md, "Local store"). PDXK's panels and cis table, request
// counts per experiment and intron, the "show all tested introns" toggle, an sQTL-only gene, an
// untested gene, the sQTL CSV, the trans tables and the QTL-versus-GWAS panel, the collection page
// and its experiment pages, the landing page and About, GTEx against the DCM GWAS, and dead code.
// Each page opens inside one experiment (src/lib/scope.ts), so sQTL and trans checks go to those
// experiments. The numbers are the genome-wide TOPCHeF store's (scripts/store_reference.py's
// docstring says how to recompute them).
//   node bench/smoke_gene_page.mjs
import { chromium } from '../node_modules/playwright/index.mjs'
import { readFileSync } from 'node:fs'
import { spawnSync } from 'node:child_process'
import { fileURLToPath } from 'node:url'

const BASE = process.env.SMOKE_BASE ?? 'http://localhost:4173'
const UI = fileURLToPath(new URL('..', import.meta.url))
const READY = '.plot-host:not(.invisible) svg g[aria-label="dot"] > *'
const ROUNDING = /^Values are rounded\.$/   // the bounds are in its tooltip (components/rounding-note.tsx)
const CIS_COLS = ['position', 'rsid', 'A1', 'A2', 'tss_distance', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'beta', 'beta_se', 'pip', 'cs_id']
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }
const badPaths = []
let dataTotal = 0

// the object each per-gene request reads, by name from the store's pointers: one document per scan
const store = await (await fetch(`${BASE}/data/store.json`)).json()
// TOPCHeF's collection and the DCM GWAS, not every cohort in the store
const studyIds = Object.keys(store.experiments).filter(id =>
  store.experiments[id].collection === store.experiments.topchef_ge_cis.collection || id === 'dcm_jurgens2024_biobanks')
// pages open inside one experiment, a cis scan with the DCM GWAS chosen (src/lib/scope.ts)
const coll = id => store.experiments[id].collection.replace(/\.json$/, '')
const at = (id, gwas = 'dcm_jurgens2024_biobanks') => `/c/${coll(id)}/e/${id}${gwas ? `/gwas/${gwas}` : ''}`
const PREFIX = at('topchef_ge_cis')
const SQTL = at('topchef_leafcutter_cis'), TRANS_E = at('topchef_ge_trans', null), TRANS_S = at('topchef_leafcutter_trans', null)
const GTEX = at('gtex_v8_heart_lv_ge_cis', null)
const docs = Object.fromEntries(await Promise.all(studyIds.map(async id => [id, await (await fetch(`${BASE}/data/experiments/${id}.json`)).json()])))
const cisOf = type => Object.values(docs).find(d => d.kind === 'qtl_cis' && d.phenotype_type === type)
const res = { ge: cisOf('ge').files, leafcutter: cisOf('leafcutter').files }
const cat = await (await fetch(`${BASE}/data/variant_catalogs/${cisOf('ge').catalog}.json`)).json()
const VARIANTS = new Set(cat.chromosomes.map(c => c.file))
const TRANS = new Set(Object.values(docs).filter(d => d.kind === 'qtl_trans').map(d => d.file))
const GWAS = new Set(Object.values(docs).filter(d => d.kind === 'gwas').flatMap(d => Object.values(d.files)))
const kindOf = p => {
  const name = p.split('/').pop()
  if (Object.values(res.ge).includes(name)) return 'ge'
  if (Object.values(res.leafcutter).includes(name)) return 'leafcutter'
  if (VARIANTS.has(name)) return 'variants'
  if (TRANS.has(name)) return 'trans'
  if (GWAS.has(name)) return 'gwas'
  return 'other'
}
check(store.format_version === 2 && TRANS.size === 2 && GWAS.size > 0, `store.json is format version 2 (${store.name}): ${Object.keys(docs).length} experiments, ${TRANS.size} trans objects, ${GWAS.size} GWAS objects`)

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
      if (!/^\/data\/(store\.json|(experiments|variant_catalogs|annotations|overlaps|colocs)\/[\w.-]+\.json|immutable\/[A-Za-z0-9_-]{32}\.[a-z0-9.]+)$/.test(u.pathname)) badPaths.push(`${path}: ${u.pathname}`)
    }
  })
  const done = r => { if (net.delete(r)) { inflight--; last = Date.now() } }
  ctx.on('requestfinished', done)
  ctx.on('requestfailed', done)
  await page.goto(BASE + (path.startsWith('/c/') || path === '/' || path === '/about' ? path : PREFIX + path))
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
  check(csRows === 3, `PDXK: credible-set table has ${csRows} set rows (the block's 32 records are in sets 1, 2 and 3)`)
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
  const toggle = page.locator('label:has-text("Gene track") input[type=checkbox]')
  const before = await page.locator('svg').count()
  await toggle.click()
  await page.waitForTimeout(1500)
  const after = await page.locator('svg').count()
  check(await toggle.isChecked() && after > before, `PDXK: gene track toggle adds the track (${before} -> ${after} svg elements)`)
  await toggle.click()
  await idle()
  // one range per ranged object and no header reads (the pointer says what each object is); the
  // eQTL experiment holds no trans scan and no sQTL scan, so neither is read
  const k = Object.fromEntries(['ge', 'variants', 'leafcutter', 'gwas', 'trans'].map(x => [x, count(reqs, x)]))
  check(k.ge === 1 && k.variants === 1 && k.leafcutter === 0 && k.gwas === 1 && k.trans === 0,
    `PDXK: requests eQTL results ${k.ge} (block), variants ${k.variants} (range), GWAS ${k.gwas} (window), trans ${k.trans}, sQTL results ${k.leafcutter}`)
  check(await page.getByText(/^sQTL/).count() === 0, 'PDXK: one cis scan, so no eQTL/sQTL tab switch')
  check(errors.length === 0, `PDXK: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- PDXK in the sQTL experiment, and in each trans experiment ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${SQTL}/gene/PDXK`)
  await idle()
  const sqtlText = await page.locator('text=/tested introns/').first().innerText()
  check(/tested introns/.test(sqtlText) && count(reqs, 'leafcutter') === 1 && count(reqs, 'ge') === 0,
    `PDXK (sQTL experiment): the introns span is the one sQTL read, no eQTL read: "${sqtlText.slice(0, 70)}"`)
  // PDXK's eQTL phenotype has 115 trans rows and its 4 introns 404 (search index n_trans)
  for (const [where, n, what] of [[TRANS_E, 115, 'eQTL'], [TRANS_S, 404, 'sQTL']]) {
    await page.goto(`${BASE}${where}/gene/PDXK`)
    const ok = await page.locator(`[data-trans-total="${n}"]`).waitFor({ timeout: 60_000 }).then(() => true, () => false)
    check(ok, `PDXK (${what} trans experiment): the trans table has its ${n} rows`)
  }
  check(errors.length === 0, `PDXK sQTL/trans: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- MICAL3 (chr22): requests per intron choice, show all, sQTL CSV ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${SQTL}/gene/MICAL3`)
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
  check(count(reqs, 'ge') === 0, `MICAL3: the sQTL experiment reads no eQTL block (${count(reqs, 'ge')})`)
  await page.goto(`${BASE}${TRANS_S}/gene/MICAL3`)
  const t = await page.locator('[data-trans-total="5436"]').waitFor({ timeout: 60_000 }).then(() => true, () => false)
  check(t, 'MICAL3 (sQTL trans experiment): 5,436 rows (search index n_trans over its 33 introns)')
  check(errors.length === 0, `MICAL3: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- SMARCB1: the GWAS panel with shared variants, and its trans rows in each trans experiment ----
{
  const { ctx, page, errors, idle } = await open('/gene/SMARCB1')
  await idle()
  const dots = await page.locator('.plot-host svg').nth(1).locator('g[aria-label="dot"] > *').count()
  check(dots > 0, `SMARCB1: the QTL-versus-GWAS panel draws ${dots} shared variants`)
  await page.goto(`${BASE}${TRANS_E}/gene/SMARCB1`)
  await page.waitForSelector('[data-trans-total="895"]', { timeout: 60_000 })
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 30_000 }), page.locator('[data-trans-total] button:has-text("CSV")').first().click()])
  const csv = readFileSync(await dl.path(), 'utf8').trim().split('\n')
  check(csv.length === 896, `SMARCB1 (eQTL trans experiment): trans CSV has ${csv.length - 1} rows, header ${csv[0]}`)
  await page.goto(`${BASE}${TRANS_S}/gene/SMARCB1`)
  const s933 = await page.locator('[data-trans-total="933"]').waitFor({ timeout: 60_000 }).then(() => true, () => false)
  check(s933, 'SMARCB1 (sQTL trans experiment): 933 rows')
  check(errors.length === 0, `SMARCB1: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- GUSBP11: sQTL only; opens on its first significant intron, and is untested in the eQTL scan ----
{
  const { ctx, page, errors, reqs, idle } = await open(`${SQTL}/gene/GUSBP11`)
  await idle()
  const rows = await intronRows(page)
  const sel = rows.findIndex(r => r.selected)
  const header = await locusHeader(page)
  check(sel === 0 && nVariants(header) > 0, `GUSBP11: opens with the first significant intron selected (${rows.length} listed): "${header}"`)
  check(count(reqs, 'ge') === 0 && count(reqs, 'leafcutter') === 1 && count(reqs, 'variants') === 1,
    `GUSBP11: requests eQTL ${count(reqs, 'ge')}, sQTL ${count(reqs, 'leafcutter')}, variants ${count(reqs, 'variants')}`)
  await page.goto(`${BASE}${PREFIX}/gene/GUSBP11`)
  const untested = await page.getByText(/was not tested for QTL/).first().waitFor({ timeout: 30_000 }).then(() => true, () => false)
  check(untested, 'GUSBP11 (eQTL experiment): not tested in this experiment')
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

// ---- landing, collection page, experiment page, genes list, region, About ----
{
  const { ctx, page, errors, idle } = await open('/', { ready: false })
  await page.locator('a[href^="/c/"]').first().waitFor({ timeout: 30_000 })
  check((await page.locator('a[href^="/c/"]').count()) === 3, 'Landing: one row per collection, studies and GWAS (TOPCHeF, GTEx, DCM)')
  await page.goto(`${BASE}/c/${coll('topchef_ge_cis')}`)
  await page.getByText('topchef_ge_cis', { exact: true }).waitFor({ timeout: 30_000 })
  await idle()
  const rows = await page.locator('table:has(th:has-text("Significance")) tbody tr').count()
  const links = await page.locator('table:has(th:has-text("Significance")) tbody a[href*="/e/"]').count()
  check(rows === 4 && links === 4, `Collection: the TOPCHeF experiments table lists its 4 scans, each a link (${rows}, ${links})`)
  await page.goto(`${BASE}${PREFIX}`)
  const top = await page.locator('section:has(h2:has-text("Top loci")) [data-genome-track]').first()
    .waitFor({ timeout: 30_000 }).then(() => true, () => false)
  check(top, 'Experiment: the Top loci genome track renders from the index\'s lead columns')
  await page.goto(`${BASE}${PREFIX}/genes`)
  await page.locator('table tbody tr').first().waitFor({ timeout: 30_000 })
  const meta = await page.locator('h1 + *').first().innerText().catch(() => '')
  check((await page.locator('table tbody tr').count()) > 0, `Genes: the eGene list renders (${meta.trim()})`)
  await page.goto(`${BASE}${PREFIX}/region/chr21:43000000-44000000`)
  const pdxk = await page.getByText('PDXK', { exact: true }).first().waitFor({ timeout: 30_000 }).then(() => true, () => false)
  check(pdxk, 'Region chr21:43-44 Mb lists PDXK')
  await page.goto(`${BASE}/about`)
  await page.getByText('This store').first().waitFor({ timeout: 30_000 })
  const text = await page.locator('body').innerText()
  check(/Experiments\s+7/.test(text) && /EiFob05aCWgVU_B_Ae0cypnQut3cxUP1/.test(text), 'About: the store\'s 7 experiments and its sequence collection')
  check(errors.length === 0, `Landing/Collection/Genes/Region/About: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- GTEx: no GWAS until one is chosen; against DCM, on another catalog, matched by position and alleles ----
{
  const { ctx, page, errors, idle } = await open(GTEX, { ready: false })
  const topG = await page.locator('section:has(h2:has-text("Top loci")) [data-genome-track]').first().waitFor({ timeout: 30_000 }).then(() => true, () => false)
  check(topG, 'GTEx eQTL experiment: Top loci renders')
  await page.goto(`${BASE}${GTEX}/gene/PDXK`)
  await page.waitForSelector(READY, { timeout: 60_000 })
  check((await page.getByText('Choose a GWAS to compare it with this locus.').count()) === 1, 'GTEx PDXK without a GWAS: the QTL-versus-GWAS panel asks for one')
  await page.goto(`${BASE}${GTEX}/gwas/dcm_jurgens2024_biobanks/gene/SMARCB1`)
  await page.waitForSelector(READY, { timeout: 60_000 })
  await idle()
  const dots = await page.locator('.plot-host svg').nth(1).locator('g[aria-label="dot"] > *').count()
  check(dots > 0, `GTEx SMARCB1 against DCM: the QTL-versus-GWAS panel draws ${dots} shared variants`)
  const h4 = await page.getByText('H4 · shared causal variant').count()
  check(h4 === 1, 'GTEx SMARCB1 against DCM: coloc.abf runs')
  check(errors.length === 0, `GTEx: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

await browser.close()

check(badPaths.length === 0, `every /data/ request is a pointer or an immutable object (${dataTotal} seen)${badPaths.length ? `: ${badPaths.slice(0, 3).join(', ')}` : ''}`)
const g = spawnSync('grep', ['-rn', 'manifest\\|pack-decode\\|trans-pack\\|variant-pack\\|lib/pack\'\\|slope', 'src'], { cwd: UI, encoding: 'utf8' })
check(g.status === 1 && !g.stdout.trim(), `dead v0/v1 code grep (packs, slope) finds nothing in ui/src${g.stdout.trim() ? `: ${g.stdout.trim().split('\n').slice(0, 3).join(' | ')}` : ''}`)

const failed = results.filter(r => !r[0])
console.log(failed.length ? `smoke: ${failed.length} of ${results.length} failed` : `smoke: all ${results.length} checks passed`)
process.exit(failed.length ? 1 : 0)
