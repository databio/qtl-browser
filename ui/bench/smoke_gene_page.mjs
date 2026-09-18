// Plan 5 step 10 smoke check against the local preview: the gene page on packs. Adapted from Plan 2's
// handoff/smoke.mjs (FLNC panels and cis table), plus request counts per tab and intron, the
// "show all tested introns" toggle, double credible-set membership, chrX, the sQTL CSV, and dead code.
import { chromium } from '/home/nsheff/Dropbox/workspaces/assistant/tasks/qtl-browser/ui/node_modules/playwright/index.mjs'
import { readFileSync } from 'node:fs'
import { spawnSync } from 'node:child_process'

const BASE = 'http://localhost:4173'
const REPO = '/home/nsheff/Dropbox/workspaces/assistant/tasks/qtl-browser'
const READY = '.plot-host:not(.invisible) svg g[aria-label="dot"] > *'
const FORBIDDEN = /\/(gene_detail|gwas_dcm|cis_eqtl_nominal|cis_sqtl_nominal)\//
// the note is built from manifest.precision now, so match its shape, not one set of numbers
const ROUNDING = /are rounded \(p within [\d.]+%\)\. Exact values: Zenodo\./
const CIS_COLS = ['position', 'rsid', 'A1', 'A2', 'tss_distance', 'af', 'ma_samples', 'ma_count', 'pval_nominal', 'slope', 'slope_se', 'pip', 'cs_id']
const LINC_INTRON = 'chr2:10844505:10844753:clu_48622_+:ENSG00000271952.2'
// packs live flat under immutable/ as <stem>.<chr>.<sha16>.<ext> (SPEC section 3)
const PACK = { eqtl: /\/immutable\/eqtl\./, sqtl: /\/immutable\/sqtl\./, variants: /\/immutable\/variants\./, gwas: /\/immutable\/gwas\./ }
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }
const forbiddenSeen = []
let dataTotal = 0

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
      if (FORBIDDEN.test(u.pathname)) forbiddenSeen.push(`${path}: ${u.pathname}`)
    }
  })
  const done = r => { if (net.delete(r)) { inflight--; last = Date.now() } }
  ctx.on('requestfinished', done)
  ctx.on('requestfailed', done)
  await page.goto(BASE + path)
  if (ready) await page.waitForSelector(READY, { timeout: 60_000 })
  // no request in flight and none started or ended for `quiet` ms
  const idle = async (quiet = 2000) => {
    const t0 = Date.now()
    while (Date.now() - t0 < 30_000) { if (inflight === 0 && Date.now() - last >= quiet) return true; await page.waitForTimeout(100) }
    return false
  }
  return { ctx, page, errors, reqs, idle }
}

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
/** "intron <id> · N variants" in the sQTL locus header. */
const locusHeader = page => page.evaluate(() => [...document.querySelectorAll('span')].map(s => s.textContent.trim())
  .find(t => /^intron .* · [\d,]+ variants$/.test(t)) ?? null)
const nVariants = h => (h ? Number(/· ([\d,]+) variants/.exec(h)[1].replace(/,/g, '')) : null)
async function downloadCsv(page) {
  const [dl] = await Promise.all([page.waitForEvent('download', { timeout: 30_000 }), page.locator('button:has-text("CSV")').first().click()])
  return readFileSync(await dl.path(), 'utf8').trim().split('\n')
}

// ---- FLNC, Plan 2's list ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/FLNC')
  await page.waitForFunction(() => document.querySelectorAll('.plot-host svg').length >= 2, null, { timeout: 30_000 })
  check((await page.locator('.plot-host svg').count()) >= 2, 'FLNC: locus plot and LocusCompare both drawn')
  const csRows = await page.locator('table:has(th:has-text("Top PIP")) tbody tr').count()
  check(csRows >= 1, `FLNC: credible-set table has ${csRows} set row(s)`)
  await cisTable(page).locator('tbody tr').first().waitFor({ timeout: 30_000 })
  const total = await pagerText(page)
  check(/of 6,050/.test(total), `FLNC: cis table pager reads "${total}" (FLNC has 6,050 variants)`)
  const asc = await cisRows(page)
  await cisTable(page).locator('th button:text-is("p")').click()
  const desc = await cisRows(page)
  check(asc[0][0] !== desc[0][0], `FLNC: sorting by p flips the first row (${asc[0][1]} -> ${desc[0][1]})`)
  await cisTable(page).locator('th button:text-is("p")').click()
  await page.selectOption('select[title="Nominal p-value threshold"]', '1e-5')
  const filtered = await waitPager(page, t => !/of 6,050/.test(t))
  check(!/of 6,050/.test(filtered), `FLNC: p <= 1e-5 filter narrows the table: "${filtered}"`)
  await page.selectOption('select[title="Nominal p-value threshold"]', '')
  await waitPager(page, t => /of 6,050/.test(t))
  await page.fill('input[placeholder="rsID or position"]', 'rs73238147')
  const searched = await waitPager(page, t => /of 1$/.test(t) || /^1–1 /.test(t))
  check(/of 1$/.test(searched), `FLNC: rsID search rs73238147 finds one row: "${searched}"`)
  await page.fill('input[placeholder="rsID or position"]', '')
  await waitPager(page, t => /of 6,050/.test(t))
  await page.locator('button:text-is("Next")').first().click()
  const paged = await waitPager(page, t => t.startsWith('11–'))
  check(paged.startsWith('11–20'), `FLNC: paging Next shows "${paged}"`)
  const csv = await downloadCsv(page)
  check(csv.length === 6051 && csv[0] === CIS_COLS.join(','), `FLNC: CSV export has ${csv.length - 1} rows, header ${csv[0]}`)
  check(await page.getByText(ROUNDING).first().isVisible(), 'FLNC: rounding note by the CSV button')
  // the toggle is controlled and reaches the header through the parent's state, so its checked
  // state can lag the click by a render: click once, then poll, and count any extra click needed
  const toggle = page.locator('label:has-text("Gene track") input[type=checkbox]')
  const setToggle = async want => {
    for (let clicks = 1; clicks <= 3; clicks++) {
      await toggle.click()
      for (let i = 0; i < 20; i++) { if ((await toggle.isChecked()) === want) return clicks; await page.waitForTimeout(100) }
    }
    return null
  }
  const before = await page.locator('svg').count()
  const on = await setToggle(true)
  await page.waitForTimeout(1500)
  const after = await page.locator('svg').count()
  check(on === 1 && after > before, `FLNC: gene track toggle adds the track (${before} -> ${after} svg elements; clicks needed ${on})`)
  const off = await setToggle(false)
  check(off === 1, `FLNC: gene track toggle turns off again (clicks needed ${off})`)
  await idle()
  const kinds = { block: reqs.filter(p => PACK.eqtl.test(p)).length, variants: reqs.filter(p => PACK.variants.test(p)).length,
    gwas: reqs.filter(p => PACK.gwas.test(p)).length, intron: reqs.filter(p => PACK.sqtl.test(p)).length }
  check(kinds.block === 1 && kinds.variants === 1 && kinds.gwas === 1 && kinds.intron === 0,
    `FLNC: pack requests eQTL block ${kinds.block}, variants ${kinds.variants}, GWAS ${kinds.gwas}, intron ${kinds.intron}`)
  const d = await drawnCount(page)
  await page.getByText(/^sQTL/).first().click()
  await waitDrawn(page, d)
  const sqtlText = await page.locator('text=/tested introns/').first().innerText()
  check(/tested introns/.test(sqtlText), `FLNC: sQTL tab renders from the pack details: "${sqtlText.slice(0, 80)}"`)
  check(errors.length === 0, `FLNC: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- SYNPO2L: requests per tab switch and intron choice ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/SYNPO2L')
  await idle()
  const n0 = reqs.length
  let d = await drawnCount(page)
  await page.getByText(/^sQTL/).first().click()
  await waitDrawn(page, d); await idle()
  const n1 = reqs.length
  const introns = n => reqs.slice(n).every(p => PACK.sqtl.test(p))
  check(n1 - n0 === 1 && introns(n0), `SYNPO2L: eQTL tab -> sQTL tab sends ${n1 - n0} /data/ request(s): ${reqs.slice(n0).join(', ')}`)
  const rows = await intronRows(page)
  const first = rows.findIndex(r => r.selected), other = rows.findIndex(r => !r.selected)
  d = await drawnCount(page)
  await intronTable(page).locator('tbody tr').nth(other).locator('td').first().click()
  await waitDrawn(page, d); await idle()
  const n2 = reqs.length
  check(n2 - n1 === 1 && introns(n1) && (await intronRows(page))[other].selected, `SYNPO2L: choosing another intron (${rows[other].text[1]}) sends ${n2 - n1}: ${reqs.slice(n1).join(', ')}`)
  d = await drawnCount(page)
  await intronTable(page).locator('tbody tr').nth(first).locator('td').first().click()
  await waitDrawn(page, d); await idle()
  const n3 = reqs.length
  check(n3 === n2, `SYNPO2L: back to the first intron (${rows[first].text[1]}) sends ${n3 - n2}`)
  d = await drawnCount(page)
  await page.getByText(/^eQTL$/).first().click()
  await waitDrawn(page, d); await idle()
  const n4 = reqs.length
  check(n4 === n3, `SYNPO2L: back to the eQTL tab sends ${n4 - n3}`)
  check(errors.length === 0, `SYNPO2L: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- CAMK2D: show all tested introns, a non-significant intron's locus, sQTL CSV ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/CAMK2D?tab=sqtl')
  const label = (await page.locator('label:has-text("Show all")').first().innerText()).trim()
  const N = Number(/Show all ([\d,]+) tested/.exec(label)[1].replace(/,/g, ''))
  const toggle = showAll(page)
  const onAtStart = await toggle.isChecked()
  const sig = await intronRows(page)
  await toggle.check()
  await page.waitForTimeout(300)
  const all = await intronRows(page)
  check(!onAtStart && all.length === N && N > sig.length, `CAMK2D: toggle off at start; "${label}" lists ${all.length} introns (${sig.length} significant)`)
  const key = r => `${r.text[0]} ${r.text[1]}`
  const sigKeys = new Set(sig.map(key))
  const idx = all.findIndex(r => !sigKeys.has(key(r)))
  const d = await drawnCount(page)
  const before = reqs.length
  await intronTable(page).locator('tbody tr').nth(idx).locator('td').first().click()
  await waitDrawn(page, d); await idle()
  const header = await locusHeader(page)
  const selected = (await intronRows(page))[idx].selected
  check(selected && nVariants(header) > 0 && reqs.length - before === 1,
    `CAMK2D: non-significant intron ${key(all[idx])} (perm p ${all[idx].text[4]}) draws its locus: "${header}", ${reqs.length - before} request`)
  await cisTable(page).locator('tbody tr').first().waitFor({ timeout: 30_000 })
  const csv = await downloadCsv(page)
  const want = ['phenotype_id', ...CIS_COLS].join(',')
  check(csv[0] === want && csv.length - 1 === nVariants(header), `CAMK2D: sQTL CSV has ${csv.length - 1} rows, header ${csv[0]}`)
  check(await page.getByText(ROUNDING).first().isVisible(), 'CAMK2D: rounding note by the sQTL CSV button')
  check(errors.length === 0, `CAMK2D: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- AL031282.2: sQTL only, opens on its first significant intron ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/AL031282.2')
  await idle()
  const rows = await intronRows(page)
  const on = await showAll(page).isChecked()
  const sel = rows.findIndex(r => r.selected)
  const header = await locusHeader(page)
  check(!on && sel === 0 && nVariants(header) > 0, `AL031282.2: opens on its sQTL tab with the first significant intron selected (row ${sel + 1} of ${rows.length} listed, toggle ${on ? 'on' : 'off'}): "${header}"`)
  const k = kind => reqs.filter(p => PACK[kind].test(p)).length
  check(k('eqtl') === 1 && k('variants') === 1 && k('gwas') === 1 && k('sqtl') === 1,
    `AL031282.2: pack requests eQTL block ${k('eqtl')}, variants ${k('variants')}, GWAS ${k('gwas')}, intron ${k('sqtl')}`)
  check(errors.length === 0, `AL031282.2: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- LINC01954: the intron whose credible sets share 12 variants ----
{
  const { ctx, page, errors, idle } = await open('/gene/LINC01954?tab=sqtl')
  const toggle = showAll(page)
  if (!(await toggle.isChecked())) { await toggle.check(); await page.waitForTimeout(300) }
  const all = await intronRows(page)
  const idx = all.findIndex(r => r.text[1].startsWith('10,844,505–10,844,753'))
  if (idx < 0) check(false, `LINC01954: intron ${LINC_INTRON} is not listed`)
  else {
    if (!all[idx].selected) {
      const d = await drawnCount(page)
      await intronTable(page).locator('tbody tr').nth(idx).locator('td').first().click()
      await waitDrawn(page, d)
    }
    await idle()
    const groups = page.locator('table:has(th:has-text("Top PIP")) > tbody > tr[aria-expanded]')
    await groups.first().waitFor({ timeout: 30_000 })
    const ng = await groups.count()
    for (let i = 0; i < ng; i++) if ((await groups.nth(i).getAttribute('aria-expanded')) !== 'true') await groups.nth(i).locator('td').first().click()
    await page.waitForTimeout(300)
    const members = await page.locator('table:has(th:has-text("Top PIP"))').first().evaluate(t => {
      const out = []
      let set = null
      for (const tr of t.tBodies[0].rows) {
        if (tr.hasAttribute('aria-expanded')) set = tr.cells[1].textContent.trim()
        // :scope keeps the match inside the expanded row's own member table (not its header row)
        else for (const r of tr.querySelectorAll(':scope table > tbody > tr')) out.push([set, [...r.cells].slice(0, 3).map(c => c.textContent.trim()).join(' ')])
      }
      return out
    })
    const sets = new Map()
    for (const [s, v] of members) sets.set(v, new Set([...(sets.get(v) ?? []), s]))
    const doubles = [...sets.entries()].filter(([, s]) => s.size > 1)
    check(doubles.length === 12, `LINC01954 ${LINC_INTRON}: ${ng} sets, ${members.length} member rows; ${doubles.length} variants listed under two sets (SPEC: 12), e.g. ${doubles[0]?.[0]} in ${[...(doubles[0]?.[1] ?? [])].join(' and ')}`)
  }
  check(errors.length === 0, `LINC01954: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- chrX: no GWAS request, LocusCompare empty state ----
{
  const { ctx, page, errors, reqs, idle } = await open('/gene/AC244197.2')
  await idle()
  const empty = await page.getByText('No variants in this window are present in the DCM GWAS.').count()
  const gwas = reqs.filter(p => PACK.gwas.test(p)).length
  check(empty === 1 && gwas === 0, `AC244197.2 (chrX): LocusCompare empty state shown ${empty}, GWAS requests ${gwas}`)
  check(errors.length === 0, `AC244197.2: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

// ---- About: the corrected definitions and the GWAS count from manifest.packs ----
{
  const { ctx, page, errors, idle } = await open('/about', { ready: false })
  await page.getByText('DCM GWAS variants').first().waitFor({ timeout: 30_000 })
  await idle()
  const text = await page.locator('body').innerText()
  check(/listed under both in the credible-set table/.test(text) && /its per-variant nominal statistics/.test(text) && /12,504,079/.test(text),
    'About: double-membership and all-intron sentences, DCM GWAS variants 12,504,079')
  check(errors.length === 0, `About: no console errors${errors.length ? ` (${errors.slice(0, 3).join(' | ')})` : ''}`)
  await ctx.close()
}

await browser.close()

check(forbiddenSeen.length === 0, `no gene page requested gene_detail/, gwas_dcm/, or cis_*_nominal/ (${dataTotal} /data/ requests seen)${forbiddenSeen.length ? `: ${forbiddenSeen.slice(0, 3).join(', ')}` : ''}`)
const g = spawnSync('grep', ['-rn', 'gene_detail\\|gwas_dcm/\\|nominalFile\\|geneDetail\\|packEnabled\\|VITE_PACK', 'ui/src'], { cwd: REPO, encoding: 'utf8' })
check(g.status === 1 && !g.stdout.trim(), `dead code grep finds nothing in ui/src${g.stdout.trim() ? `: ${g.stdout.trim().split('\n').slice(0, 3).join(' | ')}` : ''}`)

const failed = results.filter(r => !r[0])
console.log(failed.length ? `smoke: ${failed.length} of ${results.length} failed` : `smoke: all ${results.length} checks passed`)
process.exit(failed.length ? 1 : 0)
