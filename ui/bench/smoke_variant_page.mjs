// Plan 9 step 11 smoke check against the local preview: the variant page on packs.
// Checks the three lists, the trans table, the cis scan, the request counts the plan gates on,
// and that no parquet is read on any variant page.
import { chromium } from '/home/nsheff/Dropbox/workspaces/assistant/tasks/qtl-browser/ui/node_modules/playwright/index.mjs'

const BASE = 'http://localhost:4173'
const LISTS = '[data-variant-lists="1"]'
const results = []
const check = (ok, msg) => { results.push([ok, msg]); console.log(`${ok ? 'PASS' : 'FAIL'} ${msg}`) }

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

/** The data requests the plan gates on: not the manifest, and not the three files boot fetches
 *  whole (search_index, the GWAS index, and the variant index). */
const gated = reqs => reqs.filter(p => !/manifest\.json$|search_index|variant_index|gwas_index/.test(p))
// packs live flat under immutable/ as <stem>.<chr>.<sha16>.<ext> (SPEC section 3)
const kind = p => /\/immutable\/hits\./.test(p) ? 'hits' : /rsid_index/.test(p) ? 'rsid_index'
  : /\/immutable\/variants\./.test(p) ? 'variants' : /\/immutable\/eqtl\./.test(p) ? 'eqtl'
  : /\/immutable\/sqtl\./.test(p) ? 'sqtl' : /\.parquet/.test(p) ? 'parquet' : 'other'
const kinds = reqs => gated(reqs).map(kind).sort().join(',')
const table = (page, head) => page.locator(`table:has(th:has-text("${head}"))`)
const rowsOf = loc => loc.locator('tbody tr').evaluateAll(trs => trs.map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim())))

// ---- rs10824026: the paper's variant, by rsID ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs10824026')
  await page.waitForSelector(`${LISTS} [data-trans-total]`, { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === 'hits,rsid_index,variants',
    `rs10824026: 3 gated cold requests, one each of rsid_index, variants, hits (got ${kinds(reqs) || 'none'})`)
  check(!gated(reqs).some(p => /\.parquet/.test(p)), `rs10824026: no parquet request (${gated(reqs).filter(p => /\.parquet/.test(p)).join(', ') || 'none'})`)
  const head = await page.locator('h1').first().innerText()
  check(head.trim() === 'rs10824026', `rs10824026: header title is "${head.trim()}"`)
  const pos = await page.getByText('chr10:73,661,450 (GRCh38)').count()
  check(pos === 1, 'rs10824026: position row reads chr10:73,661,450 (GRCh38)')
  // the hits frame holds exactly one row for this variant, a trans sQTL: it leads nothing and is in
  // no credible set, which is what the parquet page showed too (genes, splice_phenotypes and
  // credible_sets all have 0 rows at chr10:73,661,450)
  const noLead = await page.getByText('Not the lead variant for any gene or splice phenotype.').count()
  const noCs = await page.getByText('Not in any credible set.').count()
  check(noLead === 1 && noCs === 1, `rs10824026: leads and credible sets are both empty, as the parquet page had them (${noLead}, ${noCs})`)
  const transTotal = await page.locator('[data-trans-total]').first().getAttribute('data-trans-total')
  check(transTotal === '1', `rs10824026: trans table reports ${transTotal} row (expected 1)`)

  // the scan: two more requests, one eQTL span and one sQTL span
  const before = reqs.length
  const label = await page.locator('[data-scan-button]').innerText()
  check(/^Scan cis windows \(861 KB \+ splicing\)$/.test(label.trim()), `rs10824026: scan button reads "${label.trim()}" (the eQTL span; splicing is only known once it arrives)`)
  await page.locator('[data-scan-button]').click()
  await page.waitForSelector('[data-scan-ready="1"]', { timeout: 120_000 })
  await idle()
  const scanReqs = reqs.slice(before)
  check(kinds(scanReqs) === 'eqtl,sqtl', `rs10824026 scan: 2 requests, one eQTL span and one sQTL span (got ${kinds(scanReqs) || 'none'})`)
  const eTitle = await page.getByRole('heading', { name: /^Expression \(/ }).innerText()
  check(eTitle.trim() === 'Expression (42)', `rs10824026 scan: ${eTitle.trim()} (the raw chr10 file has 42 genes at this position)`)
  const sTitle = await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()
  const toggle = await page.locator('button:has-text("Show all")').innerText()
  check(/^Show all 239 tested introns$/.test(toggle.trim()),
    `rs10824026 scan: splicing shows ${sTitle.trim()} with "${toggle.trim()}" (the raw file has 239 introns)`)
  await page.locator('button:has-text("Show all")').click()
  const sAll = await page.getByRole('heading', { name: /^Splicing \(/ }).innerText()
  check(sAll.trim() === 'Splicing (239)', `rs10824026 scan: the toggle shows ${sAll.trim()}`)
  check(errors.length === 0, `rs10824026: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- the same variant by chr:pos: one request fewer ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/chr10:73661450')
  await page.waitForSelector(`${LISTS} [data-trans-total]`, { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === 'hits,variants', `chr10:73661450: 2 gated cold requests, variants then hits (got ${kinds(reqs) || 'none'})`)
  const head = await page.locator('h1').first().innerText()
  check(head.trim() === 'rs10824026', `chr10:73661450: resolves to ${head.trim()}`)
  check(errors.length === 0, `chr10:73661450: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs141809548: trans-only, outside every cis window ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs141809548')
  await page.waitForSelector(`${LISTS} [data-trans-total="9"]`, { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === 'hits,rsid_index,variants', `rs141809548: 3 gated cold requests (got ${kinds(reqs) || 'none'})`)
  check((await page.locator('[data-scan-button]').count()) === 0, 'rs141809548: no scan button (outside every cis window)')
  const outside = await page.getByText('outside every cis window', { exact: false }).count()
  check(outside === 3, `rs141809548: the outside-cis message stands in for all three cis sections (${outside} of 3)`)
  const alleles = await page.getByText('ATGTCT / A', { exact: true }).count()
  check(alleles === 1, 'rs141809548: alleles read ATGTCT / A (an indel the trans-only section still reports)')
  check(errors.length === 0, `rs141809548: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- a chrX variant ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs1204407')
  await page.waitForSelector(`${LISTS} [data-trans-total]`, { timeout: 60_000 })
  await idle()
  check(kinds(reqs) === 'hits,rsid_index,variants', `rs1204407 (chrX): 3 gated cold requests (got ${kinds(reqs) || 'none'})`)
  const pos = await page.getByText('chrX:100,649,875 (GRCh38)').count()
  check(pos === 1, 'rs1204407: position row reads chrX:100,649,875 (GRCh38)')
  const leads = await rowsOf(table(page, 'Perm p'))
  check(leads.some(r => r.includes('TSPAN6')), `rs1204407: lead of TSPAN6 (${leads.length} lead row(s))`)
  check(errors.length === 0, `rs1204407: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- rs1: a real rsID the study does not hold ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/rs1')
  await page.waitForSelector('text=is not among the variants tested', { timeout: 60_000 })
  await idle()
  // rs1 is below the first rsID block's first record (rs3), so SPEC section 14 step 1 answers
  // "not held" from the startup file alone and sends no request at all
  check(gated(reqs).length === 0, `rs1: answered from the startup file with no request (got ${kinds(reqs) || 'none'})`)
  check((await page.getByText('Look it up in dbSNP').count()) === 1, 'rs1: offers the dbSNP link')
  check(errors.length === 0, `rs1: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

// ---- a malformed id: no request at all ----
{
  const { ctx, page, errors, reqs, idle } = await open('/variant/not-a-variant')
  await page.waitForSelector('text=is not among the variants tested', { timeout: 60_000 })
  await idle()
  check(gated(reqs).length === 0, `not-a-variant: no data request (got ${gated(reqs).join(', ') || 'none'})`)
  check(errors.length === 0, `not-a-variant: no console errors (${errors.slice(0, 2).join(' | ')})`)
  await ctx.close()
}

await browser.close()
const failed = results.filter(([ok]) => !ok)
console.log(`\n${results.length - failed.length}/${results.length} checks passed`)
if (failed.length) process.exit(1)
