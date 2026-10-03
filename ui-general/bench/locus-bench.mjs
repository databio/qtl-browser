#!/usr/bin/env node
/**
 * Gene-page cost harness: loads fixed gene pages in headless Chromium (Playwright) and records
 * every request until the locus plot is drawn, then until the network is idle.
 *
 *   node bench/locus-bench.mjs --target live|preview|rehearsal --label baseline
 *        [--runs 3] [--scenario cold,warm,nav,pair] [--pages FLNC-eqtl,...] [--note "network ..."]
 *        [--nav-mode popstate|click] [--query pack=1] [--latency 90] [--manifest immutable/manifest.<sha16>.json]
 *
 * See bench/README.md for the targets, scenarios, the readiness selector contract, and outputs.
 */
import { chromium } from 'playwright'
import { createRequire } from 'node:module'
import { execSync } from 'node:child_process'
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs'
import os from 'node:os'
import { dirname, join, relative } from 'node:path'
import { fileURLToPath } from 'node:url'

const HERE = dirname(fileURLToPath(import.meta.url))
const UI = join(HERE, '..')
const REPO = join(UI, '..')
const require = createRequire(import.meta.url)

// The locus scatter's host loses `invisible` once drawn and holds Plot's dot marks. LocusCompare
// also uses .plot-host but mounts only after the scatter is ready, so this cannot fire early.
const READY_SELECTOR = '.plot-host:not(.invisible) svg g[aria-label="dot"] > *'
const READY_TIMEOUT = 120_000
const IDLE_QUIET = 2_000       // no new request (and none in flight) for this long
const IDLE_CAP = 30_000        // after ready
const VIEWPORT = { width: 1280, height: 900 }
const CATEGORIES = ['data', 'app', 'duckdb-cdn', 'duckdb-ext', 'other']
const SCENARIOS = ['cold', 'warm', 'nav', 'pair']
const PREVIEW_HELP = 'cd ui && VITE_DATA_BASE= npm run build && npm run preview'
const USAGE = 'usage: node bench/locus-bench.mjs --target live|preview|rehearsal --label <name> [--runs 3] [--scenario cold,warm,nav,pair] [--pages FLNC-eqtl,...] [--note "..."] [--nav-mode popstate|click] [--query pack=1] [--latency <ms>] [--manifest <key>] [--self-check-page <id>]'
// Any manifest name: `manifest.json`, or a staged copy `immutable/manifest.<sha16>.json` (rehearsal).
// qtlb v1 (SPEC.md section 2): store.json and the pointer documents count as `manifest` too.
const MANIFEST_RE = /\/(manifest(\.[0-9a-f]{16})?\.json|store\.json|(experiments|variant_catalogs|annotations)\/[\w.-]+\.json)(\?|$)/
// Data requests by what they read; gated plan metrics count some kinds separately. Every file the
// browser reads at a byte offset lives flat under immutable/ as <stem>.<sha16>.<ext> (SPEC section
// 3), so a kind is matched by its stem. Kind names are unchanged, so older results stay comparable;
// the parquet kinds are kept for the same reason, and should all read 0 now.
const DATA_KINDS = [['eqtl_pack', /\/immutable\/eqtl\./], ['variants', /\/immutable\/variants\./], ['gene_detail', /\/gene_detail\//],
  ['eqtl_nominal', /\/cis_eqtl_nominal\//], ['sqtl_nominal', /\/cis_sqtl_nominal\//], ['gwas', /\/gwas_dcm\//], ['trans', /\/trans_pairs\//],
  ['gwas_pack', /\/immutable\/gwas\./], ['sqtl_pack', /\/immutable\/sqtl\./], ['trans_pack', /\/immutable\/trans\./], ['gwas_index', /\/immutable\/gwas_index\./],
  ['hits', /\/immutable\/hits\./], ['rsid_index', /\/immutable\/rsid_index\./], ['variant_index', /\/immutable\/variant_index\./],
  ['manifest', MANIFEST_RE], ['search_index', /\/immutable\/search_index\./],
  // qtlb v1 objects are `immutable/<sha512t24u>.<ext>`: the kind is the extension. Results files
  // (.qbe) do not say their phenotype type in the name, so eQTL and sQTL reads share `results`;
  // `.arrow.zst` is the search index and the annotation's gene and exon tables; `.qbt`, `.qbg` and
  // `.qgi` count under the v0 names trans_pack, gwas_pack and gwas_index.
  ['variants', /\/immutable\/[\w-]{32}\.qbv(\?|$)/], ['results', /\/immutable\/[\w-]{32}\.qbe(\?|$)/],
  ['hits', /\/immutable\/[\w-]{32}\.qbh(\?|$)/], ['rsid_index', /\/immutable\/[\w-]{32}\.qbr(\?|$)/],
  ['variant_index', /\/immutable\/[\w-]{32}\.qbx(\?|$)/], ['arrow_object', /\/immutable\/[\w-]{32}\.arrow\.zst(\?|$)/],
  ['trans_pack', /\/immutable\/[\w-]{32}\.qbt(\?|$)/], ['gwas_pack', /\/immutable\/[\w-]{32}\.qbg(\?|$)/],
  ['gwas_index', /\/immutable\/[\w-]{32}\.qgi(\?|$)/],
  // last: every parquet a page still reads. Plan 9 removed the last of them, so this is a tripwire
  // that should stay 0.
  ['parquet_other', /\.parquet(\?|$)/]]
const KIND_NAMES = [...new Set(DATA_KINDS.map(k => k[0])), 'other']
const dataKind = url => (DATA_KINDS.find(([, re]) => re.test(url)) ?? ['other'])[0]
/** Run-wide options set once in main: extra query parameters for every page, and CDP latency. */
const OPTS = { query: null, latency: 0 }

/** A page path with the run's extra query parameters (for example pack=1). */
function withQuery(path, query) {
  if (!query) return path
  const u = new URL(path, 'http://x')
  for (const [k, v] of new URLSearchParams(query)) u.searchParams.set(k, v)
  return u.pathname + u.search
}

function die(msg, code = 1) {
  console.error(`locus-bench: ${msg}`)
  process.exit(code)
}
const log = msg => console.log(`[${new Date().toISOString().slice(11, 19)}] ${msg}`)
const sleep = ms => new Promise(r => setTimeout(r, ms))
const withTimeout = (p, ms) => Promise.race([p, sleep(ms).then(() => undefined)])

/** The data host production bundles read, from ui/.env.production. */
function productionDataBase() {
  const m = /^VITE_DATA_BASE=(.*)$/m.exec(readFileSync(join(UI, '.env.production'), 'utf8'))
  return m ? m[1].trim().replace(/\/$/, '') : null
}

const DATA_BASE = productionDataBase()
const TARGETS = {
  live: { name: 'live', origin: 'https://topchef.databio.org', dataBase: DATA_BASE },
  preview: { name: 'preview', origin: 'http://localhost:4173', dataBase: 'http://localhost:4173/data' },
  // the deploy rehearsal: a local preview of a production bundle reading the real bucket, pointed
  // at a staged manifest copy with --manifest before manifest.json itself is switched over
  rehearsal: { name: 'rehearsal', origin: 'http://localhost:4173', dataBase: DATA_BASE },
}

function parseArgs(argv) {
  const a = { runs: 3, scenario: 'cold,warm,nav', pages: null, target: null, label: null, note: null, navMode: 'popstate', query: null, latency: 0, manifest: null,
    selfCheckPage: 'FLNC-eqtl' }
  for (let i = 0; i < argv.length; i++) {
    const k = argv[i]
    const v = () => { const x = argv[++i]; if (x === undefined) die(`${k} needs a value\n${USAGE}`); return x }
    if (k === '--target') a.target = v()
    else if (k === '--label') a.label = v()
    else if (k === '--runs') a.runs = Number(v())
    else if (k === '--scenario') a.scenario = v()
    else if (k === '--pages') a.pages = v()
    else if (k === '--note') a.note = v()
    else if (k === '--nav-mode') a.navMode = v()
    else if (k === '--query') a.query = v()
    else if (k === '--latency') a.latency = Number(v())
    else if (k === '--manifest') a.manifest = v().replace(/^\/+/, '')
    else if (k === '--self-check-page') a.selfCheckPage = v()
    else if (k === '--help' || k === '-h') { console.log(USAGE); process.exit(0) }
    else die(`unknown argument ${k}\n${USAGE}`)
  }
  if (!TARGETS[a.target]) die(`--target must be one of ${Object.keys(TARGETS).join(', ')}\n${USAGE}`)
  if (!a.label || !/^[\w.-]+$/.test(a.label)) die(`--label is required (letters, digits, . _ -)\n${USAGE}`)
  if (!Number.isInteger(a.runs) || a.runs < 1) die('--runs must be a positive integer')
  a.scenarios = a.scenario.split(',').map(s => s.trim()).filter(Boolean)
  for (const s of a.scenarios) if (!SCENARIOS.includes(s)) die(`unknown scenario ${s}`)
  if (!['popstate', 'click'].includes(a.navMode)) die('--nav-mode must be popstate or click')
  if (!Number.isFinite(a.latency) || a.latency < 0) die('--latency must be a number of ms >= 0')
  if (a.target === 'rehearsal' && !a.manifest) die('--target rehearsal needs --manifest <key>, the staged manifest copy (immutable/manifest.<sha16>.json)')
  if (a.manifest && a.target !== 'rehearsal') die('--manifest applies to --target rehearsal only')
  // a pair is one measured navigation between two named pages, so the ids cannot be defaulted
  if (a.scenarios.includes('pair')) {
    if (a.scenarios.length > 1) die('--scenario pair runs on its own, not alongside cold, warm or nav')
    if ((a.pages ?? '').split(',').filter(s => s.trim()).length !== 2) die('--scenario pair needs exactly two --pages ids, for example --pages MYOZ1-eqtl,SYNPO2L-sqtl')
  }
  return a
}

// ---- request classification and recording ----------------------------------------------------

function classifier(target) {
  const dataHosts = new Set([new URL(target.dataBase).host, DATA_BASE && new URL(DATA_BASE).host].filter(Boolean))
  const origin = new URL(target.origin).origin
  return url => {
    let u
    try { u = new URL(url) } catch { return 'other' }
    if (u.protocol === 'blob:' || u.protocol === 'data:') return 'local'
    if (u.origin === origin && u.pathname.startsWith('/data/')) return 'data'
    if (dataHosts.has(u.host) && !(u.origin === origin && !u.pathname.startsWith('/data/'))) return 'data'
    if (u.origin === origin) return 'app'
    if (u.host === 'cdn.jsdelivr.net') return 'duckdb-cdn'
    if (u.host === 'extensions.duckdb.org') return 'duckdb-ext'
    return 'other'
  }
}

/** Collects the requests of one scenario at a time from every page and worker of a context. */
class Recorder {
  constructor(context, classify) {
    this.classify = classify
    this.cur = null
    context.on('request', r => this.onRequest(r))
    context.on('requestfinished', r => this.onDone(r, false))
    context.on('requestfailed', r => this.onDone(r, true))
  }

  start(name) {
    this.cur = { name, map: new Map(), list: [], pending: [], inflight: 0, last: Date.now() }
  }

  onRequest(req) {
    const c = this.cur
    if (!c) return
    const url = req.url()
    let frame = true
    try { req.frame() } catch { frame = false }     // worker requests have no frame
    const rec = {
      url, method: req.method(), resource_type: req.resourceType(), category: this.classify(url),
      range: req.headers()['range'] ?? null, has_frame: frame, wall_start: Date.now(),
    }
    rec.kind = rec.category === 'data' ? dataKind(url) : null
    c.list.push(rec)
    if (rec.category === 'local') return             // blob:/data: URLs are not network; logged, never counted
    c.map.set(req, rec)
    c.inflight++
    c.last = Date.now()
  }

  onDone(req, failed) {
    const c = this.cur
    if (!c) return
    const rec = c.map.get(req)
    if (!rec || rec.wall_end) return
    rec.wall_end = Date.now()
    rec.failed = failed
    c.inflight--
    c.last = Date.now()
    c.pending.push((async () => {
      if (failed) rec.failure = req.failure()?.errorText ?? null
      rec.timing = req.timing()
      const resp = await withTimeout(req.response().catch(() => null), 5_000)
      if (resp) {
        rec.status = resp.status()
        const h = resp.headers()
        rec.content_length = h['content-length'] ?? null
        rec.content_range = h['content-range'] ?? null
        rec.from_service_worker = resp.fromServiceWorker()
      }
      rec.sizes = (await withTimeout(req.sizes().catch(() => null), 5_000)) ?? null
      if (!rec.range) {
        const all = await withTimeout(req.allHeaders().catch(() => null), 5_000)
        if (all?.range) rec.range = all.range
      }
    })())
  }

  /** Resolves true once nothing is in flight and no request started or ended for IDLE_QUIET ms,
   *  false at IDLE_CAP. */
  async waitIdle(cap = IDLE_CAP) {
    const c = this.cur
    const t0 = Date.now()
    while (Date.now() - t0 < cap) {
      if (c.inflight <= 0 && Date.now() - c.last >= IDLE_QUIET) return true
      await sleep(100)
    }
    return false
  }

  async stop() {
    const c = this.cur
    await Promise.all(c.pending)
    this.cur = null
    return c.list
  }
}

/** Times relative to the scenario's zero (epoch ms), bytes, and the from-cache rule. The raw
 *  fields stay in the record so the rule can be refined without re-running. */
function place(rec, zeroEpoch) {
  const t = rec.timing
  const hasT = t && t.startTime > 0
  const start = hasT ? t.startTime : rec.wall_start
  const end = hasT && t.responseEnd >= 0 ? t.startTime + t.responseEnd : (rec.wall_end ?? null)
  rec.start_ms = round1(start - zeroEpoch)
  rec.end_ms = end == null ? null : round1(end - zeroEpoch)
  const s = rec.sizes
  // Playwright derives responseBodySize as encodedDataLength minus header bytes, so a response
  // served from the HTTP cache (nothing on the wire) reads body = -headers. Transferred bytes are
  // their sum, and zero transferred bytes on a 200/206 is the from-cache mark. The plan's first
  // rule (zero body and responseStart -1) never matched in Chromium 153: cache hits report
  // responseStart ~0.1 ms. It is kept as a second condition.
  rec.bytes = s ? Math.max(0, s.responseBodySize + s.responseHeadersSize) : 0
  rec.from_cache = !!rec.from_service_worker ||
    ((rec.status === 200 || rec.status === 206) && !!s && s.responseBodySize + s.responseHeadersSize <= 0) ||
    ((rec.status === 200 || rec.status === 206) && s?.responseBodySize === 0 && t?.responseStart === -1)
}

const round1 = x => Math.round(x * 10) / 10
const sum = xs => xs.reduce((a, b) => a + b, 0)

function maxConcurrency(recs) {
  const ev = []
  for (const r of recs) {
    ev.push([r.start_ms, 1])
    ev.push([r.end_ms ?? r.start_ms, -1])
  }
  ev.sort((a, b) => a[0] - b[0] || a[1] - b[1])     // an end before a start at the same instant
  let cur = 0, max = 0
  for (const [, d] of ev) { cur += d; max = Math.max(max, cur) }
  return max
}

/** The data request that took longest in this window, and how long. A stall (Chrome's cache lock,
 *  a slow range read) shows up here and nowhere else in the summary. */
function slowestData(data) {
  let best = null
  for (const r of data) {
    if (r.end_ms == null) continue
    const ms = r.end_ms - r.start_ms
    if (!best || ms > best.ms) best = { ms: round1(ms), url: r.url }
  }
  return best
}

function windowStats(recs) {
  const net = recs.filter(r => r.category !== 'local')
  const data = net.filter(r => r.category === 'data')
  const bytes_by_category = Object.fromEntries(CATEGORIES.map(c => [c, sum(net.filter(r => r.category === c).map(r => r.bytes))]))
  const slowest = slowestData(data)
  return {
    requests_all: net.length,
    requests_data: data.length,
    data_get: data.filter(r => r.method === 'GET').length,
    data_head: data.filter(r => r.method === 'HEAD').length,
    data_other_method: data.filter(r => r.method !== 'GET' && r.method !== 'HEAD').length,
    data_206: data.filter(r => r.status === 206).length,
    data_200: data.filter(r => r.status === 200).length,
    data_other_status: data.filter(r => r.status !== 206 && r.status !== 200).length,
    data_from_cache: data.filter(r => r.from_cache).length,
    bytes_all: sum(net.map(r => r.bytes)),
    bytes_data: bytes_by_category.data,
    bytes_by_category,
    max_concurrency_data: maxConcurrency(data),
    // the longest single data request (end_ms - start_ms) and its URL
    data_max_ms: slowest?.ms ?? null,
    data_max_url: slowest?.url ?? null,
    kinds: Object.fromEntries(KIND_NAMES.map(k => {
      const q = data.filter(r => r.kind === k)
      return [k, { requests: q.length, bytes: sum(q.map(r => r.bytes)), head: q.filter(r => r.method === 'HEAD').length,
        not_206: q.filter(r => r.status !== 206).length, not_one_range: q.filter(r => !/^bytes=\d+-\d+$/.test(r.range ?? '')).length }]
    })),
    // 1 when the eQTL block and variants requests were in flight at the same time
    pack_overlap: packOverlap(data),
    // 1 when the eQTL block, variants, and GWAS pack requests were all in flight at one moment
    pack_overlap_gwas: allInFlight(data, ['eqtl_pack', 'variants', 'gwas_pack']),
    // ms from the end of the eQTL block response to the start of the first intron block request
    intron_after_block_ms: intronAfterBlock(data),
  }
}

function allInFlight(data, kinds) {
  const rs = kinds.map(k => data.find(r => r.kind === k))
  if (rs.some(r => !r || r.end_ms == null)) return null
  return Math.max(...rs.map(r => r.start_ms)) < Math.min(...rs.map(r => r.end_ms)) ? 1 : 0
}

function intronAfterBlock(data) {
  const a = data.find(r => r.kind === 'eqtl_pack'), b = data.find(r => r.kind === 'sqtl_pack')
  return a && b && a.end_ms != null ? round1(b.start_ms - a.end_ms) : null
}

function packOverlap(data) {
  const a = data.find(r => r.kind === 'eqtl_pack'), b = data.find(r => r.kind === 'variants')
  if (!a || !b || a.end_ms == null || b.end_ms == null) return null
  return a.start_ms < b.end_ms && b.start_ms < a.end_ms ? 1 : 0
}

/** Gene page marks and pack measures after the scenario's zero (page clock), in ms from zero. */
async function pageMarks(tab, zeroPerf) {
  return await tab.evaluate(z => {
    const first = name => performance.getEntriesByName(name, 'mark').find(m => m.startTime >= z) ?? null
    const hit = first('gene:hit'), detail = first('gene:detail'), drawn = first('locus:drawn')
    const measures = performance.getEntriesByType('measure').filter(m => m.name.startsWith('pack:') && m.startTime >= z)
    const total = (name, keep = () => true) => { const ms = measures.filter(m => m.name === name && keep(m)); return ms.length ? ms.reduce((a, m) => a + m.duration, 0) : null }
    // the trans table's worker wait and insert are reported apart, so the locus and GWAS totals stay comparable with earlier results
    const trans = m => m.detail?.part === 'trans', notTrans = m => !trans(m)
    return {
      summary: {
        ms_hit: hit && hit.startTime - z, ms_detail: detail && detail.startTime - z, ms_drawn: drawn && drawn.startTime - z,
        ms_hit_to_detail: hit && detail ? detail.startTime - hit.startTime : null,
        ms_hit_to_drawn: hit && drawn ? drawn.startTime - hit.startTime : null,
        pack_block_ms: total('pack:block'), pack_variants_ms: total('pack:variants'),
        pack_gwas_ms: total('pack:gwas'), pack_intron_ms: total('pack:intron'),
        pack_decode_ms: total('pack:decode'), pack_worker_wait_ms: total('pack:worker-wait', notTrans), pack_insert_ms: total('pack:insert', notTrans),
        pack_trans_ms: total('pack:trans'), pack_trans_decode_ms: total('pack:trans-decode'),
        pack_trans_wait_ms: total('pack:worker-wait', trans), pack_trans_insert_ms: total('pack:insert', trans),
      },
      marks: performance.getEntriesByType('mark').filter(m => m.startTime >= z).map(m => ({ name: m.name, t: m.startTime - z, detail: m.detail })),
      measures: measures.map(m => ({ name: m.name, t: m.startTime - z, duration: m.duration, detail: m.detail })),
    }
  }, zeroPerf ?? 0)
}

function summarizeScenario(recs, readyMs) {
  const net = recs.filter(r => r.category !== 'local')
  const beforeReady = net.filter(r => r.start_ms <= readyMs)
  const dataReady = beforeReady.filter(r => r.category === 'data')
  const ends = net.map(r => r.end_ms).filter(x => x != null)
  return {
    ms_ready: round1(readyMs),
    ms_idle: round1(Math.max(readyMs, ...ends)),
    // first data request start to the last data request end, among data requests started
    // before ready, clipped at ready
    ms_data_first_to_last_before_ready: dataReady.length
      ? round1(Math.min(readyMs, Math.max(...dataReady.map(r => r.end_ms ?? r.start_ms))) - Math.min(...dataReady.map(r => r.start_ms)))
      : null,
    to_ready: windowStats(beforeReady),
    to_idle: windowStats(net),
  }
}

// ---- scenarios -------------------------------------------------------------------------------

/** A page's readiness selector: its `ready` entry in pages.json, else the locus scatter. */
const readySelector = page => page.ready ?? READY_SELECTOR

async function waitReady(tab, selector) {
  const h = await tab.waitForFunction(
    sel => (document.querySelector(sel) ? performance.now() : 0),
    selector, { polling: 'raf', timeout: READY_TIMEOUT })
  return await h.jsonValue()
}

const symbolOf = page => decodeURIComponent(new URL(page.path, 'http://x').pathname.split('/').pop())

/** `rsid` allows an h1 that is any rsID: a variant page asked for by chr:pos renders the rsID the
 *  packs hold for that position, so its header never repeats the id in the path. */
async function headerRendered(tab, symbol, rsid = false) {
  await tab.waitForFunction(([sym, any]) => [...document.querySelectorAll('h1')]
    .some(h => h.textContent.trim() === sym || (any && /^rs\d+$/.test(h.textContent.trim()))),
  [symbol, rsid], { timeout: 10_000 })
}

/**
 * One measured scenario. `prepare` runs unrecorded; `act` runs recorded and returns the
 * scenario's zero on the page clock (0 = navigation start of the document that is current
 * once ready).
 */
async function measure({ rec, tab, name, page, target, prepare, act, after }) {
  if (prepare) await prepare()
  const load_start = os.loadavg()
  rec.start(name)
  let error = null, readyPerf = null, zeroPerf = null, timeOrigin = null, idle = null, marks = null
  try {
    zeroPerf = await act()
    // a page with `click` is two phases: wait for `pre_ready`, click, then wait for `ready`.
    // The scan button is the only user of this, so its ready time includes the page it acts on.
    if (page.click) {
      if (page.pre_ready) await waitReady(tab, page.pre_ready)
      await tab.click(page.click, { timeout: READY_TIMEOUT })
    }
    readyPerf = await waitReady(tab, readySelector(page))
    timeOrigin = await tab.evaluate(() => performance.timeOrigin)
    if (after) await after()
    idle = await rec.waitIdle()
    marks = await pageMarks(tab, zeroPerf)
  } catch (e) {
    error = String(e?.message ?? e).split('\n')[0]
    await rec.waitIdle(5_000).catch(() => {})
  }
  const recs = await rec.stop()
  const load_end = os.loadavg()
  if (timeOrigin == null) timeOrigin = await tab.evaluate(() => performance.timeOrigin).catch(() => null)
  const zeroEpoch = (timeOrigin ?? recs[0]?.wall_start ?? Date.now()) + (zeroPerf ?? 0)
  for (const r of recs) place(r, zeroEpoch)
  recs.sort((a, b) => a.start_ms - b.start_ms)
  const readyMs = readyPerf == null ? null : readyPerf - (zeroPerf ?? 0)
  const summary = readyMs == null ? null : { ...summarizeScenario(recs, readyMs), marks: marks?.summary ?? null }
  const status = error ? `ERROR ${error}` : `ready ${Math.round(readyMs)} ms, ${summary.to_ready.requests_data} data req to ready, ${summary.to_idle.requests_data} to idle${idle ? '' : ' (idle cap hit)'}` +
    (marks?.summary?.ms_hit_to_drawn != null ? `, hit->drawn ${Math.round(marks.summary.ms_hit_to_drawn)} ms, worker wait ${Math.round(marks.summary.pack_worker_wait_ms ?? 0)} ms` : '')
  log(`${target.name} ${page.id} ${name}: ${status}`)
  return { page: page.id, scenario: name, error, idle_reached: idle, ready_ms: readyMs, zero_epoch: zeroEpoch,
    loadavg_start: load_start, loadavg_end: load_end, summary, page_marks: marks?.marks ?? null, page_measures: marks?.measures ?? null, requests: recs }
}

async function newContext(browser, target) {
  const context = await browser.newContext({ viewport: VIEWPORT })
  const rec = new Recorder(context, classifier(target))
  const tab = await context.newPage()
  if (OPTS.latency) {
    // added round-trip latency on every request of the tab (the self-check reports whether the
    // DuckDB worker's range reads are delayed too)
    const cdp = await context.newCDPSession(tab)
    await cdp.send('Network.enable')
    await cdp.send('Network.emulateNetworkConditions', { offline: false, latency: OPTS.latency, downloadThroughput: -1, uploadThroughput: -1 })
  }
  return { context, rec, tab }
}

/** cold, then warm (reload) in the same context. Warm needs cold first, so cold always runs. */
async function coldWarm(browser, target, page, wantWarm) {
  const { context, rec, tab } = await newContext(browser, target)
  const url = target.origin + page.path
  const out = []
  try {
    out.push(await measure({ rec, tab, name: 'cold', page, target, act: async () => { await tab.goto(url, { waitUntil: 'commit' }); return 0 } }))
    if (wantWarm && !out[0].error) {
      out.push(await measure({ rec, tab, name: 'warm', page, target, act: async () => { await tab.reload({ waitUntil: 'commit' }); return 0 } }))
    }
  } finally {
    await context.close()
  }
  return out
}

/** In-app navigation to `page` by synthetic popstate. Throws NAV_IGNORED if the router sits still. */
const popstateAct = (tab, page) => async () => {
  const want = new URL(page.path, 'http://x')
  const t0 = await tab.evaluate(path => {
    const t = performance.now()
    history.pushState({}, '', path)
    dispatchEvent(new PopStateEvent('popstate', { state: {} }))
    return t
  }, page.path)
  const routed = await tab.waitForFunction(
    p => location.pathname === p && !document.querySelector('input[placeholder^="Filter by symbol"]'),
    want.pathname, { timeout: 5_000 }).then(() => true, () => false)
  if (!routed) throw new Error('NAV_IGNORED: React Router did not react to the synthetic popstate')
  return t0
}

/** The same navigation as a visitor makes it, from /genes: filter by symbol, click the row. */
const clickAct = (tab, page) => async () => {
  const want = new URL(page.path, 'http://x')
  const symbol = symbolOf(page)
  await tab.fill('input[placeholder^="Filter by symbol"]', symbol)
  let row = tab.locator('table tbody tr[role="link"]').filter({ has: tab.locator('td:first-child', { hasText: new RegExp(`^${symbol}$`) }) })
  if (await row.count() === 0) {
    await tab.getByText('All', { exact: true }).first().click()
    row = tab.locator('table tbody tr[role="link"]').filter({ has: tab.locator('td:first-child', { hasText: new RegExp(`^${symbol}$`) }) })
  }
  const t0 = await tab.evaluate(() => performance.now())
  await row.first().click()
  if (want.searchParams.get('tab') === 'sqtl') {
    await tab.getByText(/^sQTL/).first().click({ timeout: 60_000 })
  }
  return t0
}

/** After a navigation: the app is on the page it was sent to and its header has rendered. */
const navLanded = (tab, page) => async () => {
  const want = new URL(page.path, 'http://x')
  const variant = want.pathname.startsWith('/variant/')
  const path = await tab.evaluate(() => location.pathname + location.search)
  if (!path.startsWith(want.pathname) && !path.includes(variant ? '/variant/' : '/gene/')) throw new Error(`nav landed on ${path}`)
  await headerRendered(tab, symbolOf(page), variant)
}

/** Load /genes and let it settle, unrecorded: the engine and search index are up before the act. */
const atGenes = (tab, rec, target) => async () => {
  await tab.goto(`${target.origin}${withQuery('/genes', OPTS.query)}`, { waitUntil: 'commit' })
  await tab.waitForSelector('table tbody tr[role="link"]', { timeout: READY_TIMEOUT })
  rec.start('nav-prepare')                       // settle: the index load must not leak into the act
  await rec.waitIdle()
  await rec.stop()
}

/** Engine already running on /genes, then an in-app navigation to the gene. */
async function nav(browser, target, page, mode) {
  const { context, rec, tab } = await newContext(browser, target)
  try {
    let r = await measure({ rec, tab, name: 'nav', page, target, prepare: atGenes(tab, rec, target),
      act: mode === 'click' ? clickAct(tab, page) : popstateAct(tab, page), after: navLanded(tab, page) })
    if (r.error?.startsWith('NAV_IGNORED')) {
      log(`${page.id} nav: synthetic popstate ignored, falling back to filter + row click`)
      await context.close()
      return await nav(browser, target, page, 'click')
    }
    r.nav_mode = mode
    return r
  } finally {
    await context.close().catch(() => {})
  }
}

/**
 * Two gene pages in one context: the first is loaded cold and thrown away, then the navigation to
 * the second is the measured one. Both pages usually sit on one chromosome, so the second sends new
 * ranges to pack URLs the browser has just cached -- the shape behind the 20 s Chrome cache-lock
 * stall logged on 2026-09-04. Recorded as scenario `pair` under the id `<first>><second>`.
 */
async function pair(browser, target, first, second, mode) {
  const { context, rec, tab } = await newContext(browser, target)
  const page = { ...second, id: `${first.id}>${second.id}` }
  try {
    const prepare = async () => {
      rec.start('pair-prepare')
      await tab.goto(target.origin + first.path, { waitUntil: 'commit' })
      await waitReady(tab, readySelector(first))
      await rec.waitIdle()
      await rec.stop()
      // the click fallback navigates from the genes list, which keeps the first page's cached files
      if (mode === 'click') await atGenes(tab, rec, target)()
    }
    let r = await measure({ rec, tab, name: 'pair', page, target, prepare,
      act: mode === 'click' ? clickAct(tab, second) : popstateAct(tab, second), after: navLanded(tab, second) })
    if (r.error?.startsWith('NAV_IGNORED')) {
      log(`${page.id} pair: synthetic popstate ignored, falling back to filter + row click`)
      await context.close()
      return await pair(browser, target, first, second, 'click')
    }
    r.nav_mode = mode
    return r
  } finally {
    await context.close().catch(() => {})
  }
}

// ---- aggregation and output ------------------------------------------------------------------

function leafAgg(objs, f) {
  const o = objs.find(x => x != null)
  if (o === undefined) return null
  if (typeof o !== 'object') {
    const vals = objs.filter(v => typeof v === 'number')
    return vals.length ? f(vals) : null
  }
  return Object.fromEntries(Object.keys(o).map(k => [k, leafAgg(objs.map(x => x?.[k]), f)]))
}
const median = xs => { const s = [...xs].sort((a, b) => a - b); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2 }

/** `pagesByScenario` maps each scenario to the page ids it recorded; `pair` uses one synthetic id. */
function aggregate(runs, pagesByScenario) {
  const out = {}
  for (const [s, ps] of Object.entries(pagesByScenario)) {
    for (const p of ps) {
      out[p.id] ??= {}
      const rs = runs.filter(r => r.page === p.id && r.scenario === s)
      const ok = rs.filter(r => !r.error && r.summary)
      const a = {
        runs: rs.length, ok: ok.length, errors: rs.filter(r => r.error).map(r => r.error),
        nav_mode: s === 'nav' || s === 'pair' ? [...new Set(rs.map(r => r.nav_mode).filter(Boolean))].join(',') || null : undefined,
        median: leafAgg(ok.map(r => r.summary), median),
        min: leafAgg(ok.map(r => r.summary), xs => Math.min(...xs)),
        max: leafAgg(ok.map(r => r.summary), xs => Math.max(...xs)),
        // a URL is not a number, so leafAgg cannot carry it: take it from the slowest run instead
        data_max_url: {},
      }
      for (const w of ['to_ready', 'to_idle']) {
        let best = null
        for (const r of ok) {
          const x = r.summary[w]
          if (x?.data_max_ms != null && (!best || x.data_max_ms > best.data_max_ms)) best = x
        }
        a.data_max_url[w] = best?.data_max_url ?? null
        for (const k of ['median', 'min', 'max']) delete a[k]?.[w]?.data_max_url
      }
      out[p.id][s] = a
    }
  }
  return out
}

const fmtInt = x => (x == null ? '–' : Math.round(x).toLocaleString('en-US'))
const fmtKB = x => (x == null ? '–' : (x / 1000).toLocaleString('en-US', { maximumFractionDigits: 1, minimumFractionDigits: 1 }))
function cell(agg, path, fmt) {
  const get = o => path.split('.').reduce((v, k) => v?.[k], o)
  const md = get(agg.median), lo = get(agg.min), hi = get(agg.max)
  if (md == null) return '–'
  return lo != null && hi != null && lo !== hi ? `${fmt(md)} (${fmt(lo)}–${fmt(hi)})` : fmt(md)
}

function markdown(meta, results, pagesByScenario, scenarios) {
  const L = []
  L.push(`# Gene page cost: ${meta.label} on ${meta.target}`, '')
  L.push(`- Recorded ${meta.time} by \`ui/bench/locus-bench.mjs\`, ${meta.runs} run(s) per page and scenario; cells are median (min–max) when runs differ.`)
  L.push(`- Target ${meta.target_origin}, data from ${meta.data_base}.`)
  L.push(`- App bundle \`${meta.app_bundle}\`; data manifest built ${meta.data_manifest.built}, pipeline commit \`${meta.data_manifest.pipeline_commit}\`; repo HEAD \`${meta.repo_head}\`${meta.repo_dirty ? ' (working tree has uncommitted changes)' : ''}.`)
  L.push(`- Playwright ${meta.playwright}, ${meta.browser}; headless, ${VIEWPORT.width}x${VIEWPORT.height}, no throttling.`)
  L.push(`- Machine ${meta.hostname}, ${meta.cpus} CPUs; network ${meta.network.interface ?? '?'}${meta.network.wireless ? ' (wireless)' : ''}${meta.network.note ? `: ${meta.network.note}` : ''}.`)
  L.push(`- Timings (ms) were taken on a shared, loaded machine: the 1-minute load average ranged ${meta.loadavg_1m.min.toFixed(1)}–${meta.loadavg_1m.max.toFixed(1)} across scenario starts and ends. Request counts and bytes are the primary baseline; treat ms as indicative.`)
  const custom = Object.entries(meta.ready_selectors ?? {}).filter(([, sel]) => sel !== meta.ready_selector)
  if (custom.length) L.push(`- Ready selectors other than the locus scatter: ${custom.map(([id, sel]) => `${id} \`${sel}\``).join('; ')}.`)
  L.push('- Windows: `ready` counts requests started before the locus scatter was drawn (or the page\'s own ready selector matched); `idle` counts everything until no request was in flight or started for 2 s (capped 30 s after ready). Bytes are response body plus headers as transferred. `nav` is measured from just before the in-app navigation, with the engine and search index already loaded on /genes. `pair` is the same navigation, but from a gene page loaded cold first, so the second page re-reads pack URLs the browser has already cached.', '')
  const cols = [
    ['data req', 'requests_data', fmtInt], ['all req', 'requests_all', fmtInt],
    ['data GET', 'data_get', fmtInt], ['data HEAD', 'data_head', fmtInt], ['data 206', 'data_206', fmtInt],
    ['data 200', 'data_200', fmtInt], ['data cached', 'data_from_cache', fmtInt],
    ['data KB', 'bytes_data', fmtKB], ['all KB', 'bytes_all', fmtKB],
    ['app KB', 'bytes_by_category.app', fmtKB], ['duckdb-cdn KB', 'bytes_by_category.duckdb-cdn', fmtKB],
    ['duckdb-ext KB', 'bytes_by_category.duckdb-ext', fmtKB], ['other KB', 'bytes_by_category.other', fmtKB],
    ['max data concurrency', 'max_concurrency_data', fmtInt],
    ['longest data req ms', 'data_max_ms', fmtInt],
  ]
  for (const s of scenarios) {
    const pages = pagesByScenario[s]
    L.push(`## ${s}`, '')
    L.push(`| page | window | ${cols.map(c => c[0]).join(' | ')} | ms to ready | ms to idle | ms first→last data before ready |`)
    L.push(`|---|---|${cols.map(() => '---:').join('|')}|---:|---:|---:|`)
    for (const p of pages) {
      const a = results[p.id][s]
      if (!a.ok) { L.push(`| ${p.id} | – | ${cols.map(() => '–').join(' | ')} | failed: ${a.errors.join('; ')} | | |`); continue }
      for (const w of ['to_ready', 'to_idle']) {
        const timeCols = w === 'to_ready'
          ? [cell(a, 'ms_ready', fmtInt), '', cell(a, 'ms_data_first_to_last_before_ready', fmtInt)]
          : ['', cell(a, 'ms_idle', fmtInt), '']
        L.push(`| ${w === 'to_ready' ? p.id : ''} | ${w === 'to_ready' ? 'ready' : 'idle'} | ${cols.map(c => cell(a, `${w}.${c[1]}`, c[2])).join(' | ')} | ${timeCols.join(' | ')} |`)
      }
    }
    L.push('')
    // data requests by kind (whole scenario), and the gene page marks
    const kcols = [
      ...KIND_NAMES.map(k => [`${k} req`, `to_idle.kinds.${k}.requests`, fmtInt]),
      ['pack+variants KB', null, null], ['block/variants overlap', 'to_idle.pack_overlap', fmtInt],
      ['block/variants/gwas overlap', 'to_idle.pack_overlap_gwas', fmtInt], ['ms block end→intron start', 'to_idle.intron_after_block_ms', fmtInt],
      ['ms nav→detail', 'marks.ms_detail', fmtInt], ['ms nav→drawn', 'marks.ms_drawn', fmtInt],
      ['ms hit→detail', 'marks.ms_hit_to_detail', fmtInt], ['ms hit→drawn', 'marks.ms_hit_to_drawn', fmtInt],
      ['block fetch ms', 'marks.pack_block_ms', fmtInt], ['variants fetch ms', 'marks.pack_variants_ms', fmtInt],
      ['gwas fetch ms', 'marks.pack_gwas_ms', fmtInt], ['intron fetch ms', 'marks.pack_intron_ms', fmtInt],
      ['decode ms', 'marks.pack_decode_ms', fmtInt], ['worker wait ms', 'marks.pack_worker_wait_ms', fmtInt],
      ['insert ms', 'marks.pack_insert_ms', fmtInt],
      ['trans fetch ms', 'marks.pack_trans_ms', fmtInt], ['trans decode ms', 'marks.pack_trans_decode_ms', fmtInt],
      ['trans worker wait ms', 'marks.pack_trans_wait_ms', fmtInt], ['trans insert ms', 'marks.pack_trans_insert_ms', fmtInt],
    ]
    // the slowest data request's file name, next to the ms the first table reports
    const slowest = a => { const u = a.data_max_url?.to_idle; return u ? `\`${u.split('/').pop().split('?')[0]}\`` : '–' }
    L.push(`| page | ${kcols.map(c => c[0]).join(' | ')} | slowest data request |`, `|---|${kcols.map(() => '---:').join('|')}|---|`)
    for (const p of pages) {
      const a = results[p.id][s]
      if (!a.ok) continue
      const kb = o => (o ? (o.to_idle.kinds.eqtl_pack.bytes + o.to_idle.kinds.variants.bytes) : null)
      L.push(`| ${p.id} | ${kcols.map(c => (c[1] ? cell(a, c[1], c[2]) : fmtKB(kb(a.median)))).join(' | ')} | ${slowest(a)} |`)
    }
    L.push('')
  }
  return L.join('\n')
}

function networkInfo(note) {
  const info = { interface: null, wireless: null, note }
  try {
    const m = /\bdev (\S+)/.exec(execSync('ip route get 1.1.1.1', { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }))
    if (m) { info.interface = m[1]; info.wireless = existsSync(`/sys/class/net/${m[1]}/wireless`) }
  } catch { /* not linux, or no route */ }
  return info
}

function git(cmd) {
  try { return execSync(`git ${cmd}`, { cwd: REPO, encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }).trim() } catch { return null }
}

// ---- main ------------------------------------------------------------------------------------

async function main() {
  const args = parseArgs(process.argv.slice(2))
  const target = TARGETS[args.target]
  const allPages = JSON.parse(readFileSync(join(HERE, 'pages.json'), 'utf8'))
  const ids = args.pages ? args.pages.split(',').map(s => s.trim()) : allPages.map(p => p.id)
  for (const id of ids) if (!allPages.some(p => p.id === id)) die(`unknown page id ${id}`)
  OPTS.query = args.query
  OPTS.latency = args.latency
  // --pages keeps its order, so `pair` navigates from the first id to the second
  const byId = new Map(allPages.map(p => [p.id, p]))
  const pages = ids.map(id => ({ ...byId.get(id), path: withQuery(byId.get(id).path, args.query) }))
  const isPair = args.scenarios.includes('pair')
  const pairPage = isPair ? { id: `${pages[0].id}>${pages[1].id}` } : null
  const pagesByScenario = Object.fromEntries(args.scenarios.map(s => [s, s === 'pair' ? [pairPage] : pages]))

  // what the target serves; rehearsal reads the staged copy rather than manifest.json
  // qtlb v1 serves store.json where v0 served manifest.json (--manifest still names either)
  const manifestKey = args.manifest ?? 'store.json'
  const manifestPath = new URL(`${target.dataBase}/${manifestKey}`).pathname
  const manRes = await fetch(`${target.dataBase}/${manifestKey}`).catch(e => ({ ok: false, status: String(e.cause?.code ?? e.message) }))
  if (!manRes.ok || manRes.status !== 200) {
    if (target.name === 'preview') die(`${target.dataBase}/${manifestKey} did not answer 200 (${manRes.status}). Start the preview first: ${PREVIEW_HELP}`)
    die(`${target.dataBase}/${manifestKey}: ${manRes.status}`)
  }
  const manifest = await manRes.json()
  const indexHtml = await (await fetch(`${target.origin}/`)).text()
  const bundle = /assets\/index-[^"'\s]+\.js/.exec(indexHtml)?.[0] ?? null

  const browser = await chromium.launch({ headless: true })
  const iso = new Date().toISOString()
  const meta = {
    time: iso, label: args.label, target: target.name, target_origin: target.origin, data_base: target.dataBase,
    runs: args.runs, scenarios: args.scenarios, pages: pages.map(p => p.id), nav_mode_requested: args.navMode, manifest_key: manifestKey,
    playwright: require('playwright/package.json').version, browser: `${browser.browserType().name()} ${browser.version()}`,
    app_bundle: bundle,
    data_manifest: { built: manifest.built ?? null, pipeline_commit: manifest.pipeline_commit ?? null },
    repo_head: git('rev-parse HEAD'), repo_dirty: !!git('status --porcelain'),
    hostname: os.hostname(), cpus: os.cpus().length, cpu_model: os.cpus()[0]?.model ?? null,
    network: networkInfo(args.note),
    ready_selector: READY_SELECTOR, ready_selectors: Object.fromEntries(pages.map(p => [p.id, readySelector(p)])), idle_quiet_ms: IDLE_QUIET, idle_cap_ms: IDLE_CAP,
    query: args.query, latency_ms: args.latency,
  }
  log(`${target.name}: bundle ${bundle}, data built ${meta.data_manifest.built}, ${meta.browser}, playwright ${meta.playwright}`)

  const runs = []
  let selfChecked = false
  const selfCheck = r => {
    const ok = !r.error && r.requests.some(q => q.category === 'data' && q.range && q.status === 206)
    if (!ok) {
      die(`worker requests not captured: a cold FLNC load recorded no data request with a Range header and status 206${r.error ? ` (${r.error})` : ''}. Switch the harness to puppeteer-core with per-worker CDP Network events (see bench/README.md).`, 2)
    }
    const xhr = r.requests.filter(q => q.category === 'data' && q.status === 206 && q.resource_type === 'xhr').length
    log(`self-check passed: ${r.requests.filter(q => q.category === 'data' && q.status === 206).length} ranged data responses (206), ${xhr} of them DuckDB worker XHRs`)
    if (args.latency) {
      // request sent to response headers, from Playwright's timing: at least the added latency when emulation applies
      const wait = q => (q.timing && q.timing.responseStart > 0 && q.timing.requestStart >= 0 ? q.timing.responseStart - q.timing.requestStart : null)
      const xhrs = r.requests.filter(q => q.category === 'data' && q.resource_type === 'xhr' && wait(q) != null)
      const fetches = r.requests.filter(q => q.category === 'data' && q.resource_type !== 'xhr' && wait(q) != null)
      const delayed = qs => qs.filter(q => wait(q) >= 0.9 * args.latency).length
      log(`latency ${args.latency} ms: ${delayed(xhrs)} of ${xhrs.length} worker XHRs and ${delayed(fetches)} of ${fetches.length} page requests waited at least 90% of it for headers`)
      if (xhrs.length && delayed(xhrs) < xhrs.length) die('latency emulation does not reach the DuckDB worker requests; the parquet path would be measured without it', 2)
    }
    selfChecked = true
  }
  const isLocal = q => ['localhost', '127.0.0.1'].includes(new URL(q.url).hostname)
  const guardTarget = r => {
    if (target.name === 'preview') {
      const remote = r.requests.filter(q => q.category === 'data' && !isLocal(q))
      if (remote.length) die(`preview run read data from ${new URL(remote[0].url).host}, not the local server; build the bundle for local data: ${PREVIEW_HELP}`)
    }
    if (target.name === 'rehearsal') {
      // the mirror image of the preview guard: a rehearsal that reads local files proves nothing
      const local = r.requests.filter(q => q.category === 'data' && isLocal(q))
      if (local.length) die(`rehearsal run read data from ${new URL(local[0].url).host}, not the bucket; build the bundle for ${target.dataBase} (npm run build, no VITE_DATA_BASE override)`)
      const wrong = r.requests.filter(q => q.category === 'data' && MANIFEST_RE.test(q.url) && new URL(q.url).pathname !== manifestPath)
      if (wrong.length) die(`rehearsal run fetched ${new URL(wrong[0].url).pathname}, not ${manifestKey}; rebuild with VITE_MANIFEST=${manifestKey}`)
    }
    if (target.name === 'live') return
    const missing = r.requests.filter(q => q.category === 'data' && q.status === 404)
    if (missing.length) die(`${target.name} data responses 404: ${[...new Set(missing.map(q => new URL(q.url).pathname))].join(', ')}`)
  }

  const flncPage = allPages.find(p => p.id === args.selfCheckPage)
  if (!flncPage) die(`unknown --self-check-page ${args.selfCheckPage}`)
  const flnc = { ...flncPage, path: withQuery(flncPage.path, args.query) }
  if (pages[0]?.id !== args.selfCheckPage || !(args.scenarios.includes('cold') || args.scenarios.includes('warm'))) {
    const [r] = await coldWarm(browser, target, flnc, false)
    selfCheck(r)
    guardTarget(r)
  }

  try {
    for (let run = 1; run <= args.runs; run++) {
      if (isPair) {
        const r = await pair(browser, target, pages[0], pages[1], args.navMode)
        guardTarget(r)
        runs.push({ run, ...r })
        continue
      }
      for (const p of pages) {
        if (args.scenarios.includes('cold') || args.scenarios.includes('warm')) {
          const rs = await coldWarm(browser, target, p, args.scenarios.includes('warm'))
          for (const r of rs) {
            if (!selfChecked && p.id === args.selfCheckPage && r.scenario === 'cold') selfCheck(r)
            guardTarget(r)
            if (args.scenarios.includes(r.scenario)) runs.push({ run, ...r })
          }
        }
        if (args.scenarios.includes('nav')) {
          const r = await nav(browser, target, p, args.navMode)
          guardTarget(r)
          runs.push({ run, ...r })
        }
      }
    }
  } finally {
    await browser.close()
  }

  const loads = runs.flatMap(r => [r.loadavg_start[0], r.loadavg_end[0]])
  meta.loadavg_1m = { min: Math.min(...loads), max: Math.max(...loads) }
  const results = aggregate(runs, pagesByScenario)

  const stamp = iso.replace(/[:.]/g, '-')
  const outDir = join(HERE, 'out')
  const resDir = join(HERE, 'results')
  mkdirSync(outDir, { recursive: true })
  mkdirSync(resDir, { recursive: true })
  const rawPath = join(outDir, `${args.label}-${target.name}-${stamp}.json`)
  writeFileSync(rawPath, JSON.stringify({ meta, runs }, null, 1))
  const base = join(resDir, `${args.label}-${target.name}`)
  writeFileSync(`${base}.json`, JSON.stringify({ meta, raw_log: relative(UI, rawPath), results }, null, 2) + '\n')
  writeFileSync(`${base}.md`, markdown(meta, results, pagesByScenario, args.scenarios) + '\n')
  log(`wrote ${relative(UI, rawPath)}, ${relative(UI, base)}.json, ${relative(UI, base)}.md`)
  const failed = runs.filter(r => r.error)
  if (failed.length) die(`${failed.length} scenario run(s) failed: ${failed.map(r => `${r.page} ${r.scenario} run ${r.run}: ${r.error}`).join('; ')}`)
}

main().catch(e => die(e?.stack ?? String(e)))
