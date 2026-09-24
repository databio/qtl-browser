#!/usr/bin/env node
/**
 * Compare two locus-bench summaries: one markdown table per scenario, each metric as
 * `before -> after (change %)` for the page ids present in both (medians).
 *
 *   node bench/compare.mjs results/baseline-live.json results/<other>.json
 */
import { readFileSync } from 'node:fs'
import { resolve } from 'node:path'

const [a, b] = process.argv.slice(2)
if (!a || !b) {
  console.error('usage: node bench/compare.mjs <before.json> <after.json>')
  process.exit(1)
}
const load = p => JSON.parse(readFileSync(resolve(process.cwd(), p), 'utf8'))
const before = load(a)
const after = load(b)

const kb = x => `${(x / 1000).toLocaleString('en-US', { maximumFractionDigits: 1 })} KB`
const int = x => Math.round(x).toLocaleString('en-US')
const METRICS = [
  ['data requests to ready', 'to_ready.requests_data', int],
  ['all requests to ready', 'to_ready.requests_all', int],
  ['data bytes to ready', 'to_ready.bytes_data', kb],
  ['all bytes to ready', 'to_ready.bytes_all', kb],
  ['max data concurrency to ready', 'to_ready.max_concurrency_data', int],
  ['data requests to idle', 'to_idle.requests_data', int],
  ['all requests to idle', 'to_idle.requests_all', int],
  ['data bytes to idle', 'to_idle.bytes_data', kb],
  ['all bytes to idle', 'to_idle.bytes_all', kb],
  ['max data concurrency to idle', 'to_idle.max_concurrency_data', int],
  // the longest single data request in each window: where a stall shows up
  ['longest data request to ready (ms)', 'to_ready.data_max_ms', int],
  ['longest data request to idle (ms)', 'to_idle.data_max_ms', int],
  ['ms to ready', 'ms_ready', int],
  ['ms to idle', 'ms_idle', int],
  ['ms first to last data request before ready', 'ms_data_first_to_last_before_ready', int],
  // the gene page's trans tables: parquet reads before Plan 8, one trans pack frame after
  ['trans parquet requests to idle', 'to_idle.kinds.trans.requests', int],
  ['trans parquet bytes to idle', 'to_idle.kinds.trans.bytes', kb],
  ['trans pack requests to idle', 'to_idle.kinds.trans_pack.requests', int],
  ['trans pack bytes to idle', 'to_idle.kinds.trans_pack.bytes', kb],
  ['trans fetch ms (pack:trans)', 'marks.pack_trans_ms', int],
  ['trans decode ms (pack:trans-decode)', 'marks.pack_trans_decode_ms', int],
  ['trans insert ms', 'marks.pack_trans_insert_ms', int],
]
const get = (o, path) => path.split('.').reduce((v, k) => v?.[k], o)

function change(x, y, fmt) {
  if (x == null && y == null) return '–'
  // a metric one of the two runs did not record
  if (x == null || y == null) return `${x == null ? '–' : fmt(x)} -> ${y == null ? '–' : fmt(y)}`
  const pct = x === 0 ? (y === 0 ? '0%' : 'new') : `${y >= x ? '+' : ''}${(((y - x) / x) * 100).toFixed(0)}%`
  return `${fmt(x)} -> ${fmt(y)} (${pct})`
}

const out = []
out.push(`Before: ${before.meta.label} on ${before.meta.target} (${before.meta.time}, bundle \`${before.meta.app_bundle}\`)`)
out.push(`After: ${after.meta.label} on ${after.meta.target} (${after.meta.time}, bundle \`${after.meta.app_bundle}\`)`, '')
const pages = Object.keys(before.results).filter(id => id in after.results)
const scenarios = ['cold', 'warm', 'nav', 'pair'].filter(s => pages.some(id => before.results[id][s] && after.results[id][s]))
for (const s of scenarios) {
  const ids = pages.filter(id => before.results[id][s]?.ok && after.results[id][s]?.ok)
  if (!ids.length) continue
  out.push(`### ${s}`, '')
  out.push(`| metric | ${ids.join(' | ')} |`)
  out.push(`|---|${ids.map(() => '---').join('|')}|`)
  for (const [label, path, fmt] of METRICS) {
    if (ids.every(id => get(before.results[id][s].median, path) == null && get(after.results[id][s].median, path) == null)) continue
    out.push(`| ${label} | ${ids.map(id => change(get(before.results[id][s].median, path), get(after.results[id][s].median, path), fmt)).join(' | ')} |`)
  }
  // the slowest request's URL is not a median, so it is listed rather than compared
  const url = (r, id) => r.results[id][s].data_max_url?.to_idle?.split('/').pop()?.split('?')[0] ?? '–'
  if (ids.some(id => url(before, id) !== '–' || url(after, id) !== '–')) {
    out.push(`| slowest data request to idle | ${ids.map(id => `${url(before, id)} -> ${url(after, id)}`).join(' | ')} |`)
  }
  out.push('')
}
console.log(out.join('\n'))
