import type { Row } from './db'

/** SPEC.md section 9: an export of the packs' rounded per-variant values says so in its own file
 *  name, so a downloaded file still says what it holds. */
export const roundedCsvName = (stem: string) => `${stem}.rounded.csv`

export function downloadCSV(name: string, data: Row[], columns?: string[]) {
  if (!data.length) return
  const cols = columns ?? Object.keys(data[0])
  const esc = (v: unknown) => {
    if (v == null) return ''
    const s = String(v)
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
  }
  const lines = [cols.join(','), ...data.map(r => cols.map(c => esc(r[c])).join(','))]
  const blob = new Blob([lines.join('\n')], { type: 'text/csv' })
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url; a.download = name; a.click()
  URL.revokeObjectURL(url)
}
