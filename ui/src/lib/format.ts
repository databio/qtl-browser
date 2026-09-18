export function fmtP(p: unknown): string {
  if (p == null || Number.isNaN(p)) return ''
  const x = Number(p)
  if (x === 0) return '0'
  if (x >= 0.001) return x.toPrecision(2)
  return x.toExponential(1)
}

export function fmtNum(x: unknown, digits = 3): string {
  if (x == null || Number.isNaN(x)) return ''
  return Number(x).toFixed(digits)
}

export function fmtInt(x: unknown): string {
  if (x == null) return ''
  return Number(x).toLocaleString('en-US')
}

export function fmtSlopeSE(slope: unknown, se: unknown): string {
  if (slope == null) return ''
  return `${fmtNum(slope)} ± ${fmtNum(se)}`
}

export function rsFromNumber(n: unknown): string {
  return n == null ? '' : `rs${n}`
}

export function neglog10(p: unknown): number {
  const x = Number(p)
  return x > 0 ? -Math.log10(x) : 320
}

export function fmtBp(x: unknown): string {
  const n = Number(x)
  if (Math.abs(n) >= 1e6) return `${(n / 1e6).toFixed(2)} Mb`
  if (Math.abs(n) >= 1e3) return `${(n / 1e3).toFixed(1)} kb`
  return `${n} bp`
}

/** leafcutter phenotype id "chr7:128849578:128849976:clu_77648_+:ENSG…" → "chr7:128,849,578–128,849,976 (+)".
 *  The cluster label is only meaningful next to sibling introns, so it stays in the phenotype list. */
export function fmtPhenotype(id: string): string {
  const m = /^(chr[^:]+):(\d+):(\d+):(clu_\d+)_([+-?])(?::|$)/.exec(id)
  if (!m) return id
  return `${m[1]}:${Number(m[2]).toLocaleString('en-US')}–${Number(m[3]).toLocaleString('en-US')} (${m[5]})`
}

/** Byte counts for the cis-scan button: the span a click will download. */
export function fmtBytes(n: number): string {
  if (n >= 1e6) return `${(n / 1e6).toFixed(1)} MB`
  if (n >= 1e3) return `${Math.round(n / 1e3)} KB`
  return `${n} B`
}
