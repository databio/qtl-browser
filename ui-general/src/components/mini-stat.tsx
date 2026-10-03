import type { ReactNode } from 'react'

/** Compact stat tile: tinted fill, no border. */
export function MiniStat({ label, value, hint }: { label: string; value: ReactNode; hint?: ReactNode }) {
  return (
    <div className="min-w-0 rounded-lg bg-base-200/50 px-3 py-2">
      <div className="text-xs text-base-content/55">{label}</div>
      <div className="mt-0.5 truncate text-sm font-semibold tabular-nums">{value ?? '—'}</div>
      {hint && <div className="truncate text-xs text-base-content/55">{hint}</div>}
    </div>
  )
}
