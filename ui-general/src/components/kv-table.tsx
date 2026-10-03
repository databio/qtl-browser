import type { ReactNode } from 'react'

export type KvRow = { label: string; value: ReactNode }

/**
 * Two-column key/value metadata table: bordered, zebra-striped.
 *
 * `labelWidth` is the label column's preferred width (auto table layout, so it is a preference, not
 * a cap). The default suits the page's usual short labels; a table with longer labels and short
 * numeric values wants more, because side by side in a grid there is no slack for the column to
 * grow into and every label wraps instead.
 */
export function KvTable({ title, rows, align = 'left', labelWidth = 'w-44' }: {
  title?: ReactNode; rows: KvRow[]; align?: 'left' | 'right'; labelWidth?: string
}) {
  if (!rows.length) return null
  return (
    <div className="space-y-2">
      {title && <h3 className="text-sm font-semibold uppercase tracking-wide text-base-content/65">{title}</h3>}
      <div className="overflow-x-auto rounded-lg border border-base-300 bg-base-100">
        <table className="table table-sm w-full text-xs">
          <tbody>
            {rows.map(({ label, value }, i) => (
              <tr key={label} className={i % 2 === 1 ? 'bg-base-200' : ''}>
                <td className={`${labelWidth} align-top font-medium text-base-content/60`}>{label}</td>
                <td className={`text-base-content/70 ${align === 'right' ? 'text-right tabular-nums' : ''}`}>{value}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  )
}
