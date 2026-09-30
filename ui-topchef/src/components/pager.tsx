export const PAGE_SIZES = [10, 25, 50, 100]

export function Pager({ total, offset, pageSize, onPage, onPageSize, pageSizes = PAGE_SIZES }: {
  total: number; offset: number; pageSize: number; onPage: (o: number) => void; onPageSize: (n: number) => void; pageSizes?: number[]
}) {
  const from = total === 0 ? 0 : offset + 1
  const to = Math.min(offset + pageSize, total)
  return (
    <div className="mt-2 flex shrink-0 items-center justify-between gap-3 text-sm text-base-content/60">
      <div className="flex items-center gap-2">
        <span className="whitespace-nowrap tabular-nums">{from.toLocaleString()}–{to.toLocaleString()} of {total.toLocaleString()}</span>
        <select className="select select-bordered select-sm h-8 rounded-lg" value={pageSize} onChange={e => onPageSize(Number(e.target.value))} title="Rows per page">
          {pageSizes.map(n => <option key={n} value={n}>{n} / page</option>)}
        </select>
      </div>
      <div className="flex items-center gap-1.5">
        <button className="btn btn-sm h-8 rounded-lg border-base-300 font-medium" disabled={offset === 0} onClick={() => onPage(Math.max(0, offset - pageSize))}>Prev</button>
        <button className="btn btn-sm h-8 rounded-lg border-base-300 font-medium" disabled={to >= total} onClick={() => onPage(offset + pageSize)}>Next</button>
      </div>
    </div>
  )
}
