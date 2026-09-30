import { CompareSkeleton, LocusSkeleton } from '@/components/plot-skeleton'

export function Loading({ label = 'Loading…' }: { label?: string }) {
  return <div className="flex items-center gap-2 p-4 text-sm text-base-content/60"><span className="loading loading-spinner loading-sm" /> {label}</div>
}

/** Empty state in the same bordered box a table would occupy, so the text is not floating. */
export function Empty({ label }: { label: string }) {
  return <div className="rounded-lg border border-base-300 px-4 py-6 text-center text-sm text-base-content/60">{label}</div>
}

/** What a page says where a section's data has no v1 object yet (trans, GWAS, coloc). */
export const NOT_YET = 'Not available in this data release yet.'

/** The plain empty state for a section whose data this release does not carry. */
export function Unavailable({ what }: { what: string }) {
  return <Empty label={`${what}: not available in this data release yet.`} />
}

export function ErrorState({ message }: { message: string }) {
  return <div className="alert alert-error text-sm">{message}</div>
}

export function Bar({ className = '' }: { className?: string }) {
  return <div className={`animate-pulse rounded bg-base-300/60 ${className}`} aria-hidden />
}

export type SkelCol = { w?: string; align?: 'right' }

/** Table placeholder matching the `table table-sm` it precedes. */
export function TableSkeleton({ columns, rows = 8 }: { columns: SkelCol[]; rows?: number }) {
  return (
    <div className="overflow-hidden rounded-lg border border-base-300" aria-busy="true">
      <table className="table table-sm">
        <thead><tr>{columns.map((c, i) => <th key={i}><Bar className={`h-3 w-12 ${c.align === 'right' ? 'ml-auto' : ''}`} /></th>)}</tr></thead>
        <tbody>
          {Array.from({ length: rows }).map((_, r) => (
            <tr key={r}>{columns.map((c, i) => <td key={i}><Bar className={`h-4 ${c.w ?? 'w-24'} ${c.align === 'right' ? 'ml-auto' : ''}`} /></td>)}</tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

/** Key/value table placeholder matching KvTable: bordered box, zebra rows, label + value bars. */
export function KvSkeleton({ rows = 6 }: { rows?: number }) {
  return (
    <div className="overflow-hidden rounded-lg border border-base-300" aria-busy="true">
      {Array.from({ length: rows }).map((_, i) => (
        <div key={i} className={`flex items-center gap-4 px-3 py-2 ${i % 2 === 1 ? 'bg-base-200' : ''}`}>
          <Bar className="h-3 w-28 shrink-0" /><Bar className={`h-3 ${i % 3 === 0 ? 'w-40' : i % 3 === 1 ? 'w-24' : 'w-32'}`} />
        </div>
      ))}
    </div>
  )
}

/** The body of a gene-page tab: two kv tables side by side and, when `plot`, the locus row
 *  (wide scatter + square LocusCompare) at the plots' real height. */
export function TabSkeleton({ plot = false, kvRows = 9, chr }: { plot?: boolean; kvRows?: number; chr?: string }) {
  return (
    <div className="space-y-8" aria-busy="true">
      <div className="grid items-start gap-4 md:grid-cols-2"><KvSkeleton rows={kvRows} /><KvSkeleton rows={kvRows} /></div>
      {plot && (
        <div className="space-y-3">
          <h2 className="text-sm font-medium text-base-content/65">Locus</h2>
          <div className="flex flex-col gap-6 md:flex-row md:gap-2">
            <div className="min-w-0 flex-1"><LocusSkeleton chr={chr} /></div>
            <div className="size-[290px] shrink-0"><CompareSkeleton /></div>
          </div>
        </div>
      )}
    </div>
  )
}

/** A detail page before its identity resolves: breadcrumb, title + id, a chip row, then the
 *  tab body (with the locus row on gene pages). */
export function DetailSkeleton({ plot = false, kvRows = 9 }: { plot?: boolean; kvRows?: number }) {
  return (
    <div aria-busy="true">
      <div className="mb-8 mt-2 space-y-3">
        <Bar className="h-3 w-24" />
        <div className="flex items-end justify-between gap-4">
          <div className="flex items-baseline gap-2.5"><Bar className="h-7 w-32" /><Bar className="h-3 w-36" /></div>
          <Bar className="h-8 w-52 rounded-[0.625rem]" />
        </div>
        <div className="flex gap-1.5"><Bar className="h-5 w-14 rounded-full" /><Bar className="h-5 w-24 rounded-full" /></div>
      </div>
      <TabSkeleton plot={plot} kvRows={kvRows} />
    </div>
  )
}
