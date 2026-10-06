import { ArrowDown, ArrowUp, ChevronsUpDown } from 'lucide-react'
import type { ReactNode } from 'react'

export type SortState = { by: string; order: 'asc' | 'desc' }

/** Three-step cycle: default direction → reversed → cleared. */
export function toggleSort(prev: SortState, key: string, def: 'asc' | 'desc' = 'desc'): SortState {
  if (prev.by !== key) return { by: key, order: def }
  if (prev.order === def) return { by: key, order: def === 'asc' ? 'desc' : 'asc' }
  return { by: '', order: 'asc' }
}

export function SortableTh({ sortKey, label, sort, onSort, defaultOrder = 'desc', className = '', align = 'left' }: {
  sortKey: string; label: ReactNode; sort: SortState; onSort: (s: SortState) => void
  defaultOrder?: 'asc' | 'desc'; className?: string; align?: 'left' | 'right' | 'center'
}) {
  const active = sort.by === sortKey
  const Icon = !active ? ChevronsUpDown : sort.order === 'asc' ? ArrowUp : ArrowDown
  const justify = align === 'right' ? 'justify-end' : align === 'center' ? 'justify-center' : 'justify-start'
  return (
    <th className={className}>
      <button type="button" className={`group flex w-full items-center gap-1 ${justify} ${active ? 'text-base-content' : ''}`}
        onClick={() => onSort(toggleSort(sort, sortKey, defaultOrder))}>
        {label}
        <Icon className={`size-3 shrink-0 ${active ? 'opacity-60' : 'opacity-0 group-hover:opacity-30'}`} />
      </button>
    </th>
  )
}
