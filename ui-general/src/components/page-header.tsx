import type { ReactNode } from 'react'
import { Link } from 'react-router'
import { ChevronLeft, ChevronRight } from 'lucide-react'

/** Title row (+ optional muted count/meta and right-aligned actions), description below. */
export function PageHeader({ title, meta, description, actions, back, crumbs, size = 'lg', className = 'mb-8 mt-2' }: {
  title: ReactNode
  meta?: ReactNode
  description?: ReactNode
  actions?: ReactNode
  back?: { to: string; label: string }
  /** Trail above the title: linked ancestors, then the current page's label. */
  crumbs?: { to?: string; label: string }[]
  size?: 'lg' | 'md'
  /** Outer margins; replaces the default rather than adding to it. */
  className?: string
}) {
  return (
    <div className={className}>
      {crumbs && (
        <nav className="mb-4 flex min-w-0 items-center gap-1.5 text-sm">
          {crumbs.map((c, i) => (
            <span key={i} className="flex min-w-0 items-center gap-1.5">
              {i > 0 && <ChevronRight className="size-3.5 shrink-0 text-base-content/30" />}
              {c.to ? <Link to={c.to} className="truncate text-base-content/60 transition-colors hover:text-base-content">{c.label}</Link>
                : <span className="truncate font-medium">{c.label}</span>}
            </span>
          ))}
        </nav>
      )}
      {back && (
        <Link to={back.to} className="mb-5 inline-flex items-center gap-1 text-sm text-base-content/50 transition-colors hover:text-base-content">
          <ChevronLeft className="size-4" /> {back.label}
        </Link>
      )}
      <div className="flex items-end justify-between gap-4">
        <div className="flex min-w-0 flex-wrap items-baseline gap-x-2.5 gap-y-1">
          <h1 className={`font-light tracking-tight ${size === 'lg' ? 'text-[1.6875rem]' : 'text-xl'}`}>{title}</h1>
          {meta && <span className="text-sm tabular-nums text-base-content/40">{meta}</span>}
        </div>
        {actions && <div className="flex shrink-0 items-center gap-2">{actions}</div>}
      </div>
      {description && <div className="text-sm text-base-content/55">{description}</div>}
    </div>
  )
}
