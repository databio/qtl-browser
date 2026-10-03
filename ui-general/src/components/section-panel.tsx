import type { ReactNode } from 'react'

/** Titled content section: muted small heading, optional description, optional action. */
export function SectionPanel({ title, description, action, children, className = '', gap = 'space-y-3' }: {
  title?: ReactNode; description?: ReactNode; action?: ReactNode; children: ReactNode; className?: string; gap?: string
}) {
  return (
    <section className={`${gap} ${className}`}>
      {(title || action) && (
        <div className="flex items-end justify-between gap-3">
          <div className="min-w-0">
            {title && <h2 className="text-sm font-medium text-base-content/65">{title}</h2>}
            {description && <p className="text-xs text-base-content/55">{description}</p>}
          </div>
          {action}
        </div>
      )}
      {children}
    </section>
  )
}
