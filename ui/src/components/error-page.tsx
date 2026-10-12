import { Component, useEffect, type ErrorInfo, type ReactNode } from 'react'
import { RotateCw } from 'lucide-react'
import { suppressToasts } from '@/lib/toast'

/**
 * The one page-level error display: for a failure that takes the page down (the store cannot
 * open, the query engine cannot start, a page's own data cannot load, a component throws). It
 * sits in the content area under the navbar, so the site stays navigable. A failure that leaves
 * the page usable is a toast instead (lib/toast.ts).
 */
export function ErrorPage({ title, error, action }: {
  title: string
  /** the technical message, shown in a muted monospace box for the console-literate */
  error?: unknown
  /** default: reload the page */
  action?: { label: string; onClick: () => void }
}) {
  const message = error instanceof Error ? error.message : error == null ? null : String(error)
  const act = action ?? { label: 'Reload', onClick: () => window.location.reload() }
  // the sections that failed alongside have nothing to add: no toasts while this page is up
  useEffect(() => suppressToasts(), [])
  return (
    <div className="flex flex-col items-center justify-center gap-4 px-4 py-24">
      <h1 className="text-xl font-light text-base-content">{title}</h1>
      {message && (
        <div className="w-full max-w-lg rounded-lg border border-base-300 bg-base-200/50 px-4 py-3">
          <p className="break-words font-mono text-xs text-base-content/50">{message}</p>
        </div>
      )}
      <button onClick={act.onClick} className="btn btn-sm btn-primary gap-1.5 rounded-full">
        <RotateCw size={14} />
        {act.label}
      </button>
    </div>
  )
}

/** Last-resort catch for a throw during render: the error page instead of React unmounting the
 *  whole tree to a blank screen. Wraps the routes, inside the navbar. */
export class ErrorBoundary extends Component<{ children: ReactNode }, { error: Error | null }> {
  state = { error: null as Error | null }

  static getDerivedStateFromError(error: Error) {
    return { error }
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    console.error('render error:', error, info.componentStack)
  }

  render() {
    if (this.state.error) return <ErrorPage title="Something went wrong." error={this.state.error} />
    return this.props.children
  }
}
