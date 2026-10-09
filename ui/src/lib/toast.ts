import { toast } from 'sonner'

/**
 * The one section-level error display: a persistent toast (the `Toaster` in App.tsx) for a failure
 * that leaves the page usable, such as a trans table, the GWAS panel or an export. A failure that
 * takes the page down is an ErrorPage instead (components/error-page.tsx), and while one is shown
 * no toast appears: the sections that failed for the same reason would only repeat it.
 *
 * One root cause, one toast: the toast's id is the error message, so the locus plot, the compare
 * panel and the trans table all failing on the same missing object collapse into the first one's
 * title with that message under it, and a retry updates it rather than stacking.
 */
let suppressed = 0

export function toastError(title: string, error: unknown): void {
  console.error(title, error)
  if (suppressed) return
  const message = error instanceof Error ? error.message : String(error)
  toast.error(title, { id: message || title, description: message, duration: Infinity })
}

/** While a page-level error is shown: clear the toasts and hold new ones. Returns the release. */
export function suppressToasts(): () => void {
  suppressed++
  toast.dismiss()
  return () => { suppressed-- }
}
