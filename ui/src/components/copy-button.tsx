import { useEffect, useState } from 'react'
import { Check, Copy } from 'lucide-react'

/** Small icon button that copies `text` to the clipboard, then shows a check and the word
 *  "Copied" for a moment. */
export function CopyButton({ text, label = 'Copy', className = '' }: { text: string; label?: string; className?: string }) {
  const [done, setDone] = useState(false)
  useEffect(() => {
    if (!done) return
    const t = setTimeout(() => setDone(false), 1500)
    return () => clearTimeout(t)
  }, [done])
  return (
    <button type="button" title={done ? 'Copied' : `${label} ${text}`} aria-label={`${label} ${text}`}
      className={`inline-flex cursor-pointer items-center gap-1 rounded-md p-1 text-xs text-base-content/40 transition-colors hover:bg-base-200 hover:text-base-content ${className}`}
      onClick={() => navigator.clipboard.writeText(text).then(() => setDone(true)).catch(() => {})}>
      {done ? <><Check className="size-3.5 text-success" />Copied</> : <Copy className="size-3.5" />}
    </button>
  )
}
