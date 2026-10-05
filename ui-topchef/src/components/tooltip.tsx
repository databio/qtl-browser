import { useLayoutEffect, useRef, useState, type ReactNode } from 'react'
import { createPortal } from 'react-dom'

const GAP = 6
const EDGE = 8

/** Hover/focus tooltip rendered in a portal and positioned against the viewport: centered under
 *  the trigger, shifted to stay on screen, flipped above when there is no room below. Wraps at a
 *  fixed measure instead of stretching. */
export function Tooltip({ tip, children, className = '' }: { tip: string; children: ReactNode; className?: string }) {
  const ref = useRef<HTMLSpanElement>(null)
  const bubble = useRef<HTMLDivElement>(null)
  const [open, setOpen] = useState(false)
  const [pos, setPos] = useState<{ top: number; left: number; above: boolean } | null>(null)

  useLayoutEffect(() => {
    if (!open || !ref.current || !bubble.current) { setPos(null); return }
    const t = ref.current.getBoundingClientRect()
    const b = bubble.current.getBoundingClientRect()
    const vw = window.innerWidth, vh = window.innerHeight
    let left = t.left + t.width / 2 - b.width / 2
    left = Math.max(EDGE, Math.min(left, vw - EDGE - b.width))
    const below = t.bottom + GAP + b.height <= vh - EDGE
    const top = below ? t.bottom + GAP : t.top - GAP - b.height
    setPos({ top, left, above: !below })
  }, [open])

  return (
    <span ref={ref} className={`inline-flex ${className}`} tabIndex={0}
      onMouseEnter={() => setOpen(true)} onMouseLeave={() => setOpen(false)}
      onFocus={() => setOpen(true)} onBlur={() => setOpen(false)}>
      {children}
      {open && createPortal(
        <div ref={bubble} role="tooltip"
          className={`pointer-events-none fixed z-50 max-w-72 rounded-md bg-neutral px-2.5 py-1.5 text-xs leading-snug text-neutral-content shadow-md ${pos ? 'opacity-100' : 'opacity-0'}`}
          /* inline style: position is measured at runtime against the viewport */
          style={{ top: pos?.top ?? 0, left: pos?.left ?? 0 }}>
          {tip}
        </div>,
        document.body,
      )}
    </span>
  )
}
