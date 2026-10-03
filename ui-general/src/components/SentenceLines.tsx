import { useLayoutEffect, useRef, useState, type ReactNode } from 'react'

/**
 * Sentences on separate lines while every sentence fits on one line; as soon as any of them
 * would wrap, the whole thing flows as one block. CSS cannot observe wrapping, so a hidden
 * copy laid out with line breaks is measured on every resize: a sentence whose text spans more
 * than one line box means "too narrow".
 */
export default function SentenceLines({ sentences, className = '' }: { sentences: ReactNode[]; className?: string }) {
  const probe = useRef<HTMLParagraphElement>(null)
  const [block, setBlock] = useState(false)

  // layout effect: measured and decided before the first paint, so the text never shows in
  // one arrangement and then switches
  useLayoutEffect(() => {
    const el = probe.current
    if (!el) return
    const check = () => {
      const wrapped = Array.from(el.querySelectorAll<HTMLElement>('[data-sentence]')).some(s => s.getClientRects().length > 1)
      setBlock(wrapped)
    }
    check()
    const obs = new ResizeObserver(check)
    obs.observe(el)
    return () => obs.disconnect()
  }, [sentences.length])

  const items = sentences.filter(s => s != null && s !== false)
  return (
    <div className="relative">
      {/* measuring copy: same width and type, never visible, never in the accessibility tree */}
      <p ref={probe} aria-hidden className={`invisible absolute inset-x-0 top-0 ${className}`}>
        {items.map((s, i) => <span key={i}><span data-sentence>{s}</span>{i < items.length - 1 && <br />}</span>)}
      </p>
      <p className={className}>
        {items.map((s, i) => <span key={i}>{s}{i < items.length - 1 && (block ? ' ' : <br />)}</span>)}
      </p>
    </div>
  )
}
