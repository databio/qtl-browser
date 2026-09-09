import type { KeyboardEvent, MouseEvent } from 'react'
import { useHref, useNavigate } from 'react-router'

/** Classes for a row that acts as a link: pointer, hover transition, keyboard focus ring. Tables
 *  add their own hover tint on top. */
export const ROW_LINK = 'cursor-pointer transition-colors focus-visible:outline-none focus-visible:bg-base-200/60'

/** For a span around the row's identifier text (rsID, gene symbol): underlines while the text
 *  itself is hovered, the way a quiet link underlines on hover. */
export const ROW_LINK_TEXT = 'underline-offset-2 hover:underline'

/** Makes a table row behave like an anchor to an in-app path. An anchor cannot be or wrap a
 *  `<tr>`, so the row carries handlers that copy a link's behavior instead: plain click
 *  navigates in place, cmd/ctrl/shift click and middle click open a new tab, Enter on the focused
 *  row navigates. Usage: `const rowLink = useRowLink(); <tr {...rowLink('/gene/X')}>`. */
export function useRowLink() {
  const navigate = useNavigate()
  const base = useHref('/')
  const url = (to: string) => base.replace(/\/$/, '') + to
  return (to: string) => ({
    role: 'link' as const,
    tabIndex: 0,
    onClick: (e: MouseEvent<HTMLTableRowElement>) => {
      if (e.defaultPrevented) return
      if (e.metaKey || e.ctrlKey || e.shiftKey) window.open(url(to), '_blank', 'noopener')
      else navigate(to)
    },
    onAuxClick: (e: MouseEvent<HTMLTableRowElement>) => {
      if (e.button !== 1) return
      e.preventDefault()
      window.open(url(to), '_blank', 'noopener')
    },
    onKeyDown: (e: KeyboardEvent<HTMLTableRowElement>) => {
      if (e.key === 'Enter' && e.target === e.currentTarget) navigate(to)
    },
  })
}
