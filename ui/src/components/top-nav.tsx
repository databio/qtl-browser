import { NavLink } from 'react-router'
import { Moon, Sun, SunMoon } from 'lucide-react'
import { useTheme } from '@/contexts/theme-context'
import { CONTAINER } from '@/components/page'
import { SCOPE } from '@/lib/scope'

const LINK = 'rounded-lg px-3 py-1.5 text-sm transition-colors'
const IDLE = `${LINK} text-base-content/70 hover:bg-base-200 hover:text-base-content`
const ACTIVE = `${LINK} bg-base-200 font-medium text-base-content`

/** Top navbar: wordmark left, text links and theme toggle right. Shares CONTAINER with Page
 *  so the navbar and body always have the same x inset. */
export function TopNav() {
  const { mode, cycle } = useTheme()
  const ThemeIcon = mode === 'light' ? Sun : mode === 'dark' ? Moon : SunMoon
  // the same two links everywhere; where you are inside a collection is the page's crumbs. Inside
  // one, the router's basename is the collection (lib/scope.ts), so both are plain anchors: leaving
  // a collection is a full navigation
  const items = [{ to: '/', label: 'Collections', end: true }, { to: '/about', label: 'About' }]
  return (
    <header className="border-b border-base-300 bg-base-100">
      <div className={`${CONTAINER} flex h-13 items-center gap-4`}>
        <a href="/" className="text-xl font-thin tracking-wide">
          <span className="text-base-content/50">qtl</span><span className="text-base-content/30">.</span><span className="font-normal text-primary">browser</span>
        </a>
        <div className="flex-1" />
        <nav className="flex items-center gap-1">
          {items.map(({ to, label, end }) => SCOPE
            ? <a key={to} href={to} className={IDLE}>{label}</a>
            : <NavLink key={to} to={to} end={end} className={({ isActive }) => (isActive ? ACTIVE : IDLE)}>{label}</NavLink>)}
        </nav>
        <button onClick={cycle} title={`Theme: ${mode} (click to change)`}
          className="-mr-1.5 inline-flex cursor-pointer items-center justify-center rounded-md p-1.5 text-base-content/60 transition-colors hover:bg-base-200 hover:text-base-content">
          <ThemeIcon className="size-4" />
        </button>
      </div>
    </header>
  )
}
