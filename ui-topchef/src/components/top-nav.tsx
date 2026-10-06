import { Link, NavLink } from 'react-router'
import { Moon, Sun, SunMoon } from 'lucide-react'
import { useTheme } from '@/contexts/theme-context'
import { CONTAINER } from '@/components/page'

const ITEMS = [
  { to: '/', label: 'Search', end: true },
  { to: '/genes', label: 'Genes' },
  { to: '/about', label: 'About' },
]

/** Top navbar: wordmark left, text links and theme toggle right. Shares CONTAINER with Page
 *  so the navbar and body always have the same x inset. */
export function TopNav() {
  const { mode, cycle } = useTheme()
  const ThemeIcon = mode === 'light' ? Sun : mode === 'dark' ? Moon : SunMoon
  return (
    <header className="border-b border-base-300 bg-base-100">
      <div className={`${CONTAINER} flex h-13 items-center gap-4`}>
        <Link to="/" className="text-xl font-thin tracking-wide">
          <span className="text-base-content/50">topchef</span><span className="text-base-content/30">.</span><span className="font-normal text-primary">qtl</span>
        </Link>
        <div className="flex-1" />
        <nav className="flex items-center gap-1">
          {ITEMS.map(({ to, label, end }) => (
            <NavLink key={to} to={to} end={end}
              className={({ isActive }) => `rounded-lg px-3 py-1.5 text-sm transition-colors ${isActive ? 'bg-base-200 font-medium text-base-content' : 'text-base-content/70 hover:bg-base-200 hover:text-base-content'}`}>
              {label}
            </NavLink>
          ))}
        </nav>
        <button onClick={cycle} title={`Theme: ${mode} (click to change)`}
          className="-mr-1.5 inline-flex cursor-pointer items-center justify-center rounded-md p-1.5 text-base-content/60 transition-colors hover:bg-base-200 hover:text-base-content">
          <ThemeIcon className="size-4" />
        </button>
      </div>
    </header>
  )
}
