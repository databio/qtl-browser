import { createContext, useContext, useEffect, useMemo, useState, type ReactNode } from 'react'

export type Mode = 'light' | 'dark' | 'system'
const KEY = 'topchef-theme'
const Ctx = createContext<{ mode: Mode; cycle: () => void }>({ mode: 'system', cycle: () => {} })

function apply(mode: Mode) {
  const dark = mode === 'dark' || (mode === 'system' && window.matchMedia('(prefers-color-scheme: dark)').matches)
  document.documentElement.setAttribute('data-theme', dark ? 'topchef-dark' : 'topchef')
}

/** light → dark → system, persisted. Mirrors atlas' settings-context theme cycling. */
export function ThemeProvider({ children }: { children: ReactNode }) {
  const [mode, setMode] = useState<Mode>(() => (localStorage.getItem(KEY) as Mode) || 'system')
  useEffect(() => {
    apply(mode)
    localStorage.setItem(KEY, mode)
    const mq = window.matchMedia('(prefers-color-scheme: dark)')
    const onChange = () => mode === 'system' && apply(mode)
    mq.addEventListener('change', onChange)
    return () => mq.removeEventListener('change', onChange)
  }, [mode])
  const value = useMemo(() => ({
    mode,
    cycle: () => setMode(m => (m === 'light' ? 'dark' : m === 'dark' ? 'system' : 'light')),
  }), [mode])
  return <Ctx.Provider value={value}>{children}</Ctx.Provider>
}

export const useTheme = () => useContext(Ctx)
