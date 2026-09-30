import { useEffect } from 'react'

/** Read a DaisyUI theme token as a resolved color string. */
function themeColor(token: string): string {
  const probe = document.createElement('span')
  probe.style.color = `var(${token})`
  probe.style.display = 'none'
  document.body.appendChild(probe)
  const c = getComputedStyle(probe).color
  probe.remove()
  return c
}

/** The browser's theme-color follows the active theme's surface. The favicon itself is the
 *  static `public/favicon.svg` (lucide "utensils-crossed", ISC) in a neutral ink and does
 *  not change with the theme. */
export function useThemedFavicon() {
  useEffect(() => {
    const apply = () => {
      let meta = document.querySelector<HTMLMetaElement>('meta[name="theme-color"]')
      if (!meta) { meta = document.createElement('meta'); meta.name = 'theme-color'; document.head.appendChild(meta) }
      meta.content = themeColor('--color-base-100')
    }
    apply()
    const obs = new MutationObserver(apply)
    obs.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] })
    return () => obs.disconnect()
  }, [])
}
