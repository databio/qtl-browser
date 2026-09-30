/**
 * Chart colors for credible-set membership (categorical, one slot per set, fixed).
 *
 * Hues are the dataviz reference palette's blue, aqua, yellow, red, plus a purple (OKLCH
 * h320 L0.45 C0.18 light, h315 L0.65 C0.18 dark) in place of the palette's green: as an 8%
 * table-row tint the green was indistinguishable from aqua (OKLab ΔE 1.4 light / 2.5 dark;
 * the purple gives 3.0 / 5.1, and its nearest tint neighbor is no closer than any other pair).
 * No five-hue subset of the palette passes the all-pairs checks on both surfaces (best: light
 * CVD ΔE 6.9 / normal 15.6, dark CVD 6.5 / normal 11.9, run 2026-09-04 against the theme
 * surfaces #ffffff and #1f1615), so each set also gets its own marker shape as the secondary
 * encoding, and the credible-set table below the plot is the table view. Variants in no set
 * are the neutral background, a warm gray matched to each surface.
 */
import { useEffect, useState } from 'react'

export const CS_DOMAIN = ['none', '1', '2', '3', '4', '5'] as const
export const CS_COLORS = {
  light: ['#b8b3b1', '#2a78d6', '#1baf7a', '#eda100', '#7e2491', '#e34948'],
  dark: ['#5a4d4a', '#3987e5', '#199e70', '#c98500', '#b767d9', '#e66767'],
}
// the neutral background keeps the plain circle; every credible set gets its own shape
export const CS_SYMBOLS = ['circle', 'diamond2', 'square', 'triangle', 'star', 'hexagon']   // diamond2 = rotated square
/** CSS clip-paths that echo the plot symbols in the legend swatches. */
export const CS_SWATCH_CLIP: Record<string, string | undefined> = {
  circle: undefined,
  square: 'inset(10%)',
  triangle: 'polygon(50% 5%, 95% 95%, 5% 95%)',
  diamond2: 'polygon(50% 0%, 100% 50%, 50% 100%, 0% 50%)',
  star: 'polygon(50% 0%, 61% 35%, 98% 35%, 68% 57%, 79% 91%, 50% 70%, 21% 91%, 32% 57%, 2% 35%, 39% 35%)',
  hexagon: 'polygon(25% 5%, 75% 5%, 100% 50%, 75% 95%, 25% 95%, 0% 50%)',
}
export function isDark(): boolean {
  return document.documentElement.getAttribute('data-theme') === 'topchef-dark'
}

/** Current theme as state, updated when the theme toggle changes `data-theme`. */
export function useIsDark(): boolean {
  const [dark, setDark] = useState(isDark)
  useEffect(() => {
    const obs = new MutationObserver(() => setDark(isDark()))
    obs.observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] })
    return () => obs.disconnect()
  }, [])
  return dark
}

/** Row background for a credible set: the set's plot color at low alpha, so table rows and
 *  plot markers read as the same thing. `undefined` for variants in no set. */
/** Row tint alphas per credible set (index as in CS_DOMAIN; 0 is unused). The hover alpha
 *  starts from a per-color solve that matches the OKLab lightness step of a plain row's hover
 *  (base-200 at 60% over base-100: light 0.0195, dark 0.0158), then is nudged by eye, mostly
 *  because the lightness match leaves yellow, and to a lesser degree aqua, reading stronger
 *  than the rest. Solved values, sets 1-5: light 1f/22/27/1d/20, dark 28/29/27/28/27.
 *  Re-solve and re-check by eye if CS_COLORS or the theme surfaces change. */
const CS_TINT_ALPHA = {
  light: { rest: '14', hover: ['', '1f', '21', '22', '1d', '20'] },
  dark: { rest: '1f', hover: ['', '28', '28', '25', '27', '27'] },
}
export function csTint(csId: number | string | null | undefined, dark: boolean, hover = false): string | undefined {
  if (csId == null) return undefined
  const i = CS_DOMAIN.indexOf(String(csId) as (typeof CS_DOMAIN)[number])
  if (i <= 0) return undefined
  const a = CS_TINT_ALPHA[dark ? 'dark' : 'light']
  return CS_COLORS[dark ? 'dark' : 'light'][i] + (hover ? a.hover[i] : a.rest)
}
