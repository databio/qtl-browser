import type { ReactNode } from 'react'

/** One horizontal inset for the navbar and every page: viewport-proportional padding,
 *  capped width. Nothing else sets its own x padding, so edges line up on every screen. */
export const CONTAINER = 'mx-auto w-full max-w-[88rem] px-[max(1.5rem,4vw)]'

export function Page({ children }: { children: ReactNode }) {
  return <div className={`${CONTAINER} pt-6 pb-12`}>{children}</div>
}
