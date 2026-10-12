/**
 * The collection, experiment and GWAS a page belongs to, read once from the URL at load.
 *
 * `/c/<digest>` is a collection's own page. `/c/<digest>/e/<experiment id>` opens one experiment
 * of it, and `/c/<digest>/e/<id>/gwas/<gwas id>` opens it with that GWAS for coloc and
 * LocusCompare. The prefix is the router's basename (main.tsx), so every in-app link (`/gene/FLNC`)
 * stays inside it. Moving to another experiment or collection, choosing another GWAS, or leaving
 * for the landing page is a full navigation: every reader below caches per page load (genes,
 * blocks, the GWAS window), so one load never mixes two experiments or two GWAS.
 *
 * Nothing is preselected: an experiment opens with no GWAS until the visitor picks one.
 */

/** A 32-character sha512t24u digest, the name of a collection object without `.json`. */
const DIGEST = /^[A-Za-z0-9_-]{32}$/

export interface Scope {
  /** the collection's digest, as in the URL */
  digest: string
  /** the collection object's name, `<digest>.json`, as experiments name it in `store.json` */
  collection: string
  /** the experiment's id, or null on the collection's own page */
  experiment: string | null
  /** the chosen GWAS experiment's id, or null */
  gwas: string | null
  /** the URL prefix every in-app path sits under */
  base: string
}

/** The URL prefix of a collection, or of one of its experiments, with a GWAS when one is chosen. */
export const collectionPath = (digest: string, experiment?: string | null, gwas?: string | null): string =>
  `/c/${digest}${experiment ? `/e/${encodeURIComponent(experiment)}${gwas ? `/gwas/${encodeURIComponent(gwas)}` : ''}` : ''}`

/** A collection object name (`<digest>.json`) to its digest. */
export const digestOf = (collection: string): string => collection.replace(/\.json$/, '')

function parse(path: string): Scope | null {
  const m = /^\/c\/([^/]+)(?:\/e\/([^/]+)(?:\/gwas\/([^/]+))?)?(?=\/|$)/.exec(path)
  if (!m || !DIGEST.test(m[1]!)) return null
  const experiment = m[2] ? decodeURIComponent(m[2]) : null
  const gwas = m[3] ? decodeURIComponent(m[3]) : null
  return { digest: m[1]!, collection: `${m[1]}.json`, experiment, gwas, base: collectionPath(m[1]!, experiment, gwas) }
}

export const SCOPE: Scope | null = parse(window.location.pathname)

/** Reopens the current page with another GWAS (or none): a full navigation, see above. */
export function chooseGwas(gwas: string | null): void {
  if (!SCOPE?.experiment) return
  const rest = window.location.pathname.slice(SCOPE.base.length)
  window.location.assign(collectionPath(SCOPE.digest, SCOPE.experiment, gwas) + rest + window.location.search)
}
