/** Chromosome sizes from the seqcol API, cached per collection in localStorage (as in pegasus-v2f-ui). */

export type ChromSizes = { names: string[]; lengths: number[] }

const SEQCOL_API = 'https://seqcolapi.databio.org'
const STANDARD_CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY']

/** The collection's names and lengths at the seqcol API, for citing: the API's root path serves nothing. */
export const seqcolCollectionUrl = (digest: string) => `${SEQCOL_API}/collection/${digest}?level=2`
const CACHE_KEY_PREFIX = 'qtlb.chromSizes.'

/** The standard chromosomes of the sequence collection `digest` names: the store's
 *  `collection_digest`, so the track's lengths and the data's positions cite one reference. */
export async function fetchChromSizes(digest: string): Promise<ChromSizes> {
  // the digest is the key: lengths cached for one collection never answer for another
  const cacheKey = `${CACHE_KEY_PREFIX}${digest}`
  const cached = localStorage.getItem(cacheKey)
  if (cached) {
    try { return JSON.parse(cached) as ChromSizes } catch { /* re-fetch */ }
  }
  const resp = await fetch(seqcolCollectionUrl(digest))
  if (!resp.ok) throw new Error(`seqcol fetch failed: ${resp.status}`)
  const data = (await resp.json()) as { names: string[]; lengths: number[] }
  const lookup = new Map<string, number>()
  data.names.forEach((n, i) => lookup.set(n, data.lengths[i]!))
  const names: string[] = [], lengths: number[] = []
  for (const chrom of STANDARD_CHROMS) {
    const len = lookup.get(chrom)
    if (len != null) { names.push(chrom); lengths.push(len) }
  }
  const result = { names, lengths }
  localStorage.setItem(cacheKey, JSON.stringify(result))
  return result
}
