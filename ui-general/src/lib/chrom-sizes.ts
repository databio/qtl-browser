/** Chromosome sizes from the seqcol API, cached per build in localStorage (as in pegasus-v2f-ui). */

export type ChromSizes = { names: string[]; lengths: number[] }

const SEQCOL_API = 'https://seqcolapi.databio.org'
/** The seqcol collection the store's variant catalogs are anchored to (`collection_digest`), so the
 *  track's chromosome lengths and the data's positions cite one reference. `npm run store-check`
 *  fails when a store disagrees. The GRCh38 no-alt analysis set: 195 sequences, and identical to the
 *  455-sequence collection this used before on every standard chromosome. */
const SEQCOL_DIGESTS: Record<string, string> = {
  hg38: 'EiFob05aCWgVU_B_Ae0cypnQut3cxUP1',
  GRCh38: 'EiFob05aCWgVU_B_Ae0cypnQut3cxUP1',
}
export const seqcolDigest = (genomeBuild = 'GRCh38'): string | undefined => SEQCOL_DIGESTS[genomeBuild]
const STANDARD_CHROMS = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), 'chrX', 'chrY']

/** The collection this reads, for citing on the About page: the names and lengths themselves, since
 *  the API's root path serves nothing. */
export const seqcolCollectionUrl = (genomeBuild = 'GRCh38') =>
  `${SEQCOL_API}/collection/${SEQCOL_DIGESTS[genomeBuild]}?level=2`
const CACHE_KEY_PREFIX = 'topchef.chromSizes.'

export async function fetchChromSizes(genomeBuild = 'GRCh38'): Promise<ChromSizes> {
  const digest = SEQCOL_DIGESTS[genomeBuild]
  if (!digest) throw new Error(`No seqcol digest known for genome build '${genomeBuild}'`)
  // the digest is in the key: a build that cites another collection must not read lengths cached
  // from the old one, which never expire
  const cacheKey = `${CACHE_KEY_PREFIX}${genomeBuild}.${digest}`
  const cached = localStorage.getItem(cacheKey)
  if (cached) {
    try { return JSON.parse(cached) as ChromSizes } catch { /* re-fetch */ }
  }
  const resp = await fetch(`${SEQCOL_API}/collection/${digest}?level=2`)
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
