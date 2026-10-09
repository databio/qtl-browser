import { useEffect, useState } from 'react'
import { colocLists, EQTL_TYPE, SQTL_TYPE, type ColocList } from '@/lib/store'
import { toastError } from '@/lib/toast'

/**
 * The store's coloc lists (SPEC section 20, `colocLists()` in lib/store.ts), read for the Home
 * track, the gene chips and About. Nothing here names a gene: which loci colocalize is data, and
 * the store says who reported it.
 */

/** A list's QTL label: eQTL or sQTL by its experiment's phenotype type, else the type itself. */
export function colocTrait(l: ColocList): string {
  const t = l.experiment.phenotype_type
  return t === EQTL_TYPE ? 'eQTL' : t === SQTL_TYPE ? 'sQTL' : t
}

/** Gene id -> the lists that name it. */
export function colocGenes(lists: ColocList[]): Map<string, ColocList[]> {
  const out = new Map<string, ColocList[]>()
  for (const l of lists) for (const g of new Set(l.rows.map(r => r.gene_id))) out.set(g, [...(out.get(g) ?? []), l])
  return out
}

/** Whether a list of that QTL label ('eQTL', 'sQTL') names the gene; false while the lists load. */
export const colocListed = (lists: ColocList[] | null, geneId: string, trait: string): boolean =>
  (lists ?? []).some(l => colocTrait(l) === trait && l.rows.some(r => r.gene_id === geneId))

/** Distinct genes the lists of that QTL label name (gene-level and phenotype-level rows count the
 *  same way): Home's and About's "21 eGenes and 4 sGenes". */
export function colocGeneCount(lists: ColocList[] | null, trait: string): number {
  return new Set((lists ?? []).filter(l => colocTrait(l) === trait).flatMap(l => l.rows.map(r => r.gene_id))).size
}

/** The lists for the open study and its GWAS; null while loading, empty when the store has none
 *  (or they could not be read, which is a toast: the page works without them). */
export function useColocLists(): ColocList[] | null {
  const [lists, setLists] = useState<ColocList[] | null>(null)
  useEffect(() => {
    let alive = true
    colocLists()
      .then(l => { if (alive) setLists(l) })
      .catch(e => { toastError('The coloc lists could not be loaded.', e); if (alive) setLists([]) })
    return () => { alive = false }
  }, [])
  return lists
}
