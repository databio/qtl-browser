import { useEffect, useState } from 'react'
import { colocLists, phenotypeLabel, type ColocList } from '@/lib/store'
import { toastError } from '@/lib/toast'

/**
 * The store's coloc lists (SPEC section 20, `colocLists()` in lib/store.ts) for the open
 * collection and the chosen GWAS, read for the genome track and the gene chips. Nothing here names
 * a gene: which loci colocalize is data, and the store says who reported it.
 */

/** A list's QTL label: eQTL or sQTL by its experiment's phenotype type, else the type itself. */
export const colocTrait = (l: ColocList): string => phenotypeLabel(l.experiment.phenotype_type)

/** Gene id -> the lists that name it. */
export function colocGenes(lists: ColocList[]): Map<string, ColocList[]> {
  const out = new Map<string, ColocList[]>()
  for (const l of lists) for (const g of new Set(l.rows.map(r => r.gene_id))) out.set(g, [...(out.get(g) ?? []), l])
  return out
}

/** The list of that QTL label ('eQTL', 'sQTL') naming the gene, or null; null while the lists load. */
export const colocListing = (lists: ColocList[] | null, geneId: string, trait: string): ColocList | null =>
  (lists ?? []).find(l => colocTrait(l) === trait && l.rows.some(r => r.gene_id === geneId)) ?? null

/** Who reported a list, for a tooltip: its source's `reported_by`, else its id. */
export const colocSource = (l: ColocList): string => l.doc.source?.reported_by ?? l.doc.id

/** The lists for the open collection and the chosen GWAS; null while loading, empty when the store has none
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
