import { ChevronRight } from 'lucide-react'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { SectionPanel } from '@/components/section-panel'
import { TableSkeleton } from '@/components/states'
import { useIndex } from '@/contexts/store-context'
import { citationLine, collectionKind, facetLine, gwasSummary, shortName, studySummary } from '@/lib/describe'
import { collectionPath } from '@/lib/scope'
import type { GwasExperiment, IndexCollection, StoreIndex } from '@/lib/store'

/** The store's collections in two groups: the QTL studies, and the GWAS any of them can be read
 *  against. Nothing is opened until one is chosen. */
export default function Landing() {
  const ix = useIndex()
  const studies = ix?.roots.filter(c => collectionKind(c) !== 'gwas') ?? []
  const gwas = ix?.roots.filter(c => collectionKind(c) === 'gwas') ?? []
  return (
    <Page>
      <PageHeader title="Collections" />
      {!ix ? <TableSkeleton columns={[{ w: 'w-48' }, { w: 'w-32' }]} rows={3} /> : (
        <div className="space-y-10">
          <SectionPanel title={<>Studies <span className="tabular-nums text-base-content/40">{studies.length}</span></>}>
            <CollectionRows list={studies} ix={ix} />
          </SectionPanel>
          <SectionPanel title={<>GWAS <span className="tabular-nums text-base-content/40">{gwas.length}</span></>}>
            <CollectionRows list={gwas} ix={ix} />
          </SectionPanel>
        </div>
      )}
    </Page>
  )
}

export function CollectionRows({ list, ix }: { list: IndexCollection[]; ix: StoreIndex }) {
  if (!list.length) return <p className="text-xs text-base-content/45">None in this store.</p>
  return <div className="space-y-2">{list.map(c => <CollectionRow key={c.name} c={c} ix={ix} />)}</div>
}

/** One collection as a row: its name, citation and contents on the first line, its facets dimmed on
 *  the second. A full navigation, since a collection is opened per page load (lib/scope.ts). */
function CollectionRow({ c, ix }: { c: IndexCollection; ix: StoreIndex }) {
  const g = c.experiments.find((e): e is GwasExperiment => e.kind === 'gwas')
  const kind = collectionKind(c)
  const stats = kind === 'qtl' ? studySummary(c)
    : kind === 'gwas' && g ? gwasSummary(g)
    : `${c.children.length} collection${c.children.length === 1 ? '' : 's'}`
  const cite = citationLine(c.doc)
  const facets = facetLine(c.experiments)
  return (
    <a href={collectionPath(c.digest)}
      className="group flex items-center gap-3 rounded-lg border border-base-300 px-4 py-3 transition-colors hover:border-base-content/30 hover:bg-base-200/40">
      <div className="min-w-0 flex-1">
        <div className="flex items-baseline gap-3">
          <span className="truncate font-medium">{shortName(c.doc)}</span>
          {cite && <span className="truncate text-xs text-base-content/45">{cite}</span>}
          <span className="ml-auto shrink-0 text-xs tabular-nums text-base-content/55">{stats}</span>
        </div>
        {facets && <div className="mt-0.5 truncate text-xs text-base-content/50">{facets}</div>}
        {c.children.length > 0 && (
          <div className="mt-0.5 truncate text-xs text-base-content/45">
            {c.children.map(n => ix.collections.get(n)).filter(Boolean).map(k => shortName(k!.doc)).join(', ')}
          </div>
        )}
      </div>
      <ChevronRight className="size-4 shrink-0 text-base-content/25 transition-colors group-hover:text-base-content/60" />
    </a>
  )
}
