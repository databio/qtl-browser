import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { SectionPanel } from '@/components/section-panel'
import { Empty, TableSkeleton } from '@/components/states'
import { CollectionRows } from '@/routes/Landing'
import { useCollection, useIndex } from '@/contexts/store-context'
import { doiUrl, experimentLabel, experimentSize, facetRows, shortName } from '@/lib/describe'
import { KvTable, type KvRow } from '@/components/kv-table'
import { collectionPath } from '@/lib/scope'
import { hrefRowLink, ROW_LINK, ROW_LINK_TEXT } from '@/lib/row-link'
import { significanceOf } from '@/lib/significance'
import type { CollectionDoc, ExperimentDoc, IndexCollection } from '@/lib/store'

/** A collection's own page: its publication and its experiments, each a link to the experiment's
 *  pages (routes/Experiment.tsx), or a parent's children. */
export default function Collection() {
  const ix = useIndex()
  const c = useCollection()
  if (!ix) return <Page><TableSkeleton columns={[{ w: 'w-48' }, { w: 'w-32' }]} rows={4} /></Page>
  if (!c) return <NotInStore />
  return (
    <Page>
      <PageHeader title={shortName(c.doc)} crumbs={[{ href: '/', label: 'Collections' }, { label: shortName(c.doc) }]}
 />
      <div className="space-y-10">
        {c.experiments.length > 0 && <Experiments c={c} />}
        <div className="grid items-start gap-4 md:grid-cols-2">
          <KvTable title="Publication" rows={publicationRows(c.doc)} />
          {/* the data: what the experiments say they measured, then where the files came from */}
          <KvTable title="Metadata" rows={[
            ...facetRows(c.experiments),
            ...sourceRows(c.doc.source)]} />
        </div>
        {c.children.length > 0 && (
          <SectionPanel title="Collections">
            <CollectionRows list={c.children.map(n => ix.collections.get(n)).filter((k): k is IndexCollection => !!k)} ix={ix} />
          </SectionPanel>
        )}
      </div>
    </Page>
  )
}

const str = (v: unknown): string | null => (typeof v === 'string' && v ? v : typeof v === 'number' ? String(v) : null)

/** The citation as rows: title, authors, where and when, DOI. */
function publicationRows(d: CollectionDoc): KvRow[] {
  const ct = d.citation ?? {}
  const authors = [str(ct.first_author), str(ct.last_author)].filter(Boolean).join(' … ')
  const where = [str(ct.venue), str(ct.volume), str(ct.pages), str(ct.year)].filter(Boolean).join(', ')
  const doi = doiUrl(d)
  const rows: (KvRow | null)[] = [
    str(ct.title) ? { label: 'Title', value: str(ct.title) } : null,
    authors ? { label: 'Authors', value: authors } : null,
    where ? { label: 'Published in', value: where } : null,
    str(ct.posted) ? { label: 'Posted', value: str(ct.posted) } : null,
    doi ? { label: 'DOI', value: <ExternalLink icon href={doi}>{String(ct.doi)}</ExternalLink> } : null,
  ]
  rows.push(sourcePublication(d.source))
  return rows.filter((r): r is KvRow => r != null)
}

/** Keys a `source` holds when the source is itself a publication (the eQTL Catalogue's
 *  reprocessing, for GTEx): they go under Publication, not Data. */
const PUB_KEYS = ['title', 'doi', 'first_author', 'last_author', 'venue', 'volume', 'pages', 'year']

/** A source that is itself a publication, as one Publication row: its title, authors, venue and
 *  year, and its DOI as a link. Null when the source names no title or DOI. */
function sourcePublication(src: Record<string, unknown> | undefined): KvRow | null {
  if (!src || (!str(src.title) && !str(src.doi))) return null
  const authors = [str(src.first_author), str(src.last_author)].filter(Boolean).join(' … ')
  const where = [str(src.venue), str(src.volume), str(src.pages), str(src.year)].filter(Boolean).join(', ')
  return { label: 'Data processed in', value: <span className="space-y-0.5">
    {str(src.title) && <span className="block">{str(src.title)}</span>}
    {(authors || where) && <span className="block">{[authors, where].filter(Boolean).join(' · ')}</span>}
    {str(src.doi) && <ExternalLink icon href={`https://doi.org/${str(src.doi)}`}>{str(src.doi)}</ExternalLink>}
  </span> }
}

/** Where the data came from, as the collection document records it: its plain fields, a URL as a
 *  link. A source's own publication fields are under Publication instead. */
function sourceRows(src: Record<string, unknown> | undefined): KvRow[] {
  if (!src) return []
  const label = (k: string) => k.charAt(0).toUpperCase() + k.slice(1).replace(/_/g, ' ')
  const out: KvRow[] = []
  for (const [k, v] of Object.entries(src)) {
    const t = str(v)
    if (!t || PUB_KEYS.includes(k)) continue
    out.push(k === 'url' ? { label: 'URL', value: <ExternalLink icon href={t}>{t}</ExternalLink> } : { label: label(k), value: t })
  }
  return out
}

function Experiments({ c }: { c: IndexCollection }) {
  return (
    <SectionPanel title="Experiments">
      <div className="overflow-x-auto rounded-lg border border-base-300">
        <table className="table table-sm text-xs">
          <thead><tr><th>Experiment</th><th>Id</th><th>Size</th><th>Significance</th></tr></thead>
          <tbody>
            {c.experiments.map(e => (
              // the whole row opens the experiment, a full navigation (lib/scope.ts); the label stays an
              // anchor so it can be copied or opened like any link
              <tr key={e.id} className={`${ROW_LINK} hover:bg-base-200`} {...hrefRowLink(collectionPath(c.digest, e.id))}>
                <td><a className={`font-medium ${ROW_LINK_TEXT}`} href={collectionPath(c.digest, e.id)} onClick={ev => ev.stopPropagation()}>{experimentLabel(e)}</a></td>
                <td className="font-mono text-base-content/55">{e.id}</td>
                <td className="tabular-nums">{experimentSize(e)}</td>
                <td className="text-base-content/70">{rule(e)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </SectionPanel>
  )
}

const rule = (e: ExperimentDoc): string => (e.kind === 'qtl_cis' ? significanceOf(e).ruleText || 'none' : '')

/** A collection digest this store does not hold: an old link to a collection whose metadata has
 *  since changed its digest, or a mistyped one. */
export function NotInStore({ what = 'collection' }: { what?: 'collection' | 'experiment' }) {
  return (
    <Page>
      <Empty label={what === 'collection'
        ? 'This collection is not in the store. Its link may predate an edit to the collection\'s details.'
        : 'This experiment is not in this collection. Its link may predate an edit to the collection\'s details.'} />
      <p className="px-4 text-sm"><a className="link-quiet" href="/">See the collections in this store</a></p>
    </Page>
  )
}
