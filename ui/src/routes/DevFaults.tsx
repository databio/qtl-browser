import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { FAULT, FAULT_ONLY, FAULTS, type Fault } from '@/lib/faults'
import { SCOPE } from '@/lib/scope'

/**
 * `/dev/faults`, dev builds only: one link per error surface, each a full page load with
 * `?fault=` (lib/faults.ts) so the fault is in force from the first request. Two displays exist:
 * the error page (components/error-page.tsx) for a failure that takes the page down, and a
 * persistent toast (lib/toast.ts) for one that leaves it usable; "Shows" names which.
 */
interface Trigger { fault: Fault; only?: string; path: string; surface: string; today: string }

const TRIGGERS: { title: string; rows: Trigger[] }[] = [
  { title: 'The store cannot open (App.tsx: the error page in place of every route)', rows: [
    { fault: 'offline', path: '/', surface: 'store.json fetch rejects', today: 'error page with the fetch error text' },
    { fault: 'store-404', path: '/', surface: 'store.json answers 404', today: 'error page: "store.json: HTTP 404"' },
    { fault: 'store-version', path: '/', surface: 'format version the build does not read', today: 'error page asking for a reload' },
    { fault: 'store-anchor', path: '/', surface: 'the configured collection is not in store.json', today: 'error page naming the missing experiment' },
    { fault: 'store-names', path: '/', surface: 'an experiment document names a non-object', today: 'error page: "bad object name hand-edited.qbe"' },
    { fault: 'object-404', only: 'json', path: '/', surface: 'the collection document is missing', today: 'error page: HTTP 404 on the collection' },
    { fault: 'object-corrupt', only: 'json', path: '/', surface: 'the collection document is unreadable', today: 'error page: "collection ... has no title"' },
  ] },
  { title: 'The query engine (App.tsx: the error page on the page that needed it)', rows: [
    { fault: 'db', path: '/gene/PDXK', surface: 'DuckDB-WASM fails to start', today: 'error page "The query engine could not start"' },
  ] },
  { title: 'The gene page (routes/Gene.tsx)', rows: [
    { fault: 'object-404', only: 'qbe', path: '/gene/PDXK', surface: 'the results block is missing', today: 'error page (the block is the page)' },
    { fault: 'object-corrupt', only: 'qbe', path: '/gene/PDXK', surface: 'the results block does not decode (PackError)', today: 'error page' },
    { fault: 'object-corrupt', only: 'qbv', path: '/gene/PDXK', surface: 'the variants range does not decode', today: 'error page' },
    { fault: 'object-404', only: 'qbg', path: '/gene/PDXK', surface: 'the GWAS window is missing', today: 'toast; the panel says it could not load' },
    { fault: 'object-404', only: 'qbt', path: '/gene/SMARCB1', surface: 'the trans frames are missing', today: 'toast; the trans table stays empty' },
    { fault: 'object-404', only: 'arrow.zst', path: '/gene/PDXK', surface: 'the search index part is missing', today: 'error page' },
    { fault: 'range-200', path: '/gene/PDXK', surface: 'the host ignores Range', today: 'error page' },
    { fault: 'slow', path: '/gene/PDXK', surface: 'every request takes 2.5 s', today: 'the skeletons, then the page' },
  ] },
  { title: 'The variant page (routes/Variant.tsx)', rows: [
    { fault: 'object-404', only: 'qbr', path: '/variant/rs34599497', surface: 'the rsID index is missing', today: 'error page (a failed lookup, not "not found")' },
    { fault: 'object-404', only: 'qbh', path: '/variant/rs34599497', surface: 'a hits file is missing', today: 'error page' },
    { fault: 'object-corrupt', only: 'qbe', path: '/variant/rs34599497', surface: 'the lead blocks do not decode', today: 'error page (the lists are the page)' },
    { fault: 'object-404', only: 'qbe', path: '/variant/rs34599497', surface: 'the cis scan\'s spans are missing (press Scan)', today: 'toast; the scan returns to "not scanned yet"' },
  ] },
  { title: 'Home, genes, region (routes/Home.tsx, Genes.tsx, Region.tsx, components/ColocLoci.tsx)', rows: [
    { fault: 'object-404', only: 'arrow.zst', path: '/', surface: 'the whole-genome tables are missing', today: 'toasts for the coloc loci and GWAS track; the counts never fill' },
    { fault: 'object-404', only: 'qgl', path: '/', surface: 'the gene lookup is missing', today: 'toast for the coloc loci; the search box finds nothing' },
    { fault: 'object-404', only: 'arrow.zst', path: '/genes', surface: 'the gene list cannot load', today: 'error page' },
    { fault: 'object-404', only: 'arrow.zst', path: '/region/chr21:43000000-44000000', surface: 'the chromosome\'s genes cannot load', today: 'error page' },
  ] },
  { title: 'Rendering and routing', rows: [
    { fault: 'render', path: '/', surface: 'a component throws during render', today: 'the error boundary shows the error page under the navbar' },
    { fault: 'none', path: '/no/such/page', surface: 'an unknown path', today: 'the error page with a Home button' },
  ] },
]

/** Not faults: the empty states a bad or untested id reaches today. */
const NOT_FOUND: { path: string; surface: string }[] = [
  { path: '/gene/NOPE', surface: 'no gene matches' },
  { path: '/gene/FP565260.4', surface: 'annotated, not tested for QTL' },
  { path: '/gene/GUSBP11', surface: 'tested for sQTL only ("no eQTL test" chip)' },
  { path: '/variant/rs1', surface: 'an rsID the catalog does not hold' },
  { path: '/variant/not-a-variant', surface: 'a malformed variant id' },
  { path: '/variant/chrY:1000000', surface: 'a chromosome the store does not hold' },
  { path: '/variant/rs457868', surface: 'a trans-only site, outside every cis window' },
  { path: '/region/garbage', surface: 'a malformed region' },
  { path: '/region/chr21:1-2', surface: 'a region with no genes' },
  { path: '/genes?x', surface: 'the gene list; type "zzzz" in the search box for "no genes match"' },
]

// open it inside a collection (/c/<digest>/dev/faults): the paths are that collection's pages
const at = (path: string) => `${SCOPE?.base ?? ''}${path}`
const href = (t: Trigger) => t.fault === 'none' ? at(t.path) : `${at(t.path)}${t.path.includes('?') ? '&' : '?'}fault=${t.fault}${t.only ? `&only=${t.only}` : ''}`

export default function DevFaults() {
  return (
    <Page>
      <div className="mx-auto max-w-5xl space-y-8">
        <PageHeader title="Fault triggers" className="mb-2 mt-2" />
        <p className="text-sm text-base-content/70">
          Each link is a full page load with <code>?fault=</code> in force from the first request (lib/faults.ts); open
          it in a new tab and come back here. The fault holds until the next load without it.
          {FAULT !== 'none' && <> Active now: <code>{FAULT}</code>{FAULT_ONLY && <> on <code>.{FAULT_ONLY}</code></>}.</>}
        </p>
        {TRIGGERS.map(s => (
          <section key={s.title} className="space-y-2">
            <h2 className="text-sm font-medium text-base-content/70">{s.title}</h2>
            <div className="overflow-hidden rounded-lg border border-base-300">
              <table className="table table-sm">
                <thead><tr><th>Fault</th><th>What fails</th><th>Shows</th></tr></thead>
                <tbody>
                  {s.rows.map(t => (
                    <tr key={href(t)}>
                      <td className="whitespace-nowrap"><a className="link" href={href(t)} target="_blank" rel="noreferrer">{t.fault}{t.only ? ` .${t.only}` : ''}</a>
                        <div className="text-xs text-base-content/50">{t.path}</div></td>
                      <td>{t.surface}<div className="text-xs text-base-content/50">{FAULTS[t.fault]}</div></td>
                      <td className="text-base-content/70">{t.today}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </section>
        ))}
        <section className="space-y-2">
          <h2 className="text-sm font-medium text-base-content/70">Not found and untested (no fault: the empty states)</h2>
          <div className="overflow-hidden rounded-lg border border-base-300">
            <table className="table table-sm">
              <thead><tr><th>Page</th><th>What it shows</th></tr></thead>
              <tbody>
                {NOT_FOUND.map(n => (
                  <tr key={n.path}><td className="whitespace-nowrap"><a className="link" href={at(n.path)} target="_blank" rel="noreferrer">{n.path}</a></td><td>{n.surface}</td></tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      </div>
    </Page>
  )
}
