import { Link } from 'react-router'
import Search from '@/components/Search'
import { KvTable, type KvRow } from '@/components/kv-table'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { TableSkeleton } from '@/components/states'
import { GwasTopLoci, StudyTopLoci } from '@/components/TopLoci'
import { useCollection, useCollectionCrumbs, useExperiment, useIndex } from '@/contexts/store-context'
import { experimentLabel, facetRows, shortName } from '@/lib/describe'
import { fmtInt } from '@/lib/format'
import { significanceOf } from '@/lib/significance'
import { phenotypeLabel, type ExperimentDoc } from '@/lib/store'
import { NotInStore } from '@/routes/Collection'

/** An experiment's own page: the landing page of everything scoped to it. A QTL scan gets search
 *  (and, a cis scan, its Top loci); a GWAS its Top loci. Its metadata is two tables: its results, and
 *  the data's facets. */
export default function Experiment() {
  const ix = useIndex()
  const c = useCollection(), e = useExperiment()
  const crumbs = useCollectionCrumbs()
  if (!ix) return <Page><TableSkeleton columns={[{ w: 'w-48' }, { w: 'w-32' }]} rows={4} /></Page>
  if (!c) return <NotInStore />
  if (!e) return <NotInStore what="experiment" />
  return (
    <Page>
      <PageHeader title={`${shortName(c.doc)} · ${experimentLabel(e)}`}
        crumbs={[...crumbs.slice(0, -1), { label: experimentLabel(e) }]} />
      <div className="space-y-10">
        {e.kind !== 'gwas' && (
          <div className="space-y-2">
            <Search hero />
            <p className="text-xs text-base-content/55">
              Search a gene, variant, or region, or <Link className="link-quiet" to="/genes">browse the genes</Link>.
            </p>
          </div>
        )}
        {e.kind === 'qtl_cis' && <StudyTopLoci />}
        {e.kind === 'gwas' && <GwasTopLoci gwas={e} />}
        <div className="grid items-start gap-4 md:grid-cols-2">
          <KvTable title="Results" rows={resultRows(e)} />
          <KvTable title="Metadata" rows={facetRows([e])} />
        </div>
      </div>
    </Page>
  )
}

/** What the experiment holds, one row per fact. */
function resultRows(e: ExperimentDoc): KvRow[] {
  const rows: KvRow[] = [{ label: 'Experiment', value: <span className="font-mono">{e.id}</span> }]
  if (e.kind === 'gwas') {
    const cases = e.source?.n_cases, controls = e.source?.n_controls
    rows.push({ label: 'Variants', value: fmtInt(e.n_rows) })
    if (typeof cases === 'number')
      rows.push({ label: 'Cases / controls', value: `${fmtInt(cases)} / ${typeof controls === 'number' ? fmtInt(controls) : '—'}` })
    return rows
  }
  rows.push({ label: 'Phenotype', value: phenotypeLabel(e.phenotype_type) })
  rows.push({ label: 'Phenotypes', value: fmtInt(e.n_phenotypes) })
  if (e.kind === 'qtl_cis') {
    const sig = e.counts?.significant
    const rule = significanceOf(e).ruleText
    rows.push({ label: 'Significant', value: sig == null ? 'not assessed' : `${fmtInt(sig)}${rule ? ` (${rule})` : ''}` })
  } else {
    rows.push({ label: 'Associations', value: fmtInt(e.n_rows) })
  }
  rows.push({ label: 'Gene annotation', value: e.annotation })
  return rows
}
