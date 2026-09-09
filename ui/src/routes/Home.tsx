import { Link } from 'react-router'
import Search from '@/components/Search'
import ExternalLink from '@/components/ExternalLink'
import ColocLoci from '@/components/ColocLoci'
import SentenceLines from '@/components/SentenceLines'
import { Page } from '@/components/page'
import { COLOC_EQTL_GENES, COLOC_SQTL_GENES } from '@/lib/coloc'
import { PREPRINT, ZENODO } from '@/lib/links'
import { useManifest } from '@/contexts/manifest-context'
import { fmtInt } from '@/lib/format'

export default function Home() {
  const counts = useManifest()?.counts ?? null

  return (
    <Page>
      <div className="mt-8">
        <h1 className="text-4xl font-extralight tracking-tight">TOPCHeF</h1>
        <SentenceLines className="mt-2.75 text-sm text-base-content/55" sentences={[
          'Expression and splicing QTL mapped in left-ventricle tissue from failing and non-failing human hearts.',
          // the counts come from the manifest fetch; hold the line's height until then so the
          // search bar and track do not move down when it lands
          counts ? `${fmtInt(counts.egenes)} eGenes and ${fmtInt(counts.sqtl_sig_phenotypes)} sQTL introns across ${fmtInt(counts.genes_tested)} tested genes; ${COLOC_EQTL_GENES.length} eGenes and ${COLOC_SQTL_GENES.length} sGenes colocalize with dilated cardiomyopathy (DCM) risk.` : ' ',
          <>
            Search a gene, variant, or region for summary statistics, fine-mapping, and locus views.
            <span className="ml-3 inline-flex gap-x-3">
              <ExternalLink className="underline" href={PREPRINT}>Preprint</ExternalLink>
              <ExternalLink className="underline" href={ZENODO}>Zenodo</ExternalLink>
            </span>
          </>,
        ]} />
        <div className="mt-6"><Search hero autoFocus /></div>
        <p className="mt-1.5 text-xs text-base-content/55">
          Try <Link className="link-quiet" to="/gene/ENSG00000128591">FLNC</Link>,{' '}
          <Link className="link-quiet" to="/gene/ENSG00000157933">SKI</Link>,{' '}
          <Link className="link-quiet" to="/variant/rs2503715">rs2503715</Link>, or{' '}
          <Link className="link-quiet" to="/region/chr10:73000000-74500000">chr10:73,000,000-74,500,000</Link>.
        </p>
      </div>

      <div className="mt-10 pb-4">
        <ColocLoci />
      </div>
    </Page>
  )
}
