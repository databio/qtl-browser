import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { KvTable } from '@/components/kv-table'
import { seqcolCollectionUrl } from '@/lib/chrom-sizes'
import { useIndex } from '@/contexts/store-context'

const REPO = 'https://github.com/databio/qtl-browser'

export default function About() {
  const ix = useIndex()
  const st = ix?.store
  return (
    <Page>
      <div className="mx-auto max-w-4xl">
        <PageHeader title="About" className="mb-4 mt-2" />
        <div className="space-y-10">
          <div className="prose prose-sm max-w-none">
            <p>
              This site serves published QTL and GWAS summary statistics from one store. Each collection is a
              publication and the experiments it released; a collection's pages show its experiments, and any GWAS
              in the store can be set beside them.
            </p>
            <h2>Definitions</h2>
            <ul>
              <li><strong>Significance</strong>: each experiment's own rule, as its source tested it; the collection
                page lists them.</li>
              <li><strong>Lead variant</strong>: the variant with the smallest nominal p-value in the phenotype's cis window.</li>
              <li><strong>Credible sets and PIP</strong>: SuSiE 95% credible sets; PIP is the posterior inclusion probability.</li>
              <li><strong>A1 and A2</strong>: variants are stored against the store's reference, so A2 is the reference
                allele and A1 the alternate. Allele frequencies and betas are A1's.</li>
              <li><strong>Colocalization</strong>: coloc.abf (Giambartolomei et al. 2014) runs in your browser over the
                locus in view, from the z statistic and standard error of each trait, with coloc's default priors. It
                assumes one causal variant in the region and uses no LD. A QTL and a GWAS are matched on position and
                alleles.</li>
            </ul>
          </div>

          {st && (
            <KvTable title="This store" rows={[
              { label: 'Name', value: st.name },
              { label: 'Reference', value: st.collection_digest
                ? <ExternalLink icon href={seqcolCollectionUrl(st.collection_digest)}>{st.collection_digest}</ExternalLink> : 'none' },
              { label: 'Collections', value: String(ix.collections.size) },
              { label: 'Experiments', value: String(Object.keys(st.experiments).length) },
              { label: 'Variant catalogs', value: st.variant_catalogs.join(', ') },
              { label: 'Gene annotations', value: st.annotations.join(', ') },
            ]} />
          )}

          <div className="flex flex-wrap gap-x-4 text-sm">
            <ExternalLink className="underline" href={REPO}>Browser source</ExternalLink>
          </div>
        </div>
      </div>
    </Page>
  )
}
