import ExternalLink from '@/components/ExternalLink'
import { Page } from '@/components/page'
import { PageHeader } from '@/components/page-header'
import { KvTable } from '@/components/kv-table'
import { PIPELINE, PREPRINT, ZENODO } from '@/lib/links'
import { useManifest } from '@/contexts/manifest-context'
import { roundingFacts } from '@/lib/rounding'

const GWAS_PAPER = 'https://doi.org/10.1038/s41588-024-01975-5'
const CVDKP = 'https://kp4cd.org/dataset_downloads/mi'
const SEQCOL = 'https://seqcolapi.databio.org'
const REPO = 'https://github.com/databio/qtl-browser'

export default function About() {
  const m = useManifest()
  const sources = m?.sources ?? {}
  const counts = m?.counts ?? {}
  const gwas = m?.gwas_dcm ?? null
  const rounding = roundingFacts(m)
  const gwasSet = gwas?.file.includes('BiobanksOnly') ? 'biobank-only meta-analysis' : gwas?.file.includes('MTAG') ? 'MTAG analysis' : 'meta-analysis'
  return (
    <Page>
      <div className="mx-auto max-w-4xl">
        {/* mb-4: the same 16 px the prose puts under its own h2 headings below */}
        <PageHeader title="About" className="mb-4 mt-2" />
        <div className="space-y-10">
          <div className="prose prose-sm max-w-none">
            <p>
              TOPCHeF (Trans-Omics for Precision Medicine in Congestive Heart Failure) paired whole-genome and RNA sequencing
              from left-ventricle tissue of dilated cardiomyopathy, ischemic cardiomyopathy, and non-failing donors, and mapped
              cis and trans expression (eQTL) and splicing (sQTL) quantitative trait loci with TensorQTL and SuSiE. Methods,
              sample counts, and the colocalization with dilated cardiomyopathy risk are in the{' '}
              <ExternalLink href={PREPRINT}>preprint</ExternalLink>. This site serves the published summary statistics.
            </p>
            <h2>Definitions</h2>
            <ul>
              <li><strong>eGene, sQTL intron</strong>: permutation p-value below 0.05. An sGene has at least one significant intron. A Benjamini-Hochberg q-value on the beta-approximated permutation p is listed alongside.</li>
              <li><strong>Lead variant</strong>: the variant with the smallest nominal p-value in the cis window, ±1 Mb of the transcription start site.</li>
              <li><strong>Credible sets and PIP</strong>: SuSiE 95% credible sets; PIP is the posterior inclusion probability. A variant in two sets of one phenotype is listed under both in the credible-set table; the locus plot and cis table show its higher-PIP membership.</li>
              <li><strong>A1 and A2</strong>: A1 is the effect allele, the minor allele in TOPCHeF; A2 is the reference allele. Slopes are in standard-deviation units of the phenotype per A1 allele.</li>
              {rounding && <li><strong>Rounded values</strong>: per-variant p-values, slopes, standard errors, and allele
                frequencies are stored in compressed form, so a p-value shown here can differ from the source by up to{' '}
                {rounding.pPct}% and a slope by up to {rounding.slopeSe} standard errors. Gene-level results and the DCM
                GWAS values are exact, and <ExternalLink href={ZENODO}>exact per-variant values are on Zenodo</ExternalLink>.</li>}
              <li><strong>Splice phenotypes</strong>: leafcutter intron excision ratios, shown as intron coordinates and strand. Introns sharing a splice site share a cluster. Every tested intron has its permutation result and its per-variant nominal statistics.</li>
              <li><strong>Colocalized loci</strong>: the 21 eQTL and 4 sQTL genes with coloc PP.H4 above 0.8 against the DCM GWAS. PJVK and CDKN1A are not eGenes by the permutation rule; their colocalization used nominal statistics.</li>
            </ul>
            <h2>Coordinates and identifiers</h2>
            <ul>
              <li>Coordinates are GRCh38. Genes follow GENCODE v34; gene models on the locus plot are the union of each gene's transcript exons.</li>
              <li>rsIDs are assigned from dbSNP by position and alleles. Variants with no dbSNP record are shown as chr:position.</li>
              <li>Chromosome lengths for the genome track come from the <ExternalLink href={SEQCOL}>seqcol</ExternalLink> GRCh38 reference.</li>
            </ul>
            <h2>DCM GWAS comparison</h2>
            <p>
              The landing track and the QTL-versus-GWAS panel use the dilated cardiomyopathy {gwasSet} of{' '}
              <ExternalLink href={GWAS_PAPER}>Jurgens et al. 2024</ExternalLink>
              {gwas && <> ({gwas.n_cases.toLocaleString()} cases, {gwas.n_controls.toLocaleString()} controls)</>} from the{' '}
              <ExternalLink href={CVDKP}>Cardiovascular Disease Knowledge Portal</ExternalLink>. Variants are matched on GRCh38
              position and alleles in either orientation, and the GWAS effect is signed to the QTL effect allele. The landing
              track shows the strongest GWAS p-value per 5 Mb window, red where the window holds a genome-wide significant
              variant; the gene page panel plots every shared variant in the cis window.
            </p>
          </div>

          {m && (
            <KvTable title="Counts" align="right" rows={[
              { label: 'Genes tested', value: counts.genes_tested?.toLocaleString() },
              { label: 'eGenes', value: counts.egenes?.toLocaleString() },
              { label: 'Splice phenotypes tested', value: counts.splice_phenotypes_tested?.toLocaleString() },
              { label: 'Significant sQTL introns', value: `${counts.sqtl_sig_phenotypes?.toLocaleString()} in ${counts.sqtl_sig_genes?.toLocaleString()} genes` },
              { label: 'Variants in cis windows', value: counts.variants_cis?.toLocaleString() },
              { label: 'Variants seen only in trans', value: (counts.variants_trans_only ?? 0).toLocaleString() },
              { label: 'DCM GWAS variants', value: counts.gwas_variants?.toLocaleString() },
              { label: 'trans pairs', value: counts.trans_pairs?.toLocaleString() },
            ]} />
          )}

          {/* file sizes and counts stay in manifest.json next to the data; only the build date is shown */}
          <KvTable
            title={<>Data versions{m?.built && <span className="ml-1.5 font-normal normal-case tracking-normal text-base-content/50">
              (updated {String(m.built).slice(0, 10)})</span>}</>}
            rows={Object.entries(sources).map(([k, v]) => ({ label: k, value: <span><span className="font-medium text-base-content">{v.version}</span> · {v.description}</span> }))} />

          <div className="flex flex-wrap gap-x-4 text-sm">
            <ExternalLink className="underline" href={PREPRINT}>Preprint</ExternalLink>
            <ExternalLink className="underline" href={ZENODO}>Summary statistics on Zenodo</ExternalLink>
            <ExternalLink className="underline" href={PIPELINE}>QTL mapping pipeline (nf-eqtls)</ExternalLink>
            <ExternalLink className="underline" href={REPO}>Browser source</ExternalLink>
          </div>
        </div>
      </div>
    </Page>
  )
}
