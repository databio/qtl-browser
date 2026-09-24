/**
 * Genes colocalized with the Jurgens et al. 2024 DCM GWAS (PP.H4 > 0.8) in the preprint.
 * Table 2 (eQTL) and Supplemental Table 2 (sQTL) bodies are not in the public PDF; this list
 * comes from the authors. Replace with the authors' coloc table once they share it.
 */
export const COLOC_EQTL_GENES = [
  'VPREB3', 'SYNPO2L', 'SQLE', 'SMARCB1', 'SKI', 'PROM1', 'PRKCA', 'PJVK', 'MYOZ1', 'MTSS1',
  'MMP11', 'MAP3K7CL', 'LMF1', 'LINC00964', 'FLNC', 'CRIM1-DT', 'CRIM1', 'CDKN1A', 'CAMK2D',
  'ADAMTS7P3', 'ACTN2',
]
export const COLOC_SQTL_GENES = ['SYNPO2L', 'CAMK2D', 'TKT', 'LMF1']

/** The annotation `COLOC_LOCI` was read from (`annotations/gencode_v34.json` `identity_digest`). */
export const COLOC_ANNOTATION = 'vQpqU_qkV212NyRfCdED01tLMiTnNiKR'

/**
 * Where each coloc gene sits in that annotation: every gene carrying the symbol, with its gene id,
 * chromosome and TSS. The Home track places its markers from this table without a request; with any
 * other annotation it looks the symbols up in the store instead (ColocLoci.tsx), and
 * `npm run store-check` fails when a store with this annotation disagrees with the table.
 */
export const COLOC_LOCI: Record<string, { gene_id: string; chr: string; tss: number }[]> = {
  VPREB3: [{ gene_id: 'ENSG00000128218', chr: 'chr22', tss: 23754425 }],
  SYNPO2L: [{ gene_id: 'ENSG00000166317', chr: 'chr10', tss: 73663803 }],
  SQLE: [{ gene_id: 'ENSG00000104549', chr: 'chr8', tss: 124998497 }],
  SMARCB1: [{ gene_id: 'ENSG00000099956', chr: 'chr22', tss: 23786931 }],
  SKI: [{ gene_id: 'ENSG00000157933', chr: 'chr1', tss: 2228319 }],
  PROM1: [{ gene_id: 'ENSG00000007062', chr: 'chr4', tss: 16084378 }],
  PRKCA: [{ gene_id: 'ENSG00000154229', chr: 'chr17', tss: 66302613 }],
  PJVK: [{ gene_id: 'ENSG00000204311', chr: 'chr2', tss: 178451346 }],
  MYOZ1: [{ gene_id: 'ENSG00000177791', chr: 'chr10', tss: 73641474 }],
  MTSS1: [{ gene_id: 'ENSG00000170873', chr: 'chr8', tss: 124728429 }],
  MMP11: [{ gene_id: 'ENSG00000099953', chr: 'chr22', tss: 23768226 }],
  MAP3K7CL: [{ gene_id: 'ENSG00000156265', chr: 'chr21', tss: 29077471 }],
  LMF1: [{ gene_id: 'ENSG00000103227', chr: 'chr16', tss: 981318 }],
  LINC00964: [{ gene_id: 'ENSG00000249816', chr: 'chr8', tss: 124823702 }],
  FLNC: [{ gene_id: 'ENSG00000128591', chr: 'chr7', tss: 128830377 }],
  'CRIM1-DT': [{ gene_id: 'ENSG00000260025', chr: 'chr2', tss: 36355114 }],
  CRIM1: [{ gene_id: 'ENSG00000150938', chr: 'chr2', tss: 36355778 }],
  CDKN1A: [{ gene_id: 'ENSG00000124762', chr: 'chr6', tss: 36676460 }],
  CAMK2D: [{ gene_id: 'ENSG00000145349', chr: 'chr4', tss: 113761927 }],
  ADAMTS7P3: [{ gene_id: 'ENSG00000261143', chr: 'chr15', tss: 77976042 }],
  ACTN2: [{ gene_id: 'ENSG00000077522', chr: 'chr1', tss: 236664141 }],
  TKT: [{ gene_id: 'ENSG00000163931', chr: 'chr3', tss: 53256052 }],
}
