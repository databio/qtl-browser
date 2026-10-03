export const dbsnp = (rsid: string) => `https://www.ncbi.nlm.nih.gov/snp/${rsid}`
export const ucsc = (chr: string, start: number, end: number) =>
  `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${chr}%3A${start}-${end}`
export const ensemblGene = (gene_id: string) =>
  `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${gene_id}`
export const gtexGene = (symbol: string) => `https://gtexportal.org/home/gene/${symbol}`
export const openTargetsGene = (gene_id: string) => `https://platform.opentargets.org/target/${gene_id}`
export const geneCards = (symbol: string) => `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${symbol}`
export const PREPRINT = 'https://www.medrxiv.org/content/10.64898/2026.01.12.26343934v1'
export const ZENODO = 'https://zenodo.org/records/21382723'
export const PIPELINE = 'https://github.com/connor122721/nf-eqtls'
