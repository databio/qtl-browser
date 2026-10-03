/** The row shapes the pages share, and the SQL over the tables a page materializes (cis windows,
 *  trans rows). Gene lookups and lists are plain JS over the per-chromosome objects (gene-index.ts),
 *  re-exported here. */
import { lit, one, rows, type Row } from './db'
export { genesInRegion, genesInWindow, resolveGene, searchGenes } from './gene-index'

/** A gene as the pages list it (gene-index.ts `hitOf`): an annotated gene on a catalog chromosome,
 *  joined to its eQTL phenotype. `tss`, bounds, symbol and biotype come from the annotation, never
 *  from the study. */
export interface SearchHit extends Row {
  gene_id: string; symbol: string | null; chr: string; tss: number
  start: number; end: number; strand: string; biotype: string; gene_version: number | null
  /** the gene's eQTL phenotype: its search-index `ord`, block in the results file, run, and the
   *  catalog bytes of the pages covering the run (all null when the gene has no eQTL phenotype) */
  ord: number | null; blk_off: number | null; blk_len: number | null; var_start: number | null; n_var: number | null
  var_off: number | null; var_len: number | null; w_lo: number | null; w_hi: number | null
  tested: boolean; is_egene: boolean | null; n_sqtl_sig: number; n_sqtl: number
  /** any phenotype of the gene (eQTL or an intron) has a block: false means not tested at all */
  has_results: boolean
}

export interface Gene extends Row {
  gene_id: string; gene_id_version: string; symbol: string | null; chr: string
  start: number; end: number; strand: string; tss: number; biotype: string; tested: boolean
  num_var: number | null; lead_position: number | null; lead_A1: string | null; lead_A2: string | null
  lead_rsid: string | null; lead_af: number | null; lead_tss_distance: number | null
  slope: number | null; slope_se: number | null; pval_nominal: number | null; pval_perm: number | null
  pval_beta: number | null; qval: number | null; is_egene: boolean | null
  n_credible_sets: number; n_trans_pairs: number
}

export interface CredibleSetRow extends Row {
  qtl_type: string; phenotype_id: string; chr: string; position: number; A1: string; A2: string
  rsid: string | null; af: number; cs_id: number; pip: number
}

export interface CisRow extends Row {
  position: number; A1: string; A2: string; rs_number: number | null; tss_distance: number
  af: number; ma_samples: number; ma_count: number; pval_nominal: number; slope: number
  slope_se: number; pip: number | null; cs_id: number | null; phenotype_id?: string
}

export interface SplicePhenotype extends Row {
  phenotype_id: string; gene_id: string; symbol: string | null; chr: string
  intron_start: number; intron_end: number; cluster_id: string; strand: string; tss: number
  num_var: number; lead_position: number | null; lead_A1: string | null; lead_A2: string | null; lead_rsid: string | null
  lead_af: number | null; lead_tss_distance: number | null; slope: number | null; slope_se: number | null
  pval_nominal: number | null; pval_perm: number | null; pval_beta: number | null
  /** not stored in v1 (null) */
  qval: number | null; is_sqtl: boolean
  n_credible_sets: number
  /** the intron's block in the sQTL results file */
  blk_off: number; blk_len: number
}

export interface TransRow extends Row {
  qtl_type: string; phenotype_id: string; gene_id: string; symbol: string | null
  gene_chr: string; gene_tss: number; variant_chr: string; position: number; rsid: string | null; af: number
  pval: number; beta: number; beta_se: number; r2: number
}

export interface Exon extends Row { start: number; end: number }
/** Everything the gene page needs besides the locus and the introns (gene.ts): the genes row
 *  (annotation plus the eQTL block's details and lead row) and the collapsed exon model. */
export interface GeneDetail { gene: Gene; exons: Exon[] }

// ---- paged tables over materialized windows -------------------------------------------------
// The gene page holds the cis window that the locus plot draws from (gene.ts `locusTable`) as an
// in-memory table. The cis table pages off it with limit/offset, a count, and an unpaged export,
// so every interaction is a local query and nothing is re-read over HTTP. The trans queries below
// serve TransTable, which has no v1 data yet (SPEC section 10); a trans reader will build its
// table with the `TransRow` columns.

interface PagedQuery {
  table: string                 // materialized table to page off
  maxP?: number                 // filter
  search?: string               // rsID or position, prefix match
  orderBy?: string
  desc?: boolean
  limit?: number
  offset?: number
}

export interface CisQuery extends PagedQuery {
  chr: string
  qtlType: 'e' | 's'
  phenotypeId?: string          // sQTL only, for the CSV
}

/** rsID / position prefix search on a window; rsIDs are integer rs_number in the cis window
 *  and text rsid in the trans table. */
function searchWhere(s: string | undefined, rsExpr: string): string | null {
  const t = s?.trim()
  if (!t) return null
  const rs = /^rs(\d+)$/i.exec(t)
  const pos = /^(?:chr[0-9xy]+:)?([\d,]+)$/i.exec(t)
  if (rs) return `${rsExpr} LIKE ${lit(rs[1] + '%')}`
  if (pos) return `CAST(position AS VARCHAR) LIKE ${lit(pos[1].replace(/,/g, '') + '%')}`
  return 'false'
}

function cisWhere(q: CisQuery): string {
  const parts = ['true']
  if (q.maxP != null) parts.push(`pval_nominal <= ${q.maxP}`)
  const s = searchWhere(q.search, 'CAST(rs_number AS VARCHAR)')
  if (s) parts.push(s)
  return `FROM ${q.table} WHERE ${parts.join(' AND ')}`
}

const CIS_COLS = 'position, A1, A2, rs_number, tss_distance, af, ma_samples, ma_count, pval_nominal, slope, slope_se, pip, cs_id'
const CIS_SORTABLE = new Set(['position', 'pval_nominal', 'slope', 'af', 'pip', 'tss_distance', 'ma_count'])

export const cisRows = (q: CisQuery) => {
  const col = q.orderBy && CIS_SORTABLE.has(q.orderBy) ? q.orderBy : 'pval_nominal'
  return rows<CisRow>(`SELECT ${CIS_COLS} ${cisWhere(q)} ORDER BY ${col} ${q.desc ? 'DESC' : 'ASC'} NULLS LAST
                       LIMIT ${q.limit ?? 50} OFFSET ${q.offset ?? 0}`)
}

export const cisCount = async (q: CisQuery) =>
  Number((await one<{ n: number }>(`SELECT count(*) AS n ${cisWhere(q)}`))?.n ?? 0)

/** Full cis window for CSV. */
export const cisAll = (q: CisQuery) =>
  rows<CisRow>(`SELECT ${CIS_COLS} ${cisWhere({ ...q, maxP: undefined, search: undefined })} ORDER BY position`)

/** A trans table is keyed by the gene (gene page: rows are variants, one QTL type per tab) or by
 *  the variant (variant page: rows are genes, both types together). The key decides what the
 *  search box matches: variant rsID or position, or gene symbol or Ensembl ID. */
export interface TransQuery extends PagedQuery { qtlType?: 'e' | 's'; keyedBy?: 'gene' | 'variant' }

function geneSearchWhere(s: string | undefined): string | null {
  const t = s?.trim()
  if (!t) return null
  return `(symbol ILIKE ${lit(t + '%')} OR gene_id ILIKE ${lit(t + '%')})`
}

function transWhere(q: TransQuery): string {
  const parts = ['true']
  if (q.qtlType) parts.push(`qtl_type = ${lit(q.qtlType)}`)
  if (q.maxP != null) parts.push(`pval <= ${q.maxP}`)
  const s = q.keyedBy === 'variant' ? geneSearchWhere(q.search) : searchWhere(q.search, 'substr(rsid, 3)')
  if (s) parts.push(s)
  return `FROM ${q.table} WHERE ${parts.join(' AND ')}`
}

// chromosome order then position, so chr2 sorts before chr10
const chromOrder = (col: string) => `CASE WHEN ${col} = 'chrX' THEN 23 WHEN ${col} = 'chrY' THEN 24 ELSE TRY_CAST(substr(${col}, 4) AS INTEGER) END`
const TRANS_SORTABLE: Record<string, string[]> = {
  position: [chromOrder('variant_chr'), 'position'],
  gene: [chromOrder('gene_chr'), 'gene_tss'],
  type: ['qtl_type', 'pval'],
  af: ['af'], pval: ['pval'], beta: ['beta'], r2: ['r2'],
}

function transOrder(q: TransQuery): string {
  const cols = (q.orderBy ? TRANS_SORTABLE[q.orderBy] : undefined) ?? ['pval']
  return cols.map(c => `${c} ${q.desc ? 'DESC' : 'ASC'} NULLS LAST`).join(', ')
}

export const transRows = (q: TransQuery) =>
  rows<TransRow>(`SELECT * ${transWhere(q)} ORDER BY ${transOrder(q)} LIMIT ${q.limit ?? 50} OFFSET ${q.offset ?? 0}`)

export const transCount = async (q: TransQuery) =>
  Number((await one<{ n: number }>(`SELECT count(*) AS n ${transWhere(q)}`))?.n ?? 0)

/** Filtered trans rows for CSV, in the table's current order. */
export const transAll = (q: TransQuery) =>
  rows<TransRow>(`SELECT * ${transWhere(q)} ORDER BY ${transOrder(q)}`)

// ---- gene track under the locus plot --------------------------------------------------------

export interface WindowGene extends Row { gene_id: string; symbol: string | null; start: number; end: number; strand: string; tss: number; biotype: string }
