/** Every SQL string in one place. Table paths match data/derived/manifest.json. */
import { lit, one, parquet, rows, type Row } from './db'

export interface SearchHit extends Row {
  gene_id: string; symbol: string | null; chr: string; tss: number
  tested: boolean; is_egene: boolean | null; n_sqtl_sig: number
  /** partition file of the gene's nominal rows and gene_detail row (null when not tested) */
  bin: number | null
  start: number; end: number; strand: string; biotype: string
}

/** Path of the nominal cis file holding this gene: chromosome and TSS-rank bin. */
export const nominalFile = (hit: SearchHit, qtlType: 'e' | 's') =>
  `${qtlType === 'e' ? 'cis_eqtl_nominal' : 'cis_sqtl_nominal'}/chr=${hit.chr}/bin=${hit.bin}/data.parquet`

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
  num_var: number; lead_position: number; lead_A1: string; lead_A2: string; lead_rsid: string | null
  lead_af: number; lead_tss_distance: number; slope: number; slope_se: number
  pval_nominal: number; pval_perm: number; pval_beta: number; qval: number; is_sqtl: boolean
  n_credible_sets: number
}

export interface TransRow extends Row {
  qtl_type: string; phenotype_id: string; gene_id: string; symbol: string | null
  gene_chr: string; gene_tss: number; variant_chr: string; position: number; rsid: string | null; af: number
  pval: number; beta: number; beta_se: number; r2: number
}

export const searchGenes = (q: string, limit = 12) =>
  rows<SearchHit>(`
    SELECT * FROM search_index
    WHERE upper(symbol) LIKE ${lit(q.toUpperCase() + '%')} OR upper(gene_id) LIKE ${lit(q.toUpperCase() + '%')}
    ORDER BY tested DESC, is_egene DESC NULLS LAST, length(symbol), symbol LIMIT ${limit}`)

export const searchBySymbols = (symbols: string[]) =>
  rows<SearchHit>(`SELECT * FROM search_index WHERE symbol IN (${symbols.map(lit).join(',')})`)

export const resolveGene = (id: string) =>
  one<SearchHit>(`SELECT * FROM search_index WHERE gene_id = ${lit(id)} OR upper(symbol) = ${lit(id.toUpperCase())}
                  ORDER BY tested DESC LIMIT 1`)

export interface Exon extends Row { start: number; end: number }
export interface GeneDetail { gene: Gene; exons: Exon[]; splice: SplicePhenotype[] }

/** Everything the gene page needs besides the locus, in one row group of one small file:
 *  the genes row, the collapsed exon model, and every tested intron. The list columns come
 *  back as JSON text so one read serves all three. */
export const geneDetail = async (hit: SearchHit): Promise<GeneDetail | null> => {
  const r = await one<Gene & { exons_json: string | null; splice_json: string | null }>(`
    SELECT * EXCLUDE (exons, splice), to_json(exons) AS exons_json, to_json(splice) AS splice_json
    FROM ${parquet(`gene_detail/chr=${hit.chr}/bin=${hit.bin}/data.parquet`)} WHERE gene_id = ${lit(hit.gene_id)}`)
  if (!r) return null
  const { exons_json, splice_json, ...gene } = r
  const splice = (JSON.parse(splice_json ?? '[]') as SplicePhenotype[]).map(p => ({ ...p, gene_id: hit.gene_id, symbol: hit.symbol, chr: hit.chr, tss: hit.tss }))
  return { gene: gene as Gene, exons: JSON.parse(exons_json ?? '[]') as Exon[], splice }
}

// ---- paged tables over materialized windows -------------------------------------------------
// The gene page materializes two in-memory tables per gene (see db.ts `materialize`): the cis
// window that the locus plot draws from, and the gene's trans rows. The cis and trans tables
// page off those with limit/offset, a count, and an unpaged export, so every interaction is a
// local query and nothing is re-read over HTTP.

/** `CREATE TABLE ... AS` body for one gene's trans rows (one or two row groups of the
 *  per-chromosome file; 65k rows for the busiest gene). */
export const transSQL = (hit: SearchHit) =>
  `SELECT * FROM ${parquet(`trans_pairs/chr=${hit.chr}/data.parquet`)} WHERE gene_id = ${lit(hit.gene_id)}`

/** Same for one variant's trans rows, from the variant-keyed copy (position-sorted, so one
 *  position is one or two row groups). */
export const transAtSQL = (chr: string, pos: number) =>
  `SELECT * FROM ${parquet(`trans_by_variant/chr=${chr}/data.parquet`)} WHERE position = ${pos}`

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

export const genesInRegion = (chr: string, start: number, end: number) =>
  rows<SearchHit>(`SELECT * FROM search_index WHERE chr = ${lit(chr)} AND tss BETWEEN ${start} AND ${end} ORDER BY tss`)

// ---- variant page ---------------------------------------------------------------------------

/** `in_cis` false marks a variant seen only in the genome-wide trans scan, outside every cis
 *  window. Alleles are null when the variant appears only in the trans eQTL file, which has no
 *  allele columns; such rows match dbSNP by position only. */
export interface VariantRow extends Row {
  chr: string; position: number; A1: string | null; A2: string | null; rsid: string | null; rs_number: number | null
  match: string; in_cis: boolean
}

export const variantByRsid = (rsNumber: number) =>
  rows<VariantRow>(`SELECT * FROM ${parquet('variants_by_rsid.parquet')} WHERE rs_number = ${rsNumber} ORDER BY chr, position, A1`)

export const variantByPosition = (chr: string, pos: number) =>
  rows<VariantRow>(`SELECT * FROM ${parquet(`variants_by_position/chr=${chr}/data.parquet`)} WHERE position = ${pos} ORDER BY A1`)

/** Genes whose lead eQTL variant is this position. TSS is within 1 Mb by construction, which
 *  is what lets DuckDB prune genes.parquet row groups on (chr, tss). */
export const leadGenesAt = (chr: string, pos: number) =>
  rows<Gene>(`SELECT * FROM ${parquet('genes.parquet')}
    WHERE chr = ${lit(chr)} AND tss BETWEEN ${pos - 1_000_000} AND ${pos + 1_000_000} AND lead_position = ${pos} ORDER BY pval_perm`)

export const leadPhenotypesAt = (chr: string, pos: number) =>
  rows<SplicePhenotype>(`SELECT * FROM ${parquet('splice_phenotypes.parquet')}
    WHERE chr = ${lit(chr)} AND tss BETWEEN ${pos - 1_000_000} AND ${pos + 1_000_000} AND lead_position = ${pos} ORDER BY pval_perm`)

export interface CredibleSetHit extends CredibleSetRow { gene_id: string; symbol: string | null; tss: number }
export const credibleSetsAt = (chr: string, pos: number) =>
  rows<CredibleSetHit>(`SELECT * FROM ${parquet('credible_sets.parquet')}
    WHERE chr = ${lit(chr)} AND tss BETWEEN ${pos - 1_000_000} AND ${pos + 1_000_000} AND position = ${pos} ORDER BY pip DESC`)


export interface CisHit extends Row {
  gene_id: string; symbol: string | null; phenotype_id?: string; tss_distance: number
  pval_nominal: number; slope: number; slope_se: number; af: number; pip: number | null; cs_id: number | null
}
/** Every gene (or splice phenotype) whose cis window covers this position: one nominal row
 *  each. Touches every row group whose position range spans `pos` (~30 genes' worth), so it
 *  runs on demand, not on page load. */
export const cisHitsAt = async (chr: string, pos: number, qtlType: 'e' | 's') => {
  // the genes whose cis window covers this position are the tested genes with a TSS within
  // 1 Mb; their bins (from the in-memory index) name the nominal files to scan
  const table = qtlType === 'e' ? 'cis_eqtl_nominal' : 'cis_sqtl_nominal'
  const bins = await rows<{ bin: number }>(`SELECT DISTINCT bin FROM search_index
    WHERE chr = ${lit(chr)} AND bin IS NOT NULL AND tss BETWEEN ${pos - 1_000_000} AND ${pos + 1_000_000}
      ${qtlType === 's' ? 'AND n_sqtl_sig > 0' : ''} ORDER BY bin`)
  if (!bins.length) return [] as CisHit[]
  const files = bins.map(b => `'${table}/chr=${chr}/bin=${b.bin}/data.parquet'`).join(', ')
  return rows<CisHit>(`
    SELECT n.gene_id, s.symbol, ${qtlType === 's' ? 'n.phenotype_id,' : ''} n.tss_distance, n.pval_nominal, n.slope, n.slope_se, n.af, n.pip, n.cs_id
    FROM read_parquet([${files}], hive_partitioning=false) n
    LEFT JOIN search_index s USING (gene_id)
    WHERE n.position = ${pos}
    ORDER BY n.pval_nominal`)
}

// ---- colocalized loci on the landing page ---------------------------------------------------

/** Full gene rows for a handful of known genes: one OR term per (chr, tss, gene_id) so row
 *  groups prune on chr/tss statistics instead of reading the whole file. */
export const geneRowsFor = (hits: SearchHit[]) =>
  hits.length === 0 ? Promise.resolve([] as Gene[]) : rows<Gene>(`SELECT * FROM ${parquet('genes.parquet')} WHERE ${
    hits.map(h => `(chr = ${lit(h.chr)} AND tss = ${h.tss} AND gene_id = ${lit(h.gene_id)})`).join(' OR ')}`)

// ---- gene track under the locus plot --------------------------------------------------------

export interface WindowGene extends Row { gene_id: string; symbol: string | null; start: number; end: number; strand: string; tss: number; biotype: string }
/** Genes overlapping a window, from the in-memory search index (no fetch). */
export const genesInWindow = (chr: string, lo: number, hi: number) =>
  rows<WindowGene>(`SELECT gene_id, symbol, start, "end", strand, tss, biotype FROM search_index
    WHERE chr = ${lit(chr)} AND "end" >= ${lo} AND start <= ${hi} ORDER BY start`)
