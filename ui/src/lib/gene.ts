/**
 * Every gene-page read, from the qtlstore (SPEC.md sections 8 and 10). Opening a gene starts plain
 * Range requests, with byte offsets from the search index: the gene's eQTL block (its rows,
 * credible sets and details), one variants range covering the runs of every phenotype of the gene
 * (its eQTL run and every intron's), and the GWAS rows over that window (through the GWAS index,
 * one whole-object read per session). The sQTL tab adds one request, a span holding all of the
 * gene's intron blocks, which sit next to each other in the results file. The gene's phenotypes
 * and exon model come from its chromosome's search index part and exon models (gene-index.ts), which
 * the page resolving the gene has already read, so none of this waits for the query engine.
 *
 * Annotation (symbol, TSS, bounds, exons) is joined on `gene_id`; the study's blocks carry none.
 * Rows reach DuckDB as Arrow IPC streams (insertArrow). Decoding lives in store-decode.ts.
 */
import { makeVector, Table, tableToIPC, Utf8, vectorFromArray } from 'apache-arrow'
import { getDB, insertArrow, tableName } from './db'
import { exonModel, phenotypesOfGene, type PhenotypeRow } from './gene-index'
import { EQTL_TYPE, fetchDecoded, fetchObject, getStore, gwasFile, resultsFile, SQTL_TYPE, variantsFile } from './store'
import { csMembers, decodeGwasIndex, decodeGwasRange, decodeResultBlock, decodeVariantRange, gwasRange, readerColumns, sliceRun,
  type GwasColumns, type GwasIndex, type ResultBlock, type VariantRange } from './store-decode'
import { dropTable } from './db'
import type { CredibleSetRow, Gene, GeneDetail, SearchHit, SplicePhenotype } from './queries'

/** A phenotype with a block (every phenotype on a chromosome part has one; only trans-only ones do not). */
type Placed = PhenotypeRow & { blk_off: number; blk_len: number }

export interface GenePack {
  hit: SearchHit
  /** request 1: the gene's eQTL block (null when the gene has no eQTL phenotype) */
  block: Promise<ResultBlock | null>
  /** request 2: the variants range covering the runs of all the gene's phenotypes (null when none has rows) */
  variants: Promise<VariantRange | null>
  /** request 3: the GWAS rows over the gene's window as a DuckDB table (empty without a GWAS or rows
   *  there); dropped when the gene leaves the cache */
  gwas: Promise<string>
  /** the genes row and exon model: no request beyond the two above and the chromosome's exon models */
  detail: Promise<GeneDetail>
  /** the sQTL tab: +1 request for the span of the gene's intron blocks, memoized */
  splice(): Promise<SplicePhenotype[]>
  /** one intron's block, from that span */
  intron(phenotypeId: string): Promise<ResultBlock>
}

// ---- the GWAS window --------------------------------------------------------------------------

let gwasIndex: Promise<GwasIndex | null> | null = null
/** The experiment's GWAS index: one whole-object read per session, null without a GWAS. */
export function loadGwasIndex(): Promise<GwasIndex | null> {
  if (!gwasIndex) {
    const p = getStore().then(async s => {
      const g = s.experiment.gwas
      if (!g) return null
      const idx = decodeGwasIndex(await fetchObject(g.index), s.catalog.collection_digest, g.index)
      if (idx.blockRows !== g.block_rows || idx.nValues.join() !== g.n_values.join())
        throw new Error(`${g.index}: block rows or n table differ from the experiment's gwas entry`)
      return idx
    })
    gwasIndex = p
    p.catch(() => { if (gwasIndex === p) gwasIndex = null })
  }
  return gwasIndex
}

/** The column names and types the locus SQL reads from a GWAS window. `ea` is ALT (the effect
 *  allele after orientation), `nea` REF, and beta and eaf describe ALT (SPEC section 10). */
const GWAS_DDL = 'position INTEGER, ea VARCHAR, nea VARCHAR, beta FLOAT, se FLOAT, p DOUBLE, eaf FLOAT, rsid VARCHAR, n INTEGER'

/** The GWAS rows in [lo, hi] on the gene's chromosome, decoded once and inserted as a table. */
async function gwasTable(hit: SearchHit, lo: number | null, hi: number | null): Promise<string> {
  const detail = { gene_id: hit.gene_id }
  const [s, index] = await Promise.all([getStore(), loadGwasIndex()])
  const f = gwasFile(s, hit.chr)
  if (index && !f !== !index.chroms.has(hit.chr)) throw new Error(`${hit.chr}: the GWAS files and the GWAS index disagree`)
  let cols: GwasColumns | null = null
  if (index && f && lo != null && hi != null) {
    const range = gwasRange(index, hit.chr, lo, hi)
    if (range) cols = await fetchDecoded('store:gwas', detail, f.name, range.off, range.len,
      (b, what) => decodeGwasRange(b, index, hit.chr, range, lo, hi, what))
  }
  const ipc = cols?.rows ? tableToIPC(new Table({
    position: makeVector(cols.position), ea: vectorFromArray(cols.alt, new Utf8()), nea: vectorFromArray(cols.ref, new Utf8()),
    beta: makeVector(Float32Array.from(cols.beta)), se: makeVector(Float32Array.from(cols.se)), p: makeVector(cols.p),
    eaf: makeVector(Float32Array.from(cols.af)),
    rsid: vectorFromArray(Array.from(cols.rsNumber, x => (x ? `rs${x}` : null)), new Utf8()),
    n: makeVector(cols.n),
  }), 'stream') : null
  const name = tableName('gwas')
  await insertArrow(name, ipc ?? GWAS_DDL, { ...detail, part: 'gwas' })
  return name
}

// ---- rows of a block --------------------------------------------------------------------------

/** The block's row at its phenotype's group lead (matched on position and alleles), with the
 *  variant fields it needs, or null when the lead is not one of the block's rows. */
function leadRow(block: ResultBlock | null, range: VariantRange | null) {
  const lead = block?.details.group?.lead
  if (!block || !lead || block.varStart == null || !range) return null
  const run = sliceRun(range, block.varStart, block.nRows, block)
  for (let i = 0; i < run.n; i++) {
    if (run.position[i] !== lead.pos || run.ref[i] !== lead.ref || run.alt[i] !== lead.alt) continue
    const nan = (x: number) => (Number.isNaN(x) ? null : x)
    return { pval: nan(block.pval[i]), slope: nan(block.slope[i]), se: nan(block.se[i]), af: nan(run.af[i]),
      rsid: run.rsNumber[i] ? `rs${run.rsNumber[i]}` : null }
  }
  return null
}

/** The genes row the gene page prints: the annotation's fields, the eQTL block's details (group
 *  lead, permutation p, credible sets) and the lead's own nominal row. No q-value is stored in v1. */
function geneRow(hit: SearchHit, block: ResultBlock | null, range: VariantRange | null): Gene {
  const d = block?.details, g = d?.group ?? null, lead = g?.lead ?? null, row = leadRow(block, range)
  return {
    gene_id: hit.gene_id, gene_id_version: hit.gene_version != null ? `${hit.gene_id}.${hit.gene_version}` : hit.gene_id,
    symbol: hit.symbol, chr: hit.chr, start: hit.start, end: hit.end, strand: hit.strand, tss: hit.tss, biotype: hit.biotype,
    tested: hit.tested, num_var: g?.n_variants ?? block?.nRows ?? null,
    lead_position: lead?.pos ?? null, lead_A1: lead?.alt ?? null, lead_A2: lead?.ref ?? null,
    lead_rsid: row?.rsid ?? null, lead_af: row?.af ?? null, lead_tss_distance: lead ? lead.pos - hit.tss : null,
    slope: row?.slope ?? null, slope_se: row?.se ?? null, pval_nominal: row?.pval ?? null,
    pval_perm: g?.p_perm ?? null, pval_beta: g?.p_beta ?? null, qval: null, is_egene: hit.is_egene,
    n_credible_sets: d?.n_credible_sets ?? 0, n_trans_pairs: 0,
  }
}

/** An intron's row for the sQTL tab, from its search-index row, its block and the variants range.
 *  Intron coordinates come from the details' `extra`, else from the leafcutter phenotype id. */
function spliceRow(hit: SearchHit, p: Placed, block: ResultBlock, range: VariantRange | null): SplicePhenotype {
  const x = block.details.extra as { intron_start?: number; intron_end?: number; cluster_id?: string; strand?: string }
  const m = /^[^:]+:(\d+):(\d+):(clu_\d+_([+-?]))(?::|$)/.exec(p.phenotype_id)
  const g = block.details.group, lead = g?.lead ?? null, row = leadRow(block, range)
  return {
    phenotype_id: p.phenotype_id, gene_id: hit.gene_id, symbol: hit.symbol, chr: hit.chr,
    intron_start: x.intron_start ?? Number(m?.[1]), intron_end: x.intron_end ?? Number(m?.[2]),
    cluster_id: x.cluster_id ?? m?.[3] ?? '', strand: x.strand ?? m?.[4] ?? '', tss: hit.tss,
    num_var: g?.n_variants ?? block.nRows,
    lead_position: lead?.pos ?? null, lead_A1: lead?.alt ?? null, lead_A2: lead?.ref ?? null, lead_rsid: row?.rsid ?? null,
    lead_af: row?.af ?? null, lead_tss_distance: lead ? lead.pos - hit.tss : null,
    slope: row?.slope ?? null, slope_se: row?.se ?? null, pval_nominal: row?.pval ?? null,
    pval_perm: g?.p_perm ?? p.p_perm, pval_beta: g?.p_beta ?? null, qval: null, is_sqtl: p.significant,
    n_credible_sets: block.details.n_credible_sets, blk_off: p.blk_off, blk_len: p.blk_len,
  }
}

// ---- opening a gene ---------------------------------------------------------------------------

function openGene(hit: SearchHit): GenePack {
  const s = getStore()
  loadGwasIndex().catch(() => {})       // its fetch overlaps the phenotype lookup the GWAS window waits for
  const detail0 = { gene_id: hit.gene_id }
  const phen = phenotypesOfGene(hit.gene_id, hit.chr).then(ps => ps.filter((p): p is Placed => p.blk_off != null && p.blk_len != null))
  const block = s.then(st => {
    if (hit.blk_off == null || hit.blk_len == null) return null
    const f = resultsFile(st, EQTL_TYPE, hit.chr)
    return fetchDecoded('store:block', detail0, f.name, hit.blk_off, hit.blk_len,
      (b, what) => decodeResultBlock(b, f.dof, { blk_len: hit.blk_len, n_var: hit.n_var, var_start: hit.var_start }, what))
  })
  const variants = Promise.all([s, phen]).then(([st, ps]) => {
    const withRuns = ps.filter(p => p.var_off != null && p.var_len != null)
    if (!withRuns.length) return null
    const off = Math.min(...withRuns.map(p => p.var_off!))
    const len = Math.max(...withRuns.map(p => p.var_off! + p.var_len!)) - off
    const f = variantsFile(st, hit.chr)
    return fetchDecoded('store:variants', detail0, f.name, off, len, (b, what) => decodeVariantRange(b, len, what))
  })
  // the GWAS window spans every run of the gene's phenotypes, as the variants range does
  const gwas = phen.then(ps => {
    const runs = ps.filter(p => p.w_lo != null && p.w_hi != null)
    return gwasTable(hit, runs.length ? Math.min(...runs.map(p => p.w_lo!)) : null, runs.length ? Math.max(...runs.map(p => p.w_hi!)) : null)
  })
  const detail = Promise.all([block, variants, exonModel(hit.gene_id, hit.chr)])
    .then(([b, v, ex]) => ({ gene: geneRow(hit, b, v), exons: ex }))

  let introns: Promise<{ list: SplicePhenotype[]; blocks: Map<string, ResultBlock> }> | null = null
  const loadIntrons = () => {
    if (!introns) {
      const p = Promise.all([s, phen, variants]).then(async ([st, ps, range]) => {
        const rowsS = ps.filter(p => p.phenotype_type === SQTL_TYPE)
        const blocks = new Map<string, ResultBlock>()
        if (!rowsS.length) return { list: [], blocks }
        const f = resultsFile(st, SQTL_TYPE, hit.chr)
        const off = Math.min(...rowsS.map(p => p.blk_off))
        const len = Math.max(...rowsS.map(p => p.blk_off + p.blk_len)) - off
        const decoded = await fetchDecoded('store:introns', detail0, f.name, off, len, (span, what) =>
          rowsS.map(p => decodeResultBlock(span.subarray(p.blk_off - off, p.blk_off - off + p.blk_len), f.dof,
            { blk_len: p.blk_len, n_var: p.n_var, var_start: p.var_start }, `${what} ${p.phenotype_id}`)))
        rowsS.forEach((p, i) => blocks.set(p.phenotype_id, decoded[i]))
        const list = rowsS.map((p, i) => spliceRow(hit, p, decoded[i], range))
          .sort((a, b) => a.intron_start - b.intron_start || a.intron_end - b.intron_end || a.phenotype_id.localeCompare(b.phenotype_id))
        return { list, blocks }
      })
      introns = p
      p.catch(() => { if (introns === p) introns = null })
    }
    return introns
  }
  const intron = (id: string) => loadIntrons().then(x => {
    const b = x.blocks.get(id)
    if (!b) throw new Error(`${hit.gene_id}: ${id} is not one of the gene's tested introns`)
    return b
  })
  // consumers see each rejection when they await; these handlers only keep an unused one quiet
  for (const p of [block, variants, gwas, detail]) p.catch(() => {})
  return { hit, block, variants, gwas, detail, splice: () => loadIntrons().then(x => x.list), intron }
}

// the last few genes stay open, so tab switches and back navigation send no request
const CACHE_SIZE = 4
const genes = new Map<string, GenePack>()
const release = (gp: GenePack) => { gp.gwas.then(dropTable, () => {}) }

/** Opens a gene (starting its requests) or returns it from the cache. */
export function loadGene(hit: SearchHit): GenePack {
  const key = hit.gene_id
  const have = genes.get(key)
  if (have) { genes.delete(key); genes.set(key, have); return have }
  const gp = openGene(hit)
  genes.set(key, gp)
  const forget = () => { if (genes.get(key) === gp) { genes.delete(key); release(gp) } }
  for (const p of [gp.block, gp.variants, gp.gwas]) p.catch(forget)
  while (genes.size > CACHE_SIZE) {
    const [k, old] = genes.entries().next().value as [string, GenePack]
    genes.delete(k)
    release(old)
  }
  return gp
}

const blockFor = async (gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<ResultBlock> => {
  if (qtlType === 's') {
    if (!phenotypeId) throw new Error(`${gp.hit.gene_id}: an sQTL locus needs an intron`)
    return gp.intron(phenotypeId)
  }
  const b = await gp.block
  if (!b) throw new Error(`${gp.hit.gene_id}: the gene has no eQTL block`)
  return b
}

// ---- the locus table --------------------------------------------------------------------------

/** The in-band null codes as SQL NULLs. */
const rawSource = (raw: string) => `(SELECT position, A1, A2,
  CASE WHEN rs_number = 0 THEN NULL ELSE rs_number::BIGINT END AS rs_number,
  tss_distance,
  CASE WHEN isnan(af) THEN NULL ELSE af END AS af,
  CASE WHEN ma_samples < 0 THEN NULL ELSE ma_samples END AS ma_samples,
  CASE WHEN ma_count < 0 THEN NULL ELSE ma_count END AS ma_count,
  CASE WHEN isnan(pval_nominal) THEN NULL ELSE pval_nominal END AS pval_nominal,
  CASE WHEN isnan(slope) THEN NULL ELSE slope END AS slope,
  CASE WHEN isnan(slope_se) THEN NULL ELSE slope_se END AS slope_se,
  CASE WHEN cs_id < 0 THEN NULL ELSE pip END AS pip,
  CASE WHEN cs_id < 0 THEN NULL ELSE cs_id END AS cs_id
  FROM ${raw})`

/** One cis window as a table for the plots: -log10 p, credible-set class, a tooltip label, and
 *  the GWAS statistics for variants present there (matched on position and alleles in either
 *  orientation, GWAS beta re-signed to the QTL effect allele A1, which is ALT). Ordered so
 *  credible-set variants are drawn last (on top). A row with no p (a site the phenotype did not
 *  test) has no -log10 p; a results set without dof has no slopes, and the label says so. */
const locusSQL = (qtl: string, gwas: string) => `
    SELECT q.position,
           -- p underflows to 0 for a few extreme variants: place them just above the largest finite value
           CASE WHEN q.pval_nominal IS NULL THEN NULL
                ELSE coalesce(-log10(nullif(q.pval_nominal, 0)), max(-log10(nullif(q.pval_nominal, 0))) OVER () * 1.05) END AS nlp,
           q.pval_nominal = 0 AS clipped,
           q.pval_nominal, q.slope, q.slope_se, q.af, q.pip, q.cs_id, q.rs_number, q.A1, q.A2,
           q.tss_distance, q.ma_samples, q.ma_count,
           coalesce(q.cs_id::VARCHAR, 'none') AS cs,
           g.p AS gwas_p, -log10(g.p) AS gwas_nlp,
           CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END AS gwas_beta,
           coalesce('rs' || q.rs_number, q.position::VARCHAR) || '  ' || q.A1 || '/' || q.A2
             || chr(10) || CASE WHEN q.pval_nominal IS NULL THEN 'not tested' WHEN q.pval_nominal = 0 THEN 'p = 0 (underflow; drawn above the maximum)' ELSE 'p = ' || format('{:.2e}', q.pval_nominal) END
             || CASE WHEN q.slope IS NULL THEN coalesce(chr(10) || 'SE ' || format('{:.3f}', q.slope_se), '')
                     ELSE chr(10) || 'slope ' || format('{:.3f}', q.slope) || ' ± ' || format('{:.3f}', q.slope_se) END
             || coalesce(chr(10) || 'AF ' || format('{:.3f}', q.af), '')
             || CASE WHEN q.pip IS NULL THEN '' ELSE chr(10) || 'PIP ' || format('{:.3f}', q.pip) || ' (set ' || q.cs_id || ')' END
             || CASE WHEN g.p IS NULL THEN '' ELSE chr(10) || 'DCM GWAS p = ' || format('{:.2e}', g.p) || ', beta ' || format('{:+.3f}', CASE WHEN g.ea = q.A1 THEN g.beta ELSE -g.beta END) || ' (A1 as effect allele)' END AS label
    FROM ${qtl} q
    LEFT JOIN ${gwas} g
      ON g.position = q.position AND ((g.ea = q.A1 AND g.nea = q.A2) OR (g.ea = q.A2 AND g.nea = q.A1))
    -- a GWAS may list some indels in both allele orientations as separate records: keep one per
    -- QTL variant, preferring the orientation that matches the QTL alleles as written
    QUALIFY row_number() OVER (PARTITION BY q.position, q.A1, q.A2 ORDER BY (g.ea = q.A1) DESC NULLS LAST, g.p) = 1
    ORDER BY q.cs_id IS NOT NULL, q.position`

/** The locus table of the gene's eQTL rows or one intron's sQTL rows, as an in-memory DuckDB
 *  table; the caller drops it. */
export async function locusTable(gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<string> {
  const [block, range, gwas] = await Promise.all([blockFor(gp, qtlType, phenotypeId), gp.variants, gp.gwas])
  const detail = { gene_id: gp.hit.gene_id, phenotype_id: phenotypeId }
  const t0 = performance.now()
  if (block.varStart == null || !range) throw new Error(`${gp.hit.gene_id}: the ${qtlType === 'e' ? 'eQTL' : 'intron'} block has no rows`)
  const c = readerColumns(block, sliceRun(range, block.varStart, block.nRows, block, `${gp.hit.gene_id} variants range`), gp.hit.tss)
  const ipc = tableToIPC(new Table({
    position: makeVector(c.position), A1: vectorFromArray(c.a1, new Utf8()), A2: vectorFromArray(c.a2, new Utf8()),
    rs_number: makeVector(c.rsNumber), tss_distance: makeVector(c.tssDistance), af: makeVector(c.af),
    ma_samples: makeVector(c.maSamples), ma_count: makeVector(c.maCount), pval_nominal: makeVector(c.pval),
    slope: makeVector(c.slope), slope_se: makeVector(c.se), pip: makeVector(c.pip), cs_id: makeVector(c.csId),
  }), 'stream')
  performance.measure('store:decode', { start: t0, end: performance.now(), detail: { ...detail, part: 'columns' } })
  const { con } = await getDB()
  const q = tableName('q')
  try {
    await insertArrow(q, ipc, { ...detail, part: 'locus' })
    const name = tableName('locus')
    await con.query(`CREATE TABLE ${name} AS ${locusSQL(rawSource(q), gwas)}`)
    return name
  } finally {
    await con.query(`DROP TABLE IF EXISTS ${q}`).catch(() => {})
  }
}

/** Every credible-set membership of the gene's eQTL rows or one intron, in the credible-set table's
 *  row shape; a variant in two sets is listed under both. */
export async function credibleSets(gp: GenePack, qtlType: 'e' | 's', phenotypeId?: string): Promise<CredibleSetRow[]> {
  const [block, range] = await Promise.all([blockFor(gp, qtlType, phenotypeId), gp.variants])
  if (block.varStart == null || !range) return []
  const run = sliceRun(range, block.varStart, block.nRows, block, `${gp.hit.gene_id} variants range`)
  return csMembers(block)
    .map(({ row, pip, csId }) => ({
      qtl_type: qtlType, phenotype_id: phenotypeId ?? gp.hit.gene_id, chr: gp.hit.chr,
      position: run.position[row], A1: run.alt[row], A2: run.ref[row], rsid: run.rsNumber[row] ? `rs${run.rsNumber[row]}` : null,
      af: run.af[row], cs_id: csId, pip,
    }))
    .sort((a, b) => a.cs_id - b.cs_id || b.pip - a.pip)
}
