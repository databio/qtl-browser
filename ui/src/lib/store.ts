/**
 * The qtlstore (SPEC.md sections 2 and 3): `store.json` -> `experiments/<id>.json` -> its variant
 * catalog and annotation pointers -> immutable objects by name. Every data URL in the app is built here,
 * so the only objects ever fetched are ones a pointer names.
 *
 * Pointers are fetched once per session with `no-cache` (they are mutable); objects are
 * content-addressed and read whole or by HTTP range. A ranged object's 64-byte header is fetched
 * once, alongside its first range, and checked for the kind, chromosome and `seq_digest` its
 * pointer promised before any of its bytes are used.
 *
 * The reader opens one experiment (`VITE_EXPERIMENT`, default `topchef`). Trans results
 * (`results[].trans`) and the GWAS (`gwas`) are optional per experiment; `hasTrans` / `hasGwas`
 * say whether a page can show them.
 */
import { tableFromIPC } from 'apache-arrow'
import { checkFileHeader, decodeArrowObject, decodeGwasBins, DIGEST, FORMAT_VERSION, HEADER_LEN, KIND, OBJECT_NAME,
  parseFileHeader, type FileHeader, type GwasBin } from './store-decode'

/** The one data host. An empty `VITE_DATA_BASE` means "same origin" (`/data`), which is what the
 *  local dev server and preview serve. */
export const DATA_BASE: string =
  (import.meta.env.VITE_DATA_BASE as string | undefined)?.replace(/\/$/, '') || `${window.location.origin}/data`

/** The experiment the pages show. */
export const EXPERIMENT_ID: string = (import.meta.env.VITE_EXPERIMENT as string | undefined) || 'topchef'

/** The phenotype types the gene page's two tabs read (TOPCHeF's `ge` and `leafcutter`). */
export const EQTL_TYPE = 'ge'
export const SQTL_TYPE = 'leafcutter'

// ---- pointer documents ------------------------------------------------------------------------

/** The variant catalogs' pointer directory, and the `store.json` key listing their ids. */
export const CATALOG_DIR = 'variant_catalogs'

export interface StoreDoc {
  name: string; format_version: number; refget: string[]
  [CATALOG_DIR]: string[]; annotations: string[]; experiments: string[]
}

export interface CatalogChrom { name: string; seq_digest: string; length: number; count: number; n_cis: number; file: string }
export interface CatalogDoc {
  id: string; identity_digest: string; collection_digest: string; orientation: string
  attributes: string[]; page_size: number; n_sites: number
  chromosomes: CatalogChrom[]; vidx: string; rsid: string
}

export interface AnnotationDoc {
  id: string; identity_digest: string; genes: string; exons: string
  n_genes: number; n_transcripts: number; n_exons: number
  source: { file?: string; name?: string; version?: string; url?: string }
}

/** Worst-case rounding of one results set (SPEC section 8, `precision`). */
export interface ResultsPrecision { neglog10p_max_error: number; slope_se_max_rel_error: number; af_max_error: number }
/** A phenotype type's trans object (SPEC section 9). */
export interface TransSet { file: string; n_rows: number; n_phenotypes: number; precision: Record<string, number> }
export interface ResultsSet {
  phenotype_type: string
  /** the type's trans object, null when it has no trans rows */
  trans?: TransSet | null
  /** degrees of freedom for rebuilding slopes; null means no slope may be shown */
  dof: number | null
  n_phenotypes: number
  files: Record<string, string>
  precision: ResultsPrecision
}
export interface ExperimentDoc {
  id: string; catalog: string; catalog_identity: string; annotation: string
  allele_orientation_source: string | null
  significance: { column: string; op: string; threshold: number }
  search_index: string; n_phenotypes: number
  hits: Record<string, string>
  results: ResultsSet[]
  source?: { experiment_id?: string; source?: Record<string, unknown> }
  /** trans row counts, null when the experiment has no trans table */
  trans?: { rows: number } | null
  /** the experiment's GWAS (SPEC section 10), null when it has none */
  gwas?: GwasDoc | null
}

export interface GwasDoc {
  id: string; title: string; files: Record<string, string>; index: string
  n_rows: number; rows_by_chrom: Record<string, number>; block_rows: number; n_values: number[]
  /** the bin summary (SPEC section 10, "Bin summary"), null or absent when the store has none */
  bins?: { file: string; bin_bp: number; n_bins: number } | null
  source?: { file?: string; n_cases?: number; n_controls?: number }
}

/** One experiment with everything its pointers name, resolved and cross-checked. */
export interface Store {
  store: StoreDoc
  experiment: ExperimentDoc
  catalog: CatalogDoc
  annotation: AnnotationDoc
  /** catalog chromosome table by name, with its 1-based ordinal */
  chroms: Map<string, CatalogChrom & { ordinal: number }>
  /** results sets by phenotype type */
  results: Map<string, ResultsSet>
  hasTrans: boolean
  hasGwas: boolean
}

// ---- fetching ---------------------------------------------------------------------------------

async function pointer<T>(path: string): Promise<T> {
  const r = await fetch(`${DATA_BASE}/${path}`, { cache: 'no-cache' })
  if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`)
  return r.json() as Promise<T>
}

/** A full URL for an object a pointer names, `<digest>.<ext>` under `immutable/`. */
export function objectUrl(name: string): string {
  if (!OBJECT_NAME.test(name)) throw new Error(`"${name}" is not an object name (<digest>.<ext>)`)
  return `${DATA_BASE}/immutable/${name}`
}

/** A whole object. */
export async function fetchObject(name: string): Promise<Uint8Array> {
  const r = await fetch(objectUrl(name))
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`)
  return new Uint8Array(await r.arrayBuffer())
}

async function rangeFetch(name: string, offset: number, length: number): Promise<Uint8Array> {
  const range = `bytes=${offset}-${offset + length - 1}`
  const r = await fetch(objectUrl(name), { headers: { Range: range } })
  // a 200 means the server ignored Range and is sending the whole object: fail instead of downloading it
  if (r.status !== 206) throw new Error(`${name} ${range}: expected 206, got ${r.status}`)
  const buf = new Uint8Array(await r.arrayBuffer())
  if (buf.length !== length) throw new Error(`${name} ${range}: got ${buf.length} bytes, expected ${length}`)
  return buf
}

/** A byte range of an object with no header check (the caller checks what it reads). */
export const rangeObject = (name: string, offset: number, length: number) => rangeFetch(name, offset, length)

/** What a ranged object's header must say. */
export interface HeaderExpect { kind: number; chrom: string; seqDigest: string; count?: number; nCis?: number }

const headers = new Map<string, Promise<FileHeader>>()
/** The object's 64-byte header, fetched once per session and checked against `expect`. */
function checkedHeader(name: string, expect: HeaderExpect): Promise<FileHeader> {
  let p = headers.get(name)
  if (!p) {
    p = rangeFetch(name, 0, HEADER_LEN).then(b => parseFileHeader(b, name))
    headers.set(name, p)
    const q = p
    q.catch(() => { if (headers.get(name) === q) headers.delete(name) })
  }
  return p.then(h => {
    checkFileHeader(h, expect, name)
    if (expect.count !== undefined && h.count !== expect.count) throw new Error(`${name}: header count ${h.count}, pointer says ${expect.count}`)
    if (expect.nCis !== undefined && h.nCis !== expect.nCis) throw new Error(`${name}: header n_cis ${h.nCis}, pointer says ${expect.nCis}`)
    return h
  })
}

/** One range request timed as `mark` (store:block, store:variants, ...), then its decode as
 *  `decodeMark`. The object's header check runs alongside the range, and the decode waits for it. */
export async function fetchDecoded<T>(mark: string, detail: Record<string, unknown>, name: string, expect: HeaderExpect,
  offset: number, length: number, decode: (bytes: Uint8Array, what: string) => T, decodeMark = 'store:decode'): Promise<T> {
  const t0 = performance.now()
  const [, bytes] = await Promise.all([checkedHeader(name, expect), rangeFetch(name, offset, length)])
  const t1 = performance.now()
  performance.measure(mark, { start: t0, end: t1, detail })
  const out = decode(bytes, `${name} bytes ${offset}+${length}`)
  performance.measure(decodeMark, { start: t1, end: performance.now(), detail: { ...detail, part: mark.replace(/^store:/, '') } })
  return out
}

// ---- the store --------------------------------------------------------------------------------

function checkNames(doc: unknown, where: string): void {
  // every string that looks like an object name must be one (qtlstore.object_names)
  const walk = (x: unknown) => {
    if (typeof x === 'string') { if (/\.(qb[vehrx]|arrow\.zst)$/.test(x) && !OBJECT_NAME.test(x)) throw new Error(`${where}: bad object name ${x}`) }
    else if (Array.isArray(x)) x.forEach(walk)
    else if (x && typeof x === 'object') Object.values(x).forEach(walk)
  }
  walk(doc)
}

async function openStore(): Promise<Store> {
  const store = await pointer<StoreDoc>('store.json')
  if (store.format_version !== FORMAT_VERSION)
    throw new Error(`store.json: format version ${store.format_version}, and this build reads version ${FORMAT_VERSION}. The data was updated; reload the page.`)
  if (!store.experiments?.includes(EXPERIMENT_ID)) throw new Error(`store.json lists no experiment "${EXPERIMENT_ID}"`)
  const experiment = await pointer<ExperimentDoc>(`experiments/${EXPERIMENT_ID}.json`)
  if (experiment.id !== EXPERIMENT_ID) throw new Error(`experiments/${EXPERIMENT_ID}.json: id is "${experiment.id}"`)
  if (!store[CATALOG_DIR]?.includes(experiment.catalog) || !store.annotations?.includes(experiment.annotation))
    throw new Error(`experiment ${experiment.id}: variant catalog ${experiment.catalog} or annotation ${experiment.annotation} is not in store.json`)
  const [catalog, annotation] = await Promise.all([
    pointer<CatalogDoc>(`${CATALOG_DIR}/${experiment.catalog}.json`),
    pointer<AnnotationDoc>(`annotations/${experiment.annotation}.json`),
  ])
  if (experiment.catalog_identity !== catalog.identity_digest)
    throw new Error(`experiment ${experiment.id}: catalog_identity ${experiment.catalog_identity} is not catalog ${catalog.id}'s ${catalog.identity_digest}`)
  if (catalog.orientation !== 'ref_alt') throw new Error(`catalog ${catalog.id}: orientation "${catalog.orientation}", this reader reads ref_alt`)
  if (!DIGEST.test(catalog.collection_digest)) throw new Error(`catalog ${catalog.id}: collection digest ${catalog.collection_digest}`)
  checkNames(experiment, `experiments/${experiment.id}.json`)
  checkNames(catalog, `${CATALOG_DIR}/${catalog.id}.json`)
  checkNames(annotation, `annotations/${annotation.id}.json`)
  const chroms = new Map(catalog.chromosomes.map((c, i) => [c.name, { ...c, ordinal: i + 1 }]))
  const results = new Map(experiment.results.map(r => [r.phenotype_type, r]))
  for (const [c] of Object.entries(experiment.hits)) if (!chroms.has(c)) throw new Error(`experiment ${experiment.id}: hits for ${c}, not in the catalog`)
  for (const r of experiment.results) for (const c of Object.keys(r.files))
    if (!chroms.has(c)) throw new Error(`experiment ${experiment.id}: ${r.phenotype_type} results for ${c}, not in the catalog`)
  return { store, experiment, catalog, annotation, chroms, results,
    hasTrans: experiment.results.some(r => r.trans != null), hasGwas: experiment.gwas != null }
}

let opened: Promise<Store> | null = null
/** Memoized; a failed open is forgotten so the next call retries. */
export function getStore(): Promise<Store> {
  if (!opened) {
    const p = openStore()
    opened = p
    p.catch(() => { if (opened === p) opened = null })
  }
  return opened
}

function memo<T>(make: () => Promise<T>): () => Promise<T> {
  let p: Promise<T> | null = null
  return () => {
    if (!p) {
      const q = make()
      p = q
      q.catch(() => { if (p === q) p = null })
    }
    return p
  }
}

/** The experiment's search index as an Arrow IPC stream (one row per phenotype, SPEC section 8). */
export const searchIndexIPC = memo(async () => {
  const s = await getStore()
  return decodeArrowObject(await fetchObject(s.experiment.search_index), s.experiment.search_index)
})

/** The annotation's gene table as an Arrow IPC stream (SPEC section 6). */
export const genesIPC = memo(async () => {
  const s = await getStore()
  return decodeArrowObject(await fetchObject(s.annotation.genes), s.annotation.genes)
})

/** The annotation's exon table as an Arrow IPC stream: one whole-object read, the first time a
 *  gene page needs its exon model. */
export const exonsIPC = memo(async () => {
  const s = await getStore()
  return decodeArrowObject(await fetchObject(s.annotation.exons), s.annotation.exons)
})

/** The GWAS bin summary as rows (one whole-object read, for the landing track), or null when the
 *  experiment's GWAS has none. */
export const gwasBins = memo(async (): Promise<GwasBin[] | null> => {
  const s = await getStore()
  const b = s.experiment.gwas?.bins
  if (!b) return null
  return decodeGwasBins(await fetchObject(b.file), b, b.file)
})

/** The results file of one phenotype type on one chromosome, with the header it must carry. */
export function resultsFile(s: Store, phenotypeType: string, chr: string): { name: string; expect: HeaderExpect; dof: number | null } {
  const r = s.results.get(phenotypeType)
  const name = r?.files[chr]
  const c = s.chroms.get(chr)
  if (!r || !name || !c) throw new Error(`experiment ${s.experiment.id}: no ${phenotypeType} results for ${chr}`)
  return { name, expect: { kind: KIND.results, chrom: chr, seqDigest: c.seq_digest }, dof: r.dof }
}

/** A phenotype type's trans object, with the header it must carry, or null when it has none. */
export function transFile(s: Store, phenotypeType: string): { name: string; expect: HeaderExpect; dof: number | null } | null {
  const r = s.results.get(phenotypeType)
  if (!r?.trans) return null
  return { name: r.trans.file, expect: { kind: KIND.trans, chrom: 'all', seqDigest: s.catalog.collection_digest }, dof: r.dof }
}

/** A chromosome's GWAS file, with the header it must carry, or null when the GWAS has no rows there. */
export function gwasFile(s: Store, chr: string): { name: string; expect: HeaderExpect } | null {
  const g = s.experiment.gwas, c = s.chroms.get(chr)
  const name = g?.files[chr]
  if (!g || !name || !c) return null
  return { name, expect: { kind: KIND.gwas, chrom: chr, seqDigest: c.seq_digest, count: g.rows_by_chrom[chr] } }
}

/** A catalog chromosome's variants file, with the header it must carry. */
export function variantsFile(s: Store, chr: string): { name: string; expect: HeaderExpect } {
  const c = s.chroms.get(chr)
  if (!c) throw new Error(`catalog ${s.catalog.id}: no chromosome ${chr}`)
  return { name: c.file, expect: { kind: KIND.variants, chrom: chr, seqDigest: c.seq_digest, count: c.count, nCis: c.n_cis } }
}

// ---- what the About, Home and gene pages print --------------------------------------------------

/** Whole-experiment counts, from the search index and the catalog. */
export interface Counts {
  genes_tested?: number; egenes?: number
  splice_phenotypes_tested?: number; sqtl_sig_phenotypes?: number; sqtl_sig_genes?: number
  variants_cis?: number; variants_trans_only?: number
  gwas_variants?: number; trans_pairs?: number
}
export interface PrecisionKind { neglog10p_max_error: number; slope_se_max_rel_error: number }
export interface Precision { af_max_error: number; eqtl: PrecisionKind; sqtl: PrecisionKind }

/** What the pages print about the release: the resolved store plus counts, sources and precision. */
export interface StoreInfo extends Store {
  counts: Counts
  sources: Record<string, { version: string; description: string }>
  precision: Precision | null
  /** annotation release label, e.g. "v34 (GRCh38.p13)" */
  annotationVersion: string | null
}

export const getStoreInfo = memo(async (): Promise<StoreInfo> => {
  const s = await getStore()
  const t = tableFromIPC(await searchIndexIPC())
  const type = t.getChild('phenotype_type')!, sig = t.getChild('significant')!, nVar = t.getChild('n_var')!, gene = t.getChild('gene_id')!
  const c: Counts = { genes_tested: 0, egenes: 0, splice_phenotypes_tested: 0, sqtl_sig_phenotypes: 0 }
  const sigGenes = new Set<string>()
  const chr = t.getChild('chr')!
  for (let i = 0; i < t.numRows; i++) {
    const pt = type.get(i), significant = sig.get(i) === true
    if (chr.get(i) == null) continue            // trans-only phenotypes: no cis test
    if (pt === EQTL_TYPE) {
      if (nVar.get(i) != null) c.genes_tested!++
      if (significant) c.egenes!++
    } else if (pt === SQTL_TYPE) {
      c.splice_phenotypes_tested!++
      if (significant) { c.sqtl_sig_phenotypes!++; const g = gene.get(i); if (g) sigGenes.add(g) }
    }
  }
  c.sqtl_sig_genes = sigGenes.size
  c.variants_cis = s.catalog.chromosomes.reduce((a, x) => a + x.n_cis, 0)
  c.variants_trans_only = s.catalog.n_sites - c.variants_cis
  if (s.experiment.gwas) c.gwas_variants = s.experiment.gwas.n_rows
  if (s.experiment.trans) c.trans_pairs = s.experiment.trans.rows
  const e = s.results.get(EQTL_TYPE)?.precision, q = s.results.get(SQTL_TYPE)?.precision
  const precision = e && q ? {
    af_max_error: Math.max(e.af_max_error, q.af_max_error),
    eqtl: { neglog10p_max_error: e.neglog10p_max_error, slope_se_max_rel_error: e.slope_se_max_rel_error },
    sqtl: { neglog10p_max_error: q.neglog10p_max_error, slope_se_max_rel_error: q.slope_se_max_rel_error },
  } : null
  const version = s.annotation.source?.version ?? null
  const src = s.experiment.source?.source as { zenodo?: string } | undefined
  const sources: StoreInfo['sources'] = {
    gencode: { version: version ?? s.annotation.id, description: `gene annotation (${s.annotation.id})` },
    catalog: { version: s.catalog.identity_digest, description: `variant catalog ${s.catalog.id}, ${s.catalog.n_sites.toLocaleString('en-US')} sites` },
    [s.experiment.id]: { version: src?.zenodo ?? s.experiment.id, description: `experiment ${s.experiment.id} (store ${s.store.name})` },
  }
  if (s.experiment.gwas) sources.gwas = { version: s.experiment.gwas.id, description: s.experiment.gwas.title }
  return { ...s, counts: c, sources, precision, annotationVersion: version?.split(' ')[0] ?? null }
})
