/**
 * The qtlstore (SPEC.md sections 2 and 3): `store.json` -> `experiments/<id>.json` -> its variant
 * catalog and annotation pointers -> immutable objects by name. Every data URL in the app is built here,
 * so the only objects ever fetched are ones a pointer names.
 *
 * Pointers are fetched once per session with `no-cache` (they are mutable); objects are
 * content-addressed and read whole or by HTTP range. A range read sends no separate header request:
 * the pointer says what an object is, Store.validate checked every header against its pointer at
 * build time, and a header is checked wherever it arrives inside bytes read anyway (whole objects,
 * the gene lookup's directory, a hits file's frame table).
 *
 * A gene or variant page reads per-chromosome objects only (SPEC sections 6 and 8): its genes, exon
 * models and search index part, and one gene lookup bucket to find a gene's chromosome. Nothing on
 * a first gene page reads a whole-genome table; the gene list and the search box read the
 * whole-genome genes table and search index the first time they are used (two requests), and the
 * browser cache keeps every object after.
 *
 * The reader opens one experiment (`VITE_EXPERIMENT`, default `topchef`). Trans results
 * (`results[].trans`) and the GWAS (`gwas`) are optional per experiment; `hasTrans` / `hasGwas`
 * say whether a page can show them.
 */
import { tableFromIPC, type Table } from 'apache-arrow'
import { decodeArrowObject, decodeGwasBins, decodeLookupDir, DIGEST, FORMAT_VERSION, lookupBucket,
  lookupDirLen, lookupKey, OBJECT_NAME, parseFileHeader, type GwasBin } from './store-decode'

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
  id: string; identity_digest: string
  /** the whole-genome tables (SPEC section 6): read only for the gene list and search (`genes`) */
  genes: string; exons: string
  /** per chromosome: its genes rows and their collapsed exon models */
  chroms: Record<string, { genes: string; exon_models: string }>
  /** the gene lookup (kind 10): gene_id and symbol -> gene_id, name, chr, tss */
  lookup: string
  n_genes: number; n_transcripts: number; n_exons: number
  source: { file?: string; name?: string; version?: string; url?: string }
}

/** One part of the search index (SPEC section 8): its object, row count and inclusive ord runs. */
export interface IndexPart { file: string; rows: number; ords: [number, number][] }
/** What the Home and About pages print for one phenotype type (phenotypes with a cis result). */
export interface TypeCounts { phenotypes: number; with_rows: number; significant: number; significant_genes: number }

/** How much one results set was rounded (SPEC section 8, `precision`). The first three are
 *  worst-case bounds; `slope_max_error_over_se` is measured at build time over `slope_rows_compared`
 *  rows, and is null when the set has no dof (no slope is shown then). */
export interface ResultsPrecision {
  neglog10p_max_error: number; slope_se_max_rel_error: number; af_max_error: number
  slope_max_error_over_se: number | null; slope_rows_compared?: number
}
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
  /** the search index by chromosome, which the browser reads instead of `search_index` */
  search_index_parts: Record<string, IndexPart>
  /** the trans-only phenotypes' rows (`chr` null), null when there are none */
  search_index_trans_only: IndexPart | null
  counts: Record<string, TypeCounts>
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

/** One range request timed as `mark` (store:block, store:variants, ...), then its decode as
 *  `decodeMark`. No header request: the pointer already says what the object is (its name is the
 *  digest of its bytes, and Store.validate checked its header at build time), and the decoders
 *  fail on bytes that are not what the offsets promise. */
export async function fetchDecoded<T>(mark: string, detail: Record<string, unknown>, name: string,
  offset: number, length: number, decode: (bytes: Uint8Array, what: string) => T, decodeMark = 'store:decode'): Promise<T> {
  const t0 = performance.now()
  const bytes = await rangeFetch(name, offset, length)
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
    // every extension SPEC section 3 lists: qbv qbe qbg qbt qbh qbr qbx, qgi qgl, arrow.zst
    if (typeof x === 'string') { if (/\.(qb[vegthrx]|qg[il]|arrow\.zst)$/.test(x) && !OBJECT_NAME.test(x)) throw new Error(`${where}: bad object name ${x}`) }
    else if (Array.isArray(x)) x.forEach(walk)
    else if (x && typeof x === 'object') Object.values(x).forEach(walk)
  }
  walk(doc)
}

/** The catalog and annotation a previous open found this experiment naming (lib/store.ts is the only
 *  writer). Kept per data host and experiment so a dev build pointed at another store cannot inherit
 *  ids from this one. */
interface Remembered { catalog: string; annotation: string }
const REMEMBER_KEY = `qtlstore:${DATA_BASE}:${EXPERIMENT_ID}`

function remembered(): Remembered | null {
  try {
    const r = JSON.parse(localStorage.getItem(REMEMBER_KEY) ?? 'null') as Remembered | null
    return typeof r?.catalog === 'string' && typeof r?.annotation === 'string' ? r : null
  } catch { return null }   // no localStorage, or a value this build did not write: guess nothing
}

function remember(r: Remembered): void {
  try { localStorage.setItem(REMEMBER_KEY, JSON.stringify(r)) } catch { /* not worth failing an open over */ }
}

/** A speculative pointer fetch when it was for the id the experiment turned out to name, else a
 *  fresh one. A guess that 404s (the store was rebuilt) falls through to the real id. */
async function speculated<T>(spec: Promise<T> | null, guessed: string | undefined, want: string,
  fetchIt: (id: string) => Promise<T>): Promise<T> {
  if (spec && guessed === want) {
    const got = await spec.catch(() => null)
    if (got) return got
  }
  return fetchIt(want)
}

async function openStore(): Promise<Store> {
  // One wave, not two. The experiment's pointer is fetched alongside store.json (its id is a
  // build-time constant) and used only once store.json lists it. The catalog and annotation names
  // live inside the experiment document, so reading them would cost a second round trip; instead the
  // ids a previous open recorded are fetched in this same wave. A right guess -- the normal case,
  // since these ids change only when the store is rebuilt -- removes the round trip entirely, and a
  // wrong one costs only the two discarded requests. The pointer documents are served with
  // `no-cache` and no validator, so each round trip is a full origin fetch; this is the only part of
  // that the reader can do anything about.
  const guess = remembered()
  const expP = pointer<ExperimentDoc>(`experiments/${EXPERIMENT_ID}.json`)
  const catP = guess ? pointer<CatalogDoc>(`${CATALOG_DIR}/${guess.catalog}.json`) : null
  const annP = guess ? pointer<AnnotationDoc>(`annotations/${guess.annotation}.json`) : null
  for (const p of [expP, catP, annP]) p?.catch(() => {})
  const store = await pointer<StoreDoc>('store.json')
  if (store.format_version !== FORMAT_VERSION)
    throw new Error(`store.json: format version ${store.format_version}, and this build reads version ${FORMAT_VERSION}. The data was updated; reload the page.`)
  if (!store.experiments?.includes(EXPERIMENT_ID)) throw new Error(`store.json lists no experiment "${EXPERIMENT_ID}"`)
  const experiment = await expP
  if (experiment.id !== EXPERIMENT_ID) throw new Error(`experiments/${EXPERIMENT_ID}.json: id is "${experiment.id}"`)
  if (!store[CATALOG_DIR]?.includes(experiment.catalog) || !store.annotations?.includes(experiment.annotation))
    throw new Error(`experiment ${experiment.id}: variant catalog ${experiment.catalog} or annotation ${experiment.annotation} is not in store.json`)
  const [catalog, annotation] = await Promise.all([
    speculated(catP, guess?.catalog, experiment.catalog, id => pointer<CatalogDoc>(`${CATALOG_DIR}/${id}.json`)),
    speculated(annP, guess?.annotation, experiment.annotation, id => pointer<AnnotationDoc>(`annotations/${id}.json`)),
  ])
  // the two URLs can come from a remembered id, so check the documents are the ones the experiment
  // names rather than trusting the path they were fetched from
  if (catalog.id !== experiment.catalog) throw new Error(`experiment ${experiment.id}: catalog ${experiment.catalog} answered with id "${catalog.id}"`)
  if (annotation.id !== experiment.annotation) throw new Error(`experiment ${experiment.id}: annotation ${experiment.annotation} answered with id "${annotation.id}"`)
  if (experiment.catalog_identity !== catalog.identity_digest)
    throw new Error(`experiment ${experiment.id}: catalog_identity ${experiment.catalog_identity} is not catalog ${catalog.id}'s ${catalog.identity_digest}`)
  if (catalog.orientation !== 'ref_alt') throw new Error(`catalog ${catalog.id}: orientation "${catalog.orientation}", this reader reads ref_alt`)
  if (!DIGEST.test(catalog.collection_digest)) throw new Error(`catalog ${catalog.id}: collection digest ${catalog.collection_digest}`)
  // the per-chromosome objects every page reads (SPEC sections 6 and 8); a store built before them has none
  if (!annotation.chroms || !annotation.lookup || !experiment.search_index_parts || !experiment.counts || experiment.search_index_trans_only === undefined)
    throw new Error(`the store predates this build (no per-chromosome annotation or search index objects); rebuild it with the current pipeline`)
  checkNames(experiment, `experiments/${experiment.id}.json`)
  checkNames(catalog, `${CATALOG_DIR}/${catalog.id}.json`)
  checkNames(annotation, `annotations/${annotation.id}.json`)
  const chroms = new Map(catalog.chromosomes.map((c, i) => [c.name, { ...c, ordinal: i + 1 }]))
  const results = new Map(experiment.results.map(r => [r.phenotype_type, r]))
  for (const [c] of Object.entries(experiment.hits)) if (!chroms.has(c)) throw new Error(`experiment ${experiment.id}: hits for ${c}, not in the catalog`)
  for (const r of experiment.results) for (const c of Object.keys(r.files))
    if (!chroms.has(c)) throw new Error(`experiment ${experiment.id}: ${r.phenotype_type} results for ${c}, not in the catalog`)
  // only after every check: the next open guesses these two ids and skips a round trip
  remember({ catalog: experiment.catalog, annotation: experiment.annotation })
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

/** One memoized promise per key; a rejection is forgotten so the next call retries. */
function memoBy<T>(make: (key: string) => Promise<T>): (key: string) => Promise<T> {
  const m = new Map<string, Promise<T>>()
  return key => {
    let p = m.get(key)
    if (!p) {
      const q = make(key)
      m.set(key, q)
      q.catch(() => { if (m.get(key) === q) m.delete(key) })
      p = q
    }
    return p
  }
}

const arrowTable = async (name: string): Promise<Table> => tableFromIPC(decodeArrowObject(await fetchObject(name), name))

/** The whole-genome annotation genes table and experiment search index, one request each: what the
 *  gene list and the search box read, instead of every chromosome's objects (gene-index.ts). */
export const wholeGenes = memo(async () => arrowTable((await getStore()).annotation.genes))
export const wholeIndex = memo(async () => arrowTable((await getStore()).experiment.search_index))

/** One chromosome's annotation genes rows (SPEC section 6), or null when the annotation has none there. */
export const chromGenes = memoBy(async (chr): Promise<Table | null> => {
  const n = (await getStore()).annotation.chroms[chr]?.genes
  return n ? arrowTable(n) : null
})

/** One chromosome's collapsed exon models (gene_id, exon_starts, exon_ends), or null. */
export const chromExonModels = memoBy(async (chr): Promise<Table | null> => {
  const n = (await getStore()).annotation.chroms[chr]?.exon_models
  return n ? arrowTable(n) : null
})

/** One chromosome's search index rows (SPEC section 8), or null when the experiment has no part there. */
export const chromIndex = memoBy(async (chr): Promise<Table | null> => {
  const n = (await getStore()).experiment.search_index_parts[chr]?.file
  return n ? arrowTable(n) : null
})

/** The trans-only phenotypes' search index rows, or null when there are none. */
export const transOnlyIndex = memo(async (): Promise<Table | null> => {
  const p = (await getStore()).experiment.search_index_trans_only
  return p ? arrowTable(p.file) : null
})

/** The search index part holding row `ord`: a chromosome name, `null` for the trans-only part, or
 *  undefined when no part lists it. */
export function partOfOrd(s: Store, ord: number): string | null | undefined {
  const has = (p: IndexPart) => p.ords.some(([a, b]) => ord >= a && ord <= b)
  for (const [c, p] of Object.entries(s.experiment.search_index_parts)) if (has(p)) return c
  const t = s.experiment.search_index_trans_only
  return t && has(t) ? null : undefined
}

/** A gene lookup row (SPEC section 6, kind 10). */
export interface LookupRow { key: string; gene_id: string; name: string; chr: string; tss: number }

/** The lookup's header and bucket offsets: one range read of the first `lookupDirLen` bytes (the
 *  builder writes 1024 buckets; a lookup with another count costs one more read). */
const lookupDir = memo(async () => {
  const s = await getStore()
  const name = s.annotation.lookup
  let b = await rangeFetch(name, 0, lookupDirLen(1024))
  const n = parseFileHeader(b, name).count
  if (n !== 1024) b = await rangeFetch(name, 0, lookupDirLen(n))
  return { name, ...decodeLookupDir(b, s.annotation.identity_digest, name) }
})

const lookupBucketRows = memoBy(async (b): Promise<LookupRow[]> => {
  const d = await lookupDir()
  const k = Number(b), off = d.off[k], len = d.off[k + 1] - off
  if (!len) return []
  const t = tableFromIPC(decodeArrowObject(await rangeFetch(d.name, off, len), `${d.name} bucket ${k}`))
  return t.toArray().map(r => { const o = r.toJSON(); return { key: o.key, gene_id: o.gene_id, name: o.name, chr: o.chr, tss: Number(o.tss) } })
})

/** Every gene whose `gene_id` or symbol is `id` (ASCII case ignored): two small range reads the
 *  first time, one per new bucket after. */
export async function lookupGenes(id: string): Promise<LookupRow[]> {
  const key = lookupKey(id)
  const d = await lookupDir()
  return (await lookupBucketRows(String(lookupBucket(key, d.nBuckets)))).filter(r => r.key === key)
}

/** The GWAS bin summary as rows (one whole-object read, for the landing track), or null when the
 *  experiment's GWAS has none. */
export const gwasBins = memo(async (): Promise<GwasBin[] | null> => {
  const s = await getStore()
  const b = s.experiment.gwas?.bins
  if (!b) return null
  return decodeGwasBins(await fetchObject(b.file), b, b.file)
})

/** The results file of one phenotype type on one chromosome, and the dof its slopes are rebuilt with. */
export function resultsFile(s: Store, phenotypeType: string, chr: string): { name: string; dof: number | null } {
  const r = s.results.get(phenotypeType)
  const name = r?.files[chr]
  const c = s.chroms.get(chr)
  if (!r || !name || !c) throw new Error(`experiment ${s.experiment.id}: no ${phenotypeType} results for ${chr}`)
  return { name, dof: r.dof }
}

/** A phenotype type's trans object and its dof, or null when it has none. */
export function transFile(s: Store, phenotypeType: string): { name: string; dof: number | null } | null {
  const r = s.results.get(phenotypeType)
  if (!r?.trans) return null
  return { name: r.trans.file, dof: r.dof }
}

/** A chromosome's GWAS file, or null when the GWAS has no rows there. */
export function gwasFile(s: Store, chr: string): { name: string } | null {
  const g = s.experiment.gwas, c = s.chroms.get(chr)
  const name = g?.files[chr]
  if (!g || !name || !c) return null
  return { name }
}

/** A catalog chromosome's variants file. */
export function variantsFile(s: Store, chr: string): { name: string } {
  const c = s.chroms.get(chr)
  if (!c) throw new Error(`catalog ${s.catalog.id}: no chromosome ${chr}`)
  return { name: c.file }
}

// ---- what the About, Home and gene pages print --------------------------------------------------

/** Whole-experiment counts, from the experiment's `counts` and the catalog. */
export interface Counts {
  genes_tested?: number; egenes?: number
  splice_phenotypes_tested?: number; sqtl_sig_phenotypes?: number; sqtl_sig_genes?: number
  variants_cis?: number; variants_trans_only?: number
  gwas_variants?: number; trans_pairs?: number
}
export interface PrecisionKind {
  neglog10p_max_error: number; slope_se_max_rel_error: number
  /** measured slope error in units of the row's SE; null when the set has no dof */
  slope_max_error_over_se: number | null
}
export interface Precision { af_max_error: number; eqtl: PrecisionKind; sqtl: PrecisionKind }

/** What the pages print about the release: the resolved store plus counts, sources and precision. */
export interface StoreInfo extends Store {
  counts: Counts
  sources: Record<string, { version: string; description: string }>
  precision: Precision | null
  /** annotation release label, e.g. "v34 (GRCh38.p13)" */
  annotationVersion: string | null
}

/** The adapter records a Zenodo release as the directory it unpacked (`zenodo_21382723`); cite it
 *  the way Zenodo itself does, as a DOI. Anything else is passed through as written. */
function zenodoDoi(s: string | undefined): string | undefined {
  const m = /^zenodo[_.]?(\d+)$/i.exec(s ?? '')
  return m ? `10.5281/zenodo.${m[1]}` : s
}

export const getStoreInfo = memo(async (): Promise<StoreInfo> => {
  const s = await getStore()
  // the builder's per-type counts (results.type_counts): no index read
  const ce = s.experiment.counts[EQTL_TYPE], cs = s.experiment.counts[SQTL_TYPE]
  const c: Counts = { genes_tested: ce?.with_rows ?? 0, egenes: ce?.significant ?? 0,
    splice_phenotypes_tested: cs?.phenotypes ?? 0, sqtl_sig_phenotypes: cs?.significant ?? 0, sqtl_sig_genes: cs?.significant_genes ?? 0 }
  c.variants_cis = s.catalog.chromosomes.reduce((a, x) => a + x.n_cis, 0)
  c.variants_trans_only = s.catalog.n_sites - c.variants_cis
  if (s.experiment.gwas) c.gwas_variants = s.experiment.gwas.n_rows
  if (s.experiment.trans) c.trans_pairs = s.experiment.trans.rows
  const e = s.results.get(EQTL_TYPE)?.precision, q = s.results.get(SQTL_TYPE)?.precision
  const precision = e && q ? {
    af_max_error: Math.max(e.af_max_error, q.af_max_error),
    eqtl: { neglog10p_max_error: e.neglog10p_max_error, slope_se_max_rel_error: e.slope_se_max_rel_error,
            slope_max_error_over_se: e.slope_max_error_over_se ?? null },
    sqtl: { neglog10p_max_error: q.neglog10p_max_error, slope_se_max_rel_error: q.slope_se_max_rel_error,
            slope_max_error_over_se: q.slope_max_error_over_se ?? null },
  } : null
  const version = s.annotation.source?.version ?? null
  const src = s.experiment.source?.source as { zenodo?: string } | undefined
  // no GWAS row: `gwas.id` names which set of the release was used, which the About text already
  // says in prose, with the case and control counts
  const sources: StoreInfo['sources'] = {
    gencode: { version: version ?? s.annotation.id, description: `gene annotation (${s.annotation.id})` },
    reference: { version: s.catalog.collection_digest, description: 'GRCh38 reference sequence collection (seqcol)' },
    catalog: { version: s.catalog.identity_digest, description: `variant catalog ${s.catalog.id}, ${s.catalog.n_sites.toLocaleString('en-US')} sites` },
    [s.experiment.id]: { version: zenodoDoi(src?.zenodo) ?? s.experiment.id, description: `experiment ${s.experiment.id} (store ${s.store.name})` },
  }
  return { ...s, counts: c, sources, precision, annotationVersion: version?.split(' ')[0] ?? null }
})
