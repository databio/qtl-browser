/**
 * The qtlstore (SPEC.md sections 2, 3 and 8): `store.json` -> the experiments of one study and
 * their variant catalog, annotation and collection -> immutable objects by name. Every data URL in
 * the app is built here, so the only objects ever fetched are ones a pointer names.
 *
 * Pointers are fetched once per session with `no-cache` (they are mutable); objects are
 * content-addressed and read whole or by HTTP range. A range read sends no separate header request:
 * the pointer says what an object is, Store.validate checked every header against its pointer at
 * build time, and a header is checked wherever it arrives inside bytes read anyway (whole objects,
 * the gene lookup's directory, a hits file's frame table).
 *
 * One experiment is one scan (SPEC section 8): TOPCHeF is four -- eQTL and sQTL, cis and trans --
 * in one collection, and the DCM GWAS a fifth in its own. This app pins its study at build time:
 * `VITE_COLLECTION_OF` names one experiment, and every experiment in `store.json` sharing that
 * experiment's `collection` is the study; `VITE_GWAS` names the GWAS shown beside it, since the
 * two are different publications and the pairing is this deployment's choice, not a fact in the
 * store. A gene or variant page reads per-chromosome objects of each experiment it shows (SPEC
 * sections 6 and 8): its genes, exon models and search index parts, and one gene lookup bucket to
 * find a gene's chromosome. Nothing on a first gene page reads a whole-genome table; the gene list
 * and the search box read the whole-genome genes table and the cis experiments' search indexes the
 * first time they are used, and the browser cache keeps every object after.
 */
import { tableFromIPC, type Table } from 'apache-arrow'
import { checkNames, decodeArrowObject, decodeGwasBins, decodeLookupDir, DIGEST, FORMAT_VERSION, lookupBucket,
  lookupDirLen, lookupKey, OBJECT_NAME, parseFileHeader, type GwasBin } from './store-decode'
import type { SignificanceRule } from './significance'
import { beforeFetch, FAULT, faultHits } from './faults'

/** The one data host. An empty `VITE_DATA_BASE` means "same origin" (`/data`), which is what the
 *  local dev server and preview serve. */
export const DATA_BASE: string =
  (import.meta.env.VITE_DATA_BASE as string | undefined)?.replace(/\/$/, '') || `${window.location.origin}/data`

/** The experiment whose collection is this app's study: every experiment sharing its `collection`
 *  value in `store.json` is opened with it. One member rather than a digest, because a collection
 *  document's name is the digest of its bytes and changes when its title or citation is edited. */
export const COLLECTION_OF: string = (import.meta.env.VITE_COLLECTION_OF as string | undefined) || 'topchef_ge_cis'

/** The GWAS experiment shown beside the study, or '' for none. A different publication, so it is
 *  named rather than derived. */
export const GWAS_ID: string = (import.meta.env.VITE_GWAS as string | undefined) ?? 'dcm_jurgens2024_biobanks'

/** The phenotype types the gene page's two tabs read (TOPCHeF's `ge` and `leafcutter`). */
export const EQTL_TYPE = 'ge'
export const SQTL_TYPE = 'leafcutter'

// ---- pointer documents ------------------------------------------------------------------------

/** The variant catalogs' pointer directory, and the `store.json` key listing their ids. */
export const CATALOG_DIR = 'variant_catalogs'

export type ExperimentKind = 'qtl_cis' | 'qtl_trans' | 'gwas'

export interface StoreDoc {
  name: string; format_version: number; refget: string[]
  /** the seqcol digest every position in the store sits on; every catalog names it (SPEC section 2) */
  collection_digest: string | null
  [CATALOG_DIR]: string[]; annotations: string[]; overlaps?: string[]
  /** coloc lists, one per declared pair of experiments (SPEC section 20) */
  colocs?: string[]
  /** every experiment's kind and collection, generated from the documents (SPEC section 2) */
  experiments: Record<string, { kind: ExperimentKind; collection: string | null }>
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

/** A collection document (SPEC section 8): the publication a set of experiments came from,
 *  content-addressed, so every member shares one copy of the citation. */
export interface CollectionDoc {
  title: string
  citation: { doi?: string; title?: string; first_author?: string; last_author?: string; year?: number; venue?: string; [k: string]: unknown }
  source: Record<string, unknown>
  parent: string | null
}

/** One part of a search index (SPEC section 8): its object, row count and inclusive ord runs. */
export interface IndexPart { file: string; rows: number; ords: [number, number][] }
/** What the Home and About pages print for a cis scan. `significant` and `significant_genes`
 *  are **null** when the scan assessed no significance; zero would say the study found nothing. */
export interface TypeCounts {
  phenotypes: number; with_rows: number
  significant: number | null; significant_genes: number | null
}

/** How much one experiment's values were quantized (SPEC section 13): `theoretical_limit` is half
 *  a code step in each quantity's own unit, recomputable from any block header, and
 *  `measured_worst` the largest error found by decoding every row at build time. `neglog10p` is
 *  absolute, `beta_over_se` in units of the row's standard error, `se_rel` and `af` relative and
 *  absolute. */
export interface Precision {
  theoretical_limit: { neglog10p: number; beta_over_se: number; se_rel: number; af: number }
  measured_worst: { neglog10p: number; beta_over_se: number; se_rel: number; rows_compared: number }
}

interface ExperimentBase {
  id: string; kind: ExperimentKind
  /** the object name of the collection document, `<digest>.json` */
  collection: string
  catalog: string; catalog_identity: string
  allele_orientation_source: string | null
  precision: Precision | null
  source?: Record<string, unknown>
}

/** A cis QTL scan: blocks per chromosome, leads and credible sets in `hits`, a search index with
 *  one part per chromosome. */
export interface CisExperiment extends ExperimentBase {
  kind: 'qtl_cis'
  phenotype_type: string
  annotation: string; annotation_source: string | null
  /** provenance; nothing decodes through it */
  dof: number | null; dof_source: string | null
  /** which `permuted` column the study tested, and at what threshold (SPEC section 8).
   *  **Null when the source assessed no significance** -- see lib/significance.ts. */
  significance: SignificanceRule | null
  files: Record<string, string>
  hits: Record<string, string>
  search_index: string; n_phenotypes: number
  search_index_parts: Record<string, IndexPart>
  counts: TypeCounts
  unplaced: { count: number; examples: string[] }
}

/** A trans QTL scan: one object of frames, trans records in `hits`, a search index split by the
 *  phenotype's chromosome with an unplaced residue. */
export interface TransExperiment extends ExperimentBase {
  kind: 'qtl_trans'
  phenotype_type: string
  annotation: string; annotation_source: string | null
  dof: number | null; dof_source: string | null
  file: string; n_rows: number
  hits: Record<string, string>
  search_index: string; n_phenotypes: number
  search_index_parts: Record<string, IndexPart>
  /** the phenotypes nothing places on a chromosome (SPEC section 8); absent when there are none */
  search_index_unplaced?: IndexPart
  placement?: { rows: number; rows_in_source: number }
}

/** A GWAS: blocks per chromosome and the index into them, and the bin summary. No annotation, no
 *  hits, no search index (SPEC section 10). */
export interface GwasExperiment extends ExperimentBase {
  kind: 'gwas'
  title: string | null
  files: Record<string, string>; index: string
  n_rows: number; rows_by_chrom: Record<string, number>; block_rows: number; n_values: number[]
  /** the bin summary (SPEC section 10, "Bin summary"), null or absent when the store has none */
  bins?: { file: string; bin_bp: number; n_bins: number } | null
  source?: { file?: string; n_cases?: number; n_controls?: number; [k: string]: unknown }
}

export type QtlExperiment = CisExperiment | TransExperiment
export type ExperimentDoc = CisExperiment | TransExperiment | GwasExperiment

/** One study with everything its pointers name, resolved and cross-checked. */
export interface Store {
  store: StoreDoc
  catalog: CatalogDoc
  annotation: AnnotationDoc
  /** catalog chromosome table by name, with its 1-based ordinal */
  chroms: Map<string, CatalogChrom & { ordinal: number }>
  /** the study: its cis and trans scans by phenotype type, and every QTL experiment by id */
  cis: Map<string, CisExperiment>
  trans: Map<string, TransExperiment>
  experiments: Map<string, QtlExperiment>
  collectionName: string
  collection: CollectionDoc
  /** the GWAS shown beside the study, with its own collection; null when none is configured or
   *  the store holds none */
  gwas: GwasExperiment | null
  gwasCollection: CollectionDoc | null
  hasTrans: boolean
  hasGwas: boolean
}

// ---- fetching ---------------------------------------------------------------------------------

// the fault hooks (lib/faults.ts) sit behind `import.meta.env.DEV` literally, so a production
// build folds them and the module away
async function pointer<T>(path: string): Promise<T> {
  if (import.meta.env.DEV) {
    await beforeFetch()
    if (FAULT === 'store-404' && path === 'store.json') throw new Error(`${path}: HTTP 404`)
  }
  const r = await fetch(`${DATA_BASE}/${path}`, { cache: 'no-cache' })
  if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`)
  const doc = await r.json() as T
  return import.meta.env.DEV && FAULT !== 'none' ? faultPointer(path, doc) : doc
}

/** The document-level faults (lib/faults.ts): the document as a wrong build or a hand edit would leave it. */
function faultPointer<T>(path: string, doc: T): T {
  const d = doc as Record<string, unknown>
  if (path === 'store.json') {
    if (FAULT === 'store-version') return { ...d, format_version: 1 } as T
    if (FAULT === 'store-anchor') { const e = { ...(d.experiments as Record<string, unknown>) }; delete e[COLLECTION_OF]; return { ...d, experiments: e } as T }
  }
  if (FAULT === 'store-names' && path === `experiments/${COLLECTION_OF}.json`)
    return { ...d, files: { ...(d.files as Record<string, string>), chr1: 'hand-edited.qbe' } } as T
  return doc
}

/** The object-level faults: a missing object, or one whose bytes are not what its name promises. */
function faultObject(name: string, length: number | null): Uint8Array | null {
  if (!faultHits(name)) return null
  if (FAULT === 'object-404') throw new Error(`${name}: HTTP 404`)
  if (FAULT === 'object-corrupt') return new Uint8Array(length ?? 4096).fill(0x5a)
  return null
}

/** A full URL for an object a pointer names, `<digest>.<ext>` under `immutable/`. */
export function objectUrl(name: string): string {
  if (!OBJECT_NAME.test(name)) throw new Error(`"${name}" is not an object name (<digest>.<ext>)`)
  return `${DATA_BASE}/immutable/${name}`
}

/** A whole object. */
export async function fetchObject(name: string): Promise<Uint8Array> {
  if (import.meta.env.DEV) {
    await beforeFetch()
    const faulty = faultObject(name, null)
    if (faulty) return faulty
  }
  const r = await fetch(objectUrl(name))
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`)
  return new Uint8Array(await r.arrayBuffer())
}

/** A content-addressed JSON document (a collection). Cached forever by the browser, unlike a pointer. */
async function fetchJsonObject<T>(name: string): Promise<T> {
  if (import.meta.env.DEV) {
    await beforeFetch()
    if (faultObject(name, null)) return {} as T
  }
  const r = await fetch(objectUrl(name))
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`)
  return r.json() as Promise<T>
}

async function rangeFetch(name: string, offset: number, length: number): Promise<Uint8Array> {
  const range = `bytes=${offset}-${offset + length - 1}`
  if (import.meta.env.DEV) {
    await beforeFetch()
    const faulty = faultObject(name, length)
    if (faulty) return faulty
    if (FAULT === 'range-200' && faultHits(name)) throw new Error(`${name} ${range}: expected 206, got 200`)
  }
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

/** The ids a previous open found: the study's experiments, the GWAS, and their catalog and
 *  annotation. Kept per data host and anchor so a dev build pointed at another store cannot
 *  inherit ids from this one. */
interface Remembered { experiments: string[]; gwas: string | null; catalog: string; annotation: string }
const REMEMBER_KEY = `qtlstore:${DATA_BASE}:${COLLECTION_OF}:${GWAS_ID}`

function remembered(): Remembered | null {
  try {
    const r = JSON.parse(localStorage.getItem(REMEMBER_KEY) ?? 'null') as Remembered | null
    return Array.isArray(r?.experiments) && typeof r?.catalog === 'string' && typeof r?.annotation === 'string' ? r : null
  } catch { return null }   // no localStorage, or a value this build did not write: guess nothing
}

function remember(r: Remembered): void {
  try { localStorage.setItem(REMEMBER_KEY, JSON.stringify(r)) } catch { /* not worth failing an open over */ }
}

/** A pointer fetched on a guess, kept when the guess was the id the store turned out to name, else
 *  fetched fresh. A guess that 404s (the store was rebuilt) falls through to the real id. */
function speculative<T>(guessed: string[], path: (id: string) => string) {
  const started = new Map(guessed.map(id => [id, pointer<T>(path(id))]))
  for (const p of started.values()) p.catch(() => {})
  return async (id: string): Promise<T> => {
    const spec = started.get(id)
    if (spec) {
      const got = await spec.catch(() => null)
      if (got) return got
    }
    return pointer<T>(path(id))
  }
}

async function openStore(): Promise<Store> {
  // One wave where it can be. store.json is the only document whose name is known in advance;
  // the study's experiment ids, their catalog and annotation live in it and in the experiment
  // documents, which would cost two more round trips. So the ids a previous open recorded are
  // fetched in the same wave as store.json -- the normal case, since they change only when the
  // store is rebuilt -- and a wrong guess costs only the discarded requests. Pointer documents are
  // served `no-cache` with no validator, so each round trip is a full origin fetch.
  const guess = remembered()
  const expOf = speculative<ExperimentDoc>(guess ? [...guess.experiments, ...(guess.gwas ? [guess.gwas] : [])] : [], id => `experiments/${id}.json`)
  const catOf = speculative<CatalogDoc>(guess ? [guess.catalog] : [], id => `${CATALOG_DIR}/${id}.json`)
  const annOf = speculative<AnnotationDoc>(guess ? [guess.annotation] : [], id => `annotations/${id}.json`)
  const store = await pointer<StoreDoc>('store.json')
  if (store.format_version !== FORMAT_VERSION)
    throw new Error(`store.json: format version ${store.format_version}, and this build reads version ${FORMAT_VERSION}. The data was updated; reload the page.`)
  if (!store.experiments || typeof store.experiments !== 'object' || Array.isArray(store.experiments))
    throw new Error('store.json: experiments is not a map of id to {kind, collection}')
  const anchor = store.experiments[COLLECTION_OF]
  if (!anchor) throw new Error(`store.json lists no experiment "${COLLECTION_OF}"`)
  if (!anchor.collection) throw new Error(`experiment ${COLLECTION_OF} names no collection, so there is no study to open`)
  const collectionName = anchor.collection
  const ids = Object.entries(store.experiments).filter(([, e]) => e.collection === collectionName).map(([id]) => id).sort()
  for (const id of ids) {
    const k = store.experiments[id].kind
    if (k !== 'qtl_cis' && k !== 'qtl_trans') throw new Error(`experiment ${id} shares ${COLLECTION_OF}'s collection but is a ${k}; this app shows a QTL study`)
  }
  const gwasId = GWAS_ID && store.experiments[GWAS_ID] ? GWAS_ID : null
  if (GWAS_ID && !gwasId) console.warn(`store.json lists no GWAS experiment "${GWAS_ID}"; the GWAS panels will say so`)
  if (gwasId && store.experiments[gwasId].kind !== 'gwas') throw new Error(`experiment ${gwasId} is a ${store.experiments[gwasId].kind}, not a GWAS`)

  const [docs, gwasDoc, collection] = await Promise.all([
    Promise.all(ids.map(expOf)),
    gwasId ? expOf(gwasId) : Promise.resolve(null),
    fetchJsonObject<CollectionDoc>(collectionName),
  ])
  const qtl: QtlExperiment[] = []
  docs.forEach((d, i) => {
    if (d.id !== ids[i]) throw new Error(`experiments/${ids[i]}.json: id is "${d.id}"`)
    if (d.kind !== store.experiments[d.id].kind || d.collection !== collectionName)
      throw new Error(`experiments/${d.id}.json: kind or collection differ from store.json`)
    checkNames(d, `experiments/${d.id}.json`)
    qtl.push(d as QtlExperiment)
  })
  const catalogId = qtl[0]?.catalog, annotationId = qtl[0]?.annotation
  if (!catalogId || !annotationId) throw new Error(`the collection of ${COLLECTION_OF} holds no QTL experiment`)
  for (const d of qtl) if (d.catalog !== catalogId || d.annotation !== annotationId)
    throw new Error(`experiment ${d.id} is on catalog ${d.catalog} and annotation ${d.annotation}; this app shows one study on one catalog (${catalogId}, ${annotationId})`)
  let gwas: GwasExperiment | null = null
  if (gwasDoc) {
    if (gwasDoc.id !== gwasId || gwasDoc.kind !== 'gwas') throw new Error(`experiments/${gwasId}.json: id "${gwasDoc.id}" or kind ${gwasDoc.kind}`)
    if (gwasDoc.catalog !== catalogId) throw new Error(`GWAS ${gwasId} is on catalog ${gwasDoc.catalog}, the study on ${catalogId}`)
    checkNames(gwasDoc, `experiments/${gwasId}.json`)
    gwas = gwasDoc
  }
  if (!store[CATALOG_DIR]?.includes(catalogId) || !store.annotations?.includes(annotationId))
    throw new Error(`variant catalog ${catalogId} or annotation ${annotationId} is not in store.json`)
  const [catalog, annotation, gwasCollection] = await Promise.all([
    catOf(catalogId), annOf(annotationId),
    gwas ? fetchJsonObject<CollectionDoc>(gwas.collection) : Promise.resolve(null),
  ])
  // the two URLs can come from a remembered id, so check the documents are the ones the experiments
  // name rather than trusting the path they were fetched from
  if (catalog.id !== catalogId) throw new Error(`catalog ${catalogId} answered with id "${catalog.id}"`)
  if (annotation.id !== annotationId) throw new Error(`annotation ${annotationId} answered with id "${annotation.id}"`)
  for (const d of [...qtl, ...(gwas ? [gwas] : [])]) if (d.catalog_identity !== catalog.identity_digest)
    throw new Error(`experiment ${d.id}: catalog_identity ${d.catalog_identity} is not catalog ${catalog.id}'s ${catalog.identity_digest}`)
  if (catalog.orientation !== 'ref_alt') throw new Error(`catalog ${catalog.id}: orientation "${catalog.orientation}", this reader reads ref_alt`)
  if (!DIGEST.test(catalog.collection_digest)) throw new Error(`catalog ${catalog.id}: collection digest ${catalog.collection_digest}`)
  if (!annotation.chroms || !annotation.lookup) throw new Error(`annotation ${annotation.id} has no per-chromosome objects; rebuild the store with the current pipeline`)
  checkNames(catalog, `${CATALOG_DIR}/${catalog.id}.json`)
  checkNames(annotation, `annotations/${annotation.id}.json`)
  if (!collection?.title) throw new Error(`collection ${collectionName} has no title`)

  const chroms = new Map(catalog.chromosomes.map((c, i) => [c.name, { ...c, ordinal: i + 1 }]))
  const cis = new Map<string, CisExperiment>(), trans = new Map<string, TransExperiment>(), experiments = new Map<string, QtlExperiment>()
  for (const d of qtl) {
    const into = d.kind === 'qtl_cis' ? cis : trans
    if (into.has(d.phenotype_type)) throw new Error(`two ${d.kind} experiments of type ${d.phenotype_type} in one collection: ${into.get(d.phenotype_type)!.id} and ${d.id}`)
    if (!d.search_index_parts) throw new Error(`experiment ${d.id} has no search index parts; rebuild the store with the current pipeline`)
    for (const c of Object.keys(d.hits)) if (!chroms.has(c)) throw new Error(`experiment ${d.id}: hits for ${c}, not in the catalog`)
    for (const c of Object.keys(d.kind === 'qtl_cis' ? d.files : {})) if (!chroms.has(c)) throw new Error(`experiment ${d.id}: results for ${c}, not in the catalog`)
    ;(into as Map<string, QtlExperiment>).set(d.phenotype_type, d)
    experiments.set(d.id, d)
  }
  // only after every check: the next open guesses these ids and skips the round trips
  remember({ experiments: ids, gwas: gwasId, catalog: catalogId, annotation: annotationId })
  return { store, catalog, annotation, chroms, cis, trans, experiments, collectionName, collection, gwas, gwasCollection,
    hasTrans: trans.size > 0, hasGwas: gwas != null }
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

/** The whole-genome annotation genes table, one request: what the gene list and the search box
 *  read, instead of every chromosome's objects (gene-index.ts). */
export const wholeGenes = memo(async () => arrowTable((await getStore()).annotation.genes))

/** A QTL experiment's whole search index (one request), for the gene list and the search box. */
export const wholeExperimentIndex = memoBy(async (expId): Promise<Table> => {
  const e = (await getStore()).experiments.get(expId)
  if (!e) throw new Error(`no experiment ${expId} in the study`)
  return arrowTable(e.search_index)
})

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

/** One experiment's search index part for one chromosome (SPEC section 8), or null when it has none there. */
export const experimentPart = memoBy(async (key): Promise<Table | null> => {
  const [expId, chr] = key.split('\t')
  const e = (await getStore()).experiments.get(expId)
  if (!e) throw new Error(`no experiment ${expId} in the study`)
  const n = e.search_index_parts[chr]?.file
  return n ? arrowTable(n) : null
})
export const partKey = (expId: string, chr: string) => `${expId}\t${chr}`

/** A trans experiment's unplaced rows (phenotypes with no chromosome), or null when it has none. */
export const unplacedIndex = memoBy(async (expId): Promise<Table | null> => {
  const e = (await getStore()).experiments.get(expId)
  if (!e || e.kind !== 'qtl_trans' || !e.search_index_unplaced) return null
  return arrowTable(e.search_index_unplaced.file)
})

/** The part of one experiment's index holding row `ord`: a chromosome name, `'unplaced'` for a trans
 *  scan's residue, or undefined when no part lists it. Ordinals are per experiment (SPEC section 8). */
export function partOfOrd(e: QtlExperiment, ord: number): string | undefined {
  const has = (p: IndexPart) => p.ords.some(([a, b]) => ord >= a && ord <= b)
  for (const [c, p] of Object.entries(e.search_index_parts)) if (has(p)) return c
  const u = e.kind === 'qtl_trans' ? e.search_index_unplaced : undefined
  return u && has(u) ? 'unplaced' : undefined
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
 *  GWAS has none. */
export const gwasBins = memo(async (): Promise<GwasBin[] | null> => {
  const s = await getStore()
  const b = s.gwas?.bins
  if (!b) return null
  return decodeGwasBins(await fetchObject(b.file), b, b.file)
})

/** A coloc list (SPEC section 20): which phenotypes of one QTL experiment someone found to
 *  colocalize with one GWAS experiment, and who said so. Membership only; PP.H4 and the lead are
 *  computed from the two experiments' rows when a locus is opened. */
export interface ColocDoc {
  id: string; qtl: string; gwas: string; threshold: string
  source: Record<string, string>; note: string | null; object: string; n: number
}
/** One list row: a gene, and the phenotype when the list names one (null at gene level). */
export interface ColocRow { gene_id: string; phenotype_id: string | null }
export interface ColocList { doc: ColocDoc; rows: ColocRow[]; experiment: QtlExperiment }

/** The coloc lists for the open study's QTL experiments against its GWAS, with their rows; empty
 *  when the store holds none for those pairs. Every listed pointer is read to find the pair's, which
 *  is a handful today; a store with many would want the pair in `store.json`. */
export const colocLists = memo(async (): Promise<ColocList[]> => {
  const s = await getStore()
  if (!s.gwas) return []
  const gwasId = s.gwas.id
  const docs = await Promise.all((s.store.colocs ?? []).map(id => pointer<ColocDoc>(`colocs/${id}.json`)))
  return Promise.all(docs.filter(d => d.gwas === gwasId && s.experiments.has(d.qtl)).map(async doc => {
    const where = `colocs/${doc.id}.json`
    if (!OBJECT_NAME.test(doc.object) || !doc.object.endsWith('.json')) throw new Error(`${where}: bad object name ${doc.object}`)
    const rows = await fetchJsonObject<ColocRow[]>(doc.object)
    if (!Array.isArray(rows) || rows.length !== doc.n) throw new Error(`${where}: the list holds ${Array.isArray(rows) ? rows.length : 'no'} rows, the pointer says ${doc.n}`)
    return { doc, rows, experiment: s.experiments.get(doc.qtl)! }
  }))
})

/** The cis scan of one phenotype type, or null when the study has none. */
export const cisOf = (s: Store, phenotypeType: string): CisExperiment | null => s.cis.get(phenotypeType) ?? null
/** The trans scan of one phenotype type, or null. */
export const transOf = (s: Store, phenotypeType: string): TransExperiment | null => s.trans.get(phenotypeType) ?? null

/** The results file of one phenotype type's cis scan on one chromosome. */
export function resultsFile(s: Store, phenotypeType: string, chr: string): { name: string; experiment: CisExperiment } {
  const e = s.cis.get(phenotypeType)
  const name = e?.files[chr]
  const c = s.chroms.get(chr)
  if (!e || !name || !c) throw new Error(`no ${phenotypeType} cis results for ${chr}`)
  return { name, experiment: e }
}

/** A phenotype type's trans object, with the dof its r² is formed with, or null when the study
 *  has no trans scan of that type. */
export function transFile(s: Store, phenotypeType: string): { name: string; dof: number | null; experiment: TransExperiment } | null {
  const e = s.trans.get(phenotypeType)
  if (!e) return null
  return { name: e.file, dof: e.dof, experiment: e }
}

/** An experiment's hits file for one chromosome. */
export function hitsFile(s: Store, e: QtlExperiment, chr: string): string {
  const name = e.hits[chr]
  if (!name || !s.chroms.has(chr)) throw new Error(`experiment ${e.id}: no hits file for ${chr}`)
  return name
}

/** A chromosome's GWAS file, or null when the GWAS has no rows there. */
export function gwasFile(s: Store, chr: string): { name: string } | null {
  const g = s.gwas, c = s.chroms.get(chr)
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

/** Whole-study counts, from the cis scans' `counts`, the trans scans and the catalog. */
export interface Counts {
  /** the significance-derived counts are null when the scan assessed none; a page omits them
   *  rather than printing 0, which would say the study found nothing */
  genes_tested?: number; egenes?: number | null
  splice_phenotypes_tested?: number; sqtl_sig_phenotypes?: number | null; sqtl_sig_genes?: number | null
  variants_cis?: number; variants_trans_only?: number
  gwas_variants?: number; trans_pairs?: number
}

/** The worst `theoretical_limit` over a set of experiments, each quantity in its own unit; null
 *  when none of them carries one. */
export type Bounds = Precision['theoretical_limit']
export interface PrecisionSummary { cis: Bounds | null; trans: Bounds | null; gwas: Bounds | null }

/** A citation the About page can print: the collection's title and the paper it cites. */
export interface Citation { title: string; doi: string | null; venue: string | null; year: number | null; first_author: string | null }

/** What the pages print about the release: the resolved store plus counts, sources and precision. */
export interface StoreInfo extends Store {
  counts: Counts
  sources: Record<string, { version: string; description: string }>
  precision: PrecisionSummary
  /** annotation release label, e.g. "v34 (GRCh38.p13)" */
  annotationVersion: string | null
  study: Citation
  gwasStudy: Citation | null
}

function worstBounds(es: { precision: Precision | null }[]): Bounds | null {
  const ps = es.map(e => e.precision?.theoretical_limit).filter((p): p is Bounds => p != null)
  if (!ps.length) return null
  const w = (k: keyof Bounds) => Math.max(...ps.map(p => p[k]))
  return { neglog10p: w('neglog10p'), beta_over_se: w('beta_over_se'), se_rel: w('se_rel'), af: w('af') }
}

function citationOf(c: CollectionDoc): Citation {
  const ct = c.citation ?? {}
  return { title: c.title, doi: typeof ct.doi === 'string' ? ct.doi : null, venue: typeof ct.venue === 'string' ? ct.venue : null,
    year: typeof ct.year === 'number' ? ct.year : null, first_author: typeof ct.first_author === 'string' ? ct.first_author : null }
}

export const getStoreInfo = memo(async (): Promise<StoreInfo> => {
  const s = await getStore()
  // the builder's per-scan counts (results.type_counts): no index read
  const ce = s.cis.get(EQTL_TYPE)?.counts, cs = s.cis.get(SQTL_TYPE)?.counts
  const c: Counts = { genes_tested: ce?.with_rows ?? 0, egenes: ce?.significant ?? null,
    splice_phenotypes_tested: cs?.phenotypes ?? 0, sqtl_sig_phenotypes: cs?.significant ?? null,
    sqtl_sig_genes: cs?.significant_genes ?? null }
  c.variants_cis = s.catalog.chromosomes.reduce((a, x) => a + x.n_cis, 0)
  c.variants_trans_only = s.catalog.n_sites - c.variants_cis
  if (s.gwas) c.gwas_variants = s.gwas.n_rows
  if (s.trans.size) c.trans_pairs = [...s.trans.values()].reduce((a, e) => a + e.n_rows, 0)
  const precision: PrecisionSummary = { cis: worstBounds([...s.cis.values()]), trans: worstBounds([...s.trans.values()]),
    gwas: s.gwas ? worstBounds([s.gwas]) : null }
  const version = s.annotation.source?.version ?? null
  const study = citationOf(s.collection), gwasStudy = s.gwasCollection ? citationOf(s.gwasCollection) : null
  const sources: StoreInfo['sources'] = {
    gencode: { version: version ?? s.annotation.id, description: `gene annotation (${s.annotation.id})` },
    reference: { version: s.catalog.collection_digest, description: 'GRCh38 reference sequence collection (seqcol)' },
    catalog: { version: s.catalog.identity_digest, description: `variant catalog ${s.catalog.id}, ${s.catalog.n_sites.toLocaleString('en-US')} sites` },
    study: { version: study.doi ?? s.collectionName, description: `${study.title}; experiments ${[...s.experiments.keys()].join(', ')} (store ${s.store.name})` },
    ...(gwasStudy && s.gwas ? { gwas: { version: gwasStudy.doi ?? s.gwas.collection, description: `${gwasStudy.title}; experiment ${s.gwas.id}` } } : {}),
  }
  return { ...s, counts: c, sources, precision, annotationVersion: version?.split(' ')[0] ?? null, study, gwasStudy }
})
