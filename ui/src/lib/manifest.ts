/**
 * `manifest.json`: the data base URL, one fetch per session, and the shapes the app reads out of it.
 * The manifest context and the pack reader share the same promise, and every data URL in the app is
 * built here, so the file names the manifest lists are the only ones ever fetched.
 */
import { PACK_FORMAT, PACK_VERSION } from './pack-decode'

/** The one data host. An empty `VITE_DATA_BASE` means "same origin", which is what the local
 *  preview and the bench harness build against. */
export const DATA_BASE: string =
  (import.meta.env.VITE_DATA_BASE as string | undefined)?.replace(/\/$/, '') || `${window.location.origin}/data`

/** A full URL for a path the manifest names, such as `immutable/eqtl.chr1.c6c916eec090a128.qbe`. */
export const dataUrl = (path: string) => `${DATA_BASE}/${path}`

/** Rehearsal only: point a build at a staged copy of the manifest (`immutable/manifest.<sha16>.json`)
 *  so the bucket can be tested before `manifest.json` itself is switched over. */
const MANIFEST_PATH = (import.meta.env.VITE_MANIFEST as string | undefined) || 'manifest.json'

/** The `packs` block (SPEC.md section 3): every pack path and dof the gene page reads. Paths are
 *  content-addressed, `immutable/<stem>.<sha16>.<ext>`, and are only ever read from here. */
export interface PackManifest {
  format: string
  version: number
  dof: { eqtl: number; sqtl: number }
  variant_page_size: number
  variant_page_codec: string
  search_index: string
  gwas_index: string
  gwas_block_rows: number
  /** SPEC.md sections 14 and 15: the two whole-file reads the variant page needs */
  rsid_index: string
  rsid_block_records: number
  variant_index: string
  hits_frame_variants: number
  /** chromosome -> path under DATA_BASE; `gwas` lists only chromosomes with GWAS rows, `trans` only
   *  gene chromosomes with trans rows (chrM included) */
  files: { variants: Record<string, string>; eqtl: Record<string, string>; sqtl: Record<string, string>; trans: Record<string, string>
    gwas: Record<string, string>; hits: Record<string, string> }
  bytes: { variants: number; eqtl: number; sqtl: number; trans: number; gwas: number; gwas_index: number
    hits: number; rsid_index: number; variant_index: number }
  counts: {
    variants: { cis: number; trans_only: number }
    eqtl: { genes: number; rows: number; memberships: number }
    sqtl: { introns: number; rows: number; memberships: number }
    gwas: { rows: number }
    trans: { genes: number; rows: number; eqtl_rows: number; sqtl_rows: number }
    hits: { rows: number; frames: number; trans_eqtl: number; trans_sqtl: number; lead_eqtl: number; lead_sqtl: number
      cs_eqtl: number; cs_sqtl: number }
    rsid_index: { records: number }
  }
}

/** The worst rounding error of one QTL type's stored per-variant values (SPEC.md section 9). */
export interface PrecisionKind {
  /** absolute, on -log10 p */
  neglog10p_max_error: number
  /** the rebuilt slope, in units of the row's standard error */
  slope_max_error_over_se: number
  /** relative, on the standard error */
  slope_se_max_rel_error: number
}

/** `precision`: what `packcheck roundtrip` measured over every row. The About page and the cis
 *  table's note are written from it, so the site never claims more accuracy than was measured. */
export interface Precision {
  /** absolute, on allele frequency */
  af_max_error: number
  eqtl: PrecisionKind
  sqtl: PrecisionKind
  /** the report the numbers came from */
  source: string
}

/** Whole-dataset counts, all optional so a manifest from an older build still renders. */
export interface Counts {
  genes_tested?: number
  egenes?: number
  splice_phenotypes_tested?: number
  sqtl_sig_phenotypes?: number
  sqtl_sig_genes?: number
  variants_cis?: number
  variants_trans_only?: number
  gwas_variants?: number
  trans_pairs?: number
  trans_variants?: number
  rsid_match?: { exact: number; position: number; none: number }
}

/** Shape of data/derived/manifest.json as the pipeline's manifest step writes it. The `immutable`
 *  block is left out on purpose: only `pipeline/upload.py` reads it. */
export interface Manifest {
  built: string
  pipeline_commit: string | null
  significance_rule: string
  counts: Counts
  sources: Record<string, { version: string; description: string }>
  packs: PackManifest
  precision?: Precision
  assets?: Record<string, { path: string; bytes: number }>
  gwas_dcm: { file: string; n_cases: number; n_controls: number; variants: number } | null
}

let manifest: Promise<Manifest> | null = null

/** Memoized; a failed fetch is forgotten so the next call retries. `no-cache` revalidates every
 *  time, so a copy the browser kept on its own never pairs new code with an older build's files. */
export function getManifest(): Promise<Manifest> {
  if (!manifest) {
    const p = fetch(dataUrl(MANIFEST_PATH), { cache: 'no-cache' })
      .then(r => {
        if (!r.ok) throw new Error(`${MANIFEST_PATH}: HTTP ${r.status}`)
        return r.json() as Promise<Manifest>
      })
      .then(m => {
        // a reader supports exactly the format version the manifest names (SPEC section 3)
        const pk = m.packs
        if (pk?.format !== PACK_FORMAT || pk.version !== PACK_VERSION)
          throw new Error(`${MANIFEST_PATH}: packs are ${pk?.format} version ${pk?.version}, and this build reads ` +
            `${PACK_FORMAT} version ${PACK_VERSION}. The data was updated; reload the page.`)
        return m
      })
    manifest = p
    p.catch(() => { if (manifest === p) manifest = null })
  }
  return manifest
}
