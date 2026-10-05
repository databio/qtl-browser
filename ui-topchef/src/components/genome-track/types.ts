/** Chromosome info from the seqcol API. */
export type ChromInfo = {
  name: string   // "chr1" … "chrX", "chrY"
  length: number // bp
}

/** A locus positioned on the genome track. */
export type TrackLocus = {
  id: string
  chr: string
  start: number
  end: number
  label: string
  /** Series key for marker color (e.g. "eQTL", "sQTL", "both"). */
  trait?: string
}

/** Nearby loci merged into one marker when zoomed out. */
export type LocusCluster = {
  type: 'cluster'
  count: number
  loci: TrackLocus[]
  chr: string
  start: number
  end: number
  centerPixel: number
}

/** A fixed genomic window drawn as a downward bar beneath the chromosome bar. */
export type TrackBin = {
  chr: string
  start: number
  end: number
  value: number     // bar height in data units (e.g. -log10 p)
  sig?: boolean     // highlighted (e.g. holds a genome-wide significant variant)
  title?: string    // hover readout
}

/** Viewport in absolute bp (chromosomes concatenated). */
export type ViewState = { startBp: number; endBp: number }

export type TrackItem = TrackLocus | LocusCluster

export function isCluster(item: TrackItem): item is LocusCluster {
  return 'type' in item && item.type === 'cluster'
}
