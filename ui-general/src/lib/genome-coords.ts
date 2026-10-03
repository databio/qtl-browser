/** Linear-genome coordinate helpers for the genome track (ported from pegasus-v2f-ui). */
import type { ChromInfo, LocusCluster, TrackItem, TrackLocus } from '@/components/genome-track/types'

export function buildChromList(names: string[], lengths: number[]): ChromInfo[] {
  return names.map((name, i) => ({ name, length: lengths[i]! }))
}

/** Cumulative offset of each chromosome in the concatenated genome, plus total length. */
export function chromOffsets(chroms: ChromInfo[], gapBp = 50_000): { offsets: Map<string, number>; totalLength: number } {
  const offsets = new Map<string, number>()
  let pos = 0
  for (const chr of chroms) {
    offsets.set(chr.name, pos)
    pos += chr.length + gapBp
  }
  if (chroms.length > 0) pos -= gapBp
  return { offsets, totalLength: pos }
}

export function toAbsolute(chr: string, pos: number, offsets: Map<string, number>): number {
  const offset = offsets.get(chr)
  if (offset === undefined) throw new Error(`Unknown chromosome: ${chr}`)
  return offset + pos
}

export function fromAbsolute(absBp: number, chroms: ChromInfo[], offsets: Map<string, number>): { chr: string; pos: number } {
  for (let i = chroms.length - 1; i >= 0; i--) {
    const chr = chroms[i]!
    const offset = offsets.get(chr.name)!
    if (absBp >= offset) return { chr: chr.name, pos: absBp - offset }
  }
  return { chr: chroms[0]!.name, pos: 0 }
}

export function locusMidpoint(locus: TrackLocus, offsets: Map<string, number>): number {
  return (toAbsolute(locus.chr, locus.start, offsets) + toAbsolute(locus.chr, locus.end, offsets)) / 2
}

/** Greedy left-to-right merge of loci closer than `minPixelGap` on the same chromosome. */
export function clusterLoci(loci: TrackLocus[], bpToPixel: (chr: string, pos: number) => number, minPixelGap: number): TrackItem[] {
  if (loci.length === 0) return []
  const withPixel = loci.map(l => ({ locus: l, px: bpToPixel(l.chr, (l.start + l.end) / 2) }))
  withPixel.sort((a, b) => a.px - b.px)
  const result: TrackItem[] = []
  let group: typeof withPixel = [withPixel[0]!]
  for (let i = 1; i < withPixel.length; i++) {
    const last = group[group.length - 1]!
    if (withPixel[i]!.px - last.px < minPixelGap && withPixel[i]!.locus.chr === last.locus.chr) group.push(withPixel[i]!)
    else { result.push(finalizeGroup(group)); group = [withPixel[i]!] }
  }
  result.push(finalizeGroup(group))
  return result
}

function finalizeGroup(group: { locus: TrackLocus; px: number }[]): TrackItem {
  if (group.length === 1) return group[0]!.locus
  const loci = group.map(g => g.locus)
  const pixels = group.map(g => g.px)
  return {
    type: 'cluster', count: loci.length, loci, chr: loci[0]!.chr,
    start: Math.min(...loci.map(l => l.start)), end: Math.max(...loci.map(l => l.end)),
    centerPixel: (Math.min(...pixels) + Math.max(...pixels)) / 2,
  } satisfies LocusCluster
}

export function sortLociByPosition(loci: TrackLocus[], offsets: Map<string, number>): TrackLocus[] {
  return [...loci].sort((a, b) => locusMidpoint(a, offsets) - locusMidpoint(b, offsets))
}
