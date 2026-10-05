import { useMemo } from 'react'
import type { ChromInfo, ViewState } from './types'
import { buildChromList, chromOffsets, fromAbsolute, toAbsolute } from '@/lib/genome-coords'

export type TrackLayout = {
  bpToPixel: (chr: string, pos: number) => number
  pixelToBp: (px: number) => { chr: string; pos: number }
  visibleChroms: { chr: string; startPx: number; widthPx: number }[]
  chroms: ChromInfo[]
  totalLength: number
  offsets: Map<string, number>
  bpPerPixel: number
}

/** Pixel positions for the current view. */
export function useTrackLayout(containerWidth: number, view: ViewState, names: string[], lengths: number[]): TrackLayout {
  const chroms = useMemo(() => buildChromList(names, lengths), [names, lengths])
  const { offsets, totalLength } = useMemo(() => chromOffsets(chroms), [chroms])
  return useMemo(() => {
    const viewSpan = view.endBp - view.startBp
    const bpPerPixel = viewSpan / containerWidth
    const pxPerBp = containerWidth / viewSpan
    const bpToPixel = (chr: string, pos: number) => (toAbsolute(chr, pos, offsets) - view.startBp) * pxPerBp
    const pixelToBp = (px: number) => fromAbsolute(view.startBp + px / pxPerBp, chroms, offsets)
    const visibleChroms: TrackLayout['visibleChroms'] = []
    for (const chr of chroms) {
      const chrStart = offsets.get(chr.name)!
      const chrEnd = chrStart + chr.length
      const buffer = viewSpan * 0.5
      if (chrEnd < view.startBp - buffer || chrStart > view.endBp + buffer) continue
      const startPx = (chrStart - view.startBp) * pxPerBp
      visibleChroms.push({ chr: chr.name, startPx, widthPx: (chrEnd - view.startBp) * pxPerBp - startPx })
    }
    return { bpToPixel, pixelToBp, visibleChroms, chroms, totalLength, offsets, bpPerPixel }
  }, [containerWidth, view, chroms, offsets, totalLength])
}
