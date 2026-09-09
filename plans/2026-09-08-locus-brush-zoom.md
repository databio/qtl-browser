---
date: 2026-09-08
status: in-progress
model: Claude Fable 5.1
description: Transient magnifier on the locus scatter - dragging a brush shows the brushed slice in a floating detail popup that vanishes on release; LocusCompare and the gene track shade follow the brush while it exists
---

# Locus brush magnifier

## Why

The locus scatter shows a fixed ±1 Mb window. Dense loci pile hundreds of points into a few
pixels around the lead, and there is no way to look closer. Pan-zoom on the scatter itself
would fight page scrolling and lose the whole-window context. A brush that exists only while
the mouse is down, with a floating detail view of the brushed slice, gives a closer look
without changing the page: nothing moves, nothing persists, release restores everything.
Mosaic drives the detail, the LocusCompare panel, and the gene-track shade from one Selection,
and the data is already in the browser as the materialized locus table, so no fetch happens.

## Change

All in `ui/`. No pipeline or R2 change.

**`components/LocusPlot.tsx`**

- A second Selection, `brush` (`Selection.intersect()`), created with the hover Selection when
  a locus materializes and replaced on locus change.
- The overview scatter gains `vg.intervalX({ as: brush })` right after the data dots. The
  overview's own dots are not filtered by the brush (Mosaic skips a clause for the marks that
  produced it).
- React mirrors the brush's `value` event, rAF-throttled, into `brushed: [lo, hi] | null`.
- While `brushed` is set: the hover Selection is reset once, and every plot's SVG gets
  `pointer-events: none` via a Tailwind arbitrary variant on the root, so the nearest-variant
  tooltip stops chasing the cursor. d3-brush tracks the drag through window listeners, so the
  brush itself is unaffected.
- Pointer-down on the scatter records its screen rect. While `brushed` is set, a React portal
  into `document.body` renders a fixed-position popup the width of the scatter, above it when
  there is room and below otherwise, with `pointer-events: none` so it never intercepts the
  drag. It holds `LocusDetail`: the same dot encoding reading `vg.from(table, { filterBy:
  brush })` with `vg.xDomain(brush)`, so Mosaic re-queries and re-scales it as the brush moves
  without React redrawing. No hover marks; it exists only during the drag.
- On window `pointerup` or `pointercancel` while brushed, a zero-delay timeout resets the
  interval interactor (clears the brush graphic) and the Selection (empties the value). The
  delay lets d3-brush's own mouseup handling, which follows pointerup, run first so it cannot
  re-publish the final extent after the clear. The popup unmounts, the shade clears, and
  LocusCompare returns to the whole window.
- The plot element exposes vgplot's Plot instance on `.value`; the draw effect keeps the brush
  interactor from its `interactors` list for the reset.
- Header hint "drag on the plot to magnify" replaces nothing; there is no reset button since
  nothing persists.

**`components/GeneTrack.tsx`**: new `shade` prop; when set, a low-alpha rect over the brushed
interval spanning all lanes, drawn under everything else.

**`components/LocusCompare.tsx`**: new `brush` prop; the dot mark reads
`vg.from(table, { filterBy: brush })`. The "no GWAS overlap" empty state still uses the
unfiltered count, so it describes the window rather than the brush.

Nothing on the gene page or in the cis table changes.

## Decisions & ownership

| # | Decision | Owner | Note |
|---|---|---|---|
| D1 | Brush on the overview rather than pan-zoom | user-owned, surfaced | Sam chose it after the pan-zoom wheel-capture concern |
| D2 | Transient: brush and detail exist only while the mouse is down | user-owned, surfaced | Sam asked for this after trying a persistent in-flow slot (first build); the slot pushed the page down while LocusCompare sat still, which read as jarring |
| D3 | Detail floats in a portal above the scatter, not in the layout | user-owned, surfaced | with a transient view a layout shift on every drag is worse than a floating panel; fixed positioning from the pointer-down rect, no scroll tracking needed since a drag does not scroll |
| D4 | Gene track shades the brushed interval | user-owned, surfaced | yes |
| D5 | No gene track under the detail | user-owned, surfaced | no |
| D6 | Cis table does not follow the brush | user-owned, taken (reverted) | the first build filtered the table; with a brush that vanishes on release there is nothing for a table to follow, and per-frame page queries that revert are churn. Removed without asking since it falls out of D2; flagging it here |
| D7 | Hover tooltip muted while brushing by disabling pointer events on the plot SVGs | AI-owned, defended | the nearest interactor listens for pointermove on the SVG and has no button check; d3-brush does not need the SVG to receive events mid-drag. The class is applied on the first brush event, not on pointer-down, because applying it before the compatibility mousedown fires would keep d3-brush from starting |
| D8 | Release resets via the interactor instance plus a deferred Selection reset | AI-owned, defended | the interactor does not watch its Selection; resetting the Selection alone leaves the rectangle. The defer avoids d3-brush's mouseup re-publishing the extent |
| D9 | Detail plot binds its x domain to the Selection | AI-owned, defended | vgplot re-sets a Param-valued attribute on change; React only mounts and unmounts the popup |
| D10 | LocusCompare's empty state describes the whole window | AI-owned, default | a brushed slice with no GWAS overlap shows an empty square rather than a message |
| D11 | No detail on touch | AI-owned, default | pointer-down anchors and window pointerup work for touch too, but a one-finger drag also scrolls; untested |

## What this changes elsewhere

- **Nothing moves in the layout.** The popup is fixed-position in a portal. Export targets are
  unchanged and never include the popup, since it only exists mid-drag.
- **Hover is dead while brushing**, in all three plots, by design.
- **Query load.** Each brush move issues a few small queries against an in-memory table of a few
  thousand rows: detail dots and LocusCompare dots. Well under a frame each. The gene track
  re-renders per frame from React state (rAF-throttled), which is a small SVG rebuild.
- **The popup can cover the section header or tab bar** when the scatter is near the top of the
  viewport and the popup is placed above it. It is transient and pointer-transparent.
- **No pipeline, data, or R2 change.**

## Implementation log (2026-09-08)

First build: persistent brush with an in-flow detail slot above the overview, cis table
following the brush. Sam's check: the layout shift with LocusCompare staying put was jarring,
and the hover tooltip chased the cursor during the drag. Replaced the same day with the
transient popup design above; the cis table filter and its query plumbing were removed. Mosaic
internals confirmed before coding: vgplot attribute setters subscribe to Param-valued arguments,
so `vg.xDomain(brush)` re-scales the detail on each brush event; `intervalX({ as })` maps `as`
to the interactor's selection; the interval clause names the brushed plot's marks as its
clients, so the overview's own dots are not filtered by it; the nearest interactor listens on
`pointerenter pointerdown pointermove` on the SVG root with no button check. Verified with
`npx tsc -b`.

Parked the same evening at Sam's request, before a browser check of the transient version:
`BRUSH_MAGNIFIER = false` at the top of `LocusPlot.tsx` leaves out the brush interactor, the
header hint, and the popup, so the Selection never gets a value and the gene-track shade and
LocusCompare filter stay inert. All the code remains and type-checks. Flip the constant to
resume.
