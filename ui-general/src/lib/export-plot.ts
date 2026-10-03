/**
 * Export the SVG plots inside a container as one standalone SVG or PNG file. Several plots
 * (e.g. the locus scatter and the gene track beneath it) are stacked vertically. Observable
 * Plot keeps its own styles in a <style> element inside each SVG, so serialization keeps the
 * look; the font and text color it inherits from the page are written onto the root.
 */
export function collectSVG(container: HTMLElement, background: string): { svg: string; width: number; height: number } {
  const svgs = Array.from(container.querySelectorAll<SVGSVGElement>('svg'))
  const cs = getComputedStyle(container)
  const width = Math.max(...svgs.map(s => Number(s.getAttribute('width')) || s.getBoundingClientRect().width))
  let y = 0
  const parts: string[] = []
  for (const s of svgs) {
    const h = Number(s.getAttribute('height')) || s.getBoundingClientRect().height
    const w = Number(s.getAttribute('width')) || s.getBoundingClientRect().width
    const clone = s.cloneNode(true) as SVGSVGElement
    clone.removeAttribute('class')
    clone.setAttribute('x', '0'); clone.setAttribute('y', String(y))
    clone.setAttribute('width', String(w)); clone.setAttribute('height', String(h))
    clone.setAttribute('overflow', 'visible')
    parts.push(new XMLSerializer().serializeToString(clone))
    y += h
  }
  const font = cs.fontFamily || 'system-ui, sans-serif'
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${width}" height="${y}" viewBox="0 0 ${width} ${y}" ` +
    `font-family="${font.replace(/"/g, "'")}" font-size="11" color="${cs.color}" style="background:${background}">` +
    `<rect width="100%" height="100%" fill="${background}"/>${parts.join('')}</svg>`
  return { svg, width, height: y }
}

function download(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url; a.download = name; a.click()
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export function exportSVG(container: HTMLElement, name: string, background: string) {
  const { svg } = collectSVG(container, background)
  download(new Blob([svg], { type: 'image/svg+xml' }), `${name}.svg`)
}

export async function exportPNG(container: HTMLElement, name: string, background: string, scale = 2) {
  const { svg, width, height } = collectSVG(container, background)
  const img = new Image()
  const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' }))
  await new Promise<void>((resolve, reject) => { img.onload = () => resolve(); img.onerror = () => reject(new Error('svg render failed')); img.src = url })
  const canvas = document.createElement('canvas')
  canvas.width = Math.round(width * scale); canvas.height = Math.round(height * scale)
  const ctx = canvas.getContext('2d')!
  ctx.scale(scale, scale)
  ctx.drawImage(img, 0, 0)
  URL.revokeObjectURL(url)
  canvas.toBlob(b => { if (b) download(b, `${name}.png`) }, 'image/png')
}
