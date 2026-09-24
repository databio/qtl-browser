import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'
import { defineConfig, type Plugin } from 'vite'
import { fileURLToPath, URL } from 'node:url'
import { createReadStream, statSync } from 'node:fs'
import { join, normalize } from 'node:path'
import type { IncomingMessage, ServerResponse } from 'node:http'

/** The local qtlstore (store.json, pointers, immutable/) served at /data. `QTL_DATA_DIR` points it
 *  at a store anywhere on disk, e.g. a copy of a Rivanna build in /tmp. */
const DATA_DIR = process.env.QTL_DATA_DIR || fileURLToPath(new URL('../data/store', import.meta.url))

/**
 * Serve DATA_DIR at /data in `vite` and `vite preview`, with HTTP Range support, the way B2
 * (cloud2.databio.org) serves it in production. This replaces a public/ symlink, which `vite build` would copy
 * wholesale into dist/.
 */
function serveDerivedData(): Plugin {
  const handler = (req: IncomingMessage, res: ServerResponse, next: () => void) => {
    const url = req.url ?? ''
    if (!url.startsWith('/data/')) return next()
    const rel = normalize(decodeURIComponent(url.slice('/data/'.length).split('?')[0]))
    if (rel.startsWith('..')) { res.statusCode = 403; return res.end() }
    const file = join(DATA_DIR, rel)
    let size: number
    try { size = statSync(file).size } catch { res.statusCode = 404; return res.end() }
    const type = file.endsWith('.json') ? 'application/json' : 'application/octet-stream'
    res.setHeader('Accept-Ranges', 'bytes')
    res.setHeader('Content-Type', type)
    // Chromium stores a 206 only when the response carries a strong validator, so without one `warm`
    // re-downloads every range. Size and mtime identify a local file. Production (cloud2.databio.org)
    // must send one too; see the deploy notes in README.md.
    const st = statSync(file)
    res.setHeader('ETag', `"${st.size.toString(16)}-${Math.floor(st.mtimeMs).toString(16)}"`)
    res.setHeader('Last-Modified', new Date(st.mtimeMs).toUTCString())
    // the production headers (SPEC section 3, set at upload to B2): an immutable/ name
    // carries the file's content hash, so a rebuilt file is a new URL and the old one can be cached
    // forever; everything else revalidates. Harness cold runs use a fresh browser context.
    res.setHeader('Cache-Control', rel.startsWith('immutable/') ? 'public, max-age=31536000, immutable' : 'no-cache')
    const range = /^bytes=(\d*)-(\d*)$/.exec(req.headers.range ?? '')
    if (req.headers.range !== undefined && !range) {
      // multi-range or malformed: refuse rather than fall through to the whole file
      res.statusCode = 416; res.setHeader('Content-Range', `bytes */${size}`); return res.end()
    }
    if (range) {
      const start = range[1] ? Number(range[1]) : Math.max(0, size - Number(range[2]))
      const end = range[1] && range[2] ? Math.min(Number(range[2]), size - 1) : range[1] ? size - 1 : size - 1
      if (start > end || start >= size) { res.statusCode = 416; res.setHeader('Content-Range', `bytes */${size}`); return res.end() }
      res.statusCode = 206
      res.setHeader('Content-Range', `bytes ${start}-${end}/${size}`)
      res.setHeader('Content-Length', String(end - start + 1))
      if (req.method === 'HEAD') return res.end()
      return createReadStream(file, { start, end }).pipe(res)
    }
    res.statusCode = 200
    res.setHeader('Content-Length', String(size))
    if (req.method === 'HEAD') return res.end()
    createReadStream(file).pipe(res)
  }
  return {
    name: 'serve-derived-data',
    configureServer(server) { server.middlewares.use(handler) },
    configurePreviewServer(server) { server.middlewares.use(handler) },
  }
}

export default defineConfig({
  plugins: [react(), tailwindcss(), serveDerivedData()],
  resolve: { alias: { '@': fileURLToPath(new URL('./src', import.meta.url)) } },
  // duckdb-wasm ships its own workers and wasm; pre-bundling breaks the worker URLs
  optimizeDeps: { exclude: ['@duckdb/duckdb-wasm'] },
  build: { target: 'es2022' },
})
