import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { BrowserRouter } from 'react-router'
import { ThemeProvider } from '@/contexts/theme-context'
import { ManifestProvider } from '@/contexts/manifest-context'
import '@/app.css'
import App from '@/App'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <ThemeProvider>
      <ManifestProvider>
        <BrowserRouter>
          <App />
        </BrowserRouter>
      </ManifestProvider>
    </ThemeProvider>
  </StrictMode>,
)
