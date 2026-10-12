import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { BrowserRouter } from 'react-router'
import { ThemeProvider } from '@/contexts/theme-context'
import { StoreProvider } from '@/contexts/store-context'
import '@/app.css'
import App from '@/App'
import { SCOPE } from '@/lib/scope'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <ThemeProvider>
      <StoreProvider>
        <BrowserRouter basename={SCOPE?.base}>
          <App />
        </BrowserRouter>
      </StoreProvider>
    </ThemeProvider>
  </StrictMode>,
)
