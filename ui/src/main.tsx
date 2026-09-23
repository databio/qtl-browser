import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { BrowserRouter } from 'react-router'
import { ThemeProvider } from '@/contexts/theme-context'
import { StoreProvider } from '@/contexts/store-context'
import '@/app.css'
import App from '@/App'

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <ThemeProvider>
      <StoreProvider>
        <BrowserRouter>
          <App />
        </BrowserRouter>
      </StoreProvider>
    </ThemeProvider>
  </StrictMode>,
)
