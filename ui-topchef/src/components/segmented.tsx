/** Compact segmented control: bordered track, active segment filled. `nav` renders a <nav>. */
export function Segmented<T extends string>({ options, value, onChange, nav = false }: {
  options: { value: T; label: string; disabled?: boolean; title?: string }[]
  value: T; onChange: (v: T) => void; nav?: boolean
}) {
  const Wrap = nav ? 'nav' : 'div'
  return (
    <Wrap className="inline-flex h-8 items-center gap-0.5 rounded-[0.625rem] bg-base-200 p-0.5">
      {options.map(o => (
        <button key={o.value} type="button" disabled={o.disabled} title={o.title}
          onClick={() => !o.disabled && onChange(o.value)}
          className={`h-full rounded-lg border px-3 text-sm font-medium transition-colors ${
            o.disabled ? 'cursor-not-allowed border-transparent text-base-content/25'
              : value === o.value ? 'border-base-300 bg-base-100 text-base-content'
                : 'border-transparent text-base-content/60 hover:text-base-content'}`}>
          {o.label}
        </button>
      ))}
    </Wrap>
  )
}
