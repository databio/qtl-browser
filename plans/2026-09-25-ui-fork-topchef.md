---
date: 2026-09-25
status: complete
model: Claude Opus 5 (1M context)
description: Fork the browser in two: ui-topchef keeps shipping TOPCHeF, ui becomes the multi-study browser
---

# Fork the UI: `ui-topchef` keeps shipping, `ui` goes multi-study

## Where this is going

The destination is one app that treats the experiment as part of the route
(`/e/topchef/gene/FLNC`), with TOPCHeF as a default rather than an assumption, and cross-study
views reachable — the thing the qtlstore was built for and the reason `qtlstore crosscat` exists.
That is a large change and the TOPCHeF site is live at https://topchef.databio.org.

So: fork first. `ui-topchef` keeps the current app deliverable and unchanged; `ui` becomes the one
that evolves. This plan covers only the fork. Generalisation is separate work.

## What is already true

The data layer needs nothing. The store on B2 already holds two experiments built from two
adapters, exercising every axis the format varies:

| | TOPCHeF | GTEx v8 heart LV |
|---|---|---|
| catalog | `topchef_grch38`, 9,182,261 sites | `gtex_v8_heart_lv_grch38`, 9,613,073 sites |
| annotation | `gencode_v34` | `gencode_v39` |
| dof | 435 / 480 published | 367 fitted |
| trans, GWAS | both | neither |

Both catalogs anchor to seqcol collection `EiFob05aCWgVU_B_Ae0cypnQut3cxUP1`, so a site key
`(seq_digest, pos, ref, alt)` means the same thing in both — the precondition for cross-study
lookup. Adding a third study is a `config.yaml` block and an adapter run, not code.

Everything missing is in the browser: one experiment per build (`VITE_EXPERIMENT`), singleton
caches, no TypeScript cross-catalog lookup, `EQTL_TYPE`/`SQTL_TYPE` as constants, and the
TOPCHeF/DCM coupling in `ColocLoci`, `About`, `Gene`, `LocusCompare`, `Home`, the nav, and the
themes.

## The deploy contract is the dangerous part

`ui/wrangler.jsonc` is not just config in the repo. Cloudflare **Workers Builds** is configured
with **root directory `ui`**, build `npm run build`, deploy `npx wrangler deploy`, Worker name
`qtl-browser`, serving https://topchef.databio.org.

If `ui/` becomes the general browser while that setting stands, a build ships the half-finished
multi-study app over the live site, with nobody running a command. This orders the whole plan.

**Unknown, must be checked in the dashboard:** which branch Workers Builds watches. If `main`,
nothing happens until PR #2 merges and the risk is deferred. If `refget`, the risk is live on the
next push. The plan assumes the worse case.

## Steps

### 1. Repoint Cloudflare (user, dashboard) — deferred, gates the push

In the Workers Builds settings for `qtl-browser`, set **root directory** from `ui` to `ui-topchef`,
and note which branch it builds from.

Doing this *before* the move means that in the window between the two, a build fails (the directory
does not exist yet) rather than deploying the wrong app. A failed build leaves the last successful
deploy serving; that is the fail-safe direction.

### 2. Commit A — the move only

- `git mv ui ui-topchef` (87 tracked files; rename detection keeps per-file history).
- `ui-topchef/wrangler.jsonc` keeps `"name": "qtl-browser"` and its comment updated to say root
  directory `ui-topchef`.
- Update references: `README.md` (6), `pipeline/README.md` (1), `SPEC.md` §Status (1),
  `yoke.toml` (`--ignore=ui/dist`, `--ignore=ui/bench/out`).
- **No `ui/` exists after this commit.** That is deliberate: there is nothing for a stale deploy
  config to pick up.

### 3. Commit B — the general app

- Copy `ui-topchef` to `ui` (tracked files only; `node_modules` and `dist` are ignored).
- **Delete `ui/wrangler.jsonc`.** The general app has no deploy target until it earns one, and its
  absence means it cannot be deployed by accident. It gets a new Worker name, never `qtl-browser`,
  when it is ready.
- `ui/.env.production` keeps `VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser` — same store,
  same objects; the two apps differ in what they render, not what they read.
- `VITE_EXPERIMENT` default stays `topchef` in both, so `ui` runs identically on day one and
  changes only as generalisation lands.
- Add `ui/README.md` saying what this directory is for and where it is going.

### 4. Divergence

Both copies stay live: features may land in either. Nothing enforces parity, so a fix that belongs
in both is a decision someone has to make each time, and nothing catches one that was missed. The
files most likely to need porting are the shared readers (`lib/store.ts`, `lib/store-decode.ts`,
`lib/rounding.ts`) and the plot components, which the fork left byte-identical.

## Decisions & ownership

| decision | tag | note |
|---|---|---|
| Repoint Workers Builds root directory to `ui-topchef` | **user-owned, open 2026-09-25** | Sam will raise it with Nathan. The live Worker is in the databio account, so it is probably not Sam's to change; his own Workers project (topchef account, builds from `main`) has the same setting and is his. **Gate: nothing here may be pushed until both are repointed.** |
| Which branch Workers Builds watches | **open / vague** | Decides whether the risk is live now or at merge. Must be read off the dashboard; the plan assumes the worse case. |
| Move and re-create in two commits, not one | AI-owned, defended | Leaves a window where no `ui/` exists, which is the only state a stale deploy config cannot misuse. |
| New `ui/` ships with no `wrangler.jsonc` | AI-owned, defended | Cannot be deployed by accident; a new Worker name is a deliberate act later. |
| Both copies stay live; no freeze on `ui-topchef` | **user-owned, decided 2026-09-25** | Sam chose this over freezing. Consequence, stated plainly: there is no rule to fall back on, so every shared fix is a conscious decision about whether to port it, and nothing flags one that was not. |
| `ui/` starts as a full copy rather than a fresh skeleton | AI-owned, default | It runs against the real store from the first commit, so generalisation can be incremental and always-working. A skeleton would be cleaner and unusable for weeks. |
| Tell Nathan before pushing the move | **user-owned, decided 2026-09-25**: Sam will confirm when to push; the work stays local until then | He pushed `e13a4ed` to `refget` yesterday touching `ui/vite.config.ts`. A directory rename conflicts with any in-flight work of his. |
| `VITE_EXPERIMENT` default stays `topchef` in the general app | AI-owned, default | Worth challenging once an experiment picker exists. |

## What this changes elsewhere

- **The live site.** Nothing else in this repo can take down https://topchef.databio.org; this can.
  The mitigation is the push gate: the rename is local until both Workers projects are repointed.
- **Disk and installs.** A second `node_modules` (`ui/` is 360 MB, 357 MB of it `node_modules`) and
  a second `package-lock.json`. Two `npm install`s, two dependency-bump surfaces.
- **The Rivanna sync.** `yoke.toml` ignores `ui/dist` and `ui/bench/out` by path; stale ignores
  would sync a 100 MB+ `ui-topchef/dist` to the cluster.
- **Bench and check scripts.** `ui/bench/` and `ui/scripts/store-check.ts` are duplicated. Both are
  run by hand, so the cost is remembering which copy was run, not broken automation.
- **Local dev localStorage.** Both apps on `localhost:5173` share an origin, so `topchef-theme` is
  shared between them. The store and chrom-size caches are keyed by data base, experiment and
  digest, so those do not collide; the theme does, harmlessly.
- **Git history.** `git mv` keeps per-file history through rename detection; the `ui/` copy in
  commit B is 87 new files with no ancestry. Blame for the general app starts there.
- **This session's fixes** (`ColocLoci` mount, speculative pointer wave, rounding notes, seqcol
  digest, `store-check` assertion) are already in `93800db`, so both copies inherit them. That is
  why the fork is cheapest now rather than after more UI work lands.
- **PR #2** grows by two commits and a large rename. Reviewable, since commit A is a pure move.

## Not in this plan

Generalisation itself: experiment in the route, per-experiment caches, TypeScript cross-catalog
lookup, reading phenotype types from `results[]`, lifting coloc/DCM/branding into a study config.
Each is its own change against `ui/` once the fork is in.

## Implementation log

**2026-09-25.** Fork done locally, uncommitted, not pushed.

- `git mv ui ui-topchef`: 87 tracked files, all detected as pure renames.
- `cp -Rc ui-topchef ui` (APFS clone, so the second 357 MB `node_modules` costs no disk until the
  copies diverge), then removed `ui/wrangler.jsonc`, `ui/dist` and `ui/bench/out`.
- References updated: `README.md` (layout table now lists both, plus five path fixes),
  `pipeline/README.md`, `SPEC.md` §Status, `yoke.toml` (added `ui-topchef/dist` and
  `ui-topchef/bench/out` to the sync ignores; the `ui/` ones still apply to the new directory).
- `ui-topchef/wrangler.jsonc` comment records the root-directory change the dashboard still needs.
- `ui/README.md` written: what the directory is for and the five steps toward shape 3.
- `npx tsc -b` clean in both.

Both decisions taken since the plan was written are folded into the ledger above. The inversion
option (leave the live app at `ui/`, build the general one at `ui-multi/`) was raised and declined:
Sam prefers the final layout now and will coordinate the dashboard change with Nathan.

**Committed 3d0ae28** (fork) on `refget`, unpushed. `git log --follow` reaches through the rename for
`ui-topchef/` (3 commits against `ui/`'s 2 on a sample file), so per-file history is retrievable for
the deployed app even though both directories now hold the same content.

**The one thing still open is the push gate**: Workers Builds root directory `ui` -> `ui-topchef` in
the databio account, and the same setting in Sam's own topchef account, which builds from `main`.
Nothing in this branch may be pushed until both are moved.

## Inverted, 2026-10-02

Sam flipped the naming: **`ui/` stays the deployed TOPCHeF app and the general fork is
`ui-general/`**, to be renamed to `ui/` along with the deployment "in a week or so" once it has made
progress. This is the inversion raised and declined on 2026-09-25, taken now that the coloc work has
gone into the deployed app and the dashboard change has not happened.

**It removes the push gate.** Workers Builds keeps root directory `ui`, which still holds the
deployed app, so nothing in either Cloudflare account has to change before this branch is pushed.
That was the only thing blocking it.

**It moves the cost to the flip.** `ui-general/` becomes the deployed app at the rename, so it has
to carry everything `ui/` does by then or the flip is a regression.

The one divergence that existed -- the coloc.abf work (`lib/coloc-abf.ts`, `gwasCols` on
`GenePack`, the computed `ColocSection`, the `KvTable` `labelWidth`) which landed in `ui/` after the
fork -- was ported the same day, while it was four files and `diff -rq ui/src ui-general/src` was
empty afterwards. From here the two diverge again with every change to either, and nothing checks:
`diff -rq ui/src ui-general/src` is the whole audit, and it is worth running before the flip.

References updated back: `README.md`, `pipeline/README.md`, `SPEC.md` §Status, `yoke.toml`,
`.gitignore`, `ui/wrangler.jsonc` (its comment now records that the dashboard setting moves only at
the flip), and `ui-general/README.md`.

## Flipped, 2026-10-05 (PR #4, merged `a3b7bea`)

The layout the 2026-09-25 plan aimed at, reached three days after the inversion rather than the
week estimated: **`ui/` is the general browser, `ui-topchef/` is the deployed TOPCHeF app.** Git
recorded it as a single rename of `ui-general/` to `ui-topchef/` with `ui/` left untouched, because
the two trees were byte-identical, so a swap and a one-way rename produce the same end state. The
cost this section predicted -- the general app having to carry everything the deployed one did by
the flip -- came to nothing, for the same reason.

**Deployment left Workers Builds for GitHub Actions**, which is what actually unblocked the flip.
Rather than reconnecting the Cloudflare GitHub App installation that the `sanghoonio` ->
`databio` repo transfer orphaned, `.github/workflows/deploy-ui.yml` and `deploy-topchef.yml` build
and deploy each app with `cloudflare/wrangler-action`. Each is filtered to its own directory and its
own workflow file and sits in its own concurrency group, so a change to one app never redeploys the
other. Needs `CLOUDFLARE_API_TOKEN` (Workers Scripts edit, plus Workers Routes and DNS edit on the
`databio.org` zone for the custom domain) and `CLOUDFLARE_ACCOUNT_ID`.

Worker names no longer follow the directory that holds them: `ui/` deploys to `qtl-browser`
(workers.dev only) and `ui-topchef/` to `qtl-browser-topchef`, which declares
`topchef.databio.org` as a custom domain so its first deploy moves the domain off `qtl-browser` --
the Worker that served TOPCHeF before the rename and now holds the general app. Both workflows ran
green on the merge and the site answers 200.

**The one thing not verified** is that the domain actually moved. `diff -rq ui/src ui-topchef/src`
is still empty, so both Workers serve the same bundle and no external request distinguishes them; a
green `wrangler deploy` only means the command exited 0. If the domain did not move, every push
touching `ui/` deploys the general app to the Worker holding `topchef.databio.org`, and that stays
invisible until the two apps diverge. `wrangler deployments list --name qtl-browser-topchef`, or the
Worker's Domains & Routes panel, settles it.
