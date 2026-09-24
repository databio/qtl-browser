---
date: 2026-09-23
status: complete
model: Claude Fable 5.1
description: Five fixes to PR #1 (feat/pack-format) before it merges into main, including repointing deploy tooling from R2 to the Backblaze B2 host that is actually serving the packs. Superseded on 2026-09-24 by the refget branch, which contains PR #1, deletes the v0 upload tooling, and is what is deployed.
---

# PR #1 pre-merge fixes

PR #1 (`feat/pack-format`, head `9cde62f`) replaces the parquet reader with the qtlb v0 pack
format. Its browser code is live at `https://qtl-browser.nsheff.workers.dev`, reading packs from
`https://cloud2.databio.org/qtl-browser`, a Backblaze B2 bucket behind Cloudflare. `main` still
describes the parquet deployment on R2. Merging PR #1 closes that gap. These five items are what
should land on the branch first. The `refget` branch is a separate, later PR and is not touched
here.

Facts established while reviewing, which the plan relies on:

- A fresh build of PR #1 fails at `pack_eqtl`: `steps_nominal.py:167` writes to
  `data/derived/cis_eqtl_nominal/` while `steps_pack.py`, `steps_finish.py`, and `packcheck.py`
  read `data/derived/_tables/cis_eqtl_nominal/`. It never showed because the tables were already in
  place on Nathan's machine.
- PR #1 deletes `data/raw/download.py` and `data/raw/sources.yaml`; the README still documents both.
- `packcheck report` reads `manifest["tables"]`, a block PR #1 removed from the manifest.
- Nothing in any branch mentions Backblaze or `cloud2.databio.org`. `ui/.env.production`,
  `config.yaml`'s `r2:` block, `upload.py`, and the README Deploy section all target R2.
- The B2 host answers Range with 206, sets `Cache-Control: public, max-age=31536000, immutable` on
  packs and `no-cache` on the manifest, and handles CORS at Cloudflare (allow-origin echoes the
  nsheff origin, allow-headers `range`). It sends **no ETag and no Last-Modified** on any object.
  Chromium stores a 206 only when the response carries a strong validator, so pack ranges are
  re-downloaded on every visit. Objects were uploaded as B2 large files (`x-bz-content-sha1: none`),
  so their S3 ETags are not MD5s.
- The About page says gene-level results are exact. On the variant page the "Lead variant for"
  table reads perm p, slope, and SE from the hits pack, where they are quantized against frame
  maxima (`variant-pack.ts` `leadValues`).
- The locus tooltip label is built by SQL concatenation in `pack.ts:229`. A p = 0 row decodes to a
  NaN slope, the concatenation goes NULL, and `LocusPlot.tsx:195` renders `String(null)`. No p = 0
  rows exist in the current data (0 of 123,458,578 eQTL rows and 0 of 84,498,302 sQTL rows locally),
  so this is dormant.

## Branch mechanics

Work on a local branch `fix/pr1-premerge` cut from `origin/feat/pack-format`. All five items land
in the working tree first and are shown as one diff; one commit covers them all, and only after
Sam has seen the diff and said so. Nothing is pushed or merged without a separate go-ahead; when
it is, the commit goes onto `feat/pack-format` so PR #1 carries it, and the PR description gets a
short "since review" list.

Order: items 1, 2, 3 are mechanical and go first. Item 5 next. Item 4 is now one env line and one
README paragraph and goes last; the upload-tooling work it originally held is a follow-up.

## Item 1: nominal writes where the pack builders read

`pipeline/steps_nominal.py:167`: `cfg.derived / out_dir` becomes `cfg.tables / out_dir`. One line.
The README's step table already says `_tables/cis_eqtl_nominal`, so no doc change.

Do not take the refget branch's version of this file: it also imports `adapters.topchef`, which
does not exist on PR #1.

Verify: `git grep -n "cis_eqtl_nominal\|cis_sqtl_nominal" pipeline/` shows every reader and the
writer on `cfg.tables`. A local smoke build is not possible here (no raw archives on this machine);
say so in the commit message.

## Item 2: restore the download tooling

```
git checkout origin/refget -- data/raw/download.py data/raw/sources.yaml
```

`sources.yaml` on refget is byte-identical to main. `download.py` on refget is main's plus the
resume guard added after a Zenodo file grew to 13.4 GB on a `-C -` resume against a server that
ignored Range. The script has no imports from `pipeline/`, so it lifts cleanly.

Verify: `python data/raw/download.py --list` runs and prints a status per source without touching
the network beyond what `--list` already does (it stats local files only).

## Item 3: packcheck reads the moved tables and the new manifest

Apply the refget branch's `packcheck.py` diff as-is; it contains only these hunks and nothing that
depends on refget code:

```
git diff origin/feat/pack-format origin/refget -- pipeline/packcheck.py | git apply
```

Three changes: `resolve_source` and `cmd_report` glob under `cfg.tables` instead of `cfg.derived`;
`manifest_rows` comes from `manifest["packs"]["counts"][name]["rows"]` instead of the deleted
`tables` block; the report's markdown names the new source.

Verify: `git apply --check` first; then `uv run python -c "import pipeline.packcheck"`.

## Item 4: point the browser and the upload tooling at the bucket in use

First deferred on 2026-09-24 for want of credentials, then done the same day once the B2
application key landed in `databio/secrets` and the bucket facts were confirmed. 4a and 4b went in
the first commit; the profile work below went in after it.

### 4a. The browser's data host

`ui/.env.production`: `VITE_DATA_BASE=https://cloud2.databio.org/qtl-browser`. The comment above it
stays. This is the only change the browser needs; `manifest.ts` already builds every URL from it.

### 4b. One README paragraph

In the Deploy section, before the command block: the live data host is `cloud2.databio.org`, a
Backblaze B2 bucket Nathan owns and uploads to by his own process; `upload.py` and the `r2:` block
still target the R2 bucket, which serves the old parquet site and is slated to be turned off; B2
support in `upload.py` is a follow-up. Also correct the sentence "He owns the bucket, the R2 token,
and the Workers project" to say which bucket it means. No other tooling change.

Verify item 4: `cd ui && npx tsc -b && npm run build`; the built bundle contains the cloud2 host
and no `r2.dev` string (`grep -c r2.dev dist/assets/*.js` is 0).

## Item 4c: B2 profile for upload.py

Implemented 2026-09-24. The bucket facts were confirmed with the application key Nathan added to
the databio secrets store:

| | |
|---|---|
| credentials | `pass show cloud-file-service/b2_access_key_id` and `cloud-file-service/b2_secret_access_key` (the `~/.password-store` clone of `databio/secrets`) |
| S3 endpoint | `https://s3.us-west-002.backblazeb2.com` |
| bucket | `cloud-databio` |
| key prefix | `qtl-browser/` (objects are `qtl-browser/manifest.json`, `qtl-browser/immutable/...`) |
| key scope | one bucket; `ListBuckets` is refused, which is expected for a scoped key |

The prefix was the one thing the tool could not express: it wrote keys at the bucket root. The
profile's `key_prefix` is now applied inside the four bucket primitives only (`_list` strips it and
drops keys outside it, `_put` and `_delete` add it, `_public` needs nothing because the public URL
already carries it as a path), so every other function, every listing file, and the manifest keep
the unprefixed key vocabulary. `""` for the R2 profile.

### The bucket block in `config.yaml`: two profiles, one tool

R2 stays supported until Sam turns it off (planned, along with possibly the topchef Cloudflare
account). B2 is added beside it. The `r2:` block becomes:

```yaml
storage:
  immutable_prefix: immutable/          # provider-independent; common.py reads it
  cache_control: {immutable: "...", mutable: "no-cache"}
  single_put_max_bytes: 4000000000
  multipart_threshold_bytes: 8388608
  default_profile: b2
  profiles:
    r2:
      provider: r2
      bucket: qtl-browser
      key_prefix: ""
      endpoint: https://02604ac5….r2.cloudflarestorage.com
      public_url: https://pub-50e5b2fc….r2.dev
      allowed_origins: [https://qtl-browser.topchef.workers.dev, http://localhost:5173, http://localhost:4173]
      free_tier: {storage_bytes: 10000000000, class_a_ops: 1000000, class_b_ops: 10000000}
    b2:
      provider: b2
      bucket: cloud-databio
      key_prefix: qtl-browser/
      endpoint: https://s3.us-west-002.backblazeb2.com
      public_url: https://cloud2.databio.org/qtl-browser
      allowed_origins: [https://qtl-browser.nsheff.workers.dev, http://localhost:5173, http://localhost:4173]
      free_tier: {storage_bytes: <from B2 pricing page>}   # no ops figures; see below
```

`upload.py` takes `--profile r2|b2`, defaulting to `storage.default_profile`. Credentials come from
`<PROFILE>_ACCESS_KEY_ID` / `<PROFILE>_SECRET_ACCESS_KEY`, so the existing `R2_*` pair keeps
working unchanged and the B2 pair is `B2_ACCESS_KEY_ID` / `B2_SECRET_ACCESS_KEY`, which `.env` can
fill from `pass` without copying the values anywhere else:

```
B2_ACCESS_KEY_ID=$(pass show cloud-file-service/b2_access_key_id)
B2_SECRET_ACCESS_KEY=$(pass show cloud-file-service/b2_secret_access_key)
```

`common.py:41` reads `storage.immutable_prefix`. `test_upload.py` runs its fake-bucket suite once
per profile, including a prefixed one.

### ETag semantics in `stage` and `check`

Two places assume a bucket ETag is the object's MD5, which held on R2 for single-part PUTs and does
not hold for the objects already in the B2 bucket.

- `stage` (upload.py ~469) refuses any present immutable key whose ETag differs from the manifest
  MD5. Against the B2 bucket that refuses every one of the 147 present keys and exits 1. Change:
  when the bucket ETag is not a 32-character hex string (a multipart hash), compare size only, the
  rule `_plain_changed` already uses above the threshold. When it is an MD5, keep the strict check.
- `check` fails every key whose public-URL ETag differs from the MD5. The B2 host sends no ETag.
  Change: skip the ETag comparison when the header is absent and print one summary line saying
  validators are missing; fail only when an ETag is present and differs. Add a per-run check line:
  "ETag or Last-Modified present on immutable keys", reported as a warning, not a failure, because
  the fix is on the host side (see 4e).

`test_upload.py` gets one case for each: a present key with a multipart ETag is "present", and a
HEAD without ETag passes with the warning.

### `budget`, `cors`, `prune`

- `budget`: the class A / class B monthly logic is R2-shaped and stays for the `r2` profile. For
  `provider: b2` report storage against the B2 free tier and omit the ops lines. Do not invent B2
  transaction-class numbers.
- `cors`: CORS for `cloud2.databio.org` is set at Cloudflare, not on the bucket. Leave the command
  in place, but the README says it applies to a bucket served directly and is not what the current
  host uses. Whether the bucket itself should also carry a CORS rule is Nathan's call.
- `prune` and `check --retired`: the B2 bucket never held parquet, so the `replaces:` prefixes are
  already empty there. Both commands work unchanged and become no-ops; note that in the README.

### The missing validators (not a repo change, but raised in the PR)

`cloud2.databio.org` strips ETag and Last-Modified. B2 sends both natively, so whatever sits in
front (a Worker or a Cloudflare transform rule) is dropping them. Until it passes them through,
every pack range is fetched again on every page load, which undoes the caching the PR's bench
numbers assume. This is Nathan's infrastructure. `upload.py check` now counts the keys without a
validator and prints one warning; against the live host on 2026-09-24 that was 140 of 142.

### What changed, as built

- `config.yaml`: `r2:` became `storage:` with shared settings (`immutable_prefix`, `cache_control`,
  the two size limits) and `profiles: {b2, r2}`, each with `provider`, `bucket`, `key_prefix`,
  `endpoint`, `public_url`, `cors` (`edge` or `bucket`), `allowed_origins`, `free_tier`.
  `default_profile: b2`. The B2 free tier is storage only, read from B2's pricing page that day:
  first 10 GB free, Class A/B/C calls free, egress free through Cloudflare.
- `common.py` reads `storage.immutable_prefix`.
- `upload.py`: `prof(cfg)` merges the active profile over the shared settings; `--profile` sets it;
  `_env` reads `<PROFILE>_ACCESS_KEY_ID` / `_SECRET_ACCESS_KEY` and signs for the region named in a
  B2 endpoint (`auto` on R2); `_is_md5` decides whether an ETag can be compared; `stage` treats a
  present key with a multipart ETag as present on size alone; `check` compares an ETag only when
  one is sent and MD5-shaped, warns once about missing validators and once per key about missing
  expose headers (neither fails the run, since the reader uses neither); `budget` prints the ops
  sections only for a tier that prices ops; `cors` prints and refuses to apply on a `cors: edge`
  profile; `_finish` names the Cloudflare zone instead of "CORS is Sam's" on that profile;
  `prune --stale` matches only v0-shaped names under `immutable/`, so the v1 store objects sharing
  the prefix are listed under "left alone" and never deleted.
- `test_upload.py` runs every case once per profile and adds four: multipart ETag accepted by
  `stage` (and still refused on a size mismatch), prefix and profile helpers, the edge-CORS refusal,
  and three `check` variants (missing expose header warns, missing ETag warns, multipart ETag is
  not compared). The `etag` failure case uses an MD5-shaped wrong value, since a non-MD5 shape is
  no longer a failure.
- README Deploy section and `pipeline/README.md` describe the profiles, the `pass` entries, the
  size-only rule, and that `prune` and `cors` are no-ops or print-only on B2.

Verified: 17 cases x 2 profiles pass; `check` through the `b2` profile against the live host,
using the live manifest as the local one, passes with 3 warnings; `inventory --public` lists all
146 named keys. Not run: `inventory`, `budget`, and `stage --dryrun` with the real key, which need
`B2_*` in `.env` and are Sam's to run.

## Item 5: two small correctness fixes in the UI

- `ui/src/routes/About.tsx`: the "Rounded values" bullet cut to two sentences. Per-variant
  p-values, slopes, SEs, and allele frequencies are "stored in compressed form" (no format
  vocabulary), with only the p-value percentage and the slope bound in SE units quoted; gene-level
  results and the GWAS values are exact; exact per-variant values are on Zenodo. The variant page's
  lead table is not called out: its slope and SE are the lead variant's own per-variant values,
  already covered by the first sentence. Doc fix, not a format change: storing exact lead values
  would mean a new hits layout, a rebuild, and a re-upload, and that is not pre-merge work.
- `ui/src/lib/pack.ts:229`: wrap the slope fragment so a NULL slope yields text instead of
  NULL-ing the whole label: `CASE WHEN q.slope IS NULL THEN 'slope n/a (p underflow)' ELSE 'slope '
  || format(...) || ' ± ' || format(...) END`.
- `ui/src/components/LocusPlot.tsx:195`: `r.label == null ? '' : String(r.label)` so a NULL label
  can never render as the word "null" even if some other fragment goes NULL later.

Verify: `cd ui && npx tsc -b` (not `tsc -p .`, which is a no-op here) and `npm run build`.

## Out of scope, noted for later

- `locusTable` throws for a gene with no eQTL rows when opened with `?tab=eqtl`; main showed an
  empty plot. Small, but a behavior change worth its own commit after merge.
- `SPEC.md` and `PACKS.md` live in the analysis repo by design (refget branch README). A link from
  this README is enough and can ride along with the refget PR.
- Decommissioning the R2 bucket and the topchef Workers deployment.
- Everything on the `refget` branch.

## Decisions & ownership

| decision | tag | note |
|---|---|---|
| Merge PR #1 before refget, as its own PR | user-owned, surfaced | agreed in conversation 2026-09-23 |
| One commit for all five items, made only after Sam has seen the diff | user-owned, decided 2026-09-24 | Sam's standing rule for every project |
| The commit goes onto `feat/pack-format` rather than a PR into that branch | AI-owned, default | simplest for a single-commit PR; challenge if Nathan wants to review it separately |
| `cloud2.databio.org/qtl-browser` is the long-term data host | user-owned, surfaced | Nathan's host; the plan assumes it stays |
| B2 support in `upload.py` lands in this PR after all, once the key arrived | user-owned, decided 2026-09-24 | the `r2` profile stays supported alongside it |
| R2 bucket and the topchef Cloudflare account stay; the old URL keeps working | user-owned, decided 2026-09-24 (revised from "turn R2 off") | free under both tiers. The topchef Workers deployment rebuilds from `main`, so after merge it serves the pack reader against the same B2 data |
| R2 is frozen as a backup: no writes, no prune, ever | user-owned, decided 2026-09-24 | the `r2` profile stays in the tool for read-only commands; the README says so. Local copies of the parquet build under `data/derived/` can be cleaned separately, outside this PR |
| Two named profiles under `storage:` with `--profile` and per-profile env vars | AI-owned, defended | keeps R2 working untouched while B2 is added; a single block with swapped values would have broken R2 the moment B2 went in |
| `stage` size-only match on multipart ETags; `check` warns rather than fails on a missing ETag or expose header | AI-owned, defended | the alternative refused every present key and failed every check against the live host; the reader uses neither header, and a same-size collision on a content-addressed name is theoretical |
| `budget` storage-only when the tier has no ops figures | AI-owned, defended | driven by the config, not the provider string; B2's page says Class A/B/C calls are free |
| `cors: edge` profiles print and refuse to apply | AI-owned, defended | a bucket rule would never be consulted behind the Cloudflare zone; applying it would look like it did something |
| Bucket `cloud-databio`, endpoint `s3.us-west-002`, prefix `qtl-browser/`, credentials in `pass` under `cloud-file-service/` | user-owned, surfaced (settled 2026-09-24) | confirmed against the live bucket with the scoped key |
| Free-tier figure, which origins to allow | open / vague | the figure is read from B2's pricing page at follow-up time, not typed from memory; origins are Nathan's call |
| CORS stays at Cloudflare; bucket-level rule is optional | user-owned, surfaced | Nathan's infrastructure. Probed 2026-09-24: the zone reflects any origin, so `allowed_origins` on the `b2` profile is documentation and the Origin `check` sends, not an enforced list |
| `https://topchef.databio.org` is the site's name; the workers.dev host is the same deployment | user-owned, surfaced (Sam, 2026-09-24) | first entry in the `b2` profile's origins; README names it |
| ETag / Last-Modified pass-through on `cloud2.databio.org` | user-owned, surfaced | Nathan's infrastructure; goes in this PR's description because it affects visitors now |
| About-page sentence corrected rather than hits pack made exact | AI-owned, defended | format change is a rebuild and re-upload; not pre-merge |
| About-page rounding note is two plain sentences with two figures, no format terms, no variant-page exception | user-owned, decided 2026-09-24 | Sam: the old text hedged with numbers, "16-bit codes" means nothing to a user, and the variant page is per-variant by nature |
| p = 0 tooltip guard is SQL `CASE` plus a null check, no format change | AI-owned, defended | zero affected rows today; the guard costs nothing |
| No local smoke build for item 1 | silently implied, now named | this machine has no raw archives; the fix is verified by reading, and by Nathan's next build on Rivanna |

## What this changes elsewhere

- **Item 1** changes where a fresh `nominal` run writes. On a machine that already has the tables
  under `_tables/` (Nathan's), the `.done` marker means nothing re-runs and nothing moves. On a
  machine with them under `data/derived/` (Sam's, from the main-branch build), the old copies stay
  where they are as dead weight until removed by hand; nothing reads them.
- **Item 2** re-adds two files PR #1 deleted; the refget branch, which also adds them, will merge
  cleanly on top only if it is rebased after this lands, since both branches will carry the same
  content under the same paths.
- **Item 3** is confined to `packcheck`; no build step reads its output except `manifest`, which
  reads `roundtrip_genome.json` and is unaffected.
- **Item 4a** changes which host every production visitor fetches from. The topchef Workers site
  keeps its old parquet-reading bundle until it is rebuilt from `main`; once `main` carries this
  change, a Workers Builds deploy of that project would also start reading from `cloud2`. Two
  frontends reading one bucket is fine; the point is that the old R2 bucket stops being reachable
  from any current build.
- **Item 4c** makes the repo's deploy commands write to the live bucket by default. Every write
  still needs the `B2_*` key in `.env` and an explicit go-ahead; `stage` cannot overwrite an
  immutable key, and `release` cannot run before `stage`. It also relaxes one guarantee: a key whose
  bucket ETag is a multipart hash is trusted on size alone. Content-addressed names make a same-size
  collision require a SHA-256 prefix collision, so the exposure is theoretical.
- **The bucket is shared.** `cloud-databio` holds other lab projects. `_list` returns only keys
  under `qtl-browser/`, so `budget`, `prune`, and `check --retired` can neither count nor delete
  anything outside the prefix. That rule lives in one function and is the thing to protect.
- **The prefix is shared too.** Sam's first `inventory` with the real key (2026-09-24) listed 373
  keys and 7.40 GB under `qtl-browser/`, against 146 keys and 3.18 GB the live manifest names. The
  other 227 keys (4.22 GB) are a qtlb v1 store Nathan uploaded: `store.json` ("store-genome-v1e",
  format version 1), experiments `topchef` and `gtex_v8_heart_lv` (dof 367, fitted), catalogs
  `topchef_grch38` (9,182,261 sites) and `gtex_v8_heart_lv_grch38`, annotations GENCODE v34 and
  v39, and `immutable/<sha512t24u>.<ext>` objects. Its pointer layout (`variant_catalogs/`, an
  experiment JSON with `hits`, `trans`, `gwas` blocks) is ahead of the pushed `refget` branch.
  `prune --stale` would have listed every one of those objects as a stale build; it now matches
  only v0-shaped names, with a test. Storage stands at 74% of the free tier, so the next v0 build
  (about 3.2 GB) does not fit until the superseded v0 keys are pruned after its release.
- **The validator gap** on `cloud2.databio.org` is a caching regression that exists today
  regardless of this plan; `check` now reports it, and naming it in the PR makes it Nathan's to fix
  on the host.
- **Item 5** changes visible text on the About page and the tooltip for a case that does not occur
  in the current data. No data rebuild.
- **Credits and cost**: no uploads, no deletes, and no bucket writes are part of this plan. The
  only network calls are read-only `check` HEADs against the public host.

## Implementation log

2026-09-24. All five items applied in the working tree on local branch `fix/pr1-premerge` (from
`origin/feat/pack-format` at `9cde62f`), uncommitted, awaiting Sam's review of the diff.

- Item 1: `pipeline/steps_nominal.py` writes under `cfg.tables`.
- Item 2: `data/raw/download.py` and `sources.yaml` restored from the refget branch (`sources.yaml`
  identical to main; `download.py` is main's plus the resume guard).
- Item 3: the refget branch's `packcheck.py` diff applied verbatim with `git apply`.
- Item 4: `ui/.env.production` points at `cloud2.databio.org/qtl-browser`; README Deploy section
  opens with the paragraph on where the data lives and the validator gap.
- Item 5: About.tsx sentence, `pack.ts` label `CASE`, `LocusPlot.tsx` null guard.

Checks: `npx tsc -b` clean; `npm run build` clean, bundle contains the cloud2 host and no
`r2.dev`; `uv sync` then `import pipeline.packcheck, steps_nominal, steps_pack` ok;
`test_packfmt --synthetic` 28 passed; `download.py --list` runs; the SQL comment inside the label
concatenation parses in DuckDB. `npm install` touched `ui/package-lock.json` (peer flags only); that
change was reverted so the diff holds only the five items. No local smoke build (no raw archives
here).

2026-09-24, later. Item 4c (the B2 profile) implemented on the same branch after the first commit
(`70524bd`) and the link commit (`c66516e`): `config.yaml`, `common.py`, `upload.py`,
`test_upload.py`, `README.md`, `pipeline/README.md`, and this plan. Uncommitted, awaiting Sam's
review of the diff. Verification as listed under "What changed, as built".

2026-09-24, end of day. Superseded. Nathan force-pushed `refget` (`b4812ec`): it merges
`feat/pack-format` at `70524bd`, adds a v1 browser reader (`ui/src/lib/store.ts` and friends,
replacing `pack.ts`, `pack-decode.ts`, `manifest.ts`), and removes the v0 build and deploy code:
`upload.py`, `test_upload.py`, `packcheck.py`, `steps_pack*.py`, `steps_finish.py`, `steps_gwas.py`,
the `r2:` config block. The bucket was switched the same day: every v0 pack, `manifest.json`, and
the three plain JSON files were deleted from `qtl-browser/`, and `topchef.databio.org` now serves
the v1 bundle. PR #1 is an ancestor of `refget`, so merging it alone would put `main` on a reader
with no data; the route to `main` is a PR from `refget`.

What landed on PR #1 (`feat/pack-format`): `70524bd` (items 1, 2, 3, 4a, 4b, 5) and, in the
commit after this entry, the `.DS_Store` ignore line, the About-page repo link and two-sentence
rounding note, and this plan. The B2 profile for `upload.py` (item 4c) is shelved on the local
branch `wip/b2-upload-profile` (`9b8b410`), never pushed: the tool it extends no longer exists on
`refget`, and a v1 uploader would start from that code but target the v1 layout (objects first,
`store.json` last). The README and `pipeline/README.md` rewrites from item 4c are shelved with it;
`refget`'s own Deploy section already names the B2 bucket, prefix, and host and points at Nathan's
cloud-management notes for the steps.
