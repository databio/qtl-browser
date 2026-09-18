"""Stage, release, check and prune the R2 bucket the browser reads (`r2:` in pipeline/config.yaml).

    # <repo>/.env (gitignored):  R2_ACCESS_KEY_ID=...  R2_SECRET_ACCESS_KEY=...
    # token: Object Read & Write on this bucket. Exported variables override the file.

    uv run python -m pipeline.upload inventory [--public] -o FILE          # read-only bucket listing, saved for --listing
    uv run python -m pipeline.upload budget [--listing FILE] [--harness JSON] [--overlap-days 3]
    uv run python -m pipeline.upload stage [--dryrun] [--listing FILE]     # new immutable/ files + a staged manifest copy; nothing live changes
    uv run python -m pipeline.upload check --staged                        # every staged file on the public URL
    uv run python -m pipeline.upload release [--dryrun] [--listing FILE]   # changed plain-path files, then manifest.json: the live switch
    uv run python -m pipeline.upload check                                 # live manifest equals local, plus the --staged checks
    uv run python -m pipeline.upload pointer FILE                          # rollback: upload FILE as manifest.json
    uv run python -m pipeline.upload prune [--stale] [--yes] [--listing FILE]  # dry run unless --yes
    uv run python -m pipeline.upload check --retired [--listing FILE]      # replaced prefixes are empty
    uv run python -m pipeline.upload cors [--print]                        # unchanged

Every command reads data/derived/manifest.json: the bucket is a copy of what the manifest names,
never a directory walk. Uses the aws CLI against the R2 S3 endpoint; the checksum variables stop
aws-cli 2.23+ from sending CRC headers that R2 rejects.

A deploy is staged, rehearsed, then switched:

    inventory -> budget -> stage --dryrun -> stage -> check --staged -> release --dryrun
             -> release -> check -> prune -> prune --yes -> check --retired

`--listing FILE` reads a saved `inventory` instead of listing the bucket again. It is refused
together with anything that writes, so a stale listing can never drive a real upload or delete.
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .common import CHROMS, ROOT, Config, digests, log

EXPOSE = ["Content-Range", "Content-Length", "Accept-Ranges", "ETag"]
# r2.dev answers 403 (error 1010) to the default Python-urllib user agent
UA = "Mozilla/5.0 (qtl-browser upload.py)"
# small JSON files the manifest names only by its own structure, not in `assets`
PLAIN_EXTRA = ("gwas_dcm.json",)
# where a bucket key's local twin may live. `data/derived/` is the live build; the `_*` folders hold
# copies of retired bucket paths, which is what `prune` means by "there is still a local copy".
LOCAL_ROOTS = ("", "_old", "_old/plan9d", "_tables", "_retired", "_deploy")
MANIFEST = "manifest.json"


# ---- environment and the four bucket primitives -------------------------------------------------

def _dotenv() -> dict[str, str]:
    """KEY=VALUE lines from <repo>/.env (gitignored); the shell environment wins over it."""
    path = ROOT / ".env"
    out = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                out[k.strip().removeprefix("export ").strip()] = v.strip().strip("'\"")
    return out


def _env(cfg: Config) -> dict[str, str]:
    merged = {**_dotenv(), **os.environ}
    key, secret = merged.get("R2_ACCESS_KEY_ID"), merged.get("R2_SECRET_ACCESS_KEY")
    if not (key and secret):
        sys.exit("put R2_ACCESS_KEY_ID and R2_SECRET_ACCESS_KEY in .env at the repo root or export them "
                 "(Cloudflare dashboard > R2 > Manage API tokens)")
    return {
        **os.environ,
        "AWS_ACCESS_KEY_ID": key, "AWS_SECRET_ACCESS_KEY": secret, "AWS_DEFAULT_REGION": "auto",
        "AWS_REQUEST_CHECKSUM_CALCULATION": "when_required", "AWS_RESPONSE_CHECKSUM_VALIDATION": "when_required",
    }


def _aws(cfg: Config, *args: str, **kw) -> subprocess.CompletedProcess:
    env = _env(cfg)                      # resolve credentials before logging, so a missing token
    cmd = ["aws", "--endpoint-url", cfg["r2"]["endpoint"], *args]   # never looks like an attempted call
    log(" ".join(cmd))
    return subprocess.run(cmd, env=env, **kw)


def _list(cfg: Config, listing: Path | None = None) -> dict[str, dict]:
    """Bucket keys -> size, etag. From a saved listing, or one ListObjectsV2 call (Class A per
    1,000 keys; aws-cli pages)."""
    if listing:
        data = json.loads(Path(listing).read_text())
    else:
        r = _aws(cfg, "s3api", "list-objects-v2", "--bucket", cfg["r2"]["bucket"], "--output", "json",
                 capture_output=True, text=True, check=True)
        data = json.loads(r.stdout) if r.stdout.strip() else {}
    return {o["Key"]: {"size": o["Size"], "etag": o["ETag"].strip('"')} for o in data.get("Contents", [])}


def _put(cfg: Config, key: str, path: Path, md5_hex: str, cache_control: str, content_type: str) -> None:
    """Single-part PUT with Content-MD5: R2 rejects corrupted bytes, and the ETag is the MD5.
    Never multipart: each part is a Class A op and the ETag stops being the MD5."""
    if path.stat().st_size > cfg["r2"]["single_put_max_bytes"]:
        sys.exit(f"{key}: {path.stat().st_size} bytes is over single_put_max_bytes")
    md5_b64 = base64.b64encode(bytes.fromhex(md5_hex)).decode()
    r = _aws(cfg, "s3api", "put-object", "--bucket", cfg["r2"]["bucket"], "--key", key, "--body", str(path),
             "--content-type", content_type, "--cache-control", cache_control, "--content-md5", md5_b64,
             capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"{key}: {r.stderr.strip()}")
    if json.loads(r.stdout)["ETag"].strip('"') != md5_hex:
        sys.exit(f"{key}: bucket ETag does not match local MD5")


def _delete(cfg: Config, keys: list[str]) -> None:
    """s3api delete-objects, 1,000 keys per call, Quiet; exit on any Errors entry."""
    for i in range(0, len(keys), 1000):
        batch = keys[i:i + 1000]
        payload = {"Objects": [{"Key": k} for k in batch], "Quiet": True}
        r = _aws(cfg, "s3api", "delete-objects", "--bucket", cfg["r2"]["bucket"], "--delete", json.dumps(payload),
                 capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"delete-objects: {r.stderr.strip()}")
        errors = (json.loads(r.stdout) if r.stdout.strip() else {}).get("Errors") or []
        if errors:
            sys.exit(f"delete-objects: {errors}")
        log(f"deleted {len(batch):,} keys")


def _public(cfg: Config, method: str, path: str, headers: dict[str, str] | None = None) -> tuple[int, dict[str, str], bytes]:
    """HEAD or GET on the public r2.dev URL with a browser User-Agent."""
    url = cfg["r2"]["public_url"].rstrip("/") + "/" + path.lstrip("/")
    req = urllib.request.Request(url, method=method, headers={"User-Agent": UA, **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.status, {k.lower(): v for k, v in resp.headers.items()}, resp.read()
    except urllib.error.HTTPError as e:
        return e.code, {k.lower(): v for k, v in (e.headers or {}).items()}, e.read()
    except urllib.error.URLError as e:                      # DNS, TLS, timeout
        return 0, {}, str(e).encode()


# ---- manifest, keys, prefixes --------------------------------------------------------------------

def _manifest(cfg: Config) -> dict:
    path = cfg.derived / MANIFEST
    if not path.exists():
        sys.exit(f"{path} is missing; run `python -m pipeline build --step manifest`")
    return json.loads(path.read_text())


def _content_type(key: str) -> str:
    """What the bucket serves a key as. Packs, the search index and the GWAS index are opaque bytes
    the browser range-reads; never set Content-Encoding, because offsets address the stored bytes."""
    if key.endswith(".parquet"):
        return "application/vnd.apache.parquet"
    if key.endswith(".json"):
        return "application/json"
    return "application/octet-stream"


def _gb(n: int) -> str:
    return f"{n:>16,} B  {n / 2**30:>8.3f} GiB  {n / 1e9:>8.3f} GB"


def _plain_paths(cfg: Config, man: dict) -> dict[str, Path]:
    """Bucket key -> local file for every file uploaded at its own plain path, which since Plan 9 is
    the `assets` JSON plus gwas_dcm.json. `tables` entries are handled for a pre-Plan-9 manifest;
    a glob path there is a build intermediate that is not uploaded any more."""
    out: dict[str, Path] = {}
    for entry in man.get("tables", {}).values():
        if "*" not in entry["path"]:
            out[entry["path"]] = cfg.derived / entry["path"]
    for entry in man.get("assets", {}).values():
        out[entry["path"]] = cfg.derived / entry["path"]
    for name in PLAIN_EXTRA:
        if (cfg.derived / name).exists():
            out[name] = cfg.derived / name
    return dict(sorted(out.items()))


def _expand(cfg: Config, pattern: str) -> list[str]:
    """A `tables` glob from a live manifest, resolved against the local trees that carry the bucket's
    layout. Returns bucket keys."""
    found: set[str] = set()
    for root in LOCAL_ROOTS:
        base = cfg.derived / root if root else cfg.derived
        if base.is_dir():
            found |= {str(p.relative_to(base)) for p in base.glob(pattern) if p.is_file()}
    return sorted(found)


def _named_keys(cfg: Config, man: dict, warn: bool = False) -> set[str]:
    """Every bucket key a manifest names. Handles both shapes: the live pre-Plan-9 manifest with a
    `tables` block of plain paths and globs, and the current one with `immutable` and `packs`."""
    keys: set[str] = set(man.get("immutable", {}))
    packs = man.get("packs", {})
    for name, value in packs.items():
        if name == "files":
            for kind in value.values():
                keys |= set(kind.values())
        elif isinstance(value, str) and "/" in value:
            keys.add(value)
    for name, entry in man.get("tables", {}).items():
        path = entry.get("path", "")
        if "*" not in path:
            keys.add(path)
            continue
        hits = _expand(cfg, path)
        keys |= set(hits)
        if warn and entry.get("files") is not None and len(hits) != entry["files"]:
            print(f"WARNING: {name}: manifest says {entry['files']} files, {len(hits)} expanded locally from {path!r}")
    for entry in man.get("assets", {}).values():
        keys.add(entry["path"])
    keys.add(MANIFEST)
    keys |= set(PLAIN_EXTRA)
    return keys


def _retire_prefixes(cfg: Config, man: dict) -> list[str]:
    """Old bucket paths the content-addressed files take over: the top-level `replaces:` map in the
    config, formatted per chromosome, plus the per-entry `replaces` the manifest already recorded.
    Some are whole prefixes (`gene_detail/chr=chr1/`), some exact keys (`search_index.parquet`)."""
    out: set[str] = set()
    for prefixes in cfg["replaces"].values():        # `_retired` is a plain list, like any other kind
        for pre in prefixes:
            if "{chr}" in pre:
                out |= {pre.format(chr=c) for c in CHROMS}
            else:
                out.add(pre)
    for entry in man.get("immutable", {}).values():
        out |= set(entry.get("replaces") or [])
    return sorted(out)


def _under(key: str, prefixes: list[str]) -> str | None:
    """The first retire prefix a key sits under, or that it exactly equals."""
    for pre in prefixes:
        if key == pre or (pre.endswith("/") and key.startswith(pre)):
            return pre
    return None


def _local_copy(cfg: Config, key: str, size: int) -> Path | None:
    """The same relative path and size under one of LOCAL_ROOTS, so a deleted key can be restored."""
    for root in LOCAL_ROOTS:
        p = (cfg.derived / root / key) if root else (cfg.derived / key)
        if p.is_file() and p.stat().st_size == size:
            return p
    return None


def _staged_manifest_key(cfg: Config) -> tuple[str, bytes]:
    """`immutable/manifest.<sha16>.json` for the local manifest, and its bytes. Staging this copy
    first is what lets `check --staged` rehearse the exact bytes `release` will make live."""
    raw = (cfg.derived / MANIFEST).read_bytes()
    return f"{cfg['r2']['immutable_prefix']}manifest.{hashlib.sha256(raw).hexdigest()[:16]}.json", raw


def _no_listing_with_write(listing, what: str) -> None:
    if listing:
        sys.exit(f"--listing is a saved snapshot; {what} refuses it. Drop --listing, or add --dryrun.")


# ---- inventory -----------------------------------------------------------------------------------

def _write_listing(out: Path, keys: dict[str, dict], source: str) -> None:
    body = {"Source": source,
            "Generated": dt.datetime.now().isoformat(timespec="seconds"),
            "Contents": [{"Key": k, "Size": v["size"], "ETag": f'"{v["etag"]}"'} for k, v in sorted(keys.items())]}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(body, indent=1))
    total = sum(v["size"] for v in keys.values())
    print(f"{len(keys):,} keys, {total:,} B ({total / 1e9:.3f} GB) -> {out}")


def inventory(cfg: Config, out: Path, public: bool = False) -> None:
    if not public:
        _write_listing(out, _list(cfg), "list-objects-v2")
        return
    status, _, body = _public(cfg, "GET", MANIFEST)
    if status != 200:
        sys.exit(f"GET {MANIFEST} on the public URL -> {status}")
    live = json.loads(body)
    names = sorted(_named_keys(cfg, live, warn=True))
    print(f"live manifest names {len(names):,} keys (built {live.get('built')}, commit {str(live.get('pipeline_commit'))[:12]})")

    def head(key: str):
        st, h, _ = _public(cfg, "HEAD", key, {})
        return key, st, h

    keys: dict[str, dict] = {}
    missing = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for key, st, h in pool.map(head, names):
            if st != 200:
                missing.append((key, st))
                continue
            keys[key] = {"size": int(h.get("content-length", 0)), "etag": (h.get("etag") or "").strip('"')}
    for key, st in missing:
        print(f"WARNING: {key} -> {st}; dropped from the listing")
    _write_listing(out, keys, "public")
    print("note: a public listing can only see keys the live manifest names. Anything else in the "
          "bucket is invisible here; use `inventory` with the token for the whole bucket.")


# ---- budget ---------------------------------------------------------------------------------------

def _median(xs: list[float]) -> float | None:
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return None
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def _harness_requests(path: Path) -> tuple[float | None, int]:
    """Median `nav` to_idle.requests_data across the harness's pages: the data requests one in-app
    gene navigation costs, which is the per-gene part of a session's Class B bill."""
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, ValueError) as e:
        print(f"WARNING: {path}: {e}")
        return None, 0
    vals = []
    for scen in (data.get("results") or {}).values():
        req = (((scen.get("nav") or {}).get("median") or {}).get("to_idle") or {}).get("requests_data")
        if req is not None:
            vals.append(req)
    return _median(vals), len(vals)


def _plain_changed(cfg: Config, man: dict, keys: dict[str, dict]) -> tuple[list[tuple[str, Path, str]], list[str]]:
    """(key, local path, why) for every plain-path file to upload, and the keys whose local file is
    gone. At or under multipart_threshold_bytes the bucket ETag is the MD5 and is compared; above it
    the ETag is a multipart hash from the old `sync`, so only the size is."""
    thr = cfg["r2"]["multipart_threshold_bytes"]
    changed, missing = [], []
    for key, path in _plain_paths(cfg, man).items():
        if not path.exists():
            missing.append(key)
            continue
        size = path.stat().st_size
        have = keys.get(key)
        if have is None:
            changed.append((key, path, "new"))
        elif have["size"] != size:
            changed.append((key, path, f"size {have['size']:,} -> {size:,}"))
        elif size > thr:
            pass                                            # size matches and the ETag is not an MD5
        elif have["etag"] != digests(path)[1]:
            changed.append((key, path, "MD5 differs"))
    return changed, missing


def budget(cfg: Config, man: dict, keys: dict[str, dict], overlap_days: int = 3,
           harness: Path | None = None) -> dict:
    """Read-only storage and operations arithmetic for this deploy. Prints; returns the totals so
    `stage` can refuse on the same numbers."""
    tier = cfg["r2"]["free_tier"]
    imm = man.get("immutable", {})
    staged_key, staged_raw = _staged_manifest_key(cfg)

    current = sum(v["size"] for v in keys.values())
    new_keys = [k for k in imm if k not in keys]
    add = sum(imm[k]["bytes"] for k in new_keys) + (0 if staged_key in keys else len(staged_raw))
    peak = current + add
    prefixes = _retire_prefixes(cfg, man)
    retire_keys = [k for k in keys if _under(k, prefixes)]
    retire = sum(keys[k]["size"] for k in retire_keys)
    final = peak - retire

    print(f"bucket {cfg['r2']['bucket']}  ({len(keys):,} keys now, {len(imm):,} content-addressed files locally)")
    print(f"  current {_gb(current)}")
    print(f"  add     {_gb(add)}   ({len(new_keys):,} new immutable files + the staged manifest copy)")
    print(f"  peak    {_gb(peak)}")
    print(f"  retire  {_gb(retire)}   ({len(retire_keys):,} keys under {len(prefixes)} replaced prefixes)")
    print(f"  final   {_gb(final)}")
    print(f"  free tier storage {_gb(tier['storage_bytes'])}  -> final is {final / tier['storage_bytes']:.1%} of it")

    day = dt.date.today().day
    days_before = day - 1
    days_after = max(0, 30 - days_before - overlap_days)
    this_month = (days_before * current + overlap_days * peak + days_after * final) / 30
    print("\nGB-month (R2 averages each day's storage over 30 days)")
    print(f"  this month  {_gb(int(this_month))}   ({days_before} d at current + {overlap_days} d at peak "
          f"+ {days_after} d at final)")
    print(f"  later       {_gb(final)}")

    over = final >= tier["storage_bytes"] or this_month >= tier["storage_bytes"]
    if over:
        print("\nREFUSE: this deploy does not fit the free tier. Prune more, or shrink the upload set.")
    else:
        if peak > tier["storage_bytes"]:
            print(f"\nWARNING: peak is over the tier by {peak - tier['storage_bytes']:,} B. Allowed, because R2 "
                  f"averages each day's peak, but the overlap window must stay short.")
        if final > 0.9 * tier["storage_bytes"]:
            print("\nWARNING: final is above 90% of the tier; there is little room for the next build.")

    changed, missing_plain = _plain_changed(cfg, man, keys)
    stage_a = len(new_keys) + (0 if staged_key in keys else 1) + 1
    release_a = len(changed) + 1 + 1
    prune_a = 1 + -(-len(retire_keys) // 1000)
    total_a = stage_a + release_a + prune_a + 3
    print(f"\nClass A for the deploy (free tier {tier['class_a_ops']:,} a month)")
    print(f"  stage    {stage_a:>6,}   {len(new_keys):,} PUTs + the staged manifest copy + 1 list")
    print(f"  release  {release_a:>6,}   {len(changed):,} changed plain files + manifest.json + 1 list")
    print(f"  prune    {prune_a:>6,}   1 list + 1 DeleteObjects per 1,000 keys")
    print(f"  listings {3:>6,}   inventory, budget, check --retired")
    print(f"  total    {total_a:>6,}   {total_a / tier['class_a_ops']:.4%} of the month's Class A")
    if missing_plain:
        print(f"  WARNING: {len(missing_plain)} plain-path file(s) the manifest names are missing locally: {missing_plain}")

    print(f"\nClass B per visit (free tier {tier['class_b_ops']:,} a month)")
    sources = [("harness", harness)] if harness else []
    baseline = ROOT / "ui" / "bench" / "results" / "baseline-live.json"
    if baseline.exists():
        sources.append(("baseline-live", baseline))
    if not sources:
        print("  no harness JSON given and ui/bench/results/baseline-live.json is missing")
    for label, path in sources:
        R, pages = _harness_requests(path)
        if R is None:
            print(f"  {label}: no `nav` medians in {path}")
            continue
        print(f"  {label}: median nav data requests R = {R:g} over {pages} pages ({Path(path).name})")
        for genes in (1, 5):
            per = 1 + 1 + 2 + genes * R                     # manifest, search_index, 2 landing JSON, genes x R
            print(f"    {genes} gene(s): {per:g} requests a session -> {tier['class_b_ops'] / per:,.0f} sessions a month")
    print("  manifest revalidations (304) and CORS preflights count as Class B too; `check` is about one HEAD "
          f"per file ({len(imm):,}), and `inventory --public` about {len(_named_keys(cfg, man)):,}.")

    return {"current": current, "add": add, "peak": peak, "retire": retire, "retire_keys": len(retire_keys),
            "final": final, "this_month": this_month, "over": over, "new_keys": new_keys,
            "staged_key": staged_key, "changed_plain": changed, "missing_plain": missing_plain}


def cmd_budget(cfg: Config, listing: Path | None, harness: Path | None, overlap_days: int) -> None:
    man = _manifest(cfg)
    b = budget(cfg, man, _list(cfg, listing), overlap_days, harness)
    sys.exit(1 if b["over"] else 0)


# ---- stage ------------------------------------------------------------------------------------------

def stage(cfg: Config, dryrun: bool = False, listing: Path | None = None, overlap_days: int = 3) -> None:
    """Upload every new `immutable/` file plus a content-addressed copy of the manifest. Nothing the
    live site reads changes: the live manifest.json still names the old files."""
    if not dryrun:
        _no_listing_with_write(listing, "stage")
    man = _manifest(cfg)
    keys = _list(cfg, listing)
    b = budget(cfg, man, keys, overlap_days)
    if b["over"]:
        sys.exit(1)
    print()

    imm, refusals, todo = man.get("immutable", {}), [], []
    header = cfg["r2"]["cache_control"]["immutable"]
    for key in sorted(imm):
        entry, have, path = imm[key], keys.get(key), cfg.derived / key
        if not path.exists():
            refusals.append(f"{key}: no local file (run the step that writes {entry['key']})")
            continue
        if have is None:
            todo.append(key)
            print(f"  new      {key:<60} {entry['bytes']:>13,} B  {_content_type(key)}  {header!r}")
        elif have["size"] != entry["bytes"] or have["etag"] != entry["md5"]:
            refusals.append(f"{key}: already in the bucket with a different size or ETag "
                            f"(bucket {have['size']:,} B / {have['etag']}, local {entry['bytes']:,} B / {entry['md5']}); "
                            f"immutable keys are never rewritten")
        else:
            print(f"  present  {key:<60} {entry['bytes']:>13,} B")

    staged_key, staged_raw = _staged_manifest_key(cfg)
    staged_new = staged_key not in keys
    print(f"  {'new     ' if staged_new else 'present '} {staged_key:<60} {len(staged_raw):>13,} B  application/json  {header!r}")
    print(f"\nstaged manifest copy: {staged_key}")
    add = sum(imm[k]['bytes'] for k in todo) + (len(staged_raw) if staged_new else 0)
    print(f"{len(todo):,} new files + {'1' if staged_new else '0'} manifest copy, {add:,} B ({add / 1e9:.3f} GB); "
          f"{len(imm) - len(todo):,} already present")

    for r in refusals:
        print(f"WOULD REFUSE: {r}" if dryrun else f"REFUSE: {r}")
    if dryrun:
        print("\ndry run: nothing was uploaded")
        return
    if refusals:
        sys.exit(1)

    staged_path = cfg.tmp / "manifest.staged.json"
    staged_path.parent.mkdir(parents=True, exist_ok=True)
    staged_path.write_bytes(staged_raw)

    def send(key: str) -> None:
        _put(cfg, key, cfg.derived / key, imm[key]["md5"], header, _content_type(key))

    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(send, todo))
    if staged_new:
        _put(cfg, staged_key, staged_path, hashlib.md5(staged_raw).hexdigest(), header, "application/json")
    print(f"staged {len(todo):,} files; rehearse with `upload.py check --staged`, then `release --dryrun`")


# ---- release and pointer ------------------------------------------------------------------------------

def release(cfg: Config, dryrun: bool = False, listing: Path | None = None) -> None:
    """The live switch: changed plain-path files, then manifest.json last."""
    if not dryrun:
        _no_listing_with_write(listing, "release")
    man = _manifest(cfg)
    keys = _list(cfg, listing)
    imm = man.get("immutable", {})
    refusals = []

    absent = [k for k in sorted(imm) if k not in keys]
    if absent:
        refusals.append(f"{len(absent)} immutable file(s) are not in the bucket; run `stage` first "
                        f"(first: {absent[:3]})")
    staged_key, staged_raw = _staged_manifest_key(cfg)
    if staged_key not in keys:
        refusals.append(f"the staged copy of this manifest ({staged_key}) is not in the bucket; what goes "
                        f"live must be exactly what was staged and rehearsed")

    changed, missing_plain = _plain_changed(cfg, man, keys)
    for key in missing_plain:
        refusals.append(f"{key}: the manifest names it but there is no local file")
    mutable = cfg["r2"]["cache_control"]["mutable"]
    print(f"plain-path files ({len(_plain_paths(cfg, man))} named, {len(changed)} to upload)")
    for key, path, why in changed:
        print(f"  upload   {key:<40} {path.stat().st_size:>10,} B  {_content_type(key)}  {mutable!r}  ({why})")
    for key in _plain_paths(cfg, man):
        if key not in {k for k, _, _ in changed} and key not in missing_plain:
            print(f"  same     {key}")
    print(f"  last     {MANIFEST:<40} {len(staged_raw):>10,} B  application/json  {mutable!r}  "
          f"(sha256 {hashlib.sha256(staged_raw).hexdigest()[:16]})")

    for r in refusals:
        print(f"WOULD REFUSE: {r}" if dryrun else f"REFUSE: {r}")
    if dryrun:
        print("\ndry run: nothing was uploaded")
        return
    if refusals:
        sys.exit(1)

    for key, path, _ in changed:
        _put(cfg, key, path, digests(path)[1], mutable, _content_type(key))
    _put(cfg, MANIFEST, cfg.derived / MANIFEST, hashlib.md5(staged_raw).hexdigest(), mutable, "application/json")
    print(f"released: {len(changed)} plain files, then {MANIFEST}. Verify with `upload.py check`.")


def pointer(cfg: Config, path: Path) -> None:
    """Rollback's bucket half: make FILE the live manifest. FILE is usually a staged copy fetched
    back out of `immutable/`."""
    path = Path(path)
    raw = path.read_bytes()
    try:
        json.loads(raw)
    except ValueError as e:
        sys.exit(f"{path}: not JSON ({e})")
    sha = hashlib.sha256(raw).hexdigest()
    _put(cfg, MANIFEST, path, hashlib.md5(raw).hexdigest(), cfg["r2"]["cache_control"]["mutable"], "application/json")
    print(f"{MANIFEST} is now {path} ({len(raw):,} B, sha256 {sha})")


# ---- check ---------------------------------------------------------------------------------------------

def _slice(path: Path, a: int, b: int) -> bytes | None:
    """Local bytes [a, b] without reading a 200 MB pack into memory."""
    if not path.exists():
        return None
    with path.open("rb") as fh:
        fh.seek(a)
        return fh.read(b - a + 1)


def _ranges(size: int) -> list[tuple[int, int]]:
    """The first and last 16 bytes, as SPEC's readers address them; one range for a tiny file."""
    first = (0, min(15, size - 1))
    return [first] if size < 32 else [first, (size - 16, size - 1)]


def check(cfg: Config, staged: bool = False, retired: bool = False, listing: Path | None = None) -> None:
    """Public-URL checks, from a browser-like origin. One line per check; exits 1 on any failure."""
    man = _manifest(cfg)
    origin = cfg["r2"]["allowed_origins"][0]
    imm = man.get("immutable", {})
    fails: list[tuple[str, str]] = []                       # (category, message)

    def bad(category: str, msg: str) -> None:
        fails.append((category, msg))
        print(f"  FAIL {msg}")

    if retired:
        keys = _list(cfg, listing)
        prefixes = _retire_prefixes(cfg, man)
        left = sorted(k for k in keys if _under(k, prefixes))
        print(f"retired prefixes ({len(prefixes)} of them, {len(keys):,} keys in the listing)")
        if left:
            bad("bucket", f"{len(left):,} keys still under a replaced prefix "
                          f"({sum(keys[k]['size'] for k in left):,} B); first: {left[:3]}")
        else:
            print("  ok   no key sits under a replaced prefix")
        status, _, _ = _public(cfg, "HEAD", "search_index.parquet")
        if status == 404:
            print("  ok   search_index.parquet -> 404")
        else:
            bad("bucket", f"search_index.parquet -> {status}, expected 404")
        _finish(fails)
        return

    targets = dict(imm)
    if staged:
        staged_key, staged_raw = _staged_manifest_key(cfg)
        targets[staged_key] = {"bytes": len(staged_raw), "md5": hashlib.md5(staged_raw).hexdigest(),
                               "key": "manifest (staged copy)"}
    want_cc = cfg["r2"]["cache_control"]["immutable"]
    print(f"HEAD {len(targets):,} immutable keys on {cfg['r2']['public_url']}")

    def head(key: str):
        st, h, _ = _public(cfg, "HEAD", key, {"Origin": origin})
        return key, st, h

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(head, sorted(targets)))
    dirty = 0
    for key, st, h in results:
        entry, n = targets[key], len(fails)
        if st != 200:
            bad("status", f"{key} -> {st}")
            dirty += 1
            continue
        if int(h.get("content-length", -1)) != entry["bytes"]:
            bad("length", f"{key}: content-length {h.get('content-length')}, manifest says {entry['bytes']:,}")
        if (h.get("etag") or "").strip('"') != entry["md5"]:
            bad("etag", f"{key}: ETag {h.get('etag')}, manifest MD5 {entry['md5']}")
        if h.get("cache-control") != want_cc:
            bad("cache-control", f"{key}: cache-control {h.get('cache-control')!r}, want {want_cc!r}")
        dirty += len(fails) > n
    print(f"  {len(targets) - dirty:,} of {len(targets):,} keys clean (200, length, ETag, cache-control)")

    # SPEC's readers only ever range-read, so prove a range works on the index and the biggest pack
    ranged = []
    si = (man.get("packs") or {}).get("search_index")
    if si in imm:
        ranged.append(si)
    packs = {k: v for k, v in imm.items() if not v["key"].startswith("search_index")}
    if packs:
        ranged.append(max(packs, key=lambda k: packs[k]["bytes"]))
    for key in ranged:
        size = imm[key]["bytes"]
        for a, b in _ranges(size):
            st, h, body = _public(cfg, "GET", key, {"Origin": origin, "Range": f"bytes={a}-{b}"})
            label, n = f"{key} bytes={a}-{b}", len(fails)
            if st != 206:
                bad("status", f"{label} -> {st}, expected 206")
                continue
            want_cr = f"bytes {a}-{b}/{size}"
            if h.get("content-range") != want_cr:
                bad("content-range", f"{label}: content-range {h.get('content-range')!r}, want {want_cr!r}")
            local = _slice(cfg.derived / key, a, b)
            if local is not None and body != local:
                bad("body", f"{label}: {len(body)} bytes do not match the local file")
            if h.get("content-encoding"):
                bad("content-encoding", f"{label}: content-encoding {h.get('content-encoding')!r}; byte offsets "
                                        f"address the stored bytes, so there must be none")
            acao = h.get("access-control-allow-origin")
            if acao not in (origin, "*"):
                bad("cors", f"{label}: access-control-allow-origin {acao!r}, want {origin!r} or '*'")
            expose = (h.get("access-control-expose-headers") or "").lower()
            miss = [x for x in EXPOSE if x.lower() not in expose]
            if miss:
                bad("cors", f"{label}: access-control-expose-headers is missing {miss}")
            if len(fails) == n:
                print(f"  ok   {label} -> 206, {h.get('content-range')}")

    if not staged:
        raw = (cfg.derived / MANIFEST).read_bytes()
        st, h, body = _public(cfg, "GET", MANIFEST, {"Origin": origin})
        n = len(fails)
        if st != 200:
            bad("status", f"{MANIFEST} -> {st}")
        else:
            if "no-cache" not in (h.get("cache-control") or ""):
                bad("cache-control", f"{MANIFEST}: cache-control {h.get('cache-control')!r} has no no-cache")
            if h.get("access-control-allow-origin") not in (origin, "*"):
                bad("cors", f"{MANIFEST}: access-control-allow-origin {h.get('access-control-allow-origin')!r}")
            got, want = hashlib.sha256(body).hexdigest(), hashlib.sha256(raw).hexdigest()
            if got != want:
                bad("body", f"{MANIFEST}: live sha256 {got[:16]} != local {want[:16]}; the live site is not on this build")
            elif len(fails) == n:
                print(f"  ok   {MANIFEST} -> 200, sha256 {got[:16]}, cache-control {h.get('cache-control')!r}")
    _finish(fails)


def _finish(fails: list[tuple[str, str]]) -> None:
    if not fails:
        print("check: OK")
        sys.exit(0)
    if all(c == "cors" for c, _ in fails):
        print("CORS is Sam's: Admin token or dashboard (`upload.py cors --print`)")
    print(f"check: FAILED ({len(fails)} problems)")
    sys.exit(1)


# ---- prune -------------------------------------------------------------------------------------------

def prune(cfg: Config, stale: bool = False, yes: bool = False, listing: Path | None = None) -> None:
    """Delete the bucket paths the content-addressed files took over. Dry run unless --yes."""
    if yes:
        _no_listing_with_write(listing, "prune --yes")
    man = _manifest(cfg)
    keys = _list(cfg, listing)
    prefixes = _retire_prefixes(cfg, man)
    referenced = _named_keys(cfg, man)
    staged_key, staged_raw = _staged_manifest_key(cfg)
    referenced.add(staged_key)

    groups: dict[str, list[str]] = {}
    for key in sorted(keys):
        if key in referenced:
            continue
        pre = _under(key, prefixes)
        if pre:
            groups.setdefault(pre, []).append(key)
    if stale:
        imm_pre = cfg["r2"]["immutable_prefix"]
        extra = [k for k in sorted(keys) if k.startswith(imm_pre) and k not in referenced]
        if extra:
            groups[f"{imm_pre} (stale builds)"] = extra
    candidates = [k for group in groups.values() for k in group]
    left_alone = sorted(k for k in keys if k not in referenced and k not in set(candidates))

    total = sum(keys[k]["size"] for k in candidates)
    print(f"prune candidates: {len(candidates):,} keys, {total:,} B ({total / 1e9:.3f} GB)")
    no_copy = []
    for pre in sorted(groups):
        group = groups[pre]
        gb = sum(keys[k]["size"] for k in group)
        copies = {k: _local_copy(cfg, k, keys[k]["size"]) for k in group}
        no_copy += [k for k, v in copies.items() if v is None]
        print(f"\n  {pre}")
        print(f"    {len(group):,} keys, {gb:,} B ({gb / 1e9:.3f} GB); "
              f"{sum(1 for v in copies.values() if v):,}/{len(group):,} have a local copy")
        head, tail = (group, []) if len(group) <= 10 else (group[:5], group[-5:])
        for k in head:
            print(f"    {'local  ' if copies[k] else 'NO COPY'} {k}  {keys[k]['size']:,} B")
        if tail:
            print(f"    ... {len(group) - 10:,} more ...")
            for k in tail:
                print(f"    {'local  ' if copies[k] else 'NO COPY'} {k}  {keys[k]['size']:,} B")
    if left_alone:
        la = sum(keys[k]["size"] for k in left_alone)
        print(f"\n  left alone ({len(left_alone):,} keys, {la:,} B): not referenced, not replaced; Sam decides")
        for k in left_alone[:20]:
            print(f"    {k}  {keys[k]['size']:,} B")
        if len(left_alone) > 20:
            print(f"    ... {len(left_alone) - 20:,} more")

    if not yes:
        print("\ndry run: nothing was deleted. Add --yes after `check` passes on the live site.")
        return

    refusals = []
    st, _, body = _public(cfg, "GET", MANIFEST)
    live = hashlib.sha256(body).hexdigest() if st == 200 else None
    local = hashlib.sha256(staged_raw).hexdigest()
    if live != local:
        refusals.append(f"the live {MANIFEST} (status {st}, sha256 {str(live)[:16]}) is not this build "
                        f"({local[:16]}); pruning before release deletes files the live site still reads")
    for key in list(_plain_paths(cfg, man)) + [MANIFEST]:
        pre = _under(key, prefixes)
        if pre:
            refusals.append(f"{key} sits under the replaced prefix {pre}; the prefix list is wrong")
    if no_copy:
        refusals.append(f"{len(no_copy):,} candidate(s) have no local copy; fetch them before deleting, so a "
                        f"rollback stays possible (first: {no_copy[:3]})")
    if refusals:
        for r in refusals:
            print(f"REFUSE: {r}")
        sys.exit(1)
    _delete(cfg, candidates)
    print(f"deleted {len(candidates):,} keys, {total:,} B")


# ---- cors ---------------------------------------------------------------------------------------------

def cors(cfg: Config, print_only: bool = False) -> None:
    rule = {"CORSRules": [{
        "AllowedOrigins": cfg["r2"]["allowed_origins"],
        "AllowedMethods": ["GET", "HEAD"],
        "AllowedHeaders": ["Range", "If-Match", "If-None-Match", "Content-Type"],
        "ExposeHeaders": EXPOSE,
        "MaxAgeSeconds": 86400,
    }]}
    if print_only:
        # paste into dashboard > R2 > bucket > Settings > CORS policy (the dashboard takes the rules array)
        print(json.dumps(rule["CORSRules"], indent=2))
        return
    r = _aws(cfg, "s3api", "put-bucket-cors", "--bucket", cfg["r2"]["bucket"], "--cors-configuration", json.dumps(rule),
             capture_output=True, text=True)
    if r.returncode:
        if "AccessDenied" in r.stderr:
            sys.exit("PutBucketCors needs an Admin Read & Write token; an Object Read & Write token cannot set it. "
                     "Either use an admin token here or paste `upload.py cors --print` into the dashboard CORS policy.")
        sys.exit(r.stderr.strip())
    r = _aws(cfg, "s3api", "get-bucket-cors", "--bucket", cfg["r2"]["bucket"], capture_output=True, text=True)
    print(r.stdout)


# ---- cli ----------------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    i = sub.add_parser("inventory", help="read-only bucket listing, saved for --listing")
    i.add_argument("--public", action="store_true", help="no token: GET the live manifest, then HEAD every key it names")
    i.add_argument("-o", "--out", required=True, type=Path, metavar="FILE", help="where to write the listing JSON")

    b = sub.add_parser("budget", help="storage and operations arithmetic for this deploy; read-only")
    b.add_argument("--listing", type=Path, metavar="FILE", help="a saved `inventory` instead of listing the bucket")
    b.add_argument("--harness", type=Path, metavar="JSON", help="a bench harness summary, for the Class B estimate")
    b.add_argument("--overlap-days", type=int, default=3, metavar="N", help="days the old and new files both exist (default 3)")

    s = sub.add_parser("stage", help="upload new immutable/ files and a staged manifest copy; nothing live changes")
    s.add_argument("--dryrun", action="store_true")
    s.add_argument("--listing", type=Path, metavar="FILE", help="a saved `inventory` (refused without --dryrun)")
    s.add_argument("--overlap-days", type=int, default=3, metavar="N")

    r = sub.add_parser("release", help="upload changed plain-path files, then manifest.json: the live switch")
    r.add_argument("--dryrun", action="store_true")
    r.add_argument("--listing", type=Path, metavar="FILE", help="a saved `inventory` (refused without --dryrun)")

    c = sub.add_parser("check", help="public-URL checks; exits 1 on any failure")
    g = c.add_mutually_exclusive_group()
    g.add_argument("--staged", action="store_true", help="the staged files only, before release")
    g.add_argument("--retired", action="store_true", help="the replaced prefixes are empty")
    c.add_argument("--listing", type=Path, metavar="FILE", help="a saved `inventory`, for --retired")

    p = sub.add_parser("pointer", help="rollback: upload FILE as manifest.json")
    p.add_argument("file", type=Path)

    pr = sub.add_parser("prune", help="delete the replaced bucket paths; dry run unless --yes")
    pr.add_argument("--stale", action="store_true", help="also immutable/ keys no manifest references")
    pr.add_argument("--yes", action="store_true", help="actually delete")
    pr.add_argument("--listing", type=Path, metavar="FILE", help="a saved `inventory` (refused with --yes)")

    co = sub.add_parser("cors", help="apply the CORS rule from config")
    co.add_argument("--print", action="store_true", dest="print_only",
                    help="print the rules for pasting into the dashboard instead of applying")

    a = ap.parse_args()
    cfg = Config()
    if a.cmd == "inventory":
        inventory(cfg, a.out, a.public)
    elif a.cmd == "budget":
        cmd_budget(cfg, a.listing, a.harness, a.overlap_days)
    elif a.cmd == "stage":
        stage(cfg, a.dryrun, a.listing, a.overlap_days)
    elif a.cmd == "release":
        release(cfg, a.dryrun, a.listing)
    elif a.cmd == "check":
        check(cfg, a.staged, a.retired, a.listing)
    elif a.cmd == "pointer":
        pointer(cfg, a.file)
    elif a.cmd == "prune":
        prune(cfg, a.stale, a.yes, a.listing)
    elif a.cmd == "cors":
        cors(cfg, a.print_only)


if __name__ == "__main__":
    main()
