"""Put a qtlstore in the bucket the browsers read, in the one order that is safe, and check it.

    uv run python -m pipeline.upload inventory [-o FILE]                 # bucket listing under the prefix, read-only
    uv run python -m pipeline.upload plan  --store DIR [--listing FILE]  # what push would do; nothing changes
    uv run python -m pipeline.upload push  --store DIR [--yes]           # dry run unless --yes
    uv run python -m pipeline.upload check --store DIR [--sample N]      # the public URL serves what DIR holds
    uv run python -m pipeline.upload prune --store DIR [--yes]           # delete keys DIR does not name; dry run unless --yes

The store's own layout decides the order (SPEC section 2). Objects under `immutable/` are named by
their content digest, so the browser caches them forever and an object already in the bucket is
never re-sent; pointer documents and `store.json` are mutable and served `no-cache`. So a release
is: every new object, then every pointer document, then `store.json` **last**. Until that final
PUT the live store is the old one, with every old object still present, and a reader that opened it
a minute ago keeps working. There is no moment at which a pointer names an object the bucket does
not hold.

`prune` is a separate command and never part of `push`: it deletes only keys under the prefix that
the local store does not name, and only after confirming the public `store.json` is the local one,
so it can never remove an object the live store still points at. The v1 store and its GTEx
experiment go that way once v2 is live.

`--listing FILE` reuses a saved `inventory` for `plan`. It is refused with `--yes`, so a stale
listing can never drive a real upload or delete.

Profiles and the key prefix are `storage:` in pipeline/config.yaml. Credentials come from the
environment or the gitignored .env; this module never prints them, and resolves them before it
logs a command so a missing token never looks like an attempted call. Uses the aws CLI against the
profile's S3 endpoint; the two checksum variables stop aws-cli 2.23+ from sending CRC headers B2
rejects.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

from . import qtlstore as qs
from .common import ROOT, Config, log

STORE_JSON = "store.json"
IMMUTABLE_CACHE = "public, max-age=31536000, immutable"
POINTER_CACHE = "no-cache"
UA = "Mozilla/5.0 (qtl-browser upload.py)"
_MD5 = re.compile(r"^[0-9a-f]{32}$")


# ---- the storage profile ------------------------------------------------------------------------

def prof(cfg: Config) -> dict:
    """The active storage profile merged over the shared `storage:` settings, plus `name`."""
    st = cfg.cfg.get("storage")
    if not st:
        sys.exit("pipeline/config.yaml has no `storage:` section")
    name = st.get("active") or st["default_profile"]
    if name not in st["profiles"]:
        sys.exit(f"storage profile {name!r} is not one of {sorted(st['profiles'])}")
    shared = {k: v for k, v in st.items() if k not in ("profiles", "active", "default_profile")}
    return {**shared, **st["profiles"][name], "name": name}


def _dotenv() -> dict[str, str]:
    """KEY=VALUE lines from <repo>/.env (gitignored); the shell environment wins over it."""
    path = ROOT / ".env"
    out: dict[str, str] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                out[k.strip().removeprefix("export ").strip()] = v.strip().strip("'\"")
    return out


def _region(p: dict) -> str:
    """B2's endpoint names its region and wants it in the signature; anything else takes `auto`."""
    if p.get("region"):
        return p["region"]
    m = re.match(r"https://s3\.([a-z0-9-]+)\.backblazeb2\.com", p["endpoint"])
    return m.group(1) if m else "auto"


def _env(cfg: Config) -> dict[str, str]:
    p = prof(cfg)
    var_id, var_secret = f"{p['name'].upper()}_ACCESS_KEY_ID", f"{p['name'].upper()}_SECRET_ACCESS_KEY"
    merged = {**_dotenv(), **os.environ}
    key, secret = merged.get(var_id), merged.get(var_secret)
    if not (key and secret):
        sys.exit(f"put {var_id} and {var_secret} in .env at the repo root or export them")
    return {**os.environ, "AWS_ACCESS_KEY_ID": key, "AWS_SECRET_ACCESS_KEY": secret,
            "AWS_DEFAULT_REGION": _region(p),
            "AWS_REQUEST_CHECKSUM_CALCULATION": "when_required",
            "AWS_RESPONSE_CHECKSUM_VALIDATION": "when_required"}


def _aws(cfg: Config, *args: str, **kw) -> subprocess.CompletedProcess:
    env = _env(cfg)
    cmd = ["aws", "--endpoint-url", prof(cfg)["endpoint"], *args]
    return subprocess.run(cmd, env=env, **kw)


def _bkey(cfg: Config, key: str) -> str:
    return prof(cfg).get("key_prefix", "") + key


# ---- the four bucket primitives -----------------------------------------------------------------

def inventory(cfg: Config) -> dict[str, dict]:
    """Store key -> {size, etag} for everything under the profile's prefix. Keys outside it belong
    to other projects sharing the bucket and are never returned, so nothing here can touch them."""
    pre = prof(cfg).get("key_prefix", "")
    args = ["s3api", "list-objects-v2", "--bucket", prof(cfg)["bucket"], "--output", "json"]
    if pre:
        args += ["--prefix", pre]
    r = _aws(cfg, *args, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"list-objects-v2: {r.stderr.strip()}")
    data = json.loads(r.stdout) if r.stdout.strip() else {}
    out = {}
    for o in data.get("Contents", []):
        key = o["Key"][len(pre):] if o["Key"].startswith(pre) else None
        if key:
            out[key] = {"size": o["Size"], "etag": o["ETag"].strip('"')}
    return out


def load_listing(path: Path) -> dict[str, dict]:
    return json.loads(Path(path).read_text())


def put(cfg: Config, key: str, path: Path, md5_hex: str, cache_control: str) -> None:
    """One single-part PUT with Content-MD5, so the bucket refuses corrupted bytes and the ETag is
    the MD5. Never multipart: the ETag would stop being the MD5 and `check` could not compare."""
    p = prof(cfg)
    size = path.stat().st_size
    if size > p["single_put_max_bytes"]:
        sys.exit(f"{key}: {size} bytes is over single_put_max_bytes")
    r = _aws(cfg, "s3api", "put-object", "--bucket", p["bucket"], "--key", _bkey(cfg, key),
             "--body", str(path), "--content-type", content_type(key), "--cache-control", cache_control,
             "--content-md5", base64.b64encode(bytes.fromhex(md5_hex)).decode(),
             capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"{key}: {r.stderr.strip()}")
    if json.loads(r.stdout)["ETag"].strip('"') != md5_hex:
        sys.exit(f"{key}: bucket ETag does not match the local MD5")


def delete(cfg: Config, keys: list[str]) -> None:
    for i in range(0, len(keys), 1000):
        batch = keys[i:i + 1000]
        payload = {"Objects": [{"Key": _bkey(cfg, k)} for k in batch], "Quiet": True}
        r = _aws(cfg, "s3api", "delete-objects", "--bucket", prof(cfg)["bucket"], "--delete", json.dumps(payload),
                 capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"delete-objects: {r.stderr.strip()}")
        errors = (json.loads(r.stdout) if r.stdout.strip() else {}).get("Errors") or []
        if errors:
            sys.exit(f"delete-objects: {errors}")
        log(f"deleted {len(batch):,} keys")


def public(cfg: Config, method: str, path: str, headers: dict[str, str] | None = None) -> tuple[int, dict, bytes]:
    """HEAD or GET on the profile's public URL, which already carries the prefix as a path."""
    url = prof(cfg)["public_url"].rstrip("/") + "/" + path.lstrip("/")
    req = urllib.request.Request(url, method=method, headers={"User-Agent": UA, **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.status, {k.lower(): v for k, v in resp.headers.items()}, resp.read()
    except urllib.error.HTTPError as e:
        return e.code, {k.lower(): v for k, v in (e.headers or {}).items()}, e.read()
    except urllib.error.URLError as e:
        return 0, {}, str(e).encode()


# ---- what a store is, as keys -------------------------------------------------------------------

def content_type(key: str) -> str:
    """JSON documents, pointer or content-addressed, are JSON; everything else is opaque bytes the
    browser range-reads, and never gets a Content-Encoding, because offsets address the stored
    bytes."""
    return "application/json" if key.endswith(".json") else "application/octet-stream"


def cache_control(key: str) -> str:
    return IMMUTABLE_CACHE if key.startswith("immutable/") else POINTER_CACHE


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class Local:
    """A store as the bucket keys it occupies, in release order: objects, then pointer levels in
    `POINTER_DIRS` order, then `store.json`."""
    objects: dict[str, Path] = field(default_factory=dict)
    pointers: dict[str, Path] = field(default_factory=dict)
    store_json: Path | None = None

    def all(self) -> dict[str, Path]:
        return {**self.objects, **self.pointers, STORE_JSON: self.store_json}


def local_keys(store: qs.Store) -> Local:
    """Every key the store occupies. Objects are the ones a pointer names (`qs.referenced`), so a
    leftover from an earlier local build is not uploaded; pointer documents are the ones
    `store.json` lists, in its order."""
    st = json.loads((store.root / STORE_JSON).read_text())
    out = Local(store_json=store.root / STORE_JSON)
    for name in sorted(qs.referenced(store)):
        p = store.immutable / name
        if not p.exists():
            sys.exit(f"{store.root}: pointer names {name}, which immutable/ does not hold (validate first)")
        out.objects[f"immutable/{name}"] = p
    for lvl in qs.POINTER_DIRS:
        ids = st.get(lvl) or {}
        for pid in (ids if isinstance(ids, dict) else list(ids)):
            p = store.root / lvl / f"{pid}.json"
            if not p.exists():
                sys.exit(f"{store.root}: store.json lists {lvl}/{pid} but the file is missing")
            out.pointers[f"{lvl}/{pid}.json"] = p
    return out


@dataclass
class Plan:
    put_objects: list[str]
    skip_objects: list[str]
    put_pointers: list[str]
    skip_pointers: list[str]
    put_store_json: bool
    stale: list[str]
    bytes_to_put: int
    md5: dict[str, str]


def plan(local: Local, listing: dict[str, dict]) -> Plan:
    """What `push` would send, and what `prune` would remove.

    An object is skipped when the bucket holds its key with the same bytes: same MD5 when the
    ETag is one, else same size (a content-addressed name already pins the bytes, so size is a
    check on the upload rather than on identity). A pointer or `store.json` is sent whenever its
    bytes differ. `stale` is every key under the prefix the store does not name."""
    md5: dict[str, str] = {}
    put_o, skip_o, put_p, skip_p, nbytes = [], [], [], [], 0
    for key, path in local.objects.items():
        have = listing.get(key)
        if have and have["size"] == path.stat().st_size and (
                not _MD5.match(have["etag"]) or have["etag"] == md5.setdefault(key, md5_of(path))):
            skip_o.append(key)
        else:
            md5.setdefault(key, md5_of(path))
            put_o.append(key)
            nbytes += path.stat().st_size
    for key, path in local.pointers.items():
        md5[key] = md5_of(path)
        have = listing.get(key)
        if have and _MD5.match(have["etag"]) and have["etag"] == md5[key]:
            skip_p.append(key)
        else:
            put_p.append(key)
            nbytes += path.stat().st_size
    md5[STORE_JSON] = md5_of(local.store_json)
    have = listing.get(STORE_JSON)
    put_s = not (have and _MD5.match(have["etag"]) and have["etag"] == md5[STORE_JSON])
    if put_s:
        nbytes += local.store_json.stat().st_size
    stale = sorted(set(listing) - set(local.all()))
    return Plan(put_o, skip_o, put_p, skip_p, put_s, stale, nbytes, md5)


def describe(p: Plan, listing: dict[str, dict]) -> str:
    stale_bytes = sum(listing[k]["size"] for k in p.stale)
    return "\n".join([
        f"objects:   {len(p.put_objects):>6,} to put, {len(p.skip_objects):>6,} already in the bucket",
        f"pointers:  {len(p.put_pointers):>6,} to put, {len(p.skip_pointers):>6,} unchanged",
        f"store.json: {'to put' if p.put_store_json else 'unchanged'}",
        f"bytes to put: {p.bytes_to_put:,} ({p.bytes_to_put / 1e9:.2f} GB)",
        f"stale in bucket (prune candidates): {len(p.stale):,} keys, {stale_bytes:,} bytes ({stale_bytes / 1e9:.2f} GB)",
    ])


# ---- commands -----------------------------------------------------------------------------------

def cmd_inventory(cfg: Config, out: Path | None) -> int:
    inv = inventory(cfg)
    text = json.dumps(inv, indent=1)
    if out:
        out.write_text(text)
        log(f"{len(inv):,} keys -> {out}")
    else:
        print(text)
    return 0


def cmd_plan(cfg: Config, store: Path, listing: Path | None) -> int:
    inv = load_listing(listing) if listing else inventory(cfg)
    p = plan(local_keys(qs.Store(store)), inv)
    print(describe(p, inv))
    return 0


def cmd_push(cfg: Config, store: Path, yes: bool, workers: int) -> int:
    st = qs.Store(store)
    fails = st.validate()
    if fails:
        print("\n".join(fails[:20]))
        sys.exit(f"{store}: validate reports {len(fails)} failures; not uploading")
    for n in st.notes:
        log(f"validate note: {n}")
    local = local_keys(st)
    inv = inventory(cfg)
    p = plan(local, inv)
    print(describe(p, inv))
    if not yes:
        print("\ndry run: pass --yes to upload")
        return 0
    # objects first, in parallel: nothing live names them yet
    def send(key: str) -> None:
        put(cfg, key, local.objects[key], p.md5[key], IMMUTABLE_CACHE)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, _ in enumerate(ex.map(send, p.put_objects), 1):
            if i % 25 == 0 or i == len(p.put_objects):
                log(f"objects {i}/{len(p.put_objects)}")
    # then the pointer documents, in level order; still nothing live has changed
    for key in p.put_pointers:
        put(cfg, key, local.pointers[key], p.md5[key], POINTER_CACHE)
        log(f"put {key}")
    # store.json last: the switch
    if p.put_store_json:
        put(cfg, STORE_JSON, local.store_json, p.md5[STORE_JSON], POINTER_CACHE)
        log(f"put {STORE_JSON}: the live store is now {store}")
    return 0


def cmd_check(cfg: Config, store: Path, sample: int, seed: int) -> int:
    local = local_keys(qs.Store(store))
    bad = []
    status, _, body = public(cfg, "GET", STORE_JSON)
    if status != 200 or body != local.store_json.read_bytes():
        bad.append(f"{STORE_JSON}: public copy is {'missing' if status != 200 else 'not the local one'} (HTTP {status})")
    for key, path in local.pointers.items():
        status, _, body = public(cfg, "GET", key)
        if status != 200 or body != path.read_bytes():
            bad.append(f"{key}: HTTP {status}" + ("" if status != 200 else ", bytes differ"))
    rng = random.Random(seed)
    keys = list(local.objects)
    for key in rng.sample(keys, min(sample, len(keys))):
        status, h, _ = public(cfg, "HEAD", key)
        want = local.objects[key].stat().st_size
        if status != 200 or int(h.get("content-length", -1)) != want:
            bad.append(f"{key}: HTTP {status}, content-length {h.get('content-length')} for {want} bytes")
    # the reader range-reads every object, so one 206 is part of the check
    if keys:
        key = keys[0]
        status, h, body = public(cfg, "GET", key, {"Range": "bytes=0-63"})
        if status != 206 or len(body) != 64:
            bad.append(f"{key}: range request answered HTTP {status} with {len(body)} bytes; the reader needs 206")
    print("\n".join(bad) if bad else f"OK: {STORE_JSON}, {len(local.pointers)} pointers, "
                                      f"{min(sample, len(keys))} sampled objects, one range read")
    return 1 if bad else 0


def cmd_prune(cfg: Config, store: Path, yes: bool, listing: Path | None) -> int:
    if prof(cfg).get("read_only"):
        sys.exit(f"profile {prof(cfg)['name']} is read-only")
    local = local_keys(qs.Store(store))
    status, _, body = public(cfg, "GET", STORE_JSON)
    if status != 200 or body != local.store_json.read_bytes():
        sys.exit(f"the public {STORE_JSON} is not {store}'s; push first, then prune")
    inv = load_listing(listing) if listing else inventory(cfg)
    stale = sorted(set(inv) - set(local.all()))
    nbytes = sum(inv[k]["size"] for k in stale)
    print(f"{len(stale):,} stale keys, {nbytes:,} bytes ({nbytes / 1e9:.2f} GB)")
    for k in stale[:20]:
        print(f"  {k}")
    if len(stale) > 20:
        print(f"  ... {len(stale) - 20} more")
    if not yes:
        print("\ndry run: pass --yes to delete")
        return 0
    delete(cfg, stale)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--profile", help="storage profile in config.yaml; default storage.default_profile")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("inventory", help="bucket listing under the prefix (read-only)")
    s.add_argument("-o", "--out", type=Path)
    for name in ("plan", "push", "check", "prune"):
        s = sub.add_parser(name)
        s.add_argument("--store", required=True, type=Path)
        if name in ("plan", "prune"):
            s.add_argument("--listing", type=Path, help="a saved inventory instead of listing the bucket")
        if name in ("push", "prune"):
            s.add_argument("--yes", action="store_true", help="actually write; dry run otherwise")
        if name == "push":
            s.add_argument("--workers", type=int, default=4)
        if name == "check":
            s.add_argument("--sample", type=int, default=20)
            s.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    cfg = Config()
    if a.profile:
        cfg.cfg["storage"]["active"] = a.profile
        prof(cfg)
    if a.cmd == "inventory":
        return cmd_inventory(cfg, a.out)
    if a.cmd == "plan":
        return cmd_plan(cfg, a.store, a.listing)
    if a.cmd == "push":
        return cmd_push(cfg, a.store, a.yes, a.workers)
    if a.cmd == "check":
        return cmd_check(cfg, a.store, a.sample, a.seed)
    if a.cmd == "prune":
        if a.yes and a.listing:
            sys.exit("--listing is refused with --yes: a stale listing must not drive a delete")
        return cmd_prune(cfg, a.store, a.yes, a.listing)
    return 2


if __name__ == "__main__":
    sys.exit(main())
