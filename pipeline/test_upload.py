"""Tests for pipeline/upload.py against a fake bucket (SPEC.md section 3, Plan 6 part E).

    uv run python -m pipeline.test_upload

Plain asserts, as in test_packfmt.py. **Nothing here touches the network.** A dict replaces the four
functions every bucket operation goes through -- `_list`, `_put`, `_delete`, `_public` -- over a
temporary `data/derived` tree with a small manifest, so a mistake in a refusal rule fails here
instead of on Sam's bucket.
"""
from __future__ import annotations

import contextlib
import datetime as dt
import hashlib
import io
import json
import shutil
import sys
import tempfile
from pathlib import Path

from . import upload as up
from .common import Config

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


# ---- harness -------------------------------------------------------------------------------------

def run(fn, *a, **kw) -> tuple[int | None, str]:
    """Call a command, capturing its output and normalising `sys.exit("message")` to (1, output)."""
    buf = io.StringIO()
    status = None
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            fn(*a, **kw)
        except SystemExit as e:
            status = e.code
    out = buf.getvalue()
    if isinstance(status, str):
        out += "\n" + status
        status = 1
    return status, out


def quiet(fn, *a, **kw):
    """Call a function that prints, and keep the test output readable."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


class Bucket:
    """A dict pretending to be R2. `objects` is key -> body, cache-control, content-type; `puts` and
    `deleted` record what a command actually did."""

    def __init__(self):
        self.objects: dict[str, dict] = {}
        self.puts: list[str] = []
        self.deleted: list[str] = []
        self.mutate = None                   # (key, status, headers, body) -> (status, headers, body)

    def add(self, key: str, body: bytes, cache_control: str = "public, max-age=31536000, immutable",
            content_type: str = "application/octet-stream") -> None:
        self.objects[key] = {"body": body, "cache_control": cache_control, "content_type": content_type}

    # the four functions upload.py routes every bucket operation through
    def list(self, cfg, listing=None):
        if listing:
            data = json.loads(Path(listing).read_text())
            return {o["Key"]: {"size": o["Size"], "etag": o["ETag"].strip('"')} for o in data.get("Contents", [])}
        return {k: {"size": len(v["body"]), "etag": hashlib.md5(v["body"]).hexdigest()}
                for k, v in self.objects.items()}

    def put(self, cfg, key, path, md5_hex, cache_control, content_type):
        body = Path(path).read_bytes()
        assert hashlib.md5(body).hexdigest() == md5_hex, f"{key}: Content-MD5 does not match the body sent"
        assert len(body) <= cfg["r2"]["single_put_max_bytes"], f"{key}: over single_put_max_bytes"
        self.add(key, body, cache_control, content_type)
        self.puts.append(key)

    def delete(self, cfg, keys):
        for k in keys:
            self.objects.pop(k, None)
        self.deleted += list(keys)

    def public(self, cfg, method, path, headers=None):
        headers = headers or {}
        obj = self.objects.get(path)
        if obj is None:
            status, h, body = 404, {}, b""
        else:
            body = obj["body"]
            h = {"content-length": str(len(body)), "etag": f'"{hashlib.md5(body).hexdigest()}"',
                 "cache-control": obj["cache_control"], "content-type": obj["content_type"],
                 "accept-ranges": "bytes"}
            if headers.get("Origin"):
                h["access-control-allow-origin"] = headers["Origin"]
                h["access-control-expose-headers"] = ", ".join(up.EXPOSE)
            status = 200
            out = b"" if method == "HEAD" else body
            if headers.get("Range") and method == "GET":
                lo, hi = (int(x) for x in headers["Range"].removeprefix("bytes=").split("-"))
                out = body[lo:hi + 1]
                status = 206
                h["content-range"] = f"bytes {lo}-{hi}/{len(body)}"
                h["content-length"] = str(len(out))
            body = out
        if self.mutate:
            status, h, body = self.mutate(path, status, h, body)
        return status, h, body


@contextlib.contextmanager
def wired(bucket: Bucket):
    saved = {n: getattr(up, n) for n in ("_list", "_put", "_delete", "_public")}
    up._list, up._put, up._delete, up._public = bucket.list, bucket.put, bucket.delete, bucket.public
    try:
        yield bucket
    finally:
        for n, f in saved.items():
            setattr(up, n, f)


PACKS = {"eqtl/chr1": ("qbe", b"E" * 200, ["old_eqtl/chr={chr}/"]),
         "eqtl/chr2": ("qbe", b"F" * 100, ["old_eqtl/chr={chr}/"]),
         "search_index": ("arrow.zst", b"S" * 64, ["search_index.parquet"])}
ASSETS = {"coloc_loci.json": b'{"small":1}' + b" " * 21,       # 32 B: at or under the threshold -> MD5 rule
          "gwas_dcm_bins.json": b'{"big":1}' + b" " * 91}      # 100 B: over the threshold -> size rule
GWAS_DCM = b'{"variants": 3}'


def make_tree() -> tuple[Config, Path, dict]:
    """A temporary data/derived with three published files, the JSON assets, and a manifest."""
    root = Path(tempfile.mkdtemp(prefix="qtlb-upload-test-"))
    derived = root / "derived"
    cfg = Config()
    cfg.cfg = json.loads(json.dumps(cfg.cfg))                  # deep copy: the test rewrites r2 and replaces
    cfg.derived, cfg.immutable, cfg.tmp = derived, derived / "immutable", derived / "_tmp"
    cfg.tables, cfg.done_dir = derived / "_tables", derived / ".done"
    cfg.cfg["r2"]["multipart_threshold_bytes"] = 64            # so both comparison rules are exercised
    cfg.cfg["replaces"] = {"eqtl": ["old_eqtl/chr={chr}/"], "search_index": ["search_index.parquet"],
                           "_retired": ["retired_table.parquet"]}
    cfg.immutable.mkdir(parents=True)

    immutable = {}
    files = {}
    for key, (ext, body, replaces) in PACKS.items():
        sha = hashlib.sha256(body).hexdigest()
        name = f"{key.replace('/', '.')}.{sha[:16]}.{ext}"
        (cfg.immutable / name).write_bytes(body)
        bkey = f"immutable/{name}"
        immutable[bkey] = {"key": key, "bytes": len(body), "sha256": sha, "md5": hashlib.md5(body).hexdigest(),
                           "replaces": [p.format(chr=key.partition("/")[2]) for p in replaces]}
        if key.startswith("eqtl/"):
            files.setdefault("eqtl", {})[key.partition("/")[2]] = bkey
        else:
            files[key] = bkey
    for name, body in ASSETS.items():
        (derived / name).write_bytes(body)
    (derived / "gwas_dcm.json").write_bytes(GWAS_DCM)
    man = {"immutable": immutable,
           "packs": {"files": {"eqtl": files["eqtl"]}, "search_index": files["search_index"], "version": 0},
           "assets": {n: {"path": n, "bytes": len(b)} for n, b in ASSETS.items()},
           "gwas_dcm": json.loads(GWAS_DCM)}
    (derived / "manifest.json").write_text(json.dumps(man, indent=2, sort_keys=True))
    return cfg, root, man


def fill(cfg: Config, bucket: Bucket, man: dict, staged: bool = True, plain: bool = True) -> None:
    """Put the state `release` expects into the fake bucket."""
    for key, entry in man["immutable"].items():
        bucket.add(key, (cfg.derived / key).read_bytes(), cfg["r2"]["cache_control"]["immutable"],
                   up._content_type(key))
    if staged:
        skey, raw = up._staged_manifest_key(cfg)
        bucket.add(skey, raw, cfg["r2"]["cache_control"]["immutable"], "application/json")
    if plain:
        for name in list(ASSETS) + ["gwas_dcm.json"]:
            bucket.add(name, (cfg.derived / name).read_bytes(), "no-cache", "application/json")


# ---- helpers -------------------------------------------------------------------------------------

@case
def test_content_types():
    assert up._content_type("immutable/eqtl.chr1.aaaaaaaaaaaaaaaa.qbe") == "application/octet-stream"
    assert up._content_type("immutable/search_index.aaaaaaaaaaaaaaaa.arrow.zst") == "application/octet-stream"
    assert up._content_type("immutable/gwas_index.aaaaaaaaaaaaaaaa.bin") == "application/octet-stream"
    assert up._content_type("search_index.parquet") == "application/vnd.apache.parquet"
    assert up._content_type("manifest.json") == "application/json"


@case
def test_retire_prefixes():
    cfg, root, man = make_tree()
    try:
        pres = up._retire_prefixes(cfg, man)
        assert "old_eqtl/chr=chr1/" in pres and "old_eqtl/chr=chrX/" in pres, pres
        assert "search_index.parquet" in pres and "retired_table.parquet" in pres, pres
        assert up._under("old_eqtl/chr=chr1/data.parquet", pres) == "old_eqtl/chr=chr1/"
        assert up._under("search_index.parquet", pres) == "search_index.parquet"
        assert up._under("search_index.parquet.bak", pres) is None   # an exact key is not a prefix
        assert up._under("manifest.json", pres) is None
    finally:
        shutil.rmtree(root)


# ---- stage ---------------------------------------------------------------------------------------

@case
def test_stage_uploads_only_missing():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        with wired(b):
            status, out = run(up.stage, cfg)
        assert status is None, out
        skey, raw = up._staged_manifest_key(cfg)
        assert set(b.puts) == set(man["immutable"]) | {skey}, b.puts
        assert skey in out, out
        for key in man["immutable"]:
            assert b.objects[key]["cache_control"] == cfg["r2"]["cache_control"]["immutable"]
            assert b.objects[key]["content_type"] == "application/octet-stream"
        assert b.objects[skey]["content_type"] == "application/json"
        assert b.objects[skey]["body"] == raw

        b.puts.clear()                                   # re-running stages nothing
        with wired(b):
            status, out = run(up.stage, cfg)
        assert status is None and b.puts == [], (status, b.puts)
        assert "0 new files + 0 manifest copy" in out, out
    finally:
        shutil.rmtree(root)


@case
def test_stage_refuses_a_different_file_at_an_immutable_key():
    cfg, root, man = make_tree()
    try:
        key = next(iter(sorted(man["immutable"])))
        for label, body in (("size", b"X" * 7), ("etag", b"X" * man["immutable"][key]["bytes"])):
            b = Bucket()
            fill(cfg, b, man, staged=False, plain=False)
            b.add(key, body)                             # same key, different bytes
            with wired(b):
                status, out = run(up.stage, cfg)
            assert status == 1, (label, status, out)
            assert "different size or ETag" in out and key in out, (label, out)
            assert b.puts == [], (label, b.puts)         # refuses before writing anything

            b.puts.clear()                               # a dry run reports it and keeps going
            with wired(b):
                status, out = run(up.stage, cfg, dryrun=True)
            assert status is None, (label, status, out)
            assert "WOULD REFUSE" in out and "dry run: nothing was uploaded" in out, out
            assert b.puts == [], b.puts
    finally:
        shutil.rmtree(root)


# ---- release --------------------------------------------------------------------------------------

@case
def test_release_refuses_an_unstaged_build():
    cfg, root, man = make_tree()
    try:
        b = Bucket()                                     # nothing staged at all
        with wired(b):
            status, out = run(up.release, cfg)
        assert status == 1 and b.puts == [], (status, b.puts)
        assert "are not in the bucket; run `stage` first" in out, out
        assert "is not in the bucket; what goes live must be exactly what was staged" in out, out

        b = Bucket()                                     # files staged, manifest copy missing
        fill(cfg, b, man, staged=False)
        with wired(b):
            status, out = run(up.release, cfg)
        assert status == 1 and b.puts == [], (status, b.puts)
        assert "immutable file(s) are not in the bucket" not in out, out
        assert up._staged_manifest_key(cfg)[0] in out, out

        b = Bucket()                                     # a dry run lists both refusals and writes nothing
        with wired(b):
            status, out = run(up.release, cfg, dryrun=True)
        assert status is None and out.count("WOULD REFUSE") == 2, (status, out)
    finally:
        shutil.rmtree(root)


@case
def test_release_uploads_changed_plain_files_and_the_manifest_last():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        fill(cfg, b, man)
        thr = cfg["r2"]["multipart_threshold_bytes"]
        small, large = "coloc_loci.json", "gwas_dcm_bins.json"
        assert len(ASSETS[small]) <= thr < len(ASSETS[large])
        # same size, different bytes: at or under the threshold the ETag is the MD5 and is compared;
        # above it the bucket ETag is a multipart hash from the old `sync`, so only the size is
        b.add(small, b"?" * len(ASSETS[small]), "no-cache", "application/json")
        b.add(large, b"?" * len(ASSETS[large]), "no-cache", "application/json")
        with wired(b):
            status, out = run(up.release, cfg)
        assert status is None, out
        assert b.puts == [small, "manifest.json"], b.puts
        assert b.objects["manifest.json"]["body"] == (cfg.derived / "manifest.json").read_bytes()
        assert b.objects["manifest.json"]["cache_control"] == cfg["r2"]["cache_control"]["mutable"]
        assert b.objects[small]["body"] == ASSETS[small]
        assert b.objects[large]["body"] == b"?" * len(ASSETS[large])    # size rule: left alone

        b = Bucket()                                     # a missing size gets re-uploaded whatever the rule
        fill(cfg, b, man)
        b.add(large, b"?" * (len(ASSETS[large]) + 1), "no-cache", "application/json")
        with wired(b):
            status, out = run(up.release, cfg)
        assert status is None and b.puts == [large, "manifest.json"], (status, b.puts)
    finally:
        shutil.rmtree(root)


@case
def test_pointer_uploads_a_manifest():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        old = root / "rollback.json"
        old.write_text('{"built": "yesterday"}')
        with wired(b):
            status, out = run(up.pointer, cfg, old)
        assert status is None and b.puts == ["manifest.json"], (status, b.puts)
        assert b.objects["manifest.json"]["cache_control"] == "no-cache"
        assert hashlib.sha256(old.read_bytes()).hexdigest() in out, out

        notjson = root / "broken.json"
        notjson.write_text("{oops")
        b.puts.clear()
        with wired(b):
            status, out = run(up.pointer, cfg, notjson)
        assert status == 1 and b.puts == [], (status, b.puts, out)
    finally:
        shutil.rmtree(root)


# ---- prune ----------------------------------------------------------------------------------------

def _retire_keys(cfg: Config, bucket: Bucket) -> tuple[str, str]:
    a, c = "old_eqtl/chr=chr1/data.parquet", "old_eqtl/chr=chr2/data.parquet"
    bucket.add(a, b"a" * 40, "no-cache", "application/vnd.apache.parquet")
    bucket.add(c, b"c" * 50, "no-cache", "application/vnd.apache.parquet")
    keep = cfg.derived / "_old" / a                      # only the first one has a local copy
    keep.parent.mkdir(parents=True, exist_ok=True)
    keep.write_bytes(b"a" * 40)
    return a, c


@case
def test_prune_dry_run_lists_candidates_stale_and_left_alone():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        fill(cfg, b, man)
        b.add("manifest.json", (cfg.derived / "manifest.json").read_bytes(), "no-cache", "application/json")
        a, c = _retire_keys(cfg, b)
        b.add("notes.txt", b"sam's own file", "no-cache", "text/plain")
        stale = "immutable/eqtl.chr1.0000000000000000.qbe"
        b.add(stale, b"an older build")

        with wired(b):
            status, out = run(up.prune, cfg)
        assert status is None and b.deleted == [], (status, b.deleted)
        assert "756" not in out                          # sanity: this is the small tree
        assert "prune candidates: 2 keys, 90 B" in out, out
        assert f"local   {a}" in out and f"NO COPY {c}" in out, out
        assert "left alone" in out and "notes.txt" in out and stale in out.split("left alone")[1], out
        assert "dry run: nothing was deleted" in out, out

        with wired(b):                                   # --stale promotes the orphaned immutable key
            status, out = run(up.prune, cfg, stale=True)
        assert status is None and b.deleted == [], (status, b.deleted)
        assert "prune candidates: 3 keys" in out, out
        assert stale in out.split("left alone")[0], out
        assert "notes.txt" in out.split("left alone")[1], out
    finally:
        shutil.rmtree(root)


@case
def test_prune_yes_refuses_before_release_and_without_a_local_copy():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        fill(cfg, b, man)
        a, c = _retire_keys(cfg, b)
        b.add("manifest.json", b'{"built": "the old build"}', "no-cache", "application/json")
        with wired(b):
            status, out = run(up.prune, cfg, yes=True)
        assert status == 1 and b.deleted == [], (status, b.deleted)
        assert "is not this build" in out and "still reads" in out, out

        b.add("manifest.json", (cfg.derived / "manifest.json").read_bytes(), "no-cache", "application/json")
        with wired(b):                                   # released now, but chr2 has no local copy
            status, out = run(up.prune, cfg, yes=True)
        assert status == 1 and b.deleted == [], (status, b.deleted)
        assert "have no local copy" in out and "rollback stays possible" in out, out

        keep = cfg.derived / "_old" / c                  # fetch it, then the delete goes ahead
        keep.parent.mkdir(parents=True, exist_ok=True)
        keep.write_bytes(b"c" * 50)
        with wired(b):
            status, out = run(up.prune, cfg, yes=True)
        assert status is None, (status, out)
        assert sorted(b.deleted) == sorted([a, c]), b.deleted
        for key in list(man["immutable"]) + ["manifest.json", "coloc_loci.json", "gwas_dcm.json"]:
            assert key in b.objects, f"prune deleted a referenced key: {key}"
    finally:
        shutil.rmtree(root)


# ---- budget ----------------------------------------------------------------------------------------

@case
def test_budget_math_and_the_free_tier():
    cfg, root, man = make_tree()
    try:
        keys = {"old_eqtl/chr=chr1/data.parquet": {"size": 1_000, "etag": "x"},
                "old_eqtl/chr=chr2/data.parquet": {"size": 2_000, "etag": "x"},
                "keep_me.bin": {"size": 500, "etag": "x"}}
        raw = (cfg.derived / "manifest.json").read_bytes()
        current, add = 3_500, sum(e["bytes"] for e in man["immutable"].values()) + len(raw)
        b = quiet(up.budget, cfg, man, keys, overlap_days=3)
        assert (b["current"], b["add"]) == (current, add), b
        assert b["peak"] == current + add and b["retire"] == 3_000 and b["retire_keys"] == 2, b
        assert b["final"] == current + add - 3_000, b
        assert b["staged_key"] == up._staged_manifest_key(cfg)[0]
        assert len(b["new_keys"]) == len(man["immutable"]) and b["over"] is False, b

        before = dt.date.today().day - 1                 # R2 averages each day's storage over 30 days
        after = max(0, 30 - before - 3)
        want = (before * b["current"] + 3 * b["peak"] + after * b["final"]) / 30
        assert abs(b["this_month"] - want) < 1e-6, (b["this_month"], want)
        zero = quiet(up.budget, cfg, man, keys, overlap_days=0)
        assert abs(zero["this_month"] - (before * current + (30 - before) * zero["final"]) / 30) < 1e-6

        cfg.cfg["r2"]["free_tier"]["storage_bytes"] = b["final"]      # `final` at the tier is already too much
        over = quiet(up.budget, cfg, man, keys, overlap_days=3)
        assert over["over"] is True
        bucket = Bucket()
        bucket.objects = {k: {"body": b"x" * v["size"], "cache_control": "", "content_type": ""}
                          for k, v in keys.items()}
        with wired(bucket):
            status, out = run(up.cmd_budget, cfg, None, None, 3)
            assert status == 1, (status, out)
            assert "REFUSE: this deploy does not fit the free tier" in out, out
            status, out = run(up.stage, cfg)             # stage runs budget first and stops on it
            assert status == 1 and bucket.puts == [], (status, bucket.puts)

        cfg.cfg["r2"]["free_tier"]["storage_bytes"] = 10_000_000_000
        with wired(bucket):
            status, out = run(up.cmd_budget, cfg, None, None, 3)
        assert status == 0, (status, out)
    finally:
        shutil.rmtree(root)


# ---- check ------------------------------------------------------------------------------------------

def _checkable() -> tuple[Config, Path, dict, Bucket]:
    cfg, root, man = make_tree()
    b = Bucket()
    fill(cfg, b, man)
    b.add("manifest.json", (cfg.derived / "manifest.json").read_bytes(), "no-cache", "application/json")
    return cfg, root, man, b


@case
def test_check_passes_on_a_correct_bucket():
    cfg, root, man, b = _checkable()
    try:
        with wired(b):
            status, out = run(up.check, cfg)
        assert status == 0, out
        assert "check: OK" in out and "of 3 keys clean" in out, out
        assert "bytes=0-15 -> 206" in out and "bytes=184-199 -> 206" in out, out    # biggest pack is 200 B
        assert "manifest.json -> 200" in out, out

        with wired(b):                                   # --staged also covers the staged manifest copy
            status, out = run(up.check, cfg, staged=True)
        assert status == 0 and "of 4 keys clean" in out, out
        assert "manifest.json -> 200" not in out, out    # the live pointer is not checked before release
    finally:
        shutil.rmtree(root)


@case
def test_check_catches_one_broken_rule_at_a_time():
    key_of = lambda man: sorted(man["immutable"])[0]

    def drop(name):
        def m(path, status, h, body):
            h = {k: v for k, v in h.items() if k != name}
            return status, h, body
        return m

    def swap(name, value, only=None):
        def m(path, status, h, body):
            if name in h and (only is None or only(path, h)):
                h = {**h, name: value}
            return status, h, body
        return m

    rules = [
        ("status", lambda man, k: (lambda p, s, h, b: (500, h, b) if p == k else (s, h, b)), "-> 500"),
        ("length", lambda man, k: swap("content-length", "3"), "content-length 3"),
        ("etag", lambda man, k: swap("etag", '"deadbeef"'), "ETag \"deadbeef\""),
        ("cache-control", lambda man, k: swap("cache-control", "public, max-age=60"), "cache-control 'public, max-age=60'"),
        ("content-range", lambda man, k: swap("content-range", "bytes 0-1/2"), "content-range 'bytes 0-1/2'"),
        ("content-encoding", lambda man, k: (lambda p, s, h, b: (s, {**h, "content-encoding": "gzip"}, b)),
         "there must be none"),
        ("cors", lambda man, k: drop("access-control-allow-origin"), "access-control-allow-origin None"),
    ]
    for label, build, needle in rules:
        cfg, root, man, b = _checkable()
        try:
            b.mutate = build(man, key_of(man))
            with wired(b):
                status, out = run(up.check, cfg)
            assert status == 1, (label, status, out)
            assert needle in out, (label, needle, out)
            assert ("CORS is Sam's" in out) == (label == "cors"), (label, out)
        finally:
            shutil.rmtree(root)

    cfg, root, man, b = _checkable()                     # a missing expose header is CORS too
    try:
        b.mutate = swap("access-control-expose-headers", "ETag")
        with wired(b):
            status, out = run(up.check, cfg)
        assert status == 1 and "is missing ['Content-Range', 'Content-Length', 'Accept-Ranges']" in out, out
        assert "CORS is Sam's" in out, out
    finally:
        shutil.rmtree(root)


@case
def test_check_retired():
    cfg, root, man, b = _checkable()
    try:
        with wired(b):
            status, out = run(up.check, cfg, retired=True)
        assert status == 0 and "no key sits under a replaced prefix" in out, out
        assert "search_index.parquet -> 404" in out, out

        b.add("old_eqtl/chr=chr1/data.parquet", b"still here", "no-cache", "application/vnd.apache.parquet")
        b.add("search_index.parquet", b"also still here", "no-cache", "application/vnd.apache.parquet")
        with wired(b):
            status, out = run(up.check, cfg, retired=True)
        assert status == 1 and "keys still under a replaced prefix" in out, out
        assert "search_index.parquet -> 200, expected 404" in out, out
    finally:
        shutil.rmtree(root)


# ---- --listing ---------------------------------------------------------------------------------------

@case
def test_listing_is_refused_with_any_write():
    cfg, root, man = make_tree()
    try:
        b = Bucket()
        fill(cfg, b, man)
        saved = root / "listing.json"
        saved.write_text(json.dumps({"Source": "test", "Contents": [
            {"Key": k, "Size": len(v["body"]), "ETag": f'"{hashlib.md5(v["body"]).hexdigest()}"'}
            for k, v in b.objects.items()]}))
        for fn, kw in ((up.stage, {}), (up.release, {}), (up.prune, {"yes": True})):
            with wired(b):
                status, out = run(fn, cfg, listing=saved, **kw)
            assert status == 1, (fn.__name__, status, out)
            assert "--listing is a saved snapshot" in out, (fn.__name__, out)
            assert b.puts == [] and b.deleted == [], (fn.__name__, b.puts, b.deleted)

        for fn, kw in ((up.stage, {"dryrun": True}), (up.release, {"dryrun": True}), (up.prune, {})):
            with wired(b):
                status, out = run(fn, cfg, listing=saved, **kw)
            assert status is None, (fn.__name__, status, out)
            assert b.puts == [] and b.deleted == [], (fn.__name__, b.puts, b.deleted)

        keys = b.list(cfg, saved)                        # the saved shape reads back as a listing
        assert keys == b.list(cfg), (keys, b.list(cfg))
    finally:
        shutil.rmtree(root)


def main() -> int:
    for fn in CASES:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"{len(CASES)} cases passed (no network)")
    return 0


def test_all():
    """pytest entry (`uv run pytest pipeline/test_upload.py`)."""
    for fn in CASES:
        fn()


if __name__ == "__main__":
    sys.exit(main())
