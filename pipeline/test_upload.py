"""Tests for the uploader's planning (pipeline/upload.py): what it would send, in what order, and
what it would leave alone. No bucket is touched; the listing is a dict.

    uv run python -m pipeline.test_upload
"""
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from pathlib import Path

from . import qtlstore as qs, upload as up

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


def fake_store(root: Path) -> qs.Store:
    """A store-shaped directory: two referenced objects, one leftover nobody names, two pointer
    levels listed in store.json (experiments as the v2 map), and store.json."""
    st = qs.Store(root)
    a = st.put(b"object A bytes", "qbe")
    b = st.put(b"object B bytes", "json")           # a collection document
    st.put(b"leftover from an earlier build", "qbt")
    (root / "variant_catalogs").mkdir(exist_ok=True)
    (root / "experiments").mkdir(exist_ok=True)
    (root / "variant_catalogs" / "cat1.json").write_text(json.dumps({"id": "cat1"}))
    (root / "experiments" / "e1.json").write_text(json.dumps({"id": "e1", "files": {"chr1": a}, "collection": b}))
    (root / "store.json").write_text(json.dumps({"name": "t", "format_version": 2, "refget": [],
                                                 "variant_catalogs": ["cat1"], "annotations": [],
                                                 "experiments": {"e1": {"kind": "qtl_cis", "collection": b}},
                                                 "overlaps": []}))
    return st


def md5(b: bytes) -> str:
    return hashlib.md5(b).hexdigest()


@case
def local_keys_are_the_referenced_objects_and_the_listed_pointers_in_release_order():
    with tempfile.TemporaryDirectory() as d:
        st = fake_store(Path(d))
        local = up.local_keys(st)
        assert len(local.objects) == 2 and all(k.startswith("immutable/") for k in local.objects)
        assert not any("qbt" in k for k in local.objects), "a leftover nobody names is not uploaded"
        assert list(local.pointers) == ["variant_catalogs/cat1.json", "experiments/e1.json"], list(local.pointers)
        assert list(local.all())[-1] == "store.json"


@case
def plan_skips_objects_the_bucket_holds_and_sends_changed_pointers():
    with tempfile.TemporaryDirectory() as d:
        st = fake_store(Path(d))
        local = up.local_keys(st)
        a_key = next(k for k in local.objects if k.endswith(".qbe"))
        b_key = next(k for k in local.objects if k.endswith(".json"))
        listing = {
            a_key: {"size": len(b"object A bytes"), "etag": md5(b"object A bytes")},     # present, same bytes
            b_key: {"size": 3, "etag": "deadbeef"},                                     # present, wrong size
            "experiments/e1.json": {"size": 1, "etag": md5(b"old pointer")},            # changed
            "variant_catalogs/cat1.json": {"size": 1, "etag": md5(json.dumps({"id": "cat1"}).encode())},
            "store.json": {"size": 1, "etag": "0" * 32},
            "immutable/OLDOLDOLDOLDOLDOLDOLDOLDOLDOLDOLD.qbe": {"size": 99, "etag": "x"},  # v1 leftover
            "experiments/gtex.json": {"size": 5, "etag": "y"},
        }
        p = up.plan(local, listing)
        assert p.skip_objects == [a_key] and p.put_objects == [b_key], (p.skip_objects, p.put_objects)
        assert p.put_pointers == ["experiments/e1.json"] and p.skip_pointers == ["variant_catalogs/cat1.json"]
        assert p.put_store_json
        assert p.stale == ["experiments/gtex.json", "immutable/OLDOLDOLDOLDOLDOLDOLDOLDOLDOLDOLD.qbe"]
        assert p.bytes_to_put == len(b"object B bytes") + (Path(d) / "experiments" / "e1.json").stat().st_size \
            + (Path(d) / "store.json").stat().st_size
        assert "stale in bucket" in up.describe(p, listing)


@case
def an_object_in_the_bucket_with_a_multipart_etag_is_trusted_by_size():
    with tempfile.TemporaryDirectory() as d:
        st = fake_store(Path(d))
        local = up.local_keys(st)
        a_key = next(k for k in local.objects if k.endswith(".qbe"))
        p = up.plan(local, {a_key: {"size": len(b"object A bytes"), "etag": "abc-3"}})
        assert a_key in p.skip_objects


@case
def headers_follow_the_key():
    assert up.cache_control("immutable/x.qbe") == up.IMMUTABLE_CACHE
    assert up.cache_control("immutable/x.json") == up.IMMUTABLE_CACHE      # a collection document is content-addressed
    assert up.cache_control("experiments/e.json") == up.POINTER_CACHE and up.cache_control("store.json") == up.POINTER_CACHE
    assert up.content_type("immutable/x.json") == "application/json"       # a collection document
    assert up.content_type("experiments/e.json") == "application/json"
    assert up.content_type("immutable/x.arrow.zst") == "application/octet-stream"


@case
def a_pointer_naming_a_missing_object_stops_before_any_listing():
    with tempfile.TemporaryDirectory() as d:
        st = fake_store(Path(d))
        (Path(d) / "experiments" / "e2.json").write_text(json.dumps({"id": "e2", "files": {"chr1": "A" * 32 + ".qbe"}}))
        try:
            up.local_keys(st)
            raise AssertionError("a dangling pointer must stop the plan")
        except SystemExit as e:
            assert "does not hold" in str(e)


def main() -> int:
    bad = 0
    for fn in CASES:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(CASES) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
