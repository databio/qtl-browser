"""Which codec wrote each object in a store, read from the bytes.

    uv run python -m pipeline.which_codec --store data/derived/store-local-v2

A v2 store must contain no v1-encoded statistics. The claim is not worth trusting from the build
code, so this reads it back out: the cis block magic (`QGB0` v1, `QGB2` v2) and the trans frame
magic (`QTT2` v1, `QTT3` v2) say which encoder ran, per block and per frame rather than per object.

Three object kinds carry no statistics and are **codec-independent by design**, so they are reported
as such rather than as a gap: the variant catalog, the rsID index and the gene lookup hold positions,
alleles and names; the search index is Arrow; and the hits records are fixed `f32`/`u32` fields with
no quantization scales, which is why plan G-section says the kind-0 `value` needed only a
redefinition and not a re-encoding.

GWAS blocks carry **no magic at all** in either version -- the payload opens with `u32 n` -- so the
only thing distinguishing a v2 GWAS block from a v1 one is the containing file header's version
byte, which is why that byte is a parameter threaded from the codec (plan D5) rather than a module
constant. A `.qbg` is therefore reported from its header and not from its bytes: it is the one kind
where this tool cannot check the build's claim independently, and the row-length arithmetic below is
the nearest thing to a second opinion.

`Store.validate` now checks every header against store.json's `format_version`, so the mixed store
this tool was written to find is also a validation failure. It is still worth reading the magics
back: validate compares headers to a declaration, and only the magic says what actually encoded a
block.
"""
from __future__ import annotations

import argparse
import json
import re
import struct
import sys
from pathlib import Path

from . import packfmt_v1 as v1
from . import packfmt_v2 as v2
from . import qtlstore as qs

PAT = re.compile(r"^[A-Za-z0-9_-]{20,}\.[a-z0-9.]{3,12}$")
CODEC_FREE = {"qbv": "variant catalog (positions and alleles)", "qbx": "variant index",
              "qbr": "rsID index", "qgl": "gene lookup", "arrow.zst": "search index / annotation (Arrow)",
              "qbh": "hits (fixed f32/u32 fields, no scales)", "qgi": "GWAS index (offsets and the n table)"}


def names(o, out):
    if isinstance(o, str):
        if PAT.match(o):
            out.add(o)
    elif isinstance(o, dict):
        for v in o.values():
            names(v, out)
    elif isinstance(o, list):
        for v in o:
            names(v, out)


def first_frame(buf: bytes) -> bytes:
    """The payload of the first zstd frame after the 64-byte file header, ignoring the frames after
    it. Trans and GWAS objects are many frames concatenated, and the decoders refuse trailing data
    by design, so a probe that wants one frame has to stream it."""
    import zstandard
    return zstandard.ZstdDecompressor().decompressobj().decompress(buf[qs.HEADER_LEN:])


def first_frame_magic(buf: bytes) -> bytes:
    try:
        return first_frame(buf)[:4]
    except Exception:                 # noqa: BLE001 -- an unreadable frame has no magic to report
        return b"????"


def gwas_row_arithmetic(buf: bytes) -> list[int]:
    """Which codecs block 0's *length* is consistent with, read without trusting the file header.

    A GWAS block has no magic, but both versions open their payload with `u32 n, u32 heap_len` at
    the same two offsets and then lay out fixed-width rows: 21 bytes in v1, 18 in v2, after headers
    of 8 and 40 bytes. So exactly one of `HL + RB*n + heap_len == len(payload)` holds unless the
    numbers collide, and that is an answer the build code did not get to supply.

    Block 0 is read with a streaming decompressor rather than from the index, so this needs only the
    one object. Returns the versions that fit, usually a single element."""
    try:
        payload = first_frame(buf)
        n, heap_len = struct.unpack_from("<II", payload, 0)
    except Exception:                 # noqa: BLE001 -- an unreadable block fits no version
        return []
    return [c.FORMAT_VERSION for c in (v1, v2)
            if len(payload) == c.GWAS_HEADER_LEN + c.GWAS_ROW_BYTES * n + heap_len]


def scan(store: Path) -> dict:
    out: dict = {"v1": [], "v2": [], "codec_free": [], "unknown": [], "file_header_versions": {},
                 "declared": json.loads((store / "store.json").read_text()).get("format_version")}
    ns: set[str] = set()
    for lvl in qs.POINTER_DIRS:
        for p in (store / lvl).glob("*.json"):
            names(json.loads(p.read_text()), ns)
    for n in sorted(ns):
        path = store / "immutable" / n
        if not path.exists():
            continue
        ext = n.split(".", 1)[1]
        buf = path.read_bytes()
        if ext in ("qbe", "qbt", "qbg"):
            h = qs.parse_file_header(buf)
            out["file_header_versions"].setdefault(h["version"], []).append(n)
        if ext == "qbe":                       # one magic per block, not per object
            got = set()
            off = qs.HEADER_LEN
            while off < len(buf):
                magic = buf[off:off + 4]
                got.add(magic.decode("ascii", "replace"))
                blk_len = int.from_bytes(buf[off + 4:off + 8], "little")
                if blk_len <= 0:
                    break
                off += blk_len
            for m in got:
                key = "v2" if m == v2.MAGIC_BLOCK.decode() else "v1" if m == v1.MAGIC_BLOCK.decode() else "unknown"
                out[key].append(f"{n} block magic {m}")
        elif ext == "qbt":
            # One zstd frame per phenotype, concatenated, so the magic has to be read out of the
            # first frame's payload: decoding `buf[HEADER_LEN:]` whole fails on the trailing frames
            # and says nothing about the codec.
            m = first_frame_magic(buf)
            key = "v2" if m == v2.MAGIC_TRANS else "v1" if m == v1.MAGIC_TRANS else "unknown"
            out[key].append(f"{n} trans frame magic {m.decode('ascii', 'replace')}")
        elif ext == "qbg":
            fits = gwas_row_arithmetic(buf)
            claimed = h["version"]
            if fits == [claimed]:
                out[f"v{claimed}"].append(f"{n} GWAS: header version {claimed}, "
                                          f"block 0 fits {qs.codec_for(claimed).GWAS_ROW_BYTES} B/row and nothing else")
            else:
                out["unknown"].append(f"{n} GWAS: header version {claimed}, block 0 fits {fits or 'neither'} "
                                      f"row width")
        else:
            out["codec_free"].append(f"{n} {CODEC_FREE.get(ext, ext)}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--store", required=True, type=Path)
    args = ap.parse_args(argv)
    r = scan(args.store)
    want = r["declared"]
    print(f"store {args.store}, store.json format_version {want}")
    for k, label in (("v2", "v2-encoded"), ("v1", "v1-encoded"), ("unknown", "UNDETERMINED")):
        if r[k]:
            print(f"  {label}:")
            for s in r[k]:
                print(f"    {s}")
    print(f"  codec-independent by design: {len(r['codec_free'])} objects")
    print(f"  file header versions: { {k: len(v) for k, v in r['file_header_versions'].items()} }")
    # the store's own declaration is what everything is judged against, not the directory's name
    wrong = [s for v in (1, 2) if v != want for s in r[f"v{v}"]] + r["unknown"]
    stray = {v: len(f) for v, f in r["file_header_versions"].items() if v != want}
    if stray:
        print(f"  HEADERS DISAGREEING WITH store.json: {stray}")
    print(f"\n{'CLEAN' if not (wrong or stray) else f'{len(wrong)} objects not encoded v{want}'}")
    return 1 if (wrong or stray) else 0


if __name__ == "__main__":
    sys.exit(main())
