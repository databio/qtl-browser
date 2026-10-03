"""The catalog overlap index (SPEC.md section 19, kind 11, `.qbo`): which variant catalogs of a store
hold which sites.

One object per store. It holds the union of the catalogs' sites in canonical order, and for every
union site a bit per catalog. From that a reader gets membership in one chunk read, a catalog's vidx
for a shared site, and exact overlap counts -- which the pointer also records outright, so the
statistics need no object read at all.

The index is derived: rebuild it from the catalogs at any time, and a store without one is complete.
It is not part of any catalog's or experiment's identity. It is invalidated by rebuilding any
catalog it names, which `Store.validate` catches through the `identity_digest` copied into the
pointer.

    uv run python -m pipeline.overlap build --store DIR --id grch38_heart --catalogs A B ...
    uv run python -m pipeline.overlap show  --store DIR --id grch38_heart
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np

from . import catalog, packfmt_v1 as pf, qtlstore as qs

MAGIC = b"QBO1"
KIND_OVERLAP = 11
PAGE_SIZE = 512          # union sites per page; pages never span chromosomes
CHUNK = 65536            # union sites per mask chunk
MAX_CATALOGS = 64        # the mask ceiling; past this, masks must become per-catalog bitmaps
OVERLAPS = "overlaps"    # pointer level
SUSPECT_WINDOW = 50      # bp, for the indel-normalisation heuristic
SUSPECT_EXAMPLES = 20


def mask_width(k: int) -> int:
    """Bytes per membership mask: the smallest of 1, 2, 4, 8 that holds k bits."""
    if not 1 <= k <= MAX_CATALOGS:
        raise ValueError(f"overlap: {k} catalogs, the mask holds 1..{MAX_CATALOGS}")
    return next(w for w in (1, 2, 4, 8) if k <= 8 * w)


_MASK_DTYPE = {1: "<u1", 2: "<u2", 4: "<u4", 8: "<u8"}


# ---- union pages ------------------------------------------------------------------------------
# Section 5's variants page with the attribute columns dropped: the key is all an index needs, and
# af/counts/rsID describe a cohort rather than a site.

def encode_union_page(first: int, pos, ref, alt, level: int = catalog.ZSTD_LEVEL) -> bytes:
    """One page: [12-byte page header][zstd frame of u32 heap_len, u32 deltas, u8 codes, heap]."""
    n = len(pos)
    if not 0 < n <= 65535:
        raise ValueError(f"union page: {n} records")
    pos = np.asarray(pos, dtype=np.int64)
    if n > 1 and np.any(np.diff(pos) < 0):
        raise ValueError("union page: positions decrease")
    codes = np.fromiter((pf.SNP_CODES.get((r, a), 0) for r, a in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = "".join(f"{r}\t{a}\n" for r, a, c in zip(ref, alt, codes.tolist()) if c == 0).encode("ascii")
    delta = pos.copy()
    delta[1:] = np.diff(pos)
    if np.any(delta < 0) or np.any(delta > 0xFFFFFFFF):
        raise ValueError("union page: position delta out of range")
    payload = (struct.pack("<I", len(heap)) + delta.astype("<u4").tobytes() + codes.tobytes() + heap)
    stored = pf.zstd_frame(payload, level)
    return pf._PAGE_HEADER.pack(len(stored), first, n, pf.CODECS["zstd"], 0) + stored


def decode_union_page(buf: bytes, off: int) -> dict:
    """The page at `off` -> {first, n, pos, ref, alt, end}."""
    stored_len, first, n, codec, reserved = pf._PAGE_HEADER.unpack_from(buf, off)
    if n == 0 or reserved:
        raise ValueError(f"union page at {off}: n {n}, reserved {reserved}")
    end = off + pf.PAGE_HEADER_LEN + stored_len
    stored = buf[off + pf.PAGE_HEADER_LEN:end]
    p = stored if codec == pf.CODECS["raw"] else pf.zstd_unframe(stored, None, f"union page at {off}")
    (heap_len,) = struct.unpack_from("<I", p, 0)
    if len(p) != 4 + 5 * n + heap_len:
        raise ValueError(f"union page at {off}: payload length {len(p)}")
    pos = np.cumsum(np.frombuffer(p, "<u4", n, 4), dtype=np.int64)
    code = np.frombuffer(p, "u1", n, 4 + 4 * n)
    if np.any(code > 12):
        raise ValueError(f"union page at {off}: allele code above 12")
    pairs = [pf.SNP_ALLELES.get(c, (None, None)) for c in code.tolist()]
    heap = p[4 + 5 * n:].decode("ascii").split("\n")[:-1] if heap_len else []
    zero = np.flatnonzero(code == 0)
    if len(heap) != zero.size:
        raise ValueError(f"union page at {off}: heap holds {len(heap)} records for {zero.size} code-0 sites")
    for i, rec in zip(zero.tolist(), heap):
        pairs[i] = tuple(rec.split("\t"))
    return {"first": first, "n": n, "pos": pos, "ref": [x[0] for x in pairs], "alt": [x[1] for x in pairs],
            "end": end + pf.pad4(end - off)}


# ---- reading a catalog's sites in canonical order ----------------------------------------------

def catalog_sites(store: qs.Store, chrom: dict) -> list[tuple]:
    """One chromosome of a catalog as (pos, ref, alt, vidx), sorted by (pos, ref, alt).

    The pack stores the cis section first and the trans-only section after it, each sorted; the
    canonical order the index uses interleaves them, so this merges the two runs and keeps each
    site's vidx. That gap between the two orders is exactly why vidx is not a rank (section 19).
    """
    d = catalog.decode_file((store.immutable / chrom["file"]).read_bytes())
    pos, ref, alt = d["pos"].tolist(), d["ref"], d["alt"]
    out = [(pos[i], ref[i], alt[i], i) for i in range(len(pos))]
    out.sort(key=lambda x: (x[0], x[1], x[2]))
    return out


# ---- build ------------------------------------------------------------------------------------

def _pairs(ids: list[str]) -> list[tuple[str, str]]:
    return [(ids[i], ids[j]) for i in range(len(ids)) for j in range(i + 1, len(ids))]


def _conflicts_and_suspects(by_cat: dict[str, list[tuple]], ids: list[str]) -> tuple[int, list, int, list]:
    """Per chromosome: positions two catalogs both hold with no shared allele pair, and indels in
    different catalogs that could be one event written two ways.

    A suspect is never a verdict: confirming one needs the reference to left-align against, and
    section 5 does not require left-aligned indels. The count is a prompt to look."""
    at_pos: dict[int, dict[str, set]] = {}
    for cid, rows in by_cat.items():
        for p, r, a, _ in rows:
            at_pos.setdefault(p, {}).setdefault(cid, set()).add((r, a))
    conflicts, conflict_ex = 0, []
    for p, per in at_pos.items():
        if len(per) < 2:
            continue
        if any(not (per[a] & per[b]) for a, b in _pairs(sorted(per))):
            conflicts += 1
            if len(conflict_ex) < SUSPECT_EXAMPLES:
                conflict_ex.append({"pos": p, **{c: sorted(v) for c, v in per.items()}})
    # indels bucketed by length change, then a sliding window over position
    buckets: dict[int, list[tuple]] = {}
    for cid, rows in by_cat.items():
        for p, r, a, _ in rows:
            if len(r) != len(a):
                buckets.setdefault(len(a) - len(r), []).append((p, r, a, cid))
    suspects, suspect_ex = 0, []
    for rows in buckets.values():
        rows.sort()
        for i, (p, r, a, c) in enumerate(rows):
            for j in range(i + 1, len(rows)):
                q, r2, a2, c2 = rows[j]
                if q - p > SUSPECT_WINDOW:
                    break
                if c2 == c or (p, r, a) == (q, r2, a2):
                    continue
                if (r, a) in at_pos.get(q, {}).get(c2, set()) or (r2, a2) in at_pos.get(p, {}).get(c, set()):
                    continue          # the same key is in both catalogs: nothing to suspect
                suspects += 1
                if len(suspect_ex) < SUSPECT_EXAMPLES:
                    suspect_ex.append({"a": {"catalog": c, "pos": p, "ref": r, "alt": a},
                                       "b": {"catalog": c2, "pos": q, "ref": r2, "alt": a2}})
    return conflicts, conflict_ex, suspects, suspect_ex


def build(store: qs.Store, oid: str, catalog_ids: list[str], *, page_size: int = PAGE_SIZE,
          chunk: int = CHUNK, level: int = catalog.ZSTD_LEVEL, log=lambda *_: None) -> dict:
    """Build the index over `catalog_ids`, write the object and the pointer, return the pointer doc."""
    ids = sorted(set(catalog_ids))
    if len(ids) != len(catalog_ids):
        raise ValueError("overlap: repeated catalog id")
    k = len(ids)
    width = mask_width(k)
    docs = {cid: store.load(qs.CATALOGS, cid) for cid in ids}
    collections = {d["collection_digest"] for d in docs.values()}
    if len(collections) != 1:
        raise ValueError(f"overlap: catalogs span {len(collections)} sequence collections; "
                         "sites on different references are not comparable")
    collection = collections.pop()

    # chromosomes in canonical order: by seq_digest, ASCII (section 5)
    chroms: dict[str, dict] = {}
    for cid, d in docs.items():
        for c in d["chromosomes"]:
            chroms.setdefault(c["seq_digest"], {})[cid] = c
    order = sorted(chroms)

    pages, masks_all, chrom_dir = [], [], []
    prefix = np.zeros((k, 0), dtype=np.int64)
    breaks: list[list[tuple[int, int]]] = [[] for _ in range(k)]
    seen = np.zeros(k, dtype=np.int64)          # sites of each catalog placed so far (the rank)
    delta_now = np.full(k, np.iinfo(np.int64).min, dtype=np.int64)
    masks_flat: list[int] = []
    n_union = 0
    by_chrom: dict[str, dict] = {}
    conflicts, conflict_ex, suspects, suspect_ex = 0, [], 0, []

    for sd in order:
        per = chroms[sd]
        name = next(iter(per.values()))["name"]
        by_cat = {cid: catalog_sites(store, c) for cid, c in per.items()}
        log(f"overlap: {name} ({sd}) " + ", ".join(f"{cid} {len(v):,}" for cid, v in by_cat.items()))
        cf, cfx, sp, spx = _conflicts_and_suspects(by_cat, ids)
        conflicts += cf
        suspects += sp
        conflict_ex += cfx[:max(0, SUSPECT_EXAMPLES - len(conflict_ex))]
        suspect_ex += spx[:max(0, SUSPECT_EXAMPLES - len(suspect_ex))]

        merged: dict[tuple, list] = {}
        for bit, cid in enumerate(ids):
            for p, r, a, vidx in by_cat.get(cid, []):
                e = merged.setdefault((p, r, a), [0, {}])
                e[0] |= 1 << bit
                e[1][bit] = vidx
        keys = sorted(merged)
        first_union_index, n_chrom_union = n_union, len(keys)
        chrom_dir.append((sd, n_chrom_union, first_union_index))
        by_chrom[name] = {"union": n_chrom_union,
                          "sites": {cid: len(by_cat.get(cid, [])) for cid in ids if cid in by_cat},
                          "shared": {f"{a}|{b}": 0 for a, b in _pairs(ids)}}

        for off in range(0, len(keys), page_size):
            part = keys[off:off + page_size]
            pages.append((n_union + off, [x[0] for x in part], [x[1] for x in part], [x[2] for x in part]))

        bit_of = {cid: i for i, cid in enumerate(ids)}
        for j, key in enumerate(keys):
            m, vidxs = merged[key]
            masks_flat.append(m)
            u = n_union + j
            for bit, vidx in vidxs.items():
                d = vidx - seen[bit]
                if d != delta_now[bit]:
                    breaks[bit].append((u, vidx))
                    delta_now[bit] = d
                seen[bit] += 1
            for a, b in _pairs(ids):
                if m >> bit_of[a] & 1 and m >> bit_of[b] & 1:
                    by_chrom[name]["shared"][f"{a}|{b}"] += 1
        n_union += n_chrom_union

    n_chunks = max(1, -(-n_union // chunk))
    masks = np.array(masks_flat, dtype=_MASK_DTYPE[width]) if masks_flat else np.zeros(0, _MASK_DTYPE[width])
    prefix = np.zeros((k, n_chunks), dtype=np.int64)
    for c in range(n_chunks):
        if c:
            part = masks[(c - 1) * chunk:c * chunk]
            for i in range(k):
                prefix[i, c] = prefix[i, c - 1] + int(np.count_nonzero(part >> np.uint64(i) & 1)) \
                    if width == 8 else prefix[i, c - 1] + int(np.count_nonzero(part >> i & 1))

    chunk_bytes = [pf.zstd_frame(masks[c * chunk:(c + 1) * chunk].tobytes(), level) for c in range(n_chunks)]
    page_bytes = [encode_union_page(f, p, r, a, level) for f, p, r, a in pages]

    n_pages = len(page_bytes)
    dir_fixed = struct.pack("<4sIIIIBBHI", MAGIC, n_union, page_size, chunk, n_chunks, width, k, 0, len(chrom_dir))
    dir_chroms = b"".join(sd.encode("ascii") + struct.pack("<II", n, f) for sd, n, f in chrom_dir)
    dir_len = (len(dir_fixed) + len(dir_chroms) + 4 * (n_pages + 1) + 4 * n_pages + 4 * (n_chunks + 1)
               + 4 * k * n_chunks + sum(4 + 8 * len(b) for b in breaks))
    body = qs.HEADER_LEN + dir_len
    page_off, at = [], body
    for b in page_bytes:
        page_off.append(at)
        at += len(b)
    page_off.append(at)
    chunk_off, at2 = [], at
    for b in chunk_bytes:
        chunk_off.append(at2)
        at2 += len(b)
    chunk_off.append(at2)

    directory = (dir_fixed + dir_chroms
                 + np.array(page_off, "<u4").tobytes()
                 + np.array([int(p[1][0]) for p in pages], "<u4").tobytes()
                 + np.array(chunk_off, "<u4").tobytes()
                 + prefix.astype("<u4").tobytes()
                 + b"".join(struct.pack("<I", len(bp)) + np.array(bp, "<u4").reshape(-1).tobytes()
                            for bp in breaks))
    if len(directory) != dir_len:
        raise AssertionError(f"overlap: directory {len(directory)} bytes, computed {dir_len}")
    head = qs.file_header(KIND_OVERLAP, "all", k, page_size, 0, collection)
    obj = head + directory + b"".join(page_bytes) + b"".join(chunk_bytes)
    name = store.put(obj, "qbo")

    totals = {f"{a}|{b}": sum(v["shared"][f"{a}|{b}"] for v in by_chrom.values()) for a, b in _pairs(ids)}
    doc = {
        "id": oid, "object": name, "collection_digest": collection,
        "catalogs": [{"id": cid, "identity_digest": docs[cid]["identity_digest"],
                      "n_sites": docs[cid]["n_sites"]} for cid in ids],
        "mask_width": width, "chunk": chunk, "page_size": page_size, "union": n_union,
        "pairs": [{"a": a, "b": b, "shared": totals[f"{a}|{b}"],
                   "jaccard": round(totals[f"{a}|{b}"] / u, 6) if (u := docs[a]["n_sites"] + docs[b]["n_sites"]
                                                                  - totals[f"{a}|{b}"]) else 0.0,
                   "a_only": docs[a]["n_sites"] - totals[f"{a}|{b}"],
                   "b_only": docs[b]["n_sites"] - totals[f"{a}|{b}"]} for a, b in _pairs(ids)],
        "by_chrom": by_chrom,
        "conflicts": {"shared_positions_no_shared_allele_pair": conflicts, "examples": conflict_ex},
        "normalisation_suspects": {"count": suspects, "window_bp": SUSPECT_WINDOW, "examples": suspect_ex},
    }
    store.write_pointer(OVERLAPS, oid, doc)
    return doc


# ---- read -------------------------------------------------------------------------------------

def decode(buf: bytes) -> dict:
    """Header and directory; pages and mask chunks stay in `buf` and are decoded on demand."""
    h = qs.parse_file_header(buf)
    if h["kind"] != KIND_OVERLAP:
        raise ValueError(f"overlap index: kind {h['kind']}")
    magic, n_union, page_size, chunk, n_chunks, width, k, res, n_chrom = struct.unpack_from("<4sIIIIBBHI", buf, qs.HEADER_LEN)
    if magic != MAGIC or res:
        raise ValueError(f"overlap index: directory magic {magic!r}")
    if k != h["count"] or width != mask_width(k) or page_size != h["page_size"]:
        raise ValueError("overlap index: directory disagrees with the header")
    at = qs.HEADER_LEN + struct.calcsize("<4sIIIIBBHI")
    chroms = []
    for _ in range(n_chrom):
        sd = buf[at:at + 32].decode("ascii")
        n, first = struct.unpack_from("<II", buf, at + 32)
        chroms.append({"seq_digest": sd, "n_union": n, "first_union_index": first})
        at += 40
    n_pages = sum(-(-c["n_union"] // page_size) for c in chroms)
    page_off = np.frombuffer(buf, "<u4", n_pages + 1, at); at += 4 * (n_pages + 1)
    page_first = np.frombuffer(buf, "<u4", n_pages, at); at += 4 * n_pages
    chunk_off = np.frombuffer(buf, "<u4", n_chunks + 1, at); at += 4 * (n_chunks + 1)
    prefix = np.frombuffer(buf, "<u4", k * n_chunks, at).reshape(k, n_chunks); at += 4 * k * n_chunks
    breaks = []
    for _ in range(k):
        (nb,) = struct.unpack_from("<I", buf, at); at += 4
        breaks.append(np.frombuffer(buf, "<u4", 2 * nb, at).reshape(-1, 2) if nb else np.zeros((0, 2), "<u4"))
        at += 8 * nb
    if at != int(page_off[0]):
        raise ValueError(f"overlap index: directory ends at {at}, first page at {int(page_off[0])}")
    return {"header": h, "buf": buf, "n_union": int(n_union), "page_size": int(page_size), "chunk": int(chunk),
            "n_chunks": int(n_chunks), "mask_width": int(width), "k": int(k), "chroms": chroms,
            "page_off": page_off, "page_first_position": page_first, "chunk_off": chunk_off,
            "prefix": prefix, "breaks": breaks}


def masks_of_chunk(d: dict, c: int) -> np.ndarray:
    lo, hi = int(d["chunk_off"][c]), int(d["chunk_off"][c + 1])
    raw = pf.zstd_unframe(d["buf"][lo:hi], None, f"mask chunk {c}")
    return np.frombuffer(raw, _MASK_DTYPE[d["mask_width"]])


def membership(d: dict, u: int) -> int:
    """The mask at union index `u`: bit i set when catalog i holds that site."""
    if not 0 <= u < d["n_union"]:
        raise IndexError(f"union index {u} of {d['n_union']}")
    return int(masks_of_chunk(d, u // d["chunk"])[u % d["chunk"]])


def rank(d: dict, i: int, u: int) -> int:
    """How many of catalog i's sites lie strictly before union index `u`."""
    c, off = u // d["chunk"], u % d["chunk"]
    part = masks_of_chunk(d, c)[:off]
    return int(d["prefix"][i, c]) + int(np.count_nonzero(np.right_shift(part, i) & 1))


def vidx_in(d: dict, i: int, u: int) -> int | None:
    """Catalog i's vidx for the site at union index `u`, or None when it does not hold it.

    Not a rank: a catalog numbers each chromosome's cis sites before its trans-only sites, and union
    order interleaves the two, so `rank` is the site's index among catalog i's sites in canonical
    order. The breakpoints carry the correction (SPEC section 19)."""
    if not membership(d, u) >> i & 1:
        return None
    bp = d["breaks"][i]
    j = int(np.searchsorted(bp[:, 0], u, side="right")) - 1
    if j < 0:
        raise ValueError(f"overlap index: no breakpoint at or below union index {u} for catalog {i}")
    u0, v0 = int(bp[j, 0]), int(bp[j, 1])
    return v0 + rank(d, i, u) - rank(d, i, u0)


def union_sites(d: dict, first_page: int = 0, n_pages: int | None = None):
    """Decode pages, yielding (union_index, pos, ref, alt)."""
    last = len(d["page_first_position"]) if n_pages is None else first_page + n_pages
    for p in range(first_page, last):
        pg = decode_union_page(d["buf"], int(d["page_off"][p]))
        for j in range(pg["n"]):
            yield pg["first"] + j, int(pg["pos"][j]), pg["ref"][j], pg["alt"][j]


def find(d: dict, seq_digest: str, pos: int, ref: str, alt: str) -> int | None:
    """Union index of a site, or None. Two reads in a browser: the directory, then one page."""
    c = next((x for x in d["chroms"] if x["seq_digest"] == seq_digest), None)
    if c is None or not c["n_union"]:
        return None
    base = sum(-(-x["n_union"] // d["page_size"]) for x in d["chroms"]
               if x["first_union_index"] < c["first_union_index"])
    n_pages = -(-c["n_union"] // d["page_size"])
    firsts = d["page_first_position"][base:base + n_pages]
    p = base + max(0, int(np.searchsorted(firsts, pos, side="right")) - 1)
    for u, q, r, a in union_sites(d, p, 1):
        if (q, r, a) == (pos, ref, alt):
            return u
    return None


# ---- CLI --------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--store", required=True)
    b.add_argument("--id", required=True)
    b.add_argument("--catalogs", nargs="+", help="catalog ids; default every catalog in the store")
    s = sub.add_parser("show")
    s.add_argument("--store", required=True)
    s.add_argument("--id", required=True)
    a = ap.parse_args(argv)
    store = qs.Store(Path(a.store))
    if a.cmd == "build":
        ids = a.catalogs or sorted(p.stem for p in (store.root / qs.CATALOGS).glob("*.json"))
        # object, then pointer, and nothing else: store.json is written last by whoever runs the
        # build (store.sbatch's validate step), and it lists whatever pointers exist then
        doc = build(store, a.id, ids, log=lambda m: print(m, file=sys.stderr))
    else:
        doc = store.load(OVERLAPS, a.id)
    print(json.dumps({k: v for k, v in doc.items() if k != "by_chrom"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
