"""Step 5: the variant catalog builder -- one variant catalog per (genome, site list), shared by experiments.

Reads a contract `sites` table (CONTRACT.md), the refget collection and a chromosome list; writes

    variant_catalogs/<id>.json               mutable pointer: identity, collection, chromosome table, objects
      <digest>.qbv                   one per chromosome: v1 header + variant pages (vidx order)
      <digest>.qbx                   variant index: per chromosome page offsets and first positions,
                                     then every rsID block's first rs_number
      <digest>.qbr                   rsID index: (rs_number, vidx, chromosome ordinal) records

Ordering
--------
Per chromosome the file holds the **cis** section (`in_cis`) then the **trans-only** section, each
sorted by (pos, ref, alt) byte-wise, exactly as v0 split them; vidx is the row number in that order.
The variant catalog sorts rather than trusting the input order, so two adapters that emit the same sites in a
different order still get the same variant catalog. The identity digest (`qtlstore.catalog_identity`) does
not follow file order at all: it runs over the sites in canonical order (chromosomes by seq_digest,
then pos, ref, alt), so the cis/trans-only split and the chromosome table's order stay out of it.

Chromosome ordinals (the rsID index) are 1-based positions in **this variant catalog's chromosome table**,
never a fixed chr1..chrX list: a variant catalog on another genome or a subset gets its own ordinals.

Variant page, v1
----------------
The v0 page layout (v0 SPEC section 4; SPEC.md section 5) with the allele meaning fixed to ref/alt:

    u32 heap_len | u32 pos delta[n] | u32 rs_number[n] | u16 af[n] | u16 ma_samples[n]
    | u16 ma_count[n] | u8 allele[n] | u8 flags[n] | heap ("ref\\talt\\n" per allele code 0)

`allele` codes the SNP pair (ref, alt) with v0's 12 codes; `af` is the **ALT** frequency coded
rint(af * 65534) (65535 null), unchanged from v0. Flags: bit 0 `alt_is_minor` (af < 0.5; 0 when af is
null), bits 1-2 the rsID match code when the variant catalog declares `match` (0 none, 1 exact, 2 position;
3 is reserved and rejected), bits 3-7 zero. There is no "alleles not reported" bit: a contract site always has both alleles.
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from . import packfmt_v1 as pf
from . import qtlstore as qs

KIND_VARIANTS, KIND_VARIANT_INDEX, KIND_RSID = qs.KIND_VARIANTS, qs.KIND_VARIANT_INDEX, qs.KIND_RSID
EXT_VARIANTS, EXT_VIDX, EXT_RSID = "qbv", "qbx", "qbr"
PAGE_SIZE = 512
RSID_BLOCK_RECORDS = 4096
ZSTD_LEVEL = 19
ALL = qs.ALL                       # header chromosome of objects spanning a whole variant catalog; their seq_digest is the collection
FLAG_ALT_IS_MINOR, FLAG_MATCH_SHIFT, FLAG_MATCH_MASK, FLAG_RESERVED = 0x01, 1, 0x06, 0xF8
MATCH_MAX = 2                      # match codes 0 none, 1 exact, 2 position; 3 is reserved
MAGIC_VIDX = b"QVX1"
_VX = struct.Struct("<4sIIIII")     # magic, n_chrom, page_size, rsid_block_records, rsid_n, rsid_blocks
RSID_DTYPE = np.dtype([("rs_number", "<u4"), ("vidx", "<u4"), ("ordinal", "<u2"), ("pad", "<u2")])
REQUIRED = ("chr", "pos", "ref", "alt", "af", "ma_samples", "in_cis")
# attributes the variant catalog stores, in the order `attributes` lists them. The rsID text is not one: the
# variant catalog stores `rs_number`, and a sites table's `rsid` must be exactly "rs<rs_number>" (checked).
KNOWN_ATTRIBUTES = ("af", "ma_samples", "ma_count", "rs_number", "match")


# ---- pages ------------------------------------------------------------------------------------
def encode_page(first_vidx: int, pos, rs_number, af, ma_samples, ma_count, ref, alt, match_code,
                codec: str = "zstd", level: int = ZSTD_LEVEL) -> bytes:
    n = len(pos)
    if not 1 <= n <= pf.MAX_PAGE_RECORDS:
        raise ValueError(f"page: {n} records")
    pos = np.asarray(pos, dtype=np.int64)
    if pos.min() < 1 or np.any(np.diff(pos) < 0):
        raise ValueError("page: positions must be >= 1 and non-decreasing")
    deltas = np.empty(n, dtype="<u4")
    deltas[0], deltas[1:] = pos[0], np.diff(pos)
    rs = np.asarray(rs_number, dtype=np.int64)
    if rs.max() > pf.U32_MAX:
        raise ValueError("page: rs_number does not fit u32")
    rs = np.where(rs < 0, 0, rs).astype("<u4")
    a = np.asarray(af, dtype=np.float64)
    ok = ~np.isnan(a)
    if np.any((a[ok] < 0) | (a[ok] > 1)):
        raise ValueError("page: af outside [0, 1]")
    afq = np.full(n, pf.AF_NULL, dtype="<u2")
    afq[ok] = np.rint(a[ok] * pf.AF_MAXQ).astype("<u2")
    counts = []
    for x, nm in ((ma_samples, "ma_samples"), (ma_count, "ma_count")):
        c = np.asarray(x, dtype=np.int64)
        if c.max(initial=-1) >= pf.COUNT_NULL:
            raise ValueError(f"page: {nm} above 65534")
        counts.append(np.where(c < 0, pf.COUNT_NULL, c).astype("<u2"))
    codes = np.fromiter((pf.SNP_CODES.get((r, t), 0) for r, t in zip(ref, alt)), dtype=np.uint8, count=n)
    heap = []
    for i in np.flatnonzero(codes == 0):
        r, t = ref[i], alt[i]
        if not (r and t and r.isascii() and t.isascii()) or "\t" in r + t or "\n" in r + t:
            raise ValueError(f"page: allele pair {r!r}/{t!r} is not two non-empty ASCII strings")
        heap.append(f"{r}\t{t}\n")
    heap = "".join(heap).encode("ascii")
    flags = np.where(ok & (a < 0.5), FLAG_ALT_IS_MINOR, 0).astype(np.uint8)
    m = np.asarray(match_code, dtype=np.int64)
    if m.min(initial=0) < 0 or m.max(initial=0) > MATCH_MAX:
        raise ValueError(f"page: match code outside 0..{MATCH_MAX}")
    flags |= (m.astype(np.uint8) << FLAG_MATCH_SHIFT)
    payload = b"".join((struct.pack("<I", len(heap)), deltas.tobytes(), rs.tobytes(), afq.tobytes(),
                        counts[0].tobytes(), counts[1].tobytes(), codes.tobytes(), flags.tobytes(), heap))
    stored = payload if codec == "raw" else pf.zstd_frame(payload, level)
    page = pf._PAGE_HEADER.pack(len(stored), first_vidx, n, pf.CODECS[codec], 0) + stored
    return page + bytes(pf.pad4(len(page)))


def decode_file(buf: bytes) -> dict:
    """Whole v1 variants file -> header plus per-vidx arrays: pos, ref, alt, af_code, af, rs_number (0 =
    none), ma_samples, ma_count (65535 = null), flags, and page offsets/first positions."""
    h = qs.parse_file_header(buf)
    if h["kind"] != KIND_VARIANTS:
        raise ValueError(f"variants file: kind {h['kind']}")
    cols = {k: [] for k in ("pos", "rs", "af", "ms", "mc", "flags")}
    ref, alt, offs, firsts = [], [], [], []
    off, nxt = qs.HEADER_LEN, 0
    while off < len(buf):
        stored_len, first, n, codec, reserved = pf._PAGE_HEADER.unpack_from(buf, off)
        if first != nxt or n == 0 or reserved:
            raise ValueError(f"variants file: bad page header at byte {off}")
        end = off + pf.PAGE_HEADER_LEN + stored_len
        stored = buf[off + pf.PAGE_HEADER_LEN:end]
        p = stored if codec == pf.CODECS["raw"] else pf.zstd_unframe(stored, None, f"page at {off}")
        (heap_len,) = struct.unpack_from("<I", p, 0)
        if len(p) != 4 + 16 * n + heap_len:
            raise ValueError(f"variants file: page at {off} payload length")
        pos = np.cumsum(np.frombuffer(p, "<u4", n, 4), dtype=np.int64)
        code = np.frombuffer(p, "u1", n, 4 + 14 * n)
        flags = np.frombuffer(p, "u1", n, 4 + 15 * n)
        if np.any(flags & FLAG_RESERVED) or np.any(code > 12) or \
                np.any((flags & FLAG_MATCH_MASK) >> FLAG_MATCH_SHIFT > MATCH_MAX):
            raise ValueError(f"variants file: page at {off} reserved flag bits, match code 3 or allele code")
        pairs = [pf.SNP_ALLELES.get(c, (None, None)) for c in code.tolist()]
        heap = p[4 + 16 * n:].decode("ascii").split("\n")[:-1] if heap_len else []
        zero = np.flatnonzero(code == 0)
        if len(heap) != zero.size:
            raise ValueError(f"variants file: page at {off} heap holds {len(heap)} records for {zero.size}")
        for i, rec in zip(zero.tolist(), heap):
            pairs[i] = tuple(rec.split("\t"))
        cols["pos"].append(pos)
        for k, o in (("rs", 4 * n), ("af", 8 * n), ("ms", 10 * n), ("mc", 12 * n)):
            dt = "<u4" if k == "rs" else "<u2"
            cols[k].append(np.frombuffer(p, dt, n, 4 + o))
        cols["flags"].append(flags)
        ref += [x[0] for x in pairs]
        alt += [x[1] for x in pairs]
        offs.append(off)
        firsts.append(int(pos[0]))
        nxt = first + n
        off = end + pf.pad4(end - off)
    if nxt != h["count"]:
        raise ValueError(f"variants file: {nxt} records, header count {h['count']}")
    c = {k: (np.concatenate(v) if v else np.zeros(0, dtype=np.int64)) for k, v in cols.items()}
    af = np.where(c["af"] == pf.AF_NULL, np.nan, c["af"] / pf.AF_MAXQ)
    return {"header": h, "pos": c["pos"], "ref": ref, "alt": alt, "af_code": c["af"], "af": af, "rs_number": c["rs"],
            "ma_samples": c["ms"], "ma_count": c["mc"], "flags": c["flags"],
            "page_off": np.array(offs + [len(buf)], dtype=np.int64), "page_first_position": np.array(firsts, dtype=np.int64)}


def read_sites(store: qs.Store, chrom: dict):
    """The `sites` reader `Store.validate` wants: (pos, ref, alt) in file order, from the pack bytes."""
    d = decode_file((store.immutable / chrom["file"]).read_bytes())
    return d["pos"], d["ref"], d["alt"]


# ---- variant index and rsID index -------------------------------------------------------------
def encode_vidx(chroms: list[dict], page_size: int, rsid_first: np.ndarray, rsid_n: int, collection: str,
                level: int = ZSTD_LEVEL) -> bytes:
    """Kind 9: header (chromosome `all`, the collection digest), one zstd frame: `QVX1` header, then per
    chromosome **in variant catalog table order** n_cis, n_trans, n_pages_cis, n_pages_trans, page offsets
    (n_pages + 1, the last the file size), page first positions; then the rsID block samples."""
    n_blocks = -(-rsid_n // RSID_BLOCK_RECORDS)
    parts = [_VX.pack(MAGIC_VIDX, len(chroms), page_size, RSID_BLOCK_RECORDS, rsid_n, n_blocks)]
    for c in chroms:
        pc_, pt_ = -(-c["n_cis"] // page_size), -(-c["n_trans"] // page_size)
        if len(c["page_off"]) != pc_ + pt_ + 1:
            raise ValueError(f"vidx {c['name']}: {len(c['page_off']) - 1} pages for {c['n_cis']} + {c['n_trans']}")
        parts += [struct.pack("<IIII", c["n_cis"], c["n_trans"], pc_, pt_),
                  np.asarray(c["page_off"], dtype="<u4").tobytes(),
                  np.asarray(c["page_first_position"], dtype="<u4").tobytes()]
    parts.append(np.asarray(rsid_first, dtype="<u4").tobytes())
    return qs.file_header(KIND_VARIANT_INDEX, ALL, len(chroms), page_size, 0, collection) + \
        pf.zstd_frame(b"".join(parts), level)


def decode_vidx(buf: bytes, names: list[str]) -> dict:
    h = qs.parse_file_header(buf)
    p = pf.zstd_unframe(buf[qs.HEADER_LEN:], None, "variant index")
    magic, n_chrom, page_size, rbr, rsid_n, n_blocks = _VX.unpack_from(p, 0)
    if magic != MAGIC_VIDX or n_chrom != len(names) or n_chrom != h["count"]:
        raise ValueError("variant index: bad magic or chromosome count")
    off, out = _VX.size, {}
    for name in names:
        n_cis, n_tr, pc_, pt_ = struct.unpack_from("<IIII", p, off)
        off += 16
        po = np.frombuffer(p, "<u4", pc_ + pt_ + 1, off).astype(np.int64)
        off += 4 * (pc_ + pt_ + 1)
        fp = np.frombuffer(p, "<u4", pc_ + pt_, off).astype(np.int64)
        off += 4 * (pc_ + pt_)
        out[name] = {"n_cis": n_cis, "n_trans": n_tr, "page_off": po, "page_first_position": fp}
    rf = np.frombuffer(p, "<u4", n_blocks, off).astype(np.int64)
    if off + 4 * n_blocks != len(p):
        raise ValueError("variant index: trailing bytes")
    return {"page_size": page_size, "rsid_block_records": rbr, "rsid_n": rsid_n, "rsid_first": rf, "chroms": out}


def encode_rsid(rs_number, vidx, ordinal, collection: str) -> tuple[bytes, np.ndarray]:
    """Kind 8: 12-byte records (rs_number u32, vidx u32, ordinal u16, 0 u16) sorted by (rs_number,
    ordinal, vidx). rs_number may repeat -- one rsID names every allele pair dbSNP puts under it -- so a
    lookup returns a run, not a single record. Returns the bytes and each block's first rs_number."""
    rec = np.zeros(len(rs_number), dtype=RSID_DTYPE)
    rec["rs_number"], rec["vidx"], rec["ordinal"] = rs_number, vidx, ordinal
    rec = rec[np.lexsort((rec["vidx"], rec["ordinal"], rec["rs_number"]))]
    return (qs.file_header(KIND_RSID, ALL, len(rec), RSID_BLOCK_RECORDS, 0, collection) + rec.tobytes(),
            rec["rs_number"][::RSID_BLOCK_RECORDS].astype(np.int64))


def decode_rsid(buf: bytes) -> np.ndarray:
    qs.parse_file_header(buf)
    return np.frombuffer(buf, RSID_DTYPE, offset=qs.HEADER_LEN)


def rsid_lookup(read, rsid_first: np.ndarray, rsid_n: int, rs_number: int,
                block_records: int = RSID_BLOCK_RECORDS) -> list[tuple[int, int]]:
    """(ordinal, vidx) of every record with `rs_number`, the way a reader does it with range requests.
    `read(offset, length) -> bytes` reads the rsID index file; `rsid_first` and `rsid_n` come from the
    variant index.

    Start at the last block whose first rs_number is **below** `rs_number` (block 0 when none is): a run
    of one rs_number can cross a block boundary, so a block that starts with the number may have the
    run's start at the end of the block before it. Then read forward, block by block, while the records
    are at or below the number."""
    if rsid_n == 0:
        return []
    n_blocks = len(rsid_first)
    b = max(int(np.searchsorted(rsid_first, rs_number, side="left")) - 1, 0)
    out = []
    while b < n_blocks:
        lo = b * block_records
        k = min(block_records, rsid_n - lo)
        rec = np.frombuffer(read(qs.HEADER_LEN + RSID_DTYPE.itemsize * lo, RSID_DTYPE.itemsize * k), RSID_DTYPE)
        i, j = np.searchsorted(rec["rs_number"], rs_number, "left"), np.searchsorted(rec["rs_number"], rs_number, "right")
        out += [(int(o), int(v)) for o, v in zip(rec["ordinal"][i:j], rec["vidx"][i:j])]
        if j < k:                                   # the run ends inside this block
            break
        b += 1
    return out


# ---- build ------------------------------------------------------------------------------------
def collection_table(refget, collection: str, chroms: list[str]) -> list[dict]:
    """[{name, seq_digest, length}] for `chroms`, from the refget collection."""
    by_name = {r.metadata.name: r.metadata for r in refget.get_collection(collection)}
    missing = [c for c in chroms if c not in by_name]
    if missing:
        raise ValueError(f"collection {collection} has no sequence named {missing}")
    return [{"name": c, "seq_digest": by_name[c].sha512t24u, "length": by_name[c].length} for c in chroms]


def _order(t: pa.Table) -> pa.Table:
    """Cis section then trans-only, each by (pos, ref, alt) byte-wise."""
    idx = pc.sort_indices(t, sort_keys=[("in_cis", "descending"), ("pos", "ascending"), ("ref", "ascending"),
                                         ("alt", "ascending")])
    return t.take(idx)


def _check_rsid(t: pa.Table) -> None:
    """A `rsid` column is not stored; it must be recoverable from `rs_number` as "rs<rs_number>", with
    null exactly where `rs_number` is -1 (or null)."""
    if "rsid" not in t.column_names:
        return
    if "rs_number" not in t.column_names:
        raise ValueError("sites: `rsid` without `rs_number`; the variant catalog stores only rs_number")
    num = t["rs_number"].fill_null(-1)
    want = pc.if_else(pc.greater(num, 0), pc.binary_join_element_wise("rs", pc.cast(num, pa.string()), ""),
                      pa.scalar(None, pa.string()))
    have = t["rsid"]
    bad = pc.invert(pc.fill_null(pc.equal(have, want), False))
    bad = pc.and_(bad, pc.invert(pc.and_(pc.is_null(have), pc.is_null(want))))
    n_bad = pc.sum(pc.cast(bad, pa.int64())).as_py() or 0
    if n_bad:
        i = pc.index(bad, True).as_py()
        raise ValueError(f"sites: {n_bad} rows whose rsid is not \"rs<rs_number>\", e.g. "
                         f"{t['rsid'][i].as_py()!r} with rs_number {t['rs_number'][i].as_py()}")


def _column(t: pa.Table, name: str, fill, dtype) -> np.ndarray:
    if name not in t.column_names:
        return np.full(t.num_rows, fill, dtype=dtype)
    return t[name].fill_null(fill).to_numpy(zero_copy_only=False).astype(dtype)


def build(store: qs.Store, cat_id: str, sites: Path | pa.Table, refget, collection: str, chroms: list[str],
          page_size: int = PAGE_SIZE, level: int = ZSTD_LEVEL, source: dict | None = None) -> dict:
    """Write the chromosome files, the variant and rsID indexes, and `variant_catalogs/<cat_id>.json`."""
    if isinstance(sites, pa.Table):
        table = sites
    else:
        table = pq.read_table(sites, filters=[("chr", "in", list(chroms))])
    missing = [c for c in REQUIRED if c not in table.column_names]
    if missing:
        raise ValueError(f"sites: required columns missing: {missing}")
    attributes = [a for a in KNOWN_ATTRIBUTES if a in table.column_names]
    _check_rsid(table)
    seqs = collection_table(refget, collection, chroms) if refget is not None else None
    if seqs is None:
        raise ValueError("catalog: a refgetstore is required to name the sequences")
    other = sorted(set(pc.unique(table["chr"]).to_pylist()) - set(chroms))
    if other:
        raise ValueError(f"sites: chromosomes {other} are not in the configured list")

    entries, ident_parts, rs_parts = [], [], []
    for ordinal, s in enumerate(seqs, start=1):
        t = _order(table.filter(pc.equal(table["chr"], s["name"])))
        n = t.num_rows
        pos = t["pos"].to_numpy().astype(np.int64)
        ref, alt = t["ref"].to_pylist(), t["alt"].to_pylist()
        in_cis = t["in_cis"].to_numpy(zero_copy_only=False)
        n_cis = int(in_cis.sum())
        if n and (pos.min() < 1 or pos.max() > s["length"]):
            raise ValueError(f"sites {s['name']}: position outside 1..{s['length']}")
        keys = list(zip(pos.tolist(), ref, alt))
        if len(set(keys)) != n:
            raise ValueError(f"sites {s['name']}: (pos, ref, alt) is not unique")
        af = _column(t, "af", np.nan, np.float64)
        ms, mc = _column(t, "ma_samples", -1, np.int64), _column(t, "ma_count", -1, np.int64)
        rs = _column(t, "rs_number", -1, np.int64)
        if "match" in attributes:
            m = t["match"].fill_null("none").to_pylist()
            match = np.fromiter((pf.MATCH_CODES[x] for x in m), dtype=np.uint8, count=n)
        else:
            match = np.zeros(n, dtype=np.uint8)
        parts = [qs.file_header(KIND_VARIANTS, s["name"], n, page_size, n_cis, s["seq_digest"])]
        offs, firsts = [qs.HEADER_LEN], []
        for lo, hi in ((0, n_cis), (n_cis, n)):
            for a in range(lo, hi, page_size):
                b = min(hi, a + page_size)
                pg = encode_page(a, pos[a:b], rs[a:b], af[a:b], ms[a:b], mc[a:b], ref[a:b], alt[a:b], match[a:b],
                                 "zstd", level)
                parts.append(pg)
                offs.append(offs[-1] + len(pg))
                firsts.append(int(pos[a]))
        name = store.put(b"".join(parts), EXT_VARIANTS)
        entries.append({**s, "count": n, "n_cis": n_cis, "file": name, "file_digest": name.split(".")[0],
                        "_idx": {"name": s["name"], "n_cis": n_cis, "n_trans": n - n_cis, "page_off": offs,
                                 "page_first_position": firsts}})
        ident_parts.append((s["seq_digest"], pos, ref, alt))
        has = rs > 0
        rs_parts.append((rs[has], np.flatnonzero(has), np.full(int(has.sum()), ordinal)))

    ident = qs.catalog_identity(ident_parts)
    rs_all = np.concatenate([p[0] for p in rs_parts]) if rs_parts else np.zeros(0, np.int64)
    rbytes, rsid_first = encode_rsid(rs_all, np.concatenate([p[1] for p in rs_parts]),
                                     np.concatenate([p[2] for p in rs_parts]), collection)
    vbytes = encode_vidx([e["_idx"] for e in entries], page_size, rsid_first, len(rs_all), collection, level)
    for e in entries:
        del e["_idx"]
    doc = {
        "id": cat_id,
        "identity_digest": ident,
        "collection_digest": collection,
        "orientation": "ref_alt",
        "attributes": attributes,
        "page_size": page_size,
        "n_sites": sum(e["count"] for e in entries),
        "chromosomes": entries,
        "vidx": store.put(vbytes, EXT_VIDX),
        "rsid": store.put(rbytes, EXT_RSID),
        "source": source or {},
    }
    store.write_pointer("variant_catalogs", cat_id, doc)
    return doc


def load_chrom(store: qs.Store, doc: dict, chrom: str) -> dict:
    for c in doc["chromosomes"]:
        if c["name"] == chrom:
            return decode_file((store.immutable / c["file"]).read_bytes())
    raise KeyError(f"variant catalog {doc['id']}: no chromosome {chrom}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--store", required=True, type=Path)
    ap.add_argument("--id", required=True)
    args = ap.parse_args(argv)
    from .common import CHROMS, Config
    from .steps_refget import open_store
    cfg = Config()
    doc = build(qs.Store(args.store), args.id, args.sites, open_store(cfg), cfg["reference"]["collection"], CHROMS,
                source={"sites": str(args.sites)})
    print(json.dumps({k: v for k, v in doc.items() if k != "chromosomes"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
