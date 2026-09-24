"""v0 vs v1 format benchmark: bytes per table type and random-access block read speed.

Compares a v0 build (`<derived>/immutable/`, the published packs) against one experiment of a
qtlb v1 store (`pipeline/qtlstore.py`), over a chromosome set (default: every chromosome the
experiment's variant catalog holds).

Sizes
-----
Per chromosome and totalled, per table type:

    variants   v0 variants.<chr>.qbv      vs  v1 variant catalog <digest>.qbv
    eqtl       v0 eqtl.<chr>.qbe          vs  v1 results `ge` <digest>.qbe
    sqtl       v0 sqtl.<chr>.qbs          vs  v1 results `leafcutter` <digest>.qbe
    hits       v0 hits.<chr>.qbh          vs  v1 <digest>.qbh (both: leads, credible sets, trans rows)
    trans      v0 trans.<chr>.qbt         vs  v1 trans <digest>.qbt, one per phenotype type (totals only)
    gwas       v0 gwas.<chr>.qbg          vs  v1 GWAS <digest>.qbg

Plus the objects that are not per chromosome ("global"), each counted once: v0 search index,
variant index, rsID index, GWAS index; v1 experiment search index, variant catalog variant index and rsID
index, and the annotation's genes and exons tables. v0 global objects and the v1 annotation are
genome-wide whatever the chromosome set, so on a subset they are listed but not compared. The
store section counts every object of every experiment once (shared annotation/catalog objects are
not double counted).

Read speed
----------
N random phenotypes per type (eQTL genes, sQTL introns), the same phenotypes in both formats
(sampled from the v1 search index, matched to v0 by gene id / intron id; v0 intron blocks are
found through their gene's `details.splice`, outside the timing). Per phenotype the timed operation
is what a gene/intron page does after startup: pread the block, `decode_gene_block` it (packfmt_v0's or packfmt_v1's, the same code),
pread the covering variant pages, decode them. Both formats use the same block decoder and the same
minimal page decoder (the page payload layout is shared), so differences are the format's, not the
code's. Every block is read once untimed first, so timings are warm page cache (local disk, no
HTTP). Reported: median and p95 ms per phenotype (total, block, pages) and bytes read per phenotype.

Trans (`trans_bench`): a gene page's trans table, the same genes in both formats: v0 reads the
gene's one frame (its eQTL rows and every intron's sQTL rows); v1 reads, per phenotype type, the one
byte range holding the gene's frames (`results.gene_trans_ranges`; a gene's frames are contiguous) and
decodes each frame in it. `requests` counts range reads. GWAS (`gwas_bench`): the GWAS rows in a gene's [w_lo, w_hi] window
through each format's GWAS index, one range read, decode, filter.

Startup: the one-time reads before a first gene page: v0 search index + variant index + GWAS
index; v1 store/experiment/variant catalog/annotation pointers, the gene lookup's directory and one
bucket, one chromosome's genes, exon models and search index part (the catalog's first chromosome,
the largest), the trans-only part, and the GWAS index (SPEC.md sections 6 and 8). The whole-genome
search index and annotation tables are no longer read on a gene page; `sizes` lists them with the
per-chromosome objects (`split`).

    python -m pipeline.bench_store --v0 /scratch/.../derived/immutable --store /scratch/.../store \
        [--experiment topchef] [--chroms all|chr21,chr22] [--n 200] [--reps 5] [--out DIR]
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import struct
import sys
import time
from pathlib import Path

import numpy as np

from . import annotation as an
from . import catalog as cat
from . import packfmt_v0 as pf0
from . import packfmt_v1 as pf
from . import qtlstore as qs

V0_NAME = re.compile(r"^(?P<kind>[a-z_]+)\.(?:(?P<chr>chr[0-9A-Za-z]+)\.)?(?P<hash>[0-9a-f]{16})\.(?P<ext>.+)$")
V0_PER_CHROM = {"variants": "variants", "eqtl": "eqtl", "sqtl": "sqtl", "hits": "hits", "trans": "trans", "gwas": "gwas"}
V0_GLOBAL = ("search_index", "variant_index", "rsid_index", "gwas_index")
TYPE_CATEGORY = {"ge": "eqtl", "leafcutter": "sqtl"}
V1_NOT_PER_CHROM = ("trans",)           # v1 trans is one object per phenotype type, not per chromosome
CATEGORIES = ("variants", "eqtl", "sqtl", "hits", "trans", "gwas")


def _chrom_key(c: str):
    s = c[3:] if c.startswith("chr") else c
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


# ---- inputs -----------------------------------------------------------------------------------
def v0_files(v0: Path) -> tuple[dict, dict]:
    """(per chromosome {category: {chr: path}}, global {name: path}) from a v0 immutable dir."""
    per, glob_ = {k: {} for k in V0_PER_CHROM.values()}, {}
    for p in sorted(v0.iterdir()):
        m = V0_NAME.match(p.name)
        if not m:
            continue
        k, c = m["kind"], m["chr"]
        if c and k in V0_PER_CHROM:
            if c in per[V0_PER_CHROM[k]]:
                raise ValueError(f"v0: two builds of {k}/{c}")
            per[V0_PER_CHROM[k]][c] = p
        elif not c and k in V0_GLOBAL:
            glob_[k] = p
    return per, glob_


def v1_experiment(store: qs.Store, exp_id: str) -> dict:
    """Pointers and object paths of one experiment: per chromosome {category: {chr: path}}, global."""
    doc = store.load("experiments", exp_id)
    cdoc = store.load("variant_catalogs", doc["catalog"])
    adoc = store.load("annotations", doc["annotation"])
    im = store.immutable
    per = {k: {} for k in CATEGORIES}
    per["variants"] = {c["name"]: im / c["file"] for c in cdoc["chromosomes"]}
    per["hits"] = {c: im / f for c, f in (doc.get("hits") or {}).items()}
    for r in doc["results"]:
        k = TYPE_CATEGORY.get(r["phenotype_type"], r["phenotype_type"])
        per.setdefault(k, {})
        for c, f in r["files"].items():
            per[k][c] = im / f
    gw = doc.get("gwas") or {}
    per["gwas"] = {c: im / f for c, f in (gw.get("files") or {}).items()}
    glob_ = {"search_index": im / doc["search_index"], "catalog_vidx": im / cdoc["vidx"],
             "catalog_rsid": im / cdoc["rsid"], "annotation_genes": im / adoc["genes"],
             "annotation_exons": im / adoc["exons"]}
    if gw.get("index"):
        glob_["gwas_index"] = im / gw["index"]
    trans = {f"trans_{r['phenotype_type']}": im / r["trans"]["file"] for r in doc["results"] if r.get("trans")}
    # the per-chromosome objects the browser reads instead of the whole-genome ones (SPEC sections 6, 8)
    split = {f"annotation_{k} {c}": im / n for c, parts in (adoc.get("chroms") or {}).items() for k, n in parts.items()}
    split |= {f"search_index_part {c}": im / e_["file"] for c, e_ in (doc.get("search_index_parts") or {}).items()}
    if doc.get("search_index_trans_only"):
        split["search_index_part trans-only"] = im / doc["search_index_trans_only"]["file"]
    if adoc.get("lookup"):
        split["annotation_lookup"] = im / adoc["lookup"]
    return {"doc": doc, "cdoc": cdoc, "adoc": adoc, "per": per, "global": glob_, "trans": trans, "split": split}


def _size(p: Path | None):
    return None if p is None else p.stat().st_size


def _cmp(a, b) -> dict:
    out = {"v0": a, "v1": b}
    if a is not None and b is not None:
        out["delta"] = b - a
        out["pct"] = round(100.0 * (b - a) / a, 2) if a else None
    return out


# ---- sizes ------------------------------------------------------------------------------------
def sizes(v0: Path, store: qs.Store, exp_id: str, chroms: list[str]) -> dict:
    v0per, v0glob = v0_files(v0)
    e = v1_experiment(store, exp_id)
    per_chrom, totals = {}, {k: [0, 0] for k in CATEGORIES}
    missing = []
    for c in chroms:
        row = {}
        for k in CATEGORIES:
            a = _size(v0per.get(k, {}).get(c))
            b = None if k in V1_NOT_PER_CHROM else _size(e["per"].get(k, {}).get(c))
            if (a is None) != (b is None) and k not in V1_NOT_PER_CHROM:
                missing.append(f"{k}/{c}: v0 {'missing' if a is None else 'ok'}, v1 {'missing' if b is None else 'ok'}")
            row[k] = _cmp(a, b)
            totals[k][0] += a or 0
            totals[k][1] += b or 0
        per_chrom[c] = row
    # v1 trans objects span every chromosome: on a subset the whole object is counted against the
    # subset's v0 gene-chromosome files, so the trans row is comparable only genome-wide
    totals["trans"][1] = sum(p.stat().st_size for p in e["trans"].values())
    tot = {k: _cmp(a, b) for k, (a, b) in totals.items()}
    tot["all"] = _cmp(sum(t[0] for t in totals.values()), sum(t[1] for t in totals.values()))
    genome_wide = sorted(set(chroms)) != sorted(set(v0per["variants"]))
    glob = {"v0": {k: _size(p) for k, p in v0glob.items()},
            "v1": {k: _size(p) for k, p in e["global"].items()},
            "note": ("v0 global objects and the v1 annotation cover the whole genome; v1 search index, "
                     "catalog variant index and rsID index cover only the variant catalog's chromosomes"
                     + ("; on this chromosome subset they are not comparable" if genome_wide else ""))}
    sp = {k: _size(p) for k, p in e["split"].items()}
    split = {"objects": len(sp), "bytes": sum(sp.values()),
             "by_kind": {k: sum(v for n, v in sp.items() if n.split(" ")[0] == k) for k in dict.fromkeys(n.split(" ")[0] for n in sp)}}
    return {"per_chrom": per_chrom, "totals": tot, "global": glob, "split": split, "missing": missing,
            "store": store_sizes(store)}


def store_sizes(store: qs.Store) -> dict:
    """Every experiment's objects, each file counted once; `shared` = objects named by more than one."""
    owners: dict[str, set] = {}
    per_exp = {}
    for p in sorted((store.root / "experiments").glob("*.json")):
        e = v1_experiment(store, p.stem)
        files = ({f for d in e["per"].values() for f in d.values()} | set(e["global"].values()) | set(e["trans"].values())
                 | set(e["split"].values()))
        per_exp[p.stem] = sum(f.stat().st_size for f in files)
        for f in files:
            owners.setdefault(f.name, set()).add(p.stem)
    shared = {n: sorted(o) for n, o in owners.items() if len(o) > 1}
    total = sum((store.immutable / n).stat().st_size for n in owners)
    return {"experiments": per_exp, "unique_bytes": total,
            "shared_objects": {n: {"bytes": (store.immutable / n).stat().st_size, "experiments": o}
                               for n, o in shared.items()},
            "sum_of_experiments": sum(per_exp.values())}


# ---- read benchmark ---------------------------------------------------------------------------
def decode_pages(buf: bytes) -> int:
    """Minimal decode of consecutive variant pages (v0 and v1 share the payload layout: u32 heap_len,
    pos deltas, rs_number, af, ma_samples, ma_count, allele code, flags, heap). Returns the rows."""
    off, rows, total = 0, 0, len(buf)
    while off < total:
        stored_len, first, n, codec, _ = pf._PAGE_HEADER.unpack_from(buf, off)
        end = off + pf.PAGE_HEADER_LEN + stored_len
        stored = buf[off + pf.PAGE_HEADER_LEN:end]
        p = stored if codec == pf.CODECS["raw"] else pf.zstd_unframe(stored, None, "page")
        (heap_len,) = struct.unpack_from("<I", p, 0)
        np.cumsum(np.frombuffer(p, "<u4", n, 4), dtype=np.int64)
        for o, t in ((4 * n, "<u4"), (8 * n, "<u2"), (10 * n, "<u2"), (12 * n, "<u2")):
            np.frombuffer(p, t, n, 4 + o)
        code = np.frombuffer(p, "u1", n, 4 + 14 * n)
        pairs = [pf.SNP_ALLELES.get(c) for c in code.tolist()]
        if heap_len:
            p[4 + 16 * n:].decode("ascii").split("\n")
        rows += n
        del pairs
        off = end + pf.pad4(end - off)
    return rows


class Reader:
    """Open file handles, one pread per range."""

    def __init__(self):
        self.fds: dict[Path, int] = {}

    def read(self, path: Path, off: int, ln: int) -> bytes:
        fd = self.fds.get(path)
        if fd is None:
            fd = self.fds[path] = os.open(path, os.O_RDONLY)
        return os.pread(fd, ln, off)

    def close(self):
        for fd in self.fds.values():
            os.close(fd)


def _ms(t0: float) -> float:
    return (time.perf_counter() - t0) * 1000.0


def plan_reads(v0: Path, store: qs.Store, exp_id: str, chroms: list[str], n: int, seed: int) -> dict:
    """The sampled phenotypes and, per format, where their block and pages are (all untimed)."""
    v0per, v0glob = v0_files(v0)
    e = v1_experiment(store, exp_id)
    doc = e["doc"]
    idx1 = an.decode(e["global"]["search_index"].read_bytes(), "v1 search index").to_pylist()
    idx0 = an.decode(v0glob["search_index"].read_bytes(), "v0 search index").to_pylist()
    g0 = {(r["gene_id"], r["chr"]): r for r in idx0 if r["blk_off"] is not None}
    vi0 = pf0.decode_variant_index(v0glob["variant_index"].read_bytes())
    manifest = v0.parent / "manifest.json"
    dof0 = json.loads(manifest.read_text())["packs"]["dof"] if manifest.exists() else {"eqtl": 435, "sqtl": 480}
    dof1 = {r["phenotype_type"]: r["dof"] for r in doc["results"]}
    files1 = {r["phenotype_type"]: r["files"] for r in doc["results"]}
    rng = np.random.default_rng(seed)
    rdr = Reader()
    plan, skipped = {}, {}
    for ptype, cat_ in TYPE_CATEGORY.items():
        if ptype not in files1:
            continue
        pool = [r for r in idx1 if r["phenotype_type"] == ptype and r["chr"] in chroms and r["has_nominal"]]
        order = rng.permutation(len(pool))
        items, skip = [], 0
        for i in order:
            if len(items) >= n:
                break
            r = pool[int(i)]
            c = r["chr"]
            gid = (r["gene_id"] or "").split(".")[0]
            g = g0.get((gid, c))
            if g is None:
                skip += 1
                continue
            v0eq = v0per["eqtl"][c]
            if cat_ == "eqtl":
                o0 = {"path": v0eq, "off": g["blk_off"], "len": g["blk_len"], "kind": pf0.KIND_EQTL,
                      "dof": dof0["eqtl"], "var": (g["var_off"], g["var_len"])}
            else:
                d = pf0.decode_gene_block(rdr.read(v0eq, g["blk_off"], g["blk_len"]), dof0["eqtl"],
                                         kind=pf0.KIND_EQTL, expect_blk_len=g["blk_len"])["details"]
                sp = next((s for s in d.get("splice", []) if s.get("phenotype_id") == r["phenotype_id"]), None)
                if sp is None:
                    skip += 1
                    continue
                o0 = {"path": v0per["sqtl"][c], "off": sp["blk_off"], "len": sp["blk_len"], "kind": pf0.KIND_SQTL,
                      "dof": dof0["sqtl"], "var": None}
            o0["vpath"] = v0per["variants"][c]
            o0["page_off"], o0["page_size"] = vi0["chroms"][c]["page_off"], vi0["page_size"]
            o1 = {"path": store.immutable / files1[ptype][c], "off": r["blk_off"], "len": r["blk_len"],
                  "kind": None, "dof": dof1[ptype] or 1, "var": (r["var_off"], r["var_len"]),
                  "vpath": e["per"]["variants"][c], "details_version": 1}
            items.append({"id": r["phenotype_id"], "chr": c, "v0": o0, "v1": o1})
        plan[cat_] = items
        skipped[cat_] = skip
    rdr.close()
    return {"plan": plan, "skipped": skipped, "pool": {k: len(v) for k, v in plan.items()}}


def read_one(rdr: Reader, o: dict) -> dict:
    t0 = time.perf_counter()
    blk = rdr.read(o["path"], o["off"], o["len"])
    # the two decoders are the same function body (packfmt_v1 copies v0's block codec), so timings
    # compare the formats, not the code
    if o["kind"] is None:
        b = pf.decode_gene_block(blk, o["dof"], expect_blk_len=o["len"], details_version=o["details_version"])
    else:
        b = pf0.decode_gene_block(blk, o["dof"], kind=o["kind"], expect_blk_len=o["len"])
    t_blk = _ms(t0)
    t1 = time.perf_counter()
    if o["var"] is not None and o["var"][0] is not None:
        voff, vlen = o["var"]
    elif b["n_rows"]:
        voff, vlen = pf0.variant_range(o["page_off"], o["page_size"], b["var_start"], b["n_rows"])
    else:
        voff = vlen = 0
    rows = decode_pages(rdr.read(o["vpath"], voff, vlen)) if vlen else 0
    t_pg = _ms(t1)
    return {"block_ms": t_blk, "pages_ms": t_pg, "total_ms": t_blk + t_pg, "block_bytes": o["len"],
            "page_bytes": vlen, "bytes": o["len"] + vlen, "n_rows": int(b["n_rows"]), "page_rows": rows}


def _stats(xs: list[float]) -> dict:
    a = np.asarray(xs, dtype=np.float64)
    if a.size == 0:
        return {}
    return {"median": round(float(np.median(a)), 4), "p95": round(float(np.percentile(a, 95)), 4),
            "mean": round(float(a.mean()), 4), "max": round(float(a.max()), 4)}


def read_bench(v0: Path, store: qs.Store, exp_id: str, chroms: list[str], n: int, reps: int, seed: int) -> dict:
    pl = plan_reads(v0, store, exp_id, chroms, n, seed)
    rdr = Reader()
    out = {}
    for cat_, items in pl["plan"].items():
        res = {}
        for fmt in ("v0", "v1"):
            for it in items:                                  # warm: one untimed read of everything
                read_one(rdr, it[fmt])
            per = []
            for it in items:
                runs = [read_one(rdr, it[fmt]) for _ in range(reps)]
                r = dict(runs[0])
                for k in ("block_ms", "pages_ms", "total_ms"):
                    r[k] = float(np.median([x[k] for x in runs]))
                per.append(r)
            res[fmt] = {k: _stats([r[k] for r in per]) for k in
                        ("total_ms", "block_ms", "pages_ms", "bytes", "block_bytes", "page_bytes", "n_rows")}
        out[cat_] = {"n": len(items), "skipped_no_v0_match": pl["skipped"][cat_], **res}
    rdr.close()
    out["startup"] = startup(v0, store, exp_id, reps)
    out["trans"] = trans_bench(v0, store, exp_id, chroms, n, reps, seed)
    out["gwas"] = gwas_bench(v0, store, exp_id, chroms, n, reps, seed)
    out["method"] = (f"{reps} timed reads per phenotype after one untimed warm read (warm page cache, pread on "
                     "open handles, no HTTP); per-phenotype value = median of reps; same phenotypes both formats; "
                     "same block decoder (decode_gene_block, packfmt_v0 and packfmt_v1 copies) and same page decoder")
    return out


def _timed(fn, reps: int) -> tuple[float, object]:
    fn()                                                   # warm
    ts, out = [], None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fn()
        ts.append(_ms(t0))
    return float(np.median(ts)), out


def trans_bench(v0: Path, store: qs.Store, exp_id: str, chroms: list[str], n: int, reps: int, seed: int) -> dict:
    """The gene page's trans table: v0 reads the gene's one frame (its eQTL rows and every intron's sQTL
    rows); v1 reads one range per phenotype type (`results.gene_trans_ranges`) and decodes the gene's
    frames in it. Same genes, sampled among v0 genes with trans rows on `chroms`."""
    from .results import gene_trans_ranges
    v0per, v0glob = v0_files(v0)
    e = v1_experiment(store, exp_id)
    doc = e["doc"]
    if not e["trans"]:
        return {"skipped": "no v1 trans objects"}
    idx0 = [r for r in an.decode(v0glob["search_index"].read_bytes(), "v0").to_pylist()
            if r.get("trans_off") is not None and r["chr"] in chroms]
    idx1 = an.decode(e["global"]["search_index"].read_bytes(), "v1").to_pylist()
    by_gene: dict[str, list] = {}
    for r in idx1:
        if r["trans_off"] is not None and r["gene_id"]:
            by_gene.setdefault(r["gene_id"], []).append(r)
    files1 = {r["phenotype_type"]: store.immutable / r["trans"]["file"] for r in doc["results"] if r.get("trans")}
    manifest = v0.parent / "manifest.json"
    dof0 = json.loads(manifest.read_text())["packs"]["dof"] if manifest.exists() else {"eqtl": 435, "sqtl": 480}
    rng = np.random.default_rng(seed)
    rdr = Reader()
    res = {"v0": [], "v1": []}
    for i in rng.permutation(len(idx0)):
        if len(res["v0"]) >= n:
            break
        g = idx0[int(i)]
        rows1 = by_gene.get(g["gene_id"])
        if not rows1:
            continue
        v0file = v0per["trans"][g["chr"]]

        def f0():
            fr = pf0.decode_trans_frame(rdr.read(v0file, g["trans_off"], g["trans_len"]), dof0["eqtl"], dof0["sqtl"])
            return g["trans_len"], fr["n_e"] + fr["n_s"]

        ranges = gene_trans_ranges(rows1, g["gene_id"])

        def f1():
            nb = nr = 0
            for ptype, (off, ln) in ranges.items():
                buf = rdr.read(files1[ptype], off, ln)
                for r in rows1:
                    if r["phenotype_type"] == ptype:
                        a = r["trans_off"] - off
                        nr += pf.decode_trans_frame(buf[a:a + r["trans_len"]])["n"]
                nb += ln
            return nb, nr
        for fmt, fn in (("v0", f0), ("v1", f1)):
            ms, (nb, nr) = _timed(fn, reps)
            res[fmt].append({"ms": ms, "bytes": nb, "rows": nr, "requests": 1 if fmt == "v0" else len(ranges)})
    rdr.close()
    return {"n": len(res["v0"]), **{fmt: {k: _stats([x[k] for x in v]) for k in ("ms", "bytes", "rows", "requests")}
                                  for fmt, v in res.items()},
            "rows_total": {fmt: int(sum(x["rows"] for x in v)) for fmt, v in res.items()}}


def gwas_bench(v0: Path, store: qs.Store, exp_id: str, chroms: list[str], n: int, reps: int, seed: int) -> dict:
    """The gene page's GWAS panel: every GWAS row in a gene's [w_lo, w_hi] window, through each format's
    GWAS index (one range read, decode the blocks, filter). Same windows (the v1 `ge` rows' w_lo/w_hi)."""
    from .gwas import decode_index
    v0per, v0glob = v0_files(v0)
    e = v1_experiment(store, exp_id)
    g1 = e["doc"].get("gwas")
    if not g1 or "gwas_index" not in v0glob:
        return {"skipped": "no GWAS in one of the formats"}
    i0 = pf0.decode_gwas_index(v0glob["gwas_index"].read_bytes())
    i1 = decode_index(e["global"]["gwas_index"].read_bytes())
    pool = [r for r in an.decode(e["global"]["search_index"].read_bytes(), "v1").to_pylist()
            if r["phenotype_type"] == "ge" and r["w_lo"] is not None and r["chr"] in chroms
            and r["chr"] in i0["chroms"] and r["chr"] in i1["chroms"]]
    rng = np.random.default_rng(seed)
    rdr = Reader()
    res = {"v0": [], "v1": []}
    for i in rng.permutation(len(pool))[:n]:
        r = pool[int(i)]
        c, lo, hi = r["chr"], r["w_lo"], r["w_hi"]

        def read(fmt):
            if fmt == "v0":
                fp, eo = i0["chroms"][c]
                w, path, nv, dec = pf0.gwas_window(fp, eo, lo, hi), v0per["gwas"][c], i0["n_values"], pf0.decode_gwas_block
                first = pf0.FILE_HEADER_LEN
            else:
                fp, eo = i1["chroms"][c]
                w, path, nv, dec = pf.gwas_window(fp, eo, lo, hi, qs.HEADER_LEN), store.immutable / g1["files"][c], \
                    i1["n_values"], pf.decode_gwas_block
                first = qs.HEADER_LEN
            if w is None:
                return 0, 0
            start, end, a, b = w
            buf = rdr.read(path, a, b - a)
            rows = 0
            for k in range(start, end + 1):
                s0 = (first if k == 0 else int(eo[k - 1])) - a
                blk = dec(buf[s0:int(eo[k]) - a], nv)
                rows += int(np.sum((blk["position"] >= lo) & (blk["position"] <= hi)))
            return b - a, rows
        for fmt in ("v0", "v1"):
            ms, (nb, nr) = _timed(lambda: read(fmt), reps)
            res[fmt].append({"ms": ms, "bytes": nb, "rows": nr})
    rdr.close()
    same = sum(a["rows"] == b["rows"] for a, b in zip(res["v0"], res["v1"]))
    return {"n": len(res["v0"]), "windows_same_row_count": same,
            **{fmt: {k: _stats([x[k] for x in v]) for k in ("ms", "bytes", "rows")} for fmt, v in res.items()},
            "rows_total": {fmt: int(sum(x["rows"] for x in v)) for fmt, v in res.items()}}


def startup(v0: Path, store: qs.Store, exp_id: str, reps: int) -> dict:
    _, v0glob = v0_files(v0)
    e = v1_experiment(store, exp_id)
    c0 = e["cdoc"]["chromosomes"][0]["name"]         # the first (largest) chromosome: the dearest gene page
    steps = {
        "v0": [("search_index", v0glob["search_index"], lambda b: an.decode(b, "v0 search index")),
               ("variant_index", v0glob["variant_index"], pf0.decode_variant_index)]
              + ([("gwas_index", v0glob["gwas_index"], pf0.decode_gwas_index)] if "gwas_index" in v0glob else []),
        "v1": [("pointers (store, experiment, variant catalog, annotation json)", None, None),
               ("gene lookup: directory and one bucket", e["split"]["annotation_lookup"], _lookup_one),
               (f"genes {c0}", e["split"][f"annotation_genes {c0}"], lambda b: an.decode(b, "genes")),
               (f"exon models {c0}", e["split"][f"annotation_exon_models {c0}"], lambda b: an.decode(b, "exon models")),
               (f"search index part {c0}", e["split"][f"search_index_part {c0}"], lambda b: an.decode(b, "index part"))]
              + ([("search index part trans-only", e["split"]["search_index_part trans-only"], lambda b: an.decode(b, "index part"))]
                 if "search_index_part trans-only" in e["split"] else [])
              + ([("gwas_index", e["global"]["gwas_index"], _gwas_index1)] if "gwas_index" in e["global"] else []),
    }
    pointers = [store.root / "store.json", store.root / "experiments" / f"{exp_id}.json",
                store.root / "variant_catalogs" / f"{e['doc']['catalog']}.json",
                store.root / "annotations" / f"{e['doc']['annotation']}.json"]
    out = {}
    for fmt, st in steps.items():
        items, tot_ms, tot_b = [], 0.0, 0
        for name, path, fn in st:
            ts = []
            for _ in range(reps):
                t0 = time.perf_counter()
                if path is None:
                    nb = sum(len(p.read_bytes()) for p in pointers if p.exists())
                    [json.loads(p.read_text()) for p in pointers if p.exists()]
                elif fn is _lookup_one:
                    nb = _lookup_one(path)
                else:
                    buf = path.read_bytes()
                    nb = len(buf)
                    fn(buf)
                ts.append(_ms(t0))
            ms = float(np.median(ts))
            items.append({"object": name, "bytes": nb, "ms": round(ms, 3)})
            tot_ms += ms
            tot_b += nb
        out[fmt] = {"items": items, "bytes": tot_b, "ms": round(tot_ms, 3)}
    return out


def _lookup_one(path: Path) -> int:
    """What a gene page reads of the gene lookup: the header and offset table, then one bucket (the
    key FLNC's), decoded. Returns the bytes read."""
    with open(path, "rb") as f:
        head = f.read(qs.HEADER_LEN)
        n = qs.parse_file_header(head)["count"]
        offs = struct.unpack(f"<{n + 1}I", f.read(4 * (n + 1)))
        b = an.lookup_bucket(an.lookup_key("FLNC"), n)
        f.seek(offs[b])
        an.decode(f.read(offs[b + 1] - offs[b]), "lookup bucket")
    return qs.HEADER_LEN + 4 * (n + 1) + offs[b + 1] - offs[b]


def _gwas_index1(b: bytes):
    from .gwas import decode_index
    return decode_index(b)


# ---- report -----------------------------------------------------------------------------------
def _mb(x) -> str:
    return "-" if x is None else f"{x / 1e6:.2f}"


def _pct(d: dict) -> str:
    return "-" if d.get("pct") is None else f"{d['pct']:+.1f}%"


def markdown(r: dict) -> str:
    s, L = r["sizes"], []
    L += [f"# v0 vs v1 format benchmark: {r['experiment']}, {len(r['chroms'])} chromosomes",
          "", f"v0 `{r['v0']}`; store `{r['store']}`; {r['date']}.", "",
          "## Bytes per table type (MB)", "", "| table | v0 | v1 | delta | % |", "|---|---:|---:|---:|---:|"]
    labels = {"variants": "variants", "eqtl": "eQTL", "sqtl": "sQTL", "hits": "hits (both incl. trans rows)",
              "trans": "trans (v1: one object per phenotype type)", "gwas": "GWAS", "all": "**all**"}
    for k, lab in labels.items():
        d = s["totals"][k]
        L.append(f"| {lab} | {_mb(d['v0'])} | {_mb(d['v1'])} | "
                 f"{_mb(d.get('delta')) if d.get('delta') is not None else '-'} | {_pct(d)} |")
    L += ["", "### Global objects (bytes, counted once)", "", "| object | v0 | v1 |", "|---|---:|---:|"]
    for k, v in s["global"]["v0"].items():
        L.append(f"| v0 {k} | {v:,} | |")
    for k, v in s["global"]["v1"].items():
        L.append(f"| v1 {k} | | {v:,} |")
    sp = s.get("split") or {}
    if sp.get("objects"):
        L += ["", f"Per-chromosome objects the browser reads instead (SPEC.md sections 6 and 8): {sp['objects']} objects, "
              f"{sp['bytes']:,} bytes ({', '.join(f'{k} {v:,}' for k, v in sp['by_kind'].items())})."]
    L += ["", s["global"]["note"] + ".", "",
          f"Store: {s['store']['unique_bytes']:,} bytes unique over experiments "
          f"{', '.join(f'{k} {v:,}' for k, v in s['store']['experiments'].items())}; "
          f"{len(s['store']['shared_objects'])} shared objects "
          f"({sum(v['bytes'] for v in s['store']['shared_objects'].values()):,} bytes) counted once.", ""]
    L += ["### Per chromosome (MB, v0 / v1)", "", "| chr | variants | eQTL | sQTL | hits | GWAS | trans v0 |",
          "|---|---|---|---|---|---|---:|"]
    for c, row in s["per_chrom"].items():
        L.append(f"| {c} | " + " | ".join(f"{_mb(row[k]['v0'])} / {_mb(row[k]['v1'])}"
                                         for k in ("variants", "eqtl", "sqtl", "hits", "gwas"))
                 + f" | {_mb(row['trans']['v0'])} |")
    if s["missing"]:
        L += ["", "Missing: " + "; ".join(s["missing"][:20])]
    rb = r.get("reads")
    if rb:
        L += ["", "## Random-access reads (warm cache, local disk)", "", rb["method"] + ".", "",
              "| type | fmt | n | total ms median | p95 | block ms median | pages ms median | bytes/phenotype median | p95 |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
        for k in ("eqtl", "sqtl"):
            if k not in rb:
                continue
            for fmt in ("v0", "v1"):
                d = rb[k][fmt]
                L.append(f"| {k} | {fmt} | {rb[k]['n']} | {d['total_ms']['median']:.3f} | {d['total_ms']['p95']:.3f} | "
                         f"{d['block_ms']['median']:.3f} | {d['pages_ms']['median']:.3f} | "
                         f"{d['bytes']['median']:,.0f} | {d['bytes']['p95']:,.0f} |")
        for k, title, extra in (("trans", "Gene page trans table (v0: one gene frame; v1: one range per phenotype "
                                  "type holding the gene's frames)", ("requests",)),
                                 ("gwas", "Gene page GWAS window (one range read through the GWAS index)", ())):
            t = rb.get(k) or {}
            if "n" not in t:
                L += ["", f"### {title}", "", f"skipped: {t.get('skipped')}"]
                continue
            L += ["", f"### {title}", "", f"{t['n']} genes; rows read in total: v0 {t['rows_total']['v0']:,}, "
                  f"v1 {t['rows_total']['v1']:,}"
                  + (f"; windows with the same row count in both: {t['windows_same_row_count']}" if k == "gwas" else "")
                  + ".", "", "| fmt | ms median | p95 | bytes median | p95 | rows median |"
                  + (" requests median | max |" if extra else ""), "|---|---:|---:|---:|---:|---:|" + ("---:|---:|" if extra else "")]
            for fmt in ("v0", "v1"):
                d = t[fmt]
                L.append(f"| {fmt} | {d['ms']['median']:.3f} | {d['ms']['p95']:.3f} | {d['bytes']['median']:,.0f} | "
                         f"{d['bytes']['p95']:,.0f} | {d['rows']['median']:,.0f} |"
                         + (f" {d['requests']['median']:.0f} | {d['requests']['max']:.0f} |" if extra else ""))
        L += ["", "### Startup (one-time)", "", "| fmt | object | bytes | ms |", "|---|---|---:|---:|"]
        for fmt in ("v0", "v1"):
            for it in rb["startup"][fmt]["items"]:
                L.append(f"| {fmt} | {it['object']} | {it['bytes']:,} | {it['ms']:.2f} |")
            L.append(f"| {fmt} | **total** | {rb['startup'][fmt]['bytes']:,} | {rb['startup'][fmt]['ms']:.2f} |")
    return "\n".join(L) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--v0", required=True, type=Path, help="v0 immutable dir (derived/immutable)")
    ap.add_argument("--store", required=True, type=Path, help="v1 store root")
    ap.add_argument("--experiment", default="topchef")
    ap.add_argument("--chroms", default="all", help="'all' (the variant catalog's chromosomes) or comma list")
    ap.add_argument("--n", type=int, default=200, help="random phenotypes per type for the read benchmark")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-reads", action="store_true", help="sizes only")
    ap.add_argument("--out", type=Path, help="write bench_store.json and bench_store.md here")
    a = ap.parse_args(argv)
    store = qs.Store(a.store)
    cdoc = store.load("variant_catalogs", store.load("experiments", a.experiment)["catalog"])
    have = [c["name"] for c in cdoc["chromosomes"]]
    chroms = have if a.chroms == "all" else [c.strip() for c in a.chroms.split(",")]
    bad = [c for c in chroms if c not in have]
    if bad:
        ap.error(f"variant catalog {cdoc['id']} has no {bad}; it has {have}")
    chroms = sorted(chroms, key=_chrom_key)
    t0 = time.perf_counter()
    r = {"experiment": a.experiment, "v0": str(a.v0), "store": str(a.store), "chroms": chroms,
         "date": dt.datetime.now().isoformat(timespec="seconds"), "host": os.uname().nodename,
         "sizes": sizes(a.v0, store, a.experiment, chroms)}
    if not a.no_reads:
        r["reads"] = read_bench(a.v0, store, a.experiment, chroms, a.n, a.reps, a.seed)
        r["reads"]["params"] = {"n": a.n, "reps": a.reps, "seed": a.seed}
    r["elapsed_s"] = round(time.perf_counter() - t0, 1)
    md = markdown(r)
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "bench_store.json").write_text(json.dumps(r, indent=1, default=str) + "\n")
        (a.out / "bench_store.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
