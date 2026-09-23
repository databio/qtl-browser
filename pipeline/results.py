"""Step 7: the results builder -- one experiment's results objects over a built variant catalog.

One code path for every phenotype type (`ge`, `leafcutter`, ...). Reads the contract tables
`nominal`, `permuted`, `credible_sets` and `phenotypes` (CONTRACT.md) plus a variant catalog built by
`pipeline/catalog.py`, and writes

    experiments/<id>.json             mutable pointer: variant catalog, annotation, results by digest, dof, rule
      <digest>.qbe                    one per (phenotype_type, chromosome): v1 header + one block per phenotype
      <digest>.qbh                    one per chromosome: lead and credible-set hits keyed by vidx
      <digest>.arrow.zst              the experiment's search index, one row per phenotype

Blocks
------
Today's block layout (SPEC.md section 8, `packfmt_v1.encode_gene_block`): 64-byte block header, u16 -log10 p
and u16 SE codes (bit 15 the slope sign), 12-byte credible-set records, one zstd details frame. The
block's rows are the variant catalog vidx run `var_start .. var_start + n_rows - 1` covering every nominal row
and every credible-set variant of the phenotype; a vidx inside the run that the phenotype did not test
is a null row (65535, 0xFFFF). For TOPCHeF every run is contiguous and there are no null rows.

`anchor` is 0: v0 stored the gene TSS there, which is annotation. A reader takes the TSS from the
annotation object by `gene_id`.

The details JSON (`v: 1`) carries study fields only, never annotation:

    {"v": 1, "phenotype_type", "phenotype_id", "phenotype_object_id", "gene_id", "has_nominal",
     "n_nominal", "n_credible_sets", "extra": {...},
     "group": {"lead_phenotype_id", "n_variants", "p_perm", "p_beta", "significant",
               "lead": {"chr", "pos", "ref", "alt"}} | null}

`group` is the permutation row of the phenotype's group, stored once per phenotype but labelled as the
group's: a non-lead intron does not get the cluster's lead variant written onto it as if it were its own.

A phenotype with no nominal rows still gets a block (n_rows 0, or null rows covering its credible
sets) and a search-index row, so the browser can say "no per-variant data published".

Hits (`.qbh`, kind 7)
---------------------
v1 header (count = records, page size = 1024 variants per frame, byte 24 = the chromosome's variant
count), a u32 frame offset table, then one zstd frame per 1024 vidx (zero bytes when empty) of 20-byte
records sorted by (vidx, kind, ord, cs_id): `u32 vidx, u32 ord, f32 value, f32 beta, u8 kind, u8 cs_id,
u8 flags, u8 0`. kind 0 = lead variant of a group (value = p_perm, flags bit 0 = significant by the
experiment's rule, ord = the lead phenotype), kind 1 = credible-set member (value = pip), kind 2 = trans
association (value = -log10 p, beta = the ALT effect). `ord` is a u32 row number in the experiment's
search index.

Trans (`.qbt`, kind 6)
----------------------
One object per phenotype type: v1 header (chromosome `all`, the collection digest, count = frames),
then one zstd frame per phenotype with trans rows (`packfmt_v1.encode_trans_frame`), a gene's frames
back to back (`trans_frame_order`), located by the search index's `trans_off`/`trans_len`. A phenotype with trans rows and no cis result gets an index row
with `chr` and the block pointers null. SPEC.md sections 8-9.
"""
from __future__ import annotations

import argparse
import json
import math
import operator
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import catalog as cat
from . import gwas
from . import packfmt_v1 as pf
from . import qtlstore as qs
from .annotation import encode as encode_arrow, decode as decode_arrow

KIND_RESULTS, KIND_HITS, KIND_TRANS = qs.KIND_RESULTS, qs.KIND_HITS, qs.KIND_TRANS
EXT_RESULTS, EXT_HITS, EXT_INDEX, EXT_TRANS = "qbe", "qbh", "arrow.zst", "qbt"
DETAILS_VERSION = 1
ZSTD_LEVEL = 19
HIT_LEAD, HIT_CS, HIT_TRANS, HIT_DTYPE = pf.HIT_LEAD, pf.HIT_CS, pf.HIT_TRANS, pf.HIT_DTYPE
HITS_FRAME_VARIANTS = pf.HITS_FRAME_VARIANTS
OPS = {"<": operator.lt, "<=": operator.le}
SIG_COLUMNS = ("p_perm", "p_beta")              # `permuted` columns a significance rule may test
DEFAULT_RULE = {"column": "p_perm", "op": "<", "threshold": 0.05}


def significance_test(rule: dict):
    """The experiment's rule as a test on a permuted group row: `test(group) -> bool`, false for no group
    or a null value. The rule names its column (`p_perm` or `p_beta`); the test reads that column."""
    col, op = rule.get("column"), rule.get("op")
    if col not in SIG_COLUMNS or op not in OPS:
        raise ValueError(f"significance rule {rule}: column must be one of {SIG_COLUMNS}, op one of {list(OPS)}")
    thr, fn = float(rule["threshold"]), OPS[op]

    def test(g) -> bool:
        v = None if g is None else _none(getattr(g, col))
        return v is not None and not math.isnan(v) and fn(v, thr)
    return test
INDEX_SCHEMA = pa.schema([
    ("ord", pa.uint32()), ("phenotype_type", pa.string()), ("phenotype_id", pa.string()),
    ("phenotype_object_id", pa.string()), ("gene_id", pa.string()), ("chr", pa.string()),
    ("has_nominal", pa.bool_()), ("is_group_lead", pa.bool_()), ("significant", pa.bool_()),
    ("p_perm", pa.float64()), ("blk_off", pa.uint32()), ("blk_len", pa.uint32()),
    ("var_start", pa.uint32()), ("n_var", pa.uint32()), ("var_off", pa.uint32()), ("var_len", pa.uint32()),
    ("w_lo", pa.int32()), ("w_hi", pa.int32()),
    ("trans_off", pa.uint32()), ("trans_len", pa.uint32()), ("n_trans", pa.uint32())])


# ---- inputs -----------------------------------------------------------------------------------
def read_tables(tables: Path, chroms: list[str]) -> dict:
    """The four small tables, with columns the contract added defaulted when an older adapter run
    lacks them: `phenotype_object_id = phenotype_id`, `has_nominal` from the nominal table."""
    t = {n: pq.read_table(Path(tables) / f"{n}.parquet").to_pandas()
         for n in ("phenotypes", "permuted", "credible_sets")}
    for n in t:
        if "phenotype_object_id" not in t[n].columns:
            t[n]["phenotype_object_id"] = t[n]["phenotype_id"]
    cs = t["credible_sets"][t["credible_sets"]["chr"].isin(chroms)]
    # one row per (phenotype, cs, site): an adapter removes a source's repeats (CONTRACT.md), so any
    # repeat reaching here is an adapter bug, identical or not
    key = ["phenotype_type", "phenotype_id", "cs_id", "chr", "pos", "ref", "alt"]
    if cs.duplicated(key).any():
        raise ValueError(f"credible_sets: {int(cs.duplicated(key).sum())} rows repeat a (phenotype, cs_id, site)")
    t["credible_sets"] = cs
    return t


def _connect():
    """A DuckDB connection for the nominal joins; QTLB_DUCKDB_MEMORY / QTLB_DUCKDB_THREADS size it."""
    import duckdb
    con = duckdb.connect()
    con.execute(f"SET memory_limit = '{os.environ.get('QTLB_DUCKDB_MEMORY', '12GB')}'; "
                f"SET threads = {int(os.environ.get('QTLB_DUCKDB_THREADS', '4'))}; "
                "SET preserve_insertion_order = false")
    return con


def _nominal_src(tables: Path, chrom: str) -> str | None:
    p = Path(tables) / "nominal" / f"chr={chrom}"
    files = sorted(p.glob("*.parquet")) if p.exists() else []
    if not files:
        return None
    return "read_parquet([" + ", ".join(f"'{f}'" for f in files) + "], hive_partitioning = false)"


# nominal rows joined to the variant catalog (vidx) and to `phenotypes` (code = its row number); equi-joins only
_NOM_JOIN = """FROM {src} n
    LEFT JOIN sites s ON n.pos = s.pos AND n.ref = s.ref AND n.alt = s.alt
    LEFT JOIN ph p ON n.phenotype_type = p.phenotype_type AND n.phenotype_id = p.phenotype_id"""


def nominal_summary(con, tables: Path, chrom: str) -> pd.DataFrame:
    """One row per phenotype with nominal rows on `chrom`: code, lo/hi vidx, n rows. Fails on an orphan
    row (a site not in the variant catalog) or a phenotype not in `phenotypes` (CONTRACT.md: never dropped).
    Expects `sites` (pos, ref, alt, vidx) and `ph` (phenotype_type, phenotype_id, code) registered."""
    src = _nominal_src(tables, chrom)
    cols = {"code": "int64", "lo": "int64", "hi": "int64", "n": "int64"}
    if src is None:
        return pd.DataFrame({c: pd.Series(dtype=t) for c, t in cols.items()})
    j = _NOM_JOIN.format(src=src)
    bad = con.execute(f"SELECT n.pos, n.ref, n.alt {j} WHERE s.vidx IS NULL LIMIT 1").fetchall()
    if bad:
        n_bad = con.execute(f"SELECT count(*) {j} WHERE s.vidx IS NULL").fetchone()[0]
        raise ValueError(f"nominal {chrom}: {n_bad} rows name a site not in the variant catalog, e.g. {list(bad[0])}")
    unk = con.execute(f"SELECT DISTINCT n.phenotype_type, n.phenotype_id {j} WHERE p.code IS NULL LIMIT 3").fetchall()
    if unk:
        raise ValueError(f"nominal rows for phenotypes not in `phenotypes`: {[tuple(x) for x in unk]}")
    return con.execute(f"SELECT p.code, min(s.vidx) lo, max(s.vidx) hi, count(*) n {j} GROUP BY p.code").df() \
        .astype(cols)


def nominal_arrays(con, tables: Path, chrom: str) -> dict:
    """The chromosome's nominal rows as compact arrays sorted by (code, vidx): code, vidx, pvalue, beta,
    se (null -> NaN), plus `runs`, code -> (start, end). About 24 bytes a row: chr1 of TOPCHeF (58.5M
    rows) is 1.4 GB, where a pandas frame with the string columns was ~45 GB."""
    src = _nominal_src(tables, chrom)
    if src is None:
        return {"runs": {}}
    t = con.execute(f"SELECT p.code, s.vidx, n.pvalue, n.beta, n.se {_NOM_JOIN.format(src=src)}").arrow()
    if not isinstance(t, pa.Table):                    # some duckdb versions return a reader
        t = t.read_all()
    code = t["code"].to_numpy().astype(np.int32)
    vidx = t["vidx"].to_numpy().astype(np.int64)
    order = np.lexsort((vidx, code))
    out = {"code": code[order], "vidx": vidx[order]}
    del code, vidx
    for c in ("pvalue", "beta", "se"):
        out[c] = t[c].to_numpy(zero_copy_only=False).astype(np.float64)[order]
    del t, order
    kv = out["code"]
    starts = np.r_[0, np.flatnonzero(kv[1:] != kv[:-1]) + 1] if len(kv) else np.zeros(0, np.int64)
    ends = np.r_[starts[1:], len(kv)]
    out["runs"] = {int(kv[a]): (int(a), int(b)) for a, b in zip(starts, ends)}
    return out


def _vidx(sites: pd.DataFrame, df: pd.DataFrame, what: str, pos="pos", ref="ref", alt="alt") -> np.ndarray:
    """Catalog vidx of every row of `df`; an orphan is a hard error (CONTRACT.md: never dropped)."""
    k = df[[pos, ref, alt]].rename(columns={pos: "pos", ref: "ref", alt: "alt"})
    k = k.astype({"pos": np.int64})
    m = k.merge(sites, on=["pos", "ref", "alt"], how="left", validate="many_to_one")
    bad = m["vidx"].isna()
    if bad.any():
        raise ValueError(f"{what}: {int(bad.sum())} rows name a site not in the variant catalog, e.g. "
                         f"{m[bad].iloc[0][['pos', 'ref', 'alt']].tolist()}")
    return m["vidx"].to_numpy().astype(np.int64)


def _page_range(ci: dict, page_size: int, lo: int, hi: int) -> tuple[int, int]:
    """(var_off, var_len): the variant catalog file byte range of the pages holding vidx lo..hi."""
    n_cis = ci["n_cis"]
    if (lo < n_cis) != (hi < n_cis):
        raise ValueError(f"vidx run {lo}..{hi} straddles the cis/trans-only boundary {n_cis}")
    pc_ = -(-n_cis // page_size)
    page = (lambda v: v // page_size) if lo < n_cis else (lambda v: pc_ + (v - n_cis) // page_size)
    a, b = page(lo), page(hi)
    return int(ci["page_off"][a]), int(ci["page_off"][b + 1] - ci["page_off"][a])


def _clean(x):
    return pf._json_clean(x)


# ---- build ------------------------------------------------------------------------------------
def build(store: qs.Store, exp_id: str, tables: Path, cat_id: str, annot_id: str, chroms: list[str],
          ingestion: dict | None = None, level: int = ZSTD_LEVEL) -> dict:
    tables = Path(tables)
    if ingestion is None:
        ip = tables / "ingestion.json"
        ingestion = json.loads(ip.read_text()) if ip.exists() else {}
    cdoc = store.load("variant_catalogs", cat_id)
    table = {c["name"]: c for c in cdoc["chromosomes"]}
    missing = [c for c in chroms if c not in table]
    if missing:
        raise ValueError(f"variant catalog {cat_id} has no chromosome {missing}")
    chroms = [c["name"] for c in cdoc["chromosomes"] if c["name"] in chroms]      # variant catalog table order
    vidx_doc = cat.decode_vidx((store.immutable / cdoc["vidx"]).read_bytes(), [c["name"] for c in cdoc["chromosomes"]])
    page_size = vidx_doc["page_size"]
    rule = ingestion.get("significance") or DEFAULT_RULE
    sig = significance_test(rule)
    t = read_tables(tables, chroms)
    ph, perm, cs = t["phenotypes"], t["permuted"], t["credible_sets"]
    if ph.duplicated(["phenotype_type", "phenotype_id"]).any():
        raise ValueError("phenotypes: (phenotype_type, phenotype_id) is not unique")
    if perm.duplicated(["phenotype_type", "phenotype_object_id"]).any():
        raise ValueError("permuted: more than one row per (phenotype_type, phenotype_object_id)")
    groups = {(r.phenotype_type, r.phenotype_object_id): r for r in perm.itertuples(index=False)}
    types = list(dict.fromkeys(ingestion.get("phenotype_types") or sorted(ph["phenotype_type"].unique())))
    trans_src = _trans_src(tables)
    trans_keys = _trans_phenotypes(trans_src, ph)

    # Pass 1, per chromosome: where each phenotype's nominal rows lie (vidx range and count, from a
    # DuckDB aggregate), its credible sets and the permuted leads' vidx. Only small summaries are kept;
    # the nominal rows themselves are read again, one chromosome at a time, in pass 2.
    ph = ph.reset_index(drop=True)
    code_of = {k: i for i, k in enumerate(zip(ph.phenotype_type, ph.phenotype_id))}
    con = _connect()
    con.register("ph", pd.DataFrame({"phenotype_type": ph.phenotype_type, "phenotype_id": ph.phenotype_id,
                                     "code": np.arange(len(ph), dtype=np.int64)}))
    per_chrom, located = {}, {}
    for c in chroms:
        d = cat.load_chrom(store, cdoc, c)
        sites = pd.DataFrame({"pos": d["pos"].astype(np.int32), "ref": d["ref"], "alt": d["alt"],
                              "vidx": np.arange(len(d["pos"]), dtype=np.int64)})
        con.register("sites", sites)
        summ = nominal_summary(con, tables, c)
        con.unregister("sites")
        runs = {}
        for code, lo, hi, n in summ.itertuples(index=False):
            k = (ph.phenotype_type.iat[code], ph.phenotype_id.iat[code])
            if k in located:
                raise ValueError(f"nominal: {k} has rows on {located[k]} and {c}")
            located[k] = c
            runs[k] = (int(lo), int(hi), int(n))
        sites_pd = sites.astype({"pos": np.int64})
        csc = cs[cs["chr"] == c].copy()
        csc["vidx"] = _vidx(sites_pd, csc, f"credible_sets {c}")
        cs_by = {k: g.sort_values(["vidx", "cs_id"]) for k, g in csc.groupby(["phenotype_type", "phenotype_id"])}
        leads = perm[perm["lead_chr"] == c]
        lead_vidx = _vidx(sites_pd, leads, f"permuted leads {c}", "lead_pos", "lead_ref", "lead_alt")
        per_chrom[c] = {"pos": np.asarray(d["pos"], dtype=np.int64), "runs": runs, "cs": cs_by,
                        "leads": leads, "lead_vidx": lead_vidx}
        del d, sites, sites_pd, csc

    # where a phenotype with no nominal rows lives: its group's lead, else its credible sets
    for r in ph.itertuples(index=False):
        k = (r.phenotype_type, r.phenotype_id)
        if k in located:
            continue
        g = groups.get((r.phenotype_type, r.phenotype_object_id))
        if g is not None and g.lead_chr in per_chrom:
            located[k] = g.lead_chr
        else:
            for c in chroms:
                if k in per_chrom[c]["cs"]:
                    located[k] = c
                    break

    # one entry per phenotype on these chromosomes, in index order
    entries, unplaced = [], []
    for r in ph.itertuples(index=False):
        k = (r.phenotype_type, r.phenotype_id)
        c = located.get(k)
        if c is None and k in trans_keys and (r.phenotype_type, r.phenotype_object_id) not in groups:
            # trans-only: no cis result places it on a chromosome, so its index row has chr null and no
            # block, only the trans pointers (SPEC.md section 8, search index)
            entries.append({"k": k, "code": code_of[k], "chr": None, "row": r, "n_nom": 0, "cs": None, "lo": None,
                            "hi": None, "group": None,
                            "sort": (types.index(r.phenotype_type) if r.phenotype_type in types else len(types),
                                     len(chroms), 0, r.phenotype_id)})
            continue
        if c is None:
            # on a chromosome outside this build, or nothing (nominal, group, credible set) places it:
            # the second kind is counted in the pointer rather than dropped silently
            if (r.phenotype_type, r.phenotype_object_id) not in groups and not any(
                    k in per_chrom[x]["runs"] for x in chroms):
                unplaced.append(r.phenotype_id)
            continue
        pc_ = per_chrom[c]
        run = pc_["runs"].get(k)
        csr = pc_["cs"].get(k)
        los, his = [], []
        if run is not None:
            los.append(run[0])
            his.append(run[1])
        if csr is not None:
            los.append(int(csr["vidx"].min()))
            his.append(int(csr["vidx"].max()))
        lo = min(los) if los else None
        hi = max(his) if his else None
        g = groups.get((r.phenotype_type, r.phenotype_object_id))
        sort_pos = int(pc_["pos"][lo]) if lo is not None else (int(g.lead_pos) if g is not None else 0)
        entries.append({"k": k, "code": code_of[k], "chr": c, "row": r, "n_nom": 0 if run is None else run[2],
                        "cs": csr, "lo": lo, "hi": hi, "group": g,
                        "sort": (types.index(r.phenotype_type) if r.phenotype_type in types else len(types),
                                 chroms.index(c), sort_pos, r.phenotype_id)})
    entries.sort(key=lambda e: e["sort"])
    ords = {e["k"]: i for i, e in enumerate(entries)}
    by_chrom_type = {}
    for e in entries:
        by_chrom_type.setdefault((e["chr"], e["k"][0]), []).append(e)

    # Pass 2, results objects, one per (phenotype_type, chromosome): one chromosome's nominal rows in
    # memory at a time
    dofs = {k: _dof(v) for k, v in (ingestion.get("dof") or {}).items()}
    index_rows = [None] * len(entries)
    acc = {pt: {"files": {}, "n_blocks": 0, "nlp_err": 0.0, "se_err": 0.0, "slope_err": 0.0, "slope_n": 0}
           for pt in types}
    for c in chroms:
        pc_ = per_chrom[c]
        ci = vidx_doc["chroms"][c]
        nom = None
        if pc_["runs"]:
            d = cat.load_chrom(store, cdoc, c)
            con.register("sites", pd.DataFrame({"pos": d["pos"].astype(np.int32), "ref": d["ref"], "alt": d["alt"],
                                                "vidx": np.arange(len(d["pos"]), dtype=np.int64)}))
            del d
            nom = nominal_arrays(con, tables, c)
            con.unregister("sites")
        for ptype in types:
            a = acc[ptype]
            blocks, off = [], qs.HEADER_LEN
            for e in by_chrom_type.get((c, ptype), []):
                span = None if nom is None else nom["runs"].get(e["code"])
                blk, info = _block(e, nom, span, pc_["pos"], sig, level, dofs.get(ptype, (None, None))[0])
                h = pf.parse_block_header(blk)
                a["nlp_err"] = max(a["nlp_err"], h["nlp_max"] / (2 * pf.NLP_MAXQ))
                a["se_err"] = max(a["se_err"], math.expm1((h["lse_max"] - h["lse_min"]) / (2 * pf.SE_MAXQ)))
                a["slope_err"] = max(a["slope_err"], info["slope_err"])
                a["slope_n"] += info["slope_n"]
                var = (None, None) if e["lo"] is None else _page_range(ci, page_size, e["lo"], e["hi"])
                r, g = e["row"], e["group"]
                index_rows[ords[e["k"]]] = {
                    "ord": ords[e["k"]], "phenotype_type": ptype, "phenotype_id": r.phenotype_id,
                    "phenotype_object_id": r.phenotype_object_id, "gene_id": _none(r.gene_id), "chr": c,
                    "has_nominal": info["has_nominal"],
                    "is_group_lead": g is not None and g.phenotype_id == r.phenotype_id,
                    "significant": sig(g),
                    "p_perm": None if g is None else _none(g.p_perm),
                    "blk_off": off, "blk_len": len(blk), "var_start": e["lo"],
                    "n_var": None if e["lo"] is None else e["hi"] - e["lo"] + 1,
                    "var_off": var[0], "var_len": var[1],
                    "w_lo": None if e["lo"] is None else int(pc_["pos"][e["lo"]]),
                    "w_hi": None if e["lo"] is None else int(pc_["pos"][e["hi"]])}
                blocks.append(blk)
                off += len(blk)
            if off > pf.U32_MAX:
                raise ValueError(f"results {ptype} {c}: over 4 GiB")
            body = qs.file_header(KIND_RESULTS, c, len(blocks), 0, 0, table[c]["seq_digest"]) + b"".join(blocks)
            a["files"][c] = store.put(body, EXT_RESULTS)
            a["n_blocks"] += len(blocks)
            del blocks, body
        del nom
    # trans-only phenotypes: an index row with no chromosome and no block
    for e in entries:
        if e["chr"] is None:
            r = e["row"]
            index_rows[ords[e["k"]]] = {
                "ord": ords[e["k"]], "phenotype_type": r.phenotype_type, "phenotype_id": r.phenotype_id,
                "phenotype_object_id": r.phenotype_object_id, "gene_id": _none(r.gene_id), "chr": None,
                "has_nominal": False, "is_group_lead": False, "significant": False, "p_perm": None,
                "blk_off": None, "blk_len": None, "var_start": None, "n_var": None, "var_off": None, "var_len": None,
                "w_lo": None, "w_hi": None}
    for row in index_rows:
        row.update({"trans_off": None, "trans_len": None, "n_trans": None})
    trans_doc, trans_hits = build_trans(store, con, trans_src, cdoc, chroms, ords, index_rows, types, level)
    gwas_doc = gwas.build(store, tables, cdoc, chroms, level=level)
    con.close()
    results = []
    for ptype in types:
        a = acc[ptype]
        results.append({"phenotype_type": ptype, "dof": dofs.get(ptype, (None, None))[0],
                        "dof_fit": dofs.get(ptype, (None, None))[1],
                        "dof_source": "ingestion.json" if ptype in dofs else None,
                        "n_phenotypes": a["n_blocks"], "files": a["files"],
                        "precision": _precision(a, dofs.get(ptype, (None, None))[0]),
                        "trans": None if trans_doc is None else trans_doc["types"].get(ptype)})

    # hits, one object per chromosome
    hits = {}
    for c in chroms:
        pc_ = per_chrom[c]
        lv = [], []
        for v, g in zip(pc_["lead_vidx"].tolist(), pc_["leads"].itertuples(index=False)):
            k = (g.phenotype_type, g.phenotype_id)
            if k not in ords:
                raise ValueError(f"permuted: lead phenotype {k} is not in `phenotypes`")
            lv[0].append(v)
            lv[1].append((ords[k], _nan(g.p_perm), int(sig(g))))
        parts = [pf.hit_records(lv[0], [x[0] for x in lv[1]], [x[1] for x in lv[1]], HIT_LEAD,
                                flags=[x[2] for x in lv[1]])]
        for k, csr in pc_["cs"].items():
            if k not in ords:
                raise ValueError(f"credible_sets: {k} is not in `phenotypes`")
            parts.append(pf.hit_records(csr["vidx"].to_numpy(), ords[k], csr["pip"].to_numpy(dtype=np.float64), HIT_CS,
                                        cs_id=csr["cs_id"].to_numpy()))
        parts += trans_hits.get(c, [])
        recs = np.concatenate(parts)
        hits[c] = store.put(encode_hits(c, table[c]["seq_digest"], recs, table[c]["count"], level=level), EXT_HITS)

    index = pa.Table.from_pylist(index_rows, schema=INDEX_SCHEMA)
    if "has_nominal" in ph.columns:
        # the adapter's flag and the nominal rows this build actually found must agree
        declared = dict(zip(zip(ph.phenotype_type, ph.phenotype_id), ph.has_nominal.astype(bool)))
        off = [(r["phenotype_type"], r["phenotype_id"]) for r in index_rows
               if declared[(r["phenotype_type"], r["phenotype_id"])] != r["has_nominal"]]
        if off:
            raise ValueError(f"phenotypes.has_nominal disagrees with nominal for {len(off)}, e.g. {off[:3]}")
    src_annot = (ingestion.get("source") or {}).get("gene_annotation") if isinstance(ingestion.get("source"), dict) \
        else None
    doc = {
        "id": exp_id,
        "catalog": cat_id,
        "catalog_identity": cdoc["identity_digest"],
        "annotation": annot_id,
        # the source's own gene annotation when it is not the one attached (a known gap, not an error)
        "annotation_source": src_annot if src_annot and src_annot != annot_id else None,
        "allele_orientation_source": ingestion.get("allele_orientation_source"),
        "significance": rule,
        "search_index": store.put(encode_arrow(index, level), EXT_INDEX),
        "n_phenotypes": len(entries),
        # phenotypes with no nominal rows, no permuted group and no credible set: nothing says which
        # chromosome they are on, so no build places them (on a subset build this includes phenotypes
        # whose only credible sets lie on other chromosomes)
        "unplaced": {"count": len(unplaced), "examples": unplaced[:5]},
        "hits": hits,
        "results": results,
        # trans rows a subset build could not place (variant or phenotype on a chromosome outside it),
        # and what the adapter left out (variants with no source alleles); null when there is no trans table
        "trans": None if trans_doc is None else {k: v for k, v in trans_doc.items() if k != "types"},
        "trans_excluded": ingestion.get("trans_eqtl_excluded"),
        "gwas": gwas_doc,
        "source": {k: ingestion[k] for k in ("experiment_id", "source") if k in ingestion},
    }
    store.write_pointer("experiments", exp_id, doc)
    return doc


def _dof(v) -> tuple:
    """(dof, fit report) from an ingestion `dof` entry: a published int (TOPCHeF) or a
    `pipeline/dof.py` fit dict (eQTL Catalogue), whose unusable fit stores dof null."""
    if isinstance(v, dict):
        d = v.get("dof_for_manifest", v.get("dof")) if v.get("usable", True) else None
        return (None if d is None else int(d)), {k: v.get(k) for k in ("residual_log10p", "margin", "rows_usable",
                                                                      "n_samples", "implied_covariates", "reason")}
    return (None if v is None else int(v)), None


def _none(x):
    if x is None or (isinstance(x, float) and math.isnan(x)) or x is pd.NA:
        return None
    return x


def _nan(x) -> float:
    x = _none(x)
    return float("nan") if x is None else float(x)


def _block(e: dict, nom: dict | None, span: tuple | None, pos: np.ndarray, sig, level: int,
           dof: int | None = None) -> tuple[bytes, dict]:
    """One phenotype's block. `nom` is `nominal_arrays` for its chromosome and `span` its (start, end)
    there, or None when it has no nominal rows. With a `dof`, the block is decoded again and every
    rebuilt slope is compared to the source `beta` (`info["slope_err"]`, the largest
    |slope - beta| / se; `info["slope_n"]`, the rows compared)."""
    r, g, csr, lo, hi = e["row"], e["group"], e["cs"], e["lo"], e["hi"]
    extra = r.extra if hasattr(r, "extra") else None
    extra = json.loads(extra) if isinstance(extra, str) and extra else {}
    n_nom = 0 if span is None else span[1] - span[0]
    if n_nom != e["n_nom"]:
        raise ValueError(f"nominal: {e['k']} has {n_nom} rows in pass 2 and {e['n_nom']} in pass 1")
    group = None
    if g is not None:
        group = {"lead_phenotype_id": g.phenotype_id, "n_variants": _none(g.n_variants), "p_perm": _none(g.p_perm),
                 "p_beta": _none(g.p_beta), "significant": sig(g),
                 "lead": {"chr": g.lead_chr, "pos": int(g.lead_pos), "ref": g.lead_ref, "alt": g.lead_alt}}
    details = {"v": DETAILS_VERSION, "phenotype_type": r.phenotype_type, "phenotype_id": r.phenotype_id,
               "phenotype_object_id": r.phenotype_object_id, "gene_id": _none(r.gene_id),
               "has_nominal": n_nom > 0, "n_nominal": n_nom,
               "n_credible_sets": 0 if csr is None else int(csr["cs_id"].nunique()),
               "extra": extra, "group": group}
    if lo is None:
        blk = pf.encode_gene_block(details, None, None, [], [], [], None, None, None, level)
        return blk, {"has_nominal": False, "slope_err": 0.0, "slope_n": 0}
    n = hi - lo + 1
    p, beta, se = (np.full(n, np.nan) for _ in range(3))
    if span is not None:
        a, b = span
        rows = nom["vidx"][a:b] - lo
        if np.any(rows[1:] == rows[:-1]):                     # sorted by vidx within a phenotype
            raise ValueError(f"nominal: {e['k']} tests one site twice")
        p[rows] = nom["pvalue"][a:b]
        beta[rows] = nom["beta"][a:b]
        se[rows] = nom["se"][a:b]
    cs_row = cs_pip = cs_id = None
    if csr is not None:
        cs_row = (csr["vidx"].to_numpy() - lo).tolist()
        cs_pip = csr["pip"].to_numpy(dtype=np.float64)
        cs_id = csr["cs_id"].to_numpy().astype(np.int64)
        if cs_id.min() < 0 or cs_id.max() > 127:
            raise ValueError(f"credible_sets: {e['k']} cs_id outside 0..127")
    blk = pf.encode_gene_block(details, lo, 0, p, beta, se, cs_row, cs_pip, cs_id, level,
                               pos_first=int(pos[lo]), pos_last=int(pos[hi]))
    err, n_cmp = slope_error(blk, beta, se, dof)
    return blk, {"has_nominal": n_nom > 0, "slope_err": err, "slope_n": n_cmp}


def slope_error(blk: bytes, beta: np.ndarray, se: np.ndarray, dof: int | None) -> tuple[float, int]:
    """Measured slope error of one encoded block: the largest |slope rebuilt from the stored codes - source
    beta| / source se, over the rows where a slope can be rebuilt (p above 0 and not null, SE present) and
    the source has a finite beta and se > 0; and how many rows that was. (0.0, 0) with `dof` null: no
    slope is rebuilt then (SPEC.md section 8)."""
    if dof is None:
        return 0.0, 0
    d = pf.decode_gene_block(blk, dof)
    ok = np.isfinite(d["slope"]) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if not ok.any():
        return 0.0, 0
    return float(np.max(np.abs(d["slope"][ok] - beta[ok]) / se[ok])), int(ok.sum())


def _precision(a: dict, dof: int | None) -> dict:
    """The `precision` block of one results set (SPEC.md section 13). The first two are worst-case
    bounds from the block scales; `slope_max_error_over_se` is measured on every row against the source
    (null with `dof` null, when no slope is rebuilt)."""
    return {"neglog10p_max_error": a["nlp_err"], "slope_se_max_rel_error": a["se_err"],
            "slope_max_error_over_se": None if dof is None else a["slope_err"],
            "slope_rows_compared": a["slope_n"], "af_max_error": 0.5 / pf.AF_MAXQ}


# ---- trans ------------------------------------------------------------------------------------
TRANS_COLUMNS = ("phenotype_type", "phenotype_id", "chr", "pos", "ref", "alt", "beta", "pvalue")


def _trans_src(tables: Path) -> str | None:
    """The contract `trans` table as a DuckDB source: `trans.parquet` or `trans/**/*.parquet`; None when
    the experiment has no trans results."""
    f, d = Path(tables) / "trans.parquet", Path(tables) / "trans"
    if f.exists():
        return f"read_parquet('{f}')"
    files = sorted(d.rglob("*.parquet")) if d.exists() else []
    if files:
        return "read_parquet([" + ", ".join(f"'{x}'" for x in files) + "], hive_partitioning = false)"
    return None


def _trans_phenotypes(src: str | None, ph: pd.DataFrame) -> set:
    """(phenotype_type, phenotype_id) of every phenotype with trans rows; fails on one not in `phenotypes`
    (CONTRACT.md: trans rows never name an undescribed phenotype)."""
    if src is None:
        return set()
    import duckdb
    con = duckdb.connect()
    keys = {tuple(x) for x in con.execute(f"SELECT DISTINCT phenotype_type, phenotype_id FROM {src}").fetchall()}
    con.close()
    known = set(zip(ph.phenotype_type, ph.phenotype_id))
    bad = sorted(keys - known)
    if bad:
        raise ValueError(f"trans: {len(bad)} phenotypes not in `phenotypes`, e.g. {bad[:3]}")
    return keys


def build_trans(store: qs.Store, con, src: str | None, cdoc: dict, chroms: list[str], ords: dict, index_rows: list,
                types: list[str], level: int) -> tuple[dict | None, dict]:
    """One trans object per phenotype type (kind 6, `.qbt`): a v1 header with chromosome `all` and the
    collection digest, then one zstd frame per phenotype with trans rows, grouped by gene
    (`trans_frame_order`), so a gene page reads one byte range per object. Fills `trans_off`,
    `trans_len`, `n_trans` of `index_rows` and returns (the experiment's `trans` entry, variant-keyed hit
    records by chromosome for the hits files).

    Every row's variant must be a site of the variant catalog, and every phenotype in `phenotypes`. On a
    subset build, rows whose variant lies on a chromosome outside the build, or whose phenotype the build
    did not place, are skipped and counted."""
    if src is None:
        return None, {}
    names = [c["name"] for c in cdoc["chromosomes"]]
    frames = []
    for c in chroms:
        d = cat.load_chrom(store, cdoc, c)
        frames.append(pd.DataFrame({"chr": c, "pos": d["pos"].astype(np.int32), "ref": d["ref"], "alt": d["alt"],
                                    "vidx": np.arange(len(d["pos"]), dtype=np.int64),
                                    "ordinal": names.index(c) + 1, "af_code": d["af_code"].astype(np.int32),
                                    "rs": d["rs_number"].astype(np.int64)}))
        del d
    con.register("allsites_df", pd.concat(frames, ignore_index=True))
    del frames
    con.execute("CREATE OR REPLACE TABLE allsites AS SELECT * FROM allsites_df")
    con.unregister("allsites_df")
    con.register("ordtab", pd.DataFrame({"phenotype_type": [k[0] for k in ords], "phenotype_id": [k[1] for k in ords],
                                         "ord": np.fromiter(ords.values(), dtype=np.int64, count=len(ords))}))
    inb = "(" + ", ".join(f"'{c}'" for c in chroms) + ")"
    total, out_var, out_ph = con.execute(f"""SELECT count(*), count(*) FILTER (WHERE t.chr NOT IN {inb}),
        count(*) FILTER (WHERE t.chr IN {inb} AND o.ord IS NULL)
        FROM {src} t LEFT JOIN ordtab o ON o.phenotype_type = t.phenotype_type AND o.phenotype_id = t.phenotype_id""").fetchone()
    j = f"""FROM {src} t JOIN ordtab o ON o.phenotype_type = t.phenotype_type AND o.phenotype_id = t.phenotype_id
        LEFT JOIN allsites s ON s.chr = t.chr AND s.pos = t.pos AND s.ref = t.ref AND s.alt = t.alt
        WHERE t.chr IN {inb}"""
    bad = con.execute(f"SELECT t.chr, t.pos, t.ref, t.alt {j} AND s.vidx IS NULL LIMIT 1").fetchall()
    if bad:
        n_bad = con.execute(f"SELECT count(*) {j} AND s.vidx IS NULL").fetchone()[0]
        raise ValueError(f"trans: {n_bad} rows name a site not in the variant catalog, e.g. {list(bad[0])}")
    nulls = con.execute(f"""SELECT count(*) {j} AND (t.pvalue IS NULL OR isnan(t.pvalue) OR t.beta IS NULL
        OR NOT isfinite(t.beta) OR t.pvalue < 0 OR t.pvalue > 1)""").fetchone()[0]
    if nulls:
        raise ValueError(f"trans: {nulls} rows with a null or out-of-range p-value or beta")
    header_collection = cdoc["collection_digest"]
    doc = {"rows": 0, "rows_in_source": int(total), "rows_skipped_variant_outside_build": int(out_var),
           "rows_skipped_phenotype_outside_build": int(out_ph), "types": {}}
    hits_by_chrom: dict[str, list] = {}
    for ptype in types:
        t = con.execute(f"""SELECT o.ord, s.ordinal, s.vidx, s.pos, s.ref, s.alt, s.af_code, s.rs, t.pvalue, t.beta
            {j} AND t.phenotype_type = ? ORDER BY o.ord, s.ordinal, s.pos, s.ref, s.alt""", [ptype]).arrow()
        if not isinstance(t, pa.Table):
            t = t.read_all()
        n = t.num_rows
        if n == 0:
            continue
        col = {k: t[k].to_numpy(zero_copy_only=False) for k in ("ord", "ordinal", "vidx", "pos", "af_code", "rs", "pvalue", "beta")}
        ref, alt = t["ref"].to_pylist(), t["alt"].to_pylist()
        del t
        o = col["ord"]
        starts = np.r_[0, np.flatnonzero(o[1:] != o[:-1]) + 1]
        ends = np.r_[starts[1:], n]
        runs = {int(o[a]): (a, b) for a, b in zip(starts.tolist(), ends.tolist())}
        parts, off, nlp_err, beta_err = [], qs.HEADER_LEN, 0.0, 0.0
        for k in trans_frame_order(index_rows, runs):
            a, b = runs[k]
            fr = pf.encode_trans_frame(col["ordinal"][a:b], col["pos"][a:b], col["rs"][a:b],
                                       col["af_code"][a:b], ref[a:b], alt[a:b], col["pvalue"][a:b], col["beta"][a:b], level)
            row = index_rows[int(o[a])]
            row["trans_off"], row["trans_len"], row["n_trans"] = off, len(fr), b - a
            h = pf._TRANS_HEADER.unpack_from(pf.zstd_unframe(fr, None, "trans frame")[:pf.TRANS_HEADER_LEN], 0)
            nlp_err, beta_err = max(nlp_err, h[2] / (2 * pf.NLP_MAXQ)), max(beta_err, h[3] / (2 * pf.BETA_MAXQ))
            parts.append(fr)
            off += len(fr)
        if off > pf.U32_MAX:
            raise ValueError(f"trans {ptype}: over 4 GiB")
        body = qs.file_header(KIND_TRANS, qs.ALL, len(parts), 0, 0, header_collection) + b"".join(parts)
        name = store.put(body, EXT_TRANS)
        del parts, body
        doc["types"][ptype] = {"file": name, "n_rows": int(n), "n_phenotypes": int(len(starts)),
                               "precision": {"neglog10p_max_error": nlp_err, "beta_max_error": beta_err,
                                             "af_max_error": 0.5 / pf.AF_MAXQ}}
        doc["rows"] += int(n)
        with np.errstate(divide="ignore"):
            nlp = -np.log10(col["pvalue"])
        for c in chroms:
            m = col["ordinal"] == names.index(c) + 1
            if m.any():
                hits_by_chrom.setdefault(c, []).append(pf.hit_records(col["vidx"][m], o[m], nlp[m], HIT_TRANS,
                                                                      beta=col["beta"][m]))
        del col, ref, alt
    con.execute("DROP TABLE allsites")
    con.unregister("ordtab")
    return doc, hits_by_chrom


def trans_frame_order(index_rows: list, ords) -> list[int]:
    """The order of one trans object's frames (SPEC.md section 9): phenotypes with a `gene_id` grouped by
    gene, genes in search-index order (a gene's place is the smallest `ord` of its phenotypes, which is its
    eQTL phenotype's when it has one), `ord` order within a gene; then phenotypes with no gene, in `ord`
    order. `ords` are the search-index rows that have frames. A gene's frames are therefore one contiguous
    byte range, whatever the positions of its introns relative to other genes' introns."""
    first: dict[str, int] = {}
    for r in index_rows:
        g = r["gene_id"]
        if g is not None and r["ord"] < first.get(g, r["ord"] + 1):
            first[g] = r["ord"]

    def key(k: int):
        g = index_rows[k]["gene_id"]
        return (0, first[g], k) if g is not None else (1, k, k)
    return sorted(ords, key=key)


def gene_trans_ranges(index_rows: list, gene_id: str) -> dict[str, tuple[int, int]]:
    """What a gene page reads of the trans objects: per phenotype type, (offset, length) from the smallest
    `trans_off` to the largest `trans_off + trans_len` over the gene's phenotypes with frames. One range
    per object because a gene's frames are contiguous (`trans_frame_order`)."""
    out: dict[str, list[int]] = {}
    for r in index_rows:
        if r["gene_id"] == gene_id and r["trans_off"] is not None:
            a = out.setdefault(r["phenotype_type"], [r["trans_off"], r["trans_off"] + r["trans_len"]])
            a[0], a[1] = min(a[0], r["trans_off"]), max(a[1], r["trans_off"] + r["trans_len"])
    return {t: (a, b - a) for t, (a, b) in out.items()}


def check_trans_layout(index_rows: list, doc: dict, sizes: dict[str, int]) -> list[str]:
    """Failures of the trans frame layout (SPEC.md section 9) against the search index: per phenotype type,
    the frames tile the object from byte 64 to its end with no gap or overlap, and each gene's frames form
    one gap-free run. `sizes` maps each trans object name to its size in bytes."""
    fails = []
    for r in doc["results"]:
        t = r.get("trans")
        if not t:
            continue
        ptype = r["phenotype_type"]
        fr = sorted((x["trans_off"], x["trans_len"], x["gene_id"], x["phenotype_id"]) for x in index_rows
                    if x["phenotype_type"] == ptype and x["trans_off"] is not None)
        what = f"trans {ptype}"
        if len(fr) != t["n_phenotypes"]:
            fails.append(f"{what}: {len(fr)} index rows with frames, the pointer says {t['n_phenotypes']}")
        end = qs.HEADER_LEN
        for off, ln, _, pid in fr:
            if off != end:
                fails.append(f"{what}: frame of {pid} starts at {off}, the previous frame ends at {end}")
                break
            end = off + ln
        if end != sizes.get(t["file"], end):
            fails.append(f"{what}: frames end at {end}, the object is {sizes[t['file']]} bytes")
        by_gene: dict[str, list] = {}
        for off, ln, g, pid in fr:
            if g is not None:
                by_gene.setdefault(g, []).append((off, ln))
        split = [g for g, xs in by_gene.items() if any(xs[i][0] != xs[i - 1][0] + xs[i - 1][1] for i in range(1, len(xs)))]
        if split:
            fails.append(f"{what}: {len(split)} genes whose frames are not one contiguous range, e.g. {split[:3]}")
    return fails


def read_trans(store: qs.Store, doc: dict, row: dict) -> dict | None:
    """Decode one phenotype's trans frame from its search-index row, or None when it has no trans rows.
    Adds `chr` (variant chromosome names from the ordinals) and, with a dof, `se` = |beta| / t(p, dof)."""
    if row.get("trans_off") is None:
        return None
    r = next(x for x in doc["results"] if x["phenotype_type"] == row["phenotype_type"])
    with open(store.immutable / r["trans"]["file"], "rb") as f:
        f.seek(row["trans_off"])
        fr = pf.decode_trans_frame(f.read(row["trans_len"]))
    cdoc = store.load(qs.CATALOGS, doc["catalog"])
    names = [c["name"] for c in cdoc["chromosomes"]]
    fr["chr"] = [names[i - 1] for i in fr["ordinal"].tolist()]
    dof = r["dof"]
    fr["se"] = np.full(fr["n"], np.nan) if dof is None else np.abs(fr["beta"]) / pf.t_from_p(fr["p"], dof)
    return fr


# ---- hits -------------------------------------------------------------------------------------
def encode_hits(chrom: str, seq_digest: str, recs: np.ndarray, n_variants: int,
                frame_variants: int = HITS_FRAME_VARIANTS, level: int = ZSTD_LEVEL) -> bytes:
    """Kind 7: header (count = records, page size = variants per frame, the u32 at byte 24 = the
    chromosome's variant count), the frame offset table, then one zstd frame per `frame_variants`
    variant indices (zero bytes when a frame has no records). `recs` is a HIT_DTYPE array."""
    body, _ = pf.encode_hits_body(recs, n_variants, frame_variants, level)
    return qs.file_header(KIND_HITS, chrom, len(recs), frame_variants, n_variants, seq_digest) + body


def hits_frame(buf: bytes, vidx: int) -> np.ndarray:
    """The records of the frame holding `vidx`, from a whole hits file (what a reader does with two range
    requests: the header and offset table, then the frame)."""
    h = qs.parse_file_header(buf)
    offs = pf.hits_frame_table(buf, h["n_cis"], h["page_size"])
    g = vidx // h["page_size"]
    return pf.decode_hits_frame(buf[offs[g]:offs[g + 1]], g * h["page_size"], h["page_size"])


def decode_hits(buf: bytes) -> np.ndarray:
    """Every record of a hits file, frame by frame, checking the frame table and the header count."""
    h = qs.parse_file_header(buf)
    if h["kind"] != KIND_HITS:
        raise ValueError(f"hits: kind {h['kind']}")
    F, n_var = h["page_size"], h["n_cis"]
    offs = pf.hits_frame_table(buf, n_var, F)
    if offs[-1] != len(buf):
        raise ValueError(f"hits: last frame offset {offs[-1]} != file size {len(buf)}")
    parts = [pf.decode_hits_frame(buf[offs[g]:offs[g + 1]], g * F, F, f"hits frame {g}") for g in range(len(offs) - 1)]
    a = np.concatenate(parts) if parts else np.zeros(0, dtype=HIT_DTYPE)
    if len(a) != h["count"]:
        raise ValueError(f"hits: {len(a)} records, header count {h['count']}")
    return a


# ---- readers ----------------------------------------------------------------------------------
def load_index(store: qs.Store, doc: dict) -> pa.Table:
    return decode_arrow((store.immutable / doc["search_index"]).read_bytes(), "search index")


def read_block(store: qs.Store, doc: dict, row: dict) -> dict:
    """Decode one phenotype's block, given its search-index row.

    With `dof: null` (no usable fit, CONTRACT.md "Degrees of freedom") no slope can be rebuilt:
    `slope` comes back all NaN and the reader has the stored -log10 p, SE and sign bit only. Any
    stand-in dof would decode to slopes that look valid and are wrong."""
    r = next(x for x in doc["results"] if x["phenotype_type"] == row["phenotype_type"])
    buf = (store.immutable / r["files"][row["chr"]]).read_bytes()
    dof = r["dof"]
    out = pf.decode_gene_block(buf[row["blk_off"]:row["blk_off"] + row["blk_len"]], dof if dof is not None else 1,
                               expect_blk_len=row["blk_len"], details_version=DETAILS_VERSION)
    if dof is None:
        out["slope"] = np.full(len(out["slope"]), np.nan)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--tables", required=True, type=Path)
    b.add_argument("--store", required=True, type=Path)
    b.add_argument("--id", required=True)
    b.add_argument("--catalog", required=True)
    b.add_argument("--annotation", required=True)
    args = ap.parse_args(argv)
    from .common import CHROMS
    doc = build(qs.Store(args.store), args.id, args.tables, args.catalog, args.annotation, CHROMS)
    print(json.dumps({k: v for k, v in doc.items() if k not in ("hits",)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
