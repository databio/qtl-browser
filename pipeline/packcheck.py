"""Genome-wide checks of the facts the qtlb pack format (SPEC.md) rests on, and the encoding
measurements that settle its open parameters. Not a build step: it writes a report and evidence
under `data/derived/<packs.checks_dir>` (never uploaded) and leaves no `.done` marker, so it can be
re-run whenever the inputs change.

    uv run python -m pipeline packcheck dof                      # Step 3: one dof per QTL type
    uv run python -m pipeline packcheck run [--type e|s] [--chrom chr7 ...] [--workers 3] [--source auto|raw|derived]
    uv run python -m pipeline packcheck report [--md5]           # merge, pass/fail
    uv run python -m pipeline packcheck measure --chrom chr7 chr10 chr22   # page size and codec (SPEC section 4)
    uv run python -m pipeline packcheck roundtrip --chrom chr7 [--type e|s]   # check 6 on every row (or --steepest 50)

Sources: eQTL and sQTL rows come from the extracted Zenodo nominal files (`raw`). `auto` falls back
to the derived nominal tables (`derived`) when a raw file is missing or its archive is incomplete;
for sQTL that means the significant introns only, and every output records which source it used.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import re
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import yaml
from scipy.special import stdtr, stdtrit

from . import packfmt, steps_pack
from .common import CHROMS, ROOT, Config, connect, die, log, phenotype_batches, register_search_index, variants_path, variants_sql

TYPES = {
    "e": {"name": "eqtl", "raw_dir": "cis_eQTL_nominal", "archive": "cis_eQTL_nominal.tar.gz",
          "pat": "topchef_{c}_MaxPC70.cis_qtl_pairs.{c}.parquet", "derived": "cis_eqtl_nominal", "pcs": 70},
    "s": {"name": "sqtl", "raw_dir": "cis_sQTL_nominal", "archive": "cis_sQTL_nominal.tar.gz",
          "pat": "topchefSplice_{c}_MaxPC25.cis_qtl_pairs.{c}.parquet", "derived": "cis_sqtl_nominal", "pcs": 25},
}
GENE_ID_FROM_INTRON = "split_part(split_part(phenotype_id, ':', 5), '.', 1)"   # as SPLICE_PARSE

# |t| bands for reporting the full-precision SE rebuild (check 3) and the round trip (check 6)
BAND_EDGES = [0.01, 0.1, 1.0]
BAND_NAMES = ["[0, 0.01)", "[0.01, 0.1)", "[0.1, 1)", "[1, inf)"]
NB = len(BAND_NAMES)
REL_TOL = 1e-5
HIST_EDGES = 10.0 ** np.round(np.arange(-16, 16.001, 0.1), 1)   # absolute SE error histogram, 10 bins/decade
DOF_RANGE = np.arange(300, 801)


def checks_dir(cfg: Config) -> Path:
    d = cfg.derived / cfg["packs"]["checks_dir"]
    d.mkdir(parents=True, exist_ok=True)
    return d


def _clean(o):
    """JSON-safe: numpy scalars to Python, non-finite floats to strings."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (int, np.integer)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        if math.isfinite(f):
            return f
        return None if math.isnan(f) else ("inf" if f > 0 else "-inf")
    return o


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(_clean(obj), indent=1))


def _num(x) -> float:
    """Inverse of _clean for floats read back from JSON."""
    if x is None:
        return float("nan")
    if x == "inf":
        return float("inf")
    if x == "-inf":
        return float("-inf")
    return float(x)


# ---- inputs -----------------------------------------------------------------------------------

def archive_state(cfg: Config) -> dict[str, dict]:
    """Size on disk versus sources.yaml for every Zenodo archive."""
    src = yaml.safe_load((cfg.raw / "sources.yaml").read_text())
    q = next(s for s in src["sources"] if s["name"] == "topchef_qtl")
    out = {}
    for f in q["files"]:
        p = cfg.zenodo / f["file"]
        size = p.stat().st_size if p.exists() else None
        out[f["file"]] = {"expected_size": f.get("size"), "size": size, "md5_expected": f.get("md5"),
                          "complete": size is not None and size == f.get("size"),
                          "extracted": (cfg.zenodo / f["file"].removesuffix(".tar.gz")).is_dir()}
    return out


def resolve_source(cfg: Config, t: str, chrom: str, source: str) -> tuple[str, str, int, int] | None:
    """(kind, path or glob, rows, bytes) for one type and chromosome, or None if nothing is there."""
    T = TYPES[t]
    raw = cfg.raw_dir(T["raw_dir"]) / T["pat"].format(c=chrom)
    complete = archive_state(cfg)[T["archive"]]["complete"]
    if source in ("auto", "raw") and raw.exists() and complete:
        return "raw", str(raw), pq.read_metadata(raw).num_rows, raw.stat().st_size
    if source == "raw":
        return None
    files = sorted((cfg.tables / T["derived"] / f"chr={chrom}").glob("bin=*/data.parquet"))
    if not files:
        return None
    return ("derived", str(cfg.tables / T["derived"] / f"chr={chrom}" / "bin=*" / "data.parquet"),
            sum(pq.read_metadata(f).num_rows for f in files), sum(f.stat().st_size for f in files))


def source_sql(t: str, kind: str, src: str) -> str:
    """Rows in one common shape. `frn` is the row's place in the source: the raw file row number,
    or for the derived tables (bin, file row number) packed into one integer, which keeps rows of
    one phenotype comparable because a phenotype never spans bins."""
    if kind == "raw":
        return f"""SELECT phenotype_id, file_row_number::BIGINT AS frn, position::INTEGER AS position, A1, A2,
                   start_distance::INTEGER AS tss_distance, af::DOUBLE AS af, ma_samples::INTEGER AS ma_samples,
                   ma_count::INTEGER AS ma_count, pval_nominal::DOUBLE AS pval_nominal, slope::DOUBLE AS slope,
                   slope_se::DOUBLE AS slope_se
                   FROM read_parquet('{src}', file_row_number = true)"""
    pid = "phenotype_id" if t == "s" else "gene_id AS phenotype_id"
    return f"""SELECT {pid}, bin::BIGINT * 4294967296 + file_row_number AS frn, position::INTEGER AS position, A1, A2,
               tss_distance::INTEGER AS tss_distance, af::DOUBLE AS af, ma_samples::INTEGER AS ma_samples,
               ma_count::INTEGER AS ma_count, pval_nominal::DOUBLE AS pval_nominal, slope::DOUBLE AS slope,
               slope_se::DOUBLE AS slope_se
               FROM read_parquet('{src}', hive_partitioning = true, file_row_number = true)"""


# ---- Step 3: dof ----------------------------------------------------------------------------

def _sample_rows(cfg: Config, t: str, chrom: str, kind: str, src: str) -> pa.Table:
    cols = ["phenotype_id", "pval_nominal", "slope", "slope_se"]
    if kind == "raw":
        tb = pq.ParquetFile(src).read_row_group(0, columns=cols)
        ids = tb["phenotype_id"].combine_chunks()
        last = ids[len(ids) - 1]
        keep = pc.not_equal(ids, last)           # drop the last, possibly partial, phenotype
        return tb.filter(keep)
    # derived: row groups are whole phenotypes already; read them in file order to about 1M rows
    files = sorted(Path(src).parent.parent.glob("bin=*/data.parquet"), key=lambda p: int(p.parent.name[4:]))
    pcol = "phenotype_id" if t == "s" else "gene_id"
    parts, n = [], 0
    for f in files:
        pf = pq.ParquetFile(f)
        for i in range(pf.num_row_groups):
            g = pf.read_row_group(i, columns=[pcol, "pval_nominal", "slope", "slope_se"])
            g = g.rename_columns(cols)
            parts.append(g.cast(pa.schema([("phenotype_id", pa.string()), ("pval_nominal", pa.float64()),
                                           ("slope", pa.float64()), ("slope_se", pa.float64())])))
            n += g.num_rows
            if n >= 1_000_000:
                return pa.concat_tables(parts)
    return pa.concat_tables(parts)


def best_dof(p: np.ndarray, slope: np.ndarray, se: np.ndarray) -> dict:
    ok = (p > 1e-300) & (p < 0.5) & np.isfinite(slope) & (se > 0)
    p, t = p[ok], np.abs(slope[ok] / se[ok])
    if len(p) == 0:
        return {"n": 0, "dof": None}
    lp = np.log10(p)
    obj = np.empty(len(DOF_RANGE))
    with np.errstate(divide="ignore"):
        for i in range(0, len(DOF_RANGE), 50):
            d = DOF_RANGE[i:i + 50, None]
            obj[i:i + 50] = np.median(np.abs(np.log10(2 * stdtr(d, -t[None, :])) - lp[None, :]), axis=1)
    k = int(np.argmin(obj))
    order = np.argsort(obj)
    return {"n": int(len(p)), "dof": int(DOF_RANGE[k]), "obj": float(obj[k]),
            "second_dof": int(DOF_RANGE[order[1]]), "second_obj": float(obj[order[1]])}


def _dof_job(args) -> dict:
    t, chrom, kind, src = args
    cfg = Config()
    tb = _sample_rows(cfg, t, chrom, kind, src)
    ids = tb["phenotype_id"].combine_chunks()
    starts = [0] + (pc.indices_nonzero(pc.not_equal(ids.slice(0, len(ids) - 1), ids.slice(1))).to_numpy() + 1).tolist()
    ends = starts[1:] + [len(ids)]
    p = tb["pval_nominal"].to_numpy(zero_copy_only=False).astype(np.float64)
    s = tb["slope"].to_numpy(zero_copy_only=False).astype(np.float64)
    se = tb["slope_se"].to_numpy(zero_copy_only=False).astype(np.float64)
    out = []
    for k in range(0, len(starts), 10):
        a, b = starts[k], ends[k]
        r = best_dof(p[a:b], s[a:b], se[a:b])
        r.update(phenotype_id=ids[a].as_py(), rows=b - a)
        out.append(r)
    return {"type": t, "chrom": chrom, "source": kind, "phenotypes": out}


def cmd_dof(cfg: Config, args) -> None:
    out_dir = checks_dir(cfg)
    jobs = []
    for t in args.type or ["e", "s"]:
        for c in args.chrom or CHROMS:
            r = resolve_source(cfg, t, c, args.source)
            if r is None:
                log(f"packcheck dof: no {TYPES[t]['name']} source for {c}, skipping")
                continue
            jobs.append((t, c, r[0], r[1]))
    results = []
    with ProcessPoolExecutor(max_workers=args.workers or cfg["workers"]) as ex:
        for f in as_completed([ex.submit(_dof_job, j) for j in jobs]):
            r = f.result()
            results.append(r)
            ds = [x["dof"] for x in r["phenotypes"] if x["dof"] is not None]
            log(f"packcheck dof: {r['type']} {r['chrom']} ({r['source']}): {len(ds)} phenotypes, dof {sorted(set(ds))}")
    report = {"run": dt.datetime.now().isoformat(timespec="seconds"), "types": {}}
    config_dof = dict(cfg["packs"]["dof"])
    new_dof = dict(config_dof)
    for t in ("e", "s"):
        rs = [r for r in results if r["type"] == t]
        if not rs:
            continue
        name = TYPES[t]["name"]
        ph = [x for r in rs for x in r["phenotypes"] if x["dof"] is not None]
        hist: dict[int, int] = {}
        for x in ph:
            hist[x["dof"]] = hist.get(x["dof"], 0) + 1
        mode = max(hist, key=hist.get)
        per_chr = {}
        for r in sorted(rs, key=lambda r: CHROMS.index(r["chrom"])):
            h: dict[int, int] = {}
            for x in r["phenotypes"]:
                if x["dof"] is not None:
                    h[x["dof"]] = h.get(x["dof"], 0) + 1
            per_chr[r["chrom"]] = {"mode": max(h, key=h.get) if h else None, "hist": h, "source": r["source"]}
        off = [x for x in ph if x["dof"] != mode]
        report["types"][name] = {
            "sampled_phenotypes": len(ph), "skipped_no_rows": sum(1 for r in rs for x in r["phenotypes"] if x["dof"] is None),
            "sources": sorted({r["source"] for r in rs}),
            "hist": dict(sorted(hist.items())), "mode": mode, "unanimous": len(hist) == 1,
            "per_chrom": per_chr, "not_mode": off[:50],
            "median_obj_at_mode": float(np.median([x["obj"] for x in ph if x["dof"] == mode])),
            "config": config_dof.get(name), "pcs": TYPES[t]["pcs"],
            "implied_samples_minus_other_covariates": mode + 2 + TYPES[t]["pcs"],
        }
        if mode != config_dof.get(name):
            new_dof[name] = mode
    report["config_before"] = config_dof
    report["config_after"] = new_dof
    report["config_updated"] = new_dof != config_dof
    if report["config_updated"]:
        path = ROOT / "pipeline" / "config.yaml"
        text = path.read_text()
        line = f"dof: {{eqtl: {new_dof['eqtl']}, sqtl: {new_dof['sqtl']}}}"
        text2, n = re.subn(r"dof: \{eqtl: \d+, sqtl: \d+\}", line, text)
        if n != 1:
            die("packcheck dof: could not find the packs.dof line in config.yaml to update")
        path.write_text(text2)
        log(f"packcheck dof: config packs.dof changed {config_dof} -> {new_dof}")
    write_json(out_dir / "dof.json", report)
    for name, r in report["types"].items():
        log(f"packcheck dof: {name}: mode {r['mode']} over {r['sampled_phenotypes']} phenotypes, "
            f"unanimous={r['unanimous']}, hist={r['hist']}, config={r['config']}")


# ---- Step 4: per-chromosome run --------------------------------------------------------------

def _sig3(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(sign, exponent, 3-digit mantissa) of a value printed to 3 significant figures."""
    ax = np.abs(x)
    nz = np.isfinite(ax) & (ax > 0)
    e = np.zeros(x.shape)
    m = np.zeros(x.shape)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        e[nz] = np.floor(np.log10(ax[nz]))
        m[nz] = np.rint(ax[nz] / 10.0 ** (e[nz] - 2))
    roll = m >= 1000
    m[roll] = np.rint(m[roll] / 10)
    e[roll] += 1
    return np.sign(x), e, m


def _changed3(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """True where y printed to 3 significant figures differs from x (x finite)."""
    sx, ex, mx = _sig3(x)
    sy, ey, my = _sig3(y)
    return ~np.isfinite(y) | (sx != sy) | (ex != ey) | (mx != my)


def _new_tally() -> dict:
    z = lambda v=0: [v] * NB                                                     # noqa: E731
    return {
        "rows": 0, "phenotypes": 0, "phenotypes_not_in_variant_order": 0,
        "cats": {"p_null": 0, "p_zero": 0, "p_one": 0, "slope_null": 0, "se_null_or_nonpositive": 0,
                 "t_defined": 0},
        "full": {"n": z(), "max_abs": z(0.0), "max_rel": z(0.0), "n_rel_gt_tol": z(),
                 "pass_domain_rows": 0, "pass_domain_fail": 0, "max_rel_pass_domain": 0.0},
        "worst_full": [],
    }


def _band(abs_t: np.ndarray) -> np.ndarray:
    return np.searchsorted(np.array(BAND_EDGES), abs_t, side="right")


def _se_from_p(slope: np.ndarray, p: np.ndarray, dof: int) -> np.ndarray:
    """Check 3, the source's own relation at full precision: |slope| / t with
    t = -stdtrit(dof, p / 2). NaN where p is NaN, 0, or 1, or slope is NaN."""
    out = np.full(np.broadcast(slope, p).shape, np.nan)
    ok = ~np.isnan(p) & (p > 0) & (p < 1) & ~np.isnan(slope)
    if np.any(ok):
        out[ok] = np.abs(slope[ok]) / -stdtrit(dof, p[ok] / 2.0)
    return out


def _stream_pass(joined: Path, dof: int, tal: dict) -> None:
    """Step 4e: tally the full-precision SE rebuild (check 3) per |t| band over the whole batch.
    Check 6, the round trip through the 16-bit codes, is `packcheck roundtrip`."""
    cols = ["phenotype_id", "frn", "vidx", "pval_nominal", "slope", "slope_se"]

    def process(tb: pa.Table) -> None:
        n = tb.num_rows
        if n == 0:
            return
        ids = tb["phenotype_id"].combine_chunks()
        starts, ends = _slices(ids)
        frn = tb["frn"].to_numpy(zero_copy_only=False)
        vidx = pc.cast(tb["vidx"], pa.int64()).fill_null(-1).to_numpy(zero_copy_only=False)
        p = tb["pval_nominal"].to_numpy(zero_copy_only=False).astype(np.float64)
        slope = tb["slope"].to_numpy(zero_copy_only=False).astype(np.float64)
        se = tb["slope_se"].to_numpy(zero_copy_only=False).astype(np.float64)
        for a, b in zip(starts, ends):
            if b - a > 1 and not np.all(np.diff(frn[a:b]) > 0):
                tal["phenotypes_not_in_variant_order"] += 1
        tal["phenotypes"] += len(starts)
        tal["rows"] += n

        p_null = np.isnan(p)
        p_fin = ~p_null
        cats = tal["cats"]
        cats["p_null"] += int(p_null.sum())
        cats["p_zero"] += int((p == 0).sum())
        cats["p_one"] += int((p == 1).sum())
        s_fin = np.isfinite(slope)
        cats["slope_null"] += int((~s_fin).sum())
        se_ok = np.isfinite(se) & (se > 0)
        cats["se_null_or_nonpositive"] += int((~se_ok).sum())
        tdef = s_fin & se_ok
        cats["t_defined"] += int(tdef.sum())
        abs_t = np.full(n, np.nan)
        abs_t[tdef] = np.abs(slope[tdef] / se[tdef])
        band = np.full(n, -1)
        band[tdef] = _band(abs_t[tdef])

        # full precision
        full_dom = tdef & p_fin & (p > 0) & (p < 1)
        se_hat = _se_from_p(slope, p, dof)
        err = np.abs(se_hat - se)
        rel = err / se
        F = tal["full"]
        for k in range(NB):
            m = full_dom & (band == k)
            if m.any():
                F["n"][k] += int(m.sum())
                F["max_abs"][k] = max(F["max_abs"][k], float(np.nanmax(err[m])))
                F["max_rel"][k] = max(F["max_rel"][k], float(np.nanmax(rel[m])))
                F["n_rel_gt_tol"][k] += int((~(rel[m] <= REL_TOL)).sum())
        pdom = full_dom & (band >= 1)
        F["pass_domain_rows"] += int(pdom.sum())
        if pdom.any():
            F["pass_domain_fail"] += int((~(rel[pdom] <= REL_TOL)).sum())
            F["max_rel_pass_domain"] = max(F["max_rel_pass_domain"], float(np.nanmax(rel[pdom])))
            idx = np.nonzero(pdom)[0]
            top = idx[np.argsort(-rel[idx])[:20]]
            for i in top:
                tal["worst_full"].append({"phenotype_id": ids[int(i)].as_py(), "vidx": int(vidx[i]), "pval_nominal": p[i],
                                          "slope": slope[i], "slope_se": se[i], "se_hat": se_hat[i],
                                          "abs_t": abs_t[i], "rel_err": rel[i]})
            tal["worst_full"] = sorted(tal["worst_full"], key=lambda r: -r["rel_err"])[:20]

    for tb in phenotype_batches(pq.ParquetFile(joined), cols):
        process(tb)


def _worker(job) -> dict:
    t, chrom, kind, src, n_src = job
    cfg = Config()
    T = TYPES[t]
    dof = int(cfg["packs"]["dof"][T["name"]])
    out = checks_dir(cfg)
    work = out / "tmp" / f"{t}_{chrom}"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    t0 = time.time()
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"], temp_dir=work / "duckdb")
    con.execute("SET preserve_insertion_order = true")
    collation = con.execute("SELECT current_setting('default_collation')").fetchone()[0]
    if collation:
        raise RuntimeError(f"default_collation is {collation!r}; the variant order must be byte-wise")
    # the variant table is one file sorted by (chr, position, A1, A2), so one chromosome's rows
    # are contiguous in it and file_row_number still proves the file order is the sorted order
    vpos = variants_path(cfg)

    # 4a. the chromosome's cis variant list, indexed by an explicit byte-wise sort
    con.execute(f"""
        CREATE TABLE v AS
        SELECT (row_number() OVER (ORDER BY position, A1, A2) - 1)::UINTEGER AS vidx,
               position, A1, A2, rs_number, match, file_row_number AS frn
        FROM read_parquet('{vpos}', file_row_number = true) WHERE chr = '{chrom}' AND in_cis""")
    n_var, n_keys = con.execute("SELECT count(*), (SELECT count(*) FROM (SELECT DISTINCT position, A1, A2 FROM v)) FROM v").fetchone()
    file_order = con.execute("SELECT count(*) FROM (SELECT vidx, row_number() OVER (ORDER BY frn) - 1 AS r FROM v) WHERE vidx <> r").fetchone()[0]
    res = {"type": t, "qtl": T["name"], "chrom": chrom, "source": kind, "src": src, "dof": dof,
           "started": dt.datetime.now().isoformat(timespec="seconds"),
           "variants": {"n_cis": n_var, "duplicate_keys": n_var - n_keys, "file_order_equals_vidx": file_order == 0,
                        "rows_out_of_file_order": file_order, "collation": collation}}

    # 4b. raw rows joined to variant indexes, sorted by phenotype then vidx
    joined = work / "joined.parquet"
    con.execute(f"""
        COPY (
          SELECT n.phenotype_id, n.frn, v.vidx, n.position, n.tss_distance,
                 n.af, n.ma_samples, n.ma_count, n.pval_nominal, n.slope, n.slope_se
          FROM ({source_sql(t, kind, src)}) n
          LEFT JOIN v ON v.position = n.position AND v.A1 = n.A1 AND v.A2 = n.A2
          ORDER BY n.phenotype_id, v.vidx NULLS LAST, n.frn
        ) TO '{joined}' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000)""")
    n_joined = pq.read_metadata(joined).num_rows
    if n_joined != n_src:
        raise RuntimeError(f"{t} {chrom}: joined rows {n_joined} != source rows {n_src} (join duplicated or dropped rows)")
    res["rows"] = {"source": n_src, "joined": n_joined}

    # 4c. per-phenotype facts
    ph_path = out / f"phenotypes_{t}_{chrom}.parquet"
    if t == "e":
        perm = f"SELECT gene_id AS phenotype_id, num_var, is_egene AS is_sig FROM '{cfg.tables / 'genes.parquet'}' WHERE chr = '{chrom}' AND tested"
        gid = "phenotype_id"
    else:
        perm = f"SELECT phenotype_id, num_var, is_sqtl AS is_sig FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE chr = '{chrom}'"
        gid = GENE_ID_FROM_INTRON
    con.execute(f"""
        COPY (
          WITH f AS (
            SELECT phenotype_id,
              count(*)                                   AS n_rows,
              count(vidx)                                AS n_matched,
              count(DISTINCT vidx)                       AS n_distinct,
              min(vidx)                                  AS var_start,
              max(vidx)::BIGINT - min(vidx)::BIGINT + 1 - count(*) AS n_gaps,
              max(frn) - min(frn) + 1 = count(*)         AS raw_grouped,
              min(position - tss_distance)               AS anchor,
              min(position - tss_distance) = max(position - tss_distance) AS anchor_constant,
              min(position) AS pos_first, max(position) AS pos_last,
              count(*) FILTER (WHERE pval_nominal IS NULL OR isnan(pval_nominal)) AS n_p_nan,
              count(*) FILTER (WHERE pval_nominal = 0)   AS n_p_zero,
              count(*) FILTER (WHERE pval_nominal = 1)   AS n_p_one,
              count(*) FILTER (WHERE pval_nominal > 1 AND NOT isnan(pval_nominal)) AS n_p_gt1,
              count(*) FILTER (WHERE slope IS NULL OR isnan(slope))       AS n_slope_nan,
              count(*) FILTER (WHERE slope_se IS NULL OR isnan(slope_se)) AS n_se_nan,
              count(*) FILTER (WHERE slope = 0)          AS n_slope_zero,
              count(*) FILTER (WHERE af IS NULL OR isnan(af)) AS n_af_nan,
              count(*) FILTER (WHERE ma_samples IS NULL OR ma_count IS NULL) AS n_counts_null,
              max(-log10(pval_nominal)) FILTER (WHERE pval_nominal > 0 AND NOT isnan(pval_nominal)) AS max_nlp,
              min(pval_nominal) FILTER (WHERE pval_nominal > 0 AND NOT isnan(pval_nominal))         AS min_p,
              max(abs(slope)) FILTER (WHERE NOT isnan(slope))    AS max_abs_slope,
              median(abs(slope)) FILTER (WHERE NOT isnan(slope)) AS median_abs_slope
            FROM '{joined}' GROUP BY 1
          ),
          perm AS ({perm}),
          si AS (SELECT gene_id, tss FROM {register_search_index(cfg, con)})
          SELECT '{chrom}' AS chr, '{kind}' AS source, f.*,
                 (n_matched = n_rows AND n_distinct = n_rows AND n_gaps = 0) AS contiguous,
                 {gid.replace('phenotype_id', 'f.phenotype_id')} AS gene_id,
                 perm.num_var AS perm_num_var, perm.is_sig,
                 f.n_rows <> perm.num_var AS rows_ne_perm,
                 si.tss, f.anchor - si.tss AS anchor_minus_tss
          FROM f LEFT JOIN perm USING (phenotype_id)
                 LEFT JOIN si ON si.gene_id = {gid.replace('phenotype_id', 'f.phenotype_id')}
          ORDER BY f.phenotype_id
        ) TO '{ph_path}' (FORMAT PARQUET, COMPRESSION ZSTD)""")
    cur = con.execute(f"""
        SELECT count(*) AS n, sum(n_rows) AS rows,
          count(*) FILTER (WHERE n_matched <> n_rows) AS with_unmatched_rows, sum(n_rows - n_matched) AS unmatched_rows,
          count(*) FILTER (WHERE n_distinct <> n_matched) AS with_duplicate_vidx,
          count(*) FILTER (WHERE n_gaps <> 0) AS with_gaps, coalesce(sum(n_gaps) FILTER (WHERE n_gaps > 0), 0) AS gap_variants,
          count(*) FILTER (WHERE NOT contiguous) AS not_contiguous,
          count(*) FILTER (WHERE NOT raw_grouped) AS not_grouped,
          count(*) FILTER (WHERE NOT anchor_constant) AS anchor_not_constant,
          count(*) FILTER (WHERE perm_num_var IS NULL) AS no_permutation_row,
          count(*) FILTER (WHERE rows_ne_perm) AS rows_ne_perm_num_var,
          count(*) FILTER (WHERE tss IS NULL) AS no_search_index_tss,
          count(*) FILTER (WHERE anchor_minus_tss <> 0) AS anchor_ne_tss,
          max(abs(anchor_minus_tss)) AS max_abs_anchor_minus_tss,
          sum(n_p_nan) AS rows_p_nan, sum(n_p_zero) AS rows_p_zero, sum(n_p_one) AS rows_p_one, sum(n_p_gt1) AS rows_p_gt1,
          sum(n_slope_nan) AS rows_slope_nan, sum(n_se_nan) AS rows_se_nan, sum(n_slope_zero) AS rows_slope_zero,
          sum(n_af_nan) AS rows_af_nan, sum(n_counts_null) AS rows_counts_null,
          count(*) FILTER (WHERE n_p_nan > 0) AS phenotypes_with_p_nan,
          min(min_p) AS min_p, max(max_nlp) AS max_nlp,
          count(*) FILTER (WHERE max_nlp > 65.533) AS phenotypes_max_nlp_over_65_533,
          max(max_abs_slope) AS max_abs_slope
        FROM '{ph_path}'""")
    row = cur.fetchone()
    res["phenotypes"] = dict(zip([d[0] for d in cur.description], row))
    rne = con.execute(f"SELECT n_rows - perm_num_var AS d, count(*) FROM '{ph_path}' WHERE rows_ne_perm GROUP BY 1 ORDER BY 1").fetchall()
    res["phenotypes"]["rows_minus_perm_num_var_hist"] = {str(d): c for d, c in rne}

    # 4d. variant-level values
    vv_path = out / f"variant_values_{t}_{chrom}.parquet"
    con.execute(f"""
        COPY (
          SELECT vidx, count(*) AS n_phen,
            min(af) AS af, max(af) <> min(af) OR (bool_or(af IS NULL) AND bool_or(af IS NOT NULL)) AS af_differs,
            min(ma_samples) AS ma_samples, max(ma_samples) <> min(ma_samples) AS ma_samples_differs,
            min(ma_count) AS ma_count, max(ma_count) <> min(ma_count) AS ma_count_differs,
            '{kind}' AS source
          FROM '{joined}' WHERE vidx IS NOT NULL GROUP BY vidx ORDER BY vidx
        ) TO '{vv_path}' (FORMAT PARQUET, COMPRESSION ZSTD)""")
    r = con.execute(f"""SELECT count(*), count(*) FILTER (WHERE af_differs), count(*) FILTER (WHERE ma_samples_differs),
        count(*) FILTER (WHERE ma_count_differs), count(*) FILTER (WHERE af_differs OR ma_samples_differs OR ma_count_differs),
        min(af), max(af), count(*) FILTER (WHERE af > 0.5), max(ma_samples), max(ma_count)
        FROM '{vv_path}'""").fetchone()
    res["variant_values"] = dict(zip(["tested", "af_differs", "ma_samples_differs", "ma_count_differs", "any_differs",
                                      "af_min", "af_max", "af_over_0_5", "ma_samples_max", "ma_count_max"], r))
    con.close()

    # 4e. streaming SE / quantization pass
    tal = _new_tally()
    _stream_pass(joined, dof, tal)
    res["stream"] = tal
    res["seconds"] = round(time.time() - t0, 1)
    write_json(out / f"{t}_{chrom}.json", res)
    shutil.rmtree(work, ignore_errors=True)
    return res


def cross_type(cfg: Config, chrom: str) -> dict | None:
    """Step 4f: eQTL versus sQTL on one chromosome, once both types have outputs."""
    out = checks_dir(cfg)
    ve, vs = out / f"variant_values_e_{chrom}.parquet", out / f"variant_values_s_{chrom}.parquet"
    pe, ps = out / f"phenotypes_e_{chrom}.parquet", out / f"phenotypes_s_{chrom}.parquet"
    if not all(p.exists() for p in (ve, vs, pe, ps)):
        return None
    je, js = json.loads((out / f"e_{chrom}.json").read_text()), json.loads((out / f"s_{chrom}.json").read_text())
    both_raw = je["source"] == "raw" and js["source"] == "raw"
    af = "af" if both_raw else "af::FLOAT"        # derived tables store af as float32
    con = connect(cfg, temp_dir=out / "tmp" / f"cross_{chrom}")
    r = con.execute(f"""
        SELECT count(*) FILTER (WHERE e.vidx IS NOT NULL AND s.vidx IS NOT NULL) AS both,
               count(*) FILTER (WHERE e.vidx IS NOT NULL AND s.vidx IS NOT NULL AND
                   ((e.{af}) IS DISTINCT FROM (s.{af}) OR e.ma_samples IS DISTINCT FROM s.ma_samples
                    OR e.ma_count IS DISTINCT FROM s.ma_count)) AS differ,
               count(*) FILTER (WHERE e.vidx IS NOT NULL AND s.vidx IS NOT NULL AND (e.{af}) IS DISTINCT FROM (s.{af})) AS af_differ,
               count(*) FILTER (WHERE e.vidx IS NOT NULL AND s.vidx IS NOT NULL AND e.ma_samples IS DISTINCT FROM s.ma_samples) AS ma_samples_differ,
               count(*) FILTER (WHERE e.vidx IS NOT NULL AND s.vidx IS NOT NULL AND e.ma_count IS DISTINCT FROM s.ma_count) AS ma_count_differ,
               count(*) FILTER (WHERE s.vidx IS NULL) AS only_eqtl,
               count(*) FILTER (WHERE e.vidx IS NULL) AS only_sqtl
        FROM '{ve}' e FULL OUTER JOIN '{vs}' s USING (vidx)""").fetchone()
    names = ["both", "differ", "af_differ", "ma_samples_differ", "ma_count_differ", "only_eqtl", "only_sqtl"]
    res = dict(zip(names, r))
    res["neither"] = je["variants"]["n_cis"] - (res["both"] + res["only_eqtl"] + res["only_sqtl"])
    diff_examples = con.execute(f"""
        SELECT e.vidx, e.af AS af_e, s.af AS af_s, e.ma_samples AS ms_e, s.ma_samples AS ms_s, e.ma_count AS mc_e, s.ma_count AS mc_s
        FROM '{ve}' e JOIN '{vs}' s USING (vidx)
        WHERE (e.{af}) IS DISTINCT FROM (s.{af}) OR e.ma_samples IS DISTINCT FROM s.ma_samples OR e.ma_count IS DISTINCT FROM s.ma_count
        LIMIT 10""").fetchall()
    res["differ_examples"] = [list(x) for x in diff_examples]
    r = con.execute(f"""
        SELECT count(*) AS introns,
               count(*) FILTER (WHERE e.phenotype_id IS NULL) AS introns_gene_not_eqtl_tested,
               count(*) FILTER (WHERE e.phenotype_id IS NOT NULL AND s.var_start = e.var_start AND s.n_rows = e.n_rows) AS same_range_as_gene,
               count(*) FILTER (WHERE e.phenotype_id IS NOT NULL AND s.anchor = e.anchor) AS same_anchor_as_gene,
               count(*) FILTER (WHERE s.anchor_minus_tss = 0) AS intron_anchor_is_tss,
               count(DISTINCT s.gene_id) AS genes,
               count(DISTINCT s.gene_id) FILTER (WHERE e.phenotype_id IS NOT NULL AND s.var_start = e.var_start AND s.n_rows = e.n_rows) AS genes_with_an_intron_sharing_range
        FROM '{ps}' s LEFT JOIN '{pe}' e ON e.phenotype_id = s.gene_id""").fetchone()
    res.update(dict(zip(["introns", "introns_gene_not_eqtl_tested", "same_range_as_gene", "same_anchor_as_gene",
                         "intron_anchor_is_tss", "genes", "genes_with_an_intron_sharing_range"], r)))
    # genes where every intron shares the eQTL range: the sQTL tab then needs no second variants request
    res["genes_all_introns_share_range"] = con.execute(f"""
        SELECT count(*) FROM (SELECT s.gene_id, bool_and(e.phenotype_id IS NOT NULL AND s.var_start = e.var_start AND s.n_rows = e.n_rows) AS ok
                              FROM '{ps}' s LEFT JOIN '{pe}' e ON e.phenotype_id = s.gene_id GROUP BY 1) WHERE ok""").fetchone()[0]
    res.update(chrom=chrom, sources={"e": je["source"], "s": js["source"]}, af_compared_as="double" if both_raw else "float32")
    con.close()
    shutil.rmtree(out / "tmp" / f"cross_{chrom}", ignore_errors=True)
    write_json(out / f"cross_{chrom}.json", res)
    return res


def cmd_run(cfg: Config, args) -> None:
    out = checks_dir(cfg)
    jobs = []
    for t in args.type or ["e", "s"]:
        for c in args.chrom or CHROMS:
            r = resolve_source(cfg, t, c, args.source)
            if r is None:
                log(f"packcheck run: no {TYPES[t]['name']} source for {c} ({args.source}), skipping")
                continue
            kind, src, rows, nbytes = r
            jobs.append((t, c, kind, src, rows, nbytes))
    jobs.sort(key=lambda j: -j[5])      # biggest first
    log(f"packcheck run: {len(jobs)} jobs with {args.workers or cfg['workers']} workers -> {out}")
    failed = []
    with ProcessPoolExecutor(max_workers=args.workers or cfg["workers"]) as ex:
        futs = {ex.submit(_worker, j[:5]): j for j in jobs}
        for f in as_completed(futs):
            j = futs[f]
            try:
                r = f.result()
            except Exception as e:                                          # noqa: BLE001
                failed.append((j[0], j[1], repr(e)))
                log(f"packcheck run: FAILED {j[0]} {j[1]}: {e!r}")
                continue
            ph, st = r["phenotypes"], r["stream"]
            log(f"packcheck run: {r['type']} {r['chrom']} ({r['source']}) {r['rows']['joined']:,} rows, "
                f"{ph['n']:,} phenotypes, not contiguous {ph['not_contiguous']}, anchor not constant {ph['anchor_not_constant']}, "
                f"variants differing {r['variant_values']['any_differs']}, SE fails {st['full']['pass_domain_fail']}, {r['seconds']} s")
    for c in CHROMS:
        x = cross_type(cfg, c)
        if x:
            log(f"packcheck run: cross {c}: both {x['both']:,}, differ {x['differ']}, only e {x['only_eqtl']:,}, "
                f"only s {x['only_sqtl']:,}, neither {x['neither']:,}; introns sharing the gene range {x['same_range_as_gene']}/{x['introns']}")
    shutil.rmtree(out / "tmp", ignore_errors=True)
    if failed:
        die(f"packcheck run: {len(failed)} job(s) failed: {failed}")


# ---- check 6: round trip through the 16-bit codes ----------------------------------------------

def _rt_tally() -> dict:
    z = lambda v=0: [v] * NB                                                     # noqa: E731
    h = lambda: [[0] * (len(HIST_EDGES) + 1) for _ in range(NB)]                 # noqa: E731
    return {
        "rows": 0, "phenotypes": 0,
        "nulls": {"p_null": 0, "se_or_slope_null": 0, "p_zero": 0,
                  "p_null_mismatch": 0, "se_null_mismatch": 0, "slope_null_mismatch": 0},
        "n": z(), "n_nonfinite": z(), "n_fail_se": z(), "n_fail_slope": z(),
        "max_abs_se": z(0.0), "max_abs_slope": z(0.0), "max_slope_err_over_se": z(0.0),
        "max_se_err_over_bound": z(0.0), "max_slope_err_over_bound": z(0.0), "max_slope_bound": z(0.0),
        "max_abs_nlp": z(0.0), "hist_abs_se": h(), "hist_abs_slope": h(),
        "changed_3sf_p": z(), "changed_3sf_slope": z(), "changed_3sf_se": z(),
        "rows_p_reads_1": 0, "max_abs_slope_p_reads_1": 0.0,
        "max_nlp_err_over_halfstep": 0.0, "max_lse_err_over_halfstep": 0.0,
        "max_nlp_max": 0.0, "max_lse_width": 0.0,
        "worst_ratio": [], "worst_abs": [],
    }


def _rt_update(tal: dict, ids: pa.Array, starts: np.ndarray, ends: np.ndarray, p: np.ndarray, slope: np.ndarray,
               se: np.ndarray, dof: int, per: list | None) -> None:
    """Check 6 on a batch of whole phenotypes: encode with the phenotype's own scales, decode as a
    reader does (packfmt), and compare slope_se and the derived slope with the source against the
    SPEC per-row bounds."""
    n = len(p)
    nlp_q = np.empty(n, dtype=np.uint16)
    nlp_hat, se_hat, se_b, sl_b, h_nlp, h_lse = (np.empty(n) for _ in range(6))
    neg = np.empty(n, dtype=bool)
    for a, b in zip(starts.tolist(), ends.tolist()):
        q1, nm = packfmt.quantize_nlp(p[a:b])
        q2, lo, hi = packfmt.quantize_se(se[a:b], slope[a:b])
        nlp_q[a:b] = q1
        nlp_hat[a:b] = packfmt.dequantize_nlp(q1, nm)
        se_hat[a:b], neg[a:b] = packfmt.dequantize_se(q2, lo, hi)
        se_b[a:b], sl_b[a:b] = packfmt.error_bounds(p[a:b], slope[a:b], se[a:b], nm, lo, hi, dof)
        h_nlp[a:b] = nm / (2 * packfmt.NLP_MAXQ)
        h_lse[a:b] = (hi - lo) / (2 * packfmt.SE_MAXQ)
        tal["max_nlp_max"] = max(tal["max_nlp_max"], nm)
        tal["max_lse_width"] = max(tal["max_lse_width"], hi - lo)
    p_hat = packfmt.p_from_nlp(nlp_hat)
    slope_hat = packfmt.slope_from_se(se_hat, neg, p_hat, dof)
    tal["rows"] += n
    tal["phenotypes"] += len(starts)

    p_null = np.isnan(p)
    sv_null = np.isnan(se) | np.isnan(slope)
    N = tal["nulls"]
    N["p_null"] += int(p_null.sum())
    N["se_or_slope_null"] += int(sv_null.sum())
    N["p_zero"] += int((p == 0).sum())
    N["p_null_mismatch"] += int(((nlp_q == packfmt.NLP_NULL) != p_null).sum())
    N["se_null_mismatch"] += int((np.isnan(se_hat) != sv_null).sum())
    N["slope_null_mismatch"] += int((np.isnan(slope_hat) != (p_null | (p == 0) | sv_null)).sum())

    dom = ~np.isnan(sl_b)
    with np.errstate(invalid="ignore", divide="ignore"):
        abs_t = np.abs(slope / se)
        se_err = np.abs(se_hat - se)
        sl_err = np.abs(slope_hat - slope)
        se_ratio = se_err / se_b
        sl_ratio = sl_err / sl_b
        nlp_err = np.abs(nlp_hat + np.log10(np.where(p > 0, p, np.nan)))
        lse_err = np.abs(np.log(se_hat) - np.log(se))
    band = np.full(n, -1)
    band[dom] = _band(abs_t[dom])
    fin = np.isfinite(se_hat) & np.isfinite(slope_hat)
    F = packfmt.BOUND_FACTOR
    m = dom & (h_nlp > 0)
    if m.any():
        tal["max_nlp_err_over_halfstep"] = max(tal["max_nlp_err_over_halfstep"], float(np.max(nlp_err[m] / h_nlp[m])))
    m = dom & (h_lse > 0)
    if m.any():
        tal["max_lse_err_over_halfstep"] = max(tal["max_lse_err_over_halfstep"], float(np.max(lse_err[m] / h_lse[m])))
    ch = {}
    for key, x, y in (("p", p, p_hat), ("slope", slope, slope_hat), ("se", se, se_hat)):
        c = np.zeros(n, bool)
        c[dom] = _changed3(x[dom], y[dom])
        ch[key] = c
    r1 = dom & (nlp_q == 0)
    tal["rows_p_reads_1"] += int(r1.sum())
    if r1.any():
        tal["max_abs_slope_p_reads_1"] = max(tal["max_abs_slope_p_reads_1"], float(np.nanmax(sl_err[r1])))
    for k in range(NB):
        mk = dom & (band == k)
        if not mk.any():
            continue
        tal["n"][k] += int(mk.sum())
        tal["n_nonfinite"][k] += int((mk & ~fin).sum())
        tal["n_fail_se"][k] += int((mk & ~(se_err <= F * se_b)).sum())       # NaN errors fail
        tal["n_fail_slope"][k] += int((mk & ~(sl_err <= F * sl_b)).sum())
        mf = mk & fin
        if mf.any():
            tal["max_abs_se"][k] = max(tal["max_abs_se"][k], float(se_err[mf].max()))
            tal["max_abs_slope"][k] = max(tal["max_abs_slope"][k], float(sl_err[mf].max()))
            tal["max_slope_err_over_se"][k] = max(tal["max_slope_err_over_se"][k], float((sl_err[mf] / se[mf]).max()))
            tal["max_se_err_over_bound"][k] = max(tal["max_se_err_over_bound"][k], float(se_ratio[mf].max()))
            tal["max_slope_err_over_bound"][k] = max(tal["max_slope_err_over_bound"][k], float(sl_ratio[mf].max()))
            tal["max_slope_bound"][k] = max(tal["max_slope_bound"][k], float(sl_b[mf].max()))
            mp = mf & (p > 0)
            if mp.any():
                tal["max_abs_nlp"][k] = max(tal["max_abs_nlp"][k], float(nlp_err[mp].max()))
            for key, err in (("hist_abs_se", se_err), ("hist_abs_slope", sl_err)):
                hh = np.bincount(np.searchsorted(HIST_EDGES, err[mf], side="right"), minlength=len(HIST_EDGES) + 1)
                tal[key][k] = (np.array(tal[key][k]) + hh).tolist()
        for key in ("p", "slope", "se"):
            tal[f"changed_3sf_{key}"][k] += int(ch[key][mk].sum())

    idx = np.flatnonzero(dom & fin)
    if idx.size:
        def row(i):
            return {"phenotype_id": ids[int(i)].as_py(), "pval_nominal": p[i], "slope": slope[i], "slope_se": se[i],
                    "abs_t": abs_t[i], "nlp_max": h_nlp[i] * 2 * packfmt.NLP_MAXQ, "lse_width": h_lse[i] * 2 * packfmt.SE_MAXQ,
                    "se_hat": se_hat[i], "slope_hat": slope_hat[i], "se_err": se_err[i], "slope_err": sl_err[i],
                    "se_bound": se_b[i], "slope_bound": sl_b[i], "slope_ratio": sl_ratio[i]}
        for key, score in (("worst_ratio", sl_ratio), ("worst_abs", sl_err)):
            top = idx[np.argsort(-score[idx])[:10]]
            tal[key] = sorted(tal[key] + [row(i) for i in top], key=lambda r: -_num(r["slope_ratio" if key == "worst_ratio" else "slope_err"]))[:10]
    if per is not None:
        for a, b in zip(starts.tolist(), ends.tolist()):
            d = dom[a:b] & fin[a:b]
            if not d.any():
                continue
            per.append({"phenotype_id": ids[a].as_py(), "rows": b - a, "nlp_max": h_nlp[a] * 2 * packfmt.NLP_MAXQ,
                        "lse_width": h_lse[a] * 2 * packfmt.SE_MAXQ, "max_abs_se": se_err[a:b][d].max(),
                        "max_abs_slope": sl_err[a:b][d].max(), "max_slope_err_over_se": (sl_err[a:b][d] / se[a:b][d]).max(),
                        "max_se_err_over_bound": se_ratio[a:b][d].max(), "max_slope_err_over_bound": sl_ratio[a:b][d].max()})


def _rt_feed(tal: dict, tb: pa.Table, dof: int, per: list | None) -> None:
    if tb.num_rows == 0:
        return
    tb = tb.combine_chunks()
    ids = tb["phenotype_id"].combine_chunks()
    starts, ends = _slices(ids)
    col = lambda c: tb[c].to_numpy(zero_copy_only=False).astype(np.float64)          # noqa: E731
    _rt_update(tal, ids, starts, ends, col("pval_nominal"), col("slope"), col("slope_se"), dof, per)


def _rt_job(job) -> dict:
    t, chrom, src, ids = job
    cfg = Config()
    dof = int(cfg["packs"]["dof"][TYPES[t]["name"]])
    tal = _rt_tally()
    per = None if ids is None else []
    t0 = time.time()
    cols = ["phenotype_id", "pval_nominal", "slope", "slope_se"]
    if ids is None:
        for tb in phenotype_batches(pq.ParquetFile(src), cols):
            _rt_feed(tal, tb, dof, None)
    else:
        work = checks_dir(cfg) / "tmp" / f"rt_{t}_{chrom}"
        con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"], temp_dir=work)
        lst = ", ".join("'" + i.replace("'", "''") + "'" for i in ids)
        tb = con.execute(f"""SELECT phenotype_id, pval_nominal::DOUBLE AS pval_nominal, slope::DOUBLE AS slope, slope_se::DOUBLE AS slope_se
                             FROM read_parquet('{src}', file_row_number = true) WHERE phenotype_id IN ({lst})
                             ORDER BY phenotype_id, file_row_number""").fetch_arrow_table()
        con.close()
        shutil.rmtree(work, ignore_errors=True)
        _rt_feed(tal, tb, dof, per)
    return {"type": t, "chrom": chrom, "src": src, "dof": dof, "tally": tal, "phenotypes": per,
            "seconds": round(time.time() - t0, 1)}


def _merge_into(dst: dict, src: dict) -> None:
    """Sum counts, take maxima, add histograms; keys the destination does not know and `worst*`
    lists are skipped (the caller merges those)."""
    for k, v in src.items():
        if k not in dst or k.startswith("worst"):
            continue
        if isinstance(v, dict):
            _merge_into(dst[k], v)
        elif isinstance(v, list) and v and isinstance(v[0], list):
            dst[k] = [(np.array(a) + np.array(b)).tolist() for a, b in zip(dst[k], v)]
        elif isinstance(v, list):
            dst[k] = [max(_num(a), _num(b)) for a, b in zip(dst[k], v)] if "max" in k else [a + b for a, b in zip(dst[k], v)]
        elif "max" in k:
            dst[k] = max(_num(dst[k]), _num(v))
        else:
            dst[k] += v


def _rt_merge(ts: list[dict]) -> dict:
    out = _rt_tally()
    for t in ts:
        _merge_into(out, t)
        out["worst_ratio"] = sorted(out["worst_ratio"] + t["worst_ratio"], key=lambda r: -_num(r["slope_ratio"]))[:10]
        out["worst_abs"] = sorted(out["worst_abs"] + t["worst_abs"], key=lambda r: -_num(r["slope_err"]))[:10]
    return out


def _rt_overall(T: dict) -> dict:
    """All |t| bands together."""
    o = {k: sum(T[k]) for k in ("n", "n_nonfinite", "n_fail_se", "n_fail_slope", "changed_3sf_p", "changed_3sf_slope", "changed_3sf_se")}
    for k in ("max_abs_se", "max_abs_slope", "max_slope_err_over_se", "max_se_err_over_bound", "max_slope_err_over_bound",
              "max_slope_bound", "max_abs_nlp"):
        o[k] = max(_num(x) for x in T[k])
    for k in ("hist_abs_se", "hist_abs_slope"):
        o[k] = np.sum(np.array(T[k]), axis=0).tolist()
    return o


def _rt_status(T: dict) -> bool:
    o, N = _rt_overall(T), T["nulls"]
    return (o["n_fail_se"] == 0 and o["n_fail_slope"] == 0 and o["n_nonfinite"] == 0
            and N["p_null_mismatch"] == 0 and N["se_null_mismatch"] == 0 and N["slope_null_mismatch"] == 0)


def cmd_roundtrip(cfg: Config, args) -> None:
    """Check 6 alone, streaming the raw nominal files: all rows of the chosen chromosomes, or with
    --steepest N the N phenotypes per type with the largest -log10 p (the largest p steps)."""
    out = checks_dir(cfg)
    types = args.type or ["e", "s"]
    jobs, scope = [], {}
    if args.steepest:
        label = f"steepest{args.steepest}"
        con = connect(cfg, temp_dir=out / "tmp" / "rt_select")
        for t in types:
            ph = out / f"phenotypes_{t}.parquet"
            if not ph.exists():
                die(f"packcheck roundtrip: {ph} is missing; run `packcheck run` and `packcheck report` first")
            rows = con.execute(f"""SELECT chr, phenotype_id, max_nlp FROM '{ph}' WHERE n_rows > 0
                                   ORDER BY max_nlp DESC, phenotype_id LIMIT {int(args.steepest)}""").fetchall()
            by_chr: dict[str, list[str]] = {}
            for c, pid, _ in rows:
                by_chr.setdefault(c, []).append(pid)
            scope[TYPES[t]["name"]] = (f"the {len(rows)} phenotypes with the largest -log10 p genome-wide "
                                       f"(max -log10 p {_f(rows[0][2])} to {_f(rows[-1][2])})")
            for c, pids in by_chr.items():
                r = resolve_source(cfg, t, c, "raw")
                if r is None:
                    die(f"packcheck roundtrip: no raw {TYPES[t]['name']} file for {c}")
                jobs.append((t, c, r[1], pids, len(pids)))
        con.close()
    else:
        chroms = args.chrom or CHROMS
        label = "genome" if not args.chrom else "_".join(chroms)
        for t in types:
            scope[TYPES[t]["name"]] = f"every row of {', '.join(chroms) if args.chrom else 'every chromosome'}"
            for c in chroms:
                r = resolve_source(cfg, t, c, "raw")
                if r is None:
                    die(f"packcheck roundtrip: no raw {TYPES[t]['name']} file for {c}")
                jobs.append((t, c, r[1], None, r[2]))
    jobs.sort(key=lambda j: -j[4])
    log(f"packcheck roundtrip ({label}): {len(jobs)} jobs with {args.workers or cfg['workers']} workers")
    results: dict[str, list[dict]] = {t: [] for t in types}
    with ProcessPoolExecutor(max_workers=args.workers or cfg["workers"]) as ex:
        futs = {ex.submit(_rt_job, j[:4]): j for j in jobs}
        for f in as_completed(futs):
            r = f.result()
            results[r["type"]].append(r)
            o = _rt_overall(r["tally"])
            log(f"packcheck roundtrip: {r['type']} {r['chrom']}: {r['tally']['rows']:,} rows, max abs err SE {_f(o['max_abs_se'])}, "
                f"slope {_f(o['max_abs_slope'])}; max err / bound SE {_f(o['max_se_err_over_bound'])}, slope {_f(o['max_slope_err_over_bound'])}; "
                f"{r['seconds']} s")
    R = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "label": label,
         "bound_factor": packfmt.BOUND_FACTOR, "rel_f32": packfmt.REL_F32, "types": {}}
    for t in types:
        rs = results[t]
        if not rs:
            continue
        per = sorted((x for r in rs for x in (r["phenotypes"] or [])), key=lambda x: -_num(x["nlp_max"]))
        R["types"][TYPES[t]["name"]] = {
            "scope": scope[TYPES[t]["name"]], "dof": rs[0]["dof"], "chroms": sorted({r["chrom"] for r in rs}, key=CHROMS.index),
            "seconds": sum(r["seconds"] for r in rs), "tally": _rt_merge([r["tally"] for r in rs]), "phenotypes": per}
    shutil.rmtree(out / "tmp", ignore_errors=True)
    write_json(out / f"roundtrip_{label}.json", R)
    (out / "roundtrip.md").write_text(_roundtrip_md(out))
    for name, X in R["types"].items():
        o = _rt_overall(X["tally"])
        log(f"packcheck roundtrip ({label}) {name}: {'PASS' if _rt_status(X['tally']) else 'FAIL'}; {X['tally']['rows']:,} rows; "
            f"max abs err SE {_f(o['max_abs_se'])}, slope {_f(o['max_abs_slope'])} ({_f(o['max_slope_err_over_se'])} of SE); "
            f"max err / bound SE {_f(o['max_se_err_over_bound'])}, slope {_f(o['max_slope_err_over_bound'])}")
    log(f"packcheck roundtrip: wrote {out / f'roundtrip_{label}.json'} and {out / 'roundtrip.md'}")


def _roundtrip_md(out: Path) -> str:
    L = ["# packcheck roundtrip (check 6)", "",
         "Each phenotype's rows are encoded with its own scales (`packfmt.quantize_nlp`, `quantize_se`) and decoded as a reader "
         "does: slope_se = exp(lse_min + q * step), slope = sign * slope_se * t with t = max(0, -stdtrit(dof, p / 2)). "
         "Errors are absolute, against the raw Zenodo values. The per-row bound is `packfmt.error_bounds` (SPEC section 9); "
         "a row fails when its error exceeds the bound factor times its bound, or a value that should exist is not finite. "
         "Bands are by |t| = |slope| / slope_se from the source.", ""]
    for path in sorted(out.glob("roundtrip_*.json")):
        R = json.loads(path.read_text())
        L += [f"## {R['label']}", "", f"Generated {R['generated']}; bound factor {R['bound_factor']}; float32 term {R['rel_f32']:g} x |slope|.", ""]
        for name, X in R["types"].items():
            T = X["tally"]
            o = _rt_overall(T)
            bands = [(BAND_NAMES[k], {kk: (v[k] if isinstance(v, list) else v) for kk, v in T.items() if isinstance(v, list) and len(v) == NB})
                     for k in range(NB)] + [("all", o)]
            L += [f"### {name}: {'PASS' if _rt_status(T) else 'FAIL'}", "",
                  f"{X['scope']} ({', '.join(X['chroms'])}); dof {X['dof']}; {T['rows']:,} rows in {T['phenotypes']:,} phenotypes; {X['seconds'] / 60:.1f} worker-minutes.", "",
                  _table(["|t| band", "rows", "not finite", "max abs SE err", "p99.9 SE err", "max abs slope err", "p99.9 slope err",
                          "max slope err / SE", "max SE err / bound", "max slope err / bound", "fail SE", "fail slope",
                          "max abs -log10 p err", "p changed at 3 s.f.", "slope changed", "SE changed"],
                         [[nm, b["n"], b["n_nonfinite"], b["max_abs_se"], _hist_quantile(b["hist_abs_se"], 0, 0.999), b["max_abs_slope"],
                           _hist_quantile(b["hist_abs_slope"], 0, 0.999), b["max_slope_err_over_se"], b["max_se_err_over_bound"],
                           b["max_slope_err_over_bound"], b["n_fail_se"], b["n_fail_slope"], b["max_abs_nlp"],
                           b["changed_3sf_p"], b["changed_3sf_slope"], b["changed_3sf_se"]] for nm, b in bands]), "",
                  f"Nulls: p null {T['nulls']['p_null']:,}, slope or SE null {T['nulls']['se_or_slope_null']:,}, p = 0 {T['nulls']['p_zero']:,}; "
                  f"decoded null pattern mismatches: p {T['nulls']['p_null_mismatch']}, SE {T['nulls']['se_null_mismatch']}, slope {T['nulls']['slope_null_mismatch']}. "
                  f"Rows whose p reads back as exactly 1: {T['rows_p_reads_1']:,} (max abs slope err {_f(T['max_abs_slope_p_reads_1'])}). "
                  f"-log10 p error / half step max {_f(T['max_nlp_err_over_halfstep'])}; log(SE) error / half step max {_f(T['max_lse_err_over_halfstep'])}. "
                  f"Largest nlp_max {_f(T['max_nlp_max'])}; widest log(SE) range {_f(T['max_lse_width'])}. "
                  "The p99.9 columns are upper bin edges of a 10-bins-per-decade histogram.", ""]
            for key, title in (("worst_ratio", "10 rows with the largest slope error / bound"), ("worst_abs", "10 rows with the largest absolute slope error")):
                L += [f"{title}:", "", _table(["phenotype", "p", "slope", "slope_se", "|t|", "nlp_max", "slope err", "slope bound", "err / bound", "SE err", "SE bound"],
                                              [[w["phenotype_id"], _num(w["pval_nominal"]), _num(w["slope"]), _num(w["slope_se"]), _num(w["abs_t"]), _num(w["nlp_max"]),
                                                _num(w["slope_err"]), _num(w["slope_bound"]), _num(w["slope_ratio"]), _num(w["se_err"]), _num(w["se_bound"])]
                                               for w in T[key]]), ""]
            if X["phenotypes"]:
                L += ["Per phenotype:", "", _table(["phenotype", "rows", "nlp_max", "log(SE) range", "max abs SE err", "max abs slope err",
                                                    "max slope err / SE", "max SE err / bound", "max slope err / bound"],
                                                   [[x["phenotype_id"], x["rows"], _num(x["nlp_max"]), _num(x["lse_width"]), _num(x["max_abs_se"]),
                                                     _num(x["max_abs_slope"]), _num(x["max_slope_err_over_se"]), _num(x["max_se_err_over_bound"]),
                                                     _num(x["max_slope_err_over_bound"])] for x in X["phenotypes"]]), ""]
    return "\n".join(L) + "\n"


# ---- Step 5: report ---------------------------------------------------------------------------

def _merge_tallies(ts: list[dict]) -> dict:
    """Sum counts, take maxima, add histograms, keep the 20 worst rows."""
    out = _new_tally()
    for t in ts:
        _merge_into(out, t)
        out["worst_full"] = sorted(out["worst_full"] + t["worst_full"], key=lambda r: -_num(r["rel_err"]))[:20]
    return out


def _hist_quantile(h: list[int], nonfinite: int, q: float) -> float | None:
    """Upper bin edge of the q-quantile of the absolute-error histogram (non-finite counted as +inf)."""
    total = sum(h) + nonfinite
    if total == 0:
        return None
    rank, cum = q * total, 0
    for k, c in enumerate(h):
        cum += c
        if cum >= rank:
            return float(HIST_EDGES[k]) if k < len(HIST_EDGES) else float("inf")
    return float("inf")


def _f(x, digits=3) -> str:
    if x is None:
        return "–"
    if isinstance(x, str):
        return x
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}"
    x = float(x)
    if not math.isfinite(x):
        return "inf" if x > 0 else ("-inf" if x < 0 else "nan")
    if x == 0:
        return "0"
    if abs(x) >= 1e5 or abs(x) < 1e-3:
        return f"{x:.{digits - 1}e}"
    return f"{x:.{digits}g}"


def _table(head: list[str], rows: list[list]) -> str:
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    lines += ["| " + " | ".join(_f(c) if not isinstance(c, str) else c for c in r) + " |" for r in rows]
    return "\n".join(lines)


def _md5s(cfg: Config, out: Path, compute: bool) -> dict:
    cache_path = out / "archive_md5.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    state = archive_state(cfg)
    for name, st in state.items():
        p = cfg.zenodo / name
        c = cache.get(name)
        if st["complete"] and compute and not (c and c["size"] == st["size"] and c["mtime"] == p.stat().st_mtime):
            h = hashlib.md5()
            with open(p, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 24), b""):
                    h.update(chunk)
            cache[name] = c = {"size": st["size"], "mtime": p.stat().st_mtime, "md5": h.hexdigest()}
            log(f"packcheck report: md5 {name} {c['md5']}")
        if c and st["complete"] and c["size"] == st["size"] and c["mtime"] == p.stat().st_mtime:
            st["md5"] = c["md5"]
            st["md5_ok"] = c["md5"] == st["md5_expected"]
        else:
            st["md5"], st["md5_ok"] = None, None
    write_json(cache_path, cache)
    return state


def cmd_report(cfg: Config, args) -> None:
    out = checks_dir(cfg)
    runs = {t: {} for t in TYPES}
    for t in TYPES:
        for c in CHROMS:
            p = out / f"{t}_{c}.json"
            if p.exists():
                runs[t][c] = json.loads(p.read_text())
    cross = {c: json.loads((out / f"cross_{c}.json").read_text()) for c in CHROMS if (out / f"cross_{c}.json").exists()}
    dof = json.loads((out / "dof.json").read_text()) if (out / "dof.json").exists() else None
    manifest = json.loads((cfg.derived / "manifest.json").read_text())
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip() or None
    except FileNotFoundError:
        head = None
    archives = _md5s(cfg, out, args.md5)
    con = connect(cfg, temp_dir=out / "tmp" / "report")

    # concatenated per-type outputs
    for t in TYPES:
        if not runs[t]:
            continue
        con.execute(f"""COPY (SELECT * FROM read_parquet('{out}/phenotypes_{t}_chr*.parquet', union_by_name = true)
                     ORDER BY chr, phenotype_id) TO '{out}/phenotypes_{t}.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)""")
        con.execute(f"""COPY (SELECT regexp_extract(filename, 'variant_values_{t}_(chr[0-9XY]+)\\.parquet', 1) AS chr, * EXCLUDE (filename)
                     FROM read_parquet('{out}/variant_values_{t}_chr*.parquet', filename = true)
                     ORDER BY chr, vidx) TO '{out}/variant_values_{t}.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)""")

    R: dict = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "repo_head": head,
               "inputs": {"manifest_built": manifest["built"], "manifest_pipeline_commit": manifest["pipeline_commit"],
                          "archives": archives},
               "types": {}, "cross": {}, "checks": []}
    derived_rows = {}
    for t, T in TYPES.items():
        derived_rows[t] = {c: sum(pq.read_metadata(f).num_rows for f in (cfg.tables / T["derived"] / f"chr={c}").glob("bin=*/data.parquet"))
                           for c in CHROMS}
    for t, T in TYPES.items():
        rs = runs[t]
        if not rs:
            continue
        name = T["name"]
        ph = {k: sum((r["phenotypes"][k] or 0) for r in rs.values())
              for k in rs[next(iter(rs))]["phenotypes"] if k not in ("min_p", "max_nlp", "max_abs_slope", "max_abs_anchor_minus_tss", "rows_minus_perm_num_var_hist")}
        ph["min_p"] = min(_num(r["phenotypes"]["min_p"]) for r in rs.values())
        ph["max_nlp"] = max(_num(r["phenotypes"]["max_nlp"]) for r in rs.values())
        ph["max_abs_slope"] = max(_num(r["phenotypes"]["max_abs_slope"]) for r in rs.values())
        ph["max_abs_anchor_minus_tss"] = max((_num(r["phenotypes"]["max_abs_anchor_minus_tss"]) for r in rs.values()), default=0)
        hist: dict[str, int] = {}
        for r in rs.values():
            for d, n in r["phenotypes"]["rows_minus_perm_num_var_hist"].items():
                hist[d] = hist.get(d, 0) + n
        ph["rows_minus_perm_num_var_hist"] = hist
        vv = {k: sum(r["variant_values"][k] for r in rs.values()) for k in ("tested", "af_differs", "ma_samples_differs", "ma_count_differs", "any_differs", "af_over_0_5")}
        vv["af_min"] = min(r["variant_values"]["af_min"] for r in rs.values())
        vv["af_max"] = max(r["variant_values"]["af_max"] for r in rs.values())
        vv["ma_samples_max"] = max(r["variant_values"]["ma_samples_max"] for r in rs.values())
        vv["ma_count_max"] = max(r["variant_values"]["ma_count_max"] for r in rs.values())
        variants = {"n_cis": sum(r["variants"]["n_cis"] for r in rs.values()),
                    "duplicate_keys": sum(r["variants"]["duplicate_keys"] for r in rs.values()),
                    "chroms_file_order_not_vidx": [c for c, r in rs.items() if not r["variants"]["file_order_equals_vidx"]]}
        st = _merge_tallies([r["stream"] for r in rs.values()])
        q = f"read_parquet('{out}/phenotypes_{t}.parquet')"
        dist = {}
        for col in ("max_nlp", "max_abs_slope", "max_abs_slope / nullif(median_abs_slope, 0)"):
            mn, med, p99, mx = con.execute(f"SELECT min({col}), median({col}), quantile_cont({col}, 0.99), max({col}) FROM {q} WHERE n_rows > 0").fetchone()
            dist[col.split(" /")[0] if "/" not in col else "max_over_median_abs_slope"] = {"min": mn, "median": med, "p99": p99, "max": mx}
        a = con.execute(f"""SELECT count(*), count(*) FILTER (WHERE anchor_minus_tss <> 0),
                            quantile_cont(abs(anchor_minus_tss), 0.5) FILTER (WHERE anchor_minus_tss <> 0),
                            quantile_cont(abs(anchor_minus_tss), 0.9) FILTER (WHERE anchor_minus_tss <> 0),
                            max(abs(anchor_minus_tss)), count(*) FILTER (WHERE tss IS NULL) FROM {q}""").fetchone()
        dist["anchor_minus_tss"] = dict(zip(["phenotypes", "nonzero", "median_abs_nonzero", "p90_abs_nonzero", "max_abs", "no_tss"], a))
        if t == "e":
            b = con.execute(f"""SELECT count(*), count(*) FILTER (WHERE p.anchor_minus_tss <> 0), max(abs(p.anchor_minus_tss))
                               FROM {q} p JOIN '{cfg.tables / 'genes.parquet'}' g ON g.gene_id = p.phenotype_id
                               WHERE g.chr = 'chr7' AND g.bin = 7""").fetchone()
            dist["anchor_minus_tss"]["chr7_bin7"] = dict(zip(["genes", "nonzero", "max_abs"], b))
        coverage = {"chroms": sorted(rs, key=CHROMS.index), "missing": [c for c in CHROMS if c not in rs],
                    "sources": {k: sorted(c for c, r in rs.items() if r["source"] == k) for k in ("raw", "derived")},
                    "rows_by_chrom": {c: {"checked": r["rows"]["joined"], "source": r["source"], "derived_table": derived_rows[t][c]} for c, r in rs.items()}}
        R["types"][name] = {"coverage": coverage, "phenotypes": ph, "variant_values": vv, "variants": variants,
                            "stream": st, "distributions": dist,
                            "dof": (dof or {}).get("types", {}).get(name), "config_dof": cfg["packs"]["dof"][name],
                            "rows_checked": sum(r["rows"]["joined"] for r in rs.values()),
                            # 9cde62f dropped the manifest's `tables` block with the parquet tables it
                            # described; the build's own row tally now lives with the pack it wrote
                            "manifest_rows": manifest["packs"]["counts"][name]["rows"],
                            "seconds": sum(r["seconds"] for r in rs.values())}
    if cross:
        keys = ["both", "differ", "af_differ", "ma_samples_differ", "ma_count_differ", "only_eqtl", "only_sqtl", "neither",
                "introns", "introns_gene_not_eqtl_tested", "same_range_as_gene", "same_anchor_as_gene", "intron_anchor_is_tss",
                "genes", "genes_with_an_intron_sharing_range", "genes_all_introns_share_range"]
        R["cross"] = {k: sum(x[k] for x in cross.values()) for k in keys}
        R["cross"]["chroms"] = sorted(cross, key=CHROMS.index)
        R["cross"]["sources"] = sorted({f"e:{x['sources']['e']} s:{x['sources']['s']}" for x in cross.values()})
        R["cross"]["af_compared_as"] = sorted({x["af_compared_as"] for x in cross.values()})
        R["cross"]["differ_examples"] = [e for x in cross.values() for e in x["differ_examples"]][:10]

    # ---- checks
    E, S = R["types"].get("eqtl"), R["types"].get("sqtl")

    def scope(X):
        if X is None:
            return "not run"
        cov = X["coverage"]
        s = f"{len(cov['chroms'])}/{len(CHROMS)} chr"
        if cov["sources"]["derived"]:
            s += ", derived" + (" (significant introns only)" if X is S else "") + (" on " + ",".join(cov["sources"]["derived"]) if cov["sources"]["raw"] else "")
        else:
            s += ", raw"
        return s

    def add(num, check, status, counts, fallback, per):
        R["checks"].append({"num": num, "check": check, "status": status, "counts": counts, "fallback": fallback, "scope": per})

    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        p = X["phenotypes"]
        ok1 = p["not_contiguous"] == 0
        add("1", f"Contiguous run of the variant list ({nm})", "PASS" if ok1 else "FAIL",
            f"{p['not_contiguous']:,} of {p['n']:,} phenotypes not contiguous: {p['with_unmatched_rows']:,} with unmatched rows "
            f"({p['unmatched_rows']:,} rows), {p['with_duplicate_vidx']:,} with a repeated variant, {p['with_gaps']:,} with gaps ({p['gap_variants']:,} gap variants)",
            "Unmatched rows: fix variants_collect. Gaps: presence bitmap over [var_start, var_start + span).", scope(X))
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        p = X["phenotypes"]
        add("1b", f"Raw rows grouped and in variant order ({nm})", "INFO",
            f"{p['not_grouped']:,} phenotypes not grouped in the source; {X['stream']['phenotypes_not_in_variant_order']:,} not in variant order",
            "None: the builder sorts by vidx.", scope(X))
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        p = X["phenotypes"]
        add("1c", f"n_rows versus permutation num_var ({nm})", "INFO",
            f"{p['rows_ne_perm_num_var']:,} of {p['n']:,} differ (n_rows - num_var: {p['rows_minus_perm_num_var_hist']}); {p['no_permutation_row']:,} without a permutation row",
            "None: n_var comes from the nominal row count.", scope(X))
    dup = (E or S)["variants"]["duplicate_keys"] if (E or S) else None
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        v = X["variant_values"]
        ok = v["any_differs"] == 0 and X["variants"]["duplicate_keys"] == 0
        add("2", f"af, ma_samples, ma_count constant per variant within {nm}; variant keys unique", "PASS" if ok else "FAIL",
            f"{v['any_differs']:,} of {v['tested']:,} tested variants differ (af {v['af_differs']:,}, ma_samples {v['ma_samples_differs']:,}, "
            f"ma_count {v['ma_count_differs']:,}); duplicate keys {X['variants']['duplicate_keys']:,}",
            "Per-row af/ma_samples/ma_count arrays in the block. Duplicate keys: fix variants_rsid.", scope(X))
    if R["cross"]:
        x = R["cross"]
        add("2", "af, ma_samples, ma_count equal across eQTL and sQTL", "PASS" if x["differ"] == 0 else "FAIL",
            f"{x['differ']:,} of {x['both']:,} variants tested by both differ (af {x['af_differ']:,}, ma_samples {x['ma_samples_differ']:,}, "
            f"ma_count {x['ma_count_differ']:,}; af compared as {','.join(x['af_compared_as'])}); only eQTL {x['only_eqtl']:,}, only sQTL {x['only_sqtl']:,}, neither {x['neither']:,}",
            "Variants file keeps the eQTL values; sQTL blocks carry per-row arrays.",
            f"{len(x['chroms'])} chr, {'; '.join(x['sources'])}")
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        F, C = X["stream"]["full"], X["stream"]["cats"]
        d = X["dof"] or {}
        ok = bool(d.get("unanimous")) and F["pass_domain_fail"] == 0 and d.get("mode") == X["config_dof"]
        add("3", f"One dof per type rebuilds slope_se at full precision ({nm}, dof {X['config_dof']})", "PASS" if ok else "FAIL",
            f"dof search: mode {d.get('mode')}, unanimous {d.get('unanimous')} over {d.get('sampled_phenotypes')} sampled phenotypes; "
            f"{F['pass_domain_fail']:,} of {F['pass_domain_rows']:,} rows in the pass domain above rel. error {REL_TOL:g} (max {_f(F['max_rel_pass_domain'])}). "
            f"Counted, not failed: p null {C['p_null']:,}, p = 0 {C['p_zero']:,}, p = 1 {C['p_one']:,}, |t| < 0.01 {F['n'][0]:,}",
            "Several dof: per-block dof. Relation wrong: store slope per row.", scope(X))
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        p = X["phenotypes"]
        add("4", f"position - tss_distance constant per phenotype ({nm})", "PASS" if p["anchor_not_constant"] == 0 else "FAIL",
            f"{p['anchor_not_constant']:,} of {p['n']:,} phenotypes not constant",
            "Per-row i32 tss_distance in the block.", scope(X))
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        a = X["distributions"]["anchor_minus_tss"]
        extra = ""
        if X is E and "chr7_bin7" in a:
            extra = f"; chr7 bin 7: {a['chr7_bin7']['nonzero']} of {a['chr7_bin7']['genes']} (max {_f(a['chr7_bin7']['max_abs'])} bp)"
        if X is S and R["cross"]:
            x = R["cross"]
            extra = (f"; introns with their gene's eQTL (var_start, n_rows): {x['same_range_as_gene']:,} of {x['introns']:,}; "
                     f"same window start as the gene: {x['same_anchor_as_gene']:,}; gene not eQTL-tested: {x['introns_gene_not_eqtl_tested']:,}; "
                     f"genes whose every intron shares the range: {x['genes_all_introns_share_range']:,} of {x['genes']:,}")
        add("4b", f"Window start versus GENCODE TSS ({nm})", "INFO",
            f"{a['nonzero']:,} of {a['phenotypes']:,} phenotypes start away from the TSS (median {_f(a['median_abs_nonzero'])} bp, max {_f(a['max_abs'])} bp){extra}",
            "None in v0: block header anchor; intron blocks carry their own run (SPEC, sQTL results pack).", scope(X))
    for X, nm in ((E, "eQTL"), (S, "sQTL")):
        if X is None:
            continue
        p = X["phenotypes"]
        sub = p["min_p"] < 2.2250738585072014e-308
        add("5", f"Null, NaN, zero, one counts; ranges ({nm})", "INFO",
            f"p NaN/null {p['rows_p_nan']:,} ({p['phenotypes_with_p_nan']:,} phenotypes), p = 0 {p['rows_p_zero']:,}, p = 1 {p['rows_p_one']:,}, p > 1 {p['rows_p_gt1']:,}, "
            f"slope NaN {p['rows_slope_nan']:,}, slope = 0 {p['rows_slope_zero']:,}, se NaN {p['rows_se_nan']:,}, af NaN {p['rows_af_nan']:,}, counts null {p['rows_counts_null']:,}; "
            f"min p {_f(p['min_p'])}{' (SUBNORMAL)' if sub else ''}, max -log10 p {_f(p['max_nlp'])}, phenotypes over 65.533 {p['phenotypes_max_nlp_over_65_533']:,}, "
            f"max |slope| {_f(p['max_abs_slope'])}; af {_f(X['variant_values']['af_min'])} to {_f(X['variant_values']['af_max'])}, max ma_samples {X['variant_values']['ma_samples_max']}, max ma_count {X['variant_values']['ma_count_max']}",
            "Null codes in SPEC.md cover these.", scope(X))
    for path in sorted(out.glob("roundtrip_*.json")):
        RT = json.loads(path.read_text())
        for name, X in RT["types"].items():
            T, o, N = X["tally"], _rt_overall(X["tally"]), X["tally"]["nulls"]
            add("6", f"slope_se and derived slope through the 16-bit codes, within {RT['bound_factor']} x the per-row bound ({'eQTL' if name == 'eqtl' else 'sQTL'})",
                "PASS" if _rt_status(T) else "FAIL",
                f"{T['rows']:,} rows; max abs error SE {_f(o['max_abs_se'])}, slope {_f(o['max_abs_slope'])} ({_f(o['max_slope_err_over_se'])} of SE); "
                "slope error by |t| band " + ", ".join(f"{BAND_NAMES[k]}: {_f(T['max_abs_slope'][k])}" for k in range(NB))
                + f"; max error / bound SE {_f(o['max_se_err_over_bound'])}, slope {_f(o['max_slope_err_over_bound'])}; "
                f"rows over the limit {o['n_fail_se'] + o['n_fail_slope']:,}, not finite {o['n_nonfinite']:,}, "
                f"null mismatches {N['p_null_mismatch'] + N['se_null_mismatch'] + N['slope_null_mismatch']}; "
                f"p reading back as 1: {T['rows_p_reads_1']:,} rows (max slope error {_f(T['max_abs_slope_p_reads_1'])})",
                "Finer codes or a larger bound factor (SPEC section 9).", f"{RT['label']}: {X['scope']}")
    con.close()
    shutil.rmtree(out / "tmp", ignore_errors=True)
    write_json(out / "report.json", R)
    (out / "report.md").write_text(_report_md(R))
    for c in R["checks"]:
        log(f"packcheck report: {c['num']:>2} {c['status']:4s} {c['check']}")
    log(f"packcheck report: wrote {out / 'report.md'}")


def _report_md(R: dict) -> str:
    L = [f"# packcheck report", "", f"Generated {R['generated']}, repo HEAD `{R['repo_head']}`. "
         f"Derived inputs: manifest built {R['inputs']['manifest_built']}, pipeline commit `{R['inputs']['manifest_pipeline_commit']}`.", ""]
    L += ["## Checks", "", _table(["#", "Check", "Result", "Counts", "Fallback if failed", "Scope"],
                                  [[c["num"], c["check"], c["status"], c["counts"], c["fallback"], c["scope"]] for c in R["checks"]]), ""]
    L += ["## Inputs", "", _table(["archive", "expected bytes", "on disk", "complete", "md5 matches", "extracted"],
                                  [[k, v["expected_size"], v["size"], v["complete"], "not computed" if v["md5_ok"] is None else v["md5_ok"], v["extracted"]]
                                   for k, v in R["inputs"]["archives"].items()]), ""]
    for name, X in R["types"].items():
        cov = X["coverage"]
        L += [f"## {name}", "",
              f"Rows checked {X['rows_checked']:,} (manifest `packs.counts.{name}.rows` {X['manifest_rows']:,}); "
              f"chromosomes {len(cov['chroms'])}, raw {len(cov['sources']['raw'])}, derived {len(cov['sources']['derived'])}"
              + (f", missing {cov['missing']}" if cov["missing"] else "") + f"; {X['seconds'] / 60:.1f} worker-minutes.", ""]
        L += [_table(["chr", "source", "rows checked", "derived table rows"],
                     [[c, v["source"], v["checked"], v["derived_table"]] for c, v in cov["rows_by_chrom"].items()]), ""]
        d = X["dof"]
        if d:
            L += [f"dof: histogram {d['hist']}, mode {d['mode']}, unanimous {d['unanimous']}, sampled phenotypes {d['sampled_phenotypes']} "
                  f"(sources {d['sources']}); config {X['config_dof']}. {d['pcs']} expression/splicing PCs, so samples minus other covariates = "
                  f"dof + 2 + PCs = {d['implied_samples_minus_other_covariates']}.", ""]
        F = X["stream"]["full"]
        L += ["Full-precision SE rebuild (source p and slope):", "",
              _table(["|t| band", "rows", "max abs error", "max rel error", f"rows rel > {REL_TOL:g}"],
                     [[BAND_NAMES[k], F["n"][k], F["max_abs"][k], F["max_rel"][k], F["n_rel_gt_tol"][k]] for k in range(NB)]), ""]
        L += ["Check 6 (slope_se and the derived slope through the 16-bit codes) is in `roundtrip.md`.", ""]
        ds = X["distributions"]
        L += ["Per-phenotype distributions:", "",
              _table(["value", "min", "median", "p99", "max"],
                     [[k, v["min"], v["median"], v["p99"], v["max"]] for k, v in ds.items() if k != "anchor_minus_tss"]), "",
              f"Window start minus TSS: {ds['anchor_minus_tss']}", ""]
        vv = X["variant_values"]
        L += [f"Variant values: {vv['tested']:,} tested variants, af {_f(vv['af_min'])} to {_f(vv['af_max'])} ({vv['af_over_0_5']:,} above 0.5), "
              f"max ma_samples {vv['ma_samples_max']}, max ma_count {vv['ma_count_max']}. Cis variant list: {X['variants']['n_cis']:,} "
              f"(duplicate keys {X['variants']['duplicate_keys']}, chromosomes whose file order is not the byte-wise order: {X['variants']['chroms_file_order_not_vidx'] or 'none'}).", ""]
        L += ["20 worst full-precision rows (pass domain, by relative error):", "",
              _table(["phenotype", "vidx", "p", "slope", "slope_se", "se_hat", "|t|", "rel err"],
                     [[w["phenotype_id"], w["vidx"], _num(w["pval_nominal"]), _num(w["slope"]), _num(w["slope_se"]), _num(w["se_hat"]), _num(w["abs_t"]), _num(w["rel_err"])]
                      for w in X["stream"]["worst_full"]]), ""]
    if R["cross"]:
        L += ["## eQTL versus sQTL", "", _table(["measure", "value"], [[k, v if not isinstance(v, list) else str(v)] for k, v in R["cross"].items()]), ""]
    return "\n".join(L) + "\n"


# ---- encoding measurements: page size and codec (SPEC section 4) --------------------------------

PAGE_SIZES = [256, 512, 1024, 2048]
CODECS = ["raw", "zstd"]


def _stats(a) -> dict:
    a = np.asarray(a, dtype=np.float64)
    if len(a) == 0:
        return {"n": 0}
    return {"n": int(len(a)), "mean": float(a.mean()), "median": float(np.median(a)), "p90": float(np.percentile(a, 90)),
            "p99": float(np.percentile(a, 99)), "max": float(a.max()), "sum": float(a.sum())}


def _slices(ids: pa.Array) -> tuple[np.ndarray, np.ndarray]:
    n = len(ids)
    if n == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    starts = np.concatenate([[0], pc.indices_nonzero(pc.not_equal(ids.slice(0, n - 1), ids.slice(1))).to_numpy() + 1]).astype(np.int64)
    return starts, np.append(starts[1:], n)


def cmd_measure(cfg: Config, args) -> None:
    out = checks_dir(cfg)
    pk = cfg["packs"]
    level = int(pk["zstd_level"])
    dof = int(pk["dof"]["eqtl"])
    chroms = args.chrom or ["chr7", "chr10", "chr22"]
    for c in chroms:
        if not (out / f"phenotypes_e_{c}.parquet").exists():
            die(f"packcheck measure: run `packcheck run --type e --chrom {c}` first")
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"], temp_dir=out / "tmp" / "measure")
    con.execute("SET preserve_insertion_order = true")
    M: dict = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "chroms": chroms, "zstd_level": level, "dof_eqtl": dof}
    fetch = {f"{P}/{c}": [] for P in PAGE_SIZES for c in CODECS}
    file_bytes = {f"{P}/{c}": {} for P in PAGE_SIZES for c in CODECS}
    variants_info = {}
    pairs = {"rows": 0, "raw": 0, "zstd_interleaved": 0, "zstd_back_to_back": 0}
    det_json, det_z, blocks, nonpair, genes_rows = [], [], [], [], []
    t0 = time.time()
    for chrom in chroms:
        vpos = variants_sql(cfg, chrom)
        con.execute(f"""CREATE OR REPLACE TABLE v AS
            SELECT (row_number() OVER (ORDER BY position, A1, A2) - 1)::UINTEGER AS vidx, position, A1, A2, rs_number, match
            FROM {vpos} WHERE in_cis""")
        ve, vs = out / f"variant_values_e_{chrom}.parquet", out / f"variant_values_s_{chrom}.parquet"
        s_join = f"LEFT JOIN '{vs}' s USING (vidx)" if vs.exists() else ""
        s_col = lambda c: f"coalesce(e.{c}, s.{c})" if vs.exists() else f"e.{c}"          # noqa: E731
        vt = con.execute(f"""SELECT v.vidx, v.position, v.A1, v.A2, v.rs_number, v.match,
                                    {s_col('af')} AS af, {s_col('ma_samples')} AS ma_samples, {s_col('ma_count')} AS ma_count
                             FROM v LEFT JOIN '{ve}' e USING (vidx) {s_join} ORDER BY vidx""").fetch_arrow_table()
        ph = con.execute(f"""SELECT phenotype_id, var_start, n_rows, anchor, pos_first, pos_last
                             FROM '{out}/phenotypes_e_{chrom}.parquet' WHERE n_rows > 0 ORDER BY phenotype_id""").fetchall()
        A1, A2 = vt["A1"].to_pylist(), vt["A2"].to_pylist()
        snp = np.array([(a, b) in packfmt.SNP_CODES for a, b in zip(A1, A2)])
        nonsnp_len = [max(len(a), len(b)) for a, b, x in zip(A1, A2, snp) if not x]
        variants_info[chrom] = {"variants": vt.num_rows, "non_snp": int((~snp).sum()), "non_snp_share": float((~snp).mean()),
                                "heap_bytes": int(sum(len(a) + len(b) + 2 for a, b, x in zip(A1, A2, snp) if not x)),
                                "non_snp_allele_len": _stats(nonsnp_len), "genes": len(ph),
                                "no_value_variants": int(pc.sum(pc.is_null(vt["af"])).as_py() or 0)}
        cols = dict(position=vt["position"], rs_number=vt["rs_number"], af=vt["af"], ma_samples=vt["ma_samples"],
                    ma_count=vt["ma_count"], A1=A1, A2=A2, match=vt["match"].to_pylist())
        for P in PAGE_SIZES:
            for codec in CODECS:
                data, offs = packfmt.encode_variants_file(chrom, page_size=P, codec=codec, level=level, **cols)
                file_bytes[f"{P}/{codec}"][chrom] = len(data)
                fetch[f"{P}/{codec}"] += [packfmt.variant_range(offs, P, r[1], r[2])[1] for r in ph]
        log(f"packcheck measure: {chrom}: variant pages done ({time.time() - t0:.0f} s)")

        # pairs and blocks
        kind, src, _, _ = resolve_source(cfg, "e", chrom, "auto")
        rows = con.execute(f"""SELECT n.phenotype_id, v.vidx, n.pval_nominal, n.slope, n.slope_se
                               FROM ({source_sql('e', kind, src)}) n JOIN v ON v.position = n.position AND v.A1 = n.A1 AND v.A2 = n.A2
                               ORDER BY n.phenotype_id, v.vidx""").fetch_arrow_table()
        ids = rows["phenotype_id"].combine_chunks()
        starts, ends = _slices(ids)
        n = rows.num_rows
        vidx = rows["vidx"].to_numpy(zero_copy_only=False).astype(np.int64)
        p = rows["pval_nominal"].to_numpy(zero_copy_only=False).astype(np.float64)
        slope = rows["slope"].to_numpy(zero_copy_only=False).astype(np.float64)
        se = rows["slope_se"].to_numpy(zero_copy_only=False).astype(np.float64)
        span_of = {}
        for a, b in zip(starts.tolist(), ends.tolist()):
            span_of[ids[a].as_py()] = (a, b)
            q1, _ = packfmt.quantize_nlp(p[a:b])
            q2, _, _ = packfmt.quantize_se(se[a:b], slope[a:b])
            pr = np.empty(b - a, dtype=packfmt.PAIR_DTYPE)
            pr["nlp"], pr["se"] = q1, q2
            pb = pr.tobytes()
            pairs["rows"] += b - a
            pairs["raw"] += len(pb)
            pairs["zstd_interleaved"] += len(packfmt.zstd_frame(pb, level))
            pairs["zstd_back_to_back"] += len(packfmt.zstd_frame(q1.astype("<u2").tobytes() + q2.astype("<u2").tobytes(), level))
        log(f"packcheck measure: {chrom}: pairs done ({time.time() - t0:.0f} s)")

        susie = cfg.raw_dir("cis_eQTL_SuSiE") / f"topchef_{chrom}_MaxPC70.SuSiE_summary.parquet"
        cs = con.execute(f"""SELECT s.phenotype_id, v.vidx, max(s.pip) AS pip, arg_max(s.cs_id, s.pip)::TINYINT AS cs_id
                             FROM read_parquet('{susie}') s JOIN v ON v.position = s.position AND v.A1 = s.A1 AND v.A2 = s.A2
                             GROUP BY 1, 2 ORDER BY 1, 2""").fetch_arrow_table()
        cids = cs["phenotype_id"].combine_chunks()
        cst, cen = _slices(cids)
        cs_of = {cids[a].as_py(): (a, b) for a, b in zip(cst.tolist(), cen.tolist())}
        cs_v = cs["vidx"].to_numpy(zero_copy_only=False).astype(np.int64)
        cs_pip = cs["pip"].to_numpy(zero_copy_only=False)
        cs_id = cs["cs_id"].to_numpy(zero_copy_only=False)
        info = {r[0]: r[1:] for r in ph}
        _, chrom_details = steps_pack._details(cfg, con, chrom)
        for g, det in chrom_details.items():
            dj = packfmt.details_bytes(det)
            det_json.append(len(dj))
            det_z.append(len(packfmt.zstd_frame(dj, level)))
            if g in info:
                var_start, n_rows, anchor, pos_first, pos_last = info[g]
                a, b = span_of[g]
                if b - a != n_rows or not np.array_equal(vidx[a:b], np.arange(var_start, var_start + n_rows)):
                    raise RuntimeError(f"measure: {g} rows are not the contiguous run the check found")
                if g in cs_of:
                    ca, cb = cs_of[g]
                    crow = cs_v[ca:cb] - var_start
                    cargs = (crow, cs_pip[ca:cb], cs_id[ca:cb])
                else:
                    cargs = (None, None, None)
                blk = packfmt.encode_gene_block(det, var_start, anchor, p[a:b], slope[a:b], se[a:b], *cargs, level,
                                                pos_first=pos_first, pos_last=pos_last)
            else:
                n_rows = 0
                blk = packfmt.encode_gene_block(det, None, None, [], [], [], None, None, None, level)
            blocks.append(len(blk))
            nonpair.append(len(blk) - 4 * n_rows)
            genes_rows.append(n_rows)
        log(f"packcheck measure: {chrom}: blocks done ({time.time() - t0:.0f} s)")

    # genome-wide details and exact eQTL pack size (block length does not depend on the p/slope values)
    nrows = dict(con.execute(f"SELECT phenotype_id, n_rows FROM read_parquet('{out}/phenotypes_e_chr*.parquet')").fetchall())
    ncs = dict(con.execute(f"""SELECT gene_id, count(DISTINCT (position::VARCHAR || ':' || A1 || ':' || A2))
                               FROM '{cfg.tables / 'credible_sets.parquet'}' WHERE qtl_type = 'e' GROUP BY 1""").fetchall())
    gw_json, gw_z, pack_bytes, n_genes = [], [], 0, 0
    per_chr_pack: dict[str, int] = {}
    for c in CHROMS:
        _, chrom_details = steps_pack._details(cfg, con, c)
        for g, det in chrom_details.items():
            dj = packfmt.details_bytes(det)
            dz = packfmt.zstd_frame(dj, level)
            gw_json.append(len(dj))
            gw_z.append(len(dz))
            body = packfmt.BLOCK_HEADER_LEN + 4 * nrows.get(g, 0) + packfmt.CS_RECORD_LEN * (ncs.get(g, 0) if nrows.get(g, 0) else 0) + len(dz)
            per_chr_pack[c] = per_chr_pack.get(c, packfmt.FILE_HEADER_LEN) + body + packfmt.pad4(body)
            n_genes += 1
    pack_bytes = sum(per_chr_pack.values())
    con.close()
    shutil.rmtree(out / "tmp", ignore_errors=True)

    # decisions
    mean_fetch = {k: float(np.mean(v)) for k, v in fetch.items()}
    best = {c: min(PAGE_SIZES, key=lambda P: mean_fetch[f"{P}/{c}"]) for c in CODECS}
    best_mean = {c: mean_fetch[f"{best[c]}/{c}"] for c in CODECS}
    saving = 1 - best_mean["zstd"] / best_mean["raw"]
    codec = "zstd" if saving >= 0.30 else "raw"
    lo = min(mean_fetch[f"{P}/{codec}"] for P in PAGE_SIZES)
    page_size = max(P for P in PAGE_SIZES if mean_fetch[f"{P}/{codec}"] <= lo * 1.05)
    M["variant_pages"] = {
        "info": variants_info,
        "by_setting": {k: {"file_bytes": file_bytes[k], "total_bytes": sum(file_bytes[k].values()),
                           "bytes_per_variant": sum(file_bytes[k].values()) / sum(x["variants"] for x in variants_info.values()),
                           "gene_fetch": _stats(v)} for k, v in fetch.items()},
    }
    M["pairs"] = {**pairs, "bytes_per_row": {k: pairs[k] / pairs["rows"] for k in ("raw", "zstd_interleaved", "zstd_back_to_back")}}
    M["details"] = {"chroms": {"json": _stats(det_json), "zstd": _stats(det_z)}, "genome": {"json": _stats(gw_json), "zstd": _stats(gw_z), "genes": n_genes}}
    M["blocks"] = {"bytes": _stats(blocks), "non_pair_bytes": _stats(nonpair), "genes": len(blocks), "rows": int(sum(genes_rows)),
                   "genome_eqtl_pack_bytes": pack_bytes, "genome_eqtl_pack_bytes_by_chrom": per_chr_pack,
                   "genome_eqtl_rows": int(sum(nrows.values()))}
    M["decisions"] = {"best_page_size": best, "best_mean_fetch": best_mean, "zstd_saving_at_best": saving,
                      "variant_page_codec": codec, "variant_page_size": page_size,
                      "config_before": {"variant_page_codec": pk["variant_page_codec"], "variant_page_size": pk["variant_page_size"]}}
    changed = codec != pk["variant_page_codec"] or page_size != pk["variant_page_size"]
    if changed:
        path = ROOT / "pipeline" / "config.yaml"
        text = path.read_text()
        def setv(key, val, text):          # keep the comment column where it was
            return re.subn(rf"(\n  {key}: )(\S+)( +)#", lambda m: m.group(1) + str(val).ljust(len(m.group(2)) + len(m.group(3)) - 1) + " #", text)
        text, n1 = setv("variant_page_size", page_size, text)
        text, n2 = setv("variant_page_codec", codec, text)
        if n1 != 1 or n2 != 1:
            die("packcheck measure: could not update packs.variant_page_size / variant_page_codec in config.yaml")
        path.write_text(text)
        log(f"packcheck measure: config packs changed to variant_page_size {page_size}, variant_page_codec {codec}")
    M["decisions"]["config_updated"] = changed

    # inverse-t reference vectors for the TypeScript reader
    rep = json.loads((out / "report.json").read_text()) if (out / "report.json").exists() else {}
    mins = [_num(t["phenotypes"]["min_p"]) for t in rep.get("types", {}).values()]
    min_p = min(mins) if mins else 1e-300
    special = [min_p, 1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95,
               0.99, 0.999, 1 - 1e-4, 1 - 1e-6, 1 - 1e-9, 1 - 1e-12, 1 - 1e-15, np.nextafter(1.0, 0.0), 1.0]
    vec = []
    spec = np.unique([x for x in special if min_p <= x <= 1])
    k = 500 - len(spec)
    while True:                         # 500 p values per dof: the specials plus log-spaced fill
        logs = np.logspace(math.log10(min_p), 0, k)
        logs = logs[np.min(np.abs(np.log10(logs)[:, None] - np.log10(spec)[None, :]), axis=1) > 1e-6]   # drop near-duplicates of specials
        if len(logs) + len(spec) >= 500:
            break
        k += 1
    logs = logs[np.round(np.linspace(0, len(logs) - 1, 500 - len(spec))).astype(int)]
    grid = np.unique(np.concatenate([logs, spec]))
    for d in sorted(set(int(x) for x in pk["dof"].values())):
        for pv in grid:
            vec.append({"dof": d, "p": float(pv), "t": float(-stdtrit(d, pv / 2)) + 0.0})   # + 0.0 turns -0.0 into 0.0
    import scipy
    write_json(out / "se_reference.json", {"description": "t = -scipy.special.stdtrit(dof, p / 2), clamped at >= 0 by readers; slope = sign * slope_se * t",
                                           "scipy": scipy.__version__, "min_p": min_p, "n": len(vec), "vectors": vec})
    M["se_reference"] = {"n": len(vec), "min_p": min_p}
    M["seconds"] = round(time.time() - t0, 1)
    write_json(out / "measure.json", M)
    (out / "measure.md").write_text(_measure_md(M))
    log(f"packcheck measure: codec {codec} (zstd saves {saving:.1%} at best sizes {best}), page size {page_size}; "
        f"wrote {out / 'measure.md'}")


def _measure_md(M: dict) -> str:
    vp = M["variant_pages"]
    L = [f"# packcheck measure", "", f"Generated {M['generated']} on {', '.join(M['chroms'])}; zstd level {M['zstd_level']}.", "",
         "## Variant pages", "",
         _table(["chr", "cis variants", "non-SNP", "non-SNP share", "heap bytes", "non-SNP longest allele median / p90 / p99 / max", "eQTL genes"],
                [[c, x["variants"], x["non_snp"], f"{100 * x['non_snp_share']:.2f}%", x["heap_bytes"],
                  f"{_f(x['non_snp_allele_len'].get('median'))} / {_f(x['non_snp_allele_len'].get('p90'))} / {_f(x['non_snp_allele_len'].get('p99'))} / {_f(x['non_snp_allele_len'].get('max'))}",
                  x["genes"]] for c, x in vp["info"].items()]), "",
         _table(["page size / codec", "total bytes", "bytes per variant"] + [f"{c} bytes" for c in M["chroms"]] + ["gene fetch mean", "p90", "max"],
                [[k, v["total_bytes"], v["bytes_per_variant"]] + [v["file_bytes"][c] for c in M["chroms"]]
                 + [v["gene_fetch"]["mean"], v["gene_fetch"]["p90"], v["gene_fetch"]["max"]] for k, v in vp["by_setting"].items()]), ""]
    d = M["decisions"]
    L += [f"Decision: best mean gene fetch raw {_f(d['best_mean_fetch']['raw'])} B at P={d['best_page_size']['raw']}, "
          f"zstd {_f(d['best_mean_fetch']['zstd'])} B at P={d['best_page_size']['zstd']}; zstd saves {100 * d['zstd_saving_at_best']:.1f}% "
          f"(threshold 30%) -> codec **{d['variant_page_codec']}**, page size **{d['variant_page_size']}** "
          f"(config {'updated' if d['config_updated'] else 'unchanged'}).", ""]
    pr = M["pairs"]
    L += ["## Pairs", "", _table(["layout", "bytes per row", "saving vs raw"],
                                 [[k, v, f"{100 * (1 - v / 4.0):.1f}%"] for k, v in pr["bytes_per_row"].items()]),
          "", f"{pr['rows']:,} eQTL rows; a pair is u16 -log10 p code + u16 SE code (sign bit and log(SE)).", ""]
    de = M["details"]
    L += ["## Gene details", "", _table(["scope", "genes", "JSON median", "JSON p99", "JSON max", "zstd median", "zstd p99", "zstd max"],
                                        [[s, x["json"]["n"], x["json"]["median"], x["json"]["p99"], x["json"]["max"], x["zstd"]["median"], x["zstd"]["p99"], x["zstd"]["max"]]
                                         for s, x in (("chosen chromosomes", de["chroms"]), ("genome", de["genome"]))]), ""]
    b = M["blocks"]
    L += ["## Blocks", "", f"{b['genes']:,} genes on the chosen chromosomes ({b['rows']:,} rows): block bytes mean {_f(b['bytes']['mean'])}, "
          f"p90 {_f(b['bytes']['p90'])}, max {_f(b['bytes']['max'])}; bytes besides pairs mean {_f(b['non_pair_bytes']['mean'])}, "
          f"median {_f(b['non_pair_bytes']['median'])}. Genome-wide eQTL pack (exact for the v0 layout, {b['genome_eqtl_rows']:,} rows): "
          f"{b['genome_eqtl_pack_bytes']:,} bytes ({b['genome_eqtl_pack_bytes'] / 1e6:.1f} MB).", ""]
    L += [f"Inverse-t reference vectors: {M['se_reference']['n']} in se_reference.json (p from {M['se_reference']['min_p']:.3g} to 1).", ""]
    return "\n".join(L) + "\n"


def main(cfg: Config, args) -> None:
    globals()[f"cmd_{args.action}"](cfg, args)
