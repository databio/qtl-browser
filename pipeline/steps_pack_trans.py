"""Steps pack_variants_trans and pack_trans (SPEC.md sections 4 and 12), and their checks.

`pack_variants_trans` runs after `pack_eqtl`. One DuckDB pass over the trans rows gives each trans
variant its af (`_tmp/trans_variant_af.parquet`, one value per variant, checked). Then, per
chromosome, it rewrites the variants file with `packfmt.encode_variants_file`: the cis records as
decoded from the file `pack_eqtl` wrote, plus the chromosome's trans-only variants (`NOT in_cis` rows
of `_tables/variants.parquet`) as the trans-only section. Before it replaces the file it checks that the
new bytes from byte 32 to the end of the last cis page equal the old file's. That is what keeps every
`var_off`/`var_len` in `search_index` and every run in the eQTL and sQTL packs valid, so those packs
are not rebuilt.

`pack_trans` publishes `trans/<gene chr>` under `data/derived/immutable/` for chr1..chr22, chrX,
and chrM: one zstd frame per gene with trans rows, in `search_index` order (chr, tss, gene_id). Each
gene's `trans_off`/`trans_len` goes to `_tmp/pack_pointers/trans_<chr>.parquet` (stats in `.json`),
which the `search_index` step joins.

`validate` reads both with its own frame reader, written from SPEC rather than from `packfmt`, and
writes reference files for `npm run pack-check` to `data/derived/_tmp/pack_check/trans/`.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import struct
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.special import stdtrit

from . import packfmt
from .common import (CHROMS, Config, connect, log, publish_file, read_search_index, register_search_index, stage,
                     variants_sql, write_parquet)
from .steps_pack import (EXT, PackError, _file_header, _read_pages, _require, _unframe, _write_arrow, pack_file,
                         pack_paths, pointer_dir)

TRANS_CHROMS = CHROMS + ["chrM"]      # gene chromosomes that can have a trans file (trans_pairs/chr=chrY/ is empty)
EXPECT = {"variants": 9_215_026, "cis": 8_872_723, "trans_only": 342_303, "no_alleles": 32_765, "af_variants": 6_594_666,
          "rows": 15_862_525, "eqtl_rows": 2_680_117, "sqtl_rows": 13_182_408}
MIN_BETA_CODE = 16                    # beta_se = |beta| / t never comes from a tiny code (measured minimum 575)
TRANS_POINTER_SCHEMA = pa.schema([("gene_id", pa.string()), ("trans_off", pa.uint32()), ("trans_len", pa.uint32())])
REF_GENES = {"ENSG00000128591": "FLNC", "ENSG00000010282": "HHATL"}
REF_RANDOM = 20
REF_VARIANT_CHROMS = ("chr1", "chrX")
# SPEC section 9 limits for trans rows, also written to the reference files
TOL = {"nlp_half_step_divisor": 131066, "nlp_slack": 1e-9, "beta_half_step_divisor": 65534, "beta_rel_f32": 1.2e-7,
       "beta_se_rel": 0.01, "r2_abs": 0.002, "af_tol": 0.5 / 65534 + 1e-12}
MATCH = {"none": 0, "exact": 1, "position": 2}

CHR_CODE_SQL = "CASE WHEN {c} = 'chrX' THEN 23 ELSE substr({c}, 4)::INTEGER END"
INTRON_SQL = """CASE WHEN {t}.qtl_type = 's' THEN split_part({t}.phenotype_id, ':', 2)::UINTEGER END AS intron_start,
       CASE WHEN {t}.qtl_type = 's' THEN split_part({t}.phenotype_id, ':', 3)::UINTEGER END AS intron_end,
       CASE WHEN {t}.qtl_type = 's' THEN regexp_extract(split_part({t}.phenotype_id, ':', 4), '^clu_([0-9]+)_[+-]$', 1)::UINTEGER END AS cluster,
       CASE WHEN {t}.qtl_type = 's' THEN right(split_part({t}.phenotype_id, ':', 4), 1) END AS strand"""


def trans_source(cfg: Config, chrom: str) -> Path:
    """One gene chromosome's trans rows. The variant-page stage moves these files under data/derived/_tables/."""
    return cfg.tables / "trans" / f"chr={chrom}" / "data.parquet"


def trans_path(cfg: Config, chrom: str) -> Path:
    return pack_file(cfg, f"trans/{chrom}", EXT["trans"])


def variant_af_path(cfg: Config) -> Path:
    return cfg.tmp / "trans_variant_af.parquet"


def _trans_sources(cfg: Config) -> list[Path]:
    """Every trans_pairs file, which must all be gene chromosomes in TRANS_CHROMS."""
    found = set((cfg.tables / "trans").glob("chr=*/data.parquet"))
    known = [trans_source(cfg, c) for c in TRANS_CHROMS if trans_source(cfg, c).exists()]
    _require(bool(known), "trans_pairs has no files; run `build --step trans`")
    _require(found == set(known), f"trans_pairs files outside {TRANS_CHROMS}: {sorted(str(p) for p in found - set(known))}")
    return known


def _scan(files: list[Path]) -> str:
    return "read_parquet([" + ", ".join(f"'{f}'" for f in files) + "], hive_partitioning = false)"


def _vpos(cfg: Config, chrom: str | None = None) -> str:
    return variants_sql(cfg, chrom)


# ---- build: the trans-only section of each variants file --------------------------------------------
def variants_trans(cfg: Config) -> None:
    d, pk = cfg.derived, cfg["packs"]
    level, page_size, codec = int(pk["zstd_level"]), int(pk["variant_page_size"]), pk["variant_page_codec"]
    con = connect(cfg, memory_limit=pk["duckdb_memory_limit"], threads=int(pk["duckdb_threads"]), temp_dir=cfg.tmp / "pack_variants_trans")
    afp = variant_af_path(cfg)
    tmp_af = afp.with_name("trans_variant_af.tmp.parquet")
    t0 = time.time()
    con.execute(f"""COPY (SELECT variant_chr, position, min(af) AS af, count(DISTINCT af) AS n_af, count(*) AS n_rows
        FROM {_scan(_trans_sources(cfg))} GROUP BY 1, 2) TO '{tmp_af}' (FORMAT parquet, COMPRESSION zstd)""")
    n_af, multi, rows = con.execute(f"SELECT count(*), count(*) FILTER (WHERE n_af <> 1), sum(n_rows) FROM '{tmp_af}'").fetchone()
    _require(multi == 0, f"trans_variant_af: {multi:,} trans variants do not have exactly one af over their trans rows")
    os.replace(tmp_af, afp)
    log(f"pack_variants_trans: {n_af:,} trans variants (expected {EXPECT['af_variants']:,}) from {int(rows):,} trans rows, "
        f"one af each; {afp.relative_to(d)} in {time.time() - t0:.1f} s")

    pdir = pointer_dir(cfg)
    pdir.mkdir(parents=True, exist_ok=True)
    tot = dict.fromkeys(("cis", "trans_only", "no_alleles", "cis_bytes", "trans_bytes", "bytes", "old_bytes"), 0)
    for chrom in CHROMS:
        t0 = time.time()
        key = f"variants/{chrom}"
        path = pack_paths(cfg, chrom)[0]
        _require(path.exists(), f"{chrom}: {path} is missing; run `build --step pack_eqtl` first")
        old = path.read_bytes()
        h = packfmt.parse_file_header(old)
        n_cis, want_tr = con.execute(f"SELECT count(*) FILTER (WHERE in_cis), count(*) FILTER (WHERE NOT in_cis) FROM {_vpos(cfg, chrom)}").fetchone()
        # the cis pages already on disk: walk page headers from byte 32 until a page would start at vidx n_cis. n_cis comes
        # from the variant table, not the old header, and the guard below compares these bytes, not the header.
        old_cis_end = packfmt.FILE_HEADER_LEN
        while old_cis_end < len(old) and struct.unpack_from("<I", old, old_cis_end + 4)[0] < n_cis:
            old_cis_end += packfmt.PAGE_HEADER_LEN + struct.unpack_from("<I", old, old_cis_end)[0]
            old_cis_end += packfmt.pad4(old_cis_end)
        cur = packfmt.decode_variant_pages(old[packfmt.FILE_HEADER_LEN:old_cis_end])
        _require(h["kind"] == packfmt.KIND_VARIANTS and h["page_size"] == page_size and cur["codec"] == codec
                 and int(cur["vidx"][0]) == 0 and len(cur["vidx"]) == n_cis,
                 f"{chrom}: {path.name} pages before vidx {n_cis:,} hold vidx {int(cur['vidx'][0]):,}..{int(cur['vidx'][-1]):,} in pages of "
                 f"{h['page_size']} ({cur['codec']}); the variant table has {n_cis:,} cis variants, config pages of {page_size} ({codec})")
        tr = con.execute(f"""SELECT v.position, v.A1, v.A2, v.rs_number, v.match, a.af, a.position IS NOT NULL AS joined
            FROM {_vpos(cfg, chrom)} v LEFT JOIN '{afp}' a ON a.variant_chr = ? AND a.position = v.position
            WHERE NOT v.in_cis ORDER BY v.position""", [chrom]).fetch_arrow_table()
        m = tr.num_rows
        a1, a2 = tr["A1"].to_pylist(), tr["A2"].to_pylist()
        unjoined = m - int(np.sum(tr["joined"].to_numpy(zero_copy_only=False)))
        half = sum((x is None) != (y is None) for x, y in zip(a1, a2))
        no_al = sum(x is None for x in a1)
        _require(m == want_tr and unjoined == 0 and half == 0,
                 f"{chrom}: {m:,} trans-only variants (table {want_tr:,}): {unjoined:,} have no trans row to give an af, {half:,} have one null allele")
        c = slice(0, n_cis)
        new, offs = packfmt.encode_variants_file(
            chrom, cur["position"][c], cur["rs_number"][c], cur["af"][c], cur["ma_samples"][c], cur["ma_count"][c],
            cur["A1"][c], cur["A2"][c], cur["match"][c], page_size, codec, level,
            trans_only={"position": tr["position"].to_numpy(), "rs_number": tr["rs_number"].to_numpy(zero_copy_only=False),
                        "af": tr["af"].to_numpy(zero_copy_only=False), "A1": a1, "A2": a2, "match": tr["match"].to_pylist()})
        # the guard: every cis byte stays where search_index and the eQTL/sQTL packs expect it
        cis_end = int(offs[-(-n_cis // page_size)])
        _require(cis_end == old_cis_end and new[packfmt.FILE_HEADER_LEN:cis_end] == old[packfmt.FILE_HEADER_LEN:old_cis_end],
                 f"{chrom}: the cis pages would change (cis end {old_cis_end:,} -> {cis_end:,} B), so var_off/var_len and the "
                 f"eQTL and sQTL packs would no longer match; not replacing {path.name}")
        back = packfmt.decode_variants_file(new)
        _require(back["n_cis"] == n_cis and len(back["position"]) == n_cis + m and int(back["no_alleles"].sum()) == no_al,
                 f"{chrom}: the new file does not decode to {n_cis:,} + {m:,} variants with {no_al:,} lacking alleles")
        tmp = stage(cfg, key, EXT["variants"])
        tmp.write_bytes(new)
        # a new name, because the bytes changed: search_index's var_off/var_len are byte offsets into
        # this file, so it is rebuilt after this step (SPEC section 3)
        path = publish_file(cfg, tmp, key, EXT["variants"])
        st = {"chrom": chrom, "cis": n_cis, "trans_only": m, "no_alleles": no_al, "cis_bytes": cis_end,
              "trans_bytes": len(new) - cis_end, "bytes": len(new), "old_bytes": len(old), "seconds": round(time.time() - t0, 1)}
        (pdir / f"variants_trans_{chrom}.json").write_text(json.dumps(st, indent=2))
        for k in tot:
            tot[k] += st[k]
        log(f"pack_variants_trans {chrom}: n_cis {n_cis:,}, n_trans_only {m:,} ({no_al:,} without alleles), trans-only section "
            f"{st['trans_bytes']:,} B; cis bytes unchanged; file {len(old):,} -> {len(new):,} B in {st['seconds']} s")
        del cur, back, tr, old, new
    for k in ("cis", "trans_only", "no_alleles"):
        log(f"pack_variants_trans total {k}: {tot[k]:,} (expected {EXPECT[k]:,}{'' if tot[k] == EXPECT[k] else ', DIFFERS'})")
    log(f"pack_variants_trans total bytes: variants files {tot['old_bytes']:,} -> {tot['bytes']:,}; trans-only sections {tot['trans_bytes']:,} B")


# ---- build: the trans pack ------------------------------------------------------------------------------
def pack(cfg: Config) -> None:
    d, pk = cfg.derived, cfg["packs"]
    level = int(pk["zstd_level"])
    con = connect(cfg, memory_limit=pk["duckdb_memory_limit"], threads=int(pk["duckdb_threads"]), temp_dir=cfg.tmp / "pack_trans")
    con.execute(f"CREATE VIEW variants AS SELECT chr, position, rs_number FROM {_vpos(cfg)}")
    # search_index's rows are genes.parquet's in (chr, tss, gene_id) order; reading genes.parquet lets a fresh
    # build run this step before the search_index step. validate checks the frames against search_index itself.
    con.execute(f"CREATE TABLE si AS SELECT gene_id, chr, tss FROM '{cfg.tables / 'genes.parquet'}'")
    pdir = pointer_dir(cfg)
    pdir.mkdir(parents=True, exist_ok=True)
    sources = _trans_sources(cfg)
    for chrom in TRANS_CHROMS:
        src, key = trans_source(cfg, chrom), f"trans/{chrom}"
        if not src.exists():
            log(f"pack_trans {chrom}: no trans rows, no file")
            continue
        t0 = time.time()
        n_src = pq.read_metadata(src).num_rows
        tb = con.execute(f"""SELECT t.gene_id, si.chr AS gene_chr, t.qtl_type, t.variant_chr, t.position, t.af, t.pval, t.beta, v.rs_number,
                   {INTRON_SQL.format(t='t')}
            FROM read_parquet('{src}', hive_partitioning = false) t
            JOIN si USING (gene_id)
            JOIN variants v ON v.chr = t.variant_chr AND v.position = t.position
            ORDER BY si.tss, t.gene_id, t.qtl_type, intron_start, intron_end, cluster, strand,
                     {CHR_CODE_SQL.format(c='t.variant_chr')}, t.position""").fetch_arrow_table()
        _require(tb.num_rows == n_src, f"{chrom}: {tb.num_rows:,} trans rows join to a gene and to the variant table; {src.name} has {n_src:,}")
        gene_chrs = set(tb["gene_chr"].unique().to_pylist())
        _require(gene_chrs == {chrom}, f"{chrom}: the trans genes sit on {sorted(gene_chrs)}")
        gid = tb["gene_id"].to_numpy(zero_copy_only=False)
        starts = np.r_[0, np.flatnonzero(gid[1:] != gid[:-1]) + 1].astype(np.int64)
        ends = np.r_[starts[1:], len(gid)].astype(np.int64)
        rank = {g: i for i, (g,) in enumerate(con.execute("SELECT gene_id FROM si WHERE chr = ? ORDER BY tss, gene_id", [chrom]).fetchall())}
        _require(bool(np.all(np.diff([rank[g] for g in gid[starts]]) > 0)),
                 f"{chrom}: frames would not follow search_index order, or a gene's rows are not contiguous")
        qt, vc, strand = tb["qtl_type"].to_pylist(), tb["variant_chr"].to_pylist(), tb["strand"].to_pylist()
        cols = {k: tb[k].to_numpy(zero_copy_only=False) for k in ("position", "af", "pval", "beta", "rs_number", "intron_start", "intron_end", "cluster")}
        del tb
        is_s = np.array([q == "s" for q in qt], dtype=bool)
        ptr: dict[str, list] = {"gene_id": [], "trans_off": [], "trans_len": []}
        lens, nrows = [], []
        min_code, min_code_gene, max_k = packfmt.BETA_MAXQ + 1, None, 0
        tmp = stage(cfg, key, EXT["trans"])
        with open(tmp, "wb") as fh:
            fh.write(packfmt.file_header(packfmt.KIND_TRANS, chrom, len(starts), 0))
            for s, e in zip(starts.tolist(), ends.tolist()):
                frame = packfmt.encode_trans_frame(qt[s:e], vc[s:e], cols["position"][s:e], cols["rs_number"][s:e], cols["af"][s:e],
                                                   cols["pval"][s:e], cols["beta"][s:e], cols["intron_start"][s:e], cols["intron_end"][s:e],
                                                   cols["cluster"][s:e], strand[s:e], level)
                ptr["gene_id"].append(str(gid[s]))
                ptr["trans_off"].append(fh.tell())
                ptr["trans_len"].append(len(frame))
                fh.write(frame)
                lens.append(len(frame))
                nrows.append(e - s)
                ab = np.abs(cols["beta"][s:e].astype(np.float64))
                code = int(np.rint(ab.min() / ab.max() * packfmt.BETA_MAXQ))
                if code < min_code:
                    min_code, min_code_gene = code, str(gid[s])
                if is_s[s:e].any():
                    m = np.flatnonzero(is_s[s:e]) + s
                    max_k = max(max_k, len(set(zip(cols["intron_start"][m].tolist(), cols["intron_end"][m].tolist(),
                                                   cols["cluster"][m].tolist(), [strand[i] for i in m]))))
            size = fh.tell()
        _require(size == packfmt.FILE_HEADER_LEN + sum(lens) and size <= packfmt.U32_MAX, f"{chrom}: trans pack is {size:,} bytes (frames {sum(lens):,}; limit 4 GiB)")
        _require(min_code >= MIN_BETA_CODE, f"{chrom}: smallest |beta code| is {min_code} ({min_code_gene}), below {MIN_BETA_CODE}")
        publish_file(cfg, tmp, key, EXT["trans"])
        write_parquet(pa.Table.from_pydict(ptr, schema=TRANS_POINTER_SCHEMA), pdir / f"trans_{chrom}.parquet", 100_000)
        L, R = np.array(lens), np.array(nrows)
        n_e = int((~is_s).sum())
        st = {"chrom": chrom, "genes": len(lens), "rows": int(R.sum()), "eqtl_rows": n_e, "sqtl_rows": int(R.sum()) - n_e, "bytes": size,
              "source_bytes": src.stat().st_size, "median_frame": int(np.median(L)), "max_frame": int(L.max()),
              "max_frame_gene": ptr["gene_id"][int(L.argmax())], "max_rows": int(R.max()), "max_introns": max_k,
              "min_beta_code": min_code, "min_beta_code_gene": min_code_gene, "seconds": round(time.time() - t0, 1)}
        (pdir / f"trans_{chrom}.json").write_text(json.dumps(st, indent=2))
        log(f"pack_trans {chrom}: {st['genes']:,} genes, {st['rows']:,} rows ({n_e:,} eQTL), {size:,} B ({size / st['rows']:.2f} B/row; "
            f"parquet {st['source_bytes']:,} B); frame median {st['median_frame']:,} B, max {st['max_frame']:,} B ({st['max_frame_gene']}, "
            f"{st['max_rows']:,} rows max); up to {max_k} introns; smallest |beta code| {min_code}; {st['seconds']} s")
        del cols, qt, vc, strand, gid
    stats = [json.loads((pdir / f"trans_{c}.json").read_text()) for c in TRANS_CHROMS if trans_source(cfg, c).exists()]
    tot = {k: sum(s[k] for s in stats) for k in ("genes", "rows", "eqtl_rows", "sqtl_rows", "bytes", "source_bytes")}
    for k in ("rows", "eqtl_rows", "sqtl_rows"):
        log(f"pack_trans total {k}: {tot[k]:,} (expected {EXPECT[k]:,}{'' if tot[k] == EXPECT[k] else ', DIFFERS'})")
    big = max(stats, key=lambda s: s["max_frame"])
    log(f"pack_trans total: {tot['genes']:,} genes in {len(stats)} files ({len(sources)} trans_pairs files), {tot['bytes']:,} B "
        f"({tot['bytes'] / tot['rows']:.2f} B/row; trans_pairs {tot['source_bytes']:,} B); largest frame {big['max_frame']:,} B "
        f"({big['max_frame_gene']}, {big['chrom']}); smallest |beta code| {min(s['min_beta_code'] for s in stats)}; "
        f"most introns in a gene {max(s['max_introns'] for s in stats)}")


# ---- validate: an independent reader of SPEC sections 4 and 12 ---------------------------------------
_TRANS_H = struct.Struct("<4sIIHHdd")


def _up4(x: int) -> int:
    return (x + 3) & ~3


def _read_trans_frame(frame: bytes, what: str) -> dict:
    """SPEC section 12 reader rules for one frame. Returns the header, the intron table (start, end, cluster,
    strand code per row of `table`), and the row codes with absolute positions."""
    raw = _unframe(frame, None, what)
    if len(raw) < 32:
        raise PackError(f"{what}: payload of {len(raw)} bytes is shorter than its header")
    magic, n_e, n_s, k, reserved, nlp_max, beta_max = _TRANS_H.unpack_from(raw, 0)
    n = n_e + n_s
    if magic != b"QTT0" or reserved or n == 0 or (k == 0) != (n_s == 0) or k > 255 \
            or not (math.isfinite(nlp_max) and math.isfinite(beta_max) and nlp_max >= 0 and beta_max >= 0):
        raise PackError(f"{what}: header magic {magic!r}, reserved {reserved}, n_e {n_e}, n_s {n_s}, k {k}, nlp_max {nlp_max}, beta_max {beta_max}")
    table_end = 32 + 13 * k
    rows = _up4(table_end)
    codes_end = rows + 14 * n
    chr_at = _up4(codes_end)
    body_end = chr_at + n + n_s
    if len(raw) != _up4(body_end) or any(raw[table_end:rows]) or any(raw[codes_end:chr_at]) or any(raw[body_end:]):
        raise PackError(f"{what}: payload of {len(raw)} bytes (layout needs {_up4(body_end)}) or nonzero padding")

    def u32(off: int, cnt: int) -> np.ndarray:
        return np.frombuffer(raw, "<u4", cnt, off).astype(np.int64)
    table = np.stack([u32(32, k), u32(32 + 4 * k, k), u32(32 + 8 * k, k), np.frombuffer(raw, "u1", k, 32 + 12 * k).astype(np.int64)], axis=1)
    if k and (np.any(table[:, 3] > 1) or any(tuple(table[i]) >= tuple(table[i + 1]) for i in range(k - 1))):
        raise PackError(f"{what}: intron table has a strand code above 1 or is not strictly ascending")
    pcol = u32(rows, n)
    rs = u32(rows + 4 * n, n)
    af = np.frombuffer(raw, "<u2", n, rows + 8 * n).astype(np.int64)
    nq = np.frombuffer(raw, "<u2", n, rows + 10 * n).astype(np.int64)
    bq = np.frombuffer(raw, "<i2", n, rows + 12 * n).astype(np.int64)
    vchr = np.frombuffer(raw, "u1", n, chr_at).astype(np.int64)
    intron = np.frombuffer(raw, "u1", n_s, chr_at + n).astype(np.int64)
    if np.any(nq > 65533) or np.any(af > 65534) or np.any(bq == -32768) or np.any((vchr < 1) | (vchr > 23)):
        raise PackError(f"{what}: an nlp code above 65533, an af code above 65534, a beta code of -32768, or a variant_chr outside 1..23")
    if n_s and (np.any(intron >= k) or np.any(np.diff(intron) < 0) or len(set(intron.tolist())) != k):
        raise PackError(f"{what}: intron indices must be below k = {k}, never decrease, and use every intron")
    run = np.r_[np.full(n_e, -1), intron]
    same_run = np.r_[False, run[1:] == run[:-1]]
    if np.any(same_run[1:] & (vchr[1:] < vchr[:-1])):
        raise PackError(f"{what}: variant_chr decreases inside a run")
    if np.any(pcol < 1):
        raise PackError(f"{what}: a position entry is 0")
    fresh = ~(same_run & np.r_[False, vchr[1:] == vchr[:-1]])
    csum = np.cumsum(pcol)
    first = np.flatnonzero(fresh)
    pos = csum - (csum[first] - pcol[first])[np.cumsum(fresh) - 1]
    if pos.max() > 0xFFFFFFFF:
        raise PackError(f"{what}: a position above 2^32-1")
    return {"n_e": n_e, "n_s": n_s, "k": k, "nlp_max": nlp_max, "beta_max": beta_max, "table": table, "position": pos,
            "rs_number": rs, "af_code": af, "nlp_code": nq, "beta_code": bq, "variant_chr": vchr, "intron": intron}


def validate(cfg: Config, con, check) -> None:
    """SPEC sections 4 and 12 on every chromosome: the trans-only sections, the af of trans variants, every
    trans frame against its source rows, and the reference files for `npm run pack-check`."""
    d, pk = cfg.derived, cfg["packs"]
    t_start = time.time()
    dof = {"e": int(pk["dof"]["eqtl"]), "s": int(pk["dof"]["sqtl"])}
    con.execute(f"SET memory_limit = '{pk['duckdb_memory_limit']}'")
    con.execute(f"SET threads = {int(pk['duckdb_threads'])}")
    idx = register_search_index(cfg, con, "si_trans")
    cols = read_search_index(cfg).column_names
    if not all(c in cols for c in ("trans_off", "trans_len", "gene_version")):
        check(False, "search_index has trans_off, trans_len, gene_version (run `build --step pack_trans --step search_index --force`)")
        return
    try:
        sources = _trans_sources(cfg)
    except RuntimeError as e:
        check(False, f"trans sources: {e}")
        return
    out_dir = cfg.tmp / "pack_check" / "trans"
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(parents=True)

    # 1. a position joins to exactly one variant
    n, distinct = con.execute(f"SELECT count(*), count(DISTINCT (chr, position)) FROM {_vpos(cfg)}").fetchone()
    check(n == distinct == EXPECT["variants"], f"trans: (chr, position) is unique over the variant table: {distinct:,} distinct of {n:,} rows "
                                               f"(expected {EXPECT['variants']:,})")

    # 3. af of every trans variant, recomputed from the trans rows, with its rs_number from the variant table
    con.execute(f"""CREATE OR REPLACE TABLE tv AS
        SELECT a.variant_chr, a.position, a.af, a.n_af, v.rs_number, v.position IS NOT NULL AS in_table
        FROM (SELECT variant_chr, position, min(af) AS af, count(DISTINCT af) AS n_af FROM {_scan(sources)} GROUP BY 1, 2) a
        LEFT JOIN {_vpos(cfg)} v ON v.chr = a.variant_chr AND v.position = a.position""")
    n_tv, n_one, n_in = con.execute("SELECT count(*), count(*) FILTER (WHERE n_af = 1), count(*) FILTER (WHERE in_table) FROM tv").fetchone()
    check(n_tv == n_one == n_in == EXPECT["af_variants"],
          f"trans: {n_tv:,} distinct trans variants (expected {EXPECT['af_variants']:,}); {n_one:,} have exactly one af over all their "
          f"trans eQTL and sQTL rows; {n_in:,} are in the variant table")

    # 2. the variants files' two sections
    ptr_dir = pointer_dir(cfg)
    runs_e = {c: (a or 0, b or 0) for c, a, b in con.execute(
        f"SELECT chr, max(var_start + n_var), max(var_off + var_len) FROM {idx} GROUP BY chr").fetchall()}
    tot = dict.fromkeys(("cis", "trans_only", "no_alleles", "shared", "trans_bytes"), 0)
    head_bad, rec_bad, run_bad, cis_af_bad = [], [], [], []
    af_max = 0.0
    vref: dict[str, list] = {"ranges": [], "rows": []}
    for chrom in CHROMS:
        path = pack_paths(cfg, chrom)[0]
        try:
            buf = path.read_bytes()
            count, P, n_cis = _file_header(buf, 1, chrom, path.name)
            want_cis, want_tr = con.execute(f"SELECT count(*) FILTER (WHERE in_cis), count(*) FILTER (WHERE NOT in_cis) FROM {_vpos(cfg, chrom)}").fetchone()
            pages = _read_pages(buf[32:], path.name, n_cis=n_cis)
            cis_pages = -(-n_cis // P)
            cis_end = 32 + pages["page_starts"][cis_pages]
            want_sizes = [min(P, n_cis - s) for s in range(0, n_cis, P)] + [min(P, count - s) for s in range(n_cis, count, P)]
            if not (n_cis == want_cis and count - n_cis == want_tr and len(pages["position"]) == count and pages["page_sizes"] == want_sizes
                    and pages["page_firsts"][cis_pages:cis_pages + 1] in ([n_cis], [])):
                head_bad.append(f"{chrom} (n_cis {n_cis:,} vs {want_cis:,}, trans-only {count - n_cis:,} vs {want_tr:,})")
            tot["cis"] += n_cis
            tot["trans_only"] += count - n_cis
            tot["trans_bytes"] += len(buf) - cis_end
            sl = slice(n_cis, count)
            src = con.execute(f"""SELECT v.position, v.A1, v.A2, coalesce(v.rs_number, 0) AS rs_number, v.match, t.af
                FROM {_vpos(cfg, chrom)} v LEFT JOIN tv t ON t.variant_chr = ? AND t.position = v.position
                WHERE NOT v.in_cis ORDER BY v.position""", [chrom]).fetch_arrow_table()
            s_af = src["af"].to_numpy(zero_copy_only=False).astype(np.float64)
            a1, a2 = src["A1"].to_pylist(), src["A2"].to_pylist()
            code = pages["af_code"][sl].astype(np.int64)
            ok = (np.array_equal(pages["position"][sl], src["position"].to_numpy().astype(np.int64))
                  and np.array_equal(pages["rs_number"][sl].astype(np.int64), src["rs_number"].to_numpy().astype(np.int64))
                  and np.array_equal(pages["match"][sl].astype(np.int64), np.array([MATCH[x] for x in src["match"].to_pylist()], dtype=np.int64))
                  and pages["A1"][sl] == a1 and pages["A2"][sl] == a2
                  and np.array_equal(pages["no_alleles"][sl], np.array([x is None for x in a1], dtype=bool))
                  and bool(np.all(pages["ms"][sl] == 65535)) and bool(np.all(pages["mc"][sl] == 65535))
                  and not np.isnan(s_af).any() and np.array_equal(code, np.rint(s_af * 65534).astype(np.int64)))
            err = np.abs(code / 65534 - s_af)
            af_max = max(af_max, float(err.max(initial=0)))
            if not ok or err.max(initial=0) > TOL["af_tol"]:
                rec_bad.append(chrom)
            tot["no_alleles"] += int(pages["no_alleles"][sl].sum())
            # trans variants that are also cis carry the same af code in the cis section
            sh = con.execute("SELECT position, af FROM tv WHERE variant_chr = ? ORDER BY position", [chrom]).fetch_arrow_table()
            spos, saf = sh["position"].to_numpy().astype(np.int64), sh["af"].to_numpy().astype(np.float64)
            cpos = pages["position"][:n_cis]
            j = np.minimum(np.searchsorted(cpos, spos), max(n_cis - 1, 0))
            in_c = (cpos[j] == spos) if n_cis else np.zeros(spos.size, dtype=bool)
            bad = int(np.sum(pages["af_code"][j[in_c]].astype(np.int64) != np.rint(saf[in_c] * 65534).astype(np.int64)))
            tot["shared"] += int(in_c.sum())
            if bad or int(in_c.sum()) + (count - n_cis) != spos.size:
                cis_af_bad.append(f"{chrom}: {bad} differ, {int(in_c.sum()):,} cis + {count - n_cis:,} trans-only of {spos.size:,}")
            # no run or range reaches the trans-only section
            sp = pq.read_table(ptr_dir / f"sqtl_{chrom}.parquet", columns=["var_start", "n_var"])
            s_max = int(max((a + b for a, b in zip(sp["var_start"].to_pylist(), sp["n_var"].to_pylist())), default=0))
            e_max, e_bytes = runs_e.get(chrom, (0, 0))
            if e_max > n_cis or e_bytes > cis_end or s_max > n_cis:
                run_bad.append(f"{chrom} (eQTL run end {e_max:,}, sQTL run end {s_max:,}, range end {e_bytes:,}; n_cis {n_cis:,} at byte {cis_end:,})")
            # reference: the trans-only page with the most records lacking alleles
            if chrom in REF_VARIANT_CHROMS and count > n_cis:
                tp = list(range(cis_pages, len(pages["page_sizes"])))
                na = [int(pages["no_alleles"][pages["page_firsts"][i]:pages["page_firsts"][i] + pages["page_sizes"][i]].sum()) for i in tp]
                pi = tp[int(np.argmax(na))]
                f0, nn = pages["page_firsts"][pi], pages["page_sizes"][pi]
                rid = len(vref["ranges"])
                vref["ranges"].append({"id": rid, "chrom": chrom, "file": str(path.relative_to(d)), "n_cis": n_cis, "count": count, "page_size": P,
                                       "byte_start": 32 + pages["page_starts"][pi], "byte_end": 32 + pages["page_starts"][pi + 1],
                                       "first_vidx": f0, "n": nn, "no_alleles": max(na)})
                part = src.slice(f0 - n_cis, nn)
                rs = part["rs_number"].to_numpy().astype(np.int64)
                vref["rows"].append(pa.table({
                    "range": pa.array(np.full(nn, rid, dtype=np.int32)), "vidx": pa.array(np.arange(f0, f0 + nn, dtype=np.int64)),
                    "position": part["position"].cast(pa.int32()), "A1": part["A1"].cast(pa.string()), "A2": part["A2"].cast(pa.string()),
                    "no_alleles": pa.array([x is None for x in part["A1"].to_pylist()]),
                    "rs_number": pa.array(rs, mask=rs == 0, type=pa.int64()), "match": part["match"].cast(pa.string()),
                    "af": part["af"].cast(pa.float32()), "af_code": pa.array(np.rint(part["af"].to_numpy().astype(np.float64) * 65534).astype(np.int32))}))
            del pages, buf, src
        except PackError as e:
            check(False, f"{chrom} variants file decodes under SPEC with its trans-only section: {e}")
    check(not head_bad, f"variants files: header n_cis equals the in_cis rows and count - n_cis the NOT in_cis rows on every chromosome; "
                        f"cis and trans-only pages follow the page size and the first trans-only variant starts a page ({head_bad[:3]})")
    check(tot["cis"] == EXPECT["cis"] and tot["trans_only"] == EXPECT["trans_only"] and tot["no_alleles"] == EXPECT["no_alleles"],
          f"variants files hold {tot['cis']:,} cis and {tot['trans_only']:,} trans-only variants, {tot['no_alleles']:,} with flags bit 2 "
          f"(expected {EXPECT['cis']:,}, {EXPECT['trans_only']:,}, {EXPECT['no_alleles']:,}); trans-only sections {tot['trans_bytes']:,} B")
    check(not rec_bad, f"trans-only records equal the variant table (position, rs_number, match, alleles or flags bit 2, null counts), and each af "
                       f"code is rint(af * 65534) of the variant's trans rows: max af error {af_max:.3g} (limit {TOL['af_tol']:.3g}) ({rec_bad})")
    check(not run_bad, f"no eQTL run, union range, or intron run (search_index, sQTL pointers) reaches n_cis ({run_bad[:3]})")
    check(not cis_af_bad and tot["shared"] + tot["trans_only"] == EXPECT["af_variants"],
          f"the af code of all {tot['shared']:,} trans variants that are also cis equals the cis section's code; cis {tot['shared']:,} + "
          f"trans-only {tot['trans_only']:,} = {tot['shared'] + tot['trans_only']:,} trans variants ({cis_af_bad[:3]})")
    if vref["ranges"]:
        (out_dir / "variants.json").write_text(json.dumps({
            "note": "One trans-only page per range: byte_start is the page's first byte in the variants file and byte_end one past its padding. "
                    "variant_rows.arrow holds the records the page must decode to, tagged by range id: A1 and A2 null with no_alleles "
                    "(flags bit 2), rs_number null for 0, ma_samples and ma_count null, af_code = rint(af * 65534).",
            "ranges": vref["ranges"]}, indent=1))
        _write_arrow(pa.concat_tables(vref["rows"]), out_dir / "variant_rows.arrow")

    # 4. the trans pack, every gene
    si = con.execute(f"SELECT gene_id, symbol, chr, tss, gene_version, trans_off, trans_len FROM {idx} ORDER BY chr, tss, gene_id").fetch_arrow_table().to_pylist()
    ver_bad = con.execute(f"""SELECT count(*) FROM {idx} s JOIN '{cfg.tables / 'genes.parquet'}' g USING (gene_id)
        WHERE s.gene_version IS DISTINCT FROM split_part(g.gene_id_version, '.', 2)::INTEGER""").fetchone()[0]
    per_gene = con.execute(f"""SELECT gene_id, any_value(gene_chr), count(*) FILTER (WHERE qtl_type = 'e'), count(*) FILTER (WHERE qtl_type = 's')
        FROM {_scan(sources)} GROUP BY 1 ORDER BY 1""").fetchall()
    has_rows = {g for g, *_ in per_gene}
    ptr_bad = sum((g["gene_id"] in has_rows) != (g["trans_off"] is not None) or (g["trans_off"] is None) != (g["trans_len"] is None) for g in si)
    check(ver_bad == 0, f"search_index gene_version equals the number after the dot in genes.gene_id_version for all {len(si):,} genes ({ver_bad} differ)")
    check(ptr_bad == 0, f"search_index trans_off and trans_len are set for exactly the {len(has_rows):,} genes with trans rows ({ptr_bad} rows disagree)")
    rng = np.random.default_rng(20260915)
    reasons = {g: why for g, why in REF_GENES.items() if g in has_rows}
    for why, cands in (("chrM", [g for g, c, e, s in per_gene if c == "chrM"]), ("chrX", [g for g, c, e, s in per_gene if c == "chrX"]),
                       ("sQTL only", [g for g, c, e, s in per_gene if e == 0]), ("eQTL only", [g for g, c, e, s in per_gene if s == 0])):
        cands = [g for g in cands if g not in reasons]
        if cands:
            reasons[str(rng.choice(cands))] = why
    for g in rng.choice(sorted(has_rows - set(reasons)), REF_RANDOM, replace=False).tolist():
        reasons[g] = "random"

    agg = dict.fromkeys(("rows", "e", "s", "bytes", "genes"), 0)
    worst = {"nlp": 0.0, "nlp_ratio": 0.0, "beta": 0.0, "beta_ratio": 0.0, "beta_se_rel": 0.0, "r2": 0.0, "af": 0.0, "min_beta_code": 1 << 16}
    files, ref_genes, ref_rows = {}, [], []
    for chrom in TRANS_CHROMS:
        src_path, path = trans_source(cfg, chrom), trans_path(cfg, chrom)
        if not src_path.exists():
            continue
        try:
            if not path.exists():
                raise PackError(f"{path.relative_to(d)} is missing")
            buf = path.read_bytes()
            files[chrom] = str(path.relative_to(d))
            count, page, _ = _file_header(buf, 6, chrom, path.name)
            genes = [g for g in si if g["chr"] == chrom and g["trans_off"] is not None]
            off, contiguous = 32, True
            for g in genes:
                contiguous &= g["trans_off"] == off
                off += g["trans_len"]
            check(contiguous and off == len(buf) and count == len(genes) and page == 0,
                  f"{chrom}: trans pack is kind 6 with {count:,} frames for {len(genes):,} genes with trans rows, back to back in "
                  f"search_index order from byte 32 to the file end ({off:,} of {len(buf):,} B): 32 + sum(trans_len)")
            frames = [_read_trans_frame(buf[g["trans_off"]:g["trans_off"] + g["trans_len"]], f"{path.name} {g['gene_id']}") for g in genes]
            n_rows = np.array([f["n_e"] + f["n_s"] for f in frames], dtype=np.int64)
            is_s = np.concatenate([np.r_[np.zeros(f["n_e"], dtype=bool), np.ones(f["n_s"], dtype=bool)] for f in frames])
            intr = np.concatenate([np.r_[np.zeros((f["n_e"], 4), dtype=np.int64), f["table"][f["intron"]]] if f["n_s"]
                                   else np.zeros((f["n_e"], 4), dtype=np.int64) for f in frames])
            D = {k: np.concatenate([f[k] for f in frames]) for k in ("position", "rs_number", "af_code", "nlp_code", "beta_code", "variant_chr")}
            nlp_max = np.repeat([f["nlp_max"] for f in frames], n_rows)
            beta_max = np.repeat([f["beta_max"] for f in frames], n_rows)
            S = con.execute(f"""WITH r AS (
                    SELECT t.*, {INTRON_SQL.format(t='t')}, {CHR_CODE_SQL.format(c='t.variant_chr')} AS chr_code
                    FROM read_parquet('{src_path}', hive_partitioning = false) t)
                SELECT r.gene_id, r.qtl_type = 's' AS is_s, coalesce(r.intron_start, 0)::BIGINT AS istart, coalesce(r.intron_end, 0)::BIGINT AS iend,
                       coalesce(r.cluster, 0)::BIGINT AS clu, CASE WHEN r.strand = '-' THEN 1 ELSE 0 END AS strand_code, r.chr_code,
                       r.position, coalesce(tv.rs_number, 0) AS rs_number, r.af, r.pval, r.beta, r.beta_se, r.r2,
                       CASE WHEN r.qtl_type = 'e' THEN r.phenotype_id = r.gene_id
                            ELSE r.phenotype_id = s.chr || ':' || r.intron_start::VARCHAR || ':' || r.intron_end::VARCHAR || ':clu_'
                                 || r.cluster::VARCHAR || '_' || r.strand || ':' || r.gene_id || '.' || s.gene_version::VARCHAR END AS id_ok,
                       tv.position IS NOT NULL AS in_tv
                FROM r JOIN {idx} s USING (gene_id)
                LEFT JOIN tv ON tv.variant_chr = r.variant_chr AND tv.position = r.position
                ORDER BY s.tss, r.gene_id, is_s, istart, iend, clu, strand_code, r.chr_code, r.position""").fetch_arrow_table()
            sg = S["gene_id"].to_numpy(zero_copy_only=False)
            sstart = np.r_[0, np.flatnonzero(sg[1:] != sg[:-1]) + 1] if len(sg) else np.zeros(0, dtype=np.int64)
            s_counts = np.diff(np.r_[sstart, len(sg)])
            col = {k: S[k].to_numpy(zero_copy_only=False) for k in ("is_s", "istart", "iend", "clu", "strand_code", "chr_code", "position",
                                                                   "rs_number", "af", "pval", "beta", "beta_se", "r2", "id_ok", "in_tv")}
            exact = {
                "genes and row counts": [str(x) for x in sg[sstart]] == [g["gene_id"] for g in genes] and np.array_equal(s_counts, n_rows),
                "qtl_type": np.array_equal(is_s, col["is_s"].astype(bool)),
                "intron": all(np.array_equal(intr[:, i], col[k].astype(np.int64)) for i, k in enumerate(("istart", "iend", "clu", "strand_code"))),
                "variant_chr, position": np.array_equal(D["variant_chr"], col["chr_code"].astype(np.int64)) and np.array_equal(D["position"], col["position"].astype(np.int64)),
                "rs_number": np.array_equal(D["rs_number"], col["rs_number"].astype(np.int64)),
                "phenotype_id rebuilt from the frame and search_index": bool(np.all(col["id_ok"])) and bool(np.all(col["in_tv"])),
            } if len(sg) == int(n_rows.sum()) else {"row count": False}
            bad = [k for k, v in exact.items() if not v]
            check(not bad, f"{chrom}: {int(n_rows.sum()):,} decoded trans rows of {len(genes):,} genes equal the source rows in frame order: "
                           f"qtl_type, phenotype_id (rebuilt), variant_chr, position, rs_number ({bad})")
            if bad:
                continue
            dofs = np.where(is_s, float(dof["s"]), float(dof["e"]))
            nlp = D["nlp_code"] * (nlp_max / 65533)
            beta = D["beta_code"] * (beta_max / 32767)
            t = np.maximum(-stdtrit(dofs, np.power(10.0, -nlp) / 2), 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                se = np.abs(beta) / t
            r2 = t * t / (t * t + dofs)
            b_src = col["beta"].astype(np.float64)
            nlp_err = np.abs(nlp + np.log10(col["pval"].astype(np.float64)))
            nlp_lim = nlp_max / TOL["nlp_half_step_divisor"] + TOL["nlp_slack"]
            beta_err = np.abs(beta - b_src)
            beta_lim = beta_max / TOL["beta_half_step_divisor"] + TOL["beta_rel_f32"] * np.abs(b_src)
            se_src = col["beta_se"].astype(np.float64)
            se_rel = np.abs(se - se_src) / se_src
            r2_err = np.abs(r2 - col["r2"].astype(np.float64))
            af_src = col["af"].astype(np.float64)
            af_err = np.abs(D["af_code"] / 65534 - af_src)
            af_exact = np.array_equal(D["af_code"], np.rint(af_src * 65534).astype(np.int64))
            vals = {"nlp": float(nlp_err.max()), "nlp_ratio": float((nlp_err / nlp_lim).max()), "beta": float(beta_err.max()),
                    "beta_ratio": float((beta_err / beta_lim).max()), "beta_se_rel": float(np.nanmax(se_rel)) if np.all(np.isfinite(se_rel)) else math.inf,
                    "r2": float(r2_err.max()), "af": float(af_err.max()), "min_beta_code": int(np.abs(D["beta_code"]).min())}
            check(vals["nlp_ratio"] <= 1 and vals["beta_ratio"] <= 1 and vals["beta_se_rel"] <= TOL["beta_se_rel"] and vals["r2"] <= TOL["r2_abs"]
                  and af_exact and vals["af"] <= TOL["af_tol"],
                  f"{chrom}: trans values within SPEC section 9: -log10 p error {vals['nlp']:.3g} (error / limit {vals['nlp_ratio']:.4f}), beta error "
                  f"{vals['beta']:.3g} ({vals['beta_ratio']:.4f}), derived beta_se relative error {vals['beta_se_rel']:.3g} (limit {TOL['beta_se_rel']}), "
                  f"derived r2 error {vals['r2']:.3g} (limit {TOL['r2_abs']}), af code rint(af * 65534) (error {vals['af']:.3g}); "
                  f"smallest |beta code| {vals['min_beta_code']}")
            for k, v in vals.items():
                worst[k] = min(worst[k], v) if k == "min_beta_code" else max(worst[k], v)
            agg["rows"] += int(n_rows.sum())
            agg["e"] += int((~is_s).sum())
            agg["s"] += int(is_s.sum())
            agg["bytes"] += len(buf)
            agg["genes"] += len(genes)
            # reference rows for the browser decoder check
            sel = [i for i, g in enumerate(genes) if g["gene_id"] in reasons]
            if sel:
                ids = ", ".join(f"'{genes[i]['gene_id']}'" for i in sel)
                ref_rows.append(con.execute(f"""WITH r AS (
                        SELECT t.*, {INTRON_SQL.format(t='t')}, {CHR_CODE_SQL.format(c='t.variant_chr')} AS chr_code
                        FROM read_parquet('{src_path}', hive_partitioning = false) t WHERE t.gene_id IN ({ids}))
                    SELECT r.gene_id, r.symbol, r.gene_chr, r.gene_tss, r.qtl_type, r.phenotype_id, r.variant_chr, r.position,
                           tv.rs_number::BIGINT AS rs_number, CASE WHEN tv.rs_number IS NULL THEN NULL ELSE 'rs' || tv.rs_number::VARCHAR END AS rsid,
                           r.af, r.pval, r.beta, r.beta_se, r.r2
                    FROM r JOIN {idx} s USING (gene_id)
                    LEFT JOIN tv ON tv.variant_chr = r.variant_chr AND tv.position = r.position
                    ORDER BY s.tss, r.gene_id, r.qtl_type, coalesce(r.intron_start, 0), coalesce(r.intron_end, 0), coalesce(r.cluster, 0),
                             CASE WHEN r.strand = '-' THEN 1 ELSE 0 END, r.chr_code, r.position""").fetch_arrow_table())
                for i in sel:
                    g, f = genes[i], frames[i]
                    ref_genes.append({"gene_id": g["gene_id"], "symbol": g["symbol"], "chr": chrom, "file": files[chrom],
                                      "gene_version": g["gene_version"], "trans_off": g["trans_off"], "trans_len": g["trans_len"],
                                      "n_e": f["n_e"], "n_s": f["n_s"], "k": f["k"], "nlp_max": f["nlp_max"], "beta_max": f["beta_max"],
                                      "reason": reasons[g["gene_id"]]})
            del frames, D, S, col, buf
        except PackError as e:
            check(False, f"{chrom} trans pack decodes under SPEC: {e}")
    check(agg["rows"] == EXPECT["rows"] and agg["e"] == EXPECT["eqtl_rows"] and agg["s"] == EXPECT["sqtl_rows"],
          f"trans packs hold {agg['rows']:,} rows ({agg['e']:,} eQTL, {agg['s']:,} sQTL) of {agg['genes']:,} genes in {len(files)} files, "
          f"{agg['bytes']:,} B (expected {EXPECT['rows']:,}, {EXPECT['eqtl_rows']:,}, {EXPECT['sqtl_rows']:,})")
    log(f"trans round trip over every row: worst -log10 p error {worst['nlp']:.3g} (error / half step {worst['nlp_ratio']:.4f}), beta error "
        f"{worst['beta']:.3g} (error / limit {worst['beta_ratio']:.4f}), beta_se relative error {worst['beta_se_rel']:.3g}, r2 error "
        f"{worst['r2']:.3g}, af error {worst['af']:.3g}; smallest |beta code| {worst['min_beta_code']}")

    # 5. rsID text, log only (the variant-page stage fixes variants_rsid)
    pos_bad, rows_bad = con.execute(f"""SELECT count(DISTINCT (t.variant_chr, t.position)), count(*) FROM {_scan(sources)} t
        JOIN tv ON tv.variant_chr = t.variant_chr AND tv.position = t.position
        WHERE t.rsid IS DISTINCT FROM CASE WHEN tv.rs_number IS NULL THEN NULL ELSE 'rs' || tv.rs_number::VARCHAR END""").fetchone()
    log(f"trans rsID text (not a failure; the variant-page stage fixes the source): {pos_bad:,} trans positions ({rows_bad:,} rows) have a trans_pairs rsid "
        f"other than 'rs' + the variant's rs_number (measured 9,940 positions); the trans pack stores rs_number only")

    # 6. reference files
    if ref_rows:
        rows = pa.concat_tables(ref_rows).cast(pa.schema([
            ("gene_id", pa.string()), ("symbol", pa.string()), ("gene_chr", pa.string()), ("gene_tss", pa.int64()), ("qtl_type", pa.string()),
            ("phenotype_id", pa.string()), ("variant_chr", pa.string()), ("position", pa.int32()), ("rs_number", pa.int64()), ("rsid", pa.string()),
            ("af", pa.float32()), ("pval", pa.float64()), ("beta", pa.float32()), ("beta_se", pa.float32()), ("r2", pa.float32())]))
        _write_arrow(rows, out_dir / "rows.arrow")
        (out_dir / "frames.json").write_text(json.dumps({
            # key order matches manifest.json, which is written with sort_keys: `npm run pack-check`
            # compares the two with JSON.stringify, which is order-sensitive
            "dof": {"eqtl": dof["e"], "sqtl": dof["s"]}, "files": dict(sorted(files.items())), "tolerances": TOL,
            "note": "Each gene's frame is bytes trans_off..trans_off + trans_len - 1 of its file. rows.arrow holds the gene's source rows "
                    "(trans_pairs, with rs_number from the variant table and rsid = 'rs' + rs_number) in frame order: eQTL rows, then sQTL rows "
                    "by intron (start, end, cluster, strand), each run by variant chromosome (chr1..chr22, chrX) then position. Limits per row: "
                    "|-log10 p error| <= nlp_max / nlp_half_step_divisor + nlp_slack; |beta error| <= beta_max / beta_half_step_divisor + "
                    "beta_rel_f32 * |beta|; |beta_se error| / beta_se <= beta_se_rel; |r2 error| <= r2_abs; |af error| <= af_tol.",
            "genes": ref_genes}, indent=1))
        log(f"pack validate: trans reference files in {out_dir.relative_to(d)}/: {len(ref_genes)} genes "
            f"({', '.join(sorted(set(reasons.values())))}), {rows.num_rows:,} rows; {len(vref['ranges'])} trans-only page ranges")
    con.execute("DROP TABLE IF EXISTS tv")
    log(f"trans validate: {(time.time() - t_start) / 60:.1f} min")
