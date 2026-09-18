"""Steps 10-11: manifest and validation."""
import json
import math
import mmap
import subprocess

import yaml

from . import packfmt, steps_pack, steps_pack_trans, steps_pack_variant
from .common import (NAME, PACKS_METADATA_KEY, ROOT, SEARCH_INDEX_NAME, Config, addressed_files, connect,
                     digests, log, search_index_packs, search_index_path, variants_sql)


def _replaces(cfg: Config, key: str) -> list[str]:
    """The old bucket prefixes one published file takes over, from `replaces` in the config. A
    per-chromosome key fills `{chr}`; a single file takes its list as it is."""
    kind, _, chrom = key.partition("/")
    return [pre.format(chr=chrom) for pre in cfg["replaces"].get(kind, [])]


def _immutable_block(cfg: Config) -> dict:
    """One entry per published file, keyed by its bucket path: what `upload.py` needs to send it and
    what `validate` needs to prove it is the file the manifest names (SPEC section 3)."""
    out = {}
    for key, path in sorted(addressed_files(cfg).items()):
        sha, md5 = digests(path)
        if NAME.match(path.name)["sha"] != sha[:16]:
            raise ValueError(f"manifest: {path.name} does not hash to its own name (sha256 starts {sha[:16]}); "
                             f"re-run the step that writes {key}")
        out[str(path.relative_to(cfg.derived))] = {
            "key": key, "bytes": path.stat().st_size, "sha256": sha, "md5": md5, "replaces": _replaces(cfg, key)}
    return out


def _precision(cfg: Config) -> dict:
    """The largest rounding error a reader can see, per SPEC section 9. The quantized maximums come
    from the block headers of every eQTL and sQTL pack; the slope error comes from `packcheck
    roundtrip`, which compares rebuilt slopes against the source rows genome-wide."""
    rt_path = cfg.derived / cfg["packs"]["checks_dir"] / "roundtrip_genome.json"
    if not rt_path.exists():
        raise FileNotFoundError(f"manifest: {rt_path.relative_to(cfg.derived)} is missing; "
                                "run `python -m pipeline packcheck roundtrip` first")
    rt = json.loads(rt_path.read_text())
    files = addressed_files(cfg)
    out = {"af_max_error": 0.5 / packfmt.AF_MAXQ}
    for kind in ("eqtl", "sqtl"):
        nlp_err = se_err = 0.0
        for key in sorted(k for k in files if k.startswith(f"{kind}/")):
            # blocks sit back to back after the 32-byte file header; mmap so a 200 MB pack is not copied
            with open(files[key], "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                off = packfmt.FILE_HEADER_LEN
                while off < len(mm):
                    h = packfmt.parse_block_header(mm, off)
                    nlp_err = max(nlp_err, h["nlp_max"] / (2 * packfmt.NLP_MAXQ))
                    se_err = max(se_err, math.expm1((h["lse_max"] - h["lse_min"]) / (2 * packfmt.SE_MAXQ)))
                    off += h["blk_len"]
        tally = rt["types"][kind]["tally"]["max_slope_err_over_se"]
        out[kind] = {"neglog10p_max_error": nlp_err, "slope_se_max_rel_error": se_err,
                     "slope_max_error_over_se": max(float(x) for x in tally)}
    out["source"] = f"{rt_path.relative_to(cfg.derived)} (generated {rt.get('generated')})"
    return out


def _check_index_packs(cfg: Config, immutable: dict) -> dict:
    """`search_index` records the SHA-256 of every pack its byte offsets reach. A pack rebuilt after
    the index would move those offsets silently, so the manifest refuses to name both."""
    recorded = search_index_packs(cfg)
    if not recorded:
        raise ValueError(f"manifest: search_index carries no `{PACKS_METADATA_KEY.decode()}` metadata; "
                         "run `build --step search_index --force`")
    by_key = {v["key"]: v for v in immutable.values()}
    stale = sorted(k for k, sha in recorded.items() if k in by_key and by_key[k]["sha256"] != sha)
    gone = sorted(k for k in recorded if k not in by_key)
    extra = sorted(k for k in by_key if k.split("/")[0] in steps_pack.INDEX_PACK_KINDS and k not in recorded)
    if stale or gone:
        raise ValueError(f"manifest: search_index was built from older packs; rebuild it "
                         f"(`build --step search_index --force`). Changed: {stale[:5]}; missing: {gone[:5]}")
    if extra:
        raise ValueError(f"manifest: search_index points into {extra[:5]} but does not record them; "
                         "run `build --step search_index --force`")
    return recorded


def manifest(cfg: Config) -> None:
    import datetime as dt
    con = connect(cfg)
    sources = yaml.safe_load((cfg.raw / "sources.yaml").read_text())
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT).stdout.strip() or None
    except FileNotFoundError:
        commit = None
    # no `tables` block: every table is a build intermediate under `_tables/` now, and the browser
    # reads only packs, the search index, and the two small JSON assets
    sig_col, thr = cfg["sig_column"], cfg["sig_threshold"]
    V = variants_sql(cfg)
    counts = {
        "egenes": con.execute(f"SELECT count(*) FROM '{cfg.tables / 'genes.parquet'}' WHERE is_egene").fetchone()[0],
        "genes_tested": con.execute(f"SELECT count(*) FROM '{cfg.tables / 'genes.parquet'}' WHERE tested").fetchone()[0],
        "sqtl_sig_phenotypes": con.execute(f"SELECT count(*) FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0],
        "sqtl_sig_genes": con.execute(f"SELECT count(DISTINCT gene_id) FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0],
        # match breakdown for cis variants; trans-only rows include allele-less positions that can only match by position
        "rsid_match": dict(con.execute(f"SELECT match, count(*) FROM {V} WHERE in_cis GROUP BY 1").fetchall()),
        "variants_cis": con.execute(f"SELECT count(*) FROM {V} WHERE in_cis").fetchone()[0],
        "variants_trans_only": con.execute(f"SELECT count(*) FROM {V} WHERE NOT in_cis").fetchone()[0],
        "splice_phenotypes_tested": con.execute(f"SELECT count(*) FROM '{cfg.tables / 'splice_phenotypes.parquet'}'").fetchone()[0],
        "gwas_variants": json.loads((cfg.derived / "gwas_dcm.json").read_text())["variants"],
        "trans_pairs": con.execute(f"SELECT count(*) FROM read_parquet('{cfg.tables}/trans/*/*.parquet', hive_partitioning = false)").fetchone()[0],
        "trans_variants": con.execute(f"SELECT count(DISTINCT (variant_chr, position)) FROM read_parquet('{cfg.tables}/trans/*/*.parquet', hive_partitioning = false)").fetchone()[0],
    }
    pack_stats, sqtl_stats, var_stats = [], [], []
    pack_files = {"variants": {}, "eqtl": {}, "sqtl": {}}
    pack_bytes = {"variants": 0, "eqtl": 0, "sqtl": 0}
    ptr = steps_pack.pointer_dir(cfg)
    for chrom in steps_pack.CHROMS:
        vp, ep = steps_pack.pack_paths(cfg, chrom)
        sp = steps_pack.sqtl_path(cfg, chrom)
        stat, sstat, vstat = ptr / f"eqtl_{chrom}.json", ptr / f"sqtl_{chrom}.json", ptr / f"variants_trans_{chrom}.json"
        missing = [str(p) for p in (vp, ep, sp, stat, sstat, vstat) if not p.exists()]
        if missing:
            raise FileNotFoundError(f"manifest: missing pack output for {chrom}: {missing}")
        for kind, path in (("variants", vp), ("eqtl", ep), ("sqtl", sp)):
            pack_files[kind][chrom] = str(path.relative_to(cfg.derived))
            pack_bytes[kind] += path.stat().st_size
        pack_stats.append(json.loads(stat.read_text()))
        sqtl_stats.append(json.loads(sstat.read_text()))
        var_stats.append(json.loads(vstat.read_text()))
    if sum(x["cis"] for x in var_stats) != sum(x["variants"] for x in pack_stats):
        raise ValueError("manifest: pack_variants_trans's cis counts differ from pack_eqtl's; re-run pack_variants_trans")
    # trans: one file per gene chromosome with trans rows (chr1-22, chrX, chrM)
    trans_stats = []
    pack_files["trans"] = {}
    for chrom in steps_pack_trans.TRANS_CHROMS:
        if not steps_pack_trans.trans_source(cfg, chrom).exists():
            continue
        tp, tstat = steps_pack_trans.trans_path(cfg, chrom), ptr / f"trans_{chrom}.json"
        missing = [str(p) for p in (tp, tstat) if not p.exists()]
        if missing:
            raise FileNotFoundError(f"manifest: missing trans pack output for {chrom}: {missing}")
        pack_files["trans"][chrom] = str(tp.relative_to(cfg.derived))
        pack_bytes["trans"] = pack_bytes.get("trans", 0) + tp.stat().st_size
        trans_stats.append(json.loads(tstat.read_text()))
    # hits: one pack per variant chromosome, plus the rsID index and the startup file (SPEC sections 13 to 15)
    hits_stats = []
    pack_files["hits"] = {}
    for chrom in steps_pack.CHROMS:
        hp, hstat = steps_pack_variant.hits_path(cfg, chrom), ptr / f"hits_{chrom}.json"
        missing = [str(p) for p in (hp, hstat) if not p.exists()]
        if missing:
            raise FileNotFoundError(f"manifest: missing hits pack output for {chrom}: {missing}")
        pack_files["hits"][chrom] = str(hp.relative_to(cfg.derived))
        pack_bytes["hits"] = pack_bytes.get("hits", 0) + hp.stat().st_size
        hits_stats.append(json.loads(hstat.read_text()))
    rsid_index = steps_pack_variant.rsid_index_path(cfg)
    variant_index = steps_pack_variant.variant_index_path(cfg)
    vxstat_path = ptr / "variant_index.json"
    missing = [str(p) for p in (rsid_index, variant_index, vxstat_path) if not p.exists()]
    if missing:
        raise FileNotFoundError(f"manifest: missing variant-index output: {missing}")
    vxstat = json.loads(vxstat_path.read_text())
    pack_bytes["rsid_index"] = rsid_index.stat().st_size
    pack_bytes["variant_index"] = variant_index.stat().st_size

    # GWAS: one file per chromosome with rows (chr1-22) and the startup index
    gstat = json.loads((ptr / "gwas.json").read_text())
    pack_files["gwas"] = {c: str(steps_pack.gwas_path(cfg, c).relative_to(cfg.derived)) for c in steps_pack.CHROMS if steps_pack.gwas_path(cfg, c).exists()}
    if list(pack_files["gwas"]) != list(gstat["rows_by_chrom"]):
        raise FileNotFoundError(f"manifest: GWAS pack files {list(pack_files['gwas'])} differ from pack_gwas's chromosomes {list(gstat['rows_by_chrom'])}")
    gwas_index = steps_pack.gwas_index_path(cfg)
    pack_bytes["gwas"] = sum((cfg.derived / p).stat().st_size for p in pack_files["gwas"].values())
    pack_bytes["gwas_index"] = gwas_index.stat().st_size
    immutable = _immutable_block(cfg)
    _check_index_packs(cfg, immutable)
    packs = {
        "format": "qtlb", "version": 0,
        "dof": {"eqtl": int(cfg["packs"]["dof"]["eqtl"]), "sqtl": int(cfg["packs"]["dof"]["sqtl"])},
        "variant_page_size": int(cfg["packs"]["variant_page_size"]),
        "variant_page_codec": cfg["packs"]["variant_page_codec"],
        "search_index": str(search_index_path(cfg).relative_to(cfg.derived)),
        "gwas_index": str(gwas_index.relative_to(cfg.derived)), "gwas_block_rows": int(cfg["packs"]["gwas_block_rows"]),
        "rsid_index": str(rsid_index.relative_to(cfg.derived)), "rsid_block_records": int(cfg["packs"]["rsid_block_records"]),
        "variant_index": str(variant_index.relative_to(cfg.derived)), "hits_frame_variants": int(cfg["packs"]["hits_frame_variants"]),
        "files": pack_files, "bytes": pack_bytes,
        "counts": {"variants": {"cis": sum(x["cis"] for x in var_stats), "trans_only": sum(x["trans_only"] for x in var_stats)},
                   "eqtl": {"genes": sum(x["genes"] for x in pack_stats), "rows": sum(x["rows"] for x in pack_stats),
                            "memberships": sum(x["memberships"] for x in pack_stats)},
                   "sqtl": {"introns": sum(x["introns"] for x in sqtl_stats), "rows": sum(x["rows"] for x in sqtl_stats),
                            "memberships": sum(x["memberships"] for x in sqtl_stats)},
                   "gwas": {"rows": sum(gstat["rows_by_chrom"].values())},
                   "trans": {"genes": sum(x["genes"] for x in trans_stats), "rows": sum(x["rows"] for x in trans_stats),
                             "eqtl_rows": sum(x["eqtl_rows"] for x in trans_stats), "sqtl_rows": sum(x["sqtl_rows"] for x in trans_stats)},
                   "hits": {"rows": sum(x["rows"] for x in hits_stats), "frames": sum(x["frames"] for x in hits_stats),
                            **{k: sum(x["by_kind"][k] for x in hits_stats) for k in steps_pack_variant.KIND_COUNT_KEY}},
                   "rsid_index": {"records": vxstat["rsid_records"]}}}
    out = {
        "built": dt.datetime.now().isoformat(timespec="seconds"),
        "pipeline_commit": commit,
        "significance_rule": f"{sig_col} < {thr}",
        "paper_counts": cfg["paper_counts"],
        "counts": counts,
        "sources": {s["name"]: {"version": s.get("version"), "description": s.get("description")} for s in sources["sources"]},
        "packs": packs,
        # every content-addressed file: what upload.py sends, and what it replaces in the bucket
        "immutable": immutable,
        # SPEC section 9: how far a value the browser shows can sit from the source value
        "precision": _precision(cfg),
        # small JSON files the app fetches directly, without the query engine
        "assets": {n: {"path": n, "bytes": (cfg.derived / n).stat().st_size} for n in ("coloc_loci.json", "gwas_dcm_bins.json") if (cfg.derived / n).exists()},
        "gwas_dcm": json.loads((cfg.derived / "gwas_dcm.json").read_text()),
    }
    (cfg.derived / "manifest.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    log(f"manifest: eGenes={counts['egenes']}, sQTL sig={counts['sqtl_sig_phenotypes']}, "
        f"introns tested={counts['splice_phenotypes_tested']:,}, trans pairs={counts['trans_pairs']:,}")
    log(f"manifest: {len(immutable)} content-addressed files, {sum(v['bytes'] for v in immutable.values()):,} B")


def _validate_immutable(cfg: Config, check) -> None:
    """SPEC section 3: `immutable/` holds one file per logical key, each named by its own SHA-256;
    manifest.json names those same files; and nothing is left at an unhashed path."""
    mpath = cfg.derived / "manifest.json"
    if not mpath.exists():
        check(False, "manifest.json exists (run `build --step manifest`)")
        return
    man = json.loads(mpath.read_text())

    # every file in immutable/ parses as <stem>.<sha16>.<ext>, and no key has two builds
    strays, seen = [], {}
    for f in sorted(cfg.immutable.glob("*")):
        m = NAME.match(f.name)
        if not m:
            strays.append(f.name)
            continue
        seen.setdefault(m["stem"].replace(".", "/"), []).append(f.name)
    twice = {k: v for k, v in seen.items() if len(v) > 1}
    check(not strays and not twice,
          f"immutable/ holds {len(seen)} keys, one file each"
          + (f" (strays {strays[:3]}; two builds of {list(twice)[:3]})" if strays or twice else ""))

    # no unhashed copy survives: neither the old packs/ tree nor a top-level search_index
    leftovers = [str(p.relative_to(cfg.derived)) for p in
                 [cfg.derived / "packs", cfg.derived / SEARCH_INDEX_NAME, cfg.derived / "search_index.parquet"]
                 if p.exists()]
    check(not leftovers, f"no unhashed pack or index left in data/derived ({leftovers})")

    # the manifest's immutable block: every entry exists with the recorded size and SHA-256
    imm = man.get("immutable") or {}
    bad_missing, bad_bytes, bad_sha = [], [], []
    for path, rec in sorted(imm.items()):
        f = cfg.derived / path
        if not f.exists():
            bad_missing.append(path)
            continue
        if f.stat().st_size != rec["bytes"]:
            bad_bytes.append(path)
            continue
        sha, md5 = digests(f)
        if sha != rec["sha256"] or md5 != rec["md5"] or sha[:16] != NAME.match(f.name)["sha"]:
            bad_sha.append(path)
    check(bool(imm) and not (bad_missing or bad_bytes or bad_sha),
          f"manifest immutable: {len(imm)} files exist with the recorded bytes, sha256 and md5"
          + (f" (missing {bad_missing[:3]}, size {bad_bytes[:3]}, digest {bad_sha[:3]})"
             if bad_missing or bad_bytes or bad_sha else ""))

    # the paths the browser follows all resolve to an immutable entry
    packs = man.get("packs") or {}
    named = [packs.get(k) for k in ("search_index", "gwas_index", "rsid_index", "variant_index")]
    named += [p for kind in (packs.get("files") or {}).values() for p in kind.values()]
    unknown = sorted({p for p in named if p and p not in imm})
    check(named and not unknown, f"every packs path is an immutable entry ({len(named)} paths"
                                 + (f", unknown {unknown[:3]}" if unknown else "") + ")")

    # search_index names the packs it indexes, and they have not been rebuilt since
    recorded = search_index_packs(cfg)
    by_key = {v["key"]: v for v in imm.values()}
    stale = sorted(k for k, sha in recorded.items() if by_key.get(k, {}).get("sha256") != sha)
    missing_kind = sorted({k for k in by_key if k.split("/")[0] in steps_pack.INDEX_PACK_KINDS} - set(recorded))
    check(bool(recorded) and not stale and not missing_kind,
          f"search_index records the sha256 of all {len(recorded)} packs it points into"
          + (f" (stale {stale[:3]}, unrecorded {missing_kind[:3]})" if stale or missing_kind else ""))

    # no plain-path asset sits under a prefix a content-addressed file has taken over
    prefixes = [pre for rec in imm.values() for pre in rec["replaces"]] + list(cfg["replaces"].get("_retired") or [])
    clash = sorted(a["path"] for a in (man.get("assets") or {}).values()
                   if any(a["path"].startswith(pre.rstrip("*")) for pre in prefixes))
    check(not clash, f"no asset path sits under a replaced bucket prefix ({clash[:3]})")

    # the manifest's rounding maximums are at least what a round trip actually measured
    prec = man.get("precision") or {}
    rt = cfg.derived / cfg["packs"]["checks_dir"] / "roundtrip_genome.json"
    low = []
    if prec and rt.exists():
        tally = json.loads(rt.read_text())["types"]
        for kind in ("eqtl", "sqtl"):
            got = max(float(x) for x in tally[kind]["tally"]["max_slope_err_over_se"])
            if prec.get(kind, {}).get("slope_max_error_over_se", 0) < got:
                low.append(f"{kind}: {prec[kind]['slope_max_error_over_se']:.3g} < {got:.3g}")
    check(bool(prec) and not low,
          "manifest precision is at least the round trip's measured slope error" + (f" ({low})" if low else ""))


def validate(cfg: Config) -> None:
    con = connect(cfg)
    fails = []

    def check(ok: bool, msg: str):
        log(("PASS " if ok else "FAIL ") + msg)
        if not ok:
            fails.append(msg)

    # 1. row counts per nominal partition equal raw (sQTL: raw restricted to significant introns when configured)
    sig_only = cfg["sqtl_nominal"] == "significant"
    sp = cfg.tables / "splice_phenotypes.parquet"
    for qtl, raw_dir, out_dir, pat in [
        ("e", "cis_eQTL_nominal", "cis_eqtl_nominal", "topchef_{c}_MaxPC70.cis_qtl_pairs.{c}.parquet"),
        ("s", "cis_sQTL_nominal", "cis_sqtl_nominal", "topchefSplice_{c}_MaxPC25.cis_qtl_pairs.{c}.parquet"),
    ]:
        bad = []
        for c in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
            r = cfg.raw_dir(raw_dir) / pat.format(c=c)
            d = cfg.tables / out_dir / f"chr={c}" / "bin=*" / "data.parquet"
            if not (r.exists() and list(d.parent.parent.glob("bin=*/data.parquet"))):
                bad.append(f"{c}:missing")
                continue
            if qtl == "s" and sig_only:
                nr = con.execute(f"SELECT count(*) FROM '{r}' SEMI JOIN (SELECT phenotype_id FROM '{sp}' WHERE is_sqtl) USING (phenotype_id)").fetchone()[0]
            else:
                nr = con.execute(f"SELECT count(*) FROM '{r}'").fetchone()[0]
            nd = con.execute(f"SELECT count(*) FROM '{d}'").fetchone()[0]
            if nr != nd:
                bad.append(f"{c}:{nr}!={nd}")
        what = "raw significant introns" if qtl == "s" and sig_only else "raw"
        check(not bad, f"{out_dir} row counts match {what}" + (f" ({', '.join(bad)})" if bad else ""))

    # 2. every nominal gene appears in genes as tested
    n = con.execute(f"""
        SELECT count(*) FROM (SELECT DISTINCT gene_id FROM read_parquet('{cfg.tables}/cis_eqtl_nominal/*/*/*.parquet')) x
        LEFT JOIN '{cfg.tables / 'genes.parquet'}' g USING (gene_id) WHERE g.gene_id IS NULL OR NOT g.tested
    """).fetchone()[0]
    check(n == 0, f"all nominal eQTL genes present and tested in genes ({n} missing)")

    # 3. counts vs paper
    pc, tol = cfg["paper_counts"], cfg["paper_counts"]["tolerance"]
    ne = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'genes.parquet'}' WHERE is_egene").fetchone()[0]
    ns = con.execute(f"SELECT count(*) FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE is_sqtl").fetchone()[0]
    check(abs(ne - pc["egenes"]) / pc["egenes"] <= tol, f"eGenes {ne} within {tol:.0%} of paper {pc['egenes']}")
    check(abs(ns - pc["sqtl_sig_phenotypes"]) / pc["sqtl_sig_phenotypes"] <= tol, f"sQTL sig {ns} within {tol:.0%} of paper {pc['sqtl_sig_phenotypes']}")

    # 4. paper variants resolve
    for key, rs in cfg["paper_variants"].items():
        c, p = key.split(":")
        got = [r[0] for r in con.execute(f"SELECT DISTINCT rsid FROM {variants_sql(cfg, c)} WHERE position = {p}").fetchall()]
        check(rs in got, f"{key} -> {rs} (got {got})")

    # 5. row-group pruning on FLNC, in its bin file
    import pyarrow.parquet as pq
    flnc_bin = con.execute(f"SELECT bin FROM '{cfg.tables / 'genes.parquet'}' WHERE gene_id = 'ENSG00000128591'").fetchone()[0]
    flnc = cfg.tables / "cis_eqtl_nominal" / "chr=chr7" / f"bin={flnc_bin}" / "data.parquet"
    md = pq.read_metadata(flnc)
    hits = 0
    for i in range(md.num_row_groups):
        st = md.row_group(i).column(0).statistics
        if st and st.min <= "ENSG00000128591" <= st.max:
            hits += 1
    check(hits == 1, f"FLNC query touches {hits} row group(s) of {md.num_row_groups} in chr7 bin {flnc_bin} eQTL file ({md.serialized_size / 1e3:.0f} KB footer)")
    # sQTL rows are grouped per intron: one intron of CAMK2D (nine significant) touches one row group
    camk2d = "ENSG00000145349"
    camk2d_bin = con.execute(f"SELECT bin FROM '{cfg.tables / 'genes.parquet'}' WHERE gene_id = '{camk2d}'").fetchone()[0]
    introns = [r[0] for r in con.execute(f"SELECT phenotype_id FROM '{cfg.tables / 'splice_phenotypes.parquet'}' WHERE gene_id = '{camk2d}' AND is_sqtl").fetchall()]
    md = pq.read_metadata(cfg.tables / "cis_sqtl_nominal" / "chr=chr4" / f"bin={camk2d_bin}" / "data.parquet")
    pcol = [md.row_group(0).column(j).path_in_schema for j in range(md.row_group(0).num_columns)].index("phenotype_id")
    touched = {p: 0 for p in introns}
    for i in range(md.num_row_groups):
        st = md.row_group(i).column(pcol).statistics
        for p in introns:
            if st and st.min <= p <= st.max:
                touched[p] += 1
    check(all(n == 1 for n in touched.values()), f"each of CAMK2D's {len(introns)} significant introns touches one row group of {md.num_row_groups} in chr4 bin {camk2d_bin} sQTL file ({md.serialized_size / 1e3:.0f} KB footer): {sorted(set(touched.values()))}")

    # 6. rsID exact rate among cis variants (trans eQTL-only positions have no alleles and can only match by position)
    vpos = variants_sql(cfg)
    tot, ex = con.execute(f"SELECT count(*), sum(match = 'exact') FROM {vpos} WHERE in_cis").fetchone()
    check(ex / tot >= 0.90, f"rsID exact match rate among cis variants {ex / tot:.1%}")

    # 7. the cis side of the variant table is unchanged by adding trans positions: same row count as
    #    the cis-only build (2026-09-03), and no position carries both allele-bearing and allele-less rows
    check(tot == 8_872_723, f"cis variant rows {tot:,} (expected 8,872,723)")
    mixed = con.execute(f"SELECT count(*) FROM (SELECT chr, position FROM {vpos} GROUP BY 1, 2 HAVING bool_or(A1 IS NULL) AND bool_or(A1 IS NOT NULL))").fetchone()[0]
    check(mixed == 0, f"positions with both allele-bearing and allele-less rows: {mixed}")

    # 8. trans rows resolve to an rsID except where dbSNP has no record at the position
    tn, tnull = con.execute(f"SELECT count(*), sum(rsid IS NULL) FROM read_parquet('{cfg.tables}/trans/*/*.parquet', hive_partitioning=false)").fetchone()
    check(tnull / tn < 0.001, f"trans rows without rsID {tnull:,} of {tn:,} ({tnull / tn:.3%})")

    # 8.5 rsID text matches its own number, everywhere a text rsID is stored. Two separate min()s
    #     in `variants_rsid` used to take the text of one dbSNP record and the number of another
    #     (14,803 variants, 39 gene leads, 86 intron leads, 1,376 credible-set rows).
    G, S, C = (f"'{cfg.tables / n}.parquet'" for n in ("genes", "splice_phenotypes", "credible_sets"))
    text_bad = []
    for label, sql in (
        ("the variant table", f"SELECT count(*) FROM {vpos} WHERE rsid IS NOT NULL AND rsid <> 'rs' || rs_number"),
        ("genes.lead_rsid", f"""SELECT count(*) FROM {G} g JOIN {vpos} v ON v.chr = g.chr AND v.position = g.lead_position
             AND v.A1 = g.lead_A1 AND v.A2 = g.lead_A2 WHERE g.lead_rsid IS NOT NULL AND g.lead_rsid <> 'rs' || v.rs_number"""),
        ("intron lead_rsid", f"""SELECT count(*) FROM {S} s JOIN {vpos} v ON v.chr = s.chr AND v.position = s.lead_position
             AND v.A1 = s.lead_A1 AND v.A2 = s.lead_A2 WHERE s.lead_rsid IS NOT NULL AND s.lead_rsid <> 'rs' || v.rs_number"""),
        ("credible_sets.rsid", f"""SELECT count(*) FROM {C} c JOIN {vpos} v ON v.chr = c.chr AND v.position = c.position
             AND v.A1 = c.A1 AND v.A2 = c.A2 WHERE c.rsid IS NOT NULL AND c.rsid <> 'rs' || v.rs_number"""),
    ):
        n = con.execute(sql).fetchone()[0]
        if n:
            text_bad.append(f"{label}: {n:,}")
    check(not text_bad, "every stored rsID text equals 'rs' || rs_number of its variant" + (f" ({', '.join(text_bad)})" if text_bad else ""))

    # 8.6 nothing the browser could read is parquet any more: every table is under `_tables/`, and
    #     `_*` folders are never uploaded
    strays = sorted(str(f.relative_to(cfg.derived)) for f in cfg.derived.rglob("*.parquet")
                    if not any(part.startswith("_") for part in f.relative_to(cfg.derived).parts))
    check(not strays, f"no parquet outside `_*` folders in data/derived ({len(strays)} found: {strays[:3]})")

    # 8.7 content addressing (SPEC section 3): every browser-read file is named by its own bytes,
    #     the manifest agrees with what is on disk, and no unhashed copy is left behind
    _validate_immutable(cfg, check)

    # 9. binary packs on every chromosome (SPEC.md): structure, round trip, reference files
    steps_pack.validate(cfg, con, check)
    # 10. trans-only variant sections and trans packs (SPEC.md sections 4 and 12)
    steps_pack_trans.validate(cfg, con, check)
    # 11. hits packs, rsID index and variant index (SPEC.md sections 13 to 15)
    steps_pack_variant.validate(cfg, con, check)

    if fails:
        raise SystemExit(f"validate: {len(fails)} check(s) failed")
    log("validate: all checks passed")
