"""The acceptance gate for a TOPCHeF re-ingestion (pipeline/CONTRACT.md, "Acceptance gate").

    uv run python -m pipeline.adapters.verify_topchef [--chrom chr21 chr22] [--genes N] [--skip KIND]

Not a build step and not a unit test: this reads the contract tables the adapter wrote and the v0
tables and packs sitting next to them, and answers the one question the gate asks -- did any number
change. The gate is written against a genome-wide inversion that no precision check can see, so
every comparison here is on exact bits, not tolerances, and the counts are printed whether they
pass or fail.

Six checks:

  orientation  the ingestion report's swapped and dropped counts, which must both be 0 for cis
  sites        every v0 cis variant has one site row; ref is the old A2 and alt the old A1 at every
               SNP; at every indel the pair is one of the two orderings and the tie-break is
               accounted for; `af` matches the v0 variants pack's stored code exactly
  nominal      every contract nominal row against what v0 serves (the eQTL build intermediate; the
               raw Zenodo sQTL file, every tested intron): `beta` bit-identical to `slope`
               **including its sign**, `se` to `slope_se`, `pvalue` to `pval_nominal`, zero unpaired
  packs        the contract rows re-quantized through `packfmt` against the v0 pack's own codes,
               for a sample of genes and introns: the u16 nlp code, the u16 SE code (which is where
               the slope's sign lives), the block scales, and the decoded values. A sampled intron
               with a v0 block and no contract rows fails
  pack_counts  every intron with a v0 sQTL block has contract nominal rows, the same number of them
               (from the block headers, nothing decoded), and `has_nominal = true`
  counts       table row counts against the v0 tables
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from .. import packfmt_v0 as packfmt
from .. import packtool as pt
from ..common import CHROMS, Config, connect, log, pack_file, variants_sql
from ..steps_pack import _raw_sqtl, pointer_dir
from . import topchef as tc

# the genome, not `common.CHROMS`, which QTLB_CHROMS narrows: the v0 tables next to a smoke build
# are still genome-wide, so a count comparison has to know it is looking at a subset
ALL_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

RESULTS: list[tuple[bool, str]] = []


def check(ok: bool, msg: str) -> bool:
    RESULTS.append((bool(ok), msg))
    print(f"{'PASS' if ok else 'FAIL'}  {msg}", flush=True)
    return bool(ok)


def _bits(col) -> np.ndarray:
    """A float column as its raw bit pattern, so two arrays compare exactly. NaN payloads and the
    sign of zero both survive, which value equality would hide."""
    a = np.asarray(col.to_numpy(zero_copy_only=False), dtype=np.float64 if col.type == pa.float64() else np.float32)
    return a.view(np.uint64 if a.dtype == np.float64 else np.uint32)


def _diff(a: pa.ChunkedArray, b: pa.ChunkedArray) -> tuple[int, int]:
    """(rows whose bits differ, rows whose null-ness differs)."""
    na, nb = pc.is_null(a).to_numpy(zero_copy_only=False), pc.is_null(b).to_numpy(zero_copy_only=False)
    null_diff = int(np.sum(na != nb))
    both = ~na & ~nb
    return int(np.sum(_bits(a)[both] != _bits(b)[both])), null_diff


# ---- 1. orientation ----------------------------------------------------------------------------
def verify_orientation(cfg: Config) -> None:
    path = tc.report_path(cfg)
    if not path.exists():
        check(False, f"ingestion report {path.name} exists (run the adapter first)")
        return
    r = json.loads(path.read_text())
    cis, tr = r["counts"]["cis"], r["counts"]["trans_only"]
    check(cis["swapped"] == 0,
          f"orientation: {cis['swapped']:,} cis sites swapped (the gate requires 0); {cis['as_is']:,} as-is")
    check(cis["dropped"] == 0,
          f"orientation: {cis['dropped']:,} cis sites dropped (the gate requires 0)")
    tie = r.get("indel_tie_break", {})
    check(True, f"orientation: indel tie-breaks: {tie.get('disputed', 0):,} disputed, "
                f"{tie.get('a2_reads', 0):,} where the reference also reads A2 (so A2 is taken), "
                f"{tie.get('kept_refcheck_call', 0):,} where it does not (the refcheck call stands)")
    check(True, f"orientation: trans-only {tr['as_is']:,} as-is, {tr['swapped']:,} swapped, "
                f"{tr['dropped']:,} dropped {r.get('dropped_by_refcheck_class', {})}")
    ex = r.get("trans_eqtl_excluded", {})
    by_class = r.get("dropped_by_refcheck_class", {})
    check(ex.get("variants") == by_class.get("no_source_alleles", 0) == tr["dropped"],
          f"trans: {ex.get('variants', 0):,} allele-less trans eQTL variants excluded, the only trans-only drops "
          f"({tr['dropped']:,}); {ex.get('trans_eqtl_rows', 0):,} of {ex.get('trans_eqtl_rows_total', 0):,} "
          f"trans eQTL rows left out, no alleles inferred")


# ---- 2. sites ----------------------------------------------------------------------------------
def verify_sites(cfg: Config, chroms: list[str]) -> None:
    con = connect(cfg, memory_limit=cfg["duckdb_memory_limit"], threads=cfg["duckdb_threads"])
    where = "chr IN (" + ", ".join(f"'{c}'" for c in chroms) + ")"
    # `lo`/`hi` are the allele pair in a fixed order, so "same site, either orientation" is an
    # equi-join. An OR of the two orderings in the ON clause makes DuckDB fall back to a nested
    # loop over 9M x 9M rows, which is what stalled the first genome-wide gate for two hours.
    con.execute(f"""CREATE VIEW s AS SELECT *, least(ref, alt) AS lo, greatest(ref, alt) AS hi
                    FROM '{tc.sites_path(cfg)}' WHERE {where}""")
    con.execute(f"""CREATE VIEW v AS SELECT *, least(A1, A2) AS lo, greatest(A1, A2) AS hi
                    FROM {variants_sql(cfg)} WHERE {where} AND A1 IS NOT NULL""")
    con.execute(f"CREATE VIEW rc AS SELECT * FROM read_parquet('{cfg.tables / 'refcheck' / '*.parquet'}') WHERE {where}")

    # v0 variants with no alleles (trans eQTL-only positions) are excluded from `sites` by rule and
    # are not in `v`; every site must be an allele-bearing v0 variant and vice versa
    missing, extra = con.execute("""
        SELECT (SELECT count(*) FROM v ANTI JOIN s ON s.chr = v.chr AND s.pos = v.position AND s.lo = v.lo AND s.hi = v.hi),
               (SELECT count(*) FROM s ANTI JOIN v ON s.chr = v.chr AND s.pos = v.position AND s.lo = v.lo AND s.hi = v.hi)
        """).fetchone()
    check(missing == 0 and extra == 0,
          f"sites: {missing:,} v0 allele-bearing variants with no site row, {extra:,} site rows with no v0 variant")

    rows = con.execute("""
        SELECT CASE WHEN length(v.A1) = 1 AND length(v.A2) = 1 THEN 'snp' ELSE 'indel' END AS shape,
               rc.match::VARCHAR AS match, v.in_cis,
               count(*) FILTER (WHERE s.ref = v.A2 AND s.alt = v.A1) AS ref_is_a2,
               count(*) FILTER (WHERE s.ref = v.A1 AND s.alt = v.A2 AND v.A1 <> v.A2) AS ref_is_a1,
               count(*) AS n
        FROM v JOIN rc ON rc.chr = v.chr AND rc.position = v.position AND rc.A1 = v.A1 AND rc.A2 = v.A2
               JOIN s ON s.chr = v.chr AND s.pos = v.position AND s.lo = v.lo AND s.hi = v.hi
        GROUP BY 1, 2, 3 ORDER BY 3 DESC, 1, 2""").fetchall()
    snp_a1 = sum(r[4] for r in rows if r[0] == "snp" and r[2])
    check(snp_a1 == 0, f"sites: ref is the old A2 and alt the old A1 at every cis SNP "
                       f"({snp_a1:,} SNPs took A1 as ref)")
    for shape, match, in_cis, a2, a1, n in rows:
        check(True, f"sites: {'cis ' if in_cis else 'trans'} {shape:5s} refcheck={match:9s} {n:>12,}: "
                    f"ref=A2 {a2:>12,}, ref=A1 {a1:>10,}")

    # `af` against the v0 variants pack, which is where the browser reads it from
    bad = total = 0
    for chrom in chroms:
        path = pack_file(cfg, f"variants/{chrom}", "qbv")
        if not path.exists():
            check(False, f"sites: v0 variants pack for {chrom} is missing ({path.name})")
            continue
        dec = pt.variant_rows(path, section="cis")
        src = con.execute("""SELECT pos, ref, alt, af FROM s WHERE chr = ? AND in_cis
                             ORDER BY pos, alt, ref""", [chrom]).fetch_arrow_table()
        if dec.num_rows != src.num_rows:
            check(False, f"sites {chrom}: the v0 variants pack holds {dec.num_rows:,} cis records, "
                         f"sites holds {src.num_rows:,}")
            continue
        # the pack stores rint(af * 65534); compare on that code, which is the stored value
        want = np.rint(src["af"].to_numpy(zero_copy_only=False).astype(np.float64) * 65534)
        got = np.rint(dec["af"].to_numpy(zero_copy_only=False).astype(np.float64) * 65534)
        nn = ~(np.isnan(want) | np.isnan(got))
        bad += int(np.sum(want[nn] != got[nn])) + int(np.sum(np.isnan(want) != np.isnan(got)))
        total += int(nn.sum())
    check(bad == 0, f"sites: `af` equals the v0 variants pack's stored code at {total:,} cis variants "
                    f"({bad:,} differ; no site takes 1 - af)")


# ---- 3. nominal --------------------------------------------------------------------------------
V0_PHENOTYPE_TYPES = ("ge", "leafcutter")


def _v0_nominal_sql(cfg: Config, chrom: str, ptype: str) -> str | None:
    """What v0 serves for one chromosome and phenotype type, as rows keyed the way the source keys
    them: (phenotype_id, position, A1, A2) with `slope`, `slope_se`, `pval_nominal`.

    eQTL: the v0 build intermediate, which is what the eQTL packs were encoded from. sQTL: the raw
    Zenodo file, which is what the sQTL packs stream (`steps_pack._raw_sqtl`), every tested intron.
    Not `_tables/cis_sqtl_nominal`: with `sqtl_nominal: significant` that intermediate holds only
    the significant introns, and comparing against it passed a contract table that had lost the
    rest. Slope and SE are cast to the contract's float32, the same cast the adapter makes, so a
    bit comparison is exact; a sign flip still shows as a changed bit."""
    if ptype == "ge":
        files = sorted((cfg.tables / "cis_eqtl_nominal" / f"chr={chrom}").glob("bin=*/data.parquet"))
        if not files:
            return None
        src = "read_parquet([" + ", ".join(f"'{f}'" for f in files) + "], hive_partitioning = false)"
        pid = "gene_id"
    else:
        raw = _raw_sqtl(cfg, chrom)
        if not raw.exists():
            return None
        src, pid = f"'{raw}'", "phenotype_id"
    return f"""SELECT {pid} AS phenotype_id, position, A1, A2, slope::FLOAT AS slope,
                      slope_se::FLOAT AS slope_se, pval_nominal, slope AS slope_raw FROM {src}"""


def _gate_con(cfg: Config):
    """A DuckDB for the big joins: every CPU the job has, and most of its memory when Slurm says how
    much that is. The config's per-worker limits are sized for three build workers sharing a node,
    and a 47-million-row join held to them spends its time spilling."""
    threads = max(int(cfg["duckdb_threads"]), len(os.sched_getaffinity(0)))
    mem = os.environ.get("SLURM_MEM_PER_NODE")
    limit = f"{int(int(mem) * 0.4)}MB" if mem and mem.isdigit() else cfg["duckdb_memory_limit"]
    return connect(cfg, memory_limit=limit, threads=threads)


def verify_nominal(cfg: Config, chroms: list[str]) -> None:
    totals = {"rows": 0, "beta": 0, "se": 0, "pvalue": 0, "sign": 0, "nulls": 0, "unjoined": 0}
    for chrom in chroms:
        cn = tc.nominal_path(cfg, chrom)
        if not cn.exists():
            check(False, f"nominal {chrom}: no contract table ({cn})")
            continue
        # a fresh DuckDB per chromosome: one kept across all of them grew until an allocation failed
        con = _gate_con(cfg)
        con.execute(f"CREATE OR REPLACE TABLE o AS SELECT position, A1, A2, ref, alt FROM {tc.orientation_sql(cfg, chrom)} WHERE in_cis")
        for ptype in V0_PHENOTYPE_TYPES:
            v0 = _v0_nominal_sql(cfg, chrom, ptype)
            if v0 is None:
                check(False, f"nominal {chrom} {ptype}: no v0 rows to compare against")
                continue
            # The contract rows take back the source's A1/A2 from the orientation table, and then
            # meet the v0 rows on the source's own key. A full outer join, so a row on either side
            # with no partner is counted, not skipped; a contract row whose site has no orientation
            # row keeps a null A1 and so cannot pair. `in_v0`/`in_n` rather than a key column: a
            # join key is coalesced and never null. Equi-joins only.
            reader = con.execute(f"""
                SELECT n.beta, n.se, n.pvalue, v.slope, v.slope_se, v.pval_nominal, v.slope_raw,
                       v.in_v0 IS NULL OR n.in_n IS NULL AS unjoined
                FROM (SELECT n.phenotype_id, n.pos, o.A1, o.A2, n.beta, n.se, n.pvalue, true AS in_n
                      FROM (SELECT * FROM '{cn}' WHERE phenotype_type = '{ptype}') n
                      LEFT JOIN o ON o.position = n.pos AND o.ref = n.ref AND o.alt = n.alt) n
                FULL OUTER JOIN (SELECT *, true AS in_v0 FROM ({v0})) v
                  ON v.phenotype_id = n.phenotype_id AND v.position = n.pos AND v.A1 = n.A1 AND v.A2 = n.A2
                """).to_arrow_reader(4_000_000)
            got = {k: 0 for k in totals}
            for t in reader:
                got["rows"] += t.num_rows
                got["unjoined"] += int(pc.sum(pc.cast(t["unjoined"], pa.int64())).as_py() or 0)
                b_bad, b_null = _diff(t["beta"], t["slope"])
                s_bad, s_null = _diff(t["se"], t["slope_se"])
                p_bad, p_null = _diff(t["pvalue"], t["pval_nominal"])
                beta = t["beta"].to_numpy(zero_copy_only=False).astype(np.float64)
                slope = t["slope_raw"].to_numpy(zero_copy_only=False).astype(np.float64)
                nn = ~(np.isnan(beta) | np.isnan(slope))
                got["sign"] += int(np.sum(np.sign(beta[nn]) != np.sign(slope[nn])))
                got["beta"] += b_bad
                got["se"] += s_bad
                got["pvalue"] += p_bad
                got["nulls"] += b_null + s_null + p_null
            reader.close()
            for k in totals:
                totals[k] += got[k]
            log(f"nominal {chrom} {ptype}: {got['rows']:,} rows, {got['unjoined']:,} unjoined, "
                f"beta {got['beta']:,} / se {got['se']:,} / pvalue {got['pvalue']:,} bits differ, "
                f"{got['sign']:,} sign differences")
        con.close()
    check(totals["unjoined"] == 0,
          f"nominal: every contract row pairs with a v0 row and back ({totals['unjoined']:,} unpaired "
          f"of {totals['rows']:,}); zero dropped rows")
    check(totals["sign"] == 0, f"nominal: {totals['sign']:,} slopes changed sign (the gate requires 0)")
    check(totals["beta"] == 0, f"nominal: `beta` is bit-identical to the v0 `slope` at {totals['rows']:,} rows "
                               f"({totals['beta']:,} differ)")
    check(totals["se"] == 0, f"nominal: `se` is bit-identical to the v0 `slope_se` ({totals['se']:,} differ)")
    check(totals["pvalue"] == 0, f"nominal: `pvalue` is bit-identical to the v0 `pval_nominal` "
                                 f"({totals['pvalue']:,} differ)")
    check(totals["nulls"] == 0, f"nominal: the null pattern matches the v0 rows ({totals['nulls']:,} differ)")


# ---- 4. packs ----------------------------------------------------------------------------------
def _phenotype_rows(con, cfg: Config, chrom: str, ptype: str, phenotype_id: str) -> pa.Table:
    """One phenotype's contract nominal rows in the v0 pack's own row order, which is the cis
    section's (position, A1, A2). The A1/A2 come back from the orientation table, so the ordering
    does not quietly assume that nothing swapped."""
    return con.execute(f"""
        SELECT n.pos, o.A1, o.A2, n.beta, n.se, n.pvalue
        FROM (SELECT * FROM '{tc.nominal_path(cfg, chrom)}'
              WHERE phenotype_type = '{ptype}' AND phenotype_id = ?) n
        JOIN {tc.orientation_sql(cfg, chrom)} o
          ON o.position = n.pos AND o.ref = n.ref AND o.alt = n.alt
        ORDER BY n.pos, o.A1, o.A2""", [phenotype_id]).fetch_arrow_table()


def verify_packs(cfg: Config, chroms: list[str], n_genes: int) -> None:
    """Re-quantize the contract rows through `packfmt` and compare against the v0 pack's own codes.

    The u16 SE code is where the slope's sign lives (SPEC section 5, bit 15), so an equal SE code
    array is the statement the gate actually wants: no slope flipped, and nothing re-encoded to a
    different-but-valid byte.
    """
    con = connect(cfg, memory_limit=cfg["packs"]["duckdb_memory_limit"], threads=cfg["packs"]["duckdb_threads"])
    si_path = pack_file(cfg, "search_index", "arrow.zst")
    dofs = tc.dof(cfg)
    wanted = ", ".join(repr(g) for g in list(cfg["packs"]["check_genes"]) + list(cfg["packs"]["sqtl_check_genes"]))
    tot = {"blocks": 0, "rows": 0, "nlp_code": 0, "se_code": 0, "scales": 0, "value": 0, "missing": 0,
           "no_rows": 0}
    for chrom in chroms:
        eqtl = pack_file(cfg, f"eqtl/{chrom}", "qbe")
        sqtl = pack_file(cfg, f"sqtl/{chrom}", "qbs")
        if not eqtl.exists():
            check(False, f"packs: no v0 eQTL pack for {chrom}")
            continue
        # the named check genes first, then the widest windows, which are where a quantization
        # scale has the most room to come out different
        genes = con.execute(f"""SELECT gene_id FROM '{cfg.tables / 'genes.parquet'}'
            WHERE chr = ? AND bin IS NOT NULL AND tested
            ORDER BY gene_id IN ({wanted}) DESC, num_var DESC LIMIT {n_genes}""", [chrom]).fetchall()
        for (gene_id,) in genes:
            try:
                off, length = pt.find_gene_block(eqtl, gene_id, si_path)
                blk = pt.read_block(eqtl, off, length, dofs["ge"])
            except Exception as e:                            # noqa: BLE001 - a missing block is a result
                tot["missing"] += 1
                log(f"packs {chrom} {gene_id}: no block ({e})")
                continue
            rows = _phenotype_rows(con, cfg, chrom, "ge", gene_id)
            _compare_block(blk, rows, dofs["ge"], f"{chrom} ge {gene_id}", tot)
            if not sqtl.exists():
                continue
            # every intron of the same gene, reached the way a reader reaches it: through the
            # gene's details, which carry each intron's block pointer
            for sp in (blk["details"] or {}).get("splice", []):
                srows = _phenotype_rows(con, cfg, chrom, "leafcutter", sp["phenotype_id"])
                if srows.num_rows == 0:
                    # v0 serves this intron, so the contract must too: every tested intron has rows
                    tot["no_rows"] += 1
                    log(f"packs {chrom} leafcutter {sp['phenotype_id']}: a v0 block of {sp.get('n_var', '?')} rows, "
                        f"no contract rows")
                    continue
                sblk = pt.read_block(sqtl, sp["blk_off"], sp["blk_len"], dofs["leafcutter"], kind=packfmt.KIND_SQTL)
                _compare_block(sblk, srows, dofs["leafcutter"], f"{chrom} leafcutter {sp['phenotype_id']}", tot)
    check(tot["nlp_code"] == 0, f"packs: the -log10 p code array re-encodes identically over "
                                f"{tot['blocks']:,} blocks / {tot['rows']:,} rows ({tot['nlp_code']:,} blocks differ)")
    check(tot["se_code"] == 0, f"packs: the SE code array, sign bit included, re-encodes identically "
                               f"({tot['se_code']:,} blocks differ)")
    check(tot["scales"] == 0, f"packs: nlp_max, lse_min and lse_max re-encode identically "
                              f"({tot['scales']:,} blocks differ)")
    check(tot["value"] == 0, f"packs: the decoded -log10 p, SE and slope arrays are bit-identical "
                             f"({tot['value']:,} blocks differ)")
    check(tot["missing"] == 0, f"packs: every sampled phenotype has a v0 block ({tot['missing']:,} missing)")
    check(tot["no_rows"] == 0, f"packs: every sampled intron with a v0 block has contract nominal rows "
                               f"({tot['no_rows']:,} have none)")


def verify_pack_counts(cfg: Config, chroms: list[str]) -> None:
    """Every intron the v0 sQTL pack serves has contract nominal rows, and the same number of them.

    Genome-wide over the requested chromosomes, and cheap: no block is decoded. The row count per
    block comes from the pack's own block headers, and the intron each block belongs to from the
    pointer table `pack_sqtl` wrote beside it (`_tmp/pack_pointers/sqtl_<chr>.parquet`); the two
    must agree block for block before either is believed. The contract side is a count over the
    nominal table's phenotype column. Equi-joins on `phenotype_id` only."""
    con = connect(cfg)
    tot = {"introns": 0, "rows": 0, "no_rows": 0, "no_block": 0, "count_differs": 0, "ptr_bad": 0, "no_flag": 0}
    for chrom in chroms:
        sqtl = pack_file(cfg, f"sqtl/{chrom}", "qbs")
        ptr_path = pointer_dir(cfg) / f"sqtl_{chrom}.parquet"
        cn = tc.nominal_path(cfg, chrom)
        if not (sqtl.exists() and ptr_path.exists() and cn.exists()):
            check(False, f"pack counts {chrom}: missing input(s): "
                         f"{[p.name for p in (sqtl, ptr_path, cn) if not p.exists()]}")
            continue
        heads = pt.list_blocks(sqtl).select(["blk_off", "n_rows"])
        ptr = pq.read_table(ptr_path, columns=["phenotype_id", "blk_off", "n_var"])
        con.register("heads", heads)
        con.register("ptr", ptr)
        bad = con.execute("""SELECT (SELECT count(*) FROM heads) <> (SELECT count(*) FROM ptr)
                                    OR (SELECT count(DISTINCT phenotype_id) FROM ptr) <> (SELECT count(*) FROM ptr),
                                    (SELECT count(*) FROM ptr p ANTI JOIN heads h ON h.blk_off = p.blk_off),
                                    (SELECT count(*) FROM ptr p JOIN heads h ON h.blk_off = p.blk_off
                                     WHERE h.n_rows <> p.n_var)""").fetchone()
        if any(bad):
            tot["ptr_bad"] += 1
            log(f"pack counts {chrom}: the pointer table disagrees with the pack's block headers {bad}")
            continue
        r = con.execute(f"""
            WITH n AS (SELECT phenotype_id, count(*) AS n FROM '{cn}'
                       WHERE phenotype_type = 'leafcutter' GROUP BY 1),
                 f AS (SELECT phenotype_id, has_nominal FROM '{tc.phenotypes_path(cfg)}'
                       WHERE phenotype_type = 'leafcutter')
            SELECT count(p.phenotype_id), coalesce(sum(p.n_var), 0),
                   count(*) FILTER (WHERE p.phenotype_id IS NOT NULL AND n.phenotype_id IS NULL),
                   count(*) FILTER (WHERE p.phenotype_id IS NULL),
                   count(*) FILTER (WHERE p.n_var <> n.n),
                   count(*) FILTER (WHERE p.phenotype_id IS NOT NULL AND f.has_nominal IS DISTINCT FROM true)
            FROM ptr p FULL OUTER JOIN n ON n.phenotype_id = p.phenotype_id
            LEFT JOIN f ON f.phenotype_id = p.phenotype_id""").fetchone()
        con.unregister("heads")
        con.unregister("ptr")
        for k, v in zip(("introns", "rows", "no_rows", "no_block", "count_differs", "no_flag"), r):
            tot[k] += int(v)
        log(f"pack counts {chrom}: {r[0]:,} introns / {r[1]:,} rows in the v0 pack; {r[2]:,} without contract rows, "
            f"{r[3]:,} contract introns without a block, {r[4]:,} row counts differ, {r[5]:,} not flagged has_nominal")
    check(tot["ptr_bad"] == 0, f"pack counts: the pointer tables agree with the sQTL pack block headers "
                               f"({tot['ptr_bad']:,} chromosomes disagree)")
    check(tot["no_rows"] == 0, f"pack counts: every one of {tot['introns']:,} introns with a v0 sQTL block has "
                               f"contract nominal rows ({tot['no_rows']:,} have none)")
    check(tot["no_block"] == 0, f"pack counts: every intron with contract nominal rows has a v0 sQTL block "
                                f"({tot['no_block']:,} have none)")
    check(tot["count_differs"] == 0, f"pack counts: each intron's contract row count equals its v0 block's, "
                                     f"{tot['rows']:,} rows in all ({tot['count_differs']:,} introns differ)")
    check(tot["no_flag"] == 0, f"pack counts: every intron with a v0 sQTL block has has_nominal = true "
                               f"({tot['no_flag']:,} do not)")


def _compare_block(blk: dict, rows: pa.Table, dof: int, label: str, tot: dict) -> None:
    if blk["n_rows"] != rows.num_rows:
        tot["blocks"] += 1
        tot["value"] += 1
        log(f"packs {label}: the v0 block holds {blk['n_rows']:,} rows, the contract table {rows.num_rows:,}")
        return
    p = rows["pvalue"].to_numpy(zero_copy_only=False).astype(np.float64)
    se = rows["se"].to_numpy(zero_copy_only=False).astype(np.float64)
    beta = rows["beta"].to_numpy(zero_copy_only=False).astype(np.float64)
    nlp_q, nlp_max = packfmt.quantize_nlp(p)
    se_q, lse_min, lse_max = packfmt.quantize_se(se, beta)
    tot["blocks"] += 1
    tot["rows"] += rows.num_rows
    if not np.array_equal(nlp_q, blk["nlp_code"]):
        tot["nlp_code"] += 1
        log(f"packs {label}: {int(np.sum(nlp_q != blk['nlp_code'])):,} of {rows.num_rows:,} nlp codes differ")
    if not np.array_equal(se_q, blk["se_code"]):
        tot["se_code"] += 1
        d = se_q ^ blk["se_code"]
        log(f"packs {label}: {int(np.sum(se_q != blk['se_code'])):,} SE codes differ, "
            f"{int(np.sum(d & 0x8000 != 0)):,} of them in the sign bit")
    if (nlp_max, lse_min, lse_max) != (blk["nlp_max"], blk["lse_min"], blk["lse_max"]):
        tot["scales"] += 1
        log(f"packs {label}: scales {(nlp_max, lse_min, lse_max)} != {(blk['nlp_max'], blk['lse_min'], blk['lse_max'])}")
    # decode both sides the way a reader does, and compare the three arrays it hands back
    nlp = packfmt.dequantize_nlp(nlp_q, nlp_max)
    se_d, neg = packfmt.dequantize_se(se_q, lse_min, lse_max)
    slope = packfmt.slope_from_se(se_d, neg, packfmt.p_from_nlp(nlp), dof)
    bad = [name for name, a, b in (("-log10 p", nlp, blk["nlp"]), ("SE", se_d, blk["slope_se"]),
                                   ("slope", slope, blk["slope"])) if not _same(a, b)]
    if bad:
        tot["value"] += 1
        v0 = np.asarray(blk["slope"], dtype=np.float64)
        nn = ~(np.isnan(slope) | np.isnan(v0))
        flips = int(np.sum(np.sign(slope[nn]) != np.sign(v0[nn])))
        log(f"packs {label}: decoded {', '.join(bad)} differs; {flips:,} slope signs differ")


def _same(a: np.ndarray, b: np.ndarray) -> bool:
    """Bitwise equality with NaN counted as equal to NaN."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    na, nb = np.isnan(a), np.isnan(b)
    return bool(np.array_equal(na, nb) and np.array_equal(a[~na].view(np.uint64), b[~nb].view(np.uint64)))


# ---- 5. counts ---------------------------------------------------------------------------------
def verify_counts(cfg: Config, chroms: list[str]) -> None:
    """Contract table row counts against the v0 tables they replace. A `QTLB_CHROMS` subset narrows
    both sides; `permuted` and `phenotypes` have no `chr` column of their own, so the contract side
    narrows on `lead_chr`, which is the phenotype's own chromosome for both phenotype types."""
    con = connect(cfg)
    subset = sorted(chroms) != sorted(ALL_CHROMS)
    inlist = "chr IN (" + ", ".join(f"'{c}'" for c in chroms) + ")"
    v0 = lambda extra="": (f" WHERE {inlist}" if subset else "") if not extra else (
        f" WHERE {extra}" + (f" AND {inlist}" if subset else ""))
    lead = lambda pt: (f"SELECT count(*) FROM '{tc.permuted_path(cfg)}' WHERE phenotype_type = '{pt}'"
                       + (f" AND lead_{inlist}" if subset else ""))
    genes, introns = cfg.tables / "genes.parquet", cfg.tables / "splice_phenotypes.parquet"
    pairs = [
        ("permuted ge", lead("ge"), f"SELECT count(*) FROM '{genes}'{v0('tested')}"),
        ("permuted leafcutter", lead("leafcutter"), f"SELECT count(*) FROM '{introns}'{v0()}"),
        ("credible_sets", f"SELECT count(*) FROM '{tc.credible_sets_path(cfg)}'{v0()}",
         f"SELECT count(*) FROM '{cfg.tables / 'credible_sets.parquet'}'{v0()}"),
        ("phenotypes", f"SELECT count(*) FROM '{tc.phenotypes_path(cfg)}'" + (
             f" WHERE phenotype_id IN (SELECT phenotype_id FROM '{tc.permuted_path(cfg)}' WHERE lead_{inlist})"
             if subset else ""),
         f"SELECT (SELECT count(*) FROM '{genes}'{v0('tested')}) + (SELECT count(*) FROM '{introns}'{v0()})"),
    ]
    for label, a, b in pairs:
        try:
            got, want = con.execute(a).fetchone()[0], con.execute(b).fetchone()[0]
        except Exception as e:                                # noqa: BLE001 - report, do not abort
            check(False, f"counts {label}: {e}")
            continue
        check(got == want, f"counts: {label} {got:,} rows against the v0 table's {want:,}")
    # n_variants counts tested rows (p-value present); that must reproduce the source num_var
    cov = json.loads(tc.report_path(cfg).read_text()).get("phenotype_coverage", {})
    for pt, c in sorted(cov.items()):
        nv = c.get("n_variants", {})
        check(nv.get("recount_differs_from_source") == 0,
              f"counts: {pt} n_variants recount equals the source num_var for "
              f"{nv.get('recounted_from_nominal', 0):,} phenotypes ({nv.get('recount_differs_from_source')} differ; "
              f"{nv.get('nominal_rows_without_pvalue', 0):,} untested nominal rows not counted)")


# ---- main --------------------------------------------------------------------------------------
KINDS = ("orientation", "sites", "nominal", "packs", "pack_counts", "counts")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chrom", nargs="+", help="chromosomes (default every one the adapter wrote)")
    ap.add_argument("--genes", type=int, default=25, help="genes per chromosome for the pack check (default 25)")
    ap.add_argument("--skip", action="append", choices=KINDS, default=[], help="skip a check (repeatable)")
    a = ap.parse_args()
    cfg = Config()
    chroms = a.chrom or [c for c in CHROMS if tc.orientation_path(cfg, c).exists()]
    if not chroms:
        sys.exit("verify_topchef: the adapter has not written an orientation table; run it first")
    log(f"verify_topchef: {len(chroms)} chromosome(s): {', '.join(chroms)}")
    for kind, fn in (("orientation", lambda: verify_orientation(cfg)),
                     ("sites", lambda: verify_sites(cfg, chroms)),
                     ("nominal", lambda: verify_nominal(cfg, chroms)),
                     ("packs", lambda: verify_packs(cfg, chroms, a.genes)),
                     ("pack_counts", lambda: verify_pack_counts(cfg, chroms)),
                     ("counts", lambda: verify_counts(cfg, chroms))):
        if kind in a.skip:
            continue
        t0 = time.time()
        log(f"== {kind}")
        fn()
        log(f"== {kind} in {time.time() - t0:.1f} s")
    bad = [m for ok, m in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} passed, {len(bad)} failed")
    for m in bad:
        print(f"  FAIL {m}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
