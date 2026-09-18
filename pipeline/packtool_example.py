"""A worked example of packtool: rows to packs and back, on synthetic data, with every command shown.

    uv run python -m pipeline.packtool_example [DIR]

Writes three small TSV tables into DIR (default: a fresh temporary directory), packs them, inspects
the packs, and decodes rows back, printing each command line before the first lines of its output.
Needs no data download and finishes in a few seconds. Read PACKS.md alongside it.

The tables are the shapes the pack builders start from:

- variants.tsv     one chromosome's cis variants: position, A1, A2, rs_number, af, ma_samples, ma_count, match
- trans_only.tsv   the variants seen only in trans: position, A1, A2 (empty when not reported), rs_number, af, match
- eqtl_rows.tsv    nominal rows: phenotype_id, position, A1, A2, tss_distance, pval_nominal, slope, slope_se, pip, cs_id
- trans_rows.tsv   trans rows: gene_id, qtl_type, phenotype_id, variant_chr, position, rs_number, af, pval, beta
- gwas.tsv         GWAS rows: chr, position, ea, nea, rs_number, beta, se, eaf, p, n
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
from scipy.special import stdtr

from . import packtool as pt

DOF = 435


def make_tables(d: Path) -> None:
    rng = np.random.default_rng(42)
    n = 2000
    snps = list(pt.packfmt.SNP_CODES)
    pos = np.cumsum(rng.integers(1, 60, n)) + 5_000_000
    pairs = [snps[i] if rng.random() < 0.95 else ("AT", "A") for i in rng.integers(0, 12, n)]
    variants = pa.table({
        "position": pos, "A1": [p[0] for p in pairs], "A2": [p[1] for p in pairs],
        "rs_number": pa.array(rng.integers(1, 10**9, n), mask=rng.random(n) < 0.05),
        "af": rng.uniform(0.01, 0.99, n), "ma_samples": rng.integers(5, 500, n), "ma_count": rng.integers(5, 520, n),
        "match": ["exact"] * n,
    })
    pt.write_table(variants, d / "variants.tsv")
    m = 120
    pt.write_table(pa.table({
        "position": np.arange(m) * 500 + 5_200_000, "A1": [None if i % 4 == 0 else "A" for i in range(m)],
        "A2": [None if i % 4 == 0 else "G" for i in range(m)], "rs_number": rng.integers(1, 10**9, m),
        "af": rng.uniform(0.01, 0.99, m), "match": ["position" if i % 4 == 0 else "exact" for i in range(m)],
    }), d / "trans_only.tsv")

    rows = {k: [] for k in ("phenotype_id", "position", "A1", "A2", "tss_distance", "pval_nominal", "slope", "slope_se", "pip", "cs_id")}
    for gid, first, count in (("ENSG00000000001", 0, 300), ("ENSG00000000002", 200, 900)):
        se = rng.uniform(0.05, 0.3, count)
        slope = rng.normal(0, 0.2, count)
        p = 2.0 * stdtr(DOF, -np.abs(slope) / se)
        tss = int(pos[first]) + 400_000                         # the window is [tss - 1 Mb, tss + 1 Mb] in the real data
        best = np.argsort(p)[:2]
        for i in range(count):
            rows["phenotype_id"].append(gid); rows["position"].append(int(pos[first + i]))
            rows["A1"].append(pairs[first + i][0]); rows["A2"].append(pairs[first + i][1])
            rows["tss_distance"].append(int(pos[first + i]) - tss)
            rows["pval_nominal"].append(float(p[i])); rows["slope"].append(float(slope[i])); rows["slope_se"].append(float(se[i]))
            rows["pip"].append(float(rng.uniform(0.5, 1)) if i in best else None); rows["cs_id"].append(1 if i in best else None)
    pt.write_table(pa.table(rows), d / "eqtl_rows.tsv")
    details = {"ENSG00000000001": {"gene": {"gene_id": "ENSG00000000001", "symbol": "EXAMPLE1", "tss": int(pos[0]) + 400_000},
                                   "exons": [[int(pos[0]) + 400_000, int(pos[0]) + 401_000]], "splice": []}}
    (d / "details.json").write_text(json.dumps(details, indent=1))

    trans = {k: [] for k in ("gene_id", "qtl_type", "phenotype_id", "variant_chr", "position", "rs_number", "af", "pval", "beta")}
    for gid, runs in (("ENSG00000000001", [("e", "ENSG00000000001", "chr3", 8), ("s", "chr21:5400100:5400900:clu_12_+:ENSG00000000001.4", "chrX", 5)]),
                      ("ENSG00000000002", [("e", "ENSG00000000002", "chr7", 6)])):
        for qt, pid, vchr, count in runs:
            for p in np.sort(rng.choice(np.arange(10_000_000, 10_400_000), count, replace=False)).tolist():
                trans["gene_id"].append(gid); trans["qtl_type"].append(qt); trans["phenotype_id"].append(pid); trans["variant_chr"].append(vchr)
                trans["position"].append(p); trans["rs_number"].append(int(rng.integers(1, 10**9))); trans["af"].append(float(rng.uniform(0.05, 0.95)))
                trans["pval"].append(float(10.0 ** -rng.uniform(5, 30))); trans["beta"].append(float(rng.normal(0, 0.3)))
    pt.write_table(pa.table(trans), d / "trans_rows.tsv")

    m = 5000
    gpos = np.cumsum(rng.integers(1, 40, m)) + 5_000_000
    gwas = pa.table({
        "chr": ["chr21"] * m, "position": gpos, "ea": ["A"] * m, "nea": ["G"] * m,
        "rs_number": pa.array(rng.integers(1, 10**9, m), mask=rng.random(m) < 0.1),
        "beta": np.round(rng.normal(0, 0.1, m), 4), "se": np.round(rng.uniform(0.02, 0.5, m), 4), "eaf": np.round(rng.uniform(0, 1, m), 4),
        "p": np.array([float(f"{x:.3e}") for x in rng.uniform(1e-8, 1, m)]), "n": rng.choice([422920, 937963], m),
    })
    pt.write_table(gwas, d / "gwas.tsv")


def run(*args, show: int = 6) -> None:
    """Run one packtool command, printing the command line first and at most `show` lines of its output."""
    args = [str(a) for a in args]
    print(f"\n$ uv run python -m pipeline.packtool {' '.join(args)}")
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = pt.main(args)
    lines = buf.getvalue().rstrip("\n").splitlines()
    for line in lines[:show]:
        print("  " + line)
    if len(lines) > show:
        print(f"  ... ({len(lines) - show} more lines)")
    if "-o" in args:
        out = Path(args[args.index("-o") + 1])
        if out.is_file() and out.suffix in pt.TABLE_EXTS:
            text = out.read_text().rstrip("\n").splitlines() if out.suffix in (".tsv", ".csv", ".json") else pt._tsv_bytes(pt.read_table(out)).decode().splitlines()
            print(f"  -> {out.name}: {len(text) - (0 if out.suffix == '.json' else 1)} rows")
            for line in text[:show]:
                print("     " + line)
    assert rc == 0, f"command failed with exit code {rc}"


def main() -> int:
    d = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="packtool_example_"))
    d.mkdir(parents=True, exist_ok=True)
    make_tables(d)
    print(f"tables written to {d}: variants.tsv, trans_only.tsv, eqtl_rows.tsv, details.json, trans_rows.tsv, gwas.tsv")
    (d / "manifest.json").write_text(json.dumps({"packs": {"dof": {"eqtl": DOF, "sqtl": 480}}}))
    print("manifest.json holds the degrees of freedom the slope, beta_se, and r2 derivations need (packs.dof)")

    print("\n# 1. rows -> packs")
    run("pack-variants", d / "variants.tsv", "-o", d / "chr21.qbv", "--chrom", "chr21", "--trans-only", d / "trans_only.tsv", show=14)
    run("pack-results", d / "eqtl_rows.tsv", "--variants", d / "chr21.qbv", "-o", d / "chr21.qbe", "--details", d / "details.json",
        "--pointers", d / "pointers.tsv", show=14)
    run("pack-trans", d / "trans_rows.tsv", "-o", d / "chr21.qbt", "--chrom", "chr21", "--pointers", d / "trans_pointers.tsv", show=14)
    run("pack-gwas", d / "gwas.tsv", "-o", d / "gwas", show=8)

    print("\n# 2. inspect")
    run("header", d / "chr21.qbv", show=13)
    run("variants", d / "chr21.qbv", "--pages")
    run("blocks", d / "chr21.qbe", "--gene-ids")
    run("block", d / "chr21.qbe", "--gene", "EXAMPLE1", "--header", show=16)
    run("frames", d / "chr21.qbt")
    run("gwas-index", d / "gwas" / "gwas_index.bin")

    print("\n# 3. packs -> rows")
    run("block", d / "chr21.qbe", "--gene", "ENSG00000000002", "--variants", d / "chr21.qbv", "-o", d / "gene2_rows.tsv")
    run("block", d / "chr21.qbe", "--gene", "ENSG00000000001", "--cs")
    run("block", d / "chr21.qbe", "--gene", "ENSG00000000001", "--details")
    run("variants", d / "chr21.qbv", "--vidx", 510, "--n", 4)
    run("variants", d / "chr21.qbv", "--section", "trans")
    ptr = pt.read_table(d / "trans_pointers.tsv")            # in a pipeline build these offsets sit in search_index (trans_off, trans_len)
    run("trans", d / "chr21.qbt", "--off", ptr["trans_off"][0].as_py(), "--len", ptr["trans_len"][0].as_py(),
        "--gene-id", "ENSG00000000001", "--gene-version", 4, "-o", d / "gene1_trans.tsv")
    run("gwas", d / "gwas" / "chr21.qbg", "--index", d / "gwas" / "gwas_index.bin", "--lo", 5_050_000, "--hi", 5_050_400)

    print("\n# 4. check every byte of a file")
    run("check", d / "chr21.qbv", show=15)
    run("check", d / "chr21.qbe", "--variants", d / "chr21.qbv", show=14)
    run("check", d / "chr21.qbt", show=14)
    run("check", d / "gwas" / "chr21.qbg", "--index", d / "gwas" / "gwas_index.bin", show=12)
    print(f"\ndone; everything is under {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
