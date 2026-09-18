"""Build the browser-ready parquet tables.

    uv run python -m pipeline build                    # all steps, skipping ones already done
    uv run python -m pipeline build --step nominal --force
    uv run python -m pipeline validate
    uv run python -m pipeline steps                    # list steps in order
"""
import argparse
import inspect
import sys
import time

from . import (steps_extract, steps_finish, steps_gtf, steps_gwas, steps_nominal, steps_pack, steps_pack_trans,
               steps_pack_variant, steps_tables,
               steps_variants)
from .common import Config, log

STEPS = [
    ("extract", steps_extract.run),
    ("gtf", steps_gtf.run),
    ("variants_collect", steps_variants.collect),
    ("variants_rsid", steps_variants.rsid),
    ("permutation_tables", steps_tables.permutation_tables),
    ("credible_sets", steps_tables.credible_sets),
    ("trans", steps_tables.trans),
    ("nominal", steps_nominal.run),
    ("pack_sqtl", steps_pack.sqtl),
    ("pack_eqtl", steps_pack.eqtl),
    ("pack_variants_trans", steps_pack_trans.variants_trans),
    ("pack_trans", steps_pack_trans.pack),
    ("search_index", steps_pack.search_index),
    ("pack_hits", steps_pack_variant.pack_hits),
    ("pack_variant_index", steps_pack_variant.pack_variant_index),
    ("coloc_stub", steps_tables.coloc_stub),
    ("gwas_bins", steps_gwas.run),
    ("pack_gwas", steps_gwas.pack),
    ("manifest", steps_finish.manifest),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--step", action="append", help="run only this step (repeatable)")
    b.add_argument("--force", action="store_true", help="re-run even if marked done")
    sub.add_parser("validate")
    sub.add_parser("steps")
    pk = sub.add_parser("packcheck", help="genome-wide checks and measurements for the pack format (SPEC.md); not a build step")
    pk.add_argument("action", choices=["dof", "run", "report", "measure", "roundtrip"])
    pk.add_argument("--type", action="append", choices=["e", "s"], help="QTL type (repeatable; default both)")
    pk.add_argument("--chrom", nargs="+", help="chromosomes (default all)")
    pk.add_argument("--workers", type=int, help="processes (default config workers)")
    pk.add_argument("--source", choices=["auto", "raw", "derived"], default="auto",
                    help="nominal rows from the extracted Zenodo files, the derived tables, or raw when complete (default)")
    pk.add_argument("--md5", action="store_true", help="report: md5 the Zenodo archives (cached by size and mtime)")
    pk.add_argument("--steepest", type=int, help="roundtrip: only the N phenotypes per type with the largest -log10 p genome-wide")
    args = ap.parse_args()
    cfg = Config()

    if args.cmd == "steps":
        print("\n".join(n for n, _ in STEPS))
        return 0
    if args.cmd == "validate":
        steps_finish.validate(cfg)
        return 0
    if args.cmd == "packcheck":
        from . import packcheck
        packcheck.main(cfg, args)
        return 0
    names = [n for n, _ in STEPS]
    for s in args.step or []:
        if s not in names:
            sys.exit(f"unknown step {s}; choose from {names}")
    for name, fn in STEPS:
        if args.step and name not in args.step:
            continue
        if cfg.is_done(name) and not args.force:
            log(f"{name}: done, skipping (use --force)")
            continue
        t0 = time.time()
        log(f"== {name}")
        if "force" in inspect.signature(fn).parameters:
            fn(cfg, force=args.force)
        else:
            fn(cfg)
        cfg.mark_done(name)
        log(f"== {name} finished in {(time.time() - t0) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
