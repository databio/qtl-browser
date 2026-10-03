"""Build the shared input tables under `_tables/`: the extracted Zenodo files, the variant table
with rsIDs, the refget reference and allele check, and the v0 tables the TOPCHeF adapter and its
acceptance gate read. The qtlb v1 store itself is built by the adapters and store modules
(pipeline/README.md), not here.

    uv run python -m pipeline build                    # all steps, skipping ones already done
    uv run python -m pipeline build --step nominal --force
    uv run python -m pipeline steps                    # list steps in order
"""
import argparse
import inspect
import sys
import time

from . import steps_extract, steps_gtf, steps_nominal, steps_refget, steps_tables, steps_variants
from .common import Config, log

STEPS = [
    ("extract", steps_extract.run),
    ("gtf", steps_gtf.run),
    ("variants_collect", steps_variants.collect),
    ("variants_rsid", steps_variants.rsid),
    ("refget_store", steps_refget.refget_store),
    ("variants_refcheck", steps_refget.variants_refcheck),
    ("permutation_tables", steps_tables.permutation_tables),
    ("credible_sets", steps_tables.credible_sets),
    ("nominal", steps_nominal.run),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--step", action="append", help="run only this step (repeatable)")
    b.add_argument("--force", action="store_true", help="re-run even if marked done")
    sub.add_parser("steps")
    args = ap.parse_args()
    cfg = Config()

    if args.cmd == "steps":
        print("\n".join(n for n, _ in STEPS))
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
