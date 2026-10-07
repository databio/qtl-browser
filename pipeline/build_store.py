"""Build a whole qtlb store in one command: annotation, variant catalog, experiment, `store.json`.

    uv run python -m pipeline.build_store --store data/derived/store-genome-v2 --version 2
    QTLB_CHROMS=chr21,chr22 uv run python -m pipeline.build_store --store data/derived/store-local-v2 --version 2
    uv run python -m pipeline.build_store --store data/derived/store-genome-v1 --version 1 --validate

Until this existed the three builders were run as three CLIs by hand, and nothing called
`Store.write_store` -- so `data/derived/store-genome-v1` was built without a `store.json` at all and
could not be validated, which is precisely the check that would have caught a v1 object sitting in
a v2 store. One command, one `--version`, and the version reaches every object's header.

`--version` selects the codec (`qtlstore.codec_for`) and is written to `store.json`. `validate`
compares the two, so a store half-built with one codec fails rather than passing quietly.

Chromosomes come from `pipeline.common.CHROMS`, which `QTLB_CHROMS` narrows -- a subset build gets a
subset variant catalog and subset results, and the two agree because they read the same list.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

from . import annotation, catalog, qtlstore as qs, results
from .common import CHROMS, Config, log

GTF = Path("data/raw/gencode_v34/gencode.v34.annotation.gtf.gz")
SITES = Path("data/derived/_tables/topchef/sites.parquet")
TABLES = Path("data/derived/_tables/topchef")


def build(store_root: Path, version: int, *, name: str = "topchef", annot_id: str = "gencode_v34",
          cat_id: str = "topchef_grch38", exp_id: str = "topchef", gtf: Path = GTF, sites: Path = SITES,
          tables: Path = TABLES, chroms: list[str] | None = None, fresh: bool = True) -> dict:
    """Annotation, variant catalog, experiment, then `store.json`. Returns a summary.

    Order matters and is not cosmetic: `Store.write_pointer` refuses a pointer naming an object that
    is not in the store yet, the experiment names the catalog and annotation by id, and `store.json`
    lists whatever pointers exist when it is written -- so it goes last."""
    chroms = list(chroms if chroms is not None else CHROMS)
    codec = qs.codec_for(version)
    if fresh and store_root.exists():
        shutil.rmtree(store_root)
    store_root.mkdir(parents=True, exist_ok=True)
    st = qs.Store(store_root)
    cfg = Config()
    out: dict = {"store": str(store_root), "format_version": version, "chromosomes": len(chroms)}

    t0 = time.time()
    log(f"== annotation {annot_id}")
    a = annotation.build(st, annot_id, gtf, annotation.gtf_source(Path("data/raw/sources.yaml"), gtf),
                         version=version)
    out["annotation"] = {"id": annot_id, "genes": a["n_genes"], "min": round((time.time() - t0) / 60, 1)}

    t0 = time.time()
    log(f"== variant catalog {cat_id}")
    from .steps_refget import open_store
    c = catalog.build(st, cat_id, sites, open_store(cfg), cfg["reference"]["collection"], chroms,
                      source={"sites": str(sites)}, version=version)
    out["catalog"] = {"id": cat_id, "n_sites": c["n_sites"], "identity": c["identity_digest"],
                      "min": round((time.time() - t0) / 60, 1)}

    t0 = time.time()
    log(f"== experiment {exp_id} (codec v{codec.FORMAT_VERSION})")
    e = results.build(st, exp_id, tables, cat_id, annot_id, chroms, codec=codec)
    out["experiment"] = {
        "id": exp_id, "min": round((time.time() - t0) / 60, 1),
        "n_blocks": {k: v.get("n_blocks") for k, v in _result_types(e).items()},
        "precision": {k: v.get("precision") for k, v in _result_types(e).items()},
        "trans": {k: (v.get("trans") or {}).get("precision") for k, v in _result_types(e).items()},
        "gwas": {k: (e.get("gwas") or {}).get(k) for k in ("n_rows", "block_rows", "n_values")} if e.get("gwas") else None,
    }

    st.write_store(name, _refget_urls(cfg), version)
    out["bytes"] = sum(p.stat().st_size for p in (store_root / "immutable").glob("*") if p.is_file())
    return out


def _result_types(doc: dict) -> dict:
    """`{phenotype_type: entry}` across the shapes an experiment pointer uses, so the summary carries
    the precision block whichever one this experiment has."""
    r = doc.get("results")
    if isinstance(r, list):
        return {x.get("phenotype_type", str(i)): x for i, x in enumerate(r)}
    if isinstance(r, dict):
        return r
    return {}


def _refget_urls(cfg) -> list[str]:
    url = (cfg["reference"] or {}).get("store_url")
    return [url] if url else []


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--store", required=True, type=Path)
    ap.add_argument("--version", type=int, default=qs.FORMAT_VERSION, choices=qs.SUPPORTED_VERSIONS,
                    help="format version: selects the codec and is written to store.json")
    ap.add_argument("--gtf", type=Path, default=GTF)
    ap.add_argument("--sites", type=Path, default=SITES)
    ap.add_argument("--tables", type=Path, default=TABLES)
    ap.add_argument("--keep", action="store_true", help="add to an existing store instead of rebuilding it")
    ap.add_argument("--validate", action="store_true",
                    help="run Store.validate with the refgetstore and a site reader afterwards")
    args = ap.parse_args(argv)
    out = build(args.store, args.version, gtf=args.gtf, sites=args.sites, tables=args.tables,
                fresh=not args.keep)
    if args.validate:
        st = qs.Store(args.store)
        fails = st.validate(refget=Config(), sites=catalog.read_sites)
        out["validate"] = {"fails": fails, "notes": st.notes}
    print(json.dumps(out, indent=1))
    return 1 if args.validate and out["validate"]["fails"] else 0


if __name__ == "__main__":
    sys.exit(main())
