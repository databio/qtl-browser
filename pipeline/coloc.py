"""Coloc lists (SPEC.md section 20): which phenotypes of a QTL experiment someone found to colocalize
with a GWAS experiment. Optional, one per declared pair, membership only.

The pipeline never runs coloc. A list comes from coloc results as commonly written -- a CSV or TSV
with a phenotype, gene or symbol column and maybe a PP.H4 column -- named by a `colocs:` entry in
config.yaml. PP.H4 and the lead variant are not stored: the browser computes both when it opens a
locus, from the two experiments' own rows.

    uv run python -m pipeline.coloc add  --store DIR [--id ID]     # every configured list whose pair is in the store
    uv run python -m pipeline.coloc show --store DIR --id ID
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import qtlstore as qs
from .common import ROOT, Config

COLOCS = qs.COLOCS
EXT = "json"
# the column names a coloc results file is read by, first match wins (case-insensitive)
PHENOTYPE_COLUMNS = ("phenotype_id", "molecular_trait_id", "trait_id")
GENE_COLUMNS = ("gene_id", "ensembl_gene_id", "ensembl_id")
SYMBOL_COLUMNS = ("symbol", "gene_name", "gene_symbol", "hgnc_symbol")
PP_H4_COLUMNS = ("pp.h4.abf", "pp.h4", "pp_h4", "pph4", "h4")


def _column(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    return next((lower[n] for n in names if n in lower), None)


def read_results(path: Path) -> pd.DataFrame:
    """A coloc results file as strings, tab- or comma-separated by its extension."""
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, comment="#")


def resolve(df: pd.DataFrame, index: pd.DataFrame, genes: pd.DataFrame, min_pp_h4: float | None = None,
            where: str = "coloc") -> list[dict]:
    """The list's rows, `{gene_id, phenotype_id}` sorted and unique, from a results table.

    `index` is the QTL experiment's search index (`phenotype_id`, `gene_id`) and `genes` its
    annotation's genes (`gene_id`, `name`). A phenotype must be one of the experiment's; a symbol must
    name exactly one gene of the annotation; a gene id loses its version. A row naming only a gene is
    a gene-level result (`phenotype_id` null). With a PP.H4 column, rows below `min_pp_h4` are left
    out, and `min_pp_h4` is required; without one, every row is a reported result."""
    pc, gc, sc, hc = (_column(df, n) for n in (PHENOTYPE_COLUMNS, GENE_COLUMNS, SYMBOL_COLUMNS, PP_H4_COLUMNS))
    if not (pc or gc or sc):
        raise ValueError(f"{where}: no phenotype, gene or symbol column among {list(df.columns)}")
    if hc:
        if min_pp_h4 is None:
            raise ValueError(f"{where}: the results have a PP.H4 column ({hc}); set min_pp_h4 for the list")
        df = df[pd.to_numeric(df[hc], errors="raise") >= min_pp_h4]
    pheno_gene = dict(zip(index["phenotype_id"], index["gene_id"]))
    by_symbol: dict[str, list[str]] = {}
    for gid, name in zip(genes["gene_id"], genes["name"]):
        if name:
            by_symbol.setdefault(name, []).append(gid)
    tested = set(g for g in index["gene_id"] if g)
    rows, bad = set(), []
    for i, r in enumerate(df.to_dict("records")):
        ph = (r.get(pc) or "").strip() if pc else ""
        gid = (r.get(gc) or "").strip().split(".")[0] if gc else ""
        sym = (r.get(sc) or "").strip() if sc else ""
        if ph:
            if ph not in pheno_gene:
                bad.append(f"row {i}: phenotype {ph} is not in the experiment")
                continue
            gid = gid or pheno_gene[ph] or ""
        if not gid and sym:
            hits = by_symbol.get(sym, [])
            if len(hits) != 1:
                bad.append(f"row {i}: symbol {sym} names {len(hits)} genes of the annotation")
                continue
            gid = hits[0]
        if not gid:
            bad.append(f"row {i}: no gene")
            continue
        if gid not in tested:
            bad.append(f"row {i}: gene {gid} has no phenotype in the experiment")
            continue
        rows.add((gid, ph or None))
    if bad:
        raise ValueError(f"{where}: {len(bad)} rows do not resolve: " + "; ".join(bad[:5]))
    return [{"gene_id": g, "phenotype_id": p} for g, p in sorted(rows, key=lambda x: (x[0], x[1] or ""))]


def add(store: qs.Store, cid: str, entry: dict, root: Path = ROOT) -> dict:
    """Write one list from its config entry: the rows object, then the pointer."""
    for k in ("qtl", "gwas", "file", "source", "threshold"):
        if not entry.get(k):
            raise ValueError(f"colocs.{cid}: {k} is required")
    qtl, gwas = store.load("experiments", entry["qtl"]), store.load("experiments", entry["gwas"])
    if qtl.get("kind") not in qs.QTL_KINDS or gwas.get("kind") != "gwas":
        raise ValueError(f"colocs.{cid}: {entry['qtl']} is a {qtl.get('kind')}, {entry['gwas']} a {gwas.get('kind')}")
    from . import annotation, results
    index = results.load_index(store, qtl).select(["phenotype_id", "gene_id"]).to_pandas()
    genes = annotation.load(store, qtl["annotation"])["genes"].select(["gene_id", "name"]).to_pandas()
    path = Path(entry["file"])
    rows = resolve(read_results(path if path.is_absolute() else root / path), index, genes,
                   entry.get("min_pp_h4"), f"colocs.{cid}")
    if not rows:
        raise ValueError(f"colocs.{cid}: no rows")
    obj = store.put(json.dumps(rows, indent=1, sort_keys=True).encode(), EXT)
    doc = {"id": cid, "qtl": entry["qtl"], "gwas": entry["gwas"], "source": entry["source"],
           "threshold": entry["threshold"], "note": entry.get("note"), "object": obj, "n": len(rows)}
    store.write_pointer(COLOCS, cid, doc)
    return doc


def add_configured(store: qs.Store, cfg: Config, only: str | None = None, log=print) -> list[dict]:
    """Every configured list whose two experiments are in the store (or just `only`)."""
    entries = cfg.cfg.get("colocs") or {}
    if only and only not in entries:
        raise ValueError(f"colocs.{only} is not in config.yaml")
    out = []
    for cid, entry in sorted(entries.items()):
        if only and cid != only:
            continue
        have = [(store.root / "experiments" / f"{entry.get(k)}.json").exists() for k in ("qtl", "gwas")]
        if not all(have):
            if only:
                raise ValueError(f"colocs.{cid}: {entry.get('qtl')} or {entry.get('gwas')} is not in the store")
            continue
        doc = add(store, cid, entry)
        log(f"coloc list {cid}: {doc['n']} rows, {entry['qtl']} x {entry['gwas']}")
        out.append(doc)
    return out


def load_rows(store: qs.Store, doc: dict) -> list[dict]:
    return json.loads((store.immutable / doc["object"]).read_bytes())


def check(store: qs.Store, cid: str, doc: dict, docs: dict) -> list[str]:
    """Check 15 (SPEC.md section 20): the pair exists with the right kinds, the rows decode, are
    unique and sorted, number `n`, and every phenotype and gene is one the QTL experiment has. That
    last part is what makes a list stale-proof: a rebuilt experiment that lost a phenotype fails here
    rather than leaving a list that names something no longer there."""
    what = f"{COLOCS}/{cid}"
    qtl, gwas = docs.get(("experiments", doc.get("qtl"))), docs.get(("experiments", doc.get("gwas")))
    fails = []
    if qtl is None or qtl.get("kind") not in qs.QTL_KINDS:
        fails.append(f"{what}: qtl {doc.get('qtl')!r} is not a QTL experiment in the store")
    if gwas is None or gwas.get("kind") != "gwas":
        fails.append(f"{what}: gwas {doc.get('gwas')!r} is not a GWAS experiment in the store")
    for k in ("source", "threshold"):
        if not doc.get(k):
            fails.append(f"{what}: {k} missing")
    if fails:
        return fails
    try:
        rows = load_rows(store, doc)
    except Exception as e:                 # noqa: BLE001 -- any decode failure is a validation failure
        return [f"{what}: {e}"]
    keys = [(r.get("gene_id"), r.get("phenotype_id") or "") for r in rows]
    if len(rows) != doc.get("n"):
        fails.append(f"{what}: {len(rows)} rows, the pointer says {doc.get('n')}")
    if keys != sorted(set(keys)):
        fails.append(f"{what}: rows not sorted and unique")
    from . import results
    index = results.load_index(store, qtl).select(["phenotype_id", "gene_id"]).to_pydict()
    phenos = dict(zip(index["phenotype_id"], index["gene_id"]))
    genes = set(g for g in index["gene_id"] if g)
    for r in rows:
        g, p = r.get("gene_id"), r.get("phenotype_id")
        if not isinstance(g, str) or g not in genes:
            fails.append(f"{what}: gene {g} has no phenotype in {doc['qtl']}")
        elif p is not None and phenos.get(p) != g:
            fails.append(f"{what}: phenotype {p} is not a phenotype of gene {g} in {doc['qtl']}")
    return fails[:20]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a_ = sub.add_parser("add")
    a_.add_argument("--store", required=True)
    a_.add_argument("--id", help="one configured list; default every one whose experiments are in the store")
    s = sub.add_parser("show")
    s.add_argument("--store", required=True)
    s.add_argument("--id", required=True)
    a = ap.parse_args(argv)
    store = qs.Store(Path(a.store))
    if a.cmd == "add":
        # object, then pointer; store.json is rewritten by whoever finishes the build, as for overlaps
        docs = add_configured(store, Config(), a.id, log=lambda m: print(m, file=sys.stderr))
        print(json.dumps(docs, indent=1))
    else:
        doc = store.load(COLOCS, a.id)
        print(json.dumps({**doc, "rows": load_rows(store, doc)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
