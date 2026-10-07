"""Gene symbols -> unversioned ENSG, for sources that publish no gene id.

    uv run python -m pipeline.genemap --symbols aric    # the report for one adapter's symbols

`CONTRACT.md` wants `gene_id` as an unversioned ENSG. The eQTL Catalogue and TOPCHeF publish one;
affinity-proteomics studies often do not -- ARIC names a SOMAmer and a gene **symbol**, and a symbol
is not an identity: it is renamed, it is reused, and two annotation releases disagree about it.

So this resolves a symbol in a stated order, and records which rule fired for every symbol rather
than silently picking one. An unresolved symbol becomes a null `gene_id`, which the contract allows;
it means the phenotype is invisible on a gene page, so the count matters and is reported.

The order, most to least trustworthy:

1. **`gencode`** -- the symbol is a `gene_name` in the annotation the experiment will be attached
   to, and exactly one gene carries it. This is the only rule that needs no external authority.
2. **`uniprot`** -- the source published a UniProt accession for the assay, and exactly one approved
   HGNC record lists it, whose `ensembl_gene_id` the annotation has. An accession is an identity,
   so this outranks every rule below, all of which reason from a *name*.
3. **`hgnc_current`** -- several annotation genes carry the symbol, and HGNC's own record for it
   names an `ensembl_gene_id` that is one of them. HGNC breaks the tie; the annotation still decides
   the candidate set.
4. **`hgnc_prev`** / **`hgnc_alias`** -- the annotation has no such symbol, but HGNC records it as a
   previous or alias symbol of exactly one approved gene, and that gene's `ensembl_gene_id` is in
   the annotation. This is the retired-symbol case (`C10orf54` -> VSIR).
5. **`ambiguous`** -- several candidates and nothing settles it. Null, listed in the report.
6. **`absent`** -- no candidate at all. Null, listed in the report.

A previous symbol is preferred over an alias because HGNC's `prev_symbol` is a statement about this
gene's own history, while `alias_symbol` is a looser "also called", which several genes can share.
Where both fire the result is the same gene or it is `ambiguous`.

Withdrawn HGNC records are skipped: they carry symbols that were reassigned.

Why the accession rule earns its place, from ARIC: the symbol `PACAP` reads as ADCYAP1 (pituitary
adenylate cyclase activating polypeptide) to a person and to an alias table, but the assay's
accession is Q8WU39, which is **MZB1** -- a gene whose previous symbol was also PACAP. Resolving
that symbol by name would have silently attached a plasma protein to the wrong gene, on a page that
would have looked entirely normal. ARIC also carries three Excel-mangled symbols (`3-Sep`,
`10-Sep`, `11-Sep` for the septins) and four immunoglobulin assays named by a list of genes
(`IGHG1 IGHG2 IGHG3 IGHG4 IGK@ IGL@`); the accession resolves all seven to one gene.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import pyarrow as pa

from .common import Config

HOW = ("override", "gencode", "uniprot", "hgnc_current", "hgnc_prev", "hgnc_alias", "ambiguous", "absent")

# Calls a person made, because the annotation is ambiguous and no outside authority settles it
# *within that annotation*. Keyed by (symbol, uniprot); the target must exist in the annotation or
# the resolver raises, so a stale entry fails the build instead of going quiet.
#
# These are annotation-version specific by nature. GENCODE v34 carries two overlapping
# protein-coding genes named SOD2 -- ENSG00000112096 (chr6:159,669,069-159,745,186) and
# ENSG00000285441 (chr6:159,679,119-159,762,529) -- and HGNC resolves the symbol and the accession
# P04179 to ENSG00000291237, an id v34 does not have at all. So nothing external can pick between
# the two the annotation offers. ENSG00000112096 is the long-standing SOD2 id, the one GTEx and the
# literature use; ENSG00000285441 appeared later and overlaps it.
OVERRIDES: dict[tuple[str, str | None], str] = {
    ("SOD2", "P04179"): "ENSG00000112096",
}


def _split(field: str) -> list[str]:
    """One HGNC multi-valued cell. They are pipe-separated and sometimes double-quoted."""
    return [x.strip().strip('"') for x in field.split("|") if x.strip().strip('"')]


class Hgnc:
    """The HGNC complete set, reduced to what a symbol lookup needs."""

    def __init__(self, path: Path):
        self.ensg_of_symbol: dict[str, set[str]] = defaultdict(set)
        self.ensg_of_prev: dict[str, set[str]] = defaultdict(set)
        self.ensg_of_alias: dict[str, set[str]] = defaultdict(set)
        self.ensg_of_uniprot: dict[str, set[str]] = defaultdict(set)
        self.rows = 0
        self.withdrawn = 0
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                self.rows += 1
                if (r.get("status") or "").lower() != "approved":
                    self.withdrawn += 1
                    continue
                ensg = (r.get("ensembl_gene_id") or "").strip()
                if not ensg:
                    continue
                sym = (r.get("symbol") or "").strip()
                if sym:
                    self.ensg_of_symbol[sym].add(ensg)
                for s in _split(r.get("prev_symbol") or ""):
                    self.ensg_of_prev[s].add(ensg)
                for s in _split(r.get("alias_symbol") or ""):
                    self.ensg_of_alias[s].add(ensg)
                for s in _split(r.get("uniprot_ids") or ""):
                    self.ensg_of_uniprot[s].add(ensg)


class Resolver:
    """Symbols against one annotation's genes, with HGNC as the tie-break and the alias source.

    `genes` is the annotation's genes table (`annotation.parse_gtf`'s `genes`, or the stored
    annotation's): it needs `gene_id` and the symbol, which the annotation object calls `name` and
    the v0 tables call `symbol`. The annotation decides what an ENSG *is* here --
    an HGNC record naming a gene the annotation does not have resolves to nothing, because a
    `gene_id` the annotation lacks would be a dangling reference in the store.
    """

    def __init__(self, genes: pa.Table, hgnc: Hgnc | None = None, overrides: dict | None = None):
        self.hgnc = hgnc
        self.overrides = OVERRIDES if overrides is None else overrides
        self.by_symbol: dict[str, set[str]] = defaultdict(set)
        self.known: set[str] = set()
        col = "name" if "name" in genes.schema.names else "symbol"
        for gid, sym in zip(genes.column("gene_id").to_pylist(), genes.column(col).to_pylist()):
            if not gid:
                continue
            self.known.add(gid)
            if sym:
                self.by_symbol[sym].add(gid)
        stale = {k: v for k, v in self.overrides.items() if v not in self.known}
        if stale:
            raise ValueError(f"genemap overrides name genes this annotation does not have: {stale}")

    def _in_annotation(self, ensgs) -> set[str]:
        return {e for e in ensgs if e in self.known}

    def resolve(self, symbol: str | None, uniprot: str | None = None) -> tuple[str | None, str]:
        """`(gene_id, how)`; `how` is one of `HOW`. `gene_id` is null unless exactly one survived."""
        over = self.overrides.get((symbol, uniprot)) or self.overrides.get((symbol, None))
        if over:
            return over, "override"
        cands = self.by_symbol.get(symbol or "", set())
        if len(cands) == 1:
            return next(iter(cands)), "gencode"
        h = self.hgnc
        if h and uniprot:
            pick = self._in_annotation(h.ensg_of_uniprot.get(uniprot.strip(), set()))
            if len(cands) > 1:
                pick &= cands              # the annotation's candidate set still decides
            if len(pick) == 1:
                return next(iter(pick)), "uniprot"
        if len(cands) > 1:
            # several annotation genes share the symbol; HGNC's own record may name one of them
            pick = self._in_annotation(h.ensg_of_symbol.get(symbol, set())) & cands if h else set()
            return (next(iter(pick)), "hgnc_current") if len(pick) == 1 else (None, "ambiguous")
        if not symbol:
            return None, "absent"
        if h:
            for how, table in (("hgnc_prev", h.ensg_of_prev), ("hgnc_alias", h.ensg_of_alias)):
                pick = self._in_annotation(table.get(symbol, set()))
                if len(pick) == 1:
                    return next(iter(pick)), how
                if len(pick) > 1:
                    return None, "ambiguous"
        return None, "absent"

    def report(self, pairs) -> dict:
        """Counts per rule plus the unresolved keys, for `ingestion.json`. An adapter stores this so
        the number of phenotypes with no gene is a recorded fact, not a surprise later.

        `pairs` is `(symbol, uniprot)` tuples, or bare symbols. `map` is keyed by the pair, since two
        assays can share a symbol and disagree on the accession."""
        counts = dict.fromkeys(HOW, 0)
        unresolved: dict[str, list[str]] = {"ambiguous": [], "absent": []}
        resolved: dict[tuple[str | None, str | None], str | None] = {}
        for p in pairs:
            key = (p, None) if isinstance(p, str) or p is None else (p[0], p[1])
            if key in resolved:
                continue
            gid, how = self.resolve(key[0], key[1])
            counts[how] += 1
            resolved[key] = gid
            if gid is None:
                unresolved[how].append(key[0] or "")
        return {"counts": counts, "keys": len(resolved),
                "resolved": sum(1 for v in resolved.values() if v),
                "unresolved": {k: sorted(set(v)) for k, v in unresolved.items()}, "map": resolved}


def load(cfg: Config, gtf: Path | None = None, hgnc: Path | None = None) -> Resolver:
    """A resolver from the configured GENCODE GTF and HGNC file."""
    from . import annotation
    gtf = gtf or cfg.gtf
    hgnc = hgnc if hgnc is not None else (cfg.raw / "hgnc" / "hgnc_complete_set.txt")
    if not Path(hgnc).exists():
        raise SystemExit(f"no HGNC set at {hgnc}; fetch it with `./data/raw/download.py --only hgnc`")
    return Resolver(annotation.parse_gtf(gtf)["genes"], Hgnc(hgnc))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--symbols", required=True, help="`aric` for ARIC's seqid.txt, or a file with one symbol per line")
    ap.add_argument("--gtf", type=Path)
    ap.add_argument("--hgnc", type=Path)
    ap.add_argument("--examples", type=int, default=12)
    args = ap.parse_args(argv)
    cfg = Config()
    if args.symbols == "aric":
        p = cfg.raw / "aric_pqtl_zhang2022" / "seqid.txt"
        with open(p, newline="", encoding="utf-8") as fh:
            syms = [(r["entrezgenesymbol"], r["uniprot_id"]) for r in csv.DictReader(fh, delimiter="\t")]
    else:
        syms = [l.strip() for l in Path(args.symbols).read_text().splitlines() if l.strip()]
    r = load(cfg, args.gtf, args.hgnc)
    rep = r.report(syms)
    rep.pop("map")
    rep["unresolved"] = {k: v[:args.examples] for k, v in rep["unresolved"].items()}
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
