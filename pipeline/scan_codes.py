"""Sweep every results block's u16 codes: the reserved -log10 p and SE values, and the scale rule.

`Store.validate` checks that each object's bytes hash to its name and that headers and pointers
agree, but it never looks at the 2n u16 codes inside a block -- there are hundreds of millions of
them and decoding each block fully (details frame, credible sets, rebuilt slopes) to count four
things would be absurd. This sweeps the codes directly instead, which is fast because a v1 block
stores them raw: 64-byte header, then `nlp_code, se_code` interleaved, one pair per row.

What it counts, per (phenotype_type, chromosome):

    p_null      -log10 p code 65535: the phenotype was never tested against that variant, even
                though the block's vidx range is contiguous (SPEC section 8)
    p_zero      code 65534: p underflowed to 0 in the source, so no t statistic and no z exists.
                Readers must drop these rows, and a dropped row is among the *strongest* signals
                in its window -- see `ui/src/lib/coloc-abf.ts`, where it would move the posteriors
    se_null     SE code 0xFFFF: no standard error
    both        rows that are p_null and se_null at once

Reserved codes are legitimate data, not failures. The one thing here that *is* a failure is the
scale rule: a block's largest finite -log10 p code must be exactly 65533, since `nlp_max` in the
header is defined as that row's value. A violation means the codes and the header disagree and
every -log10 p in the block decodes to the wrong number, so the sweep exits non-zero on one.

The published store lives on B2 rather than on disk, so `--base` reads objects over HTTP:

    uv run python -m pipeline.scan_codes --store DIR --experiment topchef
    uv run python -m pipeline.scan_codes --base https://cloud2.databio.org/qtl-browser --experiment topchef
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

from . import packfmt_v1 as pf, qtlstore as qs, results


class Objects:
    """Object bytes by name, from a store directory or over HTTP."""

    def __init__(self, store: Path | None, base: str | None):
        self.store = Path(store) if store else None
        self.base = base.rstrip("/") if base else None

    def __call__(self, name: str) -> bytes:
        if self.store:
            return (self.store / "immutable" / name).read_bytes()
        # cloud2 refuses urllib's default User-Agent with a 403
        req = urllib.request.Request(f"{self.base}/immutable/{name}", headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req) as r:
            return r.read()

    def pointer(self, level: str, pid: str) -> dict:
        if self.store:
            return json.loads((self.store / level / f"{pid}.json").read_text())
        req = urllib.request.Request(f"{self.base}/{level}/{pid}.json", headers={"User-Agent": "curl/8"})
        with urllib.request.urlopen(req) as r:
            return json.loads(r.read())


def block_codes(buf: bytes, blk_off: int, n_var: int) -> tuple[np.ndarray, np.ndarray]:
    """One block's `(nlp_code, se_code)` arrays, as views into `buf`. The codes are the `2 * n_var`
    u16s after the block's 64-byte header, interleaved one pair per row."""
    u16 = np.frombuffer(buf, dtype="<u2", count=2 * n_var, offset=blk_off + 64)
    return u16[0::2], u16[1::2]


def scan_file(buf: bytes, blocks: list[tuple[int, int]]) -> dict:
    """Counts over every block of one results object. `blocks` is `(blk_off, n_var)` per phenotype."""
    c = dict(blocks=len(blocks), rows=0, p_null=0, p_zero=0, se_null=0, both=0, bad_scale=0)
    for off, n in blocks:
        nlp, se = block_codes(buf, off, n)
        pn, se_n = nlp == pf.NLP_NULL, se == pf.SE_NULL
        finite = nlp[nlp <= pf.NLP_MAXQ]
        if finite.size and int(finite.max()) != pf.NLP_MAXQ:
            c["bad_scale"] += 1
        c["rows"] += n
        c["p_null"] += int(pn.sum())
        c["p_zero"] += int((nlp == pf.NLP_ZERO).sum())
        c["se_null"] += int(se_n.sum())
        c["both"] += int((pn & se_n).sum())
    return c


def scan(objects: Objects, exp_id: str, chroms: list[str] | None = None, log=lambda m: None) -> dict:
    """Every nominal block of one experiment, keyed `(phenotype_type, chromosome)`. `chroms` limits
    it to some chromosomes, which is the difference between a spot check and a 2.5 GB read."""
    doc = objects.pointer("experiments", exp_id)
    if objects.store:
        index = results.load_index(qs.Store(objects.store), doc)
    else:
        index = results.decode_arrow(objects(doc["search_index"]), "search index")
    cols = index.to_pydict()
    # blocks live in the object for their (type, chromosome), addressed by the index row's blk_off
    want: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for ptype, chrom, off, n, nominal in zip(cols["phenotype_type"], cols["chr"], cols["blk_off"],
                                             cols["n_var"], cols["has_nominal"]):
        if nominal and n and (chroms is None or chrom in chroms):
            want.setdefault((ptype, chrom), []).append((off, n))
    files = {r["phenotype_type"]: r["files"] for r in doc["results"]}
    out = {}
    for key in sorted(want, key=lambda k: (k[0], list(files[k[0]]).index(k[1]))):
        ptype, chrom = key
        out[key] = c = scan_file(objects(files[ptype][chrom]), want[key])
        log(f"{ptype:11s} {chrom:6s} rows {c['rows']:>12,}  p_null {c['p_null']:>7,}  "
            f"p_zero {c['p_zero']:>7,}  se_null {c['se_null']:>7,}  both {c['both']:>7,}  "
            f"bad_scale {c['bad_scale']}")
    return out


def totals(scanned: dict, ptype: str | None = None) -> dict:
    keys = [k for k in scanned if ptype is None or k[0] == ptype]
    return {f: sum(scanned[k][f] for k in keys) for f in
            ("blocks", "rows", "p_null", "p_zero", "se_null", "both", "bad_scale")}


# ---- CLI --------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--store", help="a store directory")
    src.add_argument("--base", help="a store's base URL, e.g. https://cloud2.databio.org/qtl-browser")
    ap.add_argument("--experiment", required=True)
    ap.add_argument("--chrom", nargs="+", help="only these chromosomes; default every one")
    ap.add_argument("--quiet", action="store_true", help="totals only")
    a = ap.parse_args(argv)
    objects = Objects(a.store, a.base)
    scanned = scan(objects, a.experiment, a.chrom,
                   log=(lambda m: None) if a.quiet else lambda m: print(m, flush=True))
    print()
    for ptype in sorted({k[0] for k in scanned}) + [None]:
        t = totals(scanned, ptype)
        print(f"{ptype or 'ALL':11s} blocks {t['blocks']:>7,}  rows {t['rows']:>12,}  "
              f"p_null {t['p_null']:>7,}  p_zero {t['p_zero']:>7,}  se_null {t['se_null']:>7,}  "
              f"both {t['both']:>7,}  bad_scale {t['bad_scale']}")
    bad = totals(scanned)["bad_scale"]
    if bad:
        print(f"\nFAIL: {bad} blocks break the -log10 p scale rule", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
