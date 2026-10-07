"""Sweep every results block's codes: the reserved values, the scale rules, and the stated limits.

`Store.validate` checks that each object's bytes hash to its name and that headers and pointers
agree, but it never looks at the 2n u16 codes inside a block -- there are hundreds of millions of
them and decoding each block fully (details frame, credible sets, rebuilt slopes) to count four
things would be absurd. This sweeps the codes directly instead, which is fast because a block
stores them raw -- v1 as `nlp_code, se_code` interleaved after a 64-byte header, v2 as three
columns (`nlp` u16, `beta` i16, `se` u16) after a 72-byte one. The store's `format_version` picks
which.

What it counts, per (phenotype_type, chromosome):

    p_null      -log10 p code 65535: the phenotype was never tested against that variant, even
                though the block's vidx range is contiguous (SPEC section 8)
    p_zero      code 65534: p underflowed to 0 in the source, so no t statistic and no z exists.
                Readers must drop these rows, and a dropped row is among the *strongest* signals
                in its window -- see `ui/src/lib/coloc-abf.ts`, where it would move the posteriors
    se_null     SE code 0xFFFF: no standard error
    both        rows that are p_null and se_null at once
    beta_null   v2 only: beta code -32768, no effect size. v1 has no such code -- it *derives*
                beta, so a row with no effect size is not something v1 can express

Reserved codes are legitimate data, not failures. What *is* a failure is a scale rule, because each
header scale is defined as a particular row's value and a violation means every decoded number in
the block is wrong:

    bad_scale   the largest finite -log10 p code is not 65533 (`nlp_max` is that row's value)
    bad_beta    v2: the largest |beta| code is not 32766, with `beta_max` > 0
    bad_se      v2: the SE codes do not reach both 0 and 65534, with `lse_max` > `lse_min`

And one check that is not about codes at all: every `precision` block's `measured_worst` must be
within its own `theoretical_limit`. A bad scale is the first thing that would break it, which is
why it belongs in the same pass -- the limit is computed from the stored scales and the measurement
from the source rows, so agreement means the two halves of the claim were derived independently.

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


FIELDS = ("blocks", "rows", "p_null", "p_zero", "se_null", "both", "beta_null",
          "bad_scale", "bad_beta", "bad_se")


def block_codes(buf: bytes, blk_off: int, n_var: int, codec=pf) -> dict[str, np.ndarray | None]:
    """One block's code arrays, as views into `buf`.

    v1 interleaves `nlp_code, se_code` as `2n` u16s after its 64-byte header. v2 lays its three
    codes out as columns -- `nlp` u16, `beta` i16, `se` u16 -- after a 72-byte one, which is why
    this returns a dict rather than a pair: `beta` exists only in v2."""
    if codec.FORMAT_VERSION == 1:
        u16 = np.frombuffer(buf, dtype="<u2", count=2 * n_var, offset=blk_off + 64)
        return {"nlp": u16[0::2], "se": u16[1::2], "beta": None}
    o = blk_off + codec.BLOCK_HEADER_LEN
    return {"nlp": np.frombuffer(buf, "<u2", n_var, o),
            "beta": np.frombuffer(buf, "<i2", n_var, o + 2 * n_var),
            "se": np.frombuffer(buf, "<u2", n_var, o + 4 * n_var)}


def scan_file(buf: bytes, blocks: list[tuple[int, int]], codec=pf) -> dict:
    """Counts over every block of one results object. `blocks` is `(blk_off, n_var)` per phenotype."""
    c = dict.fromkeys(FIELDS, 0)
    c["blocks"] = len(blocks)
    v2 = codec.FORMAT_VERSION != 1
    for off, n in blocks:
        col = block_codes(buf, off, n, codec)
        nlp, se = col["nlp"], col["se"]
        pn, se_n = nlp == codec.NLP_NULL, se == codec.SE_NULL
        finite = nlp[nlp <= codec.NLP_MAXQ]
        if finite.size and int(finite.max()) != codec.NLP_MAXQ:
            c["bad_scale"] += 1
        c["rows"] += n
        c["p_null"] += int(pn.sum())
        c["p_zero"] += int((nlp == codec.NLP_ZERO).sum())
        c["se_null"] += int(se_n.sum())
        c["both"] += int((pn & se_n).sum())
        if not v2:
            continue
        beta = col["beta"]
        bn = beta == codec.BETA_NULL
        c["beta_null"] += int(bn.sum())
        # `beta_max` is the largest |beta| over the non-null rows, so some row must code 32766.
        # An all-null or all-zero block has beta_max = 0 and codes nothing, which is not a breach.
        live = beta[~bn]
        if live.size and int(np.abs(live.astype(np.int32)).max()) not in (0, codec.BETA_MAXQ):
            c["bad_beta"] += 1
        # se is a log ruler between lse_min and lse_max, so the extremes code 0 and 65534 -- unless
        # every live se is equal, when the span is zero and every code is 0.
        live_se = se[~se_n]
        if live_se.size:
            lo, hi = int(live_se.min()), int(live_se.max())
            if not (lo == hi == 0 or (lo == 0 and hi == codec.SE_MAXQ)):
                c["bad_se"] += 1
    return c


def check_limits(doc: dict) -> list[str]:
    """Every `precision` block's `measured_worst` within its own `theoretical_limit`, over the cis and
    trans entries of an experiment pointer.

    The limit comes from the stored header scales and the measurement from the source rows at build
    time, so this compares two independently derived halves of the same claim. v1 states no limit
    for `beta_over_se` -- it never bounded the quantity it breaches -- and a 0.0 limit there is read
    as "not stated" rather than as a breach of zero."""
    out, unstated = [], []
    for r in doc.get("results") or []:
        for scope, entry in (("cis", r), ("trans", r.get("trans") or {})):
            p = entry.get("precision") or {}
            lim, got = p.get("theoretical_limit") or {}, p.get("measured_worst") or {}
            for field, g in got.items():
                if field == "rows_compared" or field not in lim:
                    continue
                where = f"{r.get('phenotype_type')} {scope} {field}"
                if lim[field] == 0.0:
                    if g > 0.0:
                        unstated.append(f"{where}: measured {g:.6g}, no limit stated")
                    continue
                if g > lim[field]:
                    out.append(f"{where}: measured {g:.6g} exceeds the limit {lim[field]:.6g}")
    return {"breaches": out, "unstated": unstated}


def scan(objects: Objects, exp_id: str, chroms: list[str] | None = None, log=lambda m: None,
         codec=None) -> dict:
    """Every nominal block of one experiment, keyed `(phenotype_type, chromosome)`. `chroms` limits
    it to some chromosomes, which is the difference between a spot check and a 2.5 GB read.

    `codec` defaults to whatever the objects' own file headers declare, so the right code layout is
    read without being told. Reading a v2 block at the v1 layout would not raise -- it would count
    beta codes as if they were SE codes and report a store full of reserved values."""
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
        buf = objects(files[ptype][chrom])
        cd = codec or qs.codec_for(qs.parse_file_header(buf)["version"])
        out[key] = c = scan_file(buf, want[key], cd)
        log(f"v{cd.FORMAT_VERSION} {ptype:11s} {chrom:6s} rows {c['rows']:>12,}  p_null {c['p_null']:>7,}  "
            f"p_zero {c['p_zero']:>7,}  se_null {c['se_null']:>7,}  both {c['both']:>7,}  "
            f"beta_null {c['beta_null']:>7,}  bad {c['bad_scale']}/{c['bad_beta']}/{c['bad_se']}")
    return out


def totals(scanned: dict, ptype: str | None = None) -> dict:
    keys = [k for k in scanned if ptype is None or k[0] == ptype]
    return {f: sum(scanned[k][f] for k in keys) for f in FIELDS}


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
              f"both {t['both']:>7,}  beta_null {t['beta_null']:>7,}  "
              f"bad_scale {t['bad_scale']}  bad_beta {t['bad_beta']}  bad_se {t['bad_se']}")
    t = totals(scanned)
    broken = {k: t[k] for k in ("bad_scale", "bad_beta", "bad_se") if t[k]}
    lim = check_limits(objects.pointer("experiments", a.experiment))
    print()
    for b in lim["breaches"]:
        print(f"LIMIT BREACH: {b}", file=sys.stderr)
    # A measurement with no limit to compare it to is the v1 asymmetry, not a pass: the figure that
    # sits at 99.9% of v1's slope budget is one of these, and reporting it as "0 breaches" would
    # repeat exactly the silence that let it ship.
    for u in lim["unstated"]:
        print(f"NO LIMIT STATED: {u}")
    n = len(lim["breaches"])
    print(f"precision: {n or 'no'} measured value{'' if n == 1 else 's'} outside a stated limit, "
          f"{len(lim['unstated'])} with no limit to compare against")
    if broken:
        print(f"FAIL: {broken}", file=sys.stderr)
    return 1 if (broken or lim["breaches"]) else 0


if __name__ == "__main__":
    sys.exit(main())
