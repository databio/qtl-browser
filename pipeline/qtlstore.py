"""qtlstore core (qtlb v1): content-addressed objects, the v1 file header, variant catalog identity, allele
orientation, and a store validator.

A store is a directory:

    store.json                     mutable: name, format version, experiment and variant catalog ids, refget store URL(s)
    immutable/<sha512t24u>.<ext>   every pack, index, variant catalog and annotation object; name = digest of its bytes
    variant_catalogs/<id>.json             mutable pointer: the objects that make up one variant catalog
    annotations/<id>.json          mutable pointer: one gene/transcript annotation set
    experiments/<id>.json          mutable pointer: variant catalog + annotation + results objects + constants

Objects are written first and pointers last, so a reader never sees a pointer to a missing object.
The digest is refget's sha512t24u: sha512, first 24 bytes, base64url, 32 characters.
"""
from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import struct
from pathlib import Path

import numpy as np

FORMAT_VERSION = 1

# ---- digest -----------------------------------------------------------------------------------
DIGEST = re.compile(r"^[A-Za-z0-9_-]{32}$")
OBJECT_NAME = re.compile(r"^(?P<digest>[A-Za-z0-9_-]{32})\.(?P<ext>[a-z0-9]+(?:\.[a-z0-9]+)*)$")


def sha512t24u(data: bytes) -> str:
    """refget sha512t24u: base64url of the first 24 bytes of sha512 (same as gtars.refget.sha512t24u_digest)."""
    return base64.urlsafe_b64encode(hashlib.sha512(data).digest()[:24]).decode("ascii")


def sha512t24u_file(path: Path) -> str:
    h = hashlib.sha512()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return base64.urlsafe_b64encode(h.digest()[:24]).decode("ascii")


# ---- v1 file header ---------------------------------------------------------------------------
# magic, kind, version, header length, chromosome (8 ASCII), count, page size, n_cis,
# seq_digest (32 ASCII sha512t24u of the chromosome sequence), 4 reserved zero bytes
MAGIC_FILE = b"QTLB"
HEADER_LEN = 64
_HEADER = struct.Struct("<4sBBH8sIII32s4x")
assert _HEADER.size == HEADER_LEN
U32_MAX = 0xFFFFFFFF
# header kinds v1 writes (SPEC.md section 3); 3 is unused (v0's sQTL kind)
KIND_VARIANTS, KIND_RESULTS, KIND_GWAS, KIND_GWAS_INDEX, KIND_TRANS = 1, 2, 4, 5, 6
KIND_HITS, KIND_RSID, KIND_VARIANT_INDEX, KIND_GENE_LOOKUP = 7, 8, 9, 10
ALL = "all"                        # header chromosome of objects spanning a whole variant catalog (seq_digest: the collection) or annotation (its identity)


def file_header(kind: int, chrom: str, count: int, page_size: int, n_cis: int, seq_digest: str) -> bytes:
    """64 bytes. `seq_digest` names the exact sequence the file's positions index, so a pack found
    alone still says which genome it belongs to."""
    if not 1 <= kind <= 255:
        raise ValueError(f"file header: kind {kind} not in 1..255")
    try:
        name = chrom.encode("ascii")
    except UnicodeEncodeError:
        raise ValueError(f"file header: chromosome {chrom!r} is not ASCII") from None
    if not 1 <= len(name) <= 8 or b"\0" in name:
        raise ValueError(f"file header: chromosome {chrom!r} must be 1-8 ASCII bytes without NUL")
    for field, v in (("count", count), ("page size", page_size), ("n_cis", n_cis)):
        if not 0 <= v <= U32_MAX:
            raise ValueError(f"file header: {field} {v} does not fit u32")
    if not DIGEST.match(seq_digest):
        raise ValueError(f"file header: seq_digest {seq_digest!r} is not a 32-character sha512t24u")
    return _HEADER.pack(MAGIC_FILE, kind, FORMAT_VERSION, HEADER_LEN, name.ljust(8, b"\0"), count, page_size, n_cis,
                        seq_digest.encode("ascii"))


def parse_file_header(b: bytes) -> dict:
    if len(b) < HEADER_LEN:
        raise ValueError(f"file header: {len(b)} bytes, need {HEADER_LEN}")
    magic, kind, version, hlen, name, count, page_size, n_cis, sd = _HEADER.unpack_from(b, 0)
    if magic != MAGIC_FILE:
        raise ValueError(f"file header: magic {magic!r} is not {MAGIC_FILE!r}")
    if version != FORMAT_VERSION:
        raise ValueError(f"file header: version {version}, this reader supports {FORMAT_VERSION}")
    if hlen != HEADER_LEN:
        raise ValueError(f"file header: header length {hlen} != {HEADER_LEN}")
    if any(b[60:64]):
        raise ValueError("file header: reserved bytes are not zero")
    chrom = name.rstrip(b"\0")
    if not chrom or b"\0" in chrom or max(chrom) >= 0x80:
        raise ValueError(f"file header: bad chromosome field {name!r}")
    seq_digest = sd.decode("ascii", "replace")
    if not DIGEST.match(seq_digest):
        raise ValueError(f"file header: bad seq_digest field {sd!r}")
    return {"kind": kind, "version": version, "chrom": chrom.decode("ascii"), "count": count,
            "page_size": page_size, "n_cis": n_cis, "seq_digest": seq_digest}


# ---- variant catalog identity -----------------------------------------------------------------
def catalog_identity(chromosomes, load=None) -> str:
    """sha512t24u over `seq_digest\tpos\tref\talt\n` for every site in **canonical order**: chromosomes
    by `seq_digest` (ASCII, as seqcol's `sorted_sequences`), then sites by (pos, ref, alt) byte-wise.

    The identity depends only on the set of sites. File order does not enter it, so the cis/trans-only
    split, the chromosome table's order, the chromosome names, the page size and the attribute columns
    all stay out: the same sites on the same sequences give the same identity.

    `chromosomes` is an iterable of (seq_digest, pos, ref, alt), one entry per chromosome, with
    pos/ref/alt array-likes. With `load`, it is (seq_digest, key) instead and `load(key) -> (pos, ref,
    alt)` is called one chromosome at a time, in canonical order, so a caller need not hold every
    chromosome's sites at once."""
    entries = sorted(chromosomes, key=lambda e: e[0])
    digests = [e[0] for e in entries]
    if len(set(digests)) != len(digests):
        raise ValueError("catalog identity: a seq_digest appears in more than one chromosome entry")
    h = hashlib.sha512()
    for e in entries:
        seq_digest, (pos, ref, alt) = e[0], (load(e[1]) if load is not None else e[1:])
        prefix = seq_digest + "\t"
        sites = sorted(zip((int(p) for p in pos), ref, alt))
        h.update("".join(f"{prefix}{p}\t{r}\t{a}\n" for p, r, a in sites).encode("ascii"))
    return base64.urlsafe_b64encode(h.digest()[:24]).decode("ascii")


# ---- allele orientation -----------------------------------------------------------------------
def orient_to_ref(ref_base, effect_allele, other_allele, beta, af) -> dict:
    """Reorient a site to ref/alt with ref = the reference base(s) and `beta`/`af` ALT-relative, the store's
    one orientation convention (CONTRACT.md "Orientation").

    `beta` and `af` must describe `effect_allele` -- that is what the parameter name is for. The old names
    `a1`/`a2` did not say it, and a source's column letters are not a convention. So:

      * `other_allele == REF`  -> the effect allele is already ALT: nothing changes.
      * `effect_allele == REF` -> alleles swap, `beta` negates, `af` becomes `1 - af`.
      * neither matches        -> the site is dropped; it cannot be joined across studies or given a VRS id.

    **For TOPCHeF this is a no-op on 98.9% of sites**, and the plain call is the right one:
    `orient_to_ref(ref, effect_allele=A1, other_allele=A2, beta=slope, af=af)`. A2 is the reference allele
    and `af`/`slope` are A1's, so A1 is already ALT. Both halves are measured, not assumed:

      * *A2 is REF*: the refcheck reads GRCh38 at every variant and gets A2 at all 8,419,594 cis SNPs and
        A1 at none; `neither` is 0 across all 8,872,723 cis variants. The 97,433 variants called `a1` are
        all indels, where prefix matching cannot break the tie and `classify` takes the longer allele.
        EVIDENCE.md A.10 (analysis repo, `qtlb-format/docs/`); per-variant calls in
        `_tables/refcheck/<chr>.parquet`.
      * *`af` is the A1 frequency*: over the 6,793,566 biallelic SNPs this release shares with the Jurgens
        2024 DCM GWAS, `af` correlates **+0.9932** with that study's EAFREQ oriented to A1 (mean absolute
        difference 0.022); against the mirror `1 - eaf` it is -0.9932, mean absolute difference 0.668. On
        the 35,444 of those rows where the GWAS effect allele is our A2 -- so the frequency had to be
        mirrored to be compared -- the correlation is +0.9975, which is why this is not a tautology of the
        join. Strand-ambiguous pairs cannot confound it: that GWAS carries zero palindromic SNPs out of
        9,998,282, so our 1,282,725 palindromic sites never join. Cis `af` averages 0.2420 with median
        0.1318 over SNPs, a minor/ALT distribution, not a REF one. EVIDENCE.md A.11.

    `beta` rides along with `af` by construction: `steps_nominal` copies `n.af` and `n.slope` unchanged
    from one tensorQTL nominal row, and tensorQTL codes both against the same dosage allele.

    So reorienting a TOPCHeF v0 store to ref/alt is a relabelling. It is **not** a sign flip. A rebuild
    whose slopes come out negated genome-wide is wrong, not "oriented".
    """
    ref_base = np.asarray(ref_base, dtype=object)
    effect_allele = np.asarray(effect_allele, dtype=object)
    other_allele = np.asarray(other_allele, dtype=object)
    beta, af = np.asarray(beta, dtype=np.float64), np.asarray(af, dtype=np.float64)
    keep_as_is = other_allele == ref_base           # effect allele is ALT already: beta and af stand
    swap = ~keep_as_is & (effect_allele == ref_base)
    keep = keep_as_is | swap
    ref = np.where(swap, effect_allele, other_allele)
    alt = np.where(swap, other_allele, effect_allele)
    return {"ref": ref[keep], "alt": alt[keep], "beta": np.where(swap, -beta, beta)[keep],
            "af": np.where(swap, 1.0 - af, af)[keep], "keep": keep,
            "counts": {"as_is": int(keep_as_is.sum()), "swapped": int(swap.sum()), "dropped": int((~keep).sum())}}


# ---- store ------------------------------------------------------------------------------------
CATALOGS = "variant_catalogs"            # pointer level of variant catalogs (store.json key too)
POINTER_DIRS = (CATALOGS, "annotations", "experiments")


class Store:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.immutable = self.root / "immutable"
        self.notes = []                 # what the last validate() could not check, and why

    def put(self, src: Path | bytes, ext: str) -> str:
        """Copy bytes into immutable/ under their digest; returns the object name `<digest>.<ext>`."""
        data = src if isinstance(src, bytes) else Path(src).read_bytes()
        name = f"{sha512t24u(data)}.{ext}"
        if not OBJECT_NAME.match(name):
            raise ValueError(f"put: bad extension {ext!r}")
        out = self.immutable / name
        if not out.exists():
            self.immutable.mkdir(parents=True, exist_ok=True)
            tmp = out.with_name(name + ".tmp")
            tmp.write_bytes(data)
            os.replace(tmp, out)
        return name

    def write_pointer(self, level: str, pid: str, doc: dict) -> Path:
        if level not in POINTER_DIRS:
            raise ValueError(f"pointer level {level!r} not in {POINTER_DIRS}")
        if doc.get("id") != pid:
            raise ValueError(f"{level}/{pid}: document id {doc.get('id')!r} does not match")
        missing = [n for n in object_names(doc) if not (self.immutable / n).exists()]
        if missing:
            raise ValueError(f"{level}/{pid}: objects not in the store (write objects first): {missing[:3]}")
        return _write_json(self.root / level / f"{pid}.json", doc)

    def write_store(self, name: str, refget: list[str]) -> Path:
        """store.json last: it lists whatever pointers exist now."""
        ids = {lvl: sorted(p.stem for p in (self.root / lvl).glob("*.json")) for lvl in POINTER_DIRS}
        return _write_json(self.root / "store.json", {"name": name, "format_version": FORMAT_VERSION,
                                                      "refget": refget, **ids})

    def load(self, level: str, pid: str) -> dict:
        return json.loads((self.root / level / f"{pid}.json").read_text())

    def validate(self, refget=None, sites=None) -> list[str]:
        """Failures (empty = pass): every listed id exists, every object name resolves to bytes with that
        digest, every experiment's variant catalog and annotation exist and its `catalog_identity` equals the
        variant catalog's `identity_digest`, every binary object's header parses with the right kind -- variants,
        results and hits files with the chromosome and seq_digest of the variant catalog table entry, the variant
        and rsID indexes with chromosome `all` and the variant catalog's collection digest --, a GWAS bin summary
        decodes and passes `gwas.check_bins` on chromosomes that have GWAS files, every trans object's
        frames tile it and each gene's frames are one contiguous range (`results.check_trans_layout`), every variant catalog
        chromosome is a sequence the refgetstore holds at the recorded length, and every variant catalog identity
        digest recomputes from the variant catalog's own sites.

        `refget` is an open refgetstore, or a Config to open one from; `sites` reads one chromosome's sites
        back out of its pack, `sites(store, chrom) -> (pos, ref, alt)`.

        The first four checks only prove the store agrees with itself. A variant catalog whose every chromosome
        entry and every pack header carried the same wrong `seq_digest` would pass all of them, which
        defeats the anchor: the point is agreement with the reference, not internal agreement. Only the
        refgetstore can settle that, and only a site reader can settle the identity digest. Neither is
        something this module can conjure, so when one is absent the run says so in `self.notes` rather
        than passing quietly. **An empty `fails` with a non-empty `notes` is a partial pass**; a caller
        reporting one without the other is reporting a check that did not happen.

        Experiment result packs need no separate refgetstore check: their headers are compared to the
        variant catalog table, which is itself compared to the store.
        """
        self.notes = []
        store = self._refgetstore(refget)
        fails = []
        st = json.loads((self.root / "store.json").read_text())
        if st.get("format_version") != FORMAT_VERSION:
            fails.append(f"store.json: format_version {st.get('format_version')}")
        docs = {}
        for lvl in POINTER_DIRS:
            for pid in st.get(lvl, []):
                p = self.root / lvl / f"{pid}.json"
                if not p.exists():
                    fails.append(f"{lvl}/{pid}: listed in store.json, file missing")
                    continue
                docs[(lvl, pid)] = doc = json.loads(p.read_text())
                for n in object_names(doc):
                    obj = self.immutable / n
                    if not obj.exists():
                        fails.append(f"{lvl}/{pid}: object {n} missing")
                    elif sha512t24u_file(obj) != OBJECT_NAME.match(n)["digest"]:
                        fails.append(f"{lvl}/{pid}: object {n} bytes do not match its name")
        for (lvl, pid), doc in docs.items():
            if lvl == CATALOGS:
                for c in doc["chromosomes"]:
                    fails += self._check_header(f"variant_catalogs/{pid} {c['name']}", c["file"], c["name"], c["seq_digest"],
                                                KIND_VARIANTS)
                for key, kind in (("vidx", KIND_VARIANT_INDEX), ("rsid", KIND_RSID)):
                    if key in doc:
                        fails += self._check_header(f"variant_catalogs/{pid} {key}", doc[key], ALL,
                                                    doc.get("collection_digest"), kind)
                fails += self._check_sequences(pid, doc, store)
                fails += self._check_identity(pid, doc, sites)
            if lvl == "annotations":
                # the browser's per-chromosome gene models and gene lookup agree with the tables
                from . import annotation     # late, as for gwas below
                if isinstance(doc.get("lookup"), str):
                    fails += self._check_header(f"annotations/{pid} lookup", doc["lookup"], ALL,
                                                doc.get("identity_digest"), KIND_GENE_LOOKUP)
                fails += annotation.check_split(self, doc)
            if lvl != "experiments":
                continue
            cat = docs.get((CATALOGS, doc["catalog"]))
            if cat is None:
                fails.append(f"experiments/{pid}: variant catalog {doc['catalog']!r} not in store")
                continue
            if ("annotations", doc["annotation"]) not in docs:
                fails.append(f"experiments/{pid}: annotation {doc['annotation']!r} not in store")
            if doc.get("catalog_identity") != cat.get("identity_digest"):
                fails.append(f"experiments/{pid}: catalog_identity {doc.get('catalog_identity')!r} != variant catalog "
                             f"{doc['catalog']!r} identity_digest {cat.get('identity_digest')!r}")
            table = {c["name"]: c["seq_digest"] for c in cat["chromosomes"]}
            files = [(f"{r['phenotype_type']} {chrom}", chrom, n, KIND_RESULTS)
                     for r in doc["results"] for chrom, n in r["files"].items()]
            files += [(f"hits {chrom}", chrom, n, KIND_HITS) for chrom, n in (doc.get("hits") or {}).items()]
            gw = doc.get("gwas") or {}
            files += [(f"gwas {chrom}", chrom, n, KIND_GWAS) for chrom, n in (gw.get("files") or {}).items()]
            for what, chrom, n, kind in files:
                if chrom not in table:
                    fails.append(f"experiments/{pid} {what}: {chrom} not in variant catalog table")
                else:
                    fails += self._check_header(f"experiments/{pid} {what}", n, chrom, table[chrom], kind)
            # objects that span the whole variant catalog: chromosome `all`, the collection digest
            wide = [(f"trans {r['phenotype_type']}", r["trans"]["file"], KIND_TRANS)
                    for r in doc["results"] if r.get("trans")]
            wide += [("gwas index", gw["index"], KIND_GWAS_INDEX)] if gw.get("index") else []
            for what, n, kind in wide:
                fails += self._check_header(f"experiments/{pid} {what}", n, ALL, cat.get("collection_digest"), kind)
            # trans frames tile each object and each gene's frames are one contiguous range
            tfiles = [r["trans"]["file"] for r in doc["results"] if r.get("trans")]
            if tfiles and all((self.immutable / n).exists() for n in tfiles + [doc["search_index"]]):
                from . import results     # late, as for gwas below
                rows = results.load_index(self, doc).to_pylist()
                sizes = {n: (self.immutable / n).stat().st_size for n in tfiles}
                fails += [f"experiments/{pid} {x}" for x in results.check_trans_layout(rows, doc, sizes)]
            # the browser's search index parts and counts agree with the search index
            if doc.get("search_index") and (self.immutable / doc["search_index"]).exists():
                from . import results
                fails += results.check_split(self, doc, cat)
            # the GWAS bin summary (an Arrow object, no header): schema, bins, chromosomes with GWAS rows
            bins = gw.get("bins")
            if bins and (self.immutable / bins["file"]).exists():
                from . import gwas     # late: the Arrow stack is not needed to hash bytes
                what = f"experiments/{pid} gwas bins"
                try:
                    t = gwas.decode_bins((self.immutable / bins["file"]).read_bytes())
                except Exception as e:  # noqa: BLE001 -- any decode failure is a validation failure
                    fails.append(f"{what}: {e}")
                else:
                    fails += [f"{what}: {x}" for x in gwas.check_bins(t, cat, bins.get("bin_bp", 0), bins.get("n_bins", -1))]
                    extra = sorted(set(t.column("chr").to_pylist()) - set(gw.get("files") or {}))
                    if extra:
                        fails.append(f"{what}: chromosomes {extra} have bins and no GWAS file")
            # a hits file's frame table covers exactly the chromosome's variants (u32 at byte 24)
            counts = {c["name"]: c["count"] for c in cat["chromosomes"]}
            for chrom, n in (doc.get("hits") or {}).items():
                obj = self.immutable / n
                if chrom in counts and obj.exists():
                    with open(obj, "rb") as f:
                        h = parse_file_header(f.read(HEADER_LEN))
                    if h["n_cis"] != counts[chrom] or h["page_size"] < 1:
                        fails.append(f"experiments/{pid} hits {chrom}: header covers {h['n_cis']} variants in frames "
                                     f"of {h['page_size']}; the variant catalog has {counts[chrom]}")
        return fails

    def _check_header(self, what: str, name: str, chrom: str, seq_digest: str, kind: int) -> list[str]:
        obj = self.immutable / name
        if not obj.exists():
            return []          # already reported as missing
        with open(obj, "rb") as f:
            try:
                h = parse_file_header(f.read(HEADER_LEN))
            except ValueError as e:
                return [f"{what}: {e}"]
        out = []
        if h["kind"] != kind:
            out.append(f"{what}: header kind {h['kind']} != {kind}")
        if h["chrom"] != chrom:
            out.append(f"{what}: header chromosome {h['chrom']!r} != {chrom!r}")
        if h["seq_digest"] != seq_digest:
            out.append(f"{what}: header seq_digest {h['seq_digest']} != variant catalog {seq_digest}")
        return out

    @staticmethod
    def _refgetstore(refget):
        """An open refgetstore passed straight through, or one opened from a Config. Opening goes through
        `steps_refget.open_store` so a validation reads the store the build wrote, by the same rules
        (remote when `reference.store_url` is set, local otherwise). Imported late: qtlstore is the format
        core and should not drag in the pipeline's config, DuckDB and Arrow stack to hash some bytes."""
        if refget is None or hasattr(refget, "get_sequence_metadata"):
            return refget
        from .steps_refget import open_store
        return open_store(refget)

    def _check_sequences(self, pid: str, doc: dict, store) -> list[str]:
        """Every chromosome of a variant catalog names a sequence the refgetstore really holds, at the length the
        variant catalog records. This is the check that reaches outside the store: without it a `seq_digest` only
        has to match the other copies of itself, so a variant catalog built on a wrong or invented digest is
        internally consistent and wrong."""
        if store is None:
            self.notes.append(f"variant_catalogs/{pid}: sequences not checked, no refgetstore configured")
            return []
        out = []
        for c in doc["chromosomes"]:
            meta = store.get_sequence_metadata(c["seq_digest"])
            if meta is None:
                out.append(f"variant_catalogs/{pid} {c['name']}: seq_digest {c['seq_digest']} is not in the refgetstore")
            elif meta.length != c.get("length"):
                out.append(f"variant_catalogs/{pid} {c['name']}: length {c.get('length')} != refgetstore {meta.length} "
                           f"for {c['seq_digest']}")
        return out

    def _check_identity(self, pid: str, doc: dict, sites) -> list[str]:
        """The variant catalog identity digest recomputed from the sites the variant catalog actually holds.

        Recomputing it means decoding every variant page of every chromosome -- tens of millions of sites
        for a real variant catalog, far more work than the rest of validate -- and the v1 page decoder lives in the
        variant catalog builder, not here. So the reader is injected rather than assumed, and the check runs
        whenever one is available; there is no "skip identity" flag, because the failure mode of a flag is
        that it defaults to on and nobody notices. With no reader the digest is still checked for shape and
        the skip is recorded, so an identity nobody recomputed is never reported as one that passed."""
        ident = doc.get("identity_digest")
        if not isinstance(ident, str) or not DIGEST.match(ident):
            return [f"variant_catalogs/{pid}: identity_digest {ident!r} is not a 32-character sha512t24u"]
        if sites is None:
            self.notes.append(f"variant_catalogs/{pid}: identity digest not recomputed, no site reader given")
            return []
        # streamed one chromosome at a time: a variant catalog's sites do not have to fit in memory at once
        got = catalog_identity([(c["seq_digest"], c) for c in doc["chromosomes"]], load=lambda c: sites(self, c))
        if got != ident:
            return [f"variant_catalogs/{pid}: identity digest recomputes to {got}, the variant catalog says {ident}"]
        return []


def object_names(doc) -> list[str]:
    """Every string in a pointer document that names an immutable object."""
    out = []
    if isinstance(doc, dict):
        for v in doc.values():
            out += object_names(v)
    elif isinstance(doc, list):
        for v in doc:
            out += object_names(v)
    elif isinstance(doc, str) and OBJECT_NAME.match(doc):
        out.append(doc)
    return out


def _write_json(path: Path, doc: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(doc, indent=1) + "\n")
    os.replace(tmp, path)
    return path


# ---- experiment modularity: remove, gc -----------------------------------------------------------
def remove_experiment(store: Store, exp_id: str) -> dict:
    """Drop `experiments/<exp_id>.json` and rewrite `store.json` (same name and refget URLs). Only the
    pointer goes: its objects stay in `immutable/` until `gc`, and its variant catalog and annotation
    pointers stay, since another experiment may name them. Returns what was removed."""
    p = store.root / "experiments" / f"{exp_id}.json"
    if not p.exists():
        raise ValueError(f"experiments/{exp_id}.json does not exist")
    st = json.loads((store.root / "store.json").read_text())
    doc = json.loads(p.read_text())
    p.unlink()
    store.write_store(st["name"], st.get("refget", []))
    return {"experiment": exp_id, "objects_now_unreferenced_candidates": len(object_names(doc))}


def referenced(store: Store) -> set[str]:
    """Every object any pointer file on disk names (listed in store.json or not: a pointer present but
    unlisted still protects its objects)."""
    out = set()
    for lvl in POINTER_DIRS:
        for p in sorted((store.root / lvl).glob("*.json")):
            out |= set(object_names(json.loads(p.read_text())))
    return out


def gc(store: Store, dry_run: bool = False) -> dict:
    """Delete every object in `immutable/` that no pointer names, and leftover `*.tmp` files from an
    interrupted `put`. Content addressing makes this safe to repeat: a kept object's name is its digest,
    so nothing that remains is touched. Returns counts and bytes."""
    keep = referenced(store)
    gone, nbytes = [], 0
    for f in sorted(store.immutable.iterdir()) if store.immutable.exists() else []:
        if f.name in keep:
            continue
        if not (OBJECT_NAME.match(f.name) or f.name.endswith(".tmp")):
            continue                              # not ours: leave it
        gone.append(f.name)
        nbytes += f.stat().st_size
        if not dry_run:
            f.unlink()
    return {"deleted" if not dry_run else "would_delete": len(gone), "bytes": nbytes, "kept": len(keep),
            "examples": gone[:5]}


# ---- cross-catalog lookup check ---------------------------------------------------------------------
def crosscat(store: Store, a: str, b: str, per_chrom: int = 500, rs_per_chrom: int = 200, seed: int = 0) -> dict:
    """SPEC.md section 11 on two variant catalogs of one store: per chromosome, the sites both hold (same
    seq_digest, pos, ref, alt); sampled shared sites of `a` found in `b` by site key and decoded back to the
    same site; sampled sites with an rsID looked up in `b`'s rsID index landing on the same site. `pass` is
    every key lookup exact and at least 95% of rsID lookups found (an rsID can name a different allele pair
    in the other study)."""
    from .catalog import decode_file, decode_rsid
    da, db = store.load(CATALOGS, a), store.load(CATALOGS, b)
    out = {"a": a, "b": b, "same_collection": da["collection_digest"] == db["collection_digest"], "chroms": {}}
    rb = decode_rsid((store.immutable / db["rsid"]).read_bytes())
    tb = {c["seq_digest"]: (i + 1, c) for i, c in enumerate(db["chromosomes"])}
    rng = np.random.default_rng(seed)
    tot = dict(a=0, b=0, shared=0, look=0, found=0, rs_look=0, rs_found=0, rs_same_site=0)
    for ca in da["chromosomes"]:
        if ca["seq_digest"] not in tb:
            out["chroms"][ca["name"]] = "no chromosome with this seq_digest in " + b
            continue
        ordb, cb = tb[ca["seq_digest"]]
        x = decode_file((store.immutable / ca["file"]).read_bytes())
        y = decode_file((store.immutable / cb["file"]).read_bytes())
        kb = {k: i for i, k in enumerate(zip(y["pos"].tolist(), y["ref"], y["alt"]))}
        ka = list(zip(x["pos"].tolist(), x["ref"], x["alt"]))
        shared = [i for i, k in enumerate(ka) if k in kb]
        tot["a"] += len(ka)
        tot["b"] += len(kb)
        tot["shared"] += len(shared)
        for i in (rng.choice(shared, min(per_chrom, len(shared)), replace=False) if shared else []):
            j = kb[ka[i]]
            tot["look"] += 1
            tot["found"] += (int(y["pos"][j]), y["ref"][j], y["alt"][j]) == ka[i]
        cand = [i for i in (rng.choice(shared, min(rs_per_chrom, len(shared)), replace=False) if shared else [])
                if x["rs_number"][i] > 0]
        for i in cand:
            rs = int(x["rs_number"][i])
            tot["rs_look"] += 1
            lo, hi = np.searchsorted(rb["rs_number"], rs, "left"), np.searchsorted(rb["rs_number"], rs, "right")
            hits = [int(v) for v, o in zip(rb["vidx"][lo:hi], rb["ordinal"][lo:hi]) if o == ordb]
            tot["rs_found"] += bool(hits)
            tot["rs_same_site"] += any((int(y["pos"][v]), y["ref"][v], y["alt"][v]) == ka[i] for v in hits)
        out["chroms"][ca["name"]] = {"a": len(ka), "b": len(kb), "shared": len(shared)}
    out["total"] = tot
    out["pass"] = bool(out["same_collection"] and tot["look"] > 0 and tot["found"] == tot["look"]
                       and tot["rs_found"] >= 0.95 * tot["rs_look"])
    return out


def main(argv: list[str] | None = None) -> int:
    """`python -m pipeline.qtlstore <command> --store DIR ...`:

        validate                   Store.validate with the configured refgetstore and the site reader
        remove-experiment ID       drop experiments/ID.json and rewrite store.json (objects stay until gc)
        gc [--dry-run]             delete immutable objects no pointer names
        crosscat A B               the cross-catalog lookup check between variant catalogs A and B
    """
    import argparse
    ap = argparse.ArgumentParser(description="qtlstore maintenance: validate, remove-experiment, gc, crosscat")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("validate", "remove-experiment", "gc", "crosscat"):
        p = sub.add_parser(name)
        p.add_argument("--store", required=True, type=Path)
        if name == "remove-experiment":
            p.add_argument("id")
        if name == "gc":
            p.add_argument("--dry-run", action="store_true")
        if name == "crosscat":
            p.add_argument("a")
            p.add_argument("b")
    a = ap.parse_args(argv)
    st = Store(a.store)
    if a.cmd == "validate":
        from .catalog import read_sites
        from .common import Config
        fails = st.validate(refget=Config(), sites=read_sites)
        print("validate:", "PASS" if not fails else f"{len(fails)} FAILURES", "notes:", st.notes)
        for f in fails[:50]:
            print("  ", f)
        return 1 if fails else 0
    if a.cmd == "remove-experiment":
        print(json.dumps(remove_experiment(st, a.id)))
        return 0
    if a.cmd == "gc":
        print(json.dumps(gc(st, a.dry_run)))
        return 0
    r = crosscat(st, a.a, a.b)
    for name, c in r["chroms"].items():
        print(f"{name}: {c}")
    print("TOTAL", {k: f"{v:,}" for k, v in r["total"].items()})
    print("cross-catalog lookup:", "PASS" if r["pass"] else "FAIL")
    return 0 if r["pass"] else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
