"""Tests for the qtlstore core (pipeline/qtlstore.py).

    uv run python -m pipeline.test_qtlstore

Plain asserts, no test framework. Builds a two-experiment store on one variant catalog in a temp directory.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

from . import qtlstore as qs

CASES = []


def case(fn):
    CASES.append(fn)
    return fn


def raises(rule: str, fn, *args, **kwargs) -> None:
    try:
        fn(*args, **kwargs)
    except ValueError as e:
        assert rule in str(e), f"expected {rule!r} in error, got {e!r}"
        return
    raise AssertionError(f"no ValueError raised (expected {rule!r})")


CHROMS = ("chr1", "chr2")
SEQS = {"chr1": "ACGT" * 100, "chr2": "TTGCA" * 90}
SEQ = {c: qs.sha512t24u(s.encode()) for c, s in SEQS.items()}
LEN = {c: len(s) for c, s in SEQS.items()}
# the sites a variant catalog's identity digest is computed over, one entry per chromosome
SITES = {"chr1": ([1, 2, 3], "ACG", "GTA"), "chr2": ([7, 9], "TT", "AC")}
COLLECTION = qs.sha512t24u(b"test collection")        # stands in for the seqcol digest in the index headers


def refgetstore(path: Path):
    """A real (tiny) refgetstore holding the two test sequences, built the way test_refcheck builds
    one. The unit tests carry their own genome so they run anywhere."""
    from gtars.refget import RefgetStore
    path.parent.mkdir(parents=True, exist_ok=True)
    fa = path.with_name(path.name + ".fa")
    fa.write_text("".join(f">{c} test\n{s}\n" for c, s in SEQS.items()))
    store = RefgetStore.on_disk(str(path))
    store.set_quiet(True)
    store.add_sequence_collection_from_fasta(str(fa))
    return store


class _Cfg:
    """Just enough Config for `steps_refget.open_store`: a local store path, no store_url."""

    def __init__(self, store: Path):
        self._d = {"reference": {"store": str(store)}}

    def __getitem__(self, k):
        return self._d[k]


def pack(chrom: str, seq_digest: str) -> bytes:
    """A stand-in variant pack: the v1 header, then the chromosome's sites as JSON. The v1 page codec
    is not written yet (design plan steps 5 and 8), so this plus `read_sites` plays the part of the
    decoder `validate` will be handed by the variant catalog builder."""
    pos, ref, alt = SITES[chrom]
    return qs.file_header(1, chrom, len(pos), 512, len(pos), seq_digest) + json.dumps([pos, ref, alt]).encode()


def read_sites(store: qs.Store, chrom: dict):
    """`sites(store, chrom) -> (pos, ref, alt)`: the sites really in this chromosome's pack object."""
    pos, ref, alt = json.loads((store.immutable / chrom["file"]).read_bytes()[qs.HEADER_LEN:])
    return pos, ref, alt


@case
def digest_matches_refget():
    # refget's published test vector for sha512t24u("ACGT")
    assert qs.sha512t24u(b"ACGT") == "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
    assert qs.sha512t24u(b"") == "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc"


@case
def header_round_trip():
    b = qs.file_header(1, "chrX", 1234, 512, 1200, SEQ["chr1"])
    assert len(b) == 64
    h = qs.parse_file_header(b)
    assert h == {"kind": 1, "version": 1, "chrom": "chrX", "count": 1234, "page_size": 512, "n_cis": 1200,
                 "seq_digest": SEQ["chr1"]}
    raises("sha512t24u", qs.file_header, 1, "chr1", 1, 1, 0, "abc")
    raises("1-8 ASCII", qs.file_header, 1, "chromosome1", 1, 1, 0, SEQ["chr1"])
    raises("magic", qs.parse_file_header, b"XXXX" + b[4:])
    raises("version", qs.parse_file_header, b[:5] + b"\x00" + b[6:])
    raises("reserved", qs.parse_file_header, b[:60] + b"\x01\x00\x00\x00")
    raises("need 64", qs.parse_file_header, b[:32])


@case
def identity_is_canonical():
    """The identity is the set of sites: chromosomes by seq_digest, sites by (pos, ref, alt), whatever
    order they are handed in."""
    pos, ref, alt = [100, 200, 200, 300], ["A", "C", "C", "G"], ["G", "T", "TA", "GA"]
    a = qs.catalog_identity([(SEQ["chr1"], pos, ref, alt), (SEQ["chr2"], [5], ["T"], ["C"])])
    shuffled = [3, 1, 0, 2]
    b = qs.catalog_identity([(SEQ["chr2"], [5], ["T"], ["C"]),
                             (SEQ["chr1"], *([x[i] for i in shuffled] for x in (pos, ref, alt)))])
    assert a == b                                     # chromosome order and site order stay out
    lines = sorted([f"{SEQ['chr1']}\t{p}\t{r}\t{t}\n" for p, r, t in zip(pos, ref, alt)] + [f"{SEQ['chr2']}\t5\tT\tC\n"],
                   key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]), *x.split("\t")[2:]))
    assert a == qs.sha512t24u("".join(lines).encode())
    # positions sort as numbers, alleles byte-wise
    assert qs.catalog_identity([(SEQ["chr1"], [9, 10], ["a", "a"], ["b", "b"])]) == \
        qs.sha512t24u(f"{SEQ['chr1']}\t9\ta\tb\n{SEQ['chr1']}\t10\ta\tb\n".encode())
    assert a != qs.catalog_identity([(SEQ["chr2"], pos, ref, alt), (SEQ["chr1"], [5], ["T"], ["C"])])
    assert a != qs.catalog_identity([(SEQ["chr1"], pos, ref, ["G", "T", "TA", "A"]), (SEQ["chr2"], [5], ["T"], ["C"])])
    line = f"{SEQ['chr1']}\t100\tA\tG\n".encode()
    assert qs.catalog_identity([(SEQ["chr1"], [100], ["A"], ["G"])]) == qs.sha512t24u(line)
    raises("more than one", qs.catalog_identity, [(SEQ["chr1"], [1], ["A"], ["G"]), (SEQ["chr1"], [2], ["A"], ["G"])])
    # the lazy form reads one chromosome at a time, in canonical order
    seen = []
    got = qs.catalog_identity([(SEQ["chr2"], "b"), (SEQ["chr1"], "a")],
                              load=lambda k: seen.append(k) or {"a": (pos, ref, alt), "b": ([5], ["T"], ["C"])}[k])
    assert got == a and seen == sorted(seen, key=lambda k: SEQ["chr1" if k == "a" else "chr2"])


@case
def orientation():
    """`beta` and `af` describe `effect_allele`; the output is ALT-relative. So a site whose OTHER allele
    is REF passes through untouched, and a site whose EFFECT allele is REF swaps."""
    r = qs.orient_to_ref(ref_base=["A", "C", "G", "T"],
                         effect_allele=["G", "C", "C", "TA"],     # row 0 and 3: already ALT
                         other_allele=["A", "T", "A", "T"],       # row 1: effect allele is REF -> swap
                         beta=[0.5, 0.2, 0.1, -0.3], af=[0.1, 0.3, 0.4, 0.8])
    assert r["counts"] == {"as_is": 2, "swapped": 1, "dropped": 1}
    assert list(r["ref"]) == ["A", "C", "T"] and list(r["alt"]) == ["G", "T", "TA"]
    assert np.allclose(r["beta"], [0.5, -0.2, -0.3]) and np.allclose(r["af"], [0.1, 0.7, 0.8])
    assert list(r["keep"]) == [True, True, False, True]


@case
def orientation_topchef_is_a_relabelling_not_a_sign_flip():
    """The regression guard for a genome-wide silent inversion.

    TOPCHeF's A2 is the reference allele (EVIDENCE.md A.10: A2 reads at all 8,419,594 cis SNPs, A1 at none)
    and `af`/`slope` are A1's (EVIDENCE.md A.11: af correlates +0.9932 with the Jurgens 2024 DCM GWAS EAFREQ
    oriented to A1 over 6,793,566 shared SNPs, -0.9932 against the mirror). A1 is therefore already ALT,
    so the TOPCHeF call must touch nothing: every site as_is, every beta's sign kept, every af kept.

    Under the mirror reading -- `beta`/`af` given for A2 -- every row here would swap instead, negating
    every slope and mirroring every af across the genome. That failure is invisible in a p-value or an SE,
    so it has to be caught here.
    """
    a1 = ["G", "T", "C", "A", "GA"]          # effect allele: `af` and `slope` are A1's
    a2 = ["A", "C", "G", "T", "G"]           # reference allele at every one of these sites
    beta = [0.5, -0.2, 0.1, -0.3, 0.7]
    af = [0.10, 0.25, 0.40, 0.80, 0.55]
    r = qs.orient_to_ref(a2, a1, a2, beta, af)          # the plain positional call must be the right one
    assert r["counts"] == {"as_is": 5, "swapped": 0, "dropped": 0}, r["counts"]
    assert list(r["ref"]) == a2 and list(r["alt"]) == a1
    assert np.allclose(r["beta"], beta), "slopes were negated: the orientation convention is inverted"
    assert np.allclose(r["af"], af), "af was mirrored: the orientation convention is inverted"
    assert list(r["keep"]) == [True] * 5


def build_store(root: Path, seq: dict = SEQ) -> qs.Store:
    """`seq` is the chromosome -> seq_digest table. Overriding it with digests that name no real
    sequence builds the variant catalog that is internally consistent and wrong."""
    st = qs.Store(root)
    chroms = []
    for c in CHROMS:
        var = st.put(pack(c, seq[c]), "qbv")
        chroms.append({"name": c, "seq_digest": seq[c], "length": LEN[c], "count": len(SITES[c][0]), "file": var})
    ident = qs.catalog_identity((seq[c], *SITES[c]) for c in CHROMS)
    st.write_pointer("variant_catalogs", "cat1", {"id": "cat1", "collection_digest": COLLECTION, "identity_digest": ident,
                                          "attributes": ["af", "ma_samples", "rs_number"], "chromosomes": chroms,
                                          "vidx": st.put(qs.file_header(9, "all", 2, 512, 0, COLLECTION), "qbx"),
                                          "rsid": st.put(qs.file_header(8, "all", 0, 4096, 0, COLLECTION), "qbr")})
    st.write_pointer("annotations", "gencode_v39", {"id": "gencode_v39", "genes": st.put(b"genes", "arrow.zst")})
    for eid in ("topchef", "gtex_v8_heart_lv"):
        files = {c: st.put(qs.file_header(2, c, 1, 0, 0, seq[c]) + eid.encode(), "qbe") for c in CHROMS}
        hits = {c: st.put(qs.file_header(7, c, 0, 1024, len(SITES[c][0]), seq[c]) + eid.encode(), "qbh") for c in CHROMS}
        st.write_pointer("experiments", eid, {"id": eid, "catalog": "cat1", "catalog_identity": ident,
                                              "annotation": "gencode_v39", "hits": hits,
                                              "results": [{"phenotype_type": "ge", "dof": None, "files": files}]})
    st.write_store("test", ["https://refget.example/store"])
    return st


@case
def store_builds_and_validates():
    with tempfile.TemporaryDirectory() as d:
        st = build_store(Path(d))
        assert st.validate() == []
        doc = json.loads((Path(d) / "store.json").read_text())
        assert doc["experiments"] == ["gtex_v8_heart_lv", "topchef"] and doc["variant_catalogs"] == ["cat1"]
        vidx = qs.file_header(9, "all", 2, 512, 0, COLLECTION)
        assert st.put(vidx, "qbx") == st.load("variant_catalogs", "cat1")["vidx"]      # idempotent


@case
def pointer_needs_objects_first():
    with tempfile.TemporaryDirectory() as d:
        st = qs.Store(Path(d))
        raises("write objects first", st.write_pointer, "annotations", "a",
               {"id": "a", "genes": f"{'A' * 32}.arrow.zst"})
        raises("does not match", st.write_pointer, "annotations", "a", {"id": "b"})


@case
def validate_catches_damage():
    with tempfile.TemporaryDirectory() as d:
        st = build_store(Path(d))
        exp = st.load("experiments", "topchef")
        # corrupt one object's bytes
        (st.immutable / exp["results"][0]["files"]["chr1"]).write_bytes(qs.file_header(2, "chr1", 1, 0, 0, SEQ["chr1"]))
        # point a results file at the wrong genome
        wrong = st.put(qs.file_header(2, "chr2", 1, 0, 0, SEQ["chr1"]), "qbe")
        exp["results"][0]["files"]["chr2"] = wrong
        exp["results"][0]["files"]["chr3"] = wrong
        (Path(d) / "experiments" / "topchef.json").write_text(json.dumps(exp))
        fails = st.validate()
        assert any("do not match its name" in f for f in fails), fails
        assert any("header seq_digest" in f and "chr2" in f for f in fails), fails
        assert any("chr3 not in variant catalog table" in f for f in fails), fails


@case
def validate_checks_hits_indexes_and_catalog_identity():
    """Hits headers against the variant catalog table, the variant and rsID index headers against the collection
    digest, and each experiment's `catalog_identity` against its variant catalog."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st = build_store(d)
        assert st.validate() == []
        exp = st.load("experiments", "topchef")
        exp["hits"]["chr2"] = st.put(qs.file_header(7, "chr2", 0, 1024, 5, SEQ["chr1"]), "qbh")    # wrong sequence and count
        exp["hits"]["chr1"] = st.put(qs.file_header(2, "chr1", 0, 0, 0, SEQ["chr1"]), "qbh")      # wrong kind
        exp["hits"]["chr9"] = exp["hits"]["chr2"]
        exp["catalog_identity"] = qs.sha512t24u(b"another variant catalog")
        (d / "experiments" / "topchef.json").write_text(json.dumps(exp))
        doc = st.load("variant_catalogs", "cat1")
        doc["vidx"] = st.put(qs.file_header(9, "all", 2, 512, 0, SEQ["chr1"]), "qbx")            # not the collection
        doc["rsid"] = st.put(qs.file_header(8, "chr1", 0, 4096, 0, COLLECTION), "qbr")            # not `all`
        (d / "variant_catalogs" / "cat1.json").write_text(json.dumps(doc))
        fails = st.validate()
        want = ["topchef hits chr2: header seq_digest", "topchef hits chr1: header kind 2 != 7",
                "topchef hits chr2: header covers 5 variants", "topchef hits chr1: header covers 0 variants",
                "topchef hits chr9: chr9 not in variant catalog table", "topchef: catalog_identity",
                "cat1 vidx: header seq_digest", "cat1 rsid: header chromosome 'chr1' != 'all'"]
        for w in want:
            assert any(w in f for f in fails), (w, fails)
        assert len(fails) == len(want), fails
        assert not any("gtex_v8_heart_lv" in f for f in fails), fails


@case
def consistently_wrong_catalog_is_caught():
    """The failure every internal check misses: every chromosome entry *and* every pack header carry
    the same bogus seq_digest, so nothing in the store disagrees with anything else in the store and
    the identity digest, computed over those same digests, recomputes fine. Only the refgetstore can
    say the genome does not exist. Without the refgetstore check this variant catalog validates clean, which
    is what anchoring the format to refget was for."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        bogus = {c: qs.sha512t24u(f"not a sequence: {c}".encode()) for c in CHROMS}
        st = build_store(d / "store", seq=bogus)
        assert st.validate(sites=read_sites) == []            # internally consistent, and wrong
        fails = st.validate(refget=refgetstore(d / "refget"), sites=read_sites)
        assert len(fails) == len(CHROMS) and all("is not in the refgetstore" in f for f in fails), fails
        assert st.notes == []                                  # both checks ran; nothing was skipped


@case
def refgetstore_checks_length_and_opens_from_config():
    """A digest that is in the store but at another length is the subtler version of the same error.
    The store is opened through `steps_refget.open_store`, so a validation and a build resolve the
    configured store the same way."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        refgetstore(d / "refget")
        st = build_store(d / "store")
        cfg = _Cfg(d / "refget")
        assert st.validate(refget=cfg, sites=read_sites) == []
        doc = st.load("variant_catalogs", "cat1")
        doc["chromosomes"][1]["length"] = LEN["chr2"] - 1
        (d / "store" / "variant_catalogs" / "cat1.json").write_text(json.dumps(doc))
        fails = st.validate(refget=cfg, sites=read_sites)
        assert any("!= refgetstore" in f and "chr2" in f for f in fails), fails


@case
def identity_digest_must_recompute():
    """A variant catalog document can claim any identity digest; nothing else in the store contradicts it."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st = build_store(d / "claimed")
        doc = st.load("variant_catalogs", "cat1")
        doc["identity_digest"] = qs.sha512t24u(b"some other variant catalog")
        (d / "claimed" / "variant_catalogs" / "cat1.json").write_text(json.dumps(doc))
        for eid in ("topchef", "gtex_v8_heart_lv"):            # the experiments copy the claim
            exp = st.load("experiments", eid)
            exp["catalog_identity"] = doc["identity_digest"]
            (d / "claimed" / "experiments" / f"{eid}.json").write_text(json.dumps(exp))
        assert st.validate() == []                             # shape is fine; only recomputation tells
        fails = st.validate(sites=read_sites)
        assert any("identity digest recomputes" in f for f in fails), fails
        doc["identity_digest"] = "not a digest"
        (d / "claimed" / "variant_catalogs" / "cat1.json").write_text(json.dumps(doc))
        assert any("not a 32-character" in f for f in st.validate(sites=read_sites)), fails

        # the other direction: the digest stands, the sites move under it
        st = build_store(d / "moved")
        doc = st.load("variant_catalogs", "cat1")
        pos, ref, alt = SITES["chr1"]
        body = json.dumps([[pos[0], pos[1], pos[2] + 500], ref, alt]).encode()
        doc["chromosomes"][0]["file"] = st.put(qs.file_header(1, "chr1", 3, 512, 3, SEQ["chr1"]) + body, "qbv")
        (d / "moved" / "variant_catalogs" / "cat1.json").write_text(json.dumps(doc))
        fails = st.validate(sites=read_sites)
        assert any("identity digest recomputes" in f for f in fails), fails


@case
def validate_degrades_honestly():
    """With no refgetstore and no site reader the two outward checks cannot run. They have to say so:
    an empty failure list here means "nothing disagreed", not "everything was checked"."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        st = build_store(d / "store")
        assert st.validate() == []
        assert any("no refgetstore configured" in n for n in st.notes), st.notes
        assert any("identity digest not recomputed" in n for n in st.notes), st.notes
        assert st.validate(refget=refgetstore(d / "refget"), sites=read_sites) == []
        assert st.notes == [], st.notes                        # nothing left unchecked


@case
def codecs_do_not_mix():
    """The v1 store modules import packfmt_v1 and never packfmt_v0; the v0 reader tools import only packfmt_v0.
    `verify_v0` and `bench_store` compare the two formats and are the only modules allowed both."""
    import ast
    here = Path(__file__).parent

    def codecs(name: str) -> set[str]:
        tree = ast.parse((here / name).read_text())
        out = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                mods = [node.module or ""] + [a.name for a in node.names]
                out |= {m.rsplit(".", 1)[-1] for m in mods if m.rsplit(".", 1)[-1].startswith("packfmt")}
            elif isinstance(node, ast.Import):
                out |= {a.name.rsplit(".", 1)[-1] for a in node.names if "packfmt" in a.name}
        return out
    assert not (here / "packfmt.py").exists(), "pipeline/packfmt.py is gone: use packfmt_v0 or packfmt_v1"
    for m in ("qtlstore.py", "catalog.py", "annotation.py", "results.py", "gwas.py", "test_catalog.py",
              "test_results.py"):
        if (here / m).exists():
            assert codecs(m) <= {"packfmt_v1"}, f"{m} imports {codecs(m)}"
    for m in ("packtool.py", "common.py", "test_packfmt.py", "test_packtool.py", "adapters/verify_topchef.py"):
        assert codecs(m) <= {"packfmt_v0"}, f"{m} imports {codecs(m)}"


def main() -> int:
    bad = 0
    for fn in CASES:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:        # noqa: BLE001 - report every case
            bad += 1
            print(f"FAIL {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(CASES) - bad} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
