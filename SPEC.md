# qtlb format, version 2: the qtlstore

A **qtlstore** is a directory of content-addressed objects holding N experiments over M variant
catalogs, every position anchored to a refget sequence. It is read over plain HTTP range requests.
This document is the byte-level contract for what the builders write today:

| module | writes |
|---|---|
| `pipeline/qtlstore.py` | object names, the 64-byte header, variant catalog identity, orientation, experiment kinds and collections, `store.json`, `Store.validate`, store maintenance (`remove-experiment`, `gc`, `crosscat`) |
| `pipeline/catalog.py` | a variant catalog: per-chromosome variants files, variant index, rsID index |
| `pipeline/annotation.py` | a gene/exon annotation object, its per-chromosome genes and exon models, and its gene lookup |
| `pipeline/results.py` | the QTL experiments of a cohort, one per scan: results files, hits files, search index and its per-chromosome parts, trans objects |
| `pipeline/gwas.py` | a GWAS experiment: per-chromosome GWAS files, the GWAS index and the bin summary |
| `pipeline/dof.py` | the degrees-of-freedom fit a QTL experiment records |
| `pipeline/packfmt_v2.py` | the codec: the quantizers, the results block, the trans frame and the GWAS block |
| `pipeline/packfmt_v1.py` | the codec-free pieces every version shares (zstd framing, variant pages, SNP codes, hits records, the GWAS index payload), and the version 1 layouts for reading a v1 store |
| `pipeline/build_store.py` | one command from contract tables to a validated store |
| `pipeline/upload.py` | the store into the bucket the browsers read, objects first and `store.json` last |
| `pipeline/which_codec.py`, `pipeline/scan_codes.py` | which codec wrote each object, read from the bytes; a sweep of every code in every results block |
| `pipeline/CONTRACT.md` | the input tables every builder reads (the ingestion contract) |

`qtlstore.FORMAT_VERSION` is 2 and is what every builder writes unless told otherwise; `SUPPORTED_VERSIONS`
is `(1, 2)`, because the live bucket holds a version 1 store until the version 2 release. The
builders import `packfmt_v1` only for the shared pieces and take the codec as a parameter
(`test_qtlstore.codecs_do_not_mix`); `pipeline/packfmt_v0.py` is the v0 codec and nothing current
imports it.

When this document and the code disagree, one of them is wrong and both are fixed in the same
change. Section 18 lists what is still open.

Measurements behind the design choices are in `plans/2026-10-07-qtlstore-v2.md` (the triple, its
precision, the genome-wide sizes), `plans/2026-10-07-experiment-model.md` (one experiment per scan,
collections, what merging the trans and GWAS layouts cost when it was tried) and, for the pieces
that date from v0 (page size, SNP codes, the reference allele checks), the analysis repo's
`qtlb-format/docs/EVIDENCE.md`.

**Terms.** A **variant catalog** is the list of variants a study's results point to: for each variant
its sequence digest, position, `ref` and `alt`, stored as per-chromosome variants files, plus the
variant index (vidx) and the rsID index. It is not the **eQTL Catalogue**, the EBI resource one of
the experiments comes from; this document spells that one out in full every time. An **experiment**
is one scan: one cis QTL scan, one trans QTL scan, or one GWAS (section 8). A **collection** is the
publication a set of experiments came from, named by a content-addressed document (section 8).

**Status.** Builders and `Store.validate` pass on a genome-wide store built locally
(`data/derived/store-release-v2`: four TOPCHeF experiments, the DCM GWAS, two GTEx v8 heart LV
experiments, the overlap index and two coloc lists; 5.5 GB). Nothing version 2 is uploaded yet. The
browser readers in `ui-topchef/` and `ui/` read this version; `npm run store-check` matches them
against the Python decoders on both collections.

## 1. Conventions

- **Byte order:** little-endian everywhere.
- **Types:** u8, u16, u32, i32, f32, f64. Readers compute in float64.
- **Positions:** u32, 1-based, on the sequence the header's `seq_digest` names.
- **Byte offsets:** u32. Every file must be under 4 GiB; `results.py` raises when a results file
  would pass `0xFFFFFFFF`.
- **Alignment:** variant pages and result blocks start at a multiple of 4 and are zero-padded to
  one. Padding bytes are zero.
- **Strings:** allele heaps are ASCII. JSON is UTF-8, no byte-order mark.
- **zstd frames** (`packfmt_v1.zstd_frame`): exactly one frame per unit, Frame_Content_Size present,
  content checksum on, no dictionary, no trailing bytes. Level 19 in every builder.
  `packfmt_v1.zstd_unframe` enforces all four rules on read.
- **Arrow objects** (`*.arrow.zst`): an Arrow IPC **stream**, uncompressed inside, wrapped in one
  zstd frame (`annotation.encode`). The Arrow JavaScript reader cannot decode IPC buffer
  compression; the browser already has a zstd decoder.

## 2. Store layout

```
<store>/
  store.json                    mutable: name, format version, refget URLs, the ids at each pointer level,
                                and for experiments their kind and collection
  immutable/<digest>.<ext>      every object; <digest> = sha512t24u of the object's bytes. Collection
                                documents live here too, as <digest>.json
  variant_catalogs/<id>.json    mutable pointer: one variant catalog
  overlaps/<id>.json            mutable pointer: one catalog overlap index (section 19), derived
  annotations/<id>.json         mutable pointer: one gene annotation release
  experiments/<id>.json         mutable pointer: one experiment, which is one scan (section 8)
  colocs/<id>.json              mutable pointer: one coloc list for a pair of experiments (section 20)
```

Five mutable levels; everything else is content-addressed. `overlaps/` is derived: it can be
rebuilt from the catalogs, and a store without it is complete. `colocs/` is optional: a list exists
only for a pair someone declared. This mirrors a refgetstore
(`rgstore.json` -> `collections/` -> `sequences/`).

**Write order.** Objects first, then pointers, then `store.json`. `Store.write_pointer` refuses a
pointer that names an object not yet in `immutable/`. `Store.write_store` lists whatever pointer
files exist when it runs. A reader that starts from any level never meets a dangling name. An
upload follows the same order (`pipeline/upload.py`): until `store.json` is replaced the live store
is the old one with every old object still present, and a reader mid-session keeps working.

**`store.json`:**

```json
{"name": "topchef", "format_version": 2, "refget": [],
 "experiments": {
   "topchef_ge_cis":           {"kind": "qtl_cis",   "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json"},
   "topchef_ge_trans":         {"kind": "qtl_trans", "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json"},
   "topchef_leafcutter_cis":   {"kind": "qtl_cis",   "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json"},
   "topchef_leafcutter_trans": {"kind": "qtl_trans", "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json"},
   "dcm_jurgens2024_biobanks": {"kind": "gwas",      "collection": "x9WHpI55YzTxpa95Ngk6dld7aR8pQHdR.json"}},
 "variant_catalogs": ["topchef_grch38"],
 "annotations": ["gencode_v34"],
 "overlaps": [],
 "colocs": ["murray2026_ge_dcm", "murray2026_leafcutter_dcm"]}
```

`experiments` is a **map**, generated by `write_store` from the documents on disk, carrying each
experiment's `kind` and `collection` so a reader can group experiments by publication and route by
kind from the first file it fetches, without opening any pointer; `validate` checks every entry
against its document. The other id lists are sorted. `refget` is a list of refgetstore URLs,
written `[]` today **(future: filled)**. An experiment names its variant catalog in its own key
`catalog` (section 8).

## 3. Object names

`sha512t24u(bytes)` = base64url (alphabet with `-` and `_`; 24 bytes encode to 32 characters with
no padding) of the first 24 bytes of SHA-512: 32 characters from `[A-Za-z0-9_-]`. It is the refget
function (`gtars.refget.sha512t24u_digest`); `test_qtlstore.py` checks the two agree.

```
sha512t24u(b"")     = z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc
sha512t24u(b"ACGT") = aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2
```

An object's name is `<digest>.<ext>`, matching `^[A-Za-z0-9_-]{32}\.[a-z0-9]+(\.[a-z0-9]+)*$`. The
name carries no logical key (no `eqtl.chr1.` stem as in v0): what an object is comes from the
pointer that names it.

| ext | object | kind (header) |
|---|---|---|
| `qbv` | variants file, one per variant catalog chromosome | 1 |
| `qbe` | results file, one per (cis experiment, chromosome) | 2 |
| `qbg` | GWAS file, one per (GWAS experiment, chromosome) | 4 |
| `qgi` | GWAS index, one per GWAS experiment | 5 |
| `qbt` | trans object, one per trans experiment | 6 |
| `qbh` | hits file, one per (QTL experiment, chromosome) | 7 |
| `qbr` | rsID index, one per variant catalog | 8 |
| `qbx` | variant index, one per variant catalog | 9 |
| `qgl` | gene lookup, one per annotation | 10 |
| `qbo` | catalog overlap index, one per store (section 19) | 11 |
| `arrow.zst` | annotation gene and exon tables and their per-chromosome genes and exon models, experiment search index and its per-chromosome parts, GWAS bin summary | no header |
| `json` | a collection document (section 8) | no header |

Kind 3 is unused (v0's detail-less sQTL kind). Kinds 4-6 keep v0's numbers for the same content
(GWAS, GWAS index, trans); kinds were not renumbered across versions, the header's version byte
says which layout a body has.

A pointer document names objects by bare `<digest>.<ext>` strings. `qtlstore.object_names`
collects every string anywhere in a pointer that matches the name pattern; those are the objects
`Store.validate` resolves.

## 4. File header (64 bytes)

Every binary object (`qbv`, `qbe`, `qbg`, `qgi`, `qbt`, `qbh`, `qbr`, `qbx`, `qgl`, `qbo`) starts with this
header. Struct `<4sBBH8sIII32s4x>` (`qtlstore._HEADER`). The header is the same in versions 1 and 2;
only byte 5 differs.

| offset | type | field |
|---|---|---|
| 0 | 4 bytes | magic `QTLB` |
| 4 | u8 | kind, 1..255 (section 3 table) |
| 5 | u8 | version: the layout of the file's **body**, 2 for everything a current build writes, 1 in the live store until the release. A GWAS block carries no magic of its own, so this byte is the only thing in a `.qbg` that names the codec that wrote its rows, and `validate` compares every object's byte to `store.json`'s `format_version` |
| 6 | u16 | header length: 64 |
| 8 | 8 bytes | chromosome name, ASCII, 1-8 bytes, no NUL inside, zero-padded; `all` for objects that span a whole variant catalog |
| 16 | u32 | count (by kind, below) |
| 20 | u32 | page size (by kind, below) |
| 24 | u32 | `n_cis`: cis variants (kind 1); the chromosome's variant count (kind 7, hits); 0 otherwise |
| 28 | 32 bytes | `seq_digest`, ASCII sha512t24u |
| 60 | 4 bytes | reserved, zero |

| kind | count | page size | chromosome | `seq_digest` |
|---|---|---|---|---|
| 1 variants | variants, both sections | variants per page (512) | the chromosome | that sequence's digest |
| 2 results | blocks | 0 | the chromosome | that sequence's digest |
| 4 GWAS | rows | rows per block (2048) | the chromosome | that sequence's digest |
| 5 GWAS index | chromosomes | rows per block (2048) | `all` | the **collection** digest |
| 6 trans | frames (phenotypes with trans rows) | 0 | `all` | the **collection** digest |
| 7 hits | records | variants per frame (1024) | the chromosome | that sequence's digest |
| 8 rsID index | records | records per block (4096) | `all` | the **collection** digest |
| 9 variant index | chromosomes | variants per page (512) | `all` | the **collection** digest |
| 10 gene lookup | buckets (1024) | 0 | `all` | the annotation's `identity_digest` |

`seq_digest` is what lets a file found alone say which exact sequence its positions index. For the
kinds whose chromosome is `all` (5, 6, 8, 9, 11) the same 32 bytes carry the seqcol collection digest
instead, since they span sequences. The gene lookup (kind 10) indexes no sequence, so it carries the
identity digest of the annotation it was built from.

**What a reader checks.** `Store.validate` checks every header against its pointer at build time
(section 14). The browser does not send a separate request for a header: it checks a header only
where one arrives inside bytes it reads anyway (the gene lookup's directory, a hits file's frame
table, a whole-object read). A range read at a pointer-given offset trusts the pointer, which is
what the content-addressed name and `validate` guarantee, and the decoders reject bytes that are
not what the offsets promise (block magic and length, zstd framing, record counts).

**Reader rules** (`parse_file_header`): magic `QTLB`; version in `SUPPORTED_VERSIONS`; header length
64; bytes 60-63 zero; chromosome non-empty, ASCII, no inner NUL; `seq_digest` matches
`^[A-Za-z0-9_-]{32}$`. A reader that decodes bodies must then check the version against the codec it
is about to use: the layouts differ in row width, so the wrong one returns plausible nonsense rather
than failing.

**Example**, kind 1, `chr1`, count 3, page size 512, `n_cis` 2, GRCh38 chr1
(`file_header(1, "chr1", 3, 512, 2, "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO")`):

```
00: 51 54 4c 42 01 02 40 00 63 68 72 31 00 00 00 00   QTLB, kind 1, version 2, hlen 64, "chr1"
10: 03 00 00 00 00 02 00 00 02 00 00 00 59 61 36 52   count 3, page 512, n_cis 2, "Ya6R...
20: 73 37 44 48 68 44 65 67 37 59 61 4f 53 67 31 45
30: 6f 4e 69 33 55 5f 6e 51 39 53 76 4f 00 00 00 00   ...9SvO", reserved
```

## 5. Variant catalog

One variant catalog = one set of sites on one sequence collection. Experiments name a variant catalog; two
experiments may share one. `store.sbatch` today builds one variant catalog per experiment,
`<experiment>_grch38`.

### Site identity and order

A site is `(seq_digest, pos, ref, alt)`. The variant catalog builder (`catalog.build`):

- reads the contract `sites` table (CONTRACT.md): required columns `chr, pos, ref, alt, af,
  ma_samples, in_cis`;
- takes chromosomes in the configured list order, and names each from the refget collection
  (`collection_table`: `{name, seq_digest, length}`); a site on a chromosome outside the list is
  an error;
- rejects a position outside `1..length` and a repeated `(pos, ref, alt)` within a chromosome;
- **sorts** each chromosome itself: the cis section (`in_cis` true) then the trans-only section,
  each by `(pos, ref, alt)` byte-wise. The input order does not matter;
- rejects an `rsid` column that is not exactly `"rs" + rs_number` (null where `rs_number` is -1),
  or an `rsid` without `rs_number`: only `rs_number` is stored, so the text must be rebuildable.

`vidx` is the 0-based row in that order. Cis sites are `0 .. n_cis - 1`, trans-only sites
`n_cis .. count - 1`.

### Variant catalog identity digest

```
identity_digest = sha512t24u( concat, in canonical order, of
                              "<seq_digest>\t<pos>\t<ref>\t<alt>\n" for every site )
```

**Canonical order:** chromosomes by `seq_digest` (ASCII order, as seqcol's `sorted_sequences`),
then sites by `pos` as a number, then `ref`, then `alt`, byte-wise. It is not the file's vidx
order: the identity is a function of the **set** of sites only. The cis/trans-only split, the
chromosome table's order, the chromosome names (the `seq_digest` is in, the name is not),
attribute columns (`af`, counts, rsIDs) and every encoding choice (page size, codec, zstd level)
stay out, so any two builds of the same sites on the same sequences agree. ASCII, `pos` in decimal
with no padding (`qtlstore.catalog_identity`, which sorts; it rejects two entries with one
`seq_digest`). The line is exactly the input a VRS allele id needs, so a VRS index can be added
later without changing it.

Example, three chr1 sites (given in any order):

```
Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO\t10177\tA\tAC\n
Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO\t10352\tT\tTA\n
Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO\t11012\tC\tG\n
-> uTMqNy84_Lfiat_iu8W7UfLNxVOSbehb
```

### `variant_catalogs/<id>.json`

```json
{"id": "topchef_grch38",
 "identity_digest": "zWs_ao86wMqtA1CVWyEITvEVKNSPADjA",
 "collection_digest": "EiFob05aCWgVU_B_Ae0cypnQut3cxUP1",
 "orientation": "ref_alt",
 "attributes": ["af", "ma_samples", "ma_count", "rs_number", "match"],
 "page_size": 512,
 "n_sites": 239260,
 "chromosomes": [{"name": "chr21", "seq_digest": "UfFHm4y_HR5oRfvrYbboAunqVBfhFQoy", "length": 46709983,
                  "count": 118103, "n_cis": 108583, "file": "2dFzPjIIJD8vV1x-fIVhcTocFLfwaskE.qbv",
                  "file_digest": "2dFzPjIIJD8vV1x-fIVhcTocFLfwaskE"},
                 {"name": "chr22", "...": "..."}],
 "vidx": "rDOK-l02CIzxTOwxNdZAwij5lTR82yzC.qbx",
 "rsid": "x4cfe5TOBFWGftm8udXbOvJkEJubXJh4.qbr",
 "source": {"sites": "/scratch/ns5bc/qtl-browser/derived-topchef-v2/_tables/topchef/sites.parquet"}}
```

(`variant_catalogs/topchef_grch38.json` of the chr21/chr22 store
`/scratch/ns5bc/qtl-browser/store-c2122-v1e`, built 2026-09-24. `attributes` never lists `rsid`:
stores built before the rsID-text fix, such as `store-two-sqtlfix`, did.)

- `chromosomes` is the **chromosome table**. Its order is the ordinal order: ordinal = 1-based
  position in this list. There is no fixed chr1..chrX list.
- `attributes` lists which of `af, ma_samples, ma_count, rs_number, match` the sites table
  carried, in that order: exactly what the pages store. The rsID text is never listed; a reader
  forms it as `rs<rs_number>` (the builder checked the source's `rsid` agrees).

### Variants file (kind 1, `.qbv`)

`[64-byte header][page][page]...`, pages back to back in vidx order, each starting at a multiple
of 4. The cis section fills pages of P = 512 records (the last may be shorter); the trans-only
section starts a new page. No page mixes sections. A chromosome with no sites is just the header.

**Page header, 12 bytes** (`packfmt_v1._PAGE_HEADER`, `<IIHBB>`, unchanged from v0):

| offset | type | field |
|---|---|---|
| 0 | u32 | stored payload length (the zstd frame length when compressed) |
| 4 | u32 | `first_vidx` |
| 8 | u16 | `n` records, 1 to 65535 |
| 10 | u8 | codec: 0 raw, 1 zstd (the builder always writes 1) |
| 11 | u8 | reserved, 0 |

**Payload** (after decompression), columns back to back; length exactly `4 + 16n + heap_len`:

| offset | type x count | column |
|---|---|---|
| 0 | u32 | `heap_len` |
| 4 | u32 x n | position delta: first entry absolute, then `pos[i] - pos[i-1]` |
| 4 + 4n | u32 x n | `rs_number`, 0 = none |
| 4 + 8n | u16 x n | `af` code: **ALT** frequency, `rint(af * 65534)`, 65535 = null |
| 4 + 10n | u16 x n | `ma_samples`, 65535 = null |
| 4 + 12n | u16 x n | `ma_count`, 65535 = null |
| 4 + 14n | u8 x n | allele code |
| 4 + 15n | u8 x n | flags |
| 4 + 16n | heap_len bytes | allele heap |

Then zero padding of the page to a multiple of 4.

- **Allele codes** 1-12 are SNPs, `(ref, alt)` from `packfmt_v1.SNP_CODES`: AC=1, AG=2, AT=3, CA=4,
  CG=5, CT=6, GA=7, GC=8, GT=9, TA=10, TC=11, TG=12. Code 0: the heap holds `ref + "\t" + alt +
  "\n"` for that record, in record order. A pair in the table always uses its code. Codes 13-255
  are rejected.
- **Flags** (changed from v0): bit 0 `alt_is_minor` (`af < 0.5`; 0 when `af` is null); bits 1-2 the
  rsID match code when the variant catalog declares `match` (0 none, 1 exact, 2 position; 3 is reserved:
  the encoder refuses it and `decode_file` rejects it, as v0 did); bits 3-7 zero.
  v0's bit 2 "alleles not reported" is gone: a contract site always has both alleles.
- The encoder rejects: `af` outside [0, 1], counts above 65534, `rs_number` above 2^32 - 1, a heap
  allele that is empty, non-ASCII, or holds a tab or newline, decreasing positions.
- A raw page with three records, `A>AC` at 10177 (af 0.425, rs367896724, match exact),
  `T>TA` at 10352 (af 0.4375, match exact), `C>G` at 11012 (nothing known):

```
00: 3e 00 00 00 00 00 00 00 03 00 00 00 0a 00 00 00   stored 62, first_vidx 0, n 3, raw | heap_len 10
10: c1 27 00 00 af 00 00 00 94 02 00 00 94 a8 ed 15   deltas 10177, 175, 660 | rs_number...
20: 1e a4 fc 0b cb 2c 73 20 cc 6c ff 6f ff ff 64 00   ... | af 27852, 28671, null | ma_samples...
30: 78 00 ff ff 78 00 8c 00 ff ff 00 00 05 03 03 00   ... | ma_count | codes 0,0,5 | flags 3,3,0
40: 41 09 41 43 0a 54 09 54 41 0a 00 00               "A\tAC\nT\tTA\n", pad
```

### Variant index (kind 9, `.qbx`)

`[64-byte header][one zstd frame]`. Header: chromosome `all`, count = chromosomes, page size = P,
`seq_digest` = collection digest. Payload (`catalog.encode_vidx`), all u32:

```
0    4 bytes  magic "QVX1"
4    u32      n_chrom            (= header count, = chromosome table length)
8    u32      page_size P
12   u32      rsid_block_records B (4096)
16   u32      rsid_n             records in the rsID index
20   u32      rsid_blocks        ceil(rsid_n / B)
24   per chromosome, in chromosome table order:
       u32 n_cis, u32 n_trans, u32 n_pages_cis, u32 n_pages_trans
       u32 page_off[n_pages_cis + n_pages_trans + 1]     byte offset of each page in the .qbv; last = file size
       u32 page_first_position[n_pages_cis + n_pages_trans]
     u32 rsid_first[rsid_blocks]                          first rs_number of each rsID block
```

`n_pages_cis = ceil(n_cis / P)`, `n_pages_trans = ceil(n_trans / P)`. No trailing bytes. Unlike v0
(`QVX0`) there are no reserved fields, no `n_frames` and no `hits_off`, because hits are paged in
their own files (section 8). Chromosome names are not stored: they are the variant catalog JSON's table, in order.

**Page of a vidx:** below `n_cis`, page `vidx // P`; otherwise page `n_pages_cis + (vidx - n_cis) // P`.
Bytes `page_off[k] .. page_off[k+1] - 1`. **Page of a position:** within a section, the last page
whose first position is at or below it; confirm on the decoded records.

### rsID index (kind 8, `.qbr`)

`[64-byte header][records]`, uncompressed. Header: chromosome `all`, count = records, page size = B
(4096), `seq_digest` = collection digest. Records are 12 bytes (`catalog.RSID_DTYPE`):

| offset | type | field |
|---|---|---|
| 0 | u32 | `rs_number` |
| 4 | u32 | `vidx` |
| 8 | u16 | chromosome ordinal (1-based, variant catalog table) |
| 10 | u16 | zero |

Sorted by `(rs_number, ordinal, vidx)`. Only sites with `rs_number > 0` get a record.
**`rs_number` may repeat** (v0 forbade repeats): dbSNP puts several allele pairs under one rsID, and
each is its own site, so a lookup returns a run. Block `b` starts at `64 + 12 * b * B` and holds
`min(B, count - b * B)` records.

**Lookup** (`catalog.rsid_lookup`): start at the last block whose `rsid_first` is **strictly below**
the number (block 0 when there is none), then read forward block by block while the records are at
or below the number, collecting every record that equals it. The strict `<` matters: a run of one
`rs_number` can cross a block boundary, and then the block whose first record is the number starts
in the middle of the run, so starting there misses the run's first records. (Earlier text said "at
or below", which has that bug; `test_catalog.py::test_rsid_run_across_a_block_boundary` pins the
rule.) A block that ends inside the run means the next block must be read too.

## 6. Annotation

One object pair per GENCODE/Ensembl release, shared by every experiment that names it
(`annotation.build`). Built today: `gencode_v34` (TOPCHeF) and `gencode_v39` (eQTL Catalogue GTEx,
Ensembl 105); `store.sbatch` picks the adapter's `source.gene_annotation` when an experiment names
none, else `gencode_v34`.

### `annotations/<id>.json`

```json
{"id": "gencode_v34", "identity_digest": "<32>",
 "genes": "<digest>.arrow.zst", "exons": "<digest>.arrow.zst",
 "chroms": {"chr1": {"genes": "<digest>.arrow.zst", "exon_models": "<digest>.arrow.zst"}, "...": {}},
 "lookup": "<digest>.qgl",
 "n_genes": 0, "n_transcripts": 0, "n_exons": 0,
 "source": {"file": "gencode.v34.annotation.gtf.gz", "name": "...", "version": "...", "url": "...", "md5": "...", "size": 0}}
```

No timestamp: the same GTF gives a byte-identical pointer. `genes` and `exons` are the canonical
tables; `chroms` and `lookup` are derived from them (below) for readers that need one chromosome or
one gene, and are not part of the identity.

### Tables

- **genes**: `gene_id` (unversioned) string, `version` int32, `name` string, `biotype` string,
  `chr` string, `tss` int32, `strand` string, `start` int32, `end` int32. Sorted by `(chr, start,
  end, gene_id)`. `tss` = `start` on `+`, `end` on `-`. `gene_id` is unique.
- **exons**: `gene_id`, `transcript_id`, `exon_number` int32, `chr`, `start` int32, `end` int32,
  `strand`. Sorted by `(chr, gene_id, start, end, transcript_id, exon_number)`.
- From GTF `gene` and `exon` records only. `_PAR_Y` genes are dropped (they would repeat a
  `gene_id`). A gene id without an integer version suffix is an error. Both quoted and bare GTF
  attribute values are read, so `exon_number` is the real ordinal (v0 stored 0).

### Per-chromosome objects and the gene lookup

A gene page needs one gene: its row, its exon model, its neighbours for the gene track. Reading the
whole-genome tables for that costs 12 MB (the exon table alone is 11 MB for GENCODE v34), so the
builder also writes (`annotation.split`, which `build` and `add-split` both call):

- **`chroms[chr].genes`**: the `genes` rows of that chromosome, same schema and order. Every
  chromosome of `genes` has one.
- **`chroms[chr].exon_models`**: one row per one of those genes, in the same order: `gene_id`
  string, `exon_starts` list<int32>, `exon_ends` list<int32>, the gene's **collapsed exon model**:
  the union of its transcripts' exons as sorted intervals, an exon merged into the previous one when
  its start is at or before that one's end (touching exons, `start = end + 1`, stay apart). Empty
  lists for a gene with no exon records. GENCODE v34: chr1 254 KB, chr7 125 KB, all 2.6 MB. The
  transcript-level table split by chromosome would be 1.0 MB for chr1, so the per-chromosome object
  holds the model rather than the rows. Genes and models are two objects so a reader of gene
  rows alone (a region's gene list) does not download exon models.
- **`lookup`** (kind 10, `.qgl`): gene id and symbol to chromosome, for a reader that is given a gene
  name and does not yet know which chromosome's objects to read.

**Gene lookup layout.** `[64-byte header][u32 offset[B + 1]][bucket 0][bucket 1]...` with B = 1024
(header count). Bucket `b` is bytes `offset[b] .. offset[b + 1] - 1`; `offset[0] = 64 + 4 (B + 1)`,
`offset[B]` is the object size, and an empty bucket has zero bytes. A bucket is an Arrow object
(section 1: an Arrow IPC stream in one zstd frame) with rows `key` string, `gene_id` string, `name`
string, `chr` string, `tss` int32, sorted by `(key, gene_id)`.

- Every gene has a row under the key of its `gene_id`, and one under the key of its `name` when that
  is non-empty and differs. A key may have several rows (a symbol shared by several genes).
- **Key**: the string with ASCII `a`-`z` upper-cased and every other character unchanged
  (`annotation.lookup_key`), so JavaScript and Python agree without Unicode case rules.
- **Bucket**: FNV-1a 32 over the key's UTF-8 bytes (offset basis `0x811C9DC5`, prime `0x01000193`),
  modulo B (`annotation.lookup_bucket`). Pinned vectors, checked in both `test_annotation.py` and
  `npm run store-check`: `"" -> 453`, `"A" -> 716` (FNV-1a `0xC40BF6CC`), `"FLNC" -> 208`,
  `"ENSG00000128591" -> 414`, `"HLA-DRB1" -> 280`, `"Y_RNA" -> 828`, `"Ä" -> 834`.
- **Lookup**: read bytes `0 .. 64 + 4 (B + 1) - 1` once (header and offsets, 4.2 KB), then the
  key's bucket (about 2 KB), and keep the rows whose `key` equals the key.

A reader that needs every gene (a whole gene list, a search box) reads `genes` whole, one request,
rather than every chromosome's object.

### Annotation identity digest

sha512t24u over UTF-8 lines, gene lines then exon lines, each in table order:

```
G\t<gene_id>\t<version>\t<name>\t<biotype>\t<chr>\t<strand>\t<start>\t<end>\t<tss>\n
X\t<gene_id>\t<transcript_id>\t<exon_number>\t<chr>\t<strand>\t<start>\t<end>\n
```

Every stored value enters, `tss` included, so a change in the TSS rule is a new identity. The GTF's
own bytes, the pointer id, and encoding choices do not.

## 7. Orientation

The one convention everything rests on: **`ref` is the base(s) the anchored sequence reads at `pos`,
`alt` is the other allele, and `af` and `beta` describe ALT.** The variant catalog records it as
`"orientation": "ref_alt"`; the experiment says how it was established in
`allele_orientation_source`.

`qtlstore.orient_to_ref(ref_base, effect_allele, other_allele, beta, af)` is the one
implementation. `beta` and `af` must describe `effect_allele`:

| case | result |
|---|---|
| `other_allele == ref_base` | unchanged: `ref = other`, `alt = effect` |
| `effect_allele == ref_base` (and not the case above) | swap: `ref = effect`, `alt = other`, `beta -> -beta`, `af -> 1 - af` |
| neither | dropped, counted |

It returns `counts: {as_is, swapped, dropped}`.

- **TOPCHeF** (`topchef_refcheck_a2_is_ref`): A2 is REF, A1 is ALT, and `af` and `slope` are A1's
  (EVIDENCE.md A.10, A.11). So `A2 -> ref`, `A1 -> alt`, `af` and `beta` unchanged. **A relabelling,
  not a sign flip.** At the 97,433 cis indels the refcheck table calls `a1`, the adapter asks the
  refgetstore whether the reference reads A2 there and takes A2 as `ref` where it does
  (CONTRACT.md, "Corrected"). A rebuild whose slopes come out negated genome-wide is wrong.
- **eQTL Catalogue** (`eqtl_catalogue_ref_alt`): already `ref`/`alt` with ALT the effect allele;
  passed through, checked against the reference.
- A site whose alleles do not match the reference gets no `sites` row; every row naming it is
  left out and counted in `ingestion.json`.

## 8. Experiment

One experiment = **one scan**: one cis QTL scan, one trans QTL scan, or one GWAS, over one variant
catalog and, for the QTL kinds, one annotation (`results.build`, `gwas.build`). A cohort that
published eQTL and sQTL results in cis and trans is four experiments; a GWAS it was colocalized
against is a fifth, its own, because the cohort did not produce it. So TOPCHeF is `topchef_ge_cis`,
`topchef_ge_trans`, `topchef_leafcutter_cis`, `topchef_leafcutter_trans`, and beside them
`dcm_jurgens2024_biobanks`.

`kind` is a **closed** enum, `qtl_cis | qtl_trans | gwas` (`qtlstore.EXPERIMENT_KINDS`), and it
decides the document's shape and the object kinds it names. The rule for what is open and what is
closed: open where the value is a label, closed where the value selects a layout. The phenotype type
(`ge`, `leafcutter`, a protein assay) is a label, so it is open and one code path and one block
layout serve every type; `kind` picks the layout, so a reader that meets one it does not know must
refuse rather than show a partial store as complete. Adding a kind is a version bump.

Ids are `<cohort>_<phenotype_type>_<scope>` for the QTL kinds, derived by the builder from the stem
it is given, and the adapter's own id for a GWAS. They are URL-visible and are what a deployment
names, so they are hard to change later.

### Collections

A collection is the publication a set of experiments came from. It is a **field, not a pointer
level**: each experiment carries `collection`, the object name of a small content-addressed JSON
document, and experiments carrying the same name are in the same collection. Membership is derived
by matching, never declared, so there is no member list to keep in agreement, and the digest
proves that the members share one copy of the citation.

```json
{"title": "Dilated cardiomyopathy GWAS, Jurgens et al. 2024 (biobanks-only meta-analysis)",
 "citation": {},
 "source": {"file": "Jurgens_DCM_GWAS_META_BiobanksOnly.tsv.gz", "config_key": "dcm_gwas",
            "n_cases": 5022, "n_controls": 932941, "rows_with_grch38_position_and_p_above_0": 12504079},
 "parent": null}
```

`Store.write_collection(title, citation, source, parent)` writes it with sorted keys, so identical
metadata is one object. `parent` is another collection's object name, for nesting; a chain is
followed to its root and the depth is capped at 16. A cycle is impossible among well-formed objects,
since each name depends on the other's bytes; the loop guard exists only so the walk terminates on
a store whose bytes no longer match their names, which `validate` reports separately. **Members
need not agree on anything**: a collection records who published what together, not that the
results are comparable, and one publication can release a GRCh37 scan beside a GRCh38 one. A reader
putting two experiments on one axis compares their catalogs itself.

`store.json` carries every experiment's `collection` (section 2), so grouping is a `GROUP BY` over the
first file fetched and a collections page costs one fetch per distinct document, each of which
caches forever.

### `experiments/<id>.json`

Every kind has `id`, `kind`, `collection`, `catalog`, `catalog_identity`, `allele_orientation_source`,
`source` and `precision`; the two QTL kinds add `annotation`, `annotation_source`, `phenotype_type`,
`dof`, `dof_source`, a search index with its parts, and `hits`. What differs is below.

```json
{"id": "topchef_ge_cis", "kind": "qtl_cis", "phenotype_type": "ge",
 "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 "annotation": "gencode_v34", "annotation_source": null,
 "allele_orientation_source": "topchef_refcheck_a2_is_ref",
 "source": {"experiment_id": "topchef", "source": {"zenodo": "zenodo_21382723", "archives": ["..."]}},
 "dof": 435, "dof_fit": null, "dof_source": "ingestion.json",
 "significance": {"column": "p_perm", "op": "<", "threshold": 0.05},
 "files": {"chr1": "<digest>.qbe"},
 "hits": {"chr1": "<digest>.qbh"},
 "search_index": "<digest>.arrow.zst", "n_phenotypes": 19423,
 "search_index_parts": {"chr1": {"file": "<digest>.arrow.zst", "rows": 1959, "ords": [[0, 1958]]}},
 "counts": {"phenotypes": 19423, "with_rows": 19423, "significant": 10220, "significant_genes": 10220},
 "unplaced": {"count": 0, "examples": []},
 "precision": {"theoretical_limit": {"neglog10p": 1.27e-03, "beta_over_se": 1.09e-03, "se_rel": 2.34e-05, "af": 7.63e-06},
               "measured_worst": {"neglog10p": 1.27e-03, "beta_over_se": 1.04e-03, "se_rel": 2.34e-05, "rows_compared": 123458053}}}
```

```json
{"id": "topchef_ge_trans", "kind": "qtl_trans", "phenotype_type": "ge",
 "collection": "T9HbqRhi22izfYL7gVk-NUfsLP2OpRhd.json",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 "annotation": "gencode_v34", "annotation_source": null,
 "allele_orientation_source": "topchef_refcheck_a2_is_ref", "source": {"...": "..."},
 "dof": 435, "dof_source": "ingestion.json",
 "file": "<digest>.qbt", "n_rows": 2642166,
 "hits": {"chr1": "<digest>.qbh"},
 "search_index": "<digest>.arrow.zst", "n_phenotypes": 19440,
 "search_index_parts": {"chr1": {"file": "<digest>.arrow.zst", "rows": 1959, "ords": [[0, 1958]]}},
 "search_index_unplaced": {"file": "<digest>.arrow.zst", "rows": 16, "ords": [[19423, 19438]]},
 "placement": {"rows": 2642166, "rows_in_source": 2642166,
               "rows_skipped_variant_outside_build": 0, "rows_skipped_phenotype_outside_build": 0},
 "excluded": {"variants": 32765, "rows": 37951, "rows_in_source": 2680117,
              "phenotypes_with_excluded_rows": 7200, "reason": "..."},
 "precision": {"theoretical_limit": {"...": "..."}, "measured_worst": {"...": "..."}}}
```

```json
{"id": "dcm_jurgens2024_biobanks", "kind": "gwas",
 "collection": "x9WHpI55YzTxpa95Ngk6dld7aR8pQHdR.json",
 "title": "Dilated cardiomyopathy GWAS, Jurgens et al. 2024 (biobanks-only meta-analysis)",
 "catalog": "topchef_grch38", "catalog_identity": "<32>",
 "allele_orientation_source": "topchef_refcheck_a2_is_ref",
 "files": {"chr1": "<digest>.qbg"}, "index": "<digest>.qgi",
 "n_rows": 12494731, "rows_by_chrom": {"chr1": 985110}, "rows_sharing_a_site": 1593062,
 "block_rows": 2048, "n_values": [42637, "..."],
 "bins": {"file": "<digest>.arrow.zst", "bin_bp": 5000000, "n_bins": 569},
 "precision": {"theoretical_limit": {"...": "..."}, "measured_worst": {"...": "..."}},
 "orientation": {"as_is": 10401887, "swapped": 2092846, "dropped": 9346, "...": "..."},
 "source": {"file": "...", "n_cases": 5022, "n_controls": 932941, "...": "..."}}
```

- `catalog_identity` copies the variant catalog's `identity_digest` at build time.
- `annotation_source` is the source's own annotation when it differs from the one attached
  (a known gap, reported), else null. A GWAS has no annotation: its one phenotype is a trait.
- `significance` is a cis scan's rule, from the adapter's `ingestion.json`; absent there means
  `p_perm < 0.05`. `column` names **any numeric column of `permuted`**, `op` is one of `<` `<=` `>`
  `>=`, `threshold` is a number, and `label` is what a reader shows the metric as ("Permutation p",
  "Bonferroni p"). The builder refuses a rule whose column `permuted` does not carry. `p_perm` and
  `p_beta` are the conventional names for a tensorQTL-family permutation pass, not the only possible
  metrics: a study whose source ran PLINK2 `--glm` or MatrixEQTL has neither, and may declare a rule
  on whatever it does publish or compute. **`significance` may be `null`**: the source assessed no
  significance. Then every `significant` below is null too, which is not false. "Never assessed"
  and "tested and not significant" are different claims, and a reader that prints the second for
  the first is lying for the study. Every `significant` below means "the group's `column` passes
  the rule", and `sig_value` is the value that was tested, whichever column it came from. A trans
  scan carries no rule: the permutation pass is a cis quantity.
- `unplaced`, on a cis scan, counts its phenotypes that no nominal row, permuted group or credible
  set places on a chromosome.
- `files` has one entry per chromosome built, including chromosomes with no blocks (a header-only
  file). A trans scan has one object, `file`, with chromosome `all`.
- `placement` and `excluded`, on a trans scan, are that type's own: its rows written, in the source,
  and skipped by a subset build; and what its adapter left out for want of source alleles
  (section 12), null when nothing was.
- `precision` has one shape in every kind (section 13).

### Degrees of freedom

`dof` is **provenance**: recorded for a QTL scan, reported, used by nothing that decodes. Every
stored value is decoded from its own code and a null `dof` costs an experiment nothing. From
`ingestion.json` `dof[phenotype_type]`:

- an integer (TOPCHeF publishes 435 for `ge`, 480 for `leafcutter`): stored as is, `dof_fit` null;
- a `dof.fit` result (eQTL Catalogue): stored when `usable`, else **`dof: null`**. `dof_fit` keeps
  `residual_log10p, margin, rows_usable, n_samples, implied_covariates, reason`.

The fit (`dof.fit`): integer grid `n_samples - 200 .. n_samples - 1` (widened when the minimum sits
within 3 of an edge), objective = median `|log10 p - log10(2 * stdtr(dof, -|beta/se|))|` over rows
with finite beta, `se > 0` and `1e-300 < p < 0.5`. It runs on a uniform sample (100,000 rows) and on
the same count from the large-|t| tail. `residual_log10p` is the larger of the uniform residual and
the tail residual at the uniform dof. **usable** = residual at most 0.01 and neither minimum at a
grid edge. `identified` (reported, not gating) = both fits agree and the runner-up is at least 2x
worse on the tail. A source whose p, beta and se are not related by a Student t (SAIGE, REGENIE,
Firth, a score test) fits badly and stores `dof: null`, and loses nothing by it.

### Results file (kind 2, `.qbe`)

`[64-byte header][block][block]...`, one block per phenotype of the cis scan on that chromosome,
back to back, in search-index `ord` order. Header count = blocks, page size 0, `n_cis` 0,
`seq_digest` = the chromosome's.

**Block** (`packfmt_v2.encode_gene_block`): a 72-byte header, then three **columnar** code arrays,
then the credible-set records, then the details frame. The code arrays are raw, not compressed, so
a gene page range-reads one block and decodes it with no inflate step, and a reader that wants only
p can take the first `72 + 2n` bytes.

| offset | type | field |
|---|---|---|
| 0 | 4 bytes | magic `QGB2` |
| 4 | u32 | block length including padding |
| 8 | u32 | `n_rows` |
| 12 | u32 | `var_start` (variant catalog vidx of row 0; `0xFFFFFFFF` when `n_rows` is 0) |
| 16 | u32 | `flags`, reserved 0 (where a per-block widening of the beta code would be signalled) |
| 20 | u32 | `pos_first` (0 when `n_rows` is 0) |
| 24 | u32 | `pos_last` (0 when `n_rows` is 0) |
| 28 | u32 | `n_cs` |
| 32 | f64 | `nlp_max` |
| 40 | f64 | `lse_min` |
| 48 | f64 | `lse_max` |
| 56 | f64 | `beta_max` |
| 64 | u32 | `details_zlen` |
| 68 | u32 | `details_len` |
| 72 | u16 x n | `nlp` codes |
| 72 + 2n | i16 x n | `beta` codes |
| 72 + 4n | u16 x n | `se` codes |
| 72 + 6n | 12 x k | credible-set records `{u32 row, f32 pip, u8 cs_id, 3 zero bytes}`, sorted strictly by `(row, cs_id)` |
| 72 + 6n + 12k | zstd | the details frame |

Block length = `72 + 6 n_rows + 12 n_cs + details_zlen`, rounded up to 4 with zero padding.

- **Rows** are the vidx run `var_start .. var_start + n_rows - 1` covering every nominal row and
  every credible-set site of the phenotype. A vidx in the run the phenotype did not test is a
  **null row**: nlp code 65535, beta code -32768, se code 65535. A run is a contiguous vidx range
  by construction, which does not make it free of null rows: TOPCHeF's FLNC run holds one, a catalog
  site inside the window that the phenotype did not test and that is not a credible-set site
  either. A run never crosses the cis/trans-only boundary; the builder raises if it would.
- **A phenotype with no nominal rows and no credible sets** still gets a block with `n_rows` 0 and a
  search-index row, so a reader can say "no per-variant data published".
- **Every block has a details frame.**
- **Two invariants a partial reader depends on**, stated because a reader relying on them cannot
  detect their loss: `nlp` is the first column, at offset 72, and the code arrays are raw. Reordering
  the columns or framing the row array for size would break every such reader without changing a
  magic or a length, so either is a version bump.
- Codes, scales, the scale rule and the reject rules are section 13's (`decode_gene_block`).

**Details JSON** (`v: 1`), compact, key order as written:

```json
{"v": 1, "phenotype_type": "leafcutter", "phenotype_id": "...", "phenotype_object_id": "...",
 "gene_id": "ENSG...", "has_nominal": true, "n_nominal": 0, "n_credible_sets": 0,
 "extra": {"intron_start": 0, "intron_end": 0, "cluster_id": "clu_1234", "strand": "+"},
 "group": {"lead_phenotype_id": "...", "n_variants": 0, "p_perm": 0.0, "p_beta": 0.0,
           "significant": true, "lead": {"chr": "chr21", "pos": 0, "ref": "A", "alt": "G"}}}
```

Study fields only, never annotation (no symbol, TSS, biotype, bounds). `gene_id` may be null.
`extra` is the contract's `phenotypes.extra`. `group` is the permutation row of the phenotype's
group (`phenotype_object_id`), or null; it is labelled as the group's so a non-lead intron does not
carry the cluster's lead as its own. `significant` applies the experiment's rule to the rule's
`column`, and is null when the experiment declares no rule. `p_perm` and `p_beta` here are the
source's own published numbers under their own names, null when it published neither; the quantity
the rule tested is the search index's `sig_value`.
Readers reject a details object whose `v` is not 1.

### Search index (`.arrow.zst`)

One row per phenotype of the experiment, numbered `ord` from 0 in row order, so **each experiment
owns its ordinal space**: an `ord` in a hits record indexes that experiment's index and no other.
Rows are sorted by (chromosome in variant catalog table order, position of the run's first site or
else the group lead's position, `phenotype_id`), the phenotypes with no chromosome last. A cis
scan's index and a trans scan's carry different columns (`results.CIS_INDEX_SCHEMA`,
`TRANS_INDEX_SCHEMA`), the shared ones first:

| column | in | type | meaning; null when |
|---|---|---|---|
| `ord` | both | uint32 | row number; never null |
| `phenotype_type`, `phenotype_id`, `phenotype_object_id` | both | string | |
| `gene_id` | both | string | join key to the annotation; null when the source names no gene |
| `chr` | both | string | **the phenotype's chromosome**: where its cis results are, else (trans scan) its gene's chromosome from the annotation; null when neither places it |
| `has_nominal` | both | bool | the phenotype has nominal rows |
| `is_group_lead` | both | bool | this phenotype is its group's lead |
| `significant` | both | bool | the group passes the cis scan's rule (on its `column`); false with no group; **null when the scan declares no rule**, which is not the same as false |
| `sig_value` | both | float64 | the value the rule tested, whichever `permuted` column it names; null with no group or no rule |
| `blk_off`, `blk_len` | cis | uint32 | the block in `files[chr]` |
| `var_start`, `n_var` | cis | uint32 | the run; null when `n_rows` is 0 |
| `var_off`, `var_len` | cis | uint32 | variant catalog `.qbv` bytes of the pages covering the run; null likewise |
| `w_lo`, `w_hi` | cis | int32 | positions of the run's first and last sites; null likewise |
| `trans_off`, `trans_len` | trans | uint32 | the phenotype's frame in `file` (section 9; a gene's frames are adjacent) |
| `n_trans` | trans | uint32 | trans rows in that frame |

Membership says which scan a phenotype is in. A phenotype the source reports only in trans
(TOPCHeF: 17 genes, most on chrM, and 36 introns) is a row of the trans scan and absent from the cis
one, where before one shared index carried it with null cis columns. In a trans scan's index `chr`
is filled from the annotation for such a phenotype; it stays null only when the gene is not in the
annotation either (TOPCHeF: 16 `ge` phenotypes, none of `leafcutter`).

No annotation columns, and no schema metadata listing object hashes: content-addressed names
already pin every object. When `phenotypes.has_nominal` exists the builder fails if it disagrees
with the nominal rows it found.

**Parts.** The same rows split for readers that need one chromosome (`results._experiment_index`):
`search_index_parts[chr]` holds the rows whose `chr` is that chromosome, one part per chromosome
built (empty ones included), in **both** kinds, since a gene page asks under the gene's
chromosome whichever scan it reads. A trans scan files the rows with a null `chr` under
`search_index_unplaced`, a part of the same shape (absent when there are none); a cis scan has no
such rows, since every one of its phenotypes has a block on a chromosome. Each part is an Arrow
object with the index's schema and its rows in `ord` order; `rows` is its row count and `ords` its
rows' `ord` values as sorted inclusive `[first, last]` runs, so a reader holding an `ord` (from a
hits record) knows which part to read without reading any. Genome-wide TOPCHeF: the cis parts for
chr1 are 60 KB (`ge`) and 199 KB (`leafcutter`); the trans parts 15 KB and 90 KB at the median.

**`counts`**, on a cis scan, over its phenotypes (`results.type_counts`): `phenotypes`, `with_rows`
(`n_var` set), `significant`, and `significant_genes` (distinct non-null `gene_id` among the
significant). The summary line a landing page prints, so it reads no index.

### Hits file (kind 7, `.qbh`)

What the variant page lists for one variant, keyed by the variant catalog's vidx: for a cis scan
the groups it leads and the credible sets it is in, for a trans scan its trans associations. Each
QTL experiment has its own hits files, one per chromosome, so a variant page reads one per
experiment and joins them; a GWAS experiment has none. Each file is **paged by vidx**: frame `g`
holds the records of vidx `g * F .. (g + 1) * F - 1`, F = 1,024.

```
[64-byte header][u32 frame_off[n_frames + 1]][frame 0][frame 1]...
```

- Header: count = records in the file, page size = F, the u32 at byte 24 = the chromosome's variant
  count V in the variant catalog (both sections), `seq_digest` = the chromosome's.
  `n_frames = ceil(V / F)`; a chromosome with no variants has `n_frames` 0 and a table of one entry.
- `frame_off` (uncompressed, right after the header): absolute byte offsets. `frame_off[0] = 64 + 4
  (n_frames + 1)`, `frame_off[n_frames]` = the file size, never decreasing. Frame `g` is bytes
  `frame_off[g] .. frame_off[g + 1] - 1`. **A frame with no records is zero bytes** (equal offsets):
  a reader knows there is nothing without a request.
- A non-empty frame is one zstd frame that decompresses to a multiple of 20 bytes: the frame's
  records (`packfmt_v1.HIT_DTYPE`), sorted by `(vidx, kind, ord, cs_id)`:

| offset | type | field |
|---|---|---|
| 0 | u32 | `vidx` (inside the frame's range) |
| 4 | u32 | `ord`: search-index row |
| 8 | f32 | `value`: kind 0 the index's `sig_value`, the quantity the experiment's rule tested (NaN when null), kind 1 PIP, kind 2 `f32(-log10 p)` of the source p (+inf for p = 0) |
| 12 | f32 | `beta`: kind 2 the ALT effect; NaN on kinds 0 and 1 |
| 16 | u8 | kind: 0 lead of a group, 1 credible-set member, 2 trans association |
| 17 | u8 | `cs_id` (kind 1), 0 otherwise |
| 18 | u8 | flags: bit 0 significant by the experiment's rule (kind 0 only); other bits zero. One bit, so an experiment that declares no rule reads 0 here -- the search index's null `significant` is where "not assessed" is recorded |
| 19 | u8 | zero |

Kind 0: `ord` is the group's lead phenotype in the cis scan's index. Kind 1: a site in two sets of
one phenotype keeps both records. Kind 2, in a trans scan's files only: one record per trans row
whose variant is this one (section 9), `ord` the phenotype in the trans scan's index. `value` is
`-log10(pvalue)` of the contract `trans` row, computed in float64 from the source p and rounded to
the nearest f32 (+inf for p = 0; p itself would underflow f32); `beta` is the source beta rounded to
f32. These are not the trans frame's values: the frame (section 9) stores the same row's -log10 p
as a 16-bit code, which decodes to within `nlp_max / 131066` of the source's (`nlp_max` of that
phenotype's frame), so a kind 2 `value` and the frame's decoded -log10 p for one row differ by at
most `nlp_max / 131066 + |value| * 2^-24` (the second term is the f32 rounding, at most 1.8e-5 since
a finite -log10 p is at most 300), and both are +inf when p = 0. Likewise `beta` and the frame's
decoded beta differ by at most `beta_max / 65532 + |beta| * 2^-24`. A hits record carries no SE;
the frame does.

The record layout is the same in versions 1 and 2, and so are the variant catalog, the rsID and
variant indexes, the gene lookup and the Arrow objects: these hold no quantized statistics, which is
why `which_codec` reports them codec-independent.

**Variant lookup:** read bytes 0 .. `64 + 4 (n_frames + 1)` once per chromosome (the header and the
table; 2.8 KB for chr1), then frame `vidx // F`, then the records with that `vidx`
(`results.hits_frame`). The frame table sits in the hits file, not in the variant catalog's variant
index, because a variant catalog is shared by experiments and hits belong to one.

## 9. Trans results

The gene page's trans table: every trans association of one phenotype, whichever chromosome its
variants lie on. From the contract `trans` table (CONTRACT.md); `results.build_trans`.

**One object per trans experiment**, kind 6, `.qbt`: `[64-byte header][frame][frame]...`. Header
chromosome `all`, `seq_digest` = the collection digest (rows span every chromosome), count = frames,
page size 0. One frame per phenotype with trans rows, back to back with no padding. A phenotype's
frame is `trans_off .. trans_off + trans_len - 1` from its search-index row. One object rather than
one per chromosome because a phenotype's trans hits lie on a median of 16 chromosomes (TOPCHeF,
measured genome-wide), so chromosome is not a useful axis for them.

**Frame order** (`results.trans_frame_order`): first the phenotypes with a `gene_id`, grouped by gene,
genes in search-index order (a gene's place is the smallest `ord` among its phenotypes, which is its
eQTL phenotype's when it has one), `ord` order within a gene; then the phenotypes with a null
`gene_id`, in `ord` order. So **a gene's frames in each object are one contiguous byte range**: the
gene page reads, per phenotype type, the single range from the smallest `trans_off` to the largest
`trans_off + trans_len` over the gene's phenotypes (`results.gene_trans_ranges`), at most one request
per object. `ord` order alone would not do this: introns sort by position, so one gene's introns
interleave with a neighbouring or overlapping gene's. The search index needs no per-gene columns for
it; `ord` and the search-index row order are unchanged.

**Frame** (`packfmt_v2.encode_trans_frame` / `decode_trans_frame`): one zstd frame; the payload is
exactly `48 + 18 n + 12 k + heap_len` bytes:

| offset | type x count | field |
|---|---|---|
| 0 | 4 bytes | magic `QTT3` |
| 4 | u32 | `n` rows, at least 1 |
| 8 | f64 | `nlp_max`: the frame's largest finite -log10 p |
| 16 | f64 | `beta_max`: the frame's largest \|beta\| |
| 24 | f64 | `lse_min` |
| 32 | f64 | `lse_max` |
| 40 | u32 | `heap_len` |
| 44 | u32 | `n_cs` = k, credible-set records |
| 48 | u32 x n | position delta: absolute on the first row and wherever the ordinal changes, else `pos - previous pos` (0 allowed) |
| 48 + 4n | u32 x n | `rs_number`, 0 = none (the variant catalog's) |
| 48 + 8n | u16 x n | `af` code, ALT, the variant catalog's (`rint(af * 65534)`, 65535 null) |
| 48 + 10n | u16 x n | `nlp` code (section 13) |
| 48 + 12n | i16 x n | `beta` code |
| 48 + 14n | u16 x n | `se` code |
| 48 + 16n | u8 x n | chromosome ordinal of the variant (1-based, the variant catalog's table) |
| 48 + 17n | u8 x n | allele code (section 5's SNP codes, 0 = heap) |
| 48 + 18n | 12 x k | credible-set records, as in the results block |
| 48 + 18n + 12k | heap | `ref + "\t" + alt + "\n"` for each code-0 row, in row order |

Rows are sorted by (ordinal, position, ref, alt), so ordinals never decrease. `beta` is the ALT
effect and `se` is stored, so nothing is derived and a null `dof` costs nothing
(`results.read_trans`). Alleles, rsID and `af` are inline so the table needs no variants file from
another chromosome; the alternative, a `vidx` column with the identity looked up in the catalog,
was the first version 1 layout and cost a gene page up to 33 page reads. A frame starts with an
absolute position, so a builder may later write a large phenotype as several frames inside its one
span, with an index object beside the trans object, and no frame changes; readers inflate a span
frame by frame.

**No `vidx`.** A row's `(ordinal, pos, ref, alt)` names its site, so its vidx is recoverable through
the variant catalog (the variant index and a page, as in cross-catalog lookup, section 11) when a
reader needs one; the browser links a trans row by rsID or `chr:pos` and never needed it. The hits
files' kind 2 records (section 8) keep their own vidx.

**Rules the builder enforces:** every row's `(chr, pos, ref, alt)` is a site of the experiment's
variant catalog and every phenotype is in `phenotypes` (else the build fails); p is in [0, 1] and
beta finite. The same rows, keyed by variant, are the trans scan's hits records (section 8). No
threshold is applied: the object holds every trans association the source published, which for
TOPCHeF is 1,469 B and 95 rows per frame at the median only because its authors filtered first.

**Precision**: the experiment's `precision` block, section 13, folded over every frame.

## 10. GWAS

A GWAS is its own experiment, `kind: gwas` (section 8): a single-phenotype genome-wide scan, read by
position rather than by phenotype, with no annotation, no hits and no search index. From the
contract `gwas` table and `gwas.json` in a tables directory (CONTRACT.md); `gwas.build`. Which QTL
experiments it is shown beside is a deployment's choice, not a fact the store records.

**Values are quantized** like every other statistic in the store: p, beta and se as the section 13
triple against per-block scales, af as `rint(af * 65534)`, n from a table of distinct values. The
adapter orients each row against the refgetstore with `qtlstore.orient_to_ref`: rows are
**ref/alt** like everything else, `ref` the reference base(s), `beta` and `af` ALT-relative (the
source's `EA` is the effect allele; where `EA` is the reference the alleles swap, beta negates and
`af = 1 - EAFREQ`), and rows whose alleles the reference does not read are dropped and counted. For
the DCM GWAS: 12,504,079 source rows, 10,401,887 as given, 2,092,846 swapped, 9,346 dropped. The
printed digits are not recoverable from the store (section 13); measured genome-wide on the DCM
GWAS, beta and se decode to the source's four decimals on all but 2,274 of 12,494,731 rows, and p to
its four significant digits on all but 1.8%.

**A site may repeat.** The DCM GWAS lists 796,531 indels twice, once per allele order, with
different N and statistics (two measurements the meta-analysis did not merge). After orientation
both rows name one `(chr, pos, ref, alt)`; both are kept, and `gwas.rows_sharing_a_site` counts them
(1,593,062). Two more indels are listed twice with mirrored, otherwise identical values; oriented they
are the same row, so the adapter keeps one copy of each (`identical_rows_dropped`: 2). The builder
refuses identical rows. The DCM GWAS object holds 12,494,731 rows.

**GWAS file** (kind 4, `.qbg`, one per chromosome with rows): `[64-byte header][block][block]...`.
Header: the chromosome, count = rows, page size = B (2048, rows per block), `seq_digest` = the
chromosome's. Rows are sorted by (pos, ref, alt, rs_number, p, n, beta). Each block is one zstd
frame of B rows (the last may be shorter), back to back with no padding. Block payload
(`packfmt_v2.encode_gwas_block`), exactly `48 + 18 n + 12 k + heap_len` bytes; a block carries no
magic, so the file header's version byte is what names its layout:

| offset | type x count | column |
|---|---|---|
| 0 | u32 | `n` rows |
| 4 | u32 | `heap_len` |
| 8 | u32 | `n_cs` = k |
| 12 | u32 | reserved 0, so the scales are 8-aligned |
| 16 | f64 | `nlp_max` |
| 24 | f64 | `beta_max` |
| 32 | f64 | `lse_min` |
| 40 | f64 | `lse_max` |
| 48 | u32 x n | position delta: first absolute, then `pos[i] - pos[i-1]` |
| 48 + 4n | u32 x n | `rs_number`, 0 = none |
| 48 + 8n | u16 x n | af code (ALT), `rint(af * 65534)`, 65535 null |
| 48 + 10n | u16 x n | `nlp` code (section 13) |
| 48 + 12n | i16 x n | `beta` code |
| 48 + 14n | u16 x n | `se` code |
| 48 + 16n | u8 x n | `n_code`, index into the index's n table |
| 48 + 17n | u8 x n | allele code of (ref, alt), 0 = heap |
| 48 + 18n | 12 x k | credible-set records, as in the results block: a fine-mapped locus lives in the block a locus read fetches |
| 48 + 18n + 12k | heap | `ref + "\t" + alt + "\n"` per code-0 row |

The row is the trans frame's row with `n_code` where the frame has the chromosome ordinal, since a
GWAS file is one chromosome and N varies per row. The two layouts stay separate file kinds: the
framing that suits a GWAS, blocks that never span a chromosome, shatters a trans phenotype whose
hits sit on 16 chromosomes (`plans/2026-10-07-experiment-model.md`, F2).

**GWAS index** (kind 5, `.qgi`, one per GWAS experiment): `[64-byte header][one zstd frame]`.
Header: chromosome `all`, count = chromosomes with rows, page size = B, `seq_digest` = the collection
digest. Payload (`packfmt_v1.encode_gwas_index_payload`, the same in both versions), all u32 except
names:

```
u32 n_values_count, u32 n_values[n_values_count]      distinct n, ascending
u32 n_chroms                                          = header count
per chromosome, in the variant catalog's table order:
    8 bytes name (ASCII, zero-padded), u32 n_blocks,
    u32 first_position[n_blocks], u32 end_offset[n_blocks]   byte after each block in that chromosome's .qbg
```

Block k spans `end_offset[k - 1] .. end_offset[k] - 1`; block 0 starts at byte **64**.

**Window rule** (`packfmt_v1.gwas_window`, `gwas.read_window`), rows with `lo <= pos <= hi`: `end` =
the last block whose first position is `<= hi` (none: no rows); `start` = the last block whose first
position is `< lo`, else 0; fetch `(start = 0 ? 64 : end_offset[start - 1]) .. end_offset[end] - 1`
in one request, decode blocks `start..end`, keep rows in the window. The strict `<` covers a position
repeated across a block boundary.

### Bin summary

The landing page's genome track draws, per 5 Mb window, the GWAS's strongest p. Reading every GWAS block
for that is too much for a landing page, so the experiment carries a small summary: `bins` names
it, `{"file": "<digest>.arrow.zst", "bin_bp": 5000000, "n_bins": N}`, or `bins` is null when the
GWAS has no rows.

**Source.** Every GWAS gets one, and no adapter supplies it: `gwas.build_bins` computes it from the
same contract `gwas` rows the `.qbg` files hold, already oriented, grouped by `(chr, pos // bin_bp)`.
The lead's beta and allele are therefore alt-relative like every other beta in the store, and its rsID
is the catalog's dbSNP number. Rows the reference does not read were dropped before the builder saw
them, so they are not counted (measured on the DCM GWAS against v0's bins, which counted them: `min_p`
and `n_gws` equal in all 569 bins, `n_variants` lower in 167 by 9,348 rows).

**Object**: one zstd frame around an Arrow IPC stream (section 1), no header, one row per bin that has
rows, sorted by chromosome in the variant catalog's table order, then `bin_start`
(`gwas.BINS_SCHEMA`):

| column | type | meaning |
|---|---|---|
| `chr` | string | chromosome name, in the variant catalog's table |
| `bin_start`, `bin_end` | uint32 | the bin `[bin_start, bin_end)`, `bin_start = k * bin_bp`, `bin_end = bin_start + bin_bp`; 0-based like v0's, so position `pos` falls in bin `pos // bin_bp` |
| `min_p` | float64 | the bin's smallest p, to 3 significant digits (`float(f"{p:.3g}")`) |
| `lead_position` | uint32 | position of the row holding it (inside the bin) |
| `lead_rsid` | string | that row's rsID as `rs<n>` from its `rs_number`, null when it has none |
| `lead_beta` | float64 | that row's `beta` (effect of `alt`), to 3 decimals (Python `round(b, 3)`) |
| `lead_alt` | string | that row's `alt` |
| `n_gws` | uint32 | rows with p < 5e-8 |
| `n_variants` | uint32 | rows in the bin |

The DCM GWAS: 569 bins genome-wide, a 15,650-byte object. **Ties:** when several rows share the smallest
p, the lead is the first of them in the stored row order (section 10's sort), so one row supplies
every lead field. A reader uses `min_p`, `n_gws` and `n_variants`; the lead fields label the bar.

## 11. Cross-catalog lookup

"This site in another experiment": same genome means same `seq_digest`. Take the site's
`(seq_digest, pos, ref, alt)`, find the other variant catalog's chromosome whose `seq_digest` matches, find
the page that can hold `pos` in each section with its variant index, decode, and match `(pos, ref,
alt)`. A hit gives the other variant catalog's vidx; the other experiment's hits file and search index
(`var_start <= vidx < var_start + n_var`) give its results. No merged store-wide variant catalog is built.

Two variant catalogs with the same `identity_digest` hold the same set of sites. Their vidx are
interchangeable only when their chromosome tables and every chromosome's `n_cis` also agree (the
same sites split cis/trans-only differently get the same identity but different vidx). The
builders produce everything this needs. `python -m pipeline.qtlstore crosscat --store S A B` checks it
on two variant catalogs of a store: shared sites per chromosome, sampled lookups by site key decoded
back to the same site, and sampled rsID lookups through the other catalog's rsID index. A VRS-id index
object (VRS id -> variant catalog, chromosome ordinal, vidx) is reserved and not built **(future)**.

## 12. What the store leaves out

- **Trans eQTL variants without source alleles.** TOPCHeF's trans eQTL files report 32,765
  positions as bare `chr:pos`. A position is not a variant, so these get no site and none of their
  trans eQTL rows is ingested; alleles are not inferred from dbSNP, frequency or other studies
  (CONTRACT.md, decided 2026-09-24). `ingestion.json` `trans_excluded` counts, per phenotype type,
  the variants, the rows left out and the phenotypes, and the trans experiment's `excluded`
  repeats its own type's count: 37,951 of 2,680,117 trans eQTL rows, none of the 13,182,408 trans
  sQTL rows. They enter when the authors supply the alleles.
- **GWAS rows the reference does not read.** 9,346 of the DCM GWAS's 12,504,079 rows have alleles
  that match neither reference reading; they are dropped and counted in the experiment's
  `orientation`.

The Python reference decoders are `results.read_block`, `results.read_trans`, `results.hits_frame` /
`decode_hits`, `catalog.decode_file`, `decode_vidx`, `decode_rsid`, `rsid_lookup`, `gwas.read_window`,
`gwas.read_bins` and `annotation.load`.

## 13. Precision

Every statistic is stored as the same **triple**, 6 bytes a row, in every result kind: `nlp`, `beta`
and `se` each as a 16-bit code against scales the scope's header carries, where a scope is one
results block, one trans frame or one GWAS block. Nothing is derived from anything else, so no
source assumption (a Student t relating p to beta/se, a dof, a printed precision) is needed to
decode a value, and a source that breaks one stores exactly what it published.

| code | type | value | reserved |
|---|---|---|---|
| `nlp` | u16 | `code * nlp_max / 65533` | 65534 is p = 0 (underflowed in the source); 65535 null |
| `beta` | i16 | `code * beta_max / 32766`, the sign in the code | -32768 null |
| `se` | u16 | `exp(lse_min + code * (lse_max - lse_min) / 65534)` | 65535 null |

The four `f64` scales per scope are `nlp_max`, the largest finite -log10 p; `beta_max`, the largest
\|beta\|; `lse_min` and `lse_max`, the log of the smallest and largest se. There is no ceiling on
`nlp_max`: resolution degrades gracefully at extreme magnitudes (a step of 0.15 in -log10 p at
`nlp_max` 10,000), where the former limit of 300 failed a build outright.

- **`af`:** `rint(af * 65534)`, error at most 7.63e-6, in every kind.
- **Exact:** positions, alleles, `rs_number`, counts, credible-set rows and `cs_id`; `pip` and hits
  values as f32.
- **Scale rule** (`scan_codes` checks it): when `nlp_max > 0` some row holds 65533, else every finite
  code is 0; when `beta_max > 0` some row holds ±32766; when `lse_min < lse_max` some row holds 0
  and some 65534, when equal every non-null code is 0; when no row has a value the scales are 0.0.
- A null row has all three null codes; the three never disagree.

**The `precision` block**, one shape on every experiment:

```json
"precision": {
  "theoretical_limit": {"neglog10p": 1.269e-03, "beta_over_se": 1.090e-03, "se_rel": 2.345e-05, "af": 7.63e-06},
  "measured_worst":    {"neglog10p": 1.269e-03, "beta_over_se": 1.043e-03, "se_rel": 2.344e-05, "rows_compared": 123458053}}
```

`theoretical_limit` is half a quantization step in each quantity's own unit, computed from the four
header scales and nothing else, so any reader recomputes it from a stored block and it holds for
every row: `neglog10p` = `nlp_max / 131066` absolute; `beta_over_se` = `beta_max / 65532 / se_min`,
the beta step in units of the smallest stored se, since the step is absolute and the row it hurts
most is generally not the row that set `beta_max`; `se_rel` = `expm1((lse_max - lse_min) / 131068)`
relative. `measured_worst` is found by decoding every row after encoding and comparing it with the
source value, with `rows_compared` as the n. Every entry is a maximum. `measured_worst` exceeding
`theoretical_limit` is a detectable defect, since the arithmetic forbids it; `scan_codes` asserts it.
The experiment's block is the fold over its scopes, so the limit is the loosest scope's.

Genome-wide, the TOPCHeF experiments of `data/derived/store-release-v2`, every row compared:

| experiment | rows compared | `neglog10p` limit / measured | `beta_over_se` limit / measured | `se_rel` limit / measured |
|---|---:|---|---|---|
| `topchef_ge_cis` | 123,458,053 | 1.27e-03 / 1.27e-03 | 1.09e-03 / 1.04e-03 | 2.34e-05 / 2.34e-05 |
| `topchef_leafcutter_cis` | 499,596,989 | 1.56e-03 / 1.56e-03 | 1.19e-03 / 1.17e-03 | 2.18e-05 / 2.18e-05 |
| `topchef_ge_trans` | 2,642,166 | 1.05e-03 / 1.04e-03 | 1.53e-03 / 1.28e-03 | 2.27e-05 / 2.25e-05 |
| `topchef_leafcutter_trans` | 13,182,408 | 6.09e-04 / 6.09e-04 | 1.81e-03 / 1.53e-03 | 2.29e-05 / 2.29e-05 |
| `dcm_jurgens2024_biobanks` | 12,494,731 | 2.99e-04 / 2.99e-04 | 3.90e-03 / 3.89e-03 | 3.78e-05 / 3.78e-05 |

The `beta_over_se` figures are in units of the row's own standard error, the denomination of the
4.47e-03 slope budget the format has always stated. The four QTL scans sit at a quarter of it or
less; **the DCM GWAS sits at 87%**. The one assumption still inside the beta code is that
`beta_max / se_min` stays under 293 in every scope, which is where that budget is reached, and a
GWAS block is where it is closest: 2,048 consecutive positions mix a rare large-effect variant,
which sets `beta_max`, with a common small-se one, which sets `se_min`, so the ratio runs far
above any single row's \|z\| (the DCM GWAS's worst block has `beta_max / se_min` of 255 against a
largest \|z\| under 10). A larger or more deeply imputed GWAS would cross it. The remedies, in
order: a smaller block, which narrows the mix; coding the effect as z against a per-block `z_max`,
which makes the bound per row and closes the gap outright (`plans/2026-10-08-v2-compression-audit.md`,
measured size-neutral); or a per-block widening signalled in the results block's `flags`. Never a
global widening, which would spend bytes on every row of every study to serve a few blocks.

**The store is a view, not an archive.** The printed digits of a source are not recoverable from
it; exact values live in the source releases, and this block is what says how far a stored value
can be from them.

## 14. Validation

`Store.validate(refget, sites)` returns failures; `Store.notes` lists checks it could not run. **An
empty failure list with non-empty notes is a partial pass.**

1. `store.json` `format_version` is a supported version, every listed id has its pointer file, and
   **every object header declares that same version** -- a store half-built with one codec and half
   with another fails here rather than decoding at the wrong stride.
2. Every object any pointer names exists and its bytes hash to its name.
3. Every variant catalog chromosome's `.qbv` header parses, and its chromosome and `seq_digest` equal the
   table entry.
4. With a refgetstore: every variant catalog `seq_digest` is a sequence the store holds, at the recorded
   `length`. Without one: a note. (Without this check a variant catalog that carries the same wrong digest
   everywhere passes the rest.)
5. `identity_digest` is a 32-character digest; with a site reader (`catalog.read_sites`) it
   recomputes from the decoded variants files. Without one: a note.
6. Every experiment's `kind` is one of the three, its variant catalog is in the store and so is
   its annotation when it is a QTL kind, and its `catalog_identity` equals the variant catalog's
   `identity_digest`. `store.json`'s `experiments` map agrees with every document on `kind` and
   `collection`, and every `collection` names an object in the store whose `parent` chain resolves
   and terminates. **Not** checked: that experiments sharing a collection agree on anything.
7. Every results file, hits file and GWAS file names a chromosome in its variant catalog table, and
   its header has the right kind (2, 7, 4) and the table entry's chromosome and `seq_digest`.
8. The variant index, rsID index, trans objects and GWAS index have the right kind (9, 8, 6, 5),
   chromosome `all` and the variant catalog's `collection_digest`.
9. Every hits file's header covers exactly its chromosome's variant count (byte 24) in frames of a
   positive size.
10. A GWAS bin summary decodes with exactly `gwas.BINS_SCHEMA` and `n_bins` rows; its chromosomes are
    in the variant catalog table, in table order, and each has a GWAS file; bins ascend without
    repeats, each is `[k * bin_bp, k * bin_bp + bin_bp)` starting inside the chromosome, and holds its
    `lead_position`; `0 < min_p <= 1`; `0 <= n_gws <= n_variants`, `n_variants >= 1`; `lead_beta` and
    `lead_alt` are present (`gwas.check_bins`, which the builder also runs).
11. Every trans object is tiled by its frames as its search index places them (from byte 64 to the
    end, no gap or overlap, one frame per index row), and each gene's frames in it form one
    contiguous range (`results.check_trans_layout`, section 9).
12. Every annotation has `chroms` and `lookup`: one `chroms` entry per chromosome of `genes`, each
    `genes` object equal to that chromosome's rows and each `exon_models` object equal to the models
    recomputed from `exons`; the lookup has kind 10, chromosome `all` and the annotation's
    `identity_digest`, its offsets tile it, and every bucket equals what `lookup_buckets` makes of
    `genes` (`annotation.check_split`).
13. Every QTL experiment has `search_index_parts`, one part per chromosome built, each equal to that
    chromosome's index rows with `rows` and `ords` describing it; a trans scan's
    `search_index_unplaced` equals its rows with no chromosome and is present exactly when there
    are any; a cis scan has no such rows (`results.check_split`).

`python -m pipeline.annotation add-split --store DIR --id ID` and `python -m pipeline.results
add-split --store DIR --id ID` give a store built before 12 and 13 its derived objects, from the
objects it already has, by the same code the builders run; the pointer is rewritten and no other
object changes. `add-split` serves cis scans only.

Every header check also checks the kind: a variants file must be kind 1.

Beyond `validate`: `which_codec` reads the block and frame magics back out of every object and
reports which codec wrote it, since validate compares headers to a declaration and only the bytes
say what encoded a block; `scan_codes` sweeps every code in every results block of a cis scan for
the reserved values, the scale rule and `measured_worst <= theoretical_limit`.

Build-time checks outside `validate`: orphan nominal or credible-set rows (a site not in the
variant catalog) fail the build; a repeated `(phenotype, cs_id, site)` credible-set row fails; one
phenotype's nominal rows on two chromosomes fail; `cs_id` outside 0..127 fails; nominal rows for a
phenotype missing from `phenotypes` fail; a trans row whose site is not in the variant catalog or whose
phenotype is not in `phenotypes`, or with a null p or beta, fails; identical GWAS rows fail; a GWAS
row with a null p, beta, se, af or n fails.

## 15. Store maintenance

`python -m pipeline.qtlstore <command> --store DIR`:

- `validate`: section 14, with the configured refgetstore and the site reader.
- `remove-experiment ID`: deletes `experiments/ID.json` and rewrites `store.json` (same name and
  refget list). Nothing else changes: the experiment's objects stay until `gc`, and its variant
  catalog and annotation pointers stay, since another experiment may name them.
- `gc [--dry-run]`: deletes every file in `immutable/` whose name no pointer file names (a pointer on
  disk protects its objects whether or not `store.json` lists it), and leftover `*.tmp` files from
  an interrupted write. A collection document is an object an experiment names, so it is protected
  the same way. Files that are not object names are left alone.
- `crosscat A B`: section 11's check.

Because objects are named by their bytes, adding or removing an experiment never changes another
experiment's objects: building TOPCHeF alone, adding GTEx, removing GTEx and running `gc` leaves every
TOPCHeF object byte for byte as it was (pipeline/README.md, "Experiment modularity"). Removing one
experiment can no longer delete another publication's data, because a GWAS is its own experiment
and is reachable only through its own pointer.

`python -m pipeline.upload` puts a store in the bucket in section 2's order and prunes what the
store no longer names, as separate commands, each a dry run unless told otherwise.

## 16. Differences from v1 and v0

Version 2 against version 1, the store live in the bucket until the release:

| | v1 | v2 |
|---|---|---|
| what a result row stores | two of `(nlp, beta, se)`: cis `(nlp, se)` with the slope rebuilt through a dof, trans `(nlp, beta)` with se rebuilt, GWAS all three as fixed-point integers | the same 6-byte triple everywhere, nothing derived (section 13) |
| results block | `QGB0`, 64-byte header, interleaved `(nlp, se)` pairs with the slope's sign in the se field | `QGB2`, 72-byte header with `beta_max`, columnar `nlp`, `beta`, `se`, `flags` where `anchor` was |
| trans frame | `QTT2`, 32-byte header, 16-byte row, no se | `QTT3`, 48-byte header with the se scales and `n_cs`, 18-byte row, credible-set section |
| GWAS block | 8-byte header, 21-byte lossless row | 48-byte header with scales and `n_cs`, 18-byte quantized row, credible-set section |
| `dof` | decoded through | provenance |
| `precision` | five keys with inconsistent names and units, cis and trans differing, GWAS none | one block on every experiment, `theoretical_limit` and `measured_worst` in one unit per quantity |
| experiment | one document per cohort: `results[]` per phenotype type, each with a `trans` entry, and a `gwas` field | one document per scan, `kind` closed, a GWAS its own experiment |
| collections | none; `source.zenodo` | a content-addressed document named by every member's `collection` field |
| `store.json` `experiments` | a sorted list of ids | a map of id to `{kind, collection}` |
| search index | one per cohort, both types and both scopes, 21 columns, one ordinal space, a trans-only residue | one per experiment with its own columns and ordinal space; trans parts by the gene's chromosome and an unplaced residue |
| hits | one object per chromosome holding leads, credible sets and trans records together | per experiment: a cis scan's hold leads and credible sets, a trans scan's hold trans records |
| validate | headers against pointers | also every header's version against `store.json`, the experiments map, the collection chains |
| genome-wide TOPCHeF store | 2.94 GB without a GWAS | 4.53 GB with the DCM GWAS; cis blocks +49%, the price of the stored beta |

Version 1 against v0:

| | v0 | v1 |
|---|---|---|
| scope | one study | N experiments over M variant catalogs |
| naming | `<stem>.<16 hex of sha256>.<ext>`, `manifest.json` | `<sha512t24u>.<ext>`, four pointer levels + `store.json` |
| header | 32 bytes, version 0 | 64 bytes, version 1, `seq_digest` added, reserved 4 bytes at 60 |
| genome | only in `manifest.reference` | in every header, plus variant catalog table checked against a refgetstore |
| alleles | `A1`/`A2` as given | `ref`/`alt`, ref = reference base, `af`/`beta` ALT-relative |
| cis order | `(position, A1, A2)` | `(pos, ref, alt)` |
| trans-only section | one site per position; allele-less sites allowed (flags bit 2) | `(pos, ref, alt)`; every site has alleles |
| page flags | bits 0-1 match, bit 2 no alleles | bit 0 `alt_is_minor`, bits 1-2 match |
| chromosome ordinals | fixed chr1..chr22, chrX | position in the variant catalog's own table |
| rsID index | 8-byte records, `(ordinal << 27) \| vidx`, unique rs_number | 12-byte records, u32 vidx + u16 ordinal, rs_number may repeat |
| variant index | `QVX0`, 32-byte payload header, hits offsets | `QVX1`, 24-byte payload header; hits offsets live in each hits file |
| results | kind 2 eQTL with gene details, kind 3 sQTL without details | kind 2 for every phenotype type, details `v: 1` in every block |
| `anchor` | window start per gene | 0; TSS from the annotation |
| annotation | copied into details, `search_index`, several tables | one object per release, joined on `gene_id` |
| `search_index` | one row per gene, u16 `ord`, trans/GWAS pointers, pack-hash metadata | one row per phenotype, u32 `ord`, no annotation, no hash metadata |
| hits | paged frames of 1,024 variants, 6 kinds incl. trans, u16 codes | paged frames of 1,024 variants (offset table in the file), 3 kinds (lead, credible set, trans), 20-byte records, f32 values, u32 `ord` |
| dof | `manifest.packs.dof` | per results set; fitted when not published; may be null |
| significance | fixed in code | `significance` rule in the experiment |
| trans | kind 6, one file per gene chromosome, one frame per gene, no alleles | kind 6, one object per phenotype type, one frame per phenotype, a gene's frames contiguous, alleles inline (section 9) |
| GWAS | kinds 4, 5; `ea`/`nea` as given | kinds 4, 5; ref/alt oriented, `af` ALT, same lossless codes (section 10) |
| GWAS bins | `gwas_dcm_bins.json` at the bucket root, columnar JSON, fixed name, from the source rows | an Arrow object named from `gwas.bins`, computed by the builder for every GWAS from its oriented rows; the same bins and `min_p` (section 10) |
| variant catalog identity | none | `identity_digest` |
| what a first gene page reads | the whole `search_index` (1.9 MB) | one lookup bucket, its chromosome's genes, exon models and index part (sections 6, 8) |

## 17. Where v1 and v0 live

The version 1 specification is this file at commit `be0bb44` (`git show be0bb44:SPEC.md`): the
single-document experiment, the derived slope, the lossless GWAS and the `QGB0`/`QTT2` layouts.
The version 1 store is what the bucket serves until the version 2 release (`store-genome-v1f`,
TOPCHeF and GTEx v8 heart LV), and `packfmt_v1.py` still reads it: `read_block`, `read_trans` and
`read_window` pick the codec from each object's header.

The v0 spec (sections 1-15, byte layouts for GWAS, trans and paged hits included) is in the
analysis repo's git history at `qtlb-format/docs/SPEC.md`, commit `fe5a606`. Its Appendix A is now
`qtlb-format/docs/EVIDENCE.md`. The frozen v0 build is on Rivanna at
`/scratch/ns5bc/qtl-browser/derived/`: the packs in `immutable/`, `manifest.json`, and the v0 tables
in `_tables/`. The TOPCHeF acceptance gate (`pipeline/adapters/verify_topchef.py`) compares the
contract tables against the v0 tables there, through the v0 reader `packtool.py`; nothing else
current reads the packs. The code that built v0 is in this repo's git history at commit `c62bca3`,
and the v0-against-v1 comparison tools at commits before 2026-10-09.

## 18. Open

1. `store.json` `refget` URLs are written empty; the VRS-id index is reserved, not built. The
   overlap index (section 19) matches sites across catalogs *within* one store; VRS ids would match
   them across stores, and are still the thing that is missing.
2. The general app's redesign: collection- and experiment-based navigation with no default entry
   point, gene and variant pages not tied to one experiment, the GWAS chosen rather than pinned
   (`plans/2026-10-09-general-ui-brainstorm.md`). The TOPCHeF app stays pinned to its collection
   and GWAS.
3. The assumptions the format still makes, each with what breaks it and how it would present, are
   `plans/2026-10-07-qtlstore-v2.md`'s gaps G1 to G12: the `beta_max / se_min` ratio (section 13),
   phenotypes in the millions, a per-row sample size on the QTL side (declined), structural
   variants, a second assembly in one store, coloc's standardised-trait assumption, the archive
   question, and a trans source published at a loose threshold.

## 19. Catalog overlap index

One object per store (`overlap.build`), naming a set of variant catalogs and recording, for every
site any of them holds, which of them hold it. It answers three questions the store could otherwise
only answer by merging catalogs by hand: whether a site is in a given catalog, how far two catalogs
overlap, and what a site's vidx is in another catalog. `crosscat` (section 11) samples and checks;
this indexes and counts.

The index is **derived**: it can be rebuilt from the catalogs at any time, and a store without one
is complete. It is not part of any catalog's or experiment's identity.

### `overlaps/<id>.json`

```json
{"id": "grch38_heart", "object": "<digest>.qbo",
 "collection_digest": "EiFob05aCWgVU_B_Ae0cypnQut3cxUP1",
 "catalogs": [{"id": "gtex_v8_heart_lv_grch38", "identity_digest": "<32>", "n_sites": 9613073},
              {"id": "topchef_grch38", "identity_digest": "<32>", "n_sites": 9182261}],
 "mask_width": 1, "chunk": 65536, "union": 0,
 "pairs": [{"a": "gtex_v8_heart_lv_grch38", "b": "topchef_grch38",
            "shared": 0, "jaccard": 0.0, "a_only": 0, "b_only": 0}],
 "by_chrom": {"chr22": {"union": 148485, "shared": {"gtex_v8_heart_lv_grch38|topchef_grch38": 110598}}},
 "conflicts": {"shared_positions_no_shared_allele_pair": 0, "examples": []},
 "normalisation_suspects": {"count": 0, "examples": []}}
```

- `catalogs` is sorted by id and **fixes the bit order**: catalog *i* in this list is bit *i* of
  every mask. Its `identity_digest` is copied at build time, which is what makes a stale index
  detectable: rebuild a catalog and its digest changes, and `validate` fails rather than the index
  silently describing a store that no longer exists.
- Every catalog named must be in the store and carry the same `collection_digest`. An index never
  spans sequence collections: sites on different references are not comparable (section 7).
- `pairs` is every unordered pair, `by_chrom` the same per chromosome. Both are O(k²) in this
  document and in nothing else.
- `conflicts` counts positions where two catalogs both have a site but share no `(ref, alt)` pair.
  These are ordinarily real: different variants at one position. They are reported because a sudden
  rise in them is the signature of an ingestion that got allele orientation wrong.
- `normalisation_suspects` are indel pairs in different catalogs, within a short window, with the
  same length change and no shared allele pair — candidates for one event written two ways
  (section 5 does not require left-aligned indels). Confirming one needs the reference, so the count
  is a prompt to check rather than a claim that two records describe one event.

### Object (kind 11, `.qbo`)

```
[64-byte header][directory][union pages][mask chunks]
```

Header: chromosome `all`, count = catalogs (k, 1..64), page size = union sites per page (512),
`n_cis` 0, `seq_digest` = the collection digest.

**Union order** is section 5's canonical order — chromosomes by `seq_digest` in ASCII order, then
`pos`, then `ref`, then `alt`, byte-wise — the same order `catalog_identity` digests. It is
not any catalog's vidx order, and a reader must not assume the two agree.

**Directory** (`overlap.encode_directory`), all u32 unless stated:

```
0    4 bytes  magic "QBO1"
4    u32      n_union            sites in the union
8    u32      page_size P        (512)
12   u32      chunk              union sites per mask chunk (65536)
16   u32      n_chunks           ceil(n_union / chunk)
20   u8       mask_width         bytes per mask: 1, 2, 4 or 8
21   u8       k                  catalogs (= header count)
22   2 bytes  reserved, zero
24   u32      n_chrom
28   per chromosome, in ASCII seq_digest order:
       32 bytes seq_digest, u32 n_union_chrom, u32 first_union_index
     u32 page_off[n_pages + 1]            byte offset of each union page; last = mask area start
     u32 page_first_position[n_pages]
     u32 chunk_off[n_chunks + 1]          byte offset of each mask chunk; last = file size
     u32 prefix[k][n_chunks]              set bits of catalog i strictly before chunk c
     per catalog i, u32 n_breaks[i], then n_breaks[i] x {u32 union_index, u32 vidx}
```

**Union pages** never span chromosomes: each chromosome starts a new page, so a chromosome with
`n` union sites has `ceil(n / P)` pages and `page_first_position` is a position on that chromosome.
Pages hold the site key only — position delta, allele code, heap — in the layout of a
variants page (section 5) with the attribute columns omitted. No `af`, `ma_samples`, `ma_count`,
`rs_number` or `match`: those describe a cohort, not a site, and they live in the catalogs.

**Mask chunks**: one zstd frame per chunk, decompressing to `mask_width * n` bytes for the chunk's
`n` union sites. Mask bit *i* is set when catalog *i* holds that site. `mask_width` is the smallest
of 1, 2, 4, 8 that holds k bits. Masks repeat heavily — most sites sit in all catalogs or in one
recurring subset — so a chunk usually compresses far below a byte a site.

### Membership, rank and vidx

**Membership** is one chunk read: chunk `u // chunk`, mask at offset `(u % chunk) * mask_width`.

**Rank** — how many of catalog *i*'s sites precede union index `u` — is
`prefix[i][u // chunk]` plus the set bits of *i* in that chunk below `u % chunk`. One table lookup
and a scan bounded by one chunk.

**vidx is not rank.** A catalog numbers each chromosome's cis sites first, sorted by
`(pos, ref, alt)`, and then its trans-only sites, sorted the same way; union order interleaves the
two. So rank gives the site's index among that catalog's sites *in canonical order*, which equals
its vidx only on chromosomes with no trans-only site. The breakpoint table carries the correction:
`{union_index, vidx}` pairs marking the start of each run over which `vidx - rank` is constant.
A reader takes the last breakpoint at or below `u` and returns `vidx_first + (rank(i, u) -
rank(i, union_index_first))`. A catalog whose chromosomes are all cis has one breakpoint per
chromosome; TOPCHeF, with 19 trans-only sites on chr22, has a handful more.

An implementation that returns `rank` as a vidx passes every test built from a catalog with no
trans-only section, and is wrong on TOPCHeF. `test_overlap.py` pins the case.

### What `Store.validate` checks

14. Every overlap pointer names catalogs that are in the store, each with the pointer's recorded
    `identity_digest` and a `collection_digest` equal to the pointer's; the object's header has kind
    11, chromosome `all`, that collection digest, and count = k; `mask_width` is the smallest legal
    width for k; the directory's chromosome list is the union of the catalogs' chromosomes in ASCII
    `seq_digest` order; the union decodes to exactly the merge of the catalogs' sites; each
    catalog's set-bit count equals its `n_sites`; `prefix` agrees with a scan of the masks; and
    every breakpoint reproduces the vidx the catalog's own pages give for that site.

`gc` never deletes an object an overlap pointer names.

### Scale

The mask ceiling is **64 catalogs**; past that the masks should become one bitmap per catalog and
this section needs rewriting rather than stretching. The union grows sub-linearly in the number of
catalogs — on chr22 two catalogs union to 1.08x the larger of them — so it approaches the variant
set of the imputation panel the studies share, not the sum of their sites.

## 20. Coloc lists

Which phenotypes of one QTL experiment someone found to colocalize with one GWAS experiment. A list
is **optional** and exists only for a pair someone declared; the store never holds coloc for every
pair, and the pipeline never runs coloc. It is **membership only**: no PP.H4 and no lead variant,
since a reader computes both from the two experiments' own rows when it opens a locus. What a list
adds is which loci to point at without computing every candidate.

### `colocs/<id>.json`

```json
{"id": "murray2026_ge_dcm", "qtl": "topchef_ge_cis", "gwas": "dcm_jurgens2024_biobanks",
 "threshold": "PP.H4 > 0.8",
 "source": {"reported_by": "Murray et al. 2026, ...", "table": "Table 2 (eQTL) ...", "method": "coloc.abf ..."},
 "note": "PJVK and CDKN1A are not eGenes by the permutation rule; ...",
 "object": "<digest>.json", "n": 21}
```

- `qtl` names a `qtl_cis` or `qtl_trans` experiment, `gwas` a `gwas` one. They need not share a
  variant catalog; a reader joining their windows compares catalogs itself (section 11).
- `source` is free-form and says who produced the list and how: reported by a publication, or
  computed by someone, with the method. `threshold` is the rule as text; `note` is optional.
- `object` is a content-addressed JSON array of `{"gene_id", "phenotype_id"}`, sorted by
  `(gene_id, phenotype_id)` and unique. `gene_id` is always present; `phenotype_id` is one of the
  QTL experiment's phenotypes, or null for a result reported at gene level (TOPCHeF's sQTL list
  names 4 genes, not the 17 introns).

**Input** (`pipeline/coloc.py`, a `colocs:` entry in config.yaml): coloc results as commonly written,
CSV or TSV, read by column name. A phenotype (`phenotype_id`, `molecular_trait_id`, `trait_id`), a
gene (`gene_id`, `ensembl_gene_id`, `ensembl_id`, version stripped) or a symbol (`symbol`,
`gene_name`, `gene_symbol`, `hgnc_symbol`, which must name exactly one gene of the QTL experiment's
annotation) identifies a row; a PP.H4 column (`PP.H4.abf`, `PP.H4`, `pp_h4`, `PPH4`, `H4`) filters
by the entry's `min_pp_h4` and is not stored. Every row must resolve to a gene with a phenotype in
the QTL experiment, or the list is refused.

### What `Store.validate` checks

15. Every coloc pointer names a QTL experiment and a GWAS experiment in the store, has `source` and
    `threshold`; its object decodes to `n` rows, sorted and unique; every `gene_id` is the gene of a
    phenotype of the QTL experiment, and every `phenotype_id` is a phenotype of that gene there. The
    last check is what keeps a list from going stale silently: a rebuilt experiment that lost a
    phenotype fails its list.

`gc` never deletes an object a coloc pointer names. `remove-experiment` leaves a list naming the
removed experiment, which then fails check 15; delete the pointer by hand, as for a catalog.
