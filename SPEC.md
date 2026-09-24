# qtlb format, version 1: the qtlstore

A **qtlstore** is a directory of content-addressed objects holding N experiments over M variant
catalogs, every position anchored to a refget sequence. It is read over plain HTTP range requests,
like v0. This document is the byte-level contract for what the builders write today:

| module | writes |
|---|---|
| `pipeline/qtlstore.py` | object names, the 64-byte v1 header, variant catalog identity, orientation, `Store.validate`, store maintenance (`remove-experiment`, `gc`, `crosscat`) |
| `pipeline/catalog.py` | a variant catalog: per-chromosome variants files, variant index, rsID index |
| `pipeline/annotation.py` | a gene/exon annotation object, its per-chromosome genes and exon models, and its gene lookup |
| `pipeline/results.py` | an experiment: results files, hits files, search index and its per-chromosome parts, trans objects |
| `pipeline/gwas.py` | an experiment's GWAS object: per-chromosome GWAS files, the GWAS index and the bin summary |
| `pipeline/dof.py` | the degrees-of-freedom fit an experiment records |
| `pipeline/packfmt_v1.py` | the v1 codec pieces (zstd framing, variant page header, results block, quantizers, SNP codes, hits, trans and GWAS frames) |
| `pipeline/CONTRACT.md` | the input tables every builder reads (the ingestion contract) |

`pipeline/packfmt_v0.py` is the v0 codec; no v1 module imports it (`test_qtlstore.py` checks). The two
share code by copy: where v1 reuses a v0 layout, the v1 function is the v0 one byte for byte.

When this document and the code disagree, one of them is wrong and both are fixed in the same
change. Section 18 lists what is still open.

Measurements behind the design choices (page size, codec, dof, precision, the reference allele
checks) are in the analysis repo, `qtlb-format/docs/EVIDENCE.md`. They were taken on v0 files;
where v1 reuses a v0 layout unchanged (variant pages, result blocks, quantizers) they carry over.

**Terms.** A **variant catalog** is the list of variants a study's results point to: for each variant
its sequence digest, position, `ref` and `alt`, stored as per-chromosome variants files, plus the
variant index (vidx) and the rsID index. It is not the **eQTL Catalogue**, the EBI resource one of
the experiments comes from; this document spells that one out in full every time.

**Status.** Builders and `Store.validate` exist and pass on a genome-wide two-experiment store
(`/scratch/ns5bc/qtl-browser/store-genome-v1f` on Rivanna: TOPCHeF, with its trans results and the DCM
GWAS, and GTEx v8 heart LV). Marked **(future)** below: the VRS-id index and `store.json` refget URLs.
A browser reader for v1 is being written in `ui/`; this document is its reference.

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
  store.json                    mutable: name, format version, ids at each pointer level, refget URLs
  immutable/<digest>.<ext>      every object; <digest> = sha512t24u of the object's bytes
  variant_catalogs/<id>.json    mutable pointer: one variant catalog
  annotations/<id>.json         mutable pointer: one gene annotation release
  experiments/<id>.json         mutable pointer: one experiment (variant catalog + annotation + results)
```

Three mutable levels; everything else is content-addressed. This mirrors a refgetstore
(`rgstore.json` -> `collections/` -> `sequences/`) on purpose.

**Write order.** Objects first, then pointers, then `store.json`. `Store.write_pointer` refuses a
pointer that names an object not yet in `immutable/`. `Store.write_store` lists whatever pointer
files exist when it runs. A reader that starts from any level never meets a dangling name.

**`store.json`:**

```json
{"name": "store-genome-v1f", "format_version": 1, "refget": [],
 "variant_catalogs": ["gtex_v8_heart_lv_grch38", "topchef_grch38"],
 "annotations": ["gencode_v34", "gencode_v39"],
 "experiments": ["gtex_v8_heart_lv", "topchef"]}
```

(The genome-wide two-experiment store on Rivanna, `/scratch/ns5bc/qtl-browser/store-genome-v1f/`.)
The pointer level of variant catalogs is the directory `variant_catalogs/` and the `store.json` key
`variant_catalogs` (stores built before 2026-09-24 used `catalogs`). An experiment names its variant
catalog in its own key `catalog` (section 8).

Id lists are sorted. `refget` is a list of refgetstore URLs; `store.sbatch` writes `[]` today
**(future: filled)**.

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
| `qbe` | results file, one per (experiment, phenotype type, chromosome) | 2 |
| `qbg` | GWAS file, one per (experiment GWAS, chromosome) | 4 |
| `qgi` | GWAS index, one per experiment GWAS | 5 |
| `qbt` | trans object, one per (experiment, phenotype type) | 6 |
| `qbh` | hits file, one per (experiment, chromosome) | 7 |
| `qbr` | rsID index, one per variant catalog | 8 |
| `qbx` | variant index, one per variant catalog | 9 |
| `qgl` | gene lookup, one per annotation | 10 |
| `arrow.zst` | annotation gene and exon tables and their per-chromosome genes and exon models, experiment search index and its per-chromosome parts, GWAS bin summary | no header |

Kind 3 is unused in v1 (v0's detail-less sQTL kind). Kinds 4-6 keep v0's numbers for the same
content (GWAS, GWAS index, trans) with v1 layouts. The codec module is `pipeline/packfmt_v1.py`
(plan step 8's separate v1 codec); kinds were not renumbered.

A pointer document names objects by bare `<digest>.<ext>` strings. `qtlstore.object_names`
collects every string anywhere in a pointer that matches the name pattern; those are the objects
`Store.validate` resolves.

## 4. File header, v1 (64 bytes)

Every binary object (`qbv`, `qbe`, `qbg`, `qgi`, `qbt`, `qbh`, `qbr`, `qbx`, `qgl`) starts with this header.
Struct `<4sBBH8sIII32s4x>` (`qtlstore._HEADER`).

| offset | type | field |
|---|---|---|
| 0 | 4 bytes | magic `QTLB` |
| 4 | u8 | kind, 1..255 (section 3 table) |
| 5 | u8 | version: 1 |
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

`seq_digest` is the one field v1 adds: a file found alone says which exact sequence its positions
index. For the kinds whose chromosome is `all` (5, 6, 8, 9) the same 32 bytes carry the seqcol
collection digest instead, since they span sequences. The gene lookup (kind 10) indexes no sequence,
so they carry the identity digest of the annotation it was built from.

**What a reader checks.** `Store.validate` checks every header against its pointer at build time
(section 14). The browser does not send a separate request for a header: it checks a header only
where one arrives inside bytes it reads anyway (the gene lookup's directory, a hits file's frame
table, a whole-object read). A range read at a pointer-given offset trusts the pointer, which is
what the content-addressed name and `validate` guarantee, and the decoders reject bytes that are
not what the offsets promise (block magic and length, zstd framing, record counts).

**Reader rules** (`parse_file_header`): magic `QTLB`; version 1; header length 64; bytes 60-63
zero; chromosome non-empty, ASCII, no inner NUL; `seq_digest` matches `^[A-Za-z0-9_-]{32}$`.

**Example**, kind 1, `chr1`, count 3, page size 512, `n_cis` 2, GRCh38 chr1
(`file_header(1, "chr1", 3, 512, 2, "Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO")`):

```
00: 51 54 4c 42 01 01 40 00 63 68 72 31 00 00 00 00   QTLB, kind 1, v1, hlen 64, "chr1"
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
  position in this list. There is no fixed chr1..chrX list in v1.
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
(`QVX0`) there are no reserved fields, no `n_frames` and no `hits_off`, because v1 hits are not
paged (section 8). Chromosome names are not stored: they are the variant catalog JSON's table, in order.

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
  transcript-level table split by chromosome would be 1.0 MB for chr1, which is why the per-chromosome
  object holds the model rather than the rows. Genes and models are two objects so a reader of gene
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

One experiment = one cohort and tissue, possibly several phenotype types, over one variant catalog and one
annotation (`results.build`). The phenotype type (`ge`, `leafcutter`, ...) is data, not a file
kind: one code path and one block layout serve every type.

### `experiments/<id>.json`

```json
{"id": "topchef",
 "catalog": "topchef_grch38",
 "catalog_identity": "<32>",
 "annotation": "gencode_v34",
 "annotation_source": null,
 "allele_orientation_source": "topchef_refcheck_a2_is_ref",
 "significance": {"column": "p_perm", "op": "<", "threshold": 0.05},
 "search_index": "<digest>.arrow.zst",
 "search_index_parts": {"chr21": {"file": "<digest>.arrow.zst", "rows": 986, "ords": [[0, 189], [703, 1498]]}},
 "search_index_trans_only": {"file": "<digest>.arrow.zst", "rows": 53, "ords": [[686, 702], [3565, 3600]]},
 "counts": {"ge": {"phenotypes": 686, "with_rows": 686, "significant": 372, "significant_genes": 372},
            "leafcutter": {"phenotypes": 2862, "with_rows": 2862, "significant": 512, "significant_genes": 204}},
 "n_phenotypes": 0,
 "unplaced": {"count": 0, "examples": []},
 "hits": {"chr21": "<digest>.qbh"},
 "results": [
   {"phenotype_type": "ge", "dof": 435, "dof_fit": null, "dof_source": "ingestion.json",
    "n_phenotypes": 0, "files": {"chr21": "<digest>.qbe"},
    "precision": {"neglog10p_max_error": 0.0, "slope_se_max_rel_error": 0.0,
                  "slope_max_error_over_se": 0.0, "slope_rows_compared": 0, "af_max_error": 7.63e-06},
    "trans": {"file": "<digest>.qbt", "n_rows": 0, "n_phenotypes": 0,
              "precision": {"neglog10p_max_error": 0.0, "beta_max_error": 0.0, "af_max_error": 7.63e-06}}},
   {"phenotype_type": "leafcutter", "dof": 480, "...": "..."}],
 "trans": {"rows": 0, "rows_in_source": 0, "rows_skipped_variant_outside_build": 0,
           "rows_skipped_phenotype_outside_build": 0},
 "trans_excluded": {"variants": 32765, "trans_eqtl_rows": 37951, "...": "..."},
 "gwas": {"id": "dcm_jurgens2024_biobanks", "title": "...", "files": {"chr21": "<digest>.qbg"},
          "index": "<digest>.qgi", "n_rows": 0, "rows_by_chrom": {}, "rows_sharing_a_site": 0,
          "block_rows": 2048, "n_values": [],
          "bins": {"file": "<digest>.arrow.zst", "bin_bp": 5000000, "n_bins": 0},
          "orientation": {}, "source": {}},
 "source": {"experiment_id": "topchef", "source": {}}}
```

- `catalog_identity` copies the variant catalog's `identity_digest` at build time.
- `annotation_source` is the source's own annotation when it differs from the one attached
  (a known gap, reported), else null.
- `significance` comes from the adapter's `ingestion.json`; default `p_perm < 0.05`. `column` is
  `p_perm` or `p_beta` (the `permuted` column the rule tests), `op` is `<` or `<=`; the builder
  refuses anything else. Every `significant` below means "the group's `column` passes the rule".
- `unplaced` counts phenotypes that no nominal row, permuted group or credible set places on a
  chromosome.
- `results[].files` has one entry per chromosome built, including chromosomes with no blocks (a
  header-only file).
- `precision` per results set: `neglog10p_max_error` = the largest `nlp_max / 131066` over its
  blocks, `slope_se_max_rel_error` = the largest `expm1((lse_max - lse_min) / 65532)`,
  `af_max_error` = `0.5 / 65534` (worst-case bounds, section 13), and
  `slope_max_error_over_se`, **measured**: the largest |slope rebuilt from the stored codes - source
  beta| / source se over every row the builder can rebuild a slope for (`slope_rows_compared`
  rows: p above 0, SE present, finite source beta and se). Null when `dof` is null.
- `results[].trans` is the phenotype type's trans object (section 9), null when it has no trans rows.
  Top-level `trans` counts the rows (a subset build skips rows whose variant or phenotype lies outside
  it, counted); it is null when the experiment has no trans table. `trans_excluded` copies the
  adapter's count of trans rows left out because their variants have no source alleles (section 12).
- `gwas` is the experiment's GWAS object (section 10), or null.

### Degrees of freedom

`results[].dof` is what a reader uses to rebuild a slope. From `ingestion.json` `dof[phenotype_type]`:

- an integer (TOPCHeF publishes 435 for `ge`, 480 for `leafcutter`): stored as is, `dof_fit` null;
- a `dof.fit` result (eQTL Catalogue): stored when `usable`, else **`dof: null`**. `dof_fit` keeps
  `residual_log10p, margin, rows_usable, n_samples, implied_covariates, reason`.

The fit (`dof.fit`): integer grid `n_samples - 200 .. n_samples - 1` (widened when the minimum sits
within 3 of an edge), objective = median `|log10 p - log10(2 * stdtr(dof, -|beta/se|))|` over rows
with finite beta, `se > 0` and `1e-300 < p < 0.5`. It runs on a uniform sample (100,000 rows) and on
the same count from the large-|t| tail. `residual_log10p` is the larger of the uniform residual and
the tail residual at the uniform dof. **usable** = residual at most 0.01 and neither minimum at a
grid edge. `identified` (reported, not gating) = both fits agree and the runner-up is at least 2x
worse on the tail.

**`dof: null` means no slope.** A reader of a results set with `dof: null` shows the stored
-log10 p, SE and the sign bit, and no slope value (`results.read_block` returns NaN slopes). It
must not substitute any other dof: a stand-in dof decodes to slopes that look valid and are wrong.

### Results file (kind 2, `.qbe`)

`[64-byte header][block][block]...`, one block per phenotype of that type on that chromosome,
back to back, in search-index `ord` order. Header count = blocks, page size 0, `n_cis` 0,
`seq_digest` = the chromosome's.

**Block**: the v0 layout (`packfmt_v1.encode_gene_block`, a byte-for-byte copy of v0's), unchanged in bytes:

| offset | type | field |
|---|---|---|
| 0 | 4 bytes | magic `QGB0` |
| 4 | u32 | block length including padding |
| 8 | u32 | `n_rows` |
| 12 | u32 | `var_start` (variant catalog vidx of row 0; `0xFFFFFFFF` when `n_rows` is 0) |
| 16 | i32 | `anchor`: **always 0 in v1** |
| 20 | u32 | `pos_first` (0 when `n_rows` is 0) |
| 24 | u32 | `pos_last` (0 when `n_rows` is 0) |
| 28 | u32 | `n_cs` |
| 32 | f64 | `nlp_max` |
| 40 | f64 | `lse_min` |
| 48 | f64 | `lse_max` |
| 56 | u32 | `details_zlen` |
| 60 | u32 | `details_len` |

Then `n_rows` pairs `{u16 nlp code, u16 SE code}`, `n_cs` credible-set records `{u32 row, f32 pip,
u8 cs_id, 3 zero bytes}` sorted strictly by `(row, cs_id)`, one zstd details frame, zero padding to
4. Block length = `64 + 4 n_rows + 12 n_cs + details_zlen`, rounded up to 4.

- **Rows** are the vidx run `var_start .. var_start + n_rows - 1` covering every nominal row and
  every credible-set site of the phenotype. A vidx in the run the phenotype did not test is a
  **null row**: nlp code 65535, SE code 0xFFFF. (For TOPCHeF every run is contiguous and has no
  null rows.) A run never crosses the cis/trans-only boundary; the builder raises if it would.
- **A phenotype with no nominal rows and no credible sets** still gets a block with `n_rows` 0 and a
  search-index row, so a reader can say "no per-variant data published".
- **Every v1 block has a details frame** (v0's detail-less kind 3 sQTL block is gone).
- **`anchor` is 0.** v0 stored the tensorQTL window start there. v1 treats the TSS as annotation:
  a reader takes it from the annotation object by `gene_id`.
- Codes, scales, the scale rule and the reject rules are v0's (section 13; `decode_gene_block`).

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
`column` (`p_perm` or `p_beta`).
Readers reject a details object whose `v` is not 1.

### Search index (`.arrow.zst`)

One row per phenotype, all phenotype types, in `ord` order. Sorted by (phenotype type in the
ingestion's `phenotype_types` order, chromosome in variant catalog table order with trans-only
phenotypes last, position of the run's first site or else the group lead's position,
`phenotype_id`). Schema (`results.INDEX_SCHEMA`):

| column | type | meaning; null when |
|---|---|---|
| `ord` | uint32 | row number; never null |
| `phenotype_type`, `phenotype_id`, `phenotype_object_id` | string | |
| `gene_id` | string | join key to the annotation; null when the source names no gene |
| `chr` | string | the chromosome its cis results are on; **null for a trans-only phenotype** |
| `has_nominal` | bool | the phenotype has nominal rows |
| `is_group_lead` | bool | this phenotype is its group's lead |
| `significant` | bool | the group passes the experiment's rule (on its `column`); false with no group |
| `p_perm` | float64 | group's; null with no group |
| `blk_off`, `blk_len` | uint32 | the block in `results[type].files[chr]`; **null for a trans-only phenotype** |
| `var_start`, `n_var` | uint32 | the run; null when `n_rows` is 0 |
| `var_off`, `var_len` | uint32 | variant catalog `.qbv` bytes of the pages covering the run; null likewise |
| `w_lo`, `w_hi` | int32 | positions of the run's first and last sites; null likewise |
| `trans_off`, `trans_len` | uint32 | the phenotype's frame in `results[type].trans.file` (section 9; a gene's frames are adjacent); null when it has no trans rows |
| `n_trans` | uint32 | trans rows in that frame; null likewise |

A **trans-only phenotype** is one the source reports only in trans (TOPCHeF: 17 genes, most on chrM,
and 36 introns): it is in `phenotypes`, has no nominal rows, no group and no credible set, so nothing
places it on a chromosome. Its row has `chr` and `blk_*` null and the trans pointers set; there is
no block for it in any results file. Only trans-only phenotypes have `chr` null.

No annotation columns, and no schema metadata listing pack hashes (v0's `qtl_browser.packs`):
content-addressed names already pin every object. When `phenotypes.has_nominal` exists the
builder fails if it disagrees with the nominal rows it found.

**Parts.** The same rows split for readers that need one chromosome (`results.split_index`, which
`build` and `add-split` both call): `search_index_parts[chr]` holds the rows whose `chr` is that
chromosome, one part per chromosome built (empty ones included), and `search_index_trans_only`
the rows whose `chr` is null (null when there are none). Each is an Arrow object with the search
index's schema and its rows in `ord` order; `rows` is its row count and `ords` its rows' `ord`
values as sorted inclusive `[first, last]` runs, so a reader holding an `ord` (from a hits record)
knows which part to read without reading any. Because rows are ordered by phenotype type first, a
chromosome's part has one run per phenotype type. TOPCHeF: chr1 322 KB, chr7 166 KB, all 3.0 MB.

A gene page reads the part of its gene's annotation chromosome, so a phenotype placed on another
chromosome than its gene would not show there (TOPCHeF and GTEx have none).

**`counts`**, per phenotype type, over its phenotypes with a cis result (`chr` set) (`results.type_counts`):
`phenotypes`, `with_rows` (`n_var` set), `significant`, and `significant_genes` (distinct non-null
`gene_id` among the significant). The summary line a landing page prints, so it reads no index.

### Hits file (kind 7, `.qbh`)

What the variant page lists for one variant: the groups it leads, the credible sets it is in, and
its trans associations, keyed by the variant catalog's vidx. One file per chromosome, **paged by
vidx**: frame `g` holds the records of vidx `g * F .. (g + 1) * F - 1`, F = 1,024.

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
| 8 | f32 | `value`: kind 0 `p_perm` (NaN when null), kind 1 PIP, kind 2 `f32(-log10 p)` of the source p (+inf for p = 0) |
| 12 | f32 | `beta`: kind 2 the ALT effect; NaN on kinds 0 and 1 |
| 16 | u8 | kind: 0 lead of a group, 1 credible-set member, 2 trans association |
| 17 | u8 | `cs_id` (kind 1), 0 otherwise |
| 18 | u8 | flags: bit 0 significant by the experiment's rule (kind 0 only); other bits zero |
| 19 | u8 | zero |

Kind 0: `ord` is the group's lead phenotype. Kind 1: a site in two sets of one phenotype keeps both
records. Kind 2: one record per trans row whose variant is this one (section 9): `ord` is the trans
phenotype (possibly trans-only). `value` is `-log10(pvalue)` of the contract `trans` row, computed in
float64 from the source p and rounded to the nearest f32 (+inf for p = 0; p itself would underflow
f32); `beta` is the source beta rounded to f32. These are not the trans frame's values: the frame
(section 9) stores the same row's -log10 p as a 16-bit code, which decodes to within `nlp_max / 131066`
of the source's (`nlp_max` of that phenotype's frame), so a kind 2 `value` and the frame's decoded
-log10 p for one row differ by at most `nlp_max / 131066 + |value| * 2^-24` (the second term is the
f32 rounding, at most 1.8e-5 since a finite -log10 p is at most 300), and both are +inf when p = 0.
Likewise `beta` and the frame's decoded beta differ by at most `beta_max / 65534 + |beta| * 2^-24`.
A reader derives the SE as `|beta| / t(p, dof)` with its phenotype type's `dof`.

**Variant lookup:** read bytes 0 .. `64 + 4 (n_frames + 1)` once per chromosome (the header and the
table; 2.8 KB for chr1), then frame `vidx // F`, then the records with that `vidx`
(`results.hits_frame`). The frame table sits in the hits file, not in the variant catalog's variant
index as in v0, because a variant catalog is shared by experiments and hits belong to one.

*Changed 2026-09-24 from the first v1 layout* (one zstd frame per chromosome of 16-byte records with
no `beta` field and no kind 2): records are 20 bytes, the file is paged, and the header's page size
and byte 24 are no longer 0.

## 9. Trans results

The gene page's trans table: every trans association of one phenotype, whichever chromosome its
variants lie on. From the contract `trans` table (CONTRACT.md); `results.build_trans`.

**One object per (experiment, phenotype type)**, kind 6, `.qbt`: `[64-byte header][frame][frame]...`.
Header chromosome `all`, `seq_digest` = the collection digest (rows span every chromosome), count =
frames, page size 0. One frame per phenotype with trans rows, back to back with no padding. A
phenotype's frame is `trans_off .. trans_off + trans_len - 1` from its search-index row. One object
per type rather than per chromosome because a trans-only phenotype has no chromosome to be filed under.

**Frame order** (`results.trans_frame_order`): first the phenotypes with a `gene_id`, grouped by gene,
genes in search-index order (a gene's place is the smallest `ord` among its phenotypes, which is its
eQTL phenotype's when it has one), `ord` order within a gene; then the phenotypes with a null
`gene_id`, in `ord` order. So **a gene's frames in each object are one contiguous byte range**: the
gene page reads, per phenotype type, the single range from the smallest `trans_off` to the largest
`trans_off + trans_len` over the gene's phenotypes (`results.gene_trans_ranges`), at most one request
per object. `ord` order alone would not do this: introns sort by position, so one gene's introns
interleave with a neighbouring or overlapping gene's. The search index needs no per-gene columns for
it; `ord` and the search-index row order are unchanged.

**Frame** (`packfmt_v1.encode_trans_frame` / `decode_trans_frame`): one zstd frame; the payload is
exactly `32 + 16 n + heap_len` bytes:

| offset | type x count | field |
|---|---|---|
| 0 | 4 bytes | magic `QTT2` |
| 4 | u32 | `n` rows, at least 1 |
| 8 | f64 | `nlp_max`: the frame's largest finite -log10 p |
| 16 | f64 | `beta_max`: the frame's largest \|beta\| |
| 24 | u32 | `heap_len` |
| 28 | u32 | reserved, 0 |
| 32 | u32 x n | position delta: absolute on the first row and wherever the ordinal changes, else `pos - previous pos` (0 allowed) |
| 32 + 4n | u32 x n | `rs_number`, 0 = none (the variant catalog's) |
| 32 + 8n | u16 x n | `af` code, ALT, the variant catalog's (`rint(af * 65534)`, 65535 null) |
| 32 + 10n | u16 x n | -log10 p code (section 13: `code * nlp_max / 65533`, 65534 is p = 0; 65535 is not allowed) |
| 32 + 12n | i16 x n | beta code, `rint(beta / beta_max * 32767)`; -32768 not allowed; all 0 when `beta_max` is 0 |
| 32 + 14n | u8 x n | chromosome ordinal of the variant (1-based, the variant catalog's table) |
| 32 + 15n | u8 x n | allele code (section 5's SNP codes, 0 = heap) |
| 32 + 16n | heap | `ref + "\t" + alt + "\n"` for each code-0 row, in row order |

Rows are sorted by (ordinal, position, ref, alt), so ordinals never decrease. `beta` is the ALT
effect. The SE is not stored: a reader derives `se = |beta| / t(p, dof)` with the phenotype type's
`dof` (`results.read_trans`), none when `dof` is null. Alleles, rsID and `af` are inline so the
table needs no variants file from another chromosome (v0's reason, kept).

**No `vidx`.** A row's `(ordinal, pos, ref, alt)` names its site, so its vidx is recoverable through
the variant catalog (the variant index and a page, as in cross-catalog lookup, section 11) when a
reader needs one; the browser links a trans row by rsID or `chr:pos` and never needed it. The hits
files' kind 2 records (section 8) keep their own vidx.

*Changed 2026-09-24 from the first v1 layout* (`QTT1`: a u32 `vidx` column first, `32 + 20 n +
heap_len` bytes, frames in search-index `ord` order): dropping `vidx` saves 4 bytes a row (15.8M
TOPCHeF rows; TOPCHeF trans 232.4 MB to 197.6 MB, v0 170.6 MB, the rest being the inline alleles and
one frame per phenotype), and the gene grouping takes a gene page from a median of 3 and up to 33
frame reads to one range per object. The magic changed so a `QTT1` reader fails loudly.

**Rules the builder enforces:** every row's `(chr, pos, ref, alt)` is a site of the experiment's
variant catalog and every phenotype is in `phenotypes` (else the build fails); p is in [0, 1] and
beta finite. The same rows, keyed by variant, are the hits files' kind 2 records (section 8).

**Precision** (`results[].trans.precision`): `neglog10p_max_error` = the largest `nlp_max / 131066`
over frames, `beta_max_error` = the largest `beta_max / 65534`, `af_max_error` as for variant pages.

*v0 difference:* v0 had one frame per gene holding the gene's eQTL rows and all its introns' sQTL
rows, with no alleles and a fixed chr1..chrX variant chromosome code; v1 has one object per phenotype
type, one frame per phenotype (a phenotype need not have a gene), a gene's frames contiguous in each
object (so a gene page makes at most one request per type, v0's one request per gene), alleles inline,
and ordinals from the variant catalog's table.

## 10. GWAS

A GWAS shown next to an experiment (TOPCHeF: the DCM GWAS, Jurgens et al. 2024, biobanks-only).
From the contract `gwas` table and `gwas.json` in the experiment's tables directory (CONTRACT.md);
`gwas.build`. Listed in the experiment's `gwas` entry (section 8).

**Values are the source's, lossless**: beta, se and af to 4 decimals (`rint(x * 10000)`), p to 4
significant digits (mantissa and exponent), n from a table of distinct values: v0's GWAS rules and
codes, unchanged (`packfmt_v1.gwas_codes`). What changed is orientation: rows are **ref/alt** like
everything else in the store, `ref` the reference base(s), `beta` and `af` ALT-relative. The adapter
orients each row against the refgetstore with `qtlstore.orient_to_ref` (the source's `EA` is the
effect allele; where `EA` is the reference the alleles swap, beta negates and `af = 1 - EAFREQ`, both
exact at 4 decimals) and drops, with a count, rows whose alleles the reference does not read. For the
DCM GWAS: 12,504,079 source rows, 10,401,887 as given, 2,092,846 swapped, 9,346 dropped.

**A site may repeat.** The DCM GWAS lists 796,531 indels twice, once per allele order, with
different N and statistics (two measurements the meta-analysis did not merge). After orientation
both rows name one `(chr, pos, ref, alt)`; both are kept, and `gwas.rows_sharing_a_site` counts them
(1,593,062). Two more indels are listed twice with mirrored, otherwise identical values; oriented they
are the same row, so the adapter keeps one copy of each (`identical_rows_dropped`: 2). The builder
refuses identical rows. The DCM GWAS object holds 12,494,731 rows.

**GWAS file** (kind 4, `.qbg`, one per chromosome with rows): `[64-byte header][block][block]...`.
Header: the chromosome, count = rows, page size = B (2048, rows per block), `seq_digest` = the
chromosome's. Rows are sorted by (pos, ref, alt, rs_number, p, n, beta). Each block is one zstd
frame of B rows (the last may be shorter), back to back with no padding. Block payload, exactly
`8 + 21 n + heap_len` bytes (v0's GWAS block, with `af` for v0's `eaf` and the allele pair
`(ref, alt)` for v0's `(ea, nea)`):

| offset | type x count | column |
|---|---|---|
| 0 | u32 | `n` rows |
| 4 | u32 | `heap_len` |
| 8 | u32 x n | position delta: first absolute, then `pos[i] - pos[i-1]` |
| 8 + 4n | i32 x n | beta code, `rint(beta * 10000)` |
| 8 + 8n | u32 x n | `rs_number`, 0 = none |
| 8 + 12n | u16 x n | se code, `rint(se * 10000)` |
| 8 + 14n | u16 x n | af code (ALT), `rint(af * 10000)`, at most 10000 |
| 8 + 16n | u16 x n | `p_mant`, 1000..9999 |
| 8 + 18n | i8 x n | `p_exp`; p = `p_mant / 10^(-p_exp)` |
| 8 + 19n | u8 x n | `n_code`, index into the index's n table |
| 8 + 20n | u8 x n | allele code of (ref, alt), 0 = heap |
| 8 + 21n | heap | `ref + "\t" + alt + "\n"` per code-0 row |

**GWAS index** (kind 5, `.qgi`, one per experiment GWAS): `[64-byte header][one zstd frame]`.
Header: chromosome `all`, count = chromosomes with rows, page size = B, `seq_digest` = the collection
digest. Payload (`packfmt_v1.encode_gwas_index_payload`), all u32 except names:

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
for that is too much for a landing page, so the experiment carries a small summary: `gwas.bins` names
it, `{"file": "<digest>.arrow.zst", "bin_bp": 5000000, "n_bins": N}`, or `gwas.bins` is null when the
tables have no bin table. It is v0's `gwas_dcm_bins.json` (built by v0's `gwas_bins` step) as a v1
object: the same bins and the same values.

**Source.** The adapter writes the contract table `gwas_bins.parquet` next to `gwas.parquet`
(CONTRACT.md) with v0's query, unchanged (`dcm_gwas.BINS_SQL`): over **every source row with a GRCh38
position and p > 0**, grouped by `(chr, pos // bin_bp)`. These are the source rows, not the oriented
GWAS rows of the `.qbg` files: the 9,346 rows the reference does not read are counted, and the lead's
allele and beta are the source's `EA` and `BETA`, not ref/alt. That is what keeps it equal to v0.
`gwas.build_bins` keeps the chromosomes the GWAS object has files for and rounds as v0 did.

**Object**: one zstd frame around an Arrow IPC stream (section 1), no header, one row per bin that has
rows, sorted by chromosome in the variant catalog's table order, then `bin_start`
(`gwas.BINS_SCHEMA`):

| column | type | meaning |
|---|---|---|
| `chr` | string | chromosome name, in the variant catalog's table |
| `bin_start`, `bin_end` | uint32 | the bin `[bin_start, bin_end)`, `bin_start = k * bin_bp`, `bin_end = bin_start + bin_bp`; 0-based like v0's, so position `pos` falls in bin `pos // bin_bp` |
| `min_p` | float64 | the bin's smallest source p, to 3 significant digits (`float(f"{p:.3g}")`) |
| `lead_position` | uint32 | position of the row holding it (inside the bin) |
| `lead_rsid` | string | that row's source `rsID` text, null when none |
| `lead_beta` | float64 | that row's source `BETA` (effect of `EA`), to 3 decimals (Python `round(b, 3)`) |
| `lead_ea` | string | that row's source `EA`, as written |
| `n_gws` | uint32 | rows with p < 5e-8 |
| `n_variants` | uint32 | rows in the bin |

The DCM GWAS: 569 bins genome-wide, a 16,507-byte object (18 on chr21 and chr22, 1,125 bytes). **Ties:** when several rows share the
smallest p, the lead fields come from DuckDB's `arg_min` over the file as read, as in v0; the four
`arg_min` calls are separate, so on a tie they are not guaranteed to name one row (v0 behaves the
same). A reader uses `min_p`, `n_gws` and `n_variants`; the lead fields label the bar.

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

## 12. What v1 leaves out

- **Trans eQTL variants without source alleles.** TOPCHeF's trans eQTL files report 32,765
  positions as bare `chr:pos`. A position is not a variant, so these get no site and none of their
  trans eQTL rows is ingested; alleles are not inferred from dbSNP, frequency or other studies
  (CONTRACT.md, decided 2026-09-24). `ingestion.json` `trans_eqtl_excluded` counts the variants, the
  rows left out and the genes. They enter when the authors supply the alleles.
  Every other trans row is in v1: all 13,182,408 trans sQTL rows and 2,642,166 of 2,680,117 trans
  eQTL rows (the variant has alleles from another file of the release: one tested variant per
  position, joined on `chr:pos`). The experiment JSON's `trans_excluded` repeats the count.
- **GWAS rows the reference does not read.** 9,346 of the DCM GWAS's 12,504,079 rows have alleles
  that match neither reference reading; they are dropped and counted in `gwas.orientation`.

The Python reference decoders are `results.read_block`, `results.read_trans`, `results.hits_frame` /
`decode_hits`, `catalog.decode_file`, `decode_vidx`, `decode_rsid`, `rsid_lookup`, `gwas.read_window`,
`gwas.read_bins` and `annotation.load`.

## 13. Precision

Unchanged from v0. The quantizers are v0's (`packfmt_v1.quantize_nlp`, `quantize_se`, copies of v0's),
so the v0 round-trip evidence (EVIDENCE.md A.1, A.2) applies to the codes. The lossy-versus-lossless
question is closed.

- **`-log10 p`:** code 0..65533 gives `code * (nlp_max / 65533)`; 65534 is p = 0; 65535 is null.
  The writer rejects a finite -log10 p above 300. `nlp_max` is the block's largest finite -log10 p
  (0.0 when every p is 1 or null). Worst error half a step, `nlp_max / 131066`.
- **SE and slope sign:** bit 15 set when slope < 0 (0.0 and -0.0 are positive); bits 0-14 a log(SE)
  code 0..32766, `se = exp(lse_min + q * (lse_max - lse_min) / 32766)`; 0xFFFF null (SE or slope
  null); 0x7FFF invalid. Worst relative error `expm1((lse_max - lse_min) / 65532)`.
- **Scale rule:** when `nlp_max > 0` some row holds 65533, else every finite code is 0; when
  `lse_min < lse_max` some row holds 0 and some 32766; when equal every non-null code is 0; when no
  row has an SE both are 0.0.
- **Slope:** `slope = sign * se * t`, `t = max(0, -stdtrit(dof, p / 2))`, null when p is null, p is
  0, the SE code is null, **or `dof` is null**. D13 (store SE, derive slope) stands.
- **`af`:** `rint(af * 65534)`, error at most 7.63e-6.
- **Trans** (section 9): -log10 p as above per frame; beta `code * beta_max / 32767`, worst error
  `beta_max / 65534`; SE derived.
- **GWAS** (section 10): lossless at the source's printed precision.
- **Exact:** positions, alleles, `rs_number`, counts, credible-set rows and `cs_id`; `pip` and hits
  values as f32.
- Measured v0 maxima (TOPCHeF, every row; EVIDENCE.md A.2): -log10 p 1.27e-3; `slope_se` relative
  1.41e-5 (eQTL), 1.93e-5 (sQTL); slope at most 4.47e-3 of the row's `slope_se`.

**Slope error, measured per results set.** The slope is the one value a reader derives rather than
decodes, so its error has no closed-form bound in the block header: it depends on how well `beta /
se` agrees with `p` under the stored `dof` in the source itself, plus the two rounding steps. The
builder therefore measures it: every block is decoded again after encoding and every rebuilt slope
is compared to the source beta (`results.slope_error`); the experiment records the largest
`|slope - beta| / se` as `precision.slope_max_error_over_se` with the row count
(`slope_rows_compared`). Genome-wide, `store-genome-v1e`:

| experiment | type | dof | rows compared | `slope_max_error_over_se` |
|---|---|---:|---:|---:|
| TOPCHeF | `ge` | 435 | 123,458,053 | 3.52e-3 |
| TOPCHeF | `leafcutter` | 480 | 499,596,989 | 4.47e-3 |
| GTEx v8 heart LV | `ge` | 367 (fitted) | 153,275,168 | 3.30e-3 |
| GTEx v8 heart LV | `leafcutter` | 367 (fitted) | 16,489,572 | 2.91e-3 |

TOPCHeF's two maxima equal v0's round-trip numbers to every printed digit (v0 manifest
`precision.*.slope_max_error_over_se`: 3.518e-3, 4.466e-3), as they must, since the codes are
identical. For the eQTL Catalogue, with a fitted dof, the measured error is at the same level. **So a
rebuilt slope is within 0.5% of the row's SE** in every results set built so far. The check costs a
second decode of every block at build time (the TOPCHeF results phase went from 5.4 to 19.8 minutes
with the trans objects and the slope measurement together).

## 14. Validation

`Store.validate(refget, sites)` returns failures; `Store.notes` lists checks it could not run. **An
empty failure list with non-empty notes is a partial pass.**

1. `store.json` `format_version` is 1; every listed id has its pointer file.
2. Every object any pointer names exists and its bytes hash to its name.
3. Every variant catalog chromosome's `.qbv` header parses, and its chromosome and `seq_digest` equal the
   table entry.
4. With a refgetstore: every variant catalog `seq_digest` is a sequence the store holds, at the recorded
   `length`. Without one: a note. (Without this check a variant catalog that carries the same wrong digest
   everywhere passes the rest.)
5. `identity_digest` is a 32-character digest; with a site reader (`catalog.read_sites`) it
   recomputes from the decoded variants files. Without one: a note.
6. Every experiment's variant catalog and annotation are in the store, and its `catalog_identity` equals
   the variant catalog's `identity_digest`.
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
    `lead_ea` are present (`gwas.check_bins`, which the builder also runs).
11. Every trans object is tiled by its frames as the search index places them (from byte 64 to the
    end, no gap or overlap, as many frames as the pointer's `n_phenotypes`), and each gene's frames
    in it form one contiguous range (`results.check_trans_layout`, section 9).
12. Every annotation has `chroms` and `lookup`: one `chroms` entry per chromosome of `genes`, each
    `genes` object equal to that chromosome's rows and each `exon_models` object equal to the models
    recomputed from `exons`; the lookup has kind 10, chromosome `all` and the annotation's
    `identity_digest`, its offsets tile it, and every bucket equals what `lookup_buckets` makes of
    `genes` (`annotation.check_split`).
13. Every experiment with a search index has `search_index_parts`, `search_index_trans_only` and
    `counts`: one part per chromosome built, each equal to that chromosome's index rows with `rows`
    and `ords` describing it, the trans-only part likewise, and `counts` equal to `type_counts` of
    the index (`results.check_split`).

`python -m pipeline.annotation add-split --store DIR --id ID` and `python -m pipeline.results
add-split --store DIR --id ID` give a store built before 12 and 13 its derived objects, from the
objects it already has, by the same code the builders run; the pointer is rewritten and no other
object changes.

Every header check also checks the kind: a variants file must be kind 1.

Build-time checks outside `validate`: orphan nominal or credible-set rows (a site not in the
variant catalog) fail the build; a repeated `(phenotype, cs_id, site)` credible-set row fails; one
phenotype's nominal rows on two chromosomes fail; `cs_id` outside 0..127 fails; nominal rows for a
phenotype missing from `phenotypes` fail; a trans row whose site is not in the variant catalog or whose
phenotype is not in `phenotypes`, or with a null p or beta, fails; a GWAS value with more than 4
decimals (or p with more than 4 significant digits) fails; identical GWAS rows fail.

## 15. Store maintenance

`python -m pipeline.qtlstore <command> --store DIR`:

- `validate`: section 14, with the configured refgetstore and the site reader.
- `remove-experiment ID`: deletes `experiments/ID.json` and rewrites `store.json` (same name and
  refget list). Nothing else changes: the experiment's objects stay until `gc`, and its variant
  catalog and annotation pointers stay, since another experiment may name them.
- `gc [--dry-run]`: deletes every file in `immutable/` whose name no pointer file names (a pointer on
  disk protects its objects whether or not `store.json` lists it), and leftover `*.tmp` files from
  an interrupted write. Files that are not object names are left alone.
- `crosscat A B`: section 11's check.

Because objects are named by their bytes, adding or removing an experiment never changes another
experiment's objects: building TOPCHeF alone, adding GTEx, removing GTEx and running `gc` leaves every
TOPCHeF object byte for byte as it was (pipeline/README.md, "Experiment modularity").

## 16. Differences from v0

| | v0 | v1 |
|---|---|---|
| scope | one study | N experiments over M variant catalogs |
| naming | `<stem>.<16 hex of sha256>.<ext>`, `manifest.json` | `<sha512t24u>.<ext>`, three pointer levels + `store.json` |
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
| GWAS bins | `gwas_dcm_bins.json` at the bucket root, columnar JSON, fixed name | an Arrow object named from `gwas.bins`, same bins and values (section 10) |
| variant catalog identity | none | `identity_digest` |
| what a first gene page reads | the whole `search_index` (1.9 MB) | one lookup bucket, its chromosome's genes, exon models and index part (sections 6, 8) |

## 17. Where v0 lives

The v0 spec (sections 1-15, byte layouts for GWAS, trans and paged hits included) is in the
analysis repo's git history at `qtlb-format/docs/SPEC.md`, commit `fe5a606`. Its Appendix A is now
`qtlb-format/docs/EVIDENCE.md`.

The frozen v0 build is on Rivanna at `/scratch/ns5bc/qtl-browser/derived/`: the packs in
`immutable/`, `manifest.json`, and the v0 tables in `_tables/`. The live site no longer serves it.
No v1 builder reads v0 packs; TOPCHeF is rebuilt into the store from its contract tables. Only the
comparisons read that tree: `pipeline/verify_v0.py` (store.sbatch), the TOPCHeF acceptance gate
`pipeline/adapters/verify_topchef.py` (adapter.sbatch), and `pipeline/bench_store.py`, through
`packfmt_v0.py` and the v0 reader `packtool.py`. The code that built v0 (the pack steps, `packcheck`,
the R2 upload) is in this repo's git history at commit `c62bca3`.

## 18. Open

1. `store.json` `refget` URLs are written empty; the VRS-id index is reserved, not built.
2. The v0 bridges (`pipeline/verify_v0.py`, `pipeline/bench_store.py`, `pipeline/adapters/verify_topchef.py`,
   with `packfmt_v0.py` and `packtool.py`) go away with the frozen v0 build.
