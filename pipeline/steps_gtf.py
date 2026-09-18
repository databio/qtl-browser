"""Step 2: gene annotation and exons from the GENCODE v34 GTF."""
import gzip
import re

import pyarrow as pa

from .common import Config, log, write_parquet

ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_attrs(s: str) -> dict:
    return dict(ATTR.findall(s))


def run(cfg: Config) -> None:
    genes, exons = [], []
    with gzip.open(cfg.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("gene", "exon"):
                continue
            a = parse_attrs(f[8])
            gid_v = a["gene_id"]
            if gid_v.endswith("_PAR_Y"):
                continue
            start, end, strand = int(f[3]), int(f[4]), f[6]
            if f[2] == "gene":
                genes.append({
                    "gene_id": gid_v.split(".")[0], "gene_id_version": gid_v, "symbol": a.get("gene_name"),
                    "chr": f[0], "start": start, "end": end, "strand": strand,
                    "tss": start if strand == "+" else end, "biotype": a.get("gene_type"),
                })
            else:
                exons.append({
                    "gene_id": gid_v.split(".")[0], "transcript_id": a["transcript_id"],
                    "exon_number": int(a.get("exon_number", 0)), "chr": f[0], "start": start, "end": end, "strand": strand,
                })
    gt = pa.Table.from_pylist(genes)
    # exons are read per gene by the eQTL pack builder, which puts them in the gene's details
    # JSON: sort by gene so each row group's gene_id range is tight and one gene is one or two
    # small reads. GTF order (by position) leaves gene_id statistics spanning most of the file.
    et = pa.Table.from_pylist(sorted(exons, key=lambda e: (e["chr"], e["gene_id"], e["start"], e["end"])))
    write_parquet(gt, cfg.tables / "gene_annotation.parquet", row_group_size=100_000)
    write_parquet(et, cfg.tables / "exons.parquet", row_group_size=5_000, stats_columns=["gene_id", "chr", "start", "end"])
    log(f"gtf: {gt.num_rows} genes, {et.num_rows} exons")
