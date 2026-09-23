#!/bin/bash
#SBATCH --job-name=qtlb_store
#SBATCH --output=/home/ns5bc/scratch/qtlb_store_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#
# A qtlb v1 store from contract tables: annotation (step 6), then per experiment a catalog (step 5)
# and results (step 7), then Store.validate. Writes only $QTLB_STORE. Chromosomes from QTLB_CHROMS
# (default chr21,chr22). EXPERIMENTS is a space-separated list of `id:tables_dir[:annotation]`; each
# experiment gets its own variant catalog `<id>_grch38`. A tables dir holding trans.parquet or gwas.parquet
# gets the experiment's trans and GWAS objects too. An experiment with no annotation named takes the one
# its adapter recorded (ingestion.json `source.gene_annotation`: gencode_v39 for the eQTL Catalogue),
# else gencode_v34. ANNOTATION_GTFS maps each annotation id to its GTF; each one an experiment uses is
# built once. The default is the TOPCHeF smoke build, which also runs the v0 comparison.
# QTLB_CHROMS=all builds every chromosome; give it more room: `sbatch --time=6:00:00 --mem=64G`.
# Every phase logs `phase <name>: <s> s wall, <KB> KB peak RSS` (/usr/bin/time).
set -uo pipefail
export PYTHONUNBUFFERED=1
export QTLB_CHROMS="${QTLB_CHROMS:-chr21,chr22}"
[ "$QTLB_CHROMS" = all ] && unset QTLB_CHROMS
STORE="${QTLB_STORE:-/scratch/ns5bc/qtl-browser/store-smoke}"
EXPERIMENTS="${EXPERIMENTS:-topchef:/scratch/ns5bc/qtl-browser/derived/_tables/topchef}"
cd /home/ns5bc/scratch/qtl-browser-live
echo "host $(hostname) chroms ${QTLB_CHROMS:-all} store $STORE experiments $EXPERIMENTS started $(date -Is)"
mkdir -p "$STORE"
ANNOTATION_GTFS="${ANNOTATION_GTFS:-gencode_v34:data/raw/gencode_v34/gencode.v34.annotation.gtf.gz gencode_v39:/project/shefflab/brickyard/datasets_downloaded/annotations/gencode.v39.annotation.gtf.gz}"
set -e
T() { local name="$1"; shift; /usr/bin/time -f "phase $name: %e s wall, %M KB peak RSS" "$@"; }
annotation() {  # the annotation id an experiment uses: named, else its adapter's, else gencode_v34
  [ -n "$2" ] && { echo "$2"; return; }
  uv run python -c "import json,sys; print((json.load(open(sys.argv[1])).get('source') or {}).get('gene_annotation') or 'gencode_v34')" \
    "$1/ingestion.json" 2>/dev/null || echo gencode_v34
}
build_annotation() {
  [ -f "$STORE/annotations/$1.json" ] && return
  local gtf=""
  for a in $ANNOTATION_GTFS; do [ "${a%%:*}" = "$1" ] && gtf="${a#*:}"; done
  [ -n "$gtf" ] || { echo "no GTF for annotation $1 in ANNOTATION_GTFS"; exit 1; }
  T "annotation $1" uv run python -m pipeline.annotation build --gtf "$gtf" --store "$STORE" --id "$1" | head -3
  echo "annotation $1 done $(date -Is)"
}
for e in $EXPERIMENTS; do
  IFS=: read -r EXP TABLES ANNOT <<< "$e"
  ANNOT="$(annotation "$TABLES" "$ANNOT")"
  build_annotation "$ANNOT"
  T "catalog $EXP" uv run python -m pipeline.catalog --sites "$TABLES/sites.parquet" --store "$STORE" --id "${EXP}_grch38"
  echo "catalog $EXP done $(date -Is)"
  T "results $EXP" uv run python -m pipeline.results build --tables "$TABLES" --store "$STORE" --id "$EXP" \
    --catalog "${EXP}_grch38" --annotation "$ANNOT" | grep -v '"source"'
  echo "results $EXP done $(date -Is)"
done
set +e
T validate uv run python - "$STORE" <<'EOF'
import sys
from pathlib import Path
from pipeline import catalog, qtlstore as qs
from pipeline.common import Config
st = qs.Store(sys.argv[1])
st.write_store(Path(sys.argv[1]).name, [])
fails = st.validate(refget=Config(), sites=catalog.read_sites)
print("validate:", "PASS" if not fails else f"{len(fails)} FAILURES", "notes:", st.notes)
for f in fails[:20]:
    print("  ", f)
# genes each experiment names that its annotation does not have (reported, not a failure)
from pipeline import annotation, results
for p in sorted((Path(sys.argv[1]) / "experiments").glob("*.json")):
    doc = st.load("experiments", p.stem)
    have = set(annotation.load(st, doc["annotation"])["genes"].column("gene_id").to_pylist())
    genes = {g for g in results.load_index(st, doc).column("gene_id").to_pylist() if g}
    miss = sorted(genes - have)
    print(f"genes {p.stem} vs {doc['annotation']}: {len(genes):,} named, {len(miss):,} missing", miss[:5])
EOF
# CROSSCAT="catalogA catalogB": the cross-catalog lookup check between two variant catalogs
[ -n "${CROSSCAT:-}" ] && T crosscat uv run python -m pipeline.qtlstore crosscat --store "$STORE" $CROSSCAT
[ -z "${SKIP_V0:-}" ] && T verify-v0 uv run python -m pipeline.verify_v0 --store "$STORE" --id topchef --genes "${GENES:-40}"
echo "finished $(date -Is)"
