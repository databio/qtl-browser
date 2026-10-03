#!/bin/bash
#SBATCH --job-name=qtlb_bench
#SBATCH --output=/home/ns5bc/scratch/qtlb_bench_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#
# v0 vs v1 format benchmark (pipeline/bench_store.py): bytes per table type per chromosome, and
# random-access block read speed. Reads only; writes only $OUT.
#   sbatch bench_store.sbatch                                  # chr21/22 two-study store
#   QTLB_STORE=/scratch/ns5bc/qtl-browser/store-genome-v1f sbatch bench_store.sbatch
#   EXPERIMENT=gtex_v8_heart_lv QTLB_STORE=... sbatch bench_store.sbatch --no-reads   # no v0 twin: sizes only
# CHROMS: all (the catalog's chromosomes) or a comma list. N random phenotypes per type, REPS per read.
# Extra arguments go to pipeline.bench_store. OUT defaults to bench/<store>[_<experiment>]_<date>.
set -euo pipefail
export PYTHONUNBUFFERED=1
V0="${QTLB_V0:-/scratch/ns5bc/qtl-browser/derived/immutable}"
STORE="${QTLB_STORE:-/scratch/ns5bc/qtl-browser/store-two-sqtlfix}"
EXP="${EXPERIMENT:-topchef}"
SUF=""; [ "$EXP" = topchef ] || SUF="_$EXP"
OUT="${OUT:-/scratch/ns5bc/qtl-browser/bench/$(basename "$STORE")${SUF}_$(date +%F)}"
cd /home/ns5bc/scratch/qtl-browser-live
echo "host $(hostname) v0 $V0 store $STORE experiment $EXP out $OUT started $(date -Is)"
uv run python -m pipeline.bench_store --v0 "$V0" --store "$STORE" --experiment "$EXP" \
  --chroms "${CHROMS:-all}" --n "${N:-200}" --reps "${REPS:-5}" --out "$OUT" "$@"
echo "finished $(date -Is)"
