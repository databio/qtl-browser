#!/bin/bash
#SBATCH --job-name=qtlb_smoke
#SBATCH --output=/home/ns5bc/scratch/qtlb_smoke_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#
# The whole pipeline over two small chromosomes, for finding breakage fast.
#
# QTLB_DERIVED is what makes this safe to run while a real build is going: every output, step
# marker included, hangs off `derived`, so this tree is entirely its own. Never point a subset
# build at the real data/derived.
#
# This is a smoke test, not a release. `validate` will fail its paper_counts checks (egenes 10,241,
# sQTL phenotypes 13,540) because only a genome-wide build can meet them, and the chromosome
# ordinals in the variant and rsID indexes are positions in the subset list. Read the failures for
# what broke structurally, not for those two.

set -uo pipefail

export QTLB_CHROMS="${QTLB_CHROMS:-chr21,chr22}"
export QTLB_DERIVED="${QTLB_DERIVED:-/scratch/ns5bc/qtl-browser/derived-smoke}"

cd /home/ns5bc/scratch/qtl-browser-live
mkdir -p "$QTLB_DERIVED"

echo "host $(hostname)  chroms $QTLB_CHROMS  derived $QTLB_DERIVED  started $(date -Is)"
./build_sequence.sh
echo "smoke finished $(date -Is)"
