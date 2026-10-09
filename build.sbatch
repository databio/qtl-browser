#!/bin/bash
#SBATCH --job-name=qtlb_build
#SBATCH --output=/scratch/%u/qtlb_build_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#
# Build the shared input tables under `_tables/` (`python -m pipeline build`): the extracted Zenodo
# files, the variant table with rsIDs, the refget reference and allele check, and the v0 tables the
# TOPCHeF adapter and its gate read. The v1 store is built by adapter.sbatch and store.sbatch.
#
# `build` skips every step with a marker in <derived>/.done, so this is safe to resubmit after a
# failure or a timeout: it resumes at the first unfinished step. `nominal` re-encodes the 2.7 GB
# eQTL and 9.4 GB sQTL nominal archives with 3 workers at 6 GB each (config.yaml).
#
# QTLB_CHROMS and QTLB_DERIVED pass through for a subset build into its own tree; never point a
# subset build at the frozen v0 derived tree.
#
# Runs from the checkout QTLB_REPO (default /scratch/$USER/qtl-browser); with yoke that is the sync
# target, not a git checkout (see yoke.toml for why those are never the same directory). The log
# goes to /scratch/<user>/.

set -euo pipefail

cd "${QTLB_REPO:-/scratch/$USER/qtl-browser}"

echo "host $(hostname) chroms ${QTLB_CHROMS:-all} derived ${QTLB_DERIVED:-data/derived} started $(date -Is)"
uv run python -m pipeline build "$@"
echo "finished $(date -Is)"
