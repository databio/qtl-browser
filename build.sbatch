#!/bin/bash
#SBATCH --job-name=qtlb_build
#SBATCH --output=/home/ns5bc/scratch/qtlb_build_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#
# Run the qtl-browser pipeline to completion on a compute node.
#
# Why a batch job: everything from `nominal` onward re-encodes the 2.7 GB eQTL
# and 9.4 GB sQTL nominal archives and then builds every pack. That is hours of
# CPU and tens of GB of output, which does not belong on a login node.
#
# `build` with no --step walks the whole step list and skips anything with a
# marker in data/derived/.done, so this is safe to resubmit after a failure or a
# timeout: it resumes at the first unfinished step. `--force` is deliberately
# not used.
#
# Memory: config.yaml runs `nominal` with 3 workers at 6 GB each, so 64 GB
# leaves real headroom for the pack steps' DuckDB limits on top.
#
# This runs in the yoke sync target, not in a git checkout. See yoke.toml for
# why those are never the same directory.

set -euo pipefail

cd /home/ns5bc/scratch/qtl-browser-live

echo "host $(hostname)  started $(date -Is)"

# The step order lives in build_sequence.sh so smoke.sbatch cannot drift from it.
./build_sequence.sh
