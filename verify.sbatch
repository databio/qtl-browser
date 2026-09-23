#!/bin/bash
#SBATCH --job-name=qtlb_verify
#SBATCH --output=/home/ns5bc/scratch/qtlb_verify_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#
# `validate` and `packcheck report` on a compute node.
#
# Split out from build.sbatch because these two are the tail of a build that has
# already succeeded, and re-running the whole sequence to reach them wastes an
# hour. Every build step is marked done, so nothing here rebuilds anything.
#
# validate runs FIRST. It is what closes the refget anchor work, it does not
# depend on packcheck, and the login nodes have been sitting at load 50, so the
# cheapest way to get an answer is to not need the second command to succeed.

set -uo pipefail

cd /home/ns5bc/scratch/qtl-browser-live

echo "host $(hostname)  started $(date -Is)"

uv run python -m pipeline validate
rc=$?
echo "validate exit $rc  $(date -Is)"

uv run python -m pipeline packcheck report
echo "packcheck report exit $?  $(date -Is)"

exit $rc
