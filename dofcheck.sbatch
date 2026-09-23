#!/bin/bash
#SBATCH --job-name=qtlb_dof
#SBATCH --output=/home/ns5bc/scratch/qtlb_dof_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#
# `packcheck dof` then `packcheck report`.
#
# Check 3 of the report was failing for a missing input, not a bad measurement:
# it reports the dof search's result, the search lives in `packcheck dof`, and
# `dof` had never been run, so dof.json did not exist and the check reported
# "mode None, unanimous None over None sampled phenotypes". The evidence in the
# same row already showed the relation holds: 0 of 122,561,970 eQTL rows above
# 1e-05 relative error, max 5.96e-08, which is float32 rounding.
#
# `dof` rewrites packs.dof in config.yaml if the measured mode disagrees with it.
# It should not: 435 and 480 are what the data says. The report records
# config_before, config_after and config_updated, so read those. The sync is
# one-way from the laptop, so a change here never travels back.

set -uo pipefail

cd /home/ns5bc/scratch/qtl-browser-live

echo "host $(hostname)  started $(date -Is)"
cp pipeline/config.yaml /home/ns5bc/scratch/config_before_dof.yaml

uv run python -m pipeline packcheck dof
echo "dof exit $?  $(date -Is)"

diff -u /home/ns5bc/scratch/config_before_dof.yaml pipeline/config.yaml \
  && echo "config.yaml UNCHANGED" || echo "config.yaml CHANGED (see diff above)"

uv run python -m pipeline packcheck report
echo "report exit $?  $(date -Is)"
