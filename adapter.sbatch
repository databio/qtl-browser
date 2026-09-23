#!/bin/bash
#SBATCH --job-name=qtlb_adapter
#SBATCH --output=/home/ns5bc/scratch/qtlb_adapter_%j.log
#SBATCH --partition=standard
#SBATCH --account=shefflab
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#
# One adapter into its own derived tree, then its checks. Never the validated v0 tree.
#   ADAPTER     topchef | eqtl_catalogue | dcm_gwas                   (default topchef)
#               dcm_gwas writes gwas.parquet + gwas.json into the topchef tables (`--tables DIR` for others)
#   QTLB_DERIVED  the derived tree (required; refuses the validated v0 tree); the TOPCHeF
#               tree symlinks v0's `_tables` inputs, `immutable/` and `_tmp/pack_pointers`
#   STEPS       any of: adapter gate contract                         (default: all three)
#               gate = the TOPCHeF acceptance gate (verify_topchef), skipped for other adapters;
#               contract = pipeline/adapters/contract_check.py on the tables the adapter wrote
#   QTLB_CHROMS comma list, or unset / `all` for every chromosome
# Examples:
#   QTLB_DERIVED=/scratch/ns5bc/qtl-browser/derived-topchef-v2 sbatch adapter.sbatch
#   ADAPTER=eqtl_catalogue QTLB_DERIVED=/scratch/ns5bc/qtl-browser/derived-eqtlcat-all-v2 sbatch adapter.sbatch
#   STEPS=contract ADAPTER=eqtl_catalogue QTLB_DERIVED=... sbatch --mem=32G adapter.sbatch
#   STEPS=gate QTLB_CHROMS=chr21,chr22 QTLB_DERIVED=... sbatch adapter.sbatch
#   STEPS="adapter contract" QTLB_DERIVED=... sbatch adapter.sbatch --step trans    # one TOPCHeF table
#   ADAPTER=dcm_gwas STEPS="adapter contract" QTLB_DERIVED=... sbatch adapter.sbatch
# Extra arguments go to the adapter (e.g. `--experiment ID` for eqtl_catalogue, `--force` for topchef).
set -uo pipefail
export PYTHONUNBUFFERED=1
ADAPTER="${ADAPTER:-topchef}"
STEPS="${STEPS:-adapter gate contract}"
[ "${QTLB_CHROMS:-}" = all ] && unset QTLB_CHROMS
: "${QTLB_DERIVED:?set QTLB_DERIVED to the tree to write (never the validated v0 tree)}"
case "${QTLB_DERIVED%/}" in */qtl-browser/derived) echo "refusing the validated v0 tree"; exit 2;; esac
case "$ADAPTER" in
  topchef|dcm_gwas) EXP_ID=topchef ;;
  eqtl_catalogue) EXP_ID=gtex_v8_heart_lv
    for ((i = 1; i <= $#; i++)); do [ "${!i}" = --experiment ] && { j=$((i + 1)); EXP_ID="${!j}"; }; done ;;
  *) echo "unknown ADAPTER $ADAPTER"; exit 2 ;;
esac
TABLES="$QTLB_DERIVED/_tables/$EXP_ID"
cd /home/ns5bc/scratch/qtl-browser-live
mkdir -p "$QTLB_DERIVED"
echo "host $(hostname) adapter $ADAPTER steps '$STEPS' chroms ${QTLB_CHROMS:-all} derived $QTLB_DERIVED started $(date -Is)"
quiet() { grep --line-buffered -v "DeprecationWarning\|fetch_arrow_table()"; }
rc=0
if [[ " $STEPS " == *" adapter "* ]]; then
  uv run python -m "pipeline.adapters.$ADAPTER" "$@" 2>&1 | quiet
  rc=${PIPESTATUS[0]}; echo "adapter exit $rc $(date -Is)"; [ $rc -eq 0 ] || exit $rc
fi
if [[ " $STEPS " == *" gate "* && $ADAPTER == topchef ]]; then
  args=()
  [ -n "${QTLB_CHROMS:-}" ] && args=(--chrom ${QTLB_CHROMS//,/ })
  uv run python -m pipeline.adapters.verify_topchef --genes "${GENES:-20}" "${args[@]}" 2>&1 | quiet
  g=${PIPESTATUS[0]}; echo "gate exit $g $(date -Is)"; [ $g -eq 0 ] || rc=$g
fi
if [[ " $STEPS " == *" contract "* ]]; then
  uv run python -m pipeline.adapters.contract_check "$TABLES" 2>&1 | quiet
  k=${PIPESTATUS[0]}; echo "contract exit $k $(date -Is)"; [ $k -eq 0 ] || rc=$k
fi
echo "finished rc $rc $(date -Is)"
exit $rc
