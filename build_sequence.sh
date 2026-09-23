#!/usr/bin/env bash
# The full qtlb build, in the one order that works. Called by build.sbatch (genome-wide) and
# smoke.sbatch (a two-chromosome subset); keep the order here so the two cannot drift.
#
# The order is not obvious and is not what packcheck's own error messages suggest:
#
#   1. `build` runs every step up to pack_gwas, then fails at `manifest`, which needs the packcheck
#      evidence that does not exist yet. That failure is expected on a fresh tree and is swallowed;
#      the step markers mean the second `build` picks up exactly at `manifest`.
#   2. `packcheck run` writes the per-chromosome evidence.
#   3. `packcheck roundtrip` needs only `run`. Its error message asks for `report` as well, but that
#      is true only in --steepest mode, where it reads phenotypes_<t>.parquet to choose phenotypes.
#      In genome scope it streams the raw files.
#   4. `build` again: roundtrip_genome.json now exists, so `manifest` completes and writes
#      manifest.json.
#   5. `packcheck report` reads manifest.json, so it can only run now. Putting it before the
#      manifest deadlocks the build.
#   6. `validate`.
#
# packcheck is not a build step and leaves no .done marker, so steps 2, 3 and 5 re-run every time.

set -uo pipefail

run() { echo "--- $* :: $(date -Is)"; "$@" || return $?; }

uv run python -m pipeline build || echo "expected: build stopped at manifest, continuing to packcheck"
echo "packs finished $(date -Is)"

set -e
run uv run python -m pipeline packcheck run
# `dof` must run before `report`: report's check 3 reports the dof SEARCH's result, and the
# search lives here. Skip it and dof.json is missing, so check 3 FAILs with "mode None" even
# though its own evidence shows the relation holding to 5.96e-08. That is a missing input, not
# a bad measurement, and it cost an afternoon of misreading.
run uv run python -m pipeline packcheck dof
run uv run python -m pipeline packcheck roundtrip
run uv run python -m pipeline build
run uv run python -m pipeline packcheck report
run uv run python -m pipeline validate
echo "validate finished $(date -Is)"
