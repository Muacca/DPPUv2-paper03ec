#!/usr/bin/env bash
# run_paper03ec.sh
# Run all paper03-EC scripts (proofs/ and paper03ec/).
# Each script writes its own timestamped log to the paper-level data/ directory.
#
# Usage (from the script/ directory):
#   bash run_paper03ec.sh                          # run all -> ../data/
#   bash run_paper03ec.sh proofs                   # run proofs/ only
#   bash run_paper03ec.sh paper03ec                # run paper03ec/ only
#   bash run_paper03ec.sh --output-dir my_logs     # custom output dir
#   bash run_paper03ec.sh -o my_logs paper03ec     # combined
#
# Output: <output-dir>/<script_name>_<timestamp>.log  (written by each script)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# ── argument parsing ──────────────────────────────────────────────────────────

TARGET="all"
DATA_DIR="${SCRIPT_DIR}/../data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir|-o)
            DATA_DIR="$2"
            shift 2
            ;;
        all|proofs|paper03ec)
            TARGET="$1"
            shift
            ;;
        *)
            echo "Usage: bash run_paper03ec.sh [--output-dir DIR] [all|proofs|paper03ec]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

export PYTHONIOENCODING=utf-8
export PYTHONUTF8=1
export DPPU_LOG_DIR="${DATA_DIR}"
export DPPU_LOG_STDOUT=0

# ── helper ────────────────────────────────────────────────────────────────────

FAILED_SCRIPTS=()

run_script() {
    local rel_path="$1"
    local script_file="${SCRIPTS_DIR}/${rel_path}"
    local base
    base="$(basename "${rel_path}" .py)"
    echo ">>> ${rel_path}"
    local t0
    t0=$(date +%s)

    if python "${script_file}"; then
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    OK  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
    else
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    FAILED  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
        FAILED_SCRIPTS+=("${rel_path}")
    fi
}

# ── script lists ──────────────────────────────────────────────────────────────

# Proof scripts: symbolic / algebraic proofs (no scipy numerics, fast)
PROOFS_SCRIPTS=(
    proofs/ax_vt_dropout_proof.py
    proofs/c_delta_null_proof.py
    proofs/t3_flatness_null_test.py
    proofs/palatini_protection_proof.py
    proofs/gamma_scaling_proof.py
)

# paper03ec scripts: mode dictionaries, EFT, comparisons (heavier computation)
PAPER03EC_SCRIPTS=(
    paper03ec/t3_mode_dictionary.py
    paper03ec/nil3_mode_dictionary.py
    paper03ec/s3_vt_spin_masses.py
    paper03ec/nil3_false_vacuum_eft.py
    paper03ec/nil3_spin2_quintet_splitting.py
    paper03ec/ec_cubic_vertex.py
    paper03ec/squash_shear_cross_term.py
    paper03ec/t3_mx_vacuum_search.py
)

# ── dispatch ──────────────────────────────────────────────────────────────────

echo "========================================"
echo "  run_paper03ec.sh  (target: ${TARGET})"
echo "  output dir: ${DATA_DIR}"
echo "========================================"
echo ""

if [[ "${TARGET}" == "all" || "${TARGET}" == "proofs" ]]; then
    echo "--- proofs/ ---"
    for s in "${PROOFS_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

if [[ "${TARGET}" == "all" || "${TARGET}" == "paper03ec" ]]; then
    echo "--- paper03ec/ ---"
    for s in "${PAPER03EC_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

# ── summary ───────────────────────────────────────────────────────────────────

echo "========================================"
if [[ ${#FAILED_SCRIPTS[@]} -eq 0 ]]; then
    echo "  All scripts completed successfully."
else
    echo "  ${#FAILED_SCRIPTS[@]} script(s) failed:"
    for s in "${FAILED_SCRIPTS[@]}"; do
        echo "    - ${s}"
    done
fi
echo "========================================"
