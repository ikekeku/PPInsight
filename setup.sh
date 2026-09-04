#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
# PPInsight — one-shot environment setup
#
# Usage:
#   bash setup.sh            # create env, install PyRosetta, install package
#   bash setup.sh --no-rosetta   # skip the PyRosetta download (lightweight)
#   bash setup.sh --no-iterconsrank   # skip building the Iter-CONSRANK binary
# ─────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="ppinsight"
ITERCONSRANK_REPO="https://github.com/AOCD-lab/Iter-consrank.git"
ITERCONSRANK_COMMIT="943244394b8679e7d35cf697d7c13e3e0e058094"
ITERCONSRANK_DIR="${SCRIPT_DIR}/third_party/iter_consrank"

INSTALL_ROSETTA=true
INSTALL_ITERCONSRANK=true
for arg in "$@"; do
    case "$arg" in
        --no-rosetta) INSTALL_ROSETTA=false ;;
        --no-iterconsrank) INSTALL_ITERCONSRANK=false ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

# ── 1. conda environment ────────────────────────────────────────────
echo "▸ Creating conda environment '${ENV_NAME}' from ${ENV_FILE} …"
export CFLAGS="${CFLAGS:-} -Wno-error=incompatible-pointer-types"
if conda env list | grep -qw "^${ENV_NAME} "; then
    echo "  (environment already exists — updating)"
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
    conda env create -f "${ENV_FILE}"
fi

# Resolve the environment's Python
PYBIN="$(conda run -n "${ENV_NAME}" which python)"
PIPBIN="$(dirname "${PYBIN}")/pip"
ENV_BIN_DIR="$(dirname "${PYBIN}")"
echo "  Python: ${PYBIN}"

run_in_env_context() {
    env \
        -u PYTHONPATH \
        -u PIP_USER \
        PATH="${ENV_BIN_DIR}:${PATH}" \
        PYTHONNOUSERSITE=1 \
        "$@"
}

# ── 2. PyRosetta (optional, ~1.5 GB download) ──────────────────────
if $INSTALL_ROSETTA; then
    echo "▸ Installing PyRosetta via pyrosetta-installer …"
    run_in_env_context "${PYBIN}" -c "
import pyrosetta_installer
pyrosetta_installer.install_pyrosetta(skip_if_installed=True, mirror=1)
"
    echo "  ✓ PyRosetta installed"
else
    echo "▸ Skipping PyRosetta install (--no-rosetta)"
fi

# ── 3. Install the PPInsight package in editable mode ───────────────
echo "▸ Installing ppinsight in editable mode …"
run_in_env_context "${PIPBIN}" install -e "${SCRIPT_DIR}" --quiet
echo "  ✓ ppinsight installed"

# ── 4. DockQ quality module (optional, needs git for fork install) ──
# DockQ is installed from a numpy-2-compatible fork (nrontsis/DockQ).
# This step can fail if git is unavailable or the fork URL changes.
# If it fails, everything except the ``ppinsight_quality`` command and
# the ``ppinsight.quality`` module will still work normally.
echo "▸ Installing DockQ (for quality assessment) …"
if run_in_env_context "${PIPBIN}" install -e "${SCRIPT_DIR}[quality]" --quiet 2>/dev/null; then
    echo "  ✓ DockQ installed (ppinsight_quality CLI available)"
else
    echo "  ⚠ Could not install DockQ — the ppinsight_quality command"
    echo "    will not be available, but everything else works fine."
    echo "    To retry later:  pip install 'ppinsight[quality]'"
fi

# ── 5. Iter-CONSRANK (optional, reference-free consensus ranking) ──
# Iter-CONSRANK isn't a pip/conda package — it's vendored via git clone
# into third_party/iter_consrank/, pinned to a fixed commit for
# reproducibility.  Only the CONSRANK binary is built here; the
# iteration/cutoff loop that upstream's iter.sh + cut.f implement is
# reimplemented in ppinsight.consrank so orchestration, logging, and
# error handling stay consistent with the rest of the pipeline.
# This step can fail if git/a C++ compiler is unavailable.  If it fails,
# everything except the ``ppinsight consrank`` command still works fine.
if $INSTALL_ITERCONSRANK; then
    echo "▸ Building Iter-CONSRANK (for 'ppinsight consrank') …"
    iterconsrank_ok=true

    if [ -d "${ITERCONSRANK_DIR}/.git" ]; then
        git -C "${ITERCONSRANK_DIR}" fetch --quiet origin "${ITERCONSRANK_COMMIT}" \
            2>/dev/null || iterconsrank_ok=false
    else
        rm -rf "${ITERCONSRANK_DIR}"
        mkdir -p "$(dirname "${ITERCONSRANK_DIR}")"
        git clone --quiet "${ITERCONSRANK_REPO}" "${ITERCONSRANK_DIR}" \
            2>/dev/null || iterconsrank_ok=false
    fi

    if $iterconsrank_ok; then
        git -C "${ITERCONSRANK_DIR}" checkout --quiet "${ITERCONSRANK_COMMIT}" \
            2>/dev/null || iterconsrank_ok=false
    fi

    if $iterconsrank_ok; then
        # Use `conda run -n` (not run_in_env_context) so the cxx-compiler /
        # c-compiler / fortran-compiler activation hooks set CC/CXX/FC —
        # run_in_env_context only prepends PATH and would miss those.
        if (
            cd "${ITERCONSRANK_DIR}/src" \
            && conda run -n "${ENV_NAME}" bash -c '
                set -euo pipefail
                for f in *.cpp; do "${CXX:-g++}" -c "$f"; done
                "${CXX:-g++}" -o CONSRANK *.o
                cp CONSRANK ..
            '
        ); then
            echo "  ✓ Iter-CONSRANK built (ppinsight consrank CLI available)"
        else
            iterconsrank_ok=false
        fi
    fi

    if ! $iterconsrank_ok; then
        echo "  ⚠ Could not build Iter-CONSRANK — the ppinsight consrank"
        echo "    command will not be available, but everything else works fine."
        echo "    To retry later:  rerun bash setup.sh (or build manually in"
        echo "    third_party/iter_consrank/, see its README)"
    fi
else
    echo "▸ Skipping Iter-CONSRANK build (--no-iterconsrank)"
fi

# ── done ─────────────────────────────────────────────────────────────
echo ""
echo "Setup complete!  Activate the environment with:"
echo ""
echo "    conda activate ${ENV_NAME}"
echo ""
echo "Then try:"
echo "    ppinsight rosetta 2UUY_rec 2UUY_lig --n-runs 1"
echo "    ppinsight haddock 2UUY_rec 2UUY_lig --runname test"
echo "    ppinsight lightdock 2UUY_rec 2UUY_lig --generate"
echo "    ppinsight consrank data/output/lightdock_runs/2UUY_rec_vs_2UUY_lig --engine lightdock"
