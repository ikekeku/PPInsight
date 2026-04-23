#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
# PPInsight — one-shot environment setup
#
# Usage:
#   bash setup.sh            # create env, install PyRosetta, install package
#   bash setup.sh --no-rosetta   # skip the PyRosetta download (lightweight)
# ─────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="ppinsight"

INSTALL_ROSETTA=true
for arg in "$@"; do
    case "$arg" in
        --no-rosetta) INSTALL_ROSETTA=false ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

# ── 1. conda environment ────────────────────────────────────────────
echo "▸ Creating conda environment '${ENV_NAME}' from ${ENV_FILE} …"
if conda env list | grep -qw "^${ENV_NAME} "; then
    echo "  (environment already exists — updating)"
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
    conda env create -f "${ENV_FILE}"
fi

# Resolve the environment's Python
PYBIN="$(conda run -n "${ENV_NAME}" which python)"
PIPBIN="$(dirname "${PYBIN}")/pip"
echo "  Python: ${PYBIN}"

# ── 2. PyRosetta (optional, ~1.5 GB download) ──────────────────────
if $INSTALL_ROSETTA; then
    echo "▸ Installing PyRosetta via pyrosetta-installer …"
    "${PYBIN}" -c "
import pyrosetta_installer
pyrosetta_installer.install_pyrosetta(skip_if_installed=True)
"
    echo "  ✓ PyRosetta installed"
else
    echo "▸ Skipping PyRosetta install (--no-rosetta)"
fi

# ── 3. Install the PPInsight package in editable mode ───────────────
echo "▸ Installing ppinsight in editable mode …"
"${PIPBIN}" install -e "${SCRIPT_DIR}" --quiet
echo "  ✓ ppinsight installed"

# ── 4. DockQ quality module (optional, needs git for fork install) ──
# DockQ is installed from a numpy-2-compatible fork (nrontsis/DockQ).
# This step can fail if git is unavailable or the fork URL changes.
# If it fails, everything except the ``ppinsight_quality`` command and
# the ``ppinsight.quality`` module will still work normally.
echo "▸ Installing DockQ (for quality assessment) …"
if "${PIPBIN}" install -e "${SCRIPT_DIR}[quality]" --quiet 2>/dev/null; then
    echo "  ✓ DockQ installed (ppinsight_quality CLI available)"
else
    echo "  ⚠ Could not install DockQ — the ppinsight_quality command"
    echo "    will not be available, but everything else works fine."
    echo "    To retry later:  pip install 'ppinsight[quality]'"
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
