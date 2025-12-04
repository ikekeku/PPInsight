#!/usr/bin/env bash
# Helper script to run haddock on Hyak-like clusters (SLURM)
# Usage: ./hyak_run_haddock.sh /path/to/run_dir run1.cfg

set -euo pipefail

RUN_DIR=${1:-.}
CFG=${2:-run1.cfg}

if [ ! -d "$RUN_DIR" ]; then
  echo "Run directory not found: $RUN_DIR" >&2
  exit 2
fi

cd "$RUN_DIR"

# Check cns if haddock is installed in the same python env
if python -c "import sys
try:
  import haddock, pathlib
  p = pathlib.Path(haddock.__file__).resolve().parent / 'bin' / 'cns'
  print(p)
except Exception as e:
  sys.exit(0)
" >/dev/null 2>&1; then
  CNS_PATH=$(python -c "import haddock, pathlib; print(pathlib.Path(haddock.__file__).resolve().parent / 'bin' / 'cns')")
  if [ -f "$CNS_PATH" ]; then
    echo "CNS binary info:"
    file "$CNS_PATH" || true
  fi
fi

echo "Running haddock: haddock3 $CFG"
haddock3 "$CFG"
