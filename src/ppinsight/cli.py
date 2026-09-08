"""
Umbrella CLI for PPInsight.

Dispatches to sub-command modules via lazy import so that only the
requested command's dependencies are loaded (avoids pulling in PyRosetta
or LightDock for unrelated sub-commands).

Usage::

    ppinsight fetch P69905 P68871
    ppinsight fetch-native P35968 P15692 -o native_candidates.tsv
    ppinsight lightdock 2UUY_rec 2UUY_lig --steps 100
    ppinsight haddock   2UUY_rec 2UUY_lig --run
    ppinsight rosetta   2UUY_rec 2UUY_lig --n-runs 5
    ppinsight collect   examples/haddock3/run1-test -o scores.tsv
    ppinsight compare   scores.tsv --metric dockq
    ppinsight parse     table.tsv -o pairs.csv
    ppinsight batch     pairs.csv --engines lightdock haddock
    ppinsight quality   model.pdb native.pdb
    ppinsight quality   scores.tsv native.pdb -o scores_quality.tsv
    ppinsight consrank  data/output/rosetta_runs/run1 --engine rosetta
"""

import sys

_SUBCOMMANDS = {
    "fetch":     "ppinsight.protein_fetch",
    "fetch-native": "ppinsight.fetch_native",
    "lightdock": "ppinsight.pdb_to_lightdock",
    "haddock":   "ppinsight.pdb_to_haddock",
    "rosetta":   "ppinsight.pdb_to_rosetta",
    "collect":   "ppinsight.collect_scores",
    "compare":   "ppinsight.visualizer",
    "parse":     "ppinsight.parse_pairs",
    "batch":     "ppinsight.batch_dock",
    "purge":     "ppinsight.purge_runs",
    "quality":   "ppinsight.quality",
    "consrank":  "ppinsight.consrank",
}


def main():
    """Dispatch to the appropriate sub-command's ``main()``."""
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        _print_help()
        sys.exit(0)

    cmd = sys.argv[1]
    if cmd not in _SUBCOMMANDS:
        print(f"Unknown subcommand '{cmd}'.\n", file=sys.stderr)
        _print_help()
        sys.exit(2)

    # Lazy-import the module so we only load what's needed.
    # This keeps `ppinsight fetch` fast even when PyRosetta is installed
    # (PyRosetta's import alone takes several seconds).
    import importlib
    mod = importlib.import_module(_SUBCOMMANDS[cmd])

    # Rewrite sys.argv so the sub-command sees its own argument list.
    # Without this, argparse in the sub-module would choke on "ppinsight"
    # as an unrecognised argument.
    sys.argv = [f"ppinsight {cmd}"] + sys.argv[2:]
    mod.main()


def _print_help():
    print("ppinsight – unified protein-protein docking toolkit\n")
    print("Usage: ppinsight <command> [options]\n")
    print("Commands:")
    pad = max(len(k) for k in _SUBCOMMANDS) + 2
    descriptions = {
        "fetch":     "Fetch protein data from UniProt / PDB",
        "fetch-native": "Find native experimental complexes in RCSB",
        "lightdock": "Run LightDock docking pipeline",
        "haddock":   "Stage / run a HADDOCK3 docking pipeline",
        "rosetta":   "Run PyRosetta docking pipeline",
        "collect":   "Collect docking scores into unified TSV/CSV",
        "compare":   "Visualise and compare docking scores",
        "parse":     "Parse an interaction table into a flat pairs file",
        "batch":     "Batch-run docking for all pairs in a pairs file",
        "purge":     "Preview or remove failed batch-run directories",
        "quality":   "Evaluate docking quality with DockQ (CAPRI metrics)",
        "consrank":  "Reference-free consensus ranking (Iter-CONSRANK)",
    }
    for cmd in _SUBCOMMANDS:
        desc = descriptions.get(cmd, "")
        print(f"  {cmd:<{pad}} {desc}")
    print("\nRun 'ppinsight <command> --help' for command-specific options.")


if __name__ == "__main__":
    main()
