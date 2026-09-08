"""Purge safely recorded failed PPInsight batch-run directories."""

import argparse
import os
from pathlib import Path

import pandas as pd

from ppinsight.batch_dock import _is_safe_output_dir


def failed_run_candidates(
    manifest_path: str | Path,
    output_root: str | Path,
) -> pd.DataFrame:
    """Return failed manifest rows with removable directories under the root."""
    manifest_path = Path(manifest_path)
    separator = "\t" if manifest_path.suffix == ".tsv" else ","
    results = pd.read_csv(manifest_path, sep=separator)
    required = {"proteinA", "proteinB", "engine", "output_dir", "status"}
    if not required.issubset(results.columns):
        raise ValueError(
            "Results manifest is missing required columns: "
            + ", ".join(sorted(required - set(results.columns)))
        )

    failed = results[results["status"] == "failed"].copy()
    failed = failed[failed["output_dir"].notna()]
    failed = failed[failed["output_dir"].astype(str).str.len() > 0]
    failed["safe_to_remove"] = failed["output_dir"].map(
        lambda path: _is_safe_output_dir(str(path), str(output_root))
    )
    return failed.drop_duplicates(
        ["proteinA", "proteinB", "engine"], keep="last"
    )


def purge_failed_runs(
    manifest_path: str | Path,
    output_root: str | Path,
    *,
    execute: bool = False,
) -> tuple[int, int]:
    """Report or remove failed run directories recorded in a batch manifest."""
    candidates = failed_run_candidates(manifest_path, output_root)
    removable = candidates[candidates["safe_to_remove"]]
    skipped = len(candidates) - len(removable)

    removed = 0
    for row in removable.itertuples(index=False):
        path = Path(row.output_dir)
        action = "Removing" if execute else "Would remove"
        print(f"{action}: {path}")
        if execute:
            if path.is_dir():
                import shutil

                shutil.rmtree(path)
                removed += 1
        else:
            removed += 1

    for row in candidates[~candidates["safe_to_remove"]].itertuples(index=False):
        print(f"Skipping unsafe path outside output root: {row.output_dir}")
    return removed, skipped


def main(argv=None):
    """CLI entrypoint for manifest-driven failed-run cleanup."""
    parser = argparse.ArgumentParser(
        description=(
            "List or purge failed PPInsight run directories recorded in a "
            "batch-results CSV/TSV manifest."
        )
    )
    parser.add_argument("manifest", help="Batch results CSV/TSV to inspect.")
    parser.add_argument(
        "--output-root",
        default="data/output",
        help="Only remove run paths beneath this root (default: data/output).",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Actually delete paths. Without this flag the command is a dry run.",
    )
    args = parser.parse_args(argv)

    removed, skipped = purge_failed_runs(
        args.manifest,
        os.path.abspath(os.path.expanduser(args.output_root)),
        execute=args.yes,
    )
    action = "Removed" if args.yes else "Found"
    print(
        f"{action} {removed} failed run directory(s); "
        f"skipped {skipped} unsafe path(s)."
    )


if __name__ == "__main__":
    main()
