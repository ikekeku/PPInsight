"""
parse_pairs – extract protein interaction pairs from annotation tables.

Reads tables like the RTK Interactome spreadsheet and produces a flat
``pairs.csv`` / ``pairs.tsv`` with columns::

    proteinA, proteinB, label, family, references

where *label* is one of ``interaction`` or ``non-interaction``.

Usage::

    parse_pairs "RTK Interactome.tsv" -o pairs.csv
    parse_pairs table.tsv -o pairs.tsv --format tsv
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# Core parser
# ---------------------------------------------------------------------------

_REF_RE = re.compile(r"\s*\[[\d,\s]+\]\s*")


def _strip_refs(token: str) -> str:
    """Remove citation brackets like ``[6]`` or ``[5][8]``."""
    return _REF_RE.sub("", token).strip()


def _extract_refs(token: str) -> str:
    """Return just the citation numbers, e.g. '5,8' from 'EGFR [5][8]'."""
    refs = re.findall(r"\[(\d+)\]", token)
    return ",".join(refs)


def parse_interaction_table(path: str | os.PathLike) -> pd.DataFrame:
    """Parse an RTK-interactome-style TSV into a flat pairs DataFrame.

    The expected table layout (tab-separated)::

        <blank> ProteinA    Protein B   ...
        <blank> Family  Member  Interactions(5 cols)  Non-Interactions(4 cols)
        <blank> Ephrin  EPHA1   ERBB2 [6]  FGFR2 [7]  ...  <non-int cols>
        ...

    Continuation rows have empty Family/Member cells and continue
    listing partners for the same ProteinA.

    Returns a DataFrame with columns:
    ``proteinA, proteinB, label, family, references``
    """
    # Read raw TSV
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)

    # Find the header row that contains "Interactions" and "Non-Interactions".
    # The RTK interactome table has two header rows: one with broad
    # categories (ProteinA / ProteinB) and a sub-header with column roles
    # (Family / Member / Interactions / Non-Interactions).
    header_row_idx = None
    interaction_start = None
    non_interaction_start = None

    for idx, row in raw.iterrows():
        vals = [str(v).strip() for v in row.values]
        if "Interactions" in vals and "Non-Interactions" in vals:
            header_row_idx = idx
            interaction_start = vals.index("Interactions")
            non_interaction_start = vals.index("Non-Interactions")
            break

    if header_row_idx is None:
        # Fallback: try to infer from the first two rows
        # Row 0: ProteinA / Protein B
        # Row 1: Family / Member / Interactions / Non-Interactions
        header_row_idx = 1
        vals = [str(v).strip() for v in raw.iloc[1].values]
        # Find "Member" column as anchor
        if "Member" in vals:
            member_col = vals.index("Member")
            interaction_start = member_col + 1
            # The non-interaction columns start after the interaction columns
            # Count how many columns have data in interaction rows to find boundary
            # Default: 5 interaction columns, then non-interaction
            non_interaction_start = interaction_start + 5
        else:
            raise ValueError(
                "Cannot find 'Interactions'/'Non-Interactions' header row "
                f"in {path}. Is this the right file format?"
            )

    # Determine column indices
    # Family column and Member column come before interaction_start
    family_col = interaction_start - 2
    member_col = interaction_start - 1

    # How many interaction columns vs non-interaction columns?
    # interaction cols: from interaction_start to non_interaction_start - 1
    # non-interaction cols: from non_interaction_start to end
    int_cols = list(range(interaction_start, non_interaction_start))
    non_int_cols = list(range(non_interaction_start, len(raw.columns)))

    # Parse data rows (skip header rows).  Continuation rows (empty
    # Family/Member cells) inherit from the previous non-empty value,
    # which is how the RTK spreadsheet encodes multiple partners per protein.
    data_start = header_row_idx + 1
    rows: list[dict] = []
    current_family = ""
    current_member = ""

    for idx in range(data_start, len(raw)):
        row = raw.iloc[idx]
        fam = str(row.iloc[family_col]).strip() if family_col >= 0 else ""
        mem = str(row.iloc[member_col]).strip() if member_col >= 0 else ""

        if fam:
            current_family = fam
        if mem:
            current_member = mem

        if not current_member:
            continue

        # Parse interaction partners
        for col_idx in int_cols:
            if col_idx >= len(row):
                continue
            cell = str(row.iloc[col_idx]).strip()
            if not cell:
                continue
            protein_b = _strip_refs(cell)
            refs = _extract_refs(cell)
            if protein_b:
                rows.append({
                    "proteinA": current_member,
                    "proteinB": protein_b,
                    "label": "interaction",
                    "family": current_family,
                    "references": refs,
                })

        # Parse non-interaction partners
        for col_idx in non_int_cols:
            if col_idx >= len(row):
                continue
            cell = str(row.iloc[col_idx]).strip()
            if not cell:
                continue
            protein_b = _strip_refs(cell)
            refs = _extract_refs(cell)
            if protein_b:
                rows.append({
                    "proteinA": current_member,
                    "proteinB": protein_b,
                    "label": "non-interaction",
                    "family": current_family,
                    "references": refs,
                })

    if not rows:
        raise ValueError(f"No protein pairs extracted from {path}")

    df = pd.DataFrame(rows)
    # De-duplicate (same pair may appear on continuation rows)
    df = df.drop_duplicates(subset=["proteinA", "proteinB", "label"])
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI: parse an interaction table into a flat pairs file."""
    parser = argparse.ArgumentParser(
        description=(
            "Parse an RTK-interactome-style annotation table into a flat "
            "pairs file with interaction/non-interaction labels."
        ),
    )
    parser.add_argument(
        "table",
        help=(
            "Path to the TSV annotation table (RTK-interactome-style).  "
            "The parser extracts protein pairs and interaction/non-interaction "
            "labels from the table's matrix format."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default="pairs.csv",
        help=(
            "Output file (default: pairs.csv).  Extension determines format "
            "(.tsv → tab-separated, .csv → comma-separated).  This file is "
            "consumed by 'ppinsight batch' and 'ppinsight collect --pairs'."
        ),
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help=(
            "Print summary statistics (interaction/non-interaction counts, "
            "unique proteins) after writing.  Useful for verifying that the "
            "table was parsed correctly."
        ),
    )

    args = parser.parse_args(argv)

    df = parse_interaction_table(args.table)

    sep = "\t" if args.output.endswith(".tsv") else ","
    df.to_csv(args.output, sep=sep, index=False)
    print(f"Wrote {len(df)} pairs to {args.output}")

    if args.stats:
        n_int = (df["label"] == "interaction").sum()
        n_non = (df["label"] == "non-interaction").sum()
        n_proteins = len(set(df["proteinA"]) | set(df["proteinB"]))
        print(f"\n── Statistics ──")
        print(f"  Interactions:       {n_int}")
        print(f"  Non-interactions:   {n_non}")
        print(f"  Unique proteins:    {n_proteins}")
        print(f"  Unique proteinA:    {df['proteinA'].nunique()}")
        print(f"  Unique proteinB:    {df['proteinB'].nunique()}")
        print(f"  Families:           {df['family'].nunique()}")


if __name__ == "__main__":
    main()
