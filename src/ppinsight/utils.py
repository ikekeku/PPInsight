"""Shared utilities used by all PPInsight pipeline scripts."""

import glob
import os

_PDB_COORD_RECORDS = {"ATOM", "HETATM", "ANISOU"}


def _project_root() -> str:
    """Return the absolute path to the PPInsight repository root."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def resolve_input_path(path: str, search_root: str | None = None) -> str:
    """Resolve an input path by returning it if it exists or searching the repo.

    Accepts short names like ``'2UUY_rec'`` or ``'2UUY_rec.pdb'`` and returns
    the absolute path of the first matching file found under *search_root*
    (defaults to the repository root).

    Parameters
    ----------
    path : str
        Filename, basename, or full path to a PDB file.
    search_root : str | None
        Directory to search when *path* is not an existing file.
        Defaults to the PPInsight repository root.

    Returns
    -------
    str
        Absolute path to the resolved file.

    Raises
    ------
    FileNotFoundError
        If the file cannot be found.
    """
    if not path:
        raise FileNotFoundError("Empty input path")

    p = os.path.expanduser(path)
    p = os.path.abspath(p)
    if os.path.exists(p):
        return p

    # Search repo for basename with common structure-file extensions.
    # This allows users to pass short names like "2UUY_rec" from the CLI
    # and have them resolved against the repo's example / input files.
    base = os.path.basename(path)
    stem, ext = os.path.splitext(base)
    candidates = [base]
    if not ext:
        candidates.extend([base + '.pdb', base + '.ent'])
    elif ext.lower() == '.pdb':
        candidates.append(stem + '.ent')
    elif ext.lower() == '.ent':
        candidates.append(stem + '.pdb')

    seen = set()
    candidates = [c for c in candidates if not (c in seen or seen.add(c))]

    proj = os.path.abspath(search_root) if search_root else _project_root()
    for c in candidates:
        pattern = os.path.join(proj, '**', c)
        matches = glob.glob(pattern, recursive=True)
        if matches:
            found = os.path.abspath(matches[0])
            print(f"Resolved '{path}' -> '{found}'")
            return found

    raise FileNotFoundError(
        f"Could not find input file '{path}' (searched {proj})."
    )


def dbref_chains_for_accession(
    pdb_path: str | os.PathLike[str],
    accession: str | None,
) -> list[str]:
    """Return coordinate chain IDs whose DBREF mapping matches *accession*."""
    if not accession:
        return []

    accession = accession.upper()
    matched_chains = []
    seen = set()

    with open(pdb_path, encoding="utf-8") as handle:
        for line in handle:
            record = line[:6].strip()
            tokens = line.split()
            chain = None
            mapped_accession = None

            if record == "DBREF" and len(tokens) >= 7:
                chain = tokens[2].strip()
                mapped_accession = tokens[6].strip().upper()
            elif record == "DBREF2" and len(tokens) >= 4:
                chain = tokens[2].strip()
                accession_index = 4 if len(tokens) >= 5 else 3
                mapped_accession = tokens[accession_index].strip().upper()

            if not chain or mapped_accession != accession or chain in seen:
                continue

            seen.add(chain)
            matched_chains.append(chain)

    return matched_chains


def copy_pdb_selected_chains(
    input_pdb: str | os.PathLike[str],
    output_pdb: str | os.PathLike[str],
    allowed_chains: set[str],
) -> None:
    """Copy a PDB while keeping only coordinate records for *allowed_chains*."""
    previous_coord_was_written = False

    with open(input_pdb, encoding="utf-8") as src, open(
        output_pdb,
        "w",
        encoding="utf-8",
    ) as dst:
        for line in src:
            record = line[:6].strip()

            if record in _PDB_COORD_RECORDS:
                chain = line[21].strip()
                if chain not in allowed_chains:
                    previous_coord_was_written = False
                    continue
                dst.write(line)
                previous_coord_was_written = True
                continue

            if record == "TER":
                if previous_coord_was_written:
                    dst.write(line)
                previous_coord_was_written = False
                continue

            dst.write(line)
