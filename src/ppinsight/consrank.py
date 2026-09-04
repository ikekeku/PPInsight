"""
consrank.py
-----------
Reference-free consensus ranking of docking poses via Iter-CONSRANK
(https://github.com/AOCD-lab/Iter-consrank).

Unlike DockQ/CAPRI (``ppinsight quality``), CONSRANK does not need a native
structure: it ranks a pool of poses by contact-map consensus across the
pool itself, using a weighted average of contact conservation across each
pose's contacts. This lets poses from engines with incompatible native
score scales (LightDock energy, HADDOCK score, Rosetta I_sc) be compared
on common ground.

This module calls the compiled ``CONSRANK`` binary directly per iteration
and implements the iteration/cutoff loop in Python (using the reference
project's ``cut.f`` selection logic as a guide), rather than shelling out
to upstream's ``iter.sh``. This keeps orchestration, logging, and error
handling consistent with the other PPInsight engine wrappers.

The ``CONSRANK`` binary is not on ``$PATH`` — it is built by ``setup.sh``
from a pinned commit of the upstream repo vendored into
``third_party/iter_consrank/`` (see ``--no-iterconsrank`` in setup.sh).

CONSRANK needs the receptor/ligand chain IDs for a pose pool, but PPInsight's
three engines don't share a chain-ID convention:

* HADDOCK always normalizes staged inputs to a fixed chain A (receptor) /
  chain B (ligand) -- see ``_rewrite_pdb_as_single_chain`` in
  ``pdb_to_haddock.py``.
* Rosetta relabels chains to a sequential-letter split and records the
  split only in a pose-comment line dumped into the PDB
  (``ppinsight_partners ABC_DE`` -- see ``combine_proteins`` in
  ``rosetta/prepare_structure.py``).
* LightDock does not normalize or relabel at all -- receptor and ligand
  keep whatever chain IDs their original input PDBs had. The split is only
  recoverable from those original PDBs, which ``lightdock_pipeline`` copies
  unchanged into the top of the run directory.

``resolve_chains`` implements the auto-detection for all three; pass
``--rec-chains``/``--lig-chains`` to bypass it entirely.

To rank poses pooled *across* engines (the actual point of a reference-free
consensus ranker -- comparing engines whose native scores aren't on the same
scale), pass repeated ``--pool engine=run_dir`` instead of a single run_dir.
Every pool's poses are relabeled onto one shared receptor/ligand chain
scheme before CONSRANK ever sees them (see ``harmonize_pool_chains``), since
CONSRANK's CONTROL file has exactly one chain-ID namespace for the whole
pool. This only works when every pool has the same receptor chain *count*
and the same ligand chain *count* -- e.g. HADDOCK normalizes a multi-chain
receptor down to one chain while Rosetta/LightDock preserve the original
count, so pooling those together for a mixed-complex pair is refused with a
clear error rather than silently mis-comparing them.

Usage::

    python consrank.py data/output/rosetta_runs/2UUY_rec_vs_2UUY_lig --engine rosetta
    python consrank.py data/output/haddock_runs/run1 --engine haddock --cutoff 0.85
    python consrank.py --pool rosetta=data/output/rosetta_runs/run1 \\
                        --pool haddock=data/output/haddock_runs/run1
"""

import argparse
import glob
import gzip
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

import pandas as pd

from ppinsight import provenance
from ppinsight.utils import _project_root

_PDB_COORD_RECORDS = {"ATOM", "HETATM"}


class ConsrankError(RuntimeError):
    """Raised for CONSRANK setup or execution failures."""


# ---------------------------------------------------------------------------
# Binary resolution
# ---------------------------------------------------------------------------

def default_consrank_binary() -> str:
    """Return the expected path to the vendored ``CONSRANK`` binary.

    ``CONSRANK`` is never assumed to be on ``$PATH`` -- it is built by
    ``setup.sh`` under ``third_party/iter_consrank/`` (gitignored, vendored
    at a pinned commit).
    """
    return os.path.join(_project_root(), "third_party", "iter_consrank", "CONSRANK")


def consrank_binary_available(path: str | None = None) -> bool:
    """Return True if the CONSRANK binary exists and is executable."""
    path = path or default_consrank_binary()
    return os.path.isfile(path) and os.access(path, os.X_OK)


def _require_consrank_binary(path: str) -> str:
    if not consrank_binary_available(path):
        raise ConsrankError(
            f"CONSRANK binary not found or not executable at '{path}'.\n"
            "Build it with: bash setup.sh  (or rerun without --no-iterconsrank)\n"
            "It is vendored from https://github.com/AOCD-lab/Iter-consrank "
            "into third_party/iter_consrank/ and is never installed on $PATH."
        )
    return path


# ---------------------------------------------------------------------------
# CONTROL file (CONSRANK's input format)
# ---------------------------------------------------------------------------

def write_control(
    control_path: str,
    rec_chains: str,
    lig_chains: str,
    pdb_filenames: list[str],
    *,
    distance: float = 5.0,
    gen_mat: int = 0,
) -> None:
    """Write a CONSRANK ``CONTROL`` file.

    Mirrors the format produced by upstream's ``control-gen.py``:
    ``PairwiseChain1``/``PairwiseChain2`` list one chain letter per line,
    followed by the cutoff distance, the contact-matrix flag, and the pool
    of PDB filenames (relative to the CONSRANK working directory).
    """
    lines = [
        f"PairwiseChain1\t\t{len(rec_chains)}\t\t! # of chains in pairK",
        *list(rec_chains),
        f"PairwiseChain2\t\t{len(lig_chains)}\t\t! # of chains in pairL",
        *list(lig_chains),
        f"CUTOffDistance\t\t{distance}\t\t! Angstrom",
        f"GenMat\t\t\t{gen_mat}",
        f"NumberOfPDBFiles\t\t{len(pdb_filenames)}",
        *pdb_filenames,
    ]
    with open(control_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def _control_header_lines(
    rec_chains: str,
    lig_chains: str,
    *,
    distance: float,
    gen_mat: int,
) -> list[str]:
    """Return the chain/cutoff/genmat header lines shared by every iteration.

    Equivalent to upstream's ``control-gen-ith.py``, which copies every
    CONTROL line up to and including ``GenMat`` when building the next
    iteration's CONTROL file.
    """
    return [
        f"PairwiseChain1\t\t{len(rec_chains)}\t\t! # of chains in pairK",
        *list(rec_chains),
        f"PairwiseChain2\t\t{len(lig_chains)}\t\t! # of chains in pairL",
        *list(lig_chains),
        f"CUTOffDistance\t\t{distance}\t\t! Angstrom",
        f"GenMat\t\t\t{gen_mat}",
    ]


# ---------------------------------------------------------------------------
# Running CONSRANK and selecting the next iteration's pool
# ---------------------------------------------------------------------------

def run_command(cmd, cwd=None):
    """Print and execute a CONSRANK command, raising ConsrankError on failure."""
    print(">>", " ".join(map(str, cmd)))
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True,
        )
    except FileNotFoundError as exc:
        raise ConsrankError(f"'{cmd[0]}' could not be executed: {exc}") from exc

    if proc.stdout:
        print(proc.stdout, end="")
    if proc.returncode != 0:
        raise ConsrankError(
            f"CONSRANK exited with code {proc.returncode}.\n{proc.stderr}"
        )
    return proc.stdout


def _parse_consrank_scores(scores_path: str) -> list[tuple[str, float]]:
    """Parse ``Consrank_score.txt`` into ``(filename, score)`` pairs.

    Each line is whitespace-separated ``<filename> <score> ...``; only the
    first two fields are used. Lines that don't parse as ``name, float``
    are skipped rather than raising, since trailing blank lines are common.
    """
    rows: list[tuple[str, float]] = []
    with open(scores_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                score = float(parts[1])
            except ValueError:
                continue
            rows.append((parts[0], score))
    if not rows:
        raise ConsrankError(f"No parseable rows in {scores_path}")
    return rows


def _select_top_fraction(
    scored: list[tuple[str, float]],
    cutoff: float,
) -> list[tuple[str, float]]:
    """Keep the top *cutoff* fraction of poses by CONSRANK score.

    A higher CONSRANK score means a pose shares more consensus contacts
    with the pool, so this keeps the highest-scoring tail after an
    ascending numeric sort -- the same selection cut.f's percentage/count
    logic performs. (Upstream's iter.sh instead does a lexicographic
    ``sort -k 2`` before invoking cut.f's compiled selector; this reimplements
    the clearly-intended numeric ordering rather than reproducing that quirk.)
    """
    ordered = sorted(scored, key=lambda row: row[1])
    n = len(ordered)
    n_keep = max(1, min(n, round(cutoff * n)))
    return ordered[n - n_keep:]


def run_iterative_consrank(
    pose_dir: str,
    pdb_filenames: list[str],
    rec_chains: str,
    lig_chains: str,
    *,
    distance: float = 5.0,
    gen_mat: int = 0,
    cutoff: float = 0.85,
    iterations: int = 10,
    binary: str | None = None,
) -> dict:
    """Run the CONSRANK iteration/cutoff loop over poses in *pose_dir*.

    Reimplements upstream's ``iter.sh`` orchestration in Python: writes a
    CONTROL file, runs the CONSRANK binary, keeps the top *cutoff* fraction
    of poses by score (``cut.f``'s selection logic), and repeats for
    *iterations* rounds (stopping early if the pool can no longer shrink).

    Parameters
    ----------
    pose_dir : str
        Working directory containing the pose PDB files. CONSRANK reads
        and writes its CONTROL/score files here.
    pdb_filenames : list[str]
        Pose filenames (relative to *pose_dir*) to rank.
    rec_chains, lig_chains : str
        Chain-ID strings, e.g. ``"C"`` and ``"BA"``.
    binary : str | None
        Path to the CONSRANK executable (default: the vendored build under
        third_party/iter_consrank/).

    Returns
    -------
    dict
        ``{"final_scores": [(filename, score), ...], "history": [...]}``
        where each history entry records ``iteration``, ``n_input``, and
        ``n_kept``.
    """
    binary = _require_consrank_binary(binary or default_consrank_binary())
    binary = os.path.abspath(binary)

    if len(pdb_filenames) < 2:
        raise ConsrankError(
            f"CONSRANK needs at least 2 poses to rank; got {len(pdb_filenames)}."
        )

    control_path = os.path.join(pose_dir, "CONTROL")
    scores_path = os.path.join(pose_dir, "Consrank_score.txt")

    header = _control_header_lines(
        rec_chains, lig_chains, distance=distance, gen_mat=gen_mat,
    )
    current = sorted(pdb_filenames)
    history: list[dict] = []
    final_scores: list[tuple[str, float]] = []

    for i in range(1, iterations + 1):
        write_control(
            control_path, rec_chains, lig_chains, current,
            distance=distance, gen_mat=gen_mat,
        )
        run_command([binary, "-c", "CONTROL"], cwd=pose_dir)

        if not os.path.isfile(scores_path):
            raise ConsrankError(
                f"CONSRANK did not produce {scores_path} (iteration {i})."
            )
        scored = _parse_consrank_scores(scores_path)
        final_scores = scored

        kept = _select_top_fraction(scored, cutoff)
        n_input, n_kept = len(scored), len(kept)
        history.append({"iteration": i, "n_input": n_input, "n_kept": n_kept})
        print(f"  [consrank] iteration {i}: {n_input} -> {n_kept} poses")

        # Preserve the per-iteration Consrank_score.txt / CONTROL for
        # traceability, mirroring iter.sh's file-renaming step.
        shutil.copy(scores_path, os.path.join(pose_dir, f"Consrank_score-{i}.txt"))
        shutil.copy(control_path, os.path.join(pose_dir, f"CONTROL-{i}.txt"))

        if n_kept == n_input or n_kept <= 1:
            # Pool stopped shrinking (or is down to a single pose) -- further
            # iterations would be a no-op, so stop early rather than looping
            # to the requested count for nothing.
            break
        current = sorted(name for name, _ in kept)

    return {
        "final_scores": final_scores,
        "header": header,
        "history": history,
        "control_path": control_path,
        "scores_path": scores_path,
    }


# ---------------------------------------------------------------------------
# Pose discovery per engine
# ---------------------------------------------------------------------------

def discover_poses(
    run_dir: str, engine: str, *, max_poses: int | None = None,
) -> list[str]:
    """Return absolute paths to candidate pose PDBs for *engine*'s *run_dir*.

    Each engine lays out poses differently, so discovery is engine-specific:

    * ``lightdock`` -- generated conformations in every ``swarm_*/`` dir.
    * ``rosetta`` -- ``decoy_*.pdb`` files written by ``save_all_decoys``.
    * ``haddock`` -- clustered models under ``*_seletopclusts/`` (gzipped).
    """
    if engine == "lightdock":
        pattern = os.path.join(run_dir, "swarm_*", "lightdock_*.pdb")
        found = sorted(glob.glob(pattern))
    elif engine == "rosetta":
        found = sorted(glob.glob(os.path.join(run_dir, "decoy_*.pdb")))
    elif engine == "haddock":
        pattern = os.path.join(run_dir, "**", "*_seletopclusts", "cluster_*.pdb.gz")
        found = sorted(glob.glob(pattern, recursive=True))
    else:
        raise ConsrankError(f"Unsupported engine for pose discovery: '{engine}'")

    if not found:
        raise ConsrankError(
            f"No {engine} poses found under '{run_dir}'. "
            "Pass --engine explicitly if auto-detection picked the wrong engine."
        )

    if max_poses is not None and len(found) > max_poses:
        print(
            f"  [consrank] {len(found)} poses found, capping to --max-poses {max_poses}"
        )
        found = found[:max_poses]

    return [os.path.abspath(p) for p in found]


def stage_poses(
    pose_paths: list[str], pose_dir: str, *, prefix: str = "",
) -> list[str]:
    """Copy (and gunzip, if needed) *pose_paths* into *pose_dir*.

    Returns the staged basenames (what CONSRANK will see via ``os.listdir``).
    CONSRANK only reads plain ``.pdb`` files in its working directory, so
    gzipped HADDOCK models are decompressed here rather than symlinked.

    Pass *prefix* (e.g. an engine name) when staging several pools into one
    shared directory, so filenames stay traceable to their source even when
    two engines happen to use the same naming scheme.
    """
    os.makedirs(pose_dir, exist_ok=True)
    staged: list[str] = []
    seen: set[str] = set()
    for src in pose_paths:
        base = os.path.basename(src)
        if base.endswith(".gz"):
            base = base[: -len(".gz")]
        if prefix:
            base = f"{prefix}_{base}"
        # Disambiguate collisions (e.g. same decoy name from different
        # swarms) by prefixing with the parent directory name.
        if base in seen:
            base = f"{os.path.basename(os.path.dirname(src))}_{base}"
        seen.add(base)

        dst = os.path.join(pose_dir, base)
        if src.endswith(".gz"):
            with gzip.open(src, "rb") as fh_in, open(dst, "wb") as fh_out:
                shutil.copyfileobj(fh_in, fh_out)
        else:
            shutil.copy(src, dst)
        staged.append(base)
    return staged


# ---------------------------------------------------------------------------
# Chain-ID collision handling (LightDock only -- see resolve_chains)
# ---------------------------------------------------------------------------

_FRESH_CHAIN_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _relabel_pose_chains_by_block_order(
    pose_path: str,
    chain_map: list[str],
) -> None:
    """Rewrite *pose_path* in place, relabeling contiguous chain-ID blocks.

    LightDock does not renumber chains to avoid collisions: when the
    receptor and ligand's own input PDBs reuse the same letters (a common
    case -- e.g. both start their own chain lettering at 'A'), the merged
    pose keeps those letters as-is. Verified directly against a real run
    (``P00648_vs_P11540``): the pose is 5 *contiguous* blocks in file order
    -- ``A(878) B(871) C(872) A(718) B(718)`` -- receptor's 3 original
    chains, then ligand's 2, with 'A'/'B' reused across the boundary.

    So the reused letters can't distinguish receptor from ligand, but the
    block *order* still can: the Nth contiguous same-letter run in file
    order is relabeled to ``chain_map[N]``, regardless of what letter it
    originally had. Callers build ``chain_map`` so the first blocks (however
    many the receptor has) map to receptor letters and the rest to ligand
    letters -- either fresh disjoint letters for one pool
    (``relabel_colliding_lightdock_chains``), or one shared scheme across
    several pools (``harmonize_pool_chains``).
    """
    total_chains = len(chain_map)

    lines_out = []
    block_index = -1
    prev_chain = None
    with open(pose_path, encoding="utf-8") as fh:
        for line in fh:
            record = line[:6].strip()
            if record in _PDB_COORD_RECORDS:
                chain = line[21]
                if chain != prev_chain:
                    block_index += 1
                    prev_chain = chain
                if block_index >= total_chains:
                    raise ConsrankError(
                        f"'{pose_path}' has more contiguous chain blocks than "
                        f"the expected {total_chains}; refusing to guess a "
                        "relabeling. Pass --rec-chains/--lig-chains explicitly."
                    )
                line = line[:21] + chain_map[block_index] + line[22:]
            elif record == "TER":
                prev_chain = None
            lines_out.append(line)

    with open(pose_path, "w", encoding="utf-8") as fh:
        fh.writelines(lines_out)


def relabel_colliding_lightdock_chains(
    pose_dir: str,
    staged_filenames: list[str],
    rec_chains: str,
    lig_chains: str,
) -> tuple[str, str]:
    """Relabel every staged pose when receptor/ligand chain IDs collide.

    Returns the new, guaranteed-disjoint ``(rec_chains, lig_chains)`` to use
    in CONTROL. If the two sets are already disjoint, returns them unchanged
    and touches no files.
    """
    if not set(rec_chains) & set(lig_chains):
        return rec_chains, lig_chains

    total = len(rec_chains) + len(lig_chains)
    if total > len(_FRESH_CHAIN_ALPHABET):
        raise ConsrankError(
            f"Cannot relabel {total} chains -- only "
            f"{len(_FRESH_CHAIN_ALPHABET)} letters available."
        )

    print(
        f"  [consrank] receptor/ligand chain IDs collide "
        f"(rec={rec_chains}, lig={lig_chains}) -- relabeling "
        f"{len(staged_filenames)} staged poses to disjoint letters."
    )
    chain_map = list(_FRESH_CHAIN_ALPHABET[:total])
    for filename in staged_filenames:
        _relabel_pose_chains_by_block_order(os.path.join(pose_dir, filename), chain_map)
    new_rec = _FRESH_CHAIN_ALPHABET[: len(rec_chains)]
    new_lig = _FRESH_CHAIN_ALPHABET[len(rec_chains): total]
    return new_rec, new_lig


# ---------------------------------------------------------------------------
# Cross-engine pooling
# ---------------------------------------------------------------------------

@dataclass
class PoseSource:
    """One engine's contribution to a pooled CONSRANK ranking."""

    engine: str
    run_dir: str
    pose_paths: list[str]
    rec_chains: str
    lig_chains: str
    staged_filenames: list[str] = field(default_factory=list)


def parse_pool_spec(spec: str) -> tuple[str, str]:
    """Parse an ``ENGINE=RUN_DIR`` ``--pool`` argument."""
    engine, sep, path = spec.partition("=")
    if not sep or not engine or not path:
        raise ConsrankError(f"Invalid --pool value '{spec}'; expected ENGINE=RUN_DIR.")
    if engine not in _CHAIN_RESOLVERS:
        raise ConsrankError(
            f"Unknown engine '{engine}' in --pool '{spec}'. "
            f"Known engines: {', '.join(sorted(_CHAIN_RESOLVERS))}."
        )
    return engine, path


def harmonize_pool_chains(
    pose_dir: str,
    sources: list[PoseSource],
) -> tuple[str, str]:
    """Relabel every staged pose across *sources* onto one shared chain scheme.

    CONSRANK's CONTROL file has exactly one ``PairwiseChain1``/
    ``PairwiseChain2`` definition for the whole pool, so every pose --
    regardless of which engine produced it -- must end up using the same
    chain letters for "receptor" and the same letters for "ligand". This
    only works when every source has the same receptor chain *count* and
    the same ligand chain *count*: chain topology (how many chains the
    receptor was split into) isn't something a shared letter scheme can
    paper over, so a mismatch is a real incompatibility, not a labeling
    detail -- e.g. HADDOCK normalizes a multi-chain receptor down to a
    single chain while Rosetta/LightDock preserve the original chain count,
    so a mixed-complex pair docked by both isn't poolable.

    Assumes ``source.staged_filenames`` has already been populated by
    ``stage_poses``. Requires each *source*'s poses to be laid out as
    contiguous per-chain blocks in file order (receptor's chains, then the
    ligand's) -- true for all three current engines' pose outputs.
    """
    rec_counts = {len(s.rec_chains) for s in sources}
    lig_counts = {len(s.lig_chains) for s in sources}
    if len(rec_counts) > 1 or len(lig_counts) > 1:
        detail = "; ".join(
            f"{s.engine}: receptor={s.rec_chains!r} ({len(s.rec_chains)} chains), "
            f"ligand={s.lig_chains!r} ({len(s.lig_chains)} chains)"
            for s in sources
        )
        raise ConsrankError(
            "Cannot pool poses with mismatched chain topology across "
            f"engines: {detail}. This commonly happens because HADDOCK "
            "normalizes a multi-chain partner down to a single chain while "
            "LightDock/Rosetta preserve the original chain count -- these "
            "poses aren't structurally comparable at the chain level. Pool "
            "only engines whose runs share the same receptor/ligand chain "
            "counts for this pair."
        )

    n_rec, n_lig = rec_counts.pop(), lig_counts.pop()
    total = n_rec + n_lig
    if total > len(_FRESH_CHAIN_ALPHABET):
        raise ConsrankError(
            f"Cannot harmonize {total} chains -- only "
            f"{len(_FRESH_CHAIN_ALPHABET)} letters available."
        )

    chain_map = list(_FRESH_CHAIN_ALPHABET[:total])
    for source in sources:
        print(
            f"  [consrank] relabeling {len(source.staged_filenames)} "
            f"{source.engine} poses onto the shared chain scheme."
        )
        for filename in source.staged_filenames:
            _relabel_pose_chains_by_block_order(
                os.path.join(pose_dir, filename), chain_map,
            )

    return _FRESH_CHAIN_ALPHABET[:n_rec], _FRESH_CHAIN_ALPHABET[n_rec:total]


# ---------------------------------------------------------------------------
# Chain-ID resolution (see module docstring -- no shared convention exists)
# ---------------------------------------------------------------------------

def _pdb_chain_order(pdb_path: str) -> list[str]:
    """Return chain IDs in *pdb_path*, in order of first appearance."""
    order: list[str] = []
    seen: set[str] = set()
    with open(pdb_path, encoding="utf-8") as fh:
        for line in fh:
            if line[:6].strip() not in _PDB_COORD_RECORDS:
                continue
            chain = line[21].strip()
            if chain and chain not in seen:
                seen.add(chain)
                order.append(chain)
    return order


def _resolve_chains_haddock(_run_dir: str, _pose_paths: list[str]) -> tuple[str, str]:
    """HADDOCK always normalizes staged inputs to chain A (rec) / B (lig)."""
    return "A", "B"


def _resolve_chains_rosetta(_run_dir: str, pose_paths: list[str]) -> tuple[str, str]:
    """Read the ``ppinsight_partners`` pose-comment Rosetta dumps into the PDB.

    Written by ``combine_proteins`` in ``rosetta/prepare_structure.py`` as a
    bare ``ppinsight_partners <rec_chains>_<lig_chains>`` line between
    ``##Begin comments##``/``##End comments##`` markers -- not a standard
    PDB record, so a generic parser won't find it; we look for it explicitly.
    """
    for pose_path in pose_paths:
        with open(pose_path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip().startswith("ppinsight_partners"):
                    parts = line.split()
                    if len(parts) >= 2 and "_" in parts[1]:
                        rec, _, lig = parts[1].partition("_")
                        if rec and lig:
                            return rec, lig
    raise ConsrankError(
        "Could not find a 'ppinsight_partners' comment in any Rosetta decoy "
        "under this run directory to recover the receptor/ligand chain "
        "split. Pass --rec-chains/--lig-chains explicitly."
    )


def _resolve_chains_lightdock(run_dir: str, _pose_paths: list[str]) -> tuple[str, str]:
    """Recover the split from the original inputs LightDock copies unchanged.

    ``lightdock_pipeline`` copies the receptor/ligand PDBs into the top of
    the run directory before docking (see ``pdb_to_lightdock.py``), and
    LightDock does not relabel or merge chain IDs -- so those two files'
    own chain sets *are* the receptor/ligand split in every generated pose.
    """
    top_level_pdbs = [
        p for p in glob.glob(os.path.join(run_dir, "*.pdb"))
        if not os.path.basename(p).startswith("lightdock_")
    ]
    if len(top_level_pdbs) != 2:
        raise ConsrankError(
            f"Expected exactly 2 original input PDBs at the top of '{run_dir}' "
            f"(found {len(top_level_pdbs)}); cannot recover the receptor/"
            "ligand chain split automatically. Pass --rec-chains/--lig-chains."
        )
    top_level_pdbs.sort(key=os.path.getmtime)
    rec_path, lig_path = top_level_pdbs
    rec_chains = "".join(_pdb_chain_order(rec_path))
    lig_chains = "".join(_pdb_chain_order(lig_path))
    if not rec_chains or not lig_chains:
        raise ConsrankError(
            f"Could not read chain IDs from '{rec_path}' / '{lig_path}'. "
            "Pass --rec-chains/--lig-chains explicitly."
        )
    return rec_chains, lig_chains


_CHAIN_RESOLVERS = {
    "haddock": _resolve_chains_haddock,
    "rosetta": _resolve_chains_rosetta,
    "lightdock": _resolve_chains_lightdock,
}


def resolve_chains(
    run_dir: str,
    engine: str,
    pose_paths: list[str],
    *,
    rec_chains: str | None = None,
    lig_chains: str | None = None,
) -> tuple[str, str]:
    """Resolve receptor/ligand chain IDs for a pose pool.

    Explicit ``rec_chains``/``lig_chains`` always win. Otherwise dispatches
    to the engine-specific resolver in ``_CHAIN_RESOLVERS`` -- see the
    module docstring for why there is no single cross-engine convention.
    """
    if rec_chains and lig_chains:
        return rec_chains, lig_chains
    if rec_chains or lig_chains:
        raise ConsrankError(
            "--rec-chains and --lig-chains must be given together."
        )

    resolver = _CHAIN_RESOLVERS.get(engine)
    if resolver is None:
        raise ConsrankError(
            f"No chain-ID resolver for engine '{engine}'. "
            "Pass --rec-chains/--lig-chains explicitly."
        )
    return resolver(run_dir, pose_paths)


# ---------------------------------------------------------------------------
# Output directory + unified scores
# ---------------------------------------------------------------------------

def make_consrank_output_dir(
    label: str,
    base_root: str = "data/output",
    method: str = "consrank_runs",
) -> str:
    """Create an auto-incrementing consrank working directory.

    Mirrors ``make_output_dir`` in ``pdb_to_lightdock.py``.
    """
    if not os.path.isabs(base_root):
        base_root = os.path.join(_project_root(), base_root)
    method_dir = os.path.join(base_root, method)
    os.makedirs(method_dir, exist_ok=True)

    run_dir = os.path.join(method_dir, label)
    if os.path.exists(run_dir):
        idx = 1
        while True:
            candidate = os.path.join(method_dir, f"{label}_{idx}")
            if not os.path.exists(candidate):
                run_dir = candidate
                break
            idx += 1
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def scores_to_dataframe(
    final_scores: list[tuple[str, float]],
    pose_dir: str,
    *,
    pair: tuple[str, str] | None = None,
    label: str = "consrank",
    pose_engine_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build a unified scores DataFrame matching collect_scores.py's schema.

    ``pose_engine_map`` (filename -> engine name) lets a pooled ranking's
    ``model`` column show which engine each ranked pose actually came from,
    instead of the fixed *label* used for a single-engine ranking.
    """
    pA, pB = pair or ("", "")
    ordered = sorted(final_scores, key=lambda row: row[1], reverse=True)
    rows = []
    for rank, (filename, score) in enumerate(ordered, start=1):
        output_path = os.path.abspath(os.path.join(pose_dir, filename))
        model = pose_engine_map.get(filename, label) if pose_engine_map else label
        rows.append({
            "model": model,
            "score_type": "consrank_score",
            "score_value": score,
            "proteinA": pA,
            "proteinB": pB,
            "pose_id": os.path.splitext(filename)[0],
            "output_path": output_path if os.path.isfile(output_path) else "",
            "source_file": os.path.join(pose_dir, "Consrank_score.txt"),
            "pose_rank": rank,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI entrypoint: reference-free consensus ranking of pooled poses."""
    parser = argparse.ArgumentParser(
        description=(
            "Rank a pool of docking poses by contact-map consensus "
            "(Iter-CONSRANK) -- no native reference structure required."
        ),
    )
    parser.add_argument(
        "run_dir",
        nargs="?",
        default=None,
        help=(
            "Docking engine run directory containing poses to rank. "
            "Omit this and use --pool instead to rank poses pooled across "
            "several engines."
        ),
    )
    parser.add_argument(
        "--pool",
        action="append",
        metavar="ENGINE=RUN_DIR",
        default=None,
        help=(
            "Add one engine's poses to a cross-engine pool, e.g. "
            "'--pool rosetta=data/output/rosetta_runs/run1'. Repeat for "
            "each engine to include. Mutually exclusive with run_dir/"
            "--engine/--rec-chains/--lig-chains -- each pool auto-resolves "
            "its own chains, then every pose is relabeled onto one shared "
            "scheme (see the module docstring)."
        ),
    )
    parser.add_argument(
        "--engine",
        choices=["lightdock", "haddock", "rosetta"],
        default=None,
        help=(
            "Which engine produced run_dir (default: auto-detect using "
            "PPInsight's registry detectors, the same ones 'ppinsight "
            "collect' uses). Not used with --pool."
        ),
    )
    parser.add_argument(
        "--rec-chains", default=None,
        help=(
            "Receptor chain IDs, e.g. 'C'. Overrides auto-detection -- use "
            "this when the per-engine chain resolver can't be applied "
            "(see the module docstring for why chain IDs aren't uniform "
            "across engines)."
        ),
    )
    parser.add_argument(
        "--lig-chains", default=None,
        help="Ligand chain IDs, e.g. 'BA'. Must be given together with --rec-chains.",
    )
    parser.add_argument(
        "--distance", type=float, default=5.0,
        help="Contact-distance cutoff in Angstrom (default: 5.0, upstream's default).",
    )
    parser.add_argument(
        "--gen-mat", action="store_true",
        help="Also generate the full contacts matrix (CONSRANK's GenMat=1).",
    )
    parser.add_argument(
        "--cutoff", type=float, default=0.85,
        help=(
            "Fraction of poses kept each iteration (default: 0.85, matching "
            "upstream's default 'cut' file). Lower values converge faster "
            "to a smaller consensus set; 1.0 disables pruning between "
            "iterations."
        ),
    )
    parser.add_argument(
        "--iterations", type=int, default=10,
        help=(
            "Maximum number of refinement iterations (default: 10). "
            "Stops early once the pool can no longer shrink."
        ),
    )
    parser.add_argument(
        "--max-poses", type=int, default=None,
        help=(
            "Cap the number of poses pulled from run_dir before ranking "
            "(default: no cap). LightDock pools in particular can include "
            "hundreds of generated conformations across all swarms."
        ),
    )
    parser.add_argument(
        "--consrank-bin", default=None,
        help=(
            "Path to the CONSRANK binary "
            "(default: third_party/iter_consrank/CONSRANK)."
        ),
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help=(
            "Output scores TSV/CSV (default: <run_dir's consrank_runs "
            "working dir>/consrank_scores.tsv)."
        ),
    )

    args = parser.parse_args(argv)

    if bool(args.pool) == bool(args.run_dir):
        parser.error(
            "Pass exactly one of: run_dir, or one or more --pool ENGINE=RUN_DIR."
        )
    if args.pool and (args.engine or args.rec_chains or args.lig_chains):
        parser.error(
            "--pool cannot be combined with --engine/--rec-chains/--lig-chains."
        )

    try:
        consrank_bin = _require_consrank_binary(
            args.consrank_bin or default_consrank_binary()
        )
        if args.pool:
            (
                pose_dir, rec_chains, lig_chains, staged,
                pair_label, pose_engine_map, source_dirs,
            ) = _prepare_pool(args)
        else:
            (
                pose_dir, rec_chains, lig_chains, staged,
                pair_label, pose_engine_map, source_dirs,
            ) = _prepare_single(args)

        result = run_iterative_consrank(
            pose_dir, staged, rec_chains, lig_chains,
            distance=args.distance,
            gen_mat=1 if args.gen_mat else 0,
            cutoff=args.cutoff,
            iterations=args.iterations,
            binary=consrank_bin,
        )
    except ConsrankError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    pair = tuple(pair_label.split("_vs_")) if "_vs_" in pair_label else None
    scores_df = scores_to_dataframe(
        result["final_scores"], pose_dir, pair=pair, pose_engine_map=pose_engine_map,
    )

    output_path = args.output or os.path.join(pose_dir, "consrank_scores.tsv")
    sep = "\t" if output_path.endswith(".tsv") else ","
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    scores_df.to_csv(output_path, sep=sep, index=False)

    prov = provenance.extract_run_metadata("consrank", pose_dir)
    prov["engine_meta"].update({
        "rec_chains": rec_chains,
        "lig_chains": lig_chains,
        "distance": args.distance,
        "cutoff": args.cutoff,
        "history": result["history"],
        "consrank_binary": consrank_bin,
        "source_dirs": source_dirs,
    })
    run_id = provenance.make_run_id("consrank", pair)
    provenance.write_sidecar(output_path, {run_id: prov})

    print(f"\nConsensus ranking complete: {len(scores_df)} poses ranked.")
    print(f"Scores written to: {output_path}")
    print(f"Provenance sidecar: {provenance.sidecar_path(output_path)}")
    print(f"Working directory (CONTROL/score history): {pose_dir}")


def _prepare_single(args):
    """Resolve poses/chains for a plain single-run_dir CLI invocation."""
    run_dir = os.path.abspath(args.run_dir)
    if not os.path.isdir(run_dir):
        print(f"ERROR: run_dir not found: {run_dir}", file=sys.stderr)
        sys.exit(2)

    engine = args.engine
    if engine is None:
        from ppinsight import registry
        try:
            engine = registry.detect_engine(run_dir)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            print("Hint: pass --engine explicitly.", file=sys.stderr)
            sys.exit(2)
        print(f"Detected engine: {engine}")

    pose_paths = discover_poses(run_dir, engine, max_poses=args.max_poses)
    rec_chains, lig_chains = resolve_chains(
        run_dir, engine, pose_paths,
        rec_chains=args.rec_chains, lig_chains=args.lig_chains,
    )
    print(f"Receptor chains: {rec_chains}   Ligand chains: {lig_chains}")
    print(f"Pose pool: {len(pose_paths)} poses from {run_dir}")

    pair_label = os.path.basename(run_dir.rstrip(os.sep))
    pose_dir = make_consrank_output_dir(pair_label)
    staged = stage_poses(pose_paths, pose_dir)

    if engine == "lightdock":
        # LightDock reuses the receptor's/ligand's own chain letters
        # verbatim, so they collide whenever both inputs happen to label
        # chains the same way (e.g. both start at 'A') -- see
        # relabel_colliding_lightdock_chains's docstring. Other engines
        # are disjoint by construction, so this is skipped for them.
        rec_chains, lig_chains = relabel_colliding_lightdock_chains(
            pose_dir, staged, rec_chains, lig_chains,
        )

    return pose_dir, rec_chains, lig_chains, staged, pair_label, None, [run_dir]


def _prepare_pool(args):
    """Resolve poses/chains for a --pool cross-engine CLI invocation."""
    sources: list[PoseSource] = []
    for spec in args.pool:
        engine, path = parse_pool_spec(spec)
        run_dir = os.path.abspath(path)
        if not os.path.isdir(run_dir):
            print(f"ERROR: run_dir not found: {run_dir}", file=sys.stderr)
            sys.exit(2)
        pose_paths = discover_poses(run_dir, engine, max_poses=args.max_poses)
        rec_chains, lig_chains = resolve_chains(run_dir, engine, pose_paths)
        print(
            f"[{engine}] {len(pose_paths)} poses from {run_dir} "
            f"(receptor={rec_chains}, ligand={lig_chains})"
        )
        sources.append(PoseSource(engine, run_dir, pose_paths, rec_chains, lig_chains))

    pair_labels = sorted({os.path.basename(s.run_dir.rstrip(os.sep)) for s in sources})
    pair_label = pair_labels[0]
    if len(pair_labels) > 1:
        print(
            f"  [consrank] Warning: pooled run directories have different "
            f"names ({pair_labels}); using '{pair_label}' for output labeling."
        )

    pose_dir = make_consrank_output_dir(pair_label)
    for source in sources:
        source.staged_filenames = stage_poses(
            source.pose_paths, pose_dir, prefix=source.engine,
        )

    rec_chains, lig_chains = harmonize_pool_chains(pose_dir, sources)
    print(f"Harmonized chains across pool: receptor={rec_chains}   ligand={lig_chains}")

    staged = [f for s in sources for f in s.staged_filenames]
    pose_engine_map = {f: s.engine for s in sources for f in s.staged_filenames}
    source_dirs = [s.run_dir for s in sources]

    return (
        pose_dir, rec_chains, lig_chains, staged,
        pair_label, pose_engine_map, source_dirs,
    )


if __name__ == "__main__":
    main()
