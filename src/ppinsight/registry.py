"""
registry – central plugin registry for docking engines.

Every engine is represented by an :class:`EnginePlugin` dataclass that bundles
three callables:

* **runner** – execute a docking job given receptor + ligand PDB paths.
* **parser** – extract scores from a completed run directory.
* **detector** – probe a directory and decide whether this engine produced it.

The module provides a simple :func:`register` / :func:`get` /
:func:`list_engines` API.  The three built-in engines (LightDock, HADDOCK,
Rosetta) are registered lazily the first time the registry is queried, so
importing this module does **not** pull in heavy optional dependencies like
PyRosetta.

Third-party or user engines can be added at runtime::

    from ppinsight.registry import register, EnginePlugin

    register(EnginePlugin(
        name="my_docker",
        runner=my_run_fn,
        parser=my_parse_fn,
        detector=my_detect_fn,
    ))

After that, ``batch_dock --engines my_docker`` and
``collect_scores --engine my_docker`` will "just work".
"""

from __future__ import annotations

import glob
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    import pandas as pd


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

class RunnerFunc(Protocol):
    """Signature for an engine runner callable."""
    def __call__(
        self,
        rec_pdb: str,
        lig_pdb: str,
        output_root: str,
        pair_label: str,
        **kwargs,
    ) -> str | None: ...


class ParserFunc(Protocol):
    """Signature for an engine score-parser callable."""
    def __call__(
        self,
        run_dir: str,
        pair: tuple[str, str] | None = None,
        label: str = ...,
    ) -> "pd.DataFrame": ...


class DetectorFunc(Protocol):
    """Signature for a directory detector callable."""
    def __call__(self, directory: str) -> bool: ...


class EngineRunError(RuntimeError):
    """A failed engine run with structured diagnostic context for batch mode."""

    def __init__(
        self,
        message: str,
        *,
        output_dir: str | None = None,
        log_path: str | None = None,
        error_type: str | None = None,
    ):
        super().__init__(message)
        self.output_dir = output_dir
        self.log_path = log_path
        self.error_type = error_type or type(self).__name__


@dataclass(slots=True)
class EnginePlugin:
    """Encapsulates a docking-engine integration.

    Parameters
    ----------
    name : str
        Short lowercase identifier (e.g. ``"lightdock"``).
    runner : RunnerFunc | None
        Callable that executes the engine.  May be *None* for parse-only
        engines (useful for analysing existing outputs without being able to
        re-run them).
    parser : ParserFunc | None
        Callable that extracts a unified scores DataFrame from an output
        directory.  May be *None* for run-only integrations.
    detector : DetectorFunc | None
        Callable that returns *True* if *directory* was produced by this
        engine.  Used by :func:`detect_engine`.
    description : str
        One-liner shown when listing available engines.
    """
    name: str
    runner: RunnerFunc | None = None
    parser: ParserFunc | None = None
    detector: DetectorFunc | None = None
    description: str = ""


# ---------------------------------------------------------------------------
# Registry storage
# ---------------------------------------------------------------------------

_registry: dict[str, EnginePlugin] = {}
_defaults_loaded: bool = False
_LIGHTDOCK_ANM_MISMATCH_RE = re.compile(
    r"\[ANM\]\s*ERROR:.*Number of atoms in Prody.*LightDock",
    flags=re.IGNORECASE | re.DOTALL,
)


def _is_lightdock_anm_atom_mismatch(output: str | None) -> bool:
    """Return True for known LightDock ANM setup atom-count mismatch errors."""
    if not output:
        return False
    return bool(_LIGHTDOCK_ANM_MISMATCH_RE.search(output))


def _is_interactive_session() -> bool:
    """Return True when stdin/stdout are attached to an interactive terminal."""
    return sys.stdin.isatty() and sys.stdout.isatty()


def _confirm_lightdock_retry_without_anm(rec_pdb: str, lig_pdb: str) -> bool:
    """Prompt before retrying LightDock with ANM disabled."""
    print(
        "  [lightdock] Detected ANM setup atom mismatch for pair "
        f"{os.path.basename(rec_pdb)} vs {os.path.basename(lig_pdb)}."
    )
    print(
        "  [lightdock] Retry this pair with ANM disabled (rigid-body fallback)?"
    )
    try:
        reply = input("  [lightdock] Retry without ANM? [y/N]: ")
    except EOFError:
        return False
    return reply.strip().lower() in {"y", "yes"}


def _run_rosetta_batch(rec_pdb, lig_pdb, output_root, pair_label, **kw):
    """Default Rosetta runner used by batch mode."""
    from ppinsight.pdb_to_rosetta import (
        _make_output_dir,
        _rosetta_outputs_support_clustering,
    )
    from ppinsight.rosetta.pipeline import DockingPipeline

    output_dir = None
    try:
        output_dir = kw.get("output_dir") or _make_output_dir(
            rec_pdb,
            lig_pdb,
            base_root=output_root,
            method="rosetta_runs",
        )
        os.makedirs(output_dir, exist_ok=True)
        verbose = bool(kw.get("verbose", False))

        pipeline = DockingPipeline(
            rec_pdb,
            lig_pdb,
            n_runs=int(kw.get("n_runs", 10)),
            top_n=int(kw.get("top_n", 20)),
            relax=bool(kw.get("relax", True)),
            verbose=verbose,
            auto_filter=bool(kw.get("auto_filter", True)),
            pyrosetta_debug=bool(kw.get("pyrosetta_debug", False)),
        )
        pipeline.run()

        save_top = int(kw.get("save_top", 0))
        if save_top > 0:
            pipeline.save_top_structures(output_dir, top_n=save_top)

        csv_path = os.path.join(output_dir, "docking_scores.csv")
        pipeline.save_scores(csv_path)
        pipeline.save_all_decoys(output_dir)

        if kw.get("cluster", True):
            cluster_ready, _ = _rosetta_outputs_support_clustering(
                csv_path,
                output_dir,
            )
            if cluster_ready:
                from ppinsight.rosetta.analyze import cluster_and_rank

                cluster_and_rank(
                    csv_path,
                    output_dir,
                    score_col="i_sc",
                    top_n=int(kw.get("cluster_top_n", 200)),
                    rmsd_cutoff=float(kw.get("rmsd_cutoff", 4.0)),
                )

        return str(output_dir)
    except Exception as exc:
        raise EngineRunError(
            str(exc),
            output_dir=str(output_dir) if output_dir else None,
            error_type=type(exc).__name__,
        ) from exc


def _ensure_defaults() -> None:
    """Lazily register the three built-in engines the first time needed."""
    global _defaults_loaded
    if _defaults_loaded:
        return
    _defaults_loaded = True

    # We import the registration helpers inline to avoid pulling in heavy
    # packages (pyrosetta, lightdock, haddock) at import-time.
    _register_lightdock()
    _register_haddock()
    _register_rosetta()


# ---------------------------------------------------------------------------
# Built-in engine registration
# ---------------------------------------------------------------------------

def _register_lightdock() -> None:
    from ppinsight.collect_scores import _parse_lightdock

    def _detect(directory: str) -> bool:
        """LightDock outputs always contain numbered swarm_* directories."""
        return bool(glob.glob(os.path.join(directory, "swarm_*")))

    def _run(rec_pdb, lig_pdb, output_root, pair_label, **kw):
        from ppinsight.pdb_to_lightdock import (
            LightDockSimulationError,
            _cleanup_previous_outputs,
            _stage_cleaned_lightdock_inputs,
            lightdock_pipeline,
            make_output_dir,
        )
        try:
            workdir = make_output_dir(rec_pdb, lig_pdb, method="lightdock_runs",
                                      base_root=output_root)
            run_opts = {
                "steps": kw.get("steps", 10),
                "swarms": kw.get("swarms"),
                "glowworms": kw.get("glowworms"),
                "cores": kw.get("cores", 1),
                "anm": kw.get("anm", True),
                "scoring": kw.get("scoring"),
                "skip_postprocess": kw.get("skip_postprocess", False),
            }
            auto_clean = kw.get("auto_clean_pdb", False)

            try:
                lightdock_pipeline(
                    receptor_pdb=rec_pdb,
                    ligand_pdb=lig_pdb,
                    working_dir=workdir,
                    **run_opts,
                )
            except subprocess.CalledProcessError as exc:
                if (
                    run_opts.get("anm", True)
                    and _is_lightdock_anm_atom_mismatch(exc.output)
                ):
                    if not _is_interactive_session():
                        raise RuntimeError(
                            "LightDock ANM setup atom mismatch detected. "
                            "Interactive confirmation is required before "
                            "retrying without ANM. Rerun in an interactive "
                            "terminal, or use --lightdock-no-anm."
                        ) from exc

                    should_retry = _confirm_lightdock_retry_without_anm(
                        rec_pdb,
                        lig_pdb,
                    )
                    if should_retry:
                        print(
                            "  [lightdock] Retrying this pair with ANM disabled."
                        )
                        _cleanup_previous_outputs(workdir)
                        retry_opts = dict(run_opts)
                        retry_opts["anm"] = False
                        lightdock_pipeline(
                            receptor_pdb=rec_pdb,
                            ligand_pdb=lig_pdb,
                            working_dir=workdir,
                            **retry_opts,
                        )
                    else:
                        raise RuntimeError(
                            "LightDock ANM setup atom mismatch detected and retry "
                            "without ANM was not confirmed."
                        ) from exc
                elif _is_lightdock_anm_atom_mismatch(exc.output):
                    raise RuntimeError(
                        "LightDock ANM setup atom mismatch detected. "
                        "Hint: rerun in an interactive terminal to confirm "
                        "a retry without ANM, or use --lightdock-no-anm."
                    ) from exc
                else:
                    raise
            except LightDockSimulationError as exc:
                if exc.cleanable and auto_clean:
                    print(
                        "  [lightdock] Unsupported residue detected; "
                        "retrying this pair with cleaned protein-only inputs."
                    )
                    clean_rec, clean_lig, _ = _stage_cleaned_lightdock_inputs(
                        rec_pdb,
                        lig_pdb,
                        workdir,
                    )
                    _cleanup_previous_outputs(workdir)
                    lightdock_pipeline(
                        receptor_pdb=clean_rec,
                        ligand_pdb=clean_lig,
                        working_dir=workdir,
                        **run_opts,
                    )
                else:
                    raise

            return workdir
        except Exception as exc:
            log_path = None
            if "workdir" in locals():
                log_path = os.path.join(workdir, "lightdock.log")
            raise EngineRunError(
                str(exc),
                output_dir=workdir if "workdir" in locals() else None,
                log_path=log_path,
                error_type=type(exc).__name__,
            ) from exc

    _registry["lightdock"] = EnginePlugin(
        name="lightdock",
        runner=_run,
        parser=_parse_lightdock,
        detector=_detect,
        description="LightDock swarm-based docking",
    )


def _register_haddock() -> None:
    from ppinsight.collect_scores import _parse_haddock

    def _detect(directory: str) -> bool:
        """Detect HADDOCK by capri_ss.tsv OR clustfcc/clustrmsd step directories."""
        # Detect HADDOCK by capri_ss.tsv OR clustfcc/clustrmsd dirs
        if glob.glob(os.path.join(directory, "**", "capri_ss.tsv"), recursive=True):
            return True
        if (glob.glob(os.path.join(directory, "*_clustfcc"))
                or glob.glob(os.path.join(directory, "*_clustrmsd"))):
            return True
        return False

    def _run(rec_pdb, lig_pdb, output_root, pair_label, **kw):
        from ppinsight.pdb_to_haddock import haddock_pipeline
        from ppinsight.utils import _project_root
        try:
            run_opts = dict(kw)
            sampling = int(run_opts.pop("sampling", 10000))
            select_top = int(run_opts.pop("select_top", 400))
            run_dir, cfg_path, _ = haddock_pipeline(
                rec_pdb, lig_pdb,
                workspace_root=_project_root(),
                base_root=output_root,
                run_haddock=True,
                sampling=sampling,
                select_top=select_top,
                **run_opts,
            )
            return str(run_dir)
        except Exception as exc:
            run_dir = getattr(exc, "run_dir", None)
            raise EngineRunError(
                str(exc),
                output_dir=str(run_dir) if run_dir else None,
                log_path=os.path.join(str(run_dir), "log") if run_dir else None,
                error_type=type(exc).__name__,
            ) from exc

    _registry["haddock"] = EnginePlugin(
        name="haddock",
        runner=_run,
        parser=_parse_haddock,
        detector=_detect,
        description="HADDOCK integrative docking",
    )


def _register_rosetta() -> None:
    from ppinsight.collect_scores import _parse_rosetta

    def _detect(directory: str) -> bool:
        """Detect Rosetta by clustered_scores.csv, .sc files, or PPInsight CSV."""
        # Detect clustered_scores.csv, native .sc score files,
        # or PPInsight docking_scores.csv
        if os.path.isfile(os.path.join(directory, "clustered_scores.csv")):
            return True
        if os.path.isfile(os.path.join(directory, "docking_scores.csv")):
            return True
        if glob.glob(os.path.join(directory, "*.sc")):
            return True
        return False

    _registry["rosetta"] = EnginePlugin(
        name="rosetta",
        runner=_run_rosetta_batch,
        parser=_parse_rosetta,
        detector=_detect,
        description="Rosetta protein-protein docking (requires PyRosetta)",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def register(plugin: EnginePlugin) -> None:
    """Register (or replace) a docking-engine plugin.

    Parameters
    ----------
    plugin : EnginePlugin
        The plugin to register.  ``plugin.name`` is used as the key.

    Raises
    ------
    TypeError
        If *plugin* is not an :class:`EnginePlugin`.
    """
    if not isinstance(plugin, EnginePlugin):
        raise TypeError(f"Expected EnginePlugin, got {type(plugin).__name__}")
    _registry[plugin.name] = plugin


def get(name: str) -> EnginePlugin:
    """Return the plugin for *name*.

    Raises ``KeyError`` if the engine hasn't been registered.
    """
    _ensure_defaults()
    return _registry[name]


def list_engines() -> list[str]:
    """Return sorted list of registered engine names."""
    _ensure_defaults()
    return sorted(_registry)


def get_runner(name: str):
    """Shorthand: return the *runner* callable for *name*, or *None*."""
    return get(name).runner


def get_parser(name: str):
    """Shorthand: return the *parser* callable for *name*, or *None*."""
    return get(name).parser


def detect_engine(directory: str) -> str:
    """Detect which engine produced *directory*.

    Iterates registered engines in deterministic order and returns the
    first whose ``detector`` returns *True*.

    Raises ``ValueError`` if no engine matches.
    """
    _ensure_defaults()
    for name in sorted(_registry):
        plugin = _registry[name]
        if plugin.detector and plugin.detector(directory):
            return name
    raise ValueError(
        f"Cannot detect docking engine for '{directory}'. "
        f"Known engines: {', '.join(sorted(_registry))}."
    )
