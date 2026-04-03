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
from dataclasses import dataclass
from typing import Protocol, TYPE_CHECKING

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
        from ppinsight.pdb_to_lightdock import lightdock_pipeline, make_output_dir
        import sys
        try:
            workdir = make_output_dir(rec_pdb, lig_pdb, method="lightdock_runs",
                                      base_root=output_root)
            lightdock_pipeline(
                receptor_pdb=rec_pdb,
                ligand_pdb=lig_pdb,
                working_dir=workdir,
                steps=kw.get("steps", 10),
                swarms=kw.get("swarms"),
                anm=kw.get("anm", True),
                scoring=kw.get("scoring"),
            )
            return workdir
        except Exception as exc:
            print(f"  [lightdock] FAILED: {exc}", file=sys.stderr)
            return None

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
        import sys
        try:
            run_dir, cfg_path, _ = haddock_pipeline(
                rec_pdb, lig_pdb,
                workspace_root=_project_root(),
                run_haddock=False,
                **kw,
            )
            return str(run_dir)
        except Exception as exc:
            print(f"  [haddock] FAILED: {exc}", file=sys.stderr)
            return None

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
        """Detect Rosetta by clustered_scores.csv, .sc files, or legacy CSV."""
        # Detect clustered_scores.csv, native .sc score files, OR legacy PPInsight docking_scores.csv
        if os.path.isfile(os.path.join(directory, "clustered_scores.csv")):
            return True
        if os.path.isfile(os.path.join(directory, "docking_scores.csv")):
            return True
        if glob.glob(os.path.join(directory, "*.sc")):
            return True
        return False

    # No default runner — PyRosetta is heavy and invocation is complex.
    # Users wanting Rosetta in batch mode can register their own runner.
    _registry["rosetta"] = EnginePlugin(
        name="rosetta",
        runner=None,
        parser=_parse_rosetta,
        detector=_detect,
        description="Rosetta protein-protein docking (parse-only by default)",
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
