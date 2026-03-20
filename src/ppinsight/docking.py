"""Compatibility shim: re-export the DockingPipeline at top-level

Some tests and docs import `ppinsight.docking.DockingPipeline`. The
implementation lives under `ppinsight.rosetta.pipeline`. Provide a small
shim so both import paths work.

If PyRosetta is not installed, importing this module will raise
ImportError (same as importing the rosetta sub-package directly).
"""
try:
    from .rosetta.pipeline import DockingPipeline
except ImportError as _exc:
    raise ImportError(
        "DockingPipeline requires PyRosetta. Install with: "
        "python -c \"import pyrosetta_installer; "
        "pyrosetta_installer.install_pyrosetta()\""
    ) from _exc

__all__ = ["DockingPipeline"]
