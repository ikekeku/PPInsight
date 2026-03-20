"""ppinsight package init"""
__all__ = [
	"pdb_to_lightdock",
	"pdb_to_haddock",
	"pdb_to_rosetta",
	"protein_fetch",
	"visualizer",
	"collect_scores",
	"utils",
	"registry",
	"parse_pairs",
	"batch_dock",
	"quality",       # lazy-imported; depends on optional DockQ package
]

# Import modules that have no heavy external dependencies.
#
# NOT eagerly imported (and why):
#   - rosetta / pdb_to_rosetta: depend on PyRosetta (~1.5 GB optional download).
#   - quality: depends on DockQ, an optional dependency installed via
#     ``pip install ppinsight[quality]``.  The module guards its own import
#     and raises a helpful error if DockQ is missing.
#
# Users who only need HADDOCK, LightDock, or the visualizer can
# ``import ppinsight`` without PyRosetta or DockQ installed.
from . import pdb_to_lightdock
from . import pdb_to_haddock
from . import protein_fetch
from . import visualizer
from . import collect_scores
from . import utils
from . import registry
from . import parse_pairs
from . import batch_dock
