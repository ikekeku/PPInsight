"""ppinsight package init"""
__all__ = [
	"pdb_to_lightdock",
	"pdb_to_haddock",
	"protein_fetch",
	"visualizer",
]

# Import the core helpers that tests commonly use. `visualizer` is not
# imported here to avoid import-time failures in environments where its
# annotations are not supported; it can be imported lazily by callers.
from . import pdb_to_lightdock
from . import pdb_to_haddock
from . import protein_fetch
# from . import visualizer

# try:
#     subprocess.run(["rosetta" "--version"], check=True)
# except:
#     raise RuntimeError("rosetta not installed")

# pyrosetta_installer.install(skipifinstalled=True)
