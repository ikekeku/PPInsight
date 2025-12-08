"""ppinsight package init"""
__all__ = ["pdb_to_lightdock", "pdb_to_haddock", "protein_fetch"]
from . import pdb_to_lightdock
from . import pdb_to_haddock
from . import protein_fetch

# try:
#     subprocess.run(["rosetta" "--version"], check=True)
# except:
#     raise RuntimeError("rosetta not installed")

# pyrosetta_installer.install(skipifinstalled=True)
