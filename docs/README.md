# PPInsight Documentation

This directory contains design notes, operational guidance, and background
material for PPInsight.

## Contents

- [Docking-model clustering](#docking-model-clustering)
  - [LightDock: per-swarm BSAS clustering](#lightdock-per-swarm-bsas-clustering)
  - [HADDOCK: Fraction of Common Contacts (FCC)](#haddock-fraction-of-common-contacts-fcc)
  - [Rosetta: hierarchical Cα-RMSD clustering](#rosetta-hierarchical-cα-rmsd-clustering)
  - [Practical interpretation](#practical-interpretation)
- [Reference-based validation](#reference-based-validation)
- [Sources](#sources)

## Docking-model clustering

Docking engines can generate many structurally similar poses (also called
*models* or *decoys*). Clustering groups similar poses into distinct candidate
binding modes. It helps avoid treating a large number of nearly identical poses
as independent evidence.

A cluster is not proof that a pose is biologically correct. It is evidence that
the docking search repeatedly found a similar geometry. Use clusters together
with the engine's within-run score and, when available, reference-based quality
metrics such as DockQ.

> **Important:** engine-native scores and cluster sizes are meaningful within a
> single docking run for a single pair. Do not compare their raw values across
> different protein pairs or engines.

### LightDock: per-swarm BSAS clustering

PPInsight runs the standard LightDock post-processing workflow after a
successful simulation:

1. `lgd_generate_conformations.py` writes PDB models for each swarm.
2. `lgd_cluster_bsas.py` clusters the models **within each swarm**.
3. `lgd_rank.py` ranks the swarm cluster representatives globally.

The clustering log reports pairwise RMSD values as it processes the models in a
swarm. A small RMSD means two poses have similar geometry; a larger RMSD can
start a separate cluster. The exact assignment depends on the BSAS clustering
criterion rather than a universal "good" RMSD threshold.

For example, messages such as these are expected:

```text
RMSD between 23 and 15 is 0.044
Glowworm 15 goes into cluster 0
...
New cluster 1
Cluster result written to .../swarm_8/cluster.repr file
```

The first pair is highly similar and is placed in the same cluster. `New
cluster` means the pose was sufficiently dissimilar from the existing cluster
representatives to begin another cluster. Each swarm writes a `cluster.repr`
file. `lgd_rank.py` then reads those representatives and writes a global,
score-ranked list, normally `rank_by_scoring.list`.

**Where to look**

- `swarm_<N>/gso_<steps>.out`: GSO simulation output for a swarm.
- `swarm_<N>/cluster.repr`: representatives selected from that swarm.
- `rank_by_scoring.list`: global ranking of representatives across swarms.

LightDock clustering therefore reduces redundancy *within* individual swarms;
the final ranking compares the representatives from all processed swarms.

### HADDOCK: Fraction of Common Contacts (FCC)

The PPInsight-generated HADDOCK workflow uses the `[clustfcc]` module after
docking/refinement. Unlike RMSD clustering, FCC compares the intermolecular
contacts made by models. Models with a sufficiently similar contact pattern are
placed in the same cluster.

The default HADDOCK `clustfcc` settings used by PPInsight retain single-member
clusters (`min_population = 1`). HADDOCK's documentation describes an FCC
cutoff of `0.6` as the default minimum fraction of common contacts and a default
minimum population of four when users do not override it. The appropriate cutoff
is system-dependent: a tighter criterion separates more binding modes, whereas
a looser criterion merges more poses.

The generated workflow subsequently uses `[seletopclusts]` to select
score-ranked models from the highest-ranked clusters. A standard protein-protein
HADDOCK workflow may alternatively use FCC clustering after rigid-body docking
to preserve structural diversity before refinement.

**Where to look**

- `<run>/6_clustfcc/`: FCC clustering results in a standard PPInsight HADDOCK
  run.
- `<run>/7_seletopclusts/`: selected models from top clusters.
- `capri_clt.tsv`: cluster-level CAPRI/DockQ metrics when a reference structure
  was supplied and HADDOCK's `caprieval` has cluster information.

### Rosetta: hierarchical Cα-RMSD clustering

PPInsight clusters Rosetta decoys after docking unless
`--rosetta-no-cluster` is supplied. It:

1. orders decoys by interface score (`I_sc`; lower is better),
2. takes the configured number of best-scoring decoys (default: 200),
3. calculates pairwise Cα-RMSD values,
4. uses average-linkage hierarchical clustering, and
5. cuts the hierarchy at the configured RMSD cutoff (default: 4 Å).

PPInsight renumbers the clusters so cluster `0` is the largest cluster; ties
are broken by the best `I_sc`. Within each cluster, `cluster_rank = 0` identifies
the best-scoring member and `is_top_of_cluster` is `True` for that member. The
best-scoring decoy in the largest cluster is a useful representative prediction,
but it should be evaluated alongside other high-ranking clusters when the search
found multiple plausible interfaces.

**Where to look**

- Rosetta score tables include `cluster`, `cluster_size`, `cluster_rank`, and
  `is_top_of_cluster` when clustering completed.
- Use `--rosetta-cluster-top-n` to control how many top-scoring decoys enter
  the clustering calculation.
- Use `--rosetta-rmsd-cutoff` to control the Cα-RMSD distance at which the
  hierarchy is split.

### Practical interpretation

| Observation | Interpretation | Recommended follow-up |
|---|---|---|
| One large cluster with favorable internal scores | Search repeatedly converged on one pose family. | Inspect the cluster representative and validate with DockQ if a native complex is available. |
| Several well-populated clusters | More than one candidate binding mode was found. | Keep representatives from multiple top clusters for review or refinement. |
| Mostly singleton clusters | The search did not converge under the current sampling/settings. | Increase sampling, add reliable restraints, or examine input structures. |
| A high-scoring isolated model | It may be a legitimate alternative or a scoring artifact. | Do not discard it solely because its cluster is small; inspect contacts and validate independently. |

## Reference-based validation

Clustering measures agreement among predicted models; it does not measure
agreement with the true complex. When a native reference is available, evaluate
selected representatives with:

```bash
ppinsight quality <predicted_model.pdb> <native_complex.pdb>
```

DockQ combines native-contact fraction, interface RMSD, and ligand RMSD to
assess model quality. See [scoring_evaluation_framework.md](scoring_evaluation_framework.md)
for score interpretation and CAPRI/DockQ thresholds.

## Sources

- [LightDock workflow](https://lightdock.org/tutorials/0.9.3/simple_docking)
  and PPInsight's LightDock pipeline implementation.
- [HADDOCK3 analysis modules](https://www.bonvinlab.org/haddock3-user-manual/modules/analysis.html),
  especially `clustfcc`, `clustrmsd`, and `seletopclusts`.
- [HADDOCK3 protein-protein docking scenario](https://www.bonvinlab.org/haddock3-user-manual/docking_scenarios/prot-prot.html).
- Rodrigues JP *et al.* (2012), *Proteins* 80, 1810–1817: Fraction of Common
  Contacts clustering.
