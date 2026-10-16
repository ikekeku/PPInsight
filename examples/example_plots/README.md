
# Example Plots (FABRICATED DATA)

⚠️ **All plots in this directory were generated with fabricated
(synthetic) data.** No real docking simulations were run. These images
exist solely to demonstrate the visual output of each `ppinsight compare`
plot type.

## How they were generated

```bash
cd PPInsight
python examples/generate_example_plots.py
```

The script (`generate_example_plots.py`) writes a standalone fabricated
scores file at `examples/example_plots/fabricated_scores.tsv`, preserves
the same row-level traceability fields used by collected score tables
(`run_id`, `pose_id`, `output_path`), and then renders the plot gallery
from that synthetic TSV.

This file is intentionally **not** part of the real PPInsight pipeline.
It exists only so the example plot gallery is reproducible from a single
fabricated long-format scores table.

## Plot types

| Plot type | Files | What it shows |
|-----------|-------|---------------|
| **violin** | `violin_*.png` | Full score distribution per engine (faceted by pair when ≥ 2 pairs): DockQ for 1, 2 and 3 engines with 1, 2 and 3 pairs, plus HADDOCK's own score and LightDock's luciferin score |
| **strip** | `strip_*.png` | One dot per pose with a box summary, per engine (faceted by pair when ≥ 2 pairs): 1, 2 and 3 engines with 1, 2 and 3 pairs |
| **scatter** | `scatter_*.png` | Model-agreement scatter — same metric, same pairs, one engine per axis |
| **roc** | `roc_<Engine>_<score>.png` | ROC curve for one engine on its own score (requires interaction / non-interaction labels; DockQ is not allowed because non-binders have no native) |
| **quality_bar** | `quality_bar_*.png` | Percentage of poses in each CAPRI tier (from DockQ) per engine, for one protein pair |
| **cdf** | `cdf_*.png` | Step-CDF per engine for one protein pair; CAPRI threshold lines shown for DockQ metrics |
| **funnel** | `funnel_*.png` | Rank audit for one target: each engine's own score vs DockQ for the same poses (uses the separate `demo_rec` / `demo_lig` target) |
| **topn** | `topn_*.png` | Rank audit for one target: best DockQ within each engine's top N poses (same `demo_rec` / `demo_lig` target) |

All DockQ, Fnat and RMSD values in the gallery are the native-referenced
`quality_*` metrics (what `ppinsight quality` writes).  Files whose names end
in `_dockq` show `quality_dockq`; the name is kept short.

## Generated files

| File | Purpose |
|------|---------|
| `fabricated_scores.tsv` | Synthetic long-format scores table used to generate the PNG gallery |
| `*.png` | Plot examples rendered from `fabricated_scores.tsv` |

## Filename convention

```
<plot_type>_<N>eng_<N>pair_<metric>.png
```

- `1eng` / `2eng` / `3eng` — number of docking engines
- `1pair` / `2pair` / `3pair` — number of protein pairs
- Metric name when relevant (e.g. `dockq`, `score`, `i_sc`)

Exceptions:

- Single-engine violins for every metric name the engine instead of the
  count: `violin_<Engine>_1pair_<metric>.png` (e.g.
  `violin_LightDock_1pair_luciferin_score.png`).
- The scatter plot names both engines: `scatter_<EngineA>_vs_<EngineB>_<metric>.png`.
  It uses the 50 labelled binder pairs, which have a native complex and
  therefore a native-referenced DockQ.
- The ROC plot names the engine and its own score: `roc_<Engine>_<score>.png`.
- Strip and the multi-engine violins use the same number of engines and
  pairs (`1eng_1pair`, `2eng_2pair`, `3eng_3pair`).
