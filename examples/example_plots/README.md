
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
| **violin** | `violin_*.png` | Full score distribution per engine (faceted by pair when ≥ 2 pairs) |
| **ridge** | `ridge_*.png` | Overlapping KDE density curves, one row per engine |
| **scatter** | `scatter_*.png` | Model-agreement scatter — same metric, same pairs, one engine per axis |
| **difference** | `difference_*.png` | Bland-Altman histogram of per-pair score differences between two engines |
| **roc** | `roc_*.png` | ROC curve per engine (requires interaction / non-interaction labels) |
| **quality_bar** | `quality_bar_*.png` | 100 % stacked pose-level CAPRI tier bar (incorrect / acceptable / medium / high) |
| **cdf** | `cdf_*.png` | Step-CDF per engine; CAPRI threshold lines shown for DockQ metrics |

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
- Metric name when relevant (e.g. `dockq`, `score`, `prodigy_ddg`)
