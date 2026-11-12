# 4_visualisation

Figure producers of the paper. They read the tables written by `3_analysis/` and are run by the same driver, in the stages
`figs`, `s7`, `fig1` and `copy` (see `3_analysis/README.md`):

```
ARCHI_STAGE=H5 python3 3_analysis/integration/run_h3.py figs s7 fig1 copy
```

The R scripts run from `4_visualisation/figures/` with `ARCHI_PKG` pointing at the folder that holds `results_table.csv`
(`3_analysis/integration/inputs` for the human panels, `3_analysis/package` for the others) and `ARCHI_DRAW` set to `H`
(human) or `M` (mouse). Panels go to the script folder or to `figures/output/<figure>/`; `3_analysis/integration/copy_figures.py`
copies every final figure to `outputs/figures/`, except Fig. 3, which `Fig3_S3/assemble_figure3.py` writes there directly.

## Producers

| figure | producer | panel folder |
|---|---|---|
| Fig. 1 | `figure1/make_figure1_R2.py` (values from `3_analysis/integration/figure1_values.csv`), then `rsvg-convert` (PNG 300 dpi, PDF) and `sips` | `figure1/` |
| Fig. 2 | `figures/Fig2_S1/Figure_2_and_S1_restyled.R` (human) | `figures/Fig2_S1/` |
| Fig. S1 | `figures/Fig2_S1/Figure_2_and_S1_restyled.R` (mouse) | `figures/output/S1/` |
| Fig. 3 | A-C `figures/Fig3_S3/Figure_3_and_S3_restyled.R` (human); D `figures/Fig3_S3/build_Figure_3D_source_data_h3.py` and `Figure_3D_logos.R`; composite `figures/Fig3_S3/assemble_figure3.py` | `figures/Fig3_S3/` |
| Fig. S3 | `figures/Fig3_S3/Figure_3_and_S3_restyled.R` (mouse) | `figures/output/S3/` |
| Fig. 4, Fig. S4 | `figures/Fig4_S4/Figure_4_and_S4_restyled.R` (also writes `Figure_4_stats.csv`, `Figure_S4_stats.csv`) | `figures/output/Fig4/` |
| Fig. S2 | `figures/FigS2/Figure_S2_restyled.py` | `figures/output/S2/` |
| Fig. S5 | `figures/FigS5/Figure_S5_restyled.R` | `figures/output/S5/` |
| Fig. S6 | `figures/FigS6/Figure_S6_restyled.R` | `figures/output/S6/` |
| Fig. S7 | `figures/FigS7/build_Figure_S7_source_data_h3.py`, `hocomoco11_check/subtype_vs_hocomoco11.py`, `subtype_vs_named_v11.py`, `figures/FigS7/Figure_S7_motif_subtypes.R` | `figures/FigS7/` |
| Fig. S8 | `3_analysis/low_fpr_metrics/b2_figure_h3.py` | `3_analysis/low_fpr_metrics/` |
| Fig. S9 | `3_analysis/top_peaks/b3_summary_figure_w2.py` | `3_analysis/top_peaks/` |
| Fig. S10 | `3_analysis/feature_importance/scripts/b4_06_figure.py` | `3_analysis/feature_importance/` |
| Fig. S11 | `3_analysis/transfer_drivers/B5_02_figure.py` | `3_analysis/transfer_drivers/` |
| Fig. S12 | `3_analysis/random_split/b9_figure_h3.py` | `3_analysis/random_split/` |

`figures/FigS7/hocomoco11_check/cluster_vs_v11.py` compares the ablation subtypes with HOCOMOCO v11 (Sup. Table 9 sheets).
The HOCOMOCO v11 HUMAN mono models and their annotation are in `results/datasets/hocomoco11/`.

## Style

`figures/archi_style.R` and `figures/archi_style.py` hold the shared palette, fonts and the PDF + PNG writer: all text
black, Arial, no transparency, PNG rendered from the PDF at 300 dpi with `pdftoppm`. `figures/check_palette.py` reports the
colour distances of the palette under simulated colour vision deficiency
(`cd 4_visualisation/figures && python3 check_palette.py`).

## Software

R 4.4.3 with dplyr 1.2.0, ggplot2 4.0.2, patchwork 1.3.2, cowplot 1.2.0, ggrepel 0.9.8, ggExtra 0.11.0, viridis 0.6.5,
tidyr 1.3.2, ggbeeswarm 0.7.3, ggseqlogo 0.2, png 0.1.9. Python 3.9.12 with matplotlib 3.5.2, pandas 2.3.0, numpy 1.26.4,
Pillow 11.3.0, pypdf 6.9.1. Poppler 24.04 (`pdftoppm`, `pdfimages`), `rsvg-convert` 2.50.7, macOS `sips`, Arial.
