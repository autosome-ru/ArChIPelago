# 4_visualisation

Figure producers of the paper, each in the folder of its figure, next to its output. They read the tables written by
`../3_analysis/` and are run by the driver `../3_analysis/integration/run_h3.py` (stages `figs`, `b234`, `b5`, `b9` and
`s7`; see `../3_analysis/README.md`). Figure 1 is drawn by hand and is not in the repository; its numbers are in
`../3_analysis/integration/figure1_values.csv`.

```
ARCHI_STAGE=H5 python3 3_analysis/integration/run_h3.py figs s7
```

(run from the repository root). Each figure folder has a README with its producer, inputs and stage.

## Folders

| folder | figures | producer |
|---|---|---|
| `main_figures/Figure_2/` | Fig. 2, Fig. S1 | `Figure_2_and_S1.R` |
| `main_figures/Figure_3/` | Fig. 3, Fig. S3 | `Figure_3_and_S3.R` (A-C), `build_Figure_3D_source_data_h3.py` and `Figure_3D_logos.R` (D), `assemble_figure3.py` (composite) |
| `main_figures/Figure_4/` | Fig. 4, Fig. S4 | `Figure_4_and_S4.R` |
| `supplementary_figures/Figure_S2/` | Fig. S2 | `Figure_S2.py` |
| `supplementary_figures/Figure_S5/`, `supplementary_figures/Figure_S6/` | Fig. S5, S6 | `Figure_S5.R`, `Figure_S6.R` |
| `supplementary_figures/Figure_S7/` | Fig. S7 | `build_Figure_S7_source_data_h3.py`, `Figure_S7.R`; HOCOMOCO v11 comparison in `hocomoco11_check/` |
| `supplementary_figures/Figure_S8/` to `supplementary_figures/Figure_S12/` | Fig. S8-S12 | `Figure_S8.py` to `Figure_S12.py` |
| `style/` | shared style | `archi_style.R`, `archi_style.py`, `check_palette.py` |

The joint scripts in `main_figures/` write Fig. S1, S3 and S4 into `supplementary_figures/Figure_S1/`, `Figure_S3/` and
`Figure_S4/`. The R scripts of Fig. 2, 3 and 4 read the results table or the figure source data from the folder in
`ARCHI_PKG` (`3_analysis/integration/inputs` for the human panels of Fig. 2 and 3, `3_analysis/package` for the mouse panels and for Fig. 4); `ARCHI_DRAW=H`
draws the human figure only, `ARCHI_DRAW=M` the mouse one.

## Style

`style/archi_style.R` and `style/archi_style.py` hold the shared palette, fonts and the PDF and PNG writer: all text
black, Arial, no transparency, PNG at 600 dpi (rendered from the PDF with `pdftoppm` by the R writer, saved by matplotlib in the Python writer). `style/check_palette.py` reports the
colour distances of the palette under simulated colour vision deficiency (`cd 4_visualisation/style && python3 check_palette.py`).

## Software

R 4.4.3 with dplyr 1.2.0, ggplot2 4.0.2, patchwork 1.3.2, cowplot 1.2.0, ggrepel 0.9.8, ggExtra 0.11.0, viridis 0.6.5,
tidyr 1.3.2, ggbeeswarm 0.7.3, ggseqlogo 0.2, png 0.1.9. Python 3.9.12 with matplotlib 3.5.2, pandas 2.3.0, numpy 1.26.4,
Pillow 11.3.0, pypdf 6.9.1. Poppler 24.04 (`pdftoppm`, `pdfimages`), Arial.
