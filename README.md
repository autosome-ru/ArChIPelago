# ArChIPelago <img src='./Archipelago.png' width='55'>

[![tests](https://github.com/autosome-ru/ArChIPelago/actions/workflows/tests.yml/badge.svg)](https://github.com/autosome-ru/ArChIPelago/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14927303.svg)](https://doi.org/10.5281/zenodo.14927303)
[![release](https://img.shields.io/github/v/tag/autosome-ru/ArChIPelago?label=release)](https://github.com/autosome-ru/ArChIPelago/tree/v2.0.0)
[![License: WTFPL](https://img.shields.io/badge/license-WTFPL-lightgrey.svg)](http://www.wtfpl.net/)

**Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites**

ArChIPelago combines the position weight matrices (PWMs) of one transcription factor (TF) into a joint model with classic
machine learning (Random Forest, logistic regression, XGBoost, bagging) to predict TF binding sites in genomic sequences.
This repository holds the code that built the training and test data of 35 human and mouse TFs, trained and evaluated the
models and produced the tables, figures and numbers of the paper, together with the result tables of the final run. The
command-line tool [ArChIPelago-TFBS-finder](https://github.com/autosome-ru/ArChIPelago-TFBS-finder) (git submodule
`ArChIPelago-TFBS-finder/`) scans a FASTA file with the trained models.

## Results (Random Forest on monoPWMs and diPWMs, 35 TFs)

| test set | best single monoPWM | ArChIPelago | TFs improved |
|---|---|---|---|
| human, auROC (median) | 0.874 | 0.902 | 35 of 35 |
| human, auPRC (median) | 0.288 | 0.336 | 35 of 35 |
| mouse, human-trained models, auROC | 0.857 | 0.886 | 34 of 35 |
| mouse, human-trained models, auPRC | 0.230 | 0.307 | 34 of 35 |
| mouse, mouse-trained models, auROC | 0.873 | 0.897 | 35 of 35 |
| mouse, mouse-trained models, auPRC | 0.248 | 0.356 | 35 of 35 |

The baseline is the best single monoPWM of each TF, selected on the human training set by the metric of the row. Human
rows are random seed 0, mouse rows the mean of seeds 0 and 1. Every value is in `outputs/numbers/NUMBERS_R2_35TF.csv`
with its source table.

Supported TFs (35): ANDR, AP2A, CEBPB, COE1, CTCF, E2F4, ERG, ESR1, FLI1, GATA1, GATA2, GATA3, GCR, HNF4A, IRF1, IRF4,
JUND, MAFK, MAX, MYC, PPARG, PRGR, REST, RUNX1, RXRA, SOX2, SPI1, SRF, STA5A, STAT1, STAT3, TAL1, TF65, TFE2, USF2.

Chromosomes: human training chr2-7, 9-20 and 22, test chr1, 8 and 21; mouse training chr2-7 and 9-18, test chr1, 8 and
19 (`1_data_preparation/h5_splits.py`).

## Repository structure

```
ArChIPelago/
├── 1_data_preparation/        positives, GC-matched negatives, SARUS feature tables, Slim inputs, data checks, loader
├── 2_training/                Random Forest and other models, Slim and diChIPMunk models, evaluation tables
├── 3_analysis/
│   ├── server/                compute-heavy analyses (metrics, feature importance, saturation, Fig. 4 refits, runtime)
│   ├── integration/           driver run_h3.py and the producers of the integrated tables and numbers
│   ├── package/               results table, saturation, cross-species tables, figure source data
│   ├── low_fpr_metrics/  top_peaks/  random_split/  feature_importance/  transfer_drivers/
│   └── supplementary_tables/  builders of Sup. Tables 1-10
├── 4_visualisation/
│   ├── main_figures/          Figure_2, Figure_3, Figure_4, each with its producer
│   ├── supplementary_figures/ Figure_S1 to Figure_S12, each with its producer
│   └── style/                 shared palette and figure writer
├── results/                   tracked inputs of the local stages
├── outputs/                   supplementary tables and numbers of the final run
├── tests/                     pytest checks of the tables, figures and numbers
├── archi_env.py  archi_env.sh environment variables for the data outside the repository
├── ArChIPelago-TFBS-finder/   scanning tool (git submodule)
└── sarus/  seqtk/             SPRY-SARUS and seqtk (git submodules)
```

Each folder has its own README. Figure 1 of the paper is drawn by hand and is not in the repository; the numbers it
shows are in `3_analysis/integration/figure1_values.csv`.

## How to rerun

Steps 1-3 run on a compute server and need the inputs listed under "Inputs on Zenodo": `1_data_preparation/`, then
`2_training/`, then `3_analysis/server/` (see their READMEs). Their result tables are tracked in `results/server/`.

The local chain rebuilds every table, figure and number from `results/` with one command, run from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
ARCHI_STAGE=H5 \
python3 3_analysis/integration/run_h3.py all
```

It checks every input against `results/server/MIRROR_MANIFEST_md5.txt`, then runs the 11 stages of `all` in this order
and stops at the first failing command (about 3 min on a laptop): `pkg`, `figs`, `b234`, `b5`, `b9`, `sup`, `s7`,
`porig`, `v1`, `fig1`, `numbers`. Single stages run by name, for example
`python3 3_analysis/integration/run_h3.py sup`. The optional stage `gate` is not part of `all`. What each stage runs is
listed in `3_analysis/README.md`.

## Requirements

Local chain: Python 3.9 or later with pandas 2.3, numpy 1.26, scipy 1.13, matplotlib 3.5, scikit-learn, rpy2 3.5,
openpyxl, pyyaml, pypdf; R 4.4 with PRROC 1.4, dplyr, ggplot2, patchwork, cowplot, ggrepel, ggExtra, viridis, tidyr,
ggbeeswarm, ggseqlogo and png; poppler (`pdftoppm`, `pdfimages`, `pdftotext`); the Arial font. auROC and auPRC come
from PRROC through rpy2 (`1_data_preparation/prroc_metrics.py`). rpy2 must load the R that holds PRROC: set `R_HOME` to
that R. Steps 1-3 use the versions in `2_training/VERSIONS.txt`.

## Inputs on Zenodo

- Record [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303) (resolves to the latest version, v2.0.0):
  training and test sets (`train.tar.gz`, `test.tar.gz`), PWMs (`PWMs_{mono,di}_{HUMAN,MOUSE}.tar.gz`), Random Forest
  models (`Models.tar.gz`), Slim and diChIPMunk models (`Slim_diChIPMunk_models_v2.tar.gz`), RepeatMasker masks
  (`repeatmasker_hg38_mm10.tar.gz`) and the test-set scores that the stage `gate` reads (`scores_H5_human.tar.gz`). Its
  `README_v2.md` lists every archive. The MACS peak files (`macs.tar.gz`) and HOCOMOCO v11 PWMs (`hocomoco11.tar.gz`) are
  Zenodo archives as well. In these files `train` always marks the training chromosomes; `1_data_preparation/README.md`
  says how the loader reads them.
- Record v1.0.0 [10.5281/zenodo.23036475](https://doi.org/10.5281/zenodo.23036475): the four PWM archives (tracked in
  `results/datasets/zenodo_v1.0.0/`, unpacked by the driver) and `Slim.tar.gz` (md5
  `ce12823de573676ea0cee7dce5cadb79`). Steps 2 and 3 need the Slim jar, JDK 8 and `ytilib/` from it: unpack it in the
  repository root, which creates `Slim/` (not tracked).

The reference genomes (hg38, GRCm38.p6) come from the user's own copy, found through `ARCHI_GENOMES` (`archi_env.py`).

## Tests and CI

```bash
ARCHI_STAGE=H5 python3 -m pytest -q tests
```

The 54 tests recompute values of the supplementary tables, figures and `outputs/numbers/NUMBERS_R2_35TF.csv` from their
source tables, rerun the supplementary table builders in a temporary copy and check the md5 values in the provenance
files (`tests/README.md`). GitHub Actions (`.github/workflows/tests.yml`) runs them on every push with Python 3.10, R 4.4
and PRROC 1.4.

## Citation

Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine learning
on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.

## Licence

ArChIPelago is distributed under [WTFPL](http://www.wtfpl.net/). If you prefer a standard licence, treat it as CC-BY.
