# ArChIPelago <img src='./Archipelago.png' width='55'>

[![tests](https://github.com/autosome-ru/ArChIPelago/actions/workflows/tests.yml/badge.svg)](https://github.com/autosome-ru/ArChIPelago/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14927303.svg)](https://doi.org/10.5281/zenodo.14927303)
[![release](https://img.shields.io/github/v/tag/autosome-ru/ArChIPelago?label=release)](https://github.com/autosome-ru/ArChIPelago/tree/v2.0.0)
[![License: WTFPL](https://img.shields.io/badge/license-WTFPL-lightgrey.svg)](http://www.wtfpl.net/)

**Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites**

ArChIPelago combines multiple position weight matrices (PWMs) of one transcription factor (TF) into a joint model with
classic machine learning, from logistic regression to ensembles of decision trees, to predict TF binding sites in genomic
sequences.

With 637 ChIP-Seq datasets (361 human, 276 mouse) for 35 orthologous human and mouse TFs and 1,447 monoPWMs and 754 diPWMs
(human; mouse 1,178 and 617) curated for HOCOMOCO v11, ArChIPelago outperforms the best single PWM of every TF and scores
higher than Slim and diChIPMunk models for most TFs. Models trained on human data transfer to the mouse orthologs.

> Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). *Classic machine
> learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.*

This repository holds the code that built the data, trained and evaluated the models and produced every table, figure and
number of the manuscript, together with the inputs and outputs of the final run. The data and models are on Zenodo
([10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303), all versions). The command-line tool
[ArChIPelago-TFBS-finder](https://github.com/autosome-ru/ArChIPelago-TFBS-finder) (git submodule) scans a FASTA file
with the trained models, human-trained or mouse-trained.

Version 2.0.0 is the revised analysis. The first version (36 TFs, the notebooks and `Manuscript_analysis/`) is kept at
tag [v1.0.0](https://github.com/autosome-ru/ArChIPelago/tree/v1.0.0).

---

## Results (Random Forest on monoPWMs and diPWMs, 35 TFs)

| test set | best single monoPWM | ArChIPelago | TFs improved |
|---|---|---|---|
| human, auROC (median) | 0.874 | 0.902 | 35 of 35 |
| human, auPRC (median) | 0.288 | 0.336 | 35 of 35 |
| mouse, human-trained models, auROC | 0.857 | 0.886 | 34 of 35 |
| mouse, human-trained models, auPRC | 0.230 | 0.307 | 34 of 35 |
| mouse, mouse-trained models, auROC | 0.873 | 0.897 | 35 of 35 |
| mouse, mouse-trained models, auPRC | 0.248 | 0.356 | 35 of 35 |

At 50 % recall on the human test set, ArChIPelago gives 30 % fewer false positives than the best single monoPWM (median
over TFs; 33 of 35 TFs). The baseline is the best single monoPWM of each TF, selected on the human training set by the
metric of the row. Human rows are random seed 0, mouse rows the mean of seeds 0 and 1 (Table 1). Every value is in `outputs/numbers/NUMBERS_R2_35TF.csv` with its source table.

---

## Repository structure

```
ArChIPelago/
├── 1_data_preparation/     positives, GC-matched negatives, SPRY-SARUS feature tables, Slim inputs, data checks (server)
├── 2_training/             Random Forest, logistic regression, XGBoost and bagging models, evaluation tables,
│                           Slim and diChIPMunk models (server), release models
├── 3_analysis/
│   ├── server/             operational and low-FPR metrics, top peaks, random split, feature importance, saturation,
│   │                       Fig. 4 refits, runtime benchmark (server)
│   ├── integration/        local driver run_h3.py and the producers of the integrated tables and numbers
│   ├── package/            results table, saturation, cross-species tables, figure source data
│   ├── low_fpr_metrics/  top_peaks/  random_split/  feature_importance/  transfer_drivers/
│   └── supplementary_tables/   Sup. Tables 1-10
├── 4_visualisation/        producers of Fig. 1-4 and Fig. S1-S12
├── results/                tracked inputs of the local stages (mirror of the server results, dataset lists, PWM archives)
├── outputs/                figures/, supplementary_tables/, numbers/ of the final run
├── tests/                  checks of the tables, figures and numbers (pytest)
├── setup_server_layout.sh  links the repository folders into the folder layout the server scripts expect
├── Kravchenko_et_al_Supplementary_Materials.pdf   table and figure legends, Fig. S1-S12
├── ArChIPelago_code/       archipielago Python package, its tests and the CTCF demo notebook
├── ArChIPelago-TFBS-finder/  scanning tool (git submodule)
├── sarus/  seqtk/          SPRY-SARUS and seqtk (git submodules)
└── Slim/                   Slim jars (SlimDimont, TrainAndApplySlim), the ytilib Ruby library of diChIPMunk, JDK 8
```

Each numbered folder, `results/`, `outputs/` and `tests/` has its own README with the run commands, inputs and outputs.

---

## How the analysis runs

1. **Data preparation** (`1_data_preparation/`, server): for every GTRD MACS peak set the top 1,000 peaks (tags >= 10,
   ranked by -10log10 P), +-500 bp windows around the summits, removal of windows overlapping RepeatMasker repeats by more
   than 70 %, at most 10,000 positives per TF and set, 300 bp sequences around the summit. Negatives are peaks of TFs of
   other TFClass classes, without windows within 300 bp of the TF's own positives, GC-matched with BiasAway. Features are
   SPRY-SARUS best-hit log-odds scores of every PWM of the TF.
2. **Training** (`2_training/`, server): `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100)` and the
   other algorithms on the human and the mouse training sets; every matrix is standardised with the mean and standard
   deviation of the training matrix. Slim (m = 0, 1, -5) and diChIPMunk models for the comparison of Fig. 4.
3. **Server analyses** (`3_analysis/server/`), then a mirror of their result tables to `results/server/`.
4. **Local analysis and visualisation**: one command rebuilds every table, figure and number from `results/`:

   ```bash
   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
   ARCHI_STAGE=H5 ARCHI_H5_SCORES=$PWD/results/server/human/scores \
   python3 3_analysis/integration/run_h3.py all
   ```

   It checks every input against the md5 manifest of the server mirror first and stops at the first failing command
   (about 3 min on a laptop). Requirements: Python 3.9 or later with pandas 2.3, numpy 1.26, scipy 1.13, matplotlib 3.5,
   scikit-learn, openpyxl, pyyaml, pypdf; R 4.4 with dplyr, ggplot2, patchwork, cowplot, ggrepel, ggExtra, viridis,
   tidyr, ggbeeswarm, ggseqlogo and png; poppler (pdftoppm, pdfimages, pdftotext) and rsvg-convert. Stage `fig1` sets the
   resolution of Fig. 1 with the macOS tool `sips`; on Linux that command fails and the run stops there (run the other
   stages by name, see `3_analysis/README.md`).

**Training and test chromosomes.** Human: training chr2-7, 9-20 and 22, test chr1, 8 and 21. Mouse: training chr2-7 and
9-18, test chr1, 8 and 19. Each autosome is in exactly one set; the sex chromosomes are in neither.

**Server scripts** (steps 1-3) ran on the servers of the authors and keep their paths (`~/archi_work/...`);
`setup_server_layout.sh <repo> <base>` recreates that layout as symbolic links. Their software: Python 3.8.18,
scikit-learn 1.3.0, pandas 2.0.3, numpy 1.24.3, xgboost 1.7.3, pybedtools 0.9.1, bedtools 2.31.1, BiasAway 3.3.0,
SPRY-SARUS 2.0.2 (`2_training/VERSIONS.txt`). External data not in this repository: reference genomes hg38 and mm10, GTRD
MACS peak files, HOCOMOCO v11 PWMs (Zenodo archives `macs.tar.gz`, `hocomoco11.tar.gz`, `PWMs_*.tar.gz`) and the
RepeatMasker masks (`repeatmasker_hg38_mm10.tar.gz`).

---

## Supported TFs (35)

ANDR, AP2A, CEBPB, COE1, CTCF, E2F4, ERG, ESR1, FLI1, GATA1, GATA2, GATA3, GCR, HNF4A, IRF1, IRF4, JUND, MAFK, MAX, MYC,
PPARG, PRGR, REST, RUNX1, RXRA, SOX2, SPI1, SRF, STA5A, STAT1, STAT3, TAL1, TF65, TFE2, USF2.

---

## Data on Zenodo

Record [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303) (resolves to the latest version; the version of
v1.0.0 is [10.5281/zenodo.23036475](https://doi.org/10.5281/zenodo.23036475)). Version 2.0.0 holds the training and test
sets (`train.tar.gz`, `test.tar.gz`), the PWMs (`PWMs_{mono,di}_{HUMAN,MOUSE}.tar.gz`), the global sequence tables and
FASTA files, the negatives, the 210 Random Forest models (`Models.tar.gz`: 35 TFs, human-trained and mouse-trained, on
monoPWMs, diPWMs or both, each with a `.json` of feature order and training mean and standard deviation), the Slim and
diChIPMunk models, the RepeatMasker masks and the metadata. Its `README_v2.md` describes every archive. In the Zenodo
files, `train` always marks the training chromosomes. The loader `1_data_preparation/archi_data_v5.py` reads the build
folders, where the mouse files named `train` hold the test chromosomes; on the Zenodo files its `_file_split()` must
return the split unchanged for MOUSE (see `1_data_preparation/README.md`).

---

## `archipielago` package and CTCF demo

`ArChIPelago_code/archipielago/` holds general functions for training on PWM scores (configuration, FASTA and BED input,
SPRY-SARUS scanning and feature matrices, Random Forest training and evaluation, scikit-learn and PRROC scorers), with unit
tests in `ArChIPelago_code/tests/`. `ArChIPelago_code/5_CTCF_demo_pipeline.ipynb` is the demo of v1.0.0: scanning,
training and evaluation for CTCF. Besides the positives of `train.tar.gz` and `test.tar.gz` it needs negative FASTA files
(`train/CTCF_negative_train.fasta`, `test/CTCF_negative_test.fasta`), which are built from the negative ID lists of
`Archipelago_intermediate_files_negatives_v3.tar.gz` and the sequences of the `Archipelago_intermediate_files_fasta_*`
archives (paths in `ArChIPelago_code/config.yml`, conda environment `ArChIPelago_code/environment_rpy_2.yml`).

```bash
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago.git
cd ArChIPelago
conda env create -n ArChIPelago -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago
pytest ArChIPelago_code/tests
```

---

## Tests

```bash
pytest ArChIPelago_code/tests                                   # archipielago package
ARCHI_STAGE=H5 pytest tests --ignore=tests/server_checks       # tables, figures and numbers of the manuscript
```

The second suite recomputes values of the supplementary tables, the figures and `outputs/numbers/NUMBERS_R2_35TF.csv`
from their source tables, reruns the supplementary table builders in a temporary copy, and checks the md5 values recorded
in the provenance files. GitHub Actions runs both suites on every push.

---

## License

ArChIPelago is distributed under [WTFPL](http://www.wtfpl.net/). If you prefer a standard license, treat it as CC-BY.

## Citation

Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine learning
on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.
