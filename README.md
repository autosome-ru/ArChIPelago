# ArChIPelago <img src='./Archipelago.png' width='55'>

ArChIPelago (**Arr**angement of multiple **C**hIP-seq and **Hi**gh-throughput analysis using machine **l**earning for prediction of transcription factor binding sites) is a pipeline for training PWM-ensemble Random Forest classifiers that predict transcription factor binding sites from DNA sequence.

This repository reproduces the analysis from _Kravchenko et al., (2026)_.

For scanning new sequences with pre-trained models, see the companion tool:
[ArChIPelago-TFBS-finder](https://github.com/Pavel-Kravchenko/ArChIPelago-TFBS-finder)

---

## Repository Structure

```
ArChIPelago/
│
├── ArChIPelago_code/             # Main analysis notebooks and package
│   ├── config.yml                # ← EDIT THIS FIRST (set your data paths)
│   ├── archipielago/             # Python package (importable in all notebooks)
│   │   ├── __init__.py
│   │   ├── config.py             # Config loading
│   │   ├── io.py                 # FASTA/BED reading, writing, sequence extraction
│   │   ├── scanning.py           # SARUS wrapper, feature matrix construction
│   │   ├── training.py           # RandomForest training and evaluation helpers
│   │   └── scoring.py            # Scorer classes (ROC-AUC, PR-AUC)
│   ├── scorer_module.py          # Backward-compatibility shim (imports archipielago.scoring)
│   ├── tests/                    # Unit tests — run with pytest
│   │   ├── conftest.py
│   │   ├── test_io.py
│   │   ├── test_scanning.py
│   │   ├── test_training.py
│   │   └── test_scoring.py
│   ├── 0_Data_preparation and_test_train_split.ipynb
│   │                             # BED files → FASTA sequences (requires genome + GTRD data)
│   ├── 1_Scanning_with_CHIPMUNK_feature_generation_MONO_DI.ipynb
│   │                             # SARUS PWM scanning → feature matrices
│   ├── 2_ArChIPelago_and_Slim_training.ipynb
│   │                             # Model training (RF, SLIM) and evaluation
│   ├── 3_Biasaway_QC.ipynb       # Negative sequence set quality control
│   ├── 4_Plot_generation_and_analysis.ipynb
│   │                             # Publication figures and metrics
│   └── 5_CTCF_demo_pipeline.ipynb
│                                 # ← START HERE: end-to-end demo on CTCF
│
├── ArChIPelago-TFBS-finder/      # Standalone scanning tool (git submodule)
│   ├── scanning_tool.py          # CLI: scan any FASTA for 36 pre-trained TF models
│   ├── requirements.txt
│   └── tests/
│
├── hocomoco11/                   # HOCOMOCO v11 PWM database
│   ├── models/pwm/mono/          # Mononucleotide PWMs (.pwm)
│   ├── models/pwm/di/            # Dinucleotide PWMs (.dpwm)
│   └── auc/                      # Pre-computed benchmark AUC tables
│
├── sarus/                        # SARUS PWM scanner — Java (git submodule)
│   └── releases/sarus-2.2.3.jar  # Current release
│
├── Slim/                         # SLIM model tools — Java
│   ├── SlimDimont.jar
│   ├── TrainAndApplySlim.jar
│   └── jdk8u232-b09/             # Bundled JDK 8
│
├── seqtk/                        # seqtk sequence toolkit (git submodule)
├── Figures/                      # Publication figures (PDF)
├── Table_1.xlsx                  # ChIP-seq experiment identifiers
├── Table_2.xlsx                  # TF metadata (experiments, peaks, GC content)
└── Table_3.xlsx                  # ArChIPelago performance metrics
```

---

## Installation

### 1. Clone the repository with all submodules

```bash
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago.git
cd ArChIPelago
```

### 2. Create the conda environment

```bash
conda env create -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago
```

> **Note:** R 4.3.1 and rpy2 3.5.11 are included for PRROC-based ROC/PR curve computation. If R is not required, the sklearn-based scorers in `archipielago/scoring.py` can be used instead.

### 3. Verify Java

SARUS requires Java 8 or later:

```bash
java -version       # should print 1.8 or higher
# or use the bundled JDK:
Slim/jdk8u232-b09/bin/java -version
```

### 4. Edit `config.yml`

```bash
nano ArChIPelago_code/config.yml
```

Fill in the paths to your reference genomes, MACS peak directory, and RepeatMasker BED file.

---

## Data Download

### Reference genomes

```bash
# Human (hg38) — ~3 GB
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Mouse (mm10) — ~2.7 GB
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
```

Set `paths.genome_human` and `paths.genome_mouse` in `config.yml`.

### ChIP-seq peaks and pre-processed data

Download from Zenodo: **DOI: 10.5281/zenodo.14927304**

This includes:
- Pre-trained Random Forest models (`Models_sklearn13/`)
- CTCF demo FASTA sequences
- Mononucleotide PWMs (`PWMs_mono_HUMAN/`)
- Dinucleotide PWMs (`PWMs_di_HUMAN/`)

### RepeatMasker

Download `track_out.bed` from the UCSC Table Browser or ENCODE portal.
Set `paths.repeatmasker` in `config.yml`.

---

## Quick Start — CTCF Demo

To see the full pipeline on a single TF without running the full analysis:

```bash
conda activate ArChIPelago
cd ArChIPelago_code
jupyter lab 5_CTCF_demo_pipeline.ipynb
```

This notebook demonstrates:
1. Loading positive (ChIP-seq) and negative (background) FASTA sequences
2. Scanning with SARUS using HOCOMOCO v11 mono + di PWMs
3. Building a feature matrix and selecting top-1000 features
4. Training `RandomForestClassifier` with `random_state=42`
5. Evaluating with ROC-AUC and PR-AUC on held-out test chromosomes
6. Generating publication-quality ROC and PR curve figures

---

## Full Reproduction — Notebooks 0–4

Run notebooks in order. Each notebook reads `config.yml` for external paths.

**For demo/testing, use TF = ["CTCF_HUMAN"]** (or any other TF from the 36 supported).

| Step | Notebook | Input | Output |
|------|----------|-------|--------|
| 1 | `0_Data_preparation...` | GTRD BED + hg38/mm10 FASTA | Train/test FASTA, BED files |
| 2 | `1_Scanning...` | Train/test FASTA + HOCOMOCO PWMs | Per-TF feature matrices |
| 3 | `2_ArChIPelago_and_Slim_training...` | Feature matrices | Trained models, ROC/PR scores |
| 4 | `3_Biasaway_QC...` | Negative FASTA + model predictions | QC plots |
| 5 | `4_Plot_generation...` | Model scores | Publication figures (Figures/) |

> **Compute note:** Full training across all 36 TFs was performed on 100 AMD EPYC 7662 cores (~8 hours).
> For a single TF on a laptop, use `n_jobs=4` and expect ~30 min.

---

## `archipielago` Package API

All reusable functions are consolidated in `ArChIPelago_code/archipielago/`:

```python
from archipielago.config import load_config, get_path
from archipielago.io import fasta_iter, load_fasta, save_fasta, make_train_test_beds
from archipielago.scanning import run_sarus, build_feature_matrix, select_top_features
from archipielago.training import train_rf, evaluate_model, cross_validate_model
from archipielago.scoring import SklearnROCAUC, SklearnPRAUC, PRROC_PRAUC
```

See module docstrings for full parameter documentation.

---

## Running Unit Tests

```bash
conda activate ArChIPelago
cd ArChIPelago_code
pytest tests/ -v
```

Tests do not require Zenodo data downloads or external tools (SARUS, BiasAway).
PRROC-based scorer tests are skipped automatically if R/rpy2 is unavailable.

---

## Known Bugs Fixed

The following bugs were identified during the publication review and fixed in-place
(each fix marked with `# BUGFIX <ID>:` in the code):

| ID | File | Issue |
|----|------|-------|
| C1 | NB0 Cells 20, 37 | Missing comma: `"11" "12"` → `"1112"` silently excluded chr11+chr12 from all train/test sets |
| C2 | NB0 Cells 27–29, 41–42 | `np.random.choice(..., random_state=1)` is invalid (TypeError); replaced with `np.random.seed(1)` |
| C9 | NB2 Cell 8 | `/hocomoco11/` rooted at filesystem root `/`; fixed to `./hocomoco11/` |
| H11| NB0 Cell 21 | `sns.distplot` removed in seaborn ≥0.12; replaced with `sns.histplot(kde=True)` |
| H1 | NB1 Cells 11–12; NB2 Cell 13 | Cyrillic `с` (U+0441) used as variable name; replaced with ASCII `c` |
| H2 | NB2 Cell 13 | PR-di curves used `mtrx_mono` column; **all reported di PR-AUC values were the mono values**; fixed to `mtrx_di` |
| H7 | NB2 Cell 13 | `today_date = "NN"` placeholder; replaced with `date.today().strftime(...)` |
| H9 | scorer_module.py | `__attrs_post_init__` (attrs hook) → `__post_init__` (dataclass hook) |
| H10| scorer_module.py | `PRROC_ROCAUC` was missing `@dataclass` decorator; instantiation raised TypeError |
| C10| scanning_tool.py | SARUS jar hardcoded to 2.0.1; updated to 2.2.3 with 2.0.1 fallback |
| C11| scanning_tool.py | Temp FASTA file leaked on exception in `build_null_distribution`; wrapped in `try/finally` |
| H12| scanning_tool.py | Mononucleotide fallback in `dinucleotide_shuffle` was silent; now logs a warning |
| M5 | scanning_tool.py | BED export used `split('@')[0]`; changed to `rsplit('@', 1)` for ENCODE-style seq IDs |

---

## Platform

Developed and tested on Ubuntu 20.04.6 LTS (GNU/Linux 5.15.0-113-generic x86_64).

---

## License

ArChIPelago is distributed under WTFPL. If you prefer more standard licenses, treat it as CC-BY.

---

## Citation

TBA. Kravchenko et al., (2026)
