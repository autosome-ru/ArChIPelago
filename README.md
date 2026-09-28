# ArChIPelago <img src='./Archipelago.png' width='55'>

**Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites**

ArChIPelago is a computational framework that combines multiple position weight matrices (PWMs) into a joint model using classic machine learning techniques, from linear regression to ensembles of decision trees, to improve prediction of transcription factor binding sites in genomic sequences.

Using 704 ChIP-Seq datasets (414 human, 290 mouse) for 36 orthologous human-mouse transcription factors and 2,275 human and 1,816 mouse mono- and dinucleotide PWMs curated for the HOCOMOCO v11 motif collection, we show that ArChIPelago consistently outperforms the best available individual PWMs as well as sparse local inhomogeneous mixture (Slim) models. Models trained on human data transfer to the mouse orthologs.

This repository contains the analysis pipeline, the code and outputs behind every number, supplementary table and figure panel of the manuscript, and the code to reproduce them:

> Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). *Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.*

> **Want to scan your own sequences?** See the companion tool [ArChIPelago-TFBS-finder](https://github.com/autosome-ru/ArChIPelago-TFBS-finder), a command-line tool that applies the pre-trained models to any FASTA file.

---

## Overview

<table>
<tr>
<td width="55%">

**Pipeline at a glance:**

1. Extract ChIP-Seq peaks from GTRD, split by chromosome into training and test sets
2. Scan the peak regions with HOCOMOCO v11 PWMs (mono- and dinucleotide) using SPRY-SARUS
3. Build feature matrices from the log-odds best-hit PWM scores
4. Train classifiers (Random Forest, XGBoost, Logistic Regression, Bagging) on the human training set
5. Evaluate with auROC and auPRC (PRROC R package) on the human and mouse test sets
6. Compare with the best single monoPWM, Slim models and diChIPMunk diPWMs

</td>
<td width="45%">

**Key results (Random Forest on monoPWMs + diPWMs, 36 TFs):**

- Median 34.5 human PWMs per TF (up to 303)
- Human test set: median auROC **0.891** vs 0.861 for the best single monoPWM; auPRC **0.298** vs 0.255
- Median per-TF gain +0.024 auROC, +0.041 auPRC; all 36 TFs improved on both metrics
- 24 % fewer false positives at 50 % recall (median over TFs)
- Mouse test set (human-trained models): auROC 0.874 vs 0.861, auPRC 0.248 vs 0.183

</td>
</tr>
</table>

---

## Repository Structure

```
ArChIPelago/
│
├── ArChIPelago_code/                # Pipeline: data preparation, scanning, training
│   ├── config.yml                   # ← EDIT THIS: set paths to your data
│   ├── archipielago/                # Python package (used by all notebooks)
│   │   ├── config.py                # Config loading & path validation
│   │   ├── io.py                    # FASTA/BED I/O, train/test splitting
│   │   ├── scanning.py              # SPRY-SARUS wrapper, feature matrix construction
│   │   ├── training.py              # RF training, evaluation, cross-validation
│   │   └── scoring.py               # Scorer classes (sklearn + PRROC)
│   ├── tests/                       # Unit tests (pytest)
│   ├── 5_CTCF_demo_pipeline.ipynb   # ← START HERE: end-to-end CTCF demo
│   ├── 0_Data_preparation and_test_train_split.ipynb
│   ├── 1_Scanning_with_CHIPMUNK_feature_generation_MONO_DI.ipynb
│   ├── 2_ArChIPelago_and_Slim_training.ipynb
│   ├── 3_Biasaway_QC.ipynb
│   ├── 4_Plot_generation_and_analysis.ipynb
│   ├── run_simulation.py            # End-to-end pipeline check on simulated data
│   ├── scorer_module.py             # Backward-compatible import of archipielago.scoring
│   └── environment_rpy_2.yml        # Conda environment specification
│
├── Manuscript_analysis/             # Numbers, tables and figures of the manuscript
│   ├── results_table.csv            # Results table (36 TFs x 5 algorithms x 3 PWM sets)
│   ├── numbers_in_text*.json        # Every number quoted in the text
│   ├── Sup_Tables/                  # Supplementary Tables 1-6 (xlsx + csv)
│   ├── Figures/                     # Figure panels (PDF) and their source data (csv)
│   ├── scripts/                     # Table, source-data and figure builders; rebuild_tables_and_figures.sh
│   └── analysis/                    # Per-TF analyses and their outputs
│
├── Kravchenko_et_al_Supplementary_Materials.pdf   # Supplementary figures S1-S7 and table legends
│
├── ArChIPelago-TFBS-finder/         # Standalone scanning tool (git submodule)
├── sarus/                           # SPRY-SARUS PWM scanner (git submodule)
├── seqtk/                           # seqtk sequence toolkit (git submodule)
│
└── Slim/                            # Slim / diChIPMunk tools (Java)
    ├── SlimDimont.jar               # diChIPMunk: de novo diPWM discovery
    ├── TrainAndApplySlim.jar        # Slim: sparse local inhomogeneous mixture models
    ├── ytilib/                      # Ruby utilities for diChIPMunk
    └── jdk8u232-b09/                # Bundled OpenJDK 8 (required by the Slim jars)
```

---

## Quick Start: CTCF Demo

Run the full pipeline on CTCF without any external data preprocessing:

```bash
# 1. Clone
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago.git
cd ArChIPelago

# 2. Set up environment
conda env create -n ArChIPelago -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago

# 3. Download CTCF data from Zenodo (DOI: 10.5281/zenodo.14927303)
#    Extract train/, test/, PWMs_mono_HUMAN/, PWMs_di_HUMAN/
#    Edit ArChIPelago_code/config.yml with the paths

# 4. Run the demo
cd ArChIPelago_code
jupyter lab 5_CTCF_demo_pipeline.ipynb
```

The demo notebook walks through every step:

| Step | What it does |
|------|-------------|
| Load sequences | Positive (ChIP-Seq peaks) and negative (GC-matched background) FASTA |
| PWM scanning | SPRY-SARUS scan with HOCOMOCO v11 mono + di PWMs |
| Feature matrix | Build and select top-1000 features by RF importance |
| Train model | `RandomForestClassifier(n_estimators=100, max_depth=6)` |
| Evaluate | auROC, auPRC on held-out test chromosomes |
| Predict | Per-sequence binding probabilities |
| Visualize | ROC and PR curves |

---

## Installation

### Prerequisites

- **Conda** (Miniconda or Anaconda)
- **Java 8+** (for SPRY-SARUS PWM scanning; OpenJDK 8 is bundled in `Slim/jdk8u232-b09/`)
- **Git** with submodule support

### Step-by-step

```bash
# Clone with all submodules (sarus, seqtk, ArChIPelago-TFBS-finder)
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago.git
cd ArChIPelago

# Create and activate the conda environment
conda env create -n ArChIPelago -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago

# Verify Java
java -version   # should print 1.8 or higher
# or use the bundled JDK:
Slim/jdk8u232-b09/bin/java -version

# Edit configuration
nano ArChIPelago_code/config.yml
```

> **Note on R:** The environment includes R 4.3.1 and rpy2 3.5.11 for PRROC-based auROC/auPRC computation (Grau et al. 2015), as used in the manuscript. If R is not needed, the sklearn-based scorers in `archipielago/scoring.py` provide equivalent functionality.

### Configuration

Edit `ArChIPelago_code/config.yml` to set paths for your system:

```yaml
paths:
  genome_human: "/path/to/hg38.fa"                  # UCSC hg38 reference genome
  genome_mouse: "/path/to/mm10.fa"                  # UCSC mm10 reference genome
  repeatmasker: "/path/to/repeatmasker/track_out.bed"
  macs_peaks:   "/path/to/macs/"                    # GTRD MACS peak interval files
  output_dir:   "Release/TF-ML"                     # output directory of notebooks 0-2
  hocomoco11_wlogauc_mono: "/path/to/hocomoco11/wlogauc/mono"   # notebook 1, PWM selection
  hocomoco11_wlogauc_di:   "/path/to/hocomoco11/wlogauc/di"
  pcms_mono: "/path/to/pcms_for_logo_mono"          # notebook 1, logos
  pcms_di:   "/path/to/pcms_for_logo_di"

tools:
  sarus_jar:       "../sarus/releases/sarus-2.2.3.jar"
  java_bin:        "java"
  slim_dimont_jar: "../Slim/SlimDimont.jar"         # diChIPMunk
  slim_apply_jar:  "../Slim/TrainAndApplySlim.jar"  # Slim model training
  slim_java_bin:   "../Slim/jdk8u232-b09/bin/java"  # JDK 8 for the Slim jars
```

---

## Data Download

### Zenodo archive (required)

**DOI: [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303)** (resolves to the latest version of the archive)

| Archive | Contents | Required for |
|---------|----------|-------------|
| `train.tar.gz` | Training peak sets (BED) and sequences (FASTA) of the 36 TFs, human and mouse | Notebooks 2, 5 |
| `test.tar.gz` | Test peak sets (BED) and sequences (FASTA) of the 36 TFs, human and mouse | Notebooks 2, 5 |
| `PWMs_mono_HUMAN.tar.gz`, `PWMs_di_HUMAN.tar.gz` | Human monoPWMs (1,495) and diPWMs (780) of the 36 TFs | Notebooks 1, 2, 5 |
| `PWMs_mono_MOUSE.tar.gz`, `PWMs_di_MOUSE.tar.gz` | Mouse monoPWMs (1,192) and diPWMs (624) of the 36 TFs | Notebooks 1, 2; mouse-trained control |
| `hocomoco11.tar.gz` | HOCOMOCO v11 models, benchmark data and curation tables | Notebook 0 |
| `macs.tar.gz` | GTRD MACS peak intervals | Notebook 0 |
| `Archipelago_intermediate_files.tar.gz` | Global sequence tables and FASTA files of the training and test sets, result tables of notebook 2 | Skip notebook 0 |
| `Models.tar.gz` | Pre-trained Random Forest models of the 36 TFs on all their human PWMs (monoPWMs, diPWMs, both; scikit-learn 1.3), each with a `.json` feature specification (PWM order, training mean and standard deviation) | ArChIPelago-TFBS-finder |
| `Slim.tar.gz` | Slim jar files and bundled JDK 8 | Notebook 2 (Slim training) |
| `Slim_models.tar.gz` | Pre-trained Slim models | Notebook 2 (comparison) |
| `Manuscript_analysis.tar.gz` | `Manuscript_analysis/` of this repository (results table, Sup. Tables 1-6, figure panels and source data, scripts, per-TF analysis outputs) and the supplementary figures PDF | Tables and figures of the manuscript |
| `ArChIPelago_code_<commit>.tar.gz` | Snapshot of this repository at the commit named in the file name (without the submodules) | Code |

In every archive, files named `train` hold the training chromosomes and files named `test` or `control` hold the test chromosomes (see the split below).

### Reference genomes (for Notebook 0 only)

```bash
# Human hg38 (~3 GB)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz

# Mouse mm10 (~2.7 GB)
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz && gunzip mm10.fa.gz
```

### RepeatMasker (for Notebook 0 only)

Download `track_out.bed` from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) (group: Repeats, track: RepeatMasker) and set `paths.repeatmasker` in `config.yml`.

---

## Full Reproduction: Notebooks 0-4 and `Manuscript_analysis/`

Run the notebooks sequentially. Each reads `config.yml` for external paths.

| Step | Notebook | Input | Output |
|------|----------|-------|--------|
| 0 | Data preparation | GTRD BED + hg38/mm10 FASTA | Training and test FASTA and BED files |
| 1 | PWM scanning | FASTA + HOCOMOCO PWMs | Per-TF SPRY-SARUS feature matrices |
| 2 | Model training | Feature matrices | Trained RF and Slim models, auROC/auPRC scores |
| 3 | BiasAway QC | Negative sequences | GC-content quality control plots |
| 4 | Results collection | Model scores | Results table of notebook 2 (`HUMAN_MOUSE_total_100k.csv`; kept as `Manuscript_analysis/analysis/inputs/notebook4_results_table.csv`) |

The numbers, supplementary tables and figure panels of the manuscript are computed from the outputs of notebooks 0-2 by the scripts in [`Manuscript_analysis/`](Manuscript_analysis/README.md): the evaluation of every model and baseline on both test sets, the mouse-trained control, the operational metrics, the dependence on the number of PWMs, the Fig. 4 models on the best monoPWM and diPWM, the cross-species table and the motif-subtype figure. `bash Manuscript_analysis/scripts/rebuild_tables_and_figures.sh` rebuilds all tables and panels from the per-TF outputs in `Manuscript_analysis/analysis/`.

**Data preparation (Notebook 0):** ChIP-Seq peaks from GTRD were called with MACS (Zhang et al. 2008). Peak lists were filtered by `tags >= 10` and sorted by `-log10(P value)`. Putative binding regions [-150;+150] were extracted centred at the peak summit. Repeat-overlapping peaks were removed using RepeatMasker (Smit et al. 2013-2015) via pybedtools `subtract(f=0.7, N=True)`.

**Training and test sets:** following Zhou and Troyanskaya (2015), chromosomes 1, 8 and 21 (human) and 1, 8 and 19 (mouse) are the test sets; the training sets are chromosomes 2-7, 9, 10 and 13-20 (human) and 2-7, 9, 10 and 13-18 (mouse). Chromosomes 11, 12 and the sex chromosomes are used in neither set. Up to 10,000 positives per TF; negatives (1:100 class balance) sampled from peaks of unrelated TF families (TFClass; Wingender et al. 2018), GC-matched using BiasAway (Khan et al. 2021).

**Baseline:** the best single monoPWM of each TF, selected on the human training set, independently by auROC and by auPRC.

**PWM scanning:** Genomic regions were scanned with SPRY-SARUS (Kulakovskiy et al. 2016) using `--skipn --show-non-matching --output-scoring-mode score besthit`. Log-odds best-hit scores serve as features.

**Model parameters (Random Forest):** `max_depth=6, max_samples=0.8, n_estimators=100` (selected via GridSearchCV). Features were scale-transformed with `sklearn.preprocessing.StandardScaler`.

**Slim models (Notebook 2):** ArChIPelago was benchmarked against sparse local inhomogeneous mixture (Slim) models (Grau et al. 2013), trained side-by-side from extended 1001 bp genomic regions around the same peak summits. Slim models were trained using `TrainAndApplySlim.jar` with the bundled JDK 8 (`Slim/jdk8u232-b09/bin/java`). Three Slim model orders were compared: `markov_order=0` (equivalent to monoPWM), `markov_order=1` (equivalent to diPWM), and LSlim with `markov_order=-5` (limited Slim; Keilwagen and Grau 2015). Peak signal annotations were derived from the `-10*log10(pvalue)` MACS output field. Slim predictions used `max_score` for performance assessment. Additionally, diPWMs were constructed *de novo* from the positive sequences using diChIPMunk (`run_dichiphorde8.rb`). The Random Forest on the best monoPWM and the best diPWM (RF2f) was augmented with Slim and diChIPMunk features to test whether combining diverse model types improves prediction (Fig. 4, `Manuscript_analysis/analysis/slim_dichipmunk/`).

**Supported TFs (36):** ANDR, AP2A, CEBPB, COE1, CTCF, E2F4, ERG, ESR1, FLI1, GATA1, GATA2, GATA3, GCR, HNF4A, IRF1, IRF4, JUND, MAFK, MAX, MYC, P53, PPARG, PRGR, REST, RUNX1, RXRA, SOX2, SPI1, SRF, STA5A, STAT1, STAT3, TAL1, TF65, TFE2, USF2.

> **Compute requirements (Sup. Table 6):** with the Java virtual machine limited to one core of a 2.6 GHz Intel Xeon E5-4607 v2, SPRY-SARUS needs 0.15 s per monoPWM (0.19 s per diPWM) per megabase of sequence; training the Random Forest on a full human training matrix (200,000-230,000 sequences) takes 28-104 s on one core with at most 1.9 GB RAM, and prediction for 10⁶ sequences 7-20 s.

---

## `archipielago` Python Package

All reusable pipeline functions are consolidated in `ArChIPelago_code/archipielago/`:

```python
# Configuration
from archipielago.config import load_config, get_path

# I/O
from archipielago.io import fasta_iter, load_fasta, save_fasta, make_train_test_beds

# PWM scanning and feature construction
from archipielago.scanning import run_sarus, build_feature_matrix, select_top_features

# Model training and evaluation
from archipielago.training import train_rf, evaluate_model, cross_validate_model

# Scorer classes (sklearn and R/PRROC)
from archipielago.scoring import SklearnROCAUC, SklearnPRAUC, PRROC_PRAUC, PRROC_ROCAUC
```

### Minimal example

```python
from archipielago import io, scanning, training
from archipielago.config import load_config

cfg = load_config()

# Load sequences
train_pos = io.load_fasta("train/CTCF_HUMAN.PEAKS000001.macs.train.mfa")

# Scan with SPRY-SARUS
scanning.run_sarus("train.fasta", "motif.pwm", cfg['tools']['sarus_jar'],
                   "scores.txt", pwm_type="mono")

# Build features and train
X = scanning.build_feature_matrix("scan_dir/", mode="mono_di")
model = training.train_rf(X, labels, n_estimators=100, random_state=42)

# Evaluate
results = training.evaluate_model(model, X_test, y_test)
print(f"auROC: {results['roc_auc']:.4f}, auPRC: {results['pr_auc']:.4f}")
```

---

## Running Tests

```bash
conda activate ArChIPelago
cd ArChIPelago_code
pytest tests/ -v
```

The test suite covers all package modules (59 tests). Tests run without Zenodo data, external tools, or reference genomes. PRROC scorer tests are skipped automatically when R/rpy2 is unavailable.

---

## Manuscript Figures and Tables

Every figure panel and supplementary table is produced by a script in `Manuscript_analysis/`; the figure PDFs are in `Manuscript_analysis/Figures/panels/`, the plotted values in `Manuscript_analysis/Figures/source_data/`.

| Item | Script (`Manuscript_analysis/scripts/`) | Data |
|------|--------|------|
| Fig. 2 / S1 | `Figure_2_human_test.R`, `Figure_S1_mouse_test.R` | `results_table.csv` |
| Fig. S2 | `Figure_S2_RF_vs_best_monoPWM.py` | `results_table.csv` |
| Fig. 3A-C / S3 | `Figure_3ABC_human_test.R`, `Figure_S3ABC_mouse_test.R` | `results_table.csv` |
| Fig. 4 / S4 | `Figure_4_and_S4.R` | `analysis/slim_dichipmunk/` |
| Fig. S5 | `Figure_S5_saturation.R` | `analysis/saturation/` |
| Fig. S6 | `Figure_S6_cross_species.R` | Sup. Table 5 |
| Fig. S7 | `Figure_S7_motif_subtypes.R` | `analysis/motif_subtypes/` |
| Sup. Tables 1-2 | `make_sup_tables_1_2.py` | `analysis/inputs/`, `analysis/datasets/`, `analysis/mouse_transfer/` |
| Sup. Tables 3-6 | `make_sup_tables_3_to_6.py` | results table and `analysis/` outputs |

Supplementary Tables: 1, ChIP-Seq experiments; 2, PWMs and datasets per TF; 3, ArChIPelago performance; 4, operational metrics; 5, cross-species transfer and the mouse-trained control; 6, runtime and memory. Their legends and the supplementary figures are in `Kravchenko_et_al_Supplementary_Materials.pdf`.

---

## Platform

Developed and tested on:
- **Ubuntu 20.04.6 LTS** (GNU/Linux 5.15.0-113-generic x86_64)
- **Python 3.8.18**, scikit-learn 1.3.0, numpy 1.24.3, pandas 2.0.3, xgboost 1.7.3 (pinned in `ArChIPelago_code/environment_rpy_2.yml`)
- **R 4.3.1** with PRROC, rpy2 3.5.11
- **Java 8** (OpenJDK 8u232-b09, bundled)

The tables and figure panels of `Manuscript_analysis/` were built with Python 3.9 (pandas 2.3, numpy 1.26) and R 4.4.3.

---

## License

ArChIPelago is distributed under [WTFPL](http://www.wtfpl.net/). If you prefer a standard license, treat it as CC-BY.

---

## Citation

Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.

Zenodo data archive: [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303)
