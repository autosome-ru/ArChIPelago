# ArChIPelago <img src='./Archipelago.png' width='55'>

**Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites**

ArChIPelago is a computational framework that combines multiple position weight matrices (PWMs) into a joint model using classic machine learning techniques — from linear regression to ensembles of decision trees — to improve prediction of transcription factor binding sites in genomic sequences.

Using 704 ChIP-Seq datasets spanning 36 orthologous human-mouse transcription factor pairs and 1558 mono- and dinucleotide PWMs from the HOCOMOCO v11 motif collection, we show that ArChIPelago consistently outperforms the best available individual PWMs as well as sparse local inhomogeneous mixture (Slim) models. Models trained on human data generalize to mouse orthologs, demonstrating cross-species applicability.

This repository contains the full analysis pipeline, pre-trained models, and code to reproduce all results from:

> Kravchenko P., Vorontsov I.E., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). *Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites.* Bioinformatics (ECCB 2026).

> **Want to scan your own sequences?** See the companion tool [ArChIPelago-TFBS-finder](https://github.com/Pavel-Kravchenko/ArChIPelago-TFBS-finder) for a ready-to-use CLI that applies pre-trained models to any FASTA file.

---

## Overview

<table>
<tr>
<td width="55%">

**Pipeline at a glance:**

1. Extract ChIP-seq peaks from GTRD, split by chromosome into train/test sets
2. Scan 301 bp peak regions with HOCOMOCO v11 PWMs (mono- and dinucleotide) using SPRY-SARUS
3. Build feature matrices from log-odds best-hit PWM scores
4. Train ensemble classifiers (Random Forest, XGBoost, Logistic Regression, Bagging)
5. Evaluate with auROC and auPRC using the PRROC R package
6. Compare against single-PWM baselines, Slim models, and diChIPMunk

</td>
<td width="45%">

**Key results:**

- 704 ChIP-Seq datasets, 36 human-mouse TF pairs
- Median 17 mono + 17 di PWMs per TF (up to 198+99)
- Median auROC **0.891** (vs 0.869 single best PWM)
- Median auPRC improvement **+0.043** over best PWM
- Cross-species validation on mouse orthologs
- Outperforms single best PWM for 34 of 36 TFs

</td>
</tr>
</table>

---

## Repository Structure

```
ArChIPelago/
│
├── ArChIPelago_code/                # Main analysis pipeline
│   ├── config.yml                   # ← EDIT THIS: set paths to your data
│   │
│   ├── archipielago/                # Python package (used by all notebooks)
│   │   ├── config.py                # Config loading & path validation
│   │   ├── io.py                    # FASTA/BED I/O, train/test splitting
│   │   ├── scanning.py              # SPRY-SARUS wrapper, feature matrix construction
│   │   ├── training.py              # RF training, evaluation, cross-validation
│   │   └── scoring.py               # Scorer classes (sklearn + PRROC)
│   │
│   ├── tests/                       # Unit tests (pytest)
│   │   ├── conftest.py              # Shared fixtures
│   │   ├── test_config.py
│   │   ├── test_io.py
│   │   ├── test_scanning.py
│   │   ├── test_training.py
│   │   └── test_scoring.py
│   │
│   ├── 5_CTCF_demo_pipeline.ipynb   # ← START HERE: end-to-end CTCF demo
│   ├── 0_Data_preparation and_test_train_split.ipynb
│   ├── 1_Scanning_with_CHIPMUNK_feature_generation_MONO_DI.ipynb
│   ├── 2_ArChIPelago_and_Slim_training.ipynb
│   ├── 3_Biasaway_QC.ipynb
│   ├── 4_Plot_generation_and_analysis.ipynb
│   ├── run_simulation.py            # End-to-end pipeline verification script
│   ├── scorer_module.py             # Backward-compat shim → archipielago.scoring
│   └── environment_rpy_2.yml        # Conda environment specification
│
├── Figures/                         # Publication figure R scripts
│   ├── Figure 2 and S1/             # Model comparison boxplots (Fig. 2, S1)
│   ├── Figure 3 and S3/             # PWM vs ensemble scatter + delta (Fig. 3, S3)
│   └── Figure 4/                    # Slim/diChIPMunk comparison (Fig. 4, S4)
│
├── ArChIPelago-TFBS-finder/         # Standalone scanning tool (git submodule)
├── sarus/                           # SPRY-SARUS PWM scanner (git submodule)
├── seqtk/                           # seqtk sequence toolkit (git submodule)
│
├── Slim/                            # Slim model tools (Java)
│   ├── SlimDimont.jar
│   ├── TrainAndApplySlim.jar
│   └── jdk8u232-b09/                # Bundled OpenJDK 8 for Slim
│
├── Table_1.xlsx                     # ChIP-seq experiment identifiers (Sup. Table 1)
├── Table_2.xlsx                     # TF metadata: experiments, peaks, GC% (Sup. Table 2)
└── Table_3.xlsx                     # ArChIPelago performance metrics (Sup. Table 3)
```

---

## Quick Start — CTCF Demo

Run the full pipeline on CTCF without any external data preprocessing:

```bash
# 1. Clone
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago.git
cd ArChIPelago

# 2. Set up environment
conda env create -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago

# 3. Download CTCF data from Zenodo (DOI: 10.5281/zenodo.14927304)
#    Extract train/, test/, PWMs_mono_HUMAN/, PWMs_di_HUMAN/
#    Edit ArChIPelago_code/config.yml with the paths

# 4. Run the demo
cd ArChIPelago_code
jupyter lab 5_CTCF_demo_pipeline.ipynb
```

The demo notebook walks through every step:

| Step | What it does |
|------|-------------|
| Load sequences | Positive (ChIP-seq peaks) and negative (GC-matched background) FASTA |
| PWM scanning | SPRY-SARUS scan with HOCOMOCO v11 mono + di PWMs |
| Feature matrix | Build and select top-1000 features by RF importance |
| Train model | `RandomForestClassifier(n_estimators=100, max_depth=6)` |
| Evaluate | auROC, auPRC on held-out test chromosomes |
| Predict | Per-sequence binding probabilities |
| Visualize | Publication-quality ROC and PR curves |

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
conda env create -f ArChIPelago_code/environment_rpy_2.yml
conda activate ArChIPelago

# Verify Java
java -version   # should print 1.8 or higher
# or use the bundled JDK:
Slim/jdk8u232-b09/bin/java -version

# Edit configuration
nano ArChIPelago_code/config.yml
```

> **Note on R:** The environment includes R 4.3.1 and rpy2 3.5.11 for PRROC-based auROC/auPRC computation (Grau et al. 2015), as used in the publication. If R is not needed, the sklearn-based scorers in `archipielago/scoring.py` provide equivalent functionality.

### Configuration

Edit `ArChIPelago_code/config.yml` to set paths for your system:

```yaml
paths:
  genome_human: "/path/to/hg38.fa"       # UCSC hg38 reference genome
  genome_mouse: "/path/to/mm10.fa"       # UCSC mm10 reference genome
  repeatmasker: "/path/to/track_out.bed" # RepeatMasker BED (mouse only)
  macs_peaks:   "/path/to/macs/"         # GTRD MACS peak interval files
  output_dir:   "Release/TF-ML"         # Output directory

tools:
  sarus_jar:  "../sarus/releases/sarus-2.2.3.jar"
  java_bin:   "java"
```

---

## Data Download

### Zenodo archive (required)

**DOI: [10.5281/zenodo.14927304](https://zenodo.org/records/14927304)**

| Archive | Contents | Required for |
|---------|----------|-------------|
| `train.tar.gz` | Training FASTA sequences (36 TFs) | Notebooks 2, 5 |
| `test.tar.gz` | Test FASTA sequences (36 TFs) | Notebooks 2, 5 |
| `PWMs_mono_HUMAN.tar.gz` | Mononucleotide PWMs per TF | Notebooks 1, 2, 5 |
| `PWMs_di_HUMAN.tar.gz` | Dinucleotide PWMs per TF | Notebooks 1, 2, 5 |
| `hocomoco11.tar.gz` | HOCOMOCO v11 benchmark data | Notebook 1 |
| `Models.tar.gz` | Pre-trained RF models (sklearn 1.3) | Inference only |
| `macs.tar.gz` | GTRD MACS peak intervals | Notebook 0 |
| `Archipelago_intermediate_files.tar.gz` | Pre-computed feature matrices | Skip notebooks 0–1 |

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

## Full Reproduction — Notebooks 0–4

Run notebooks sequentially. Each reads `config.yml` for external paths.

| Step | Notebook | Input | Output |
|------|----------|-------|--------|
| 0 | Data preparation | GTRD BED + hg38/mm10 FASTA | Train/test FASTA and BED files |
| 1 | PWM scanning | FASTA + HOCOMOCO PWMs | Per-TF SPRY-SARUS feature matrices |
| 2 | Model training | Feature matrices | Trained RF/Slim models, auROC/auPRC scores |
| 3 | BiasAway QC | Negative sequences | GC-content quality control plots |
| 4 | Figures | Model scores | Publication figures |

**Data preparation (Notebook 0):** ChIP-Seq peaks from GTRD were called with MACS (Zhang et al. 2008). Peak lists were filtered by `tags >= 10` and sorted by `-log10(P value)`. Putative binding regions of 301 bp [-150;+150] were extracted centred at the peak summit. Repeat-overlapping peaks were removed using RepeatMasker (Smit et al. 2013-2015) via pybedtools `subtract(f=0.7, N=True)`.

**Train/test split:** Training uses chromosomes 2–7, 9–20; testing uses chromosomes 1, 8, 21 (human) and chromosomes 1, 8, 19 (mouse), following Zhou and Troyanskaya (2015). Up to 10,000 positives per TF; negatives (1:100 class balance) sampled from peaks of unrelated TF families (TFClass; Wingender et al. 2018), GC-matched using BiasAway (Khan et al. 2021). Sex chromosomes were excluded.

**PWM scanning:** Genomic regions were scanned with SPRY-SARUS (Kulakovskiy et al. 2016) using `--skipn --show-non-matching --output-scoring-mode score besthit`. Log-odds best-hit scores serve as features.

**Model parameters (Random Forest):** `max_depth=6, max_samples=0.8, n_estimators=100` (selected via GridSearchCV). Features were scale-transformed with `sklearn.preprocessing.StandardScaler`.

**Supported TFs (36):** ANDR, AP2A, CEBPB, COE1, CTCF, E2F4, ERG, ESR1, FLI1, GATA1, GATA2, GATA3, GCR, HNF4A, IRF1, IRF4, JUND, MAFK, MAX, MYC, P53, PPARG, PRGR, REST, RUNX1, RXRA, SOX2, SPI1, SRF, STA5A, STAT1, STAT3, TAL1, TF65, TFE2, USF2.

> **Compute requirements:** Full training across all 36 TFs was performed on 100 AMD EPYC 7662 cores (~8 hours). For a single TF on a laptop, use `n_jobs=4` and expect ~30 min.

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

## Publication Figures

R scripts for all publication figures are in `Figures/`. Generated PDFs and images are not tracked in git — run the scripts to regenerate them.

| Figure | Script | Description |
|--------|--------|-------------|
| Fig. 2 | `Figure 2 and S1/Figure_2_21012026_H_H.R` | Model comparison boxplots (human test set) |
| Fig. S1 | `Figure 2 and S1/Figure_2_210120267_H_M.R` | Model comparison boxplots (mouse validation) |
| Fig. 3 | `Figure 3 and S3/Figure_3_H_H_20012026.R` | PWM vs RF ensemble scatter + improvement (human) |
| Fig. S3 | `Figure 3 and S3/Figure_S3_H_M_20012026.R` | PWM vs RF ensemble scatter + improvement (mouse) |
| Fig. 4 | `Figure 4/Figure_4_H_H_7012026.R` | Slim/diChIPMunk comparison, delta (human) |
| Fig. S4 | `Figure 4/Figure_S4_H_M_7012026.R` | Slim/diChIPMunk comparison, delta (mouse) |
| Fig. 4 (abs) | `Figure 4/Figure_4_H_H_7012026_abs_values.R` | Absolute values variant |
| Fig. S4 (abs) | `Figure 4/Figure_S4_H_M_7012026_abs_values.R` | Absolute values variant |

The R scripts read from `HUMAN_MOUSE_total_100k.csv` (generated by Notebook 4). All mono+di comparisons use the mononucleotide single-PWM baseline, consistent with the analysis in the manuscript.

---

## Platform

Developed and tested on:
- **Ubuntu 20.04.6 LTS** (GNU/Linux 5.15.0-113-generic x86_64)
- **Python 3.8.18**, scikit-learn 1.3.2, numpy 1.24.4, pandas 2.0.3
- **R 4.3.1** with PRROC, ggplot2, patchwork
- **Java 8** (OpenJDK 8u232-b09, bundled)

---

## License

ArChIPelago is distributed under [WTFPL](http://www.wtfpl.net/). If you prefer a standard license, treat it as CC-BY.

---

## Citation

Kravchenko P., Vorontsov I.E., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding sites. *Bioinformatics* (ECCB 2026).

Zenodo data archive: [10.5281/zenodo.14927304](https://zenodo.org/records/14927304)
