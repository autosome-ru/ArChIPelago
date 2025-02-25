# ArChIPelago <img src='./Archipelago.png' width='55'>

This repository allows reproducing the _Kravchenko et al., biorxiv (2025)_ analysis. It contains scripts for sequences generation, model training/testing and plot generation for — ArChIPelago <https://github.com/autosome-ru/ArChIPelago> — the arrangement of multiple position weight matrices with ChIP-seq and machine learning for prediction of transcription factor binding sites.
</br>
Use the scanning tool from ArChIPelago-TFBS-finder: <https://github.com/Pavel-Kravchenko/ArChIPelago-TFBS-finder> to look for known motifs in the input sequences.
</br>
 
## Before you start

Make sure that you have installed:
<ul>
<li>Python 3.7 (or upper) https://www.python.org/</br></br>
<li>Optionally: R 4.2.1 (or upper) and RStudio https://posit.co/download/rstudio-desktop/</br></br>
</ul>

## Getting started

Please ```clone``` this directory with ```git clone https://github.com/autosome-ru/ArChIPelago/```</br></br>
Then ```cd``` in ArChIPelago: ```cd ArChIPelago```</br></br>
Install the SARUS PWM scanner into `./sarus` directory, copying the `sarus.jar`-file should be enough; see also the instructions at <https://github.com/autosome-ru/sarus></br></br>
_Note: SARUS is written in Java (hence it requires JRE)._</br></br>
Install ChIPMunk <https://autosome.org/ChIPMunk></br></br>
Download raw peaks data and PWMs from the Zenodo repository: 10.5281/zenodo.14927304. </br></br>

### Input data organization
(1) Download repeatmasker data and place it into ```Repeatmasker_data``` directory </br></br>
(2) Download and unzip the hg38 and mm10 reference genomes ```wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz; gunzip hg38.fa.gz``` ```wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz; gunzip mm10.fa.gz```. Place the files into the ```Genomic_data``` folder.</br></br>

### Steps to reproduce the ArChIPelago analysis

ArChIPelago was developed and tested on Ubuntu 20.04.6 LTS (GNU/Linux 5.15.0-113-generic x86_64). The model training process utilized 100 AMD EPYC 7662 64-Core Processor cores, with a total runtime of approximately 8 hours for training across all transcription factors (TFs).

**For DEMO tests use ["GABPA"]**

(1) Create the ArChIPelago environment by running ```conda env create -f environment_rpy2.yml``` and activate it with ```conda activate ArChIPelago```</br></br>

(2) Generate sequences and features from the train-test data splits:
-  ```0_Data_preparation and_test_train_split.ipynb``` - takes bed files from GTRD and extracts genomic sequences</br></br>
-  ```1_Scanning_with_CHIPMUNK_feature_generation_MONO_DI.ipynb``` - takes generated fasta files with genomic sequences for train/test with positive/negative examples and scanes them with PWMs using SARUS <https://github.com/autosome-ru/sarus></br></br>

(3) Train PWM-based models and classify sequences containing TFBSs:
-  ```2_ArChIPelago_and_Slim_training.ipynb``` - uses identified motif scores as features for different unified classifiers (i.e. RandomForestClassifier) models</br></br>
- To estimate the ArChIPelago performance with auROC and auPRC metrics, we complement the repository with scanning_tool.py from PRROC R package [Grau, J., Grosse, I. & Keilwagen, J. PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R. Bioinformatics 31, 2595–2597 (2015)]</br></br>

(4) Collect model performances and visualize results:
-  ```3_Biasaway_QC.ipynb``` - the script for Biasaway negative set QC</br></br>
-  ```4_Plot_generation_and_analysis.ipynb``` - the script for collection of model performances from 2_ArChIPelago_and_Slim_training.ipynb</br></br>

# Citing
TBA. Kravchenko et al., biorxiv (2025)
</br></br>

# License
ArChIPelago is distributed under WTFPL. If you prefer more standard licenses, feel free to treat WTFPL as CC-BY.

--2024--
