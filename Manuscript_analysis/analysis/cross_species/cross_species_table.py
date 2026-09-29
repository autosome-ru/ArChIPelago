#!/usr/bin/env python3
"""Per-TF diagnostic table for the human->mouse transfer of ArChIPelago models (Sup. Table 5, Fig. S6):
mouse test set = chromosomes 1, 8 and 19; mouse training set = chromosomes 2-7, 9, 10 and 13-18.

Inputs:
  * mouse performance (H>M: ArChIPelago RF mono+di seed 0 vs the best human monoPWM, train-selected by auROC /
    by auPRC) and the mouse training / test sizes (n_train_pos_M, n_test_pos_M) come from
    ../mouse_transfer/mouse_transfer_results.csv; the human values from the notebook 4 table (../inputs/notebook4_results_table.csv);
  * the motif similarity of the baseline human monoPWM to the nearest mouse monoPWM: motif_similarity.csv
    (motif_similarity.py; the baseline PWM name is asserted to be the one whose similarity was computed);
  * TF metadata: ../inputs/GTRD_experiments.xlsx, TF_table.xlsx, GTRD_metadata.csv.
Output: cross_species_table.csv in this folder.
Run with the system python3 (3.9); reads only; never overwrites an existing cross_species_table.csv.
"""
import os, re, sys, collections
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from archi_paths import INPUTS, BASIS  # noqa: E402
OUT = os.path.join(BASIS, "cross_species")
os.makedirs(OUT, exist_ok=True)
if os.path.exists(os.path.join(OUT, "cross_species_table.csv")):
    raise SystemExit("refusing to overwrite %s" % os.path.join(OUT, "cross_species_table.csv"))
SIM = pd.read_csv(os.path.join(OUT, "motif_similarity.csv")).set_index("TF").sort_index()

# 1. Results: human side from the notebook 4 table, mouse side from the mouse_transfer results
csv = pd.read_csv(os.path.join(INPUTS, "notebook4_results_table.csv"), sep="\t")
rf = (csv[(csv.Model == "RandomForestClassifier") & (csv.PWM == "mono+di")]
      .set_index("TF_name").sort_index())
assert len(rf) == 36 and not rf.index.duplicated().any()
TFS = list(rf.index)
mono = (csv[(csv.Model == "Single best mono PWM") & (csv.PWM == "mono")].drop_duplicates("TF_name")
        .set_index("TF_name").sort_index())              # the human best-monoPWM baselines
assert list(mono.index) == TFS
hm = pd.read_csv(os.path.join(BASIS, "mouse_transfer", "mouse_transfer_results.csv")).set_index("TF").sort_index()
assert list(hm.index) == TFS and hm.MM_available.all()

A = pd.DataFrame(index=TFS)
A["PWM_auROC_M"] = hm.HM_broc_test_M_auroc          # best human monoPWM (train auROC) on the mouse test set
A["PWM_auPRC_M"] = hm.HM_bprc_test_M_auprc          # best human monoPWM (train auPRC) on the mouse test set
A["ARCH_auROC_M"] = hm.HM_rf_test_M_auroc_s0
A["ARCH_auPRC_M"] = hm.HM_rf_test_M_auprc_s0
A["PWM_auROC_H"] = mono.roc_auc_test_H_PWM_mono
A["PWM_auPRC_H"] = mono.pr_auc_test_H_PWM_mono
A["ARCH_auROC_H"] = rf.roc_auc_test_H
A["ARCH_auPRC_H"] = rf.pr_auc_test_H
assert np.allclose(A.PWM_auROC_H, hm.HM_broc_test_H_auroc) and np.allclose(A.PWM_auPRC_H, hm.HM_bprc_test_H_auprc), \
    "human baselines of the mouse_transfer results must equal the notebook 4 best-monoPWM baselines"
A["dROC_M"] = A.ARCH_auROC_M - A.PWM_auROC_M
A["dPRC_M"] = A.ARCH_auPRC_M - A.PWM_auPRC_M
A["dROC_H"] = A.ARCH_auROC_H - A.PWM_auROC_H
A["dPRC_H"] = A.ARCH_auPRC_H - A.PWM_auPRC_H
A["n_PWM_csv"] = rf.Count
A["n_train_pos"] = rf.Seq_count
assert (A.n_train_pos == hm.n_train_pos_H).all()
A["n_train_pos_M"] = hm.n_train_pos_M
A["n_test_pos_M"] = hm.n_test_pos_M
below_M = sorted(A.index[(A.dROC_M < 0) | (A.dPRC_M < 0)])

# 2. Dataset and PWM counts (pipeline TF table, Sup. Table 1 experiments) and cell types from metadata
t2 = pd.read_excel(os.path.join(INPUTS, "TF_table.xlsx")).dropna(subset=["TF_name"]).set_index("TF_name").sort_index()
assert list(t2.index) == TFS
t1h = pd.read_excel(os.path.join(INPUTS, "GTRD_experiments.xlsx"), sheet_name="Sheet1").assign(sp="H")
t1m = pd.read_excel(os.path.join(INPUTS, "GTRD_experiments.xlsx"), sheet_name="Sheet2").assign(sp="M")
t1 = pd.concat([t1h, t1m], ignore_index=True)
assert (t1.groupby(["TF", "sp"]).size().unstack()["H"] == t2.HUMAN_datasets).all()
assert (t1.groupby(["TF", "sp"]).size().unstack()["M"] == t2.MOUSE_datasets).all()

# The metadata csv has unquoted commas inside antibody/treatment fields, which shifts
# columns for ~100 experiments. The cell-line field is always the one immediately
# preceding the PEAKSxxxxxx token, so parse the raw lines instead of trusting pandas.
raw = {}
for line in open(os.path.join(INPUTS, "GTRD_metadata.csv")):
    tok = line.rstrip("\n").split(",")
    if not tok[0].startswith("EXP"):
        continue
    pi = [i for i, t in enumerate(tok) if re.fullmatch(r"PEAKS\d+", t)]
    if not pi:
        continue
    i = pi[0]
    raw.setdefault(tok[0], dict(cell=tok[i - 1].strip(), tfclass=tok[2],
                                 uniprot=next((t for t in tok if t.endswith("_HUMAN") or t.endswith("_MOUSE")), "")))
t1["cell_raw"] = t1.exp.map(lambda e: raw[e]["cell"])
t1["TFclass"] = t1.exp.map(lambda e: raw[e]["tfclass"])
assert t1.cell_raw.notna().all()
meta_n = collections.Counter()
for e, d in raw.items():
    pref = d["uniprot"].split("_")[0]
    if pref in TFS:
        meta_n[(pref, "H" if d["uniprot"].endswith("_HUMAN") else "M")] += 1

# manual abbreviation of raw cell-line strings -> short label (tissue/lineage in brackets where useful)
ABBR = {
 # prostate / AR
 "VCaP": "VCaP", "VCaP cell line": "VCaP", "LNCaP": "LNCaP", "LNCAP": "LNCaP", "LNCaP cell line": "LNCaP",
 "LNCaP androgen-dependent prostate cancer cell line": "LNCaP", "LNCaP-1F5": "LNCaP-1F5", "LNCaP-1F5 cells": "LNCaP-1F5",
 "LNCaP-ARmo": "LNCaP-AR", "LNCaP-ARhi": "LNCaP-AR", "LNCaP-pcDNA3.1": "LNCaP", "C4-2B": "C4-2B",
 "prostate cancer cell line DU145": "DU145", "prostate cancer": "prostate ca.", "RWPE-1": "RWPE-1", "RWPE1": "RWPE-1",
 "caput epididymis": "epididymis", "R26-ERG/ERG": "prostate(ERG-tg)", "Pten-f/f": "prostate(Pten-/-)",
 "Mouse Prostate": "prostate", "prostate": "prostate", "ventral prostate": "prostate", "kidney": "kidney",
 # common human lines
 "HeLa-S3": "HeLa-S3", "HeLa S3": "HeLa-S3", "HelaS3": "HeLa-S3", "HeLa": "HeLa", "LoVo": "LoVo", "ECC-1": "ECC-1",
 "A549": "A549", "HCT-116": "HCT-116", "HepG2": "HepG2", "Blood monocytes": "monocytes", "GM12878": "GM12878",
 "Adipose stromal cell (ASC) adipocyte": "ASC-adipocyte", "hASC": "hASC", "HUVEC": "HUVEC", "primary HUVEC cells": "HUVEC",
 "umbilical vein endothelial cells": "HUVEC", "CD36": "CD36+ erythroid", "Pancreatic islets": "panc. islets", "HMEC": "HMEC",
 "GP5d": "GP5d", "MCF-7": "MCF-7", "MCF7": "MCF-7", "MCF-7 breast cancer cells": "MCF-7", "FB0167P": "fibroblast",
 "HFF": "HFF", "GM20000": "LCL", "Monocytes-CD14+_RO01746": "CD14+ monocytes", "GM10266": "LCL", "NHLF": "NHLF",
 "CD20+": "CD20+ B", "Pancreas_OC": "pancreas", "NHEK": "NHEK", "MM1.S": "MM1.S", "passage39": "fibroblast",
 "passage 95": "fibroblast", "fetal lung fibroblast": "fetal lung fibr.", "GM13977": "LCL", "HSMMtube": "HSMMtube",
 "Lung_OC": "lung", "WI-38": "WI-38", "BL41": "BL41", "CD34": "CD34+", "BJAB": "BJAB", "NB4": "NB4", "HMF": "HMF",
 "HCFaa": "HCFaa", "AG09309": "AG09309", "K562": "K562", "Osteobl": "osteoblast", "T-47D": "T-47D", "T47D": "T-47D",
 "T47D breast cancer cells": "T-47D", "HVMF": "HVMF", "SK-N-SH_RA": "SK-N-SH", "SK-N-SH": "SK-N-SH", "NH-A": "NH-A",
 "GM12801": "LCL", "HL-60": "HL-60", "THP-1": "THP-1", "GM06990": "LCL", "Mesenchymal stem cells-derived adipocytes": "MSC-adipocyte",
 "CD34+ cells": "CD34+", "primary breast tumor": "breast tumor", "MDA-MB-134": "MDA-MB-134", "Ishikawa breast cancer cells": "Ishikawa",
 "ECC1 breast cancer cells": "ECC-1", "H3396 cells": "H3396", "ZR751": "ZR-75-1",
 "megakaryocytes cultured from cord blood CD34-positive cells": "CB megakaryocytes", "Pro-erythroblasts": "proerythroblasts",
 "primary human PBMC derived erythroblasts": "erythroblasts", "PBDE": "erythroblasts", "CD34+ HSPC-derived proerythroblasts": "proerythroblasts",
 "Jurkat E6-1": "Jurkat", "Jurkat": "Jurkat", "Jurkat Trex (Invitrogen)": "Jurkat", "clone E6-1 T lymphocyte": "Jurkat",
 "Th1": "Th1", "Th2": "Th2", "CCRF-CEM": "CCRF-CEM", "Human Embryonal Kidney Cells": "HEK293", "Caco-2": "Caco-2",
 "CD14+ monocytes": "CD14+ monocytes", "OCI-LY3": "OCI-LY3", "OCI-LY10": "OCI-LY10", "H929": "H929", "P493-6": "P493-6",
 "H128": "H128", "H1-hESC": "H1-hESC", "Plasma cell": "plasma cell", "H2171": "H2171", "FB8470": "FB8470", "Ramos": "Ramos",
 "CA46": "CA46", "Blue1": "Blue1", "U2OS": "U2OS", "osteosarcoma U2OS": "U2OS", "Human neonatal foreskin keratinocytes": "keratinocytes",
 "H9 hESC": "H9-hESC", "H9": "H9-hESC", "IMR90": "IMR90", "lymphoblastoid cells": "LCL",
 "Human diploid foreskin fibroblasts (BJ) immortalized with hTERT": "BJ-hTERT",
 "Simpson Golabi Behmel Syndrome (SGBS) human adi- pocytes": "SGBS adipocytes", "adipocyte": "adipocyte", "AB32": "AB32",
 "PANC-1": "PANC-1", "PFSK-1": "PFSK-1", "U87": "U87", "H295R/TR SF-1 cells": "H295R",
 "K562 erythrocytic leukaemia cells (ATCC CCL-243)": "K562", "APL bone marrow blasts": "APL blasts",
 "AML-patient-derived cell lines t(8;21) Kasumi-1": "Kasumi-1", "CMK": "CMK", "LS180 cells": "LS180", "ReN-VM": "ReN-VM(NPC)",
 "KYSE70": "KYSE70", "TT": "TT", "HCC95": "HCC95", "monocyte derived macrophages": "macrophages",
 "CD133+ expanded umbilical cord blood (UCB) cells": "CB CD133+", "CD4+CD25+CD45RA+ expanded naive regulatory T cells": "Treg",
 "CD4+CD25- expanded conventional T cells": "Tconv", "MCF10A-Er-Src": "MCF10A-Src", "U-2932": "U-2932",
 "Leukemic T-cell": "T-ALL", "GM18526": "LCL", "GM12891": "LCL", "GM19193": "LCL", "GM18505": "LCL", "GM18951": "LCL",
 "GM19099": "LCL", "GM15510": "LCL", "GM12892": "LCL", "GM10847": "LCL", "RPMI-8402": "RPMI-8402(T-ALL)",
 # mouse
 "3T3-L1 preadipocyte cell line": "3T3-L1", "3T3-L1": "3T3-L1", "3T3L1": "3T3-L1", "3T3-L1 6hr of differentiation": "3T3-L1",
 "3T3-L1 adipocytes": "3T3-L1 adipocytes", "C2C12": "C2C12", "C/EBPbeta ChIP": "n/a", "C57BL6/J": "tissue(C57BL6)",
 "differentiating murine hematopoietic cells": "diff. hematopoietic", "liver": "liver", "Liver": "liver", "Mouse liver": "liver",
 "MC3T3-E1": "MC3T3-E1", "MC3T3-E1 differentiated for 15 days": "MC3T3-E1", "A-MuLV pre B-cells": "pre-B",
 "splenic B cells": "splenic B", "Bone Marrow": "bone marrow", "bone marrow": "bone marrow", "Mouse bone marrow": "bone marrow",
 "BoneMarrow": "bone marrow", "bone marrow cells": "bone marrow", "Pro-B cell": "pro-B", "bone marrow-derived dendritic cells": "BMDC",
 "Embryonic stem cells": "ESC", "Embryonic Stem Cells (ESCs)": "ESC", "embryonic stem cells (ESCs)": "ESC", "ESC": "ESC", "mESC": "ESC",
 "embryonic stem (ES) cells": "ESC", "Bruce4 embryonic stem cells": "ESC", "V6.5 embryonic stem cells": "ESC", "R1/E ESC": "ESC",
 "Differentiated mouse ES cells": "diff. ESC", "Undifferentiated mouse ES cells": "ESC", "Mouse placenta": "placenta",
 "Mouse heart": "heart", "Heart": "heart", "Mouse E14.5 heart": "E14.5 heart", "Mouse lung": "lung", "Lung": "lung",
 "Mouse spleen": "spleen", "Spleen": "spleen", "Mouse E14.5 brain": "E14.5 brain", "Mouse E14.5 limb": "E14.5 limb", "Limb": "limb",
 "Mouse intestine": "intestine", "SmIntestine": "sm. intestine", "Mouse testis": "testis", "Testis": "testis",
 "Mouse embryonic fibroblasts": "MEF", "mouse embryonic fibroblasts (MEFs)": "MEF", "mouse embryo fibroblast (MEF)": "MEF",
 "NIH-3T3": "NIH-3T3", "NIH3T3 fibroblasts": "NIH-3T3", "Tcf3(E2A) deficient pre-pro-B": "pre-pro-B", "Rag1 deficient pro-B": "pro-B",
 "Rag1-/- with TCR-beat transgene DP thymocytes of C57BL/6 mice": "DP thymocytes", "Rag1 -/- pro-B cells of C57BL/6": "pro-B",
 "splenic Plasmablasts": "plasmablasts", "pro-B cells of C57BL/6 mice": "pro-B", "Mouse cerebellum": "cerebellum", "Cerebellum": "cerebellum",
 "Mouse cortex": "cortex", "Cortex": "cortex", "WholeBrain": "brain", "Kidney": "kidney", "G1E-ER4": "G1E-ER4", "G1E ER4": "G1E-ER4",
 "asynchrnous G1E-ER4+E2 cells": "G1E-ER4", "G1E cells": "G1E", "G1E": "G1E", "Erthyroid Progenitor (G1E)": "G1E", "MEL": "MEL",
 "CH12": "CH12", "CD43- B cells extracted from spleen": "splenic B", "Rag1 -/- CD4-CD8- (DN) thymocytes": "DN thymocytes",
 "Sox2 and Oct4 (KSO)": "MEF(OSK reprogr.)", "spleen cells": "spleen", "GATA1s": "GATA1s megakaryoblasts",
 "uterus": "uterus", "Uterus": "uterus", "uterus epithelium": "uterus",
 "mouse mammary cells (express GFP-GR and Cherry-ER under a tet regulated promoter)": "mammary 3134",
 "and Afos under a tet regulated promoter)": "mammary 3134", "7438 cells": "mammary 7438",
 "Early haematopoietic cell line derived from ES cells": "HPC-7", "Megakaryo": "megakaryocytes", "Ter119+ cells": "Ter119+ erythroid",
 "Erythrobl": "erythroblasts", "Megakaryocyte Progenitor Cell Line (G1ME)": "G1ME", "G1ME": "G1ME", "ES-EP": "ES-EP", "FDCPmix": "FDCP-mix",
 "Lin- bone marrow hematopoietic progenitor cells": "Lin- HPC", "primary CD4 T cells": "CD4 T", "Fetal liver precursor derived DN1": "DN1 thymocytes",
 "Fetal liver precursor derived DN2b": "DN2b thymocytes", "Adult thymic DP (pre-positive selection)": "DP thymocytes",
 "Mouse villus": "intestinal villus", "Cdx2 KO": "intestinal villus(Cdx2-/-)", "B cells": "B cells", "Th17 Cells": "Th17",
 "STHdhQ7": "STHdh striatal", "STHdhQ111": "STHdh striatal", "p53R172H/+-mice": "tissue(p53 mut.)",
 "inguinal white adipous tissue-derived adipocytes": "WAT adipocytes", "brown adipous tissue-derived adipocytes": "BAT adipocytes",
 "Lin-SCA+CD34+ cells": "LSK HPC", "L8057": "L8057(megakaryoblast)", "megakaryocyte cells": "megakaryocytes", "BMiFLT3(15-3)": "BM iFLT3",
 "F9 embryonal carcinoma cells": "F9 EC", "RMWT_227": "RMWT_227", "IDG-SW3": "IDG-SW3(osteocyte)",
 "Day2 differentiating Teto-Brn2 ESCs": "diff. ESC", "Day2 differentiating rtta ESCs": "diff. ESC",
 "ES-derived neural progenitor cells (NPCs)": "ES-NPC", "Neuralized embryoid bodies": "neural EB", "ZHBTc4-TS cells": "TS cells",
 "EGFP-TS3.5 TS cells": "TS cells", "Macrophage": "macrophages", "resting mature peripheral primary B cells": "B cells",
 "38B9 cells": "38B9(pre-B)", "Murine erythroleukemia (MEL) cells stably transfected with Gata-1 fused to ER": "MEL",
 "control RAW macrophages": "RAW264.7", "biotinylated RevErba expressing RAW macrophages": "RAW264.7",
 "HL-1 cardiac muscle cell line": "HL-1", "total T cells": "T cells", "mammary": "mammary", "Mammary tissues": "mammary",
 "mammary gland": "mammary", "FH": "n/a(FH)", "MH": "n/a(MH)", "Bone marrow macrophages": "BMDM", "Bone marrow macrophages ikk-/-": "BMDM",
 "AtT-20": "AtT-20", "CD4+ T Cells": "CD4 T", "fetal liver erythroblast": "FL erythroblasts", "Rag2-/- thymocytes": "thymocytes",
}
t1["cell"] = t1.cell_raw.map(ABBR)
missing = sorted(set(t1.cell_raw[t1.cell.isna()]))
if missing:
    sys.exit(f"Unmapped cell-line strings: {missing}")

def celllist(tf, sp):
    v = t1[(t1.TF == tf) & (t1.sp == sp)].cell.value_counts()
    return "; ".join(f"{k}(x{n})" if n > 1 else k for k, n in v.items())

D = pd.DataFrame(index=TFS)
D["n_exp_H"] = t2.HUMAN_datasets.astype(int)
D["n_exp_M"] = t2.MOUSE_datasets.astype(int)
D["n_exp_H_metafile"] = [meta_n[(t, "H")] for t in TFS]
D["n_exp_M_metafile"] = [meta_n[(t, "M")] for t in TFS]
D["n_monoPWM_H"] = t2.HUMAN_monoPWMs.astype(int)
D["n_diPWM_H"] = t2.HUMAN_diPWMs.astype(int)
D["n_monoPWM_M"] = t2.MOUSE_monoPWMs.astype(int)
D["train_pos_T2"] = t2.HUMAN_train_set.astype(int)
D["cells_H"] = [celllist(t, "H") for t in TFS]
D["cells_M"] = [celllist(t, "M") for t in TFS]
assert (D.train_pos_T2 == A.n_train_pos).all()

# manual judgement of shared / equivalent tissue between the human and mouse experiments
OVERLAP = {
 "ANDR": ("yes", "prostate (LNCaP/VCaP vs mouse prostate)"),
 "AP2A": ("no", "HeLa/LoVo vs epididymis/prostate"),
 "CEBPB": ("partial", "HepG2 vs liver; monocytes vs hematopoietic; rest unmatched (adipocyte/myoblast)"),
 "COE1": ("yes", "GM12878 (B-LCL) vs pre-B/B cells; ASC vs 3T3-L1"),
 "CTCF": ("yes", "many shared lineages (ESC, erythroid, B, fibroblast, liver, lung)"),
 "E2F4": ("partial?", "LCL vs CH12 (B lineage); others unmatched (HeLa/K562 vs MEF/C2C12)"),
 "ERG": ("yes", "prostate (VCaP/RWPE-1 vs ERG-tg prostate); CD34+ vs hematopoietic"),
 "ESR1": ("yes", "breast/endometrium vs mammary/uterus"),
 "FLI1": ("yes", "megakaryocytes in both"),
 "GATA1": ("yes", "erythroid/megakaryocytic in both"),
 "GATA2": ("yes", "hematopoietic progenitors/K562 vs FDCP-mix/G1E"),
 "GATA3": ("yes", "T cells in both (Th1/Th2/Jurkat vs CD4 T/thymocytes)"),
 "GCR": ("no", "A549/LNCaP/ECC-1 vs mammary 3134/3T3-L1"),
 "HNF4A": ("yes", "colon lines (Caco-2/LoVo/GP5d) vs intestinal villus"),
 "IRF1": ("partial", "CD14+ monocytes vs BMDC (both LPS-stimulated myeloid)"),
 "IRF4": ("yes", "B-cell lymphoma/LCL/myeloma vs B/pro-B cells"),
 "JUND": ("partial?", "GM12878 vs CH12 (B); SK-N-SH vs striatal; majority unmatched"),
 "MAFK": ("yes", "K562 vs MEL (erythroleukemia); LoVo/CH12 unmatched"),
 "MAX": ("yes", "K562 vs MEL; H1-hESC vs ESC"),
 "MYC": ("yes", "K562 vs MEL; LCL vs CH12/pro-B; ESC"),
 "P53": ("yes", "H9-hESC vs ESC; BJ/IMR90 fibroblasts vs MEF"),
 "PPARG": ("yes", "adipocytes in both"),
 "PRGR": ("no?", "T-47D/AB32 vs uterus (both PR-responsive reproductive tissue, not the same tissue)"),
 "REST": ("yes", "H1-hESC vs ESC"),
 "RUNX1": ("yes", "hematopoietic/megakaryocytic leukaemia lines vs hematopoietic progenitors/megakaryocytes"),
 "RXRA": ("yes", "HepG2 vs liver"),
 "SOX2": ("yes", "H9-hESC/ReN-VM vs ESC/NPC"),
 "SPI1": ("yes", "macrophages/THP-1/B-lymphoma vs macrophages/B cells"),
 "SRF": ("no", "LCL/epithelial lines vs cardiomyocyte/fibroblast/myoblast"),
 "STA5A": ("yes", "T cells in both"),
 "STAT1": ("no", "K562/HeLa vs BMDM"),
 "STAT3": ("no", "HeLa/MCF10A/U-2932 vs AtT-20/ESC/CD4 T"),
 "TAL1": ("yes", "erythroid in both"),
 "TF65": ("no", "LCL/fibroblast/HUVEC vs macrophages/BMDC"),
 "TFE2": ("yes", "RPMI-8402 (T-ALL) vs thymocytes; GM12878 unmatched"),
 "USF2": ("yes", "K562 vs MEL; GM12878 vs CH12"),
}
D["celltype_overlap"] = [OVERLAP[t][0] for t in TFS]
D["celltype_overlap_note"] = [OVERLAP[t][1] for t in TFS]

# 3. Family labelling
tfclass = t1.groupby("TF").TFclass.first()
assert (t1.groupby("TF").TFclass.nunique() == 1).all()
# family label; n_paralogs = approx. number of human genes in the TFClass subfamily/genus
# that bind an essentially identical core motif (knowledge-based).
FAM = {
 "ANDR": ("NR3C steroid receptor", 4), "GCR": ("NR3C steroid receptor", 4), "PRGR": ("NR3C steroid receptor", 4),
 "ESR1": ("NR3A estrogen receptor", 2), "HNF4A": ("NR2A HNF4", 2), "RXRA": ("NR2B RXR", 3), "PPARG": ("NR1C PPAR", 3),
 "AP2A": ("AP-2", 5), "CEBPB": ("bZIP C/EBP", 6), "JUND": ("bZIP AP-1 (Jun)", 3), "MAFK": ("bZIP small MAF", 3),
 "COE1": ("EBF", 4), "CTCF": ("C2H2 ZF CTCF", 2), "REST": ("C2H2 ZF REST", 1), "E2F4": ("E2F", 8),
 "ERG": ("ETS", 27), "FLI1": ("ETS", 27), "SPI1": ("ETS (SPI)", 3), "GATA1": ("GATA", 6), "GATA2": ("GATA", 6), "GATA3": ("GATA", 6),
 "IRF1": ("IRF", 9), "IRF4": ("IRF", 9), "MAX": ("bHLH-ZIP MYC/MAX", 10), "MYC": ("bHLH-ZIP MYC/MAX", 10), "USF2": ("bHLH-ZIP USF", 3),
 "TFE2": ("bHLH E2A", 3), "TAL1": ("bHLH TAL", 3), "P53": ("p53", 3), "RUNX1": ("RUNX", 3), "SOX2": ("SOX (B1)", 3),
 "SRF": ("MADS SRF", 1), "STA5A": ("STAT", 7), "STAT1": ("STAT", 7), "STAT3": ("STAT", 7), "TF65": ("NF-kB/Rel", 5),
}
D["TFclass_id"] = tfclass.reindex(TFS)
D["family"] = [FAM[t][0] for t in TFS]
D["n_paralogs_approx"] = [FAM[t][1] for t in TFS]

# 4. Motif similarity human vs mouse monoPWMs (motif_similarity.py);
#    the baseline monoPWM must be the PWM whose similarity was computed
names = pd.read_csv(os.path.join(INPUTS, "single_PWM_features.csv"))
names = names.set_index(["TF", "feature"]).pwm_name
E = pd.DataFrame(index=TFS)
for tf in TFS:
    base = names.loc[(tf, hm.loc[tf, "HM_broc_col"])]
    assert base == SIM.loc[tf, "topH_PWM"], (tf, base, SIM.loc[tf, "topH_PWM"])
for c in ["n_monoPWM_H_found", "n_monoPWM_M_found", "topH_PWM", "topH_vs_nearestM_pcc",
          "median_H_vs_nearestM_pcc", "median_H_vs_nearestH_pcc"]:
    E[c] = SIM[c]
assert (E.n_monoPWM_M_found == D.n_monoPWM_M).all() and (E.n_monoPWM_H_found == hm.P_mono_H).all()

# 5. Assemble table
T = pd.concat([A, D, E], axis=1)
T.index.name = "TF"
T["below_baseline_M"] = T.index.isin(below_M)
order = ["below_baseline_M", "family", "TFclass_id", "n_paralogs_approx", "n_monoPWM_H", "n_diPWM_H", "n_monoPWM_M", "n_train_pos",
         "n_train_pos_M", "n_test_pos_M",
         "n_exp_H", "n_exp_M", "n_exp_H_metafile", "n_exp_M_metafile", "cells_H", "cells_M", "celltype_overlap", "celltype_overlap_note",
         "PWM_auROC_M", "PWM_auPRC_M", "ARCH_auROC_M", "ARCH_auPRC_M", "dROC_M", "dPRC_M", "dROC_H", "dPRC_H",
         "PWM_auROC_H", "PWM_auPRC_H", "ARCH_auROC_H", "ARCH_auPRC_H",
         "n_monoPWM_H_found", "n_monoPWM_M_found", "topH_PWM", "topH_vs_nearestM_pcc", "median_H_vs_nearestM_pcc", "median_H_vs_nearestH_pcc"]
T = T[order]
T.to_csv(os.path.join(OUT, "cross_species_table.csv"), float_format="%.4f")

# 6. Statistics
lines = ["| Spearman (n=36) | rho | p |", "|---|---|---|"]
def sp(x, y, label):
    r, p = stats.spearmanr(x, y)
    lines.append(f"| {label} | {r:+.3f} | {p:.3g} |")
    return r, p

for dcol, dl in [("dROC_M", "dauROC_M"), ("dPRC_M", "dauPRC_M")]:
    sp(T[dcol], T.n_train_pos, f"{dl} vs # human training positives")
    sp(T[dcol], T.n_train_pos_M, f"{dl} vs # mouse training positives")
    sp(T[dcol], T.n_exp_M, f"{dl} vs # mouse experiments (Sup. Table 1)")
    sp(T[dcol], T.n_exp_H, f"{dl} vs # human experiments (Sup. Table 1)")
    sp(T[dcol], T.n_monoPWM_H, f"{dl} vs # human monoPWMs")
    sp(T[dcol], T.PWM_auROC_M if dcol == "dROC_M" else T.PWM_auPRC_M, f"{dl} vs best-PWM {'auROC' if dcol=='dROC_M' else 'auPRC'} on mouse (ceiling)")
    sp(T[dcol], T.dROC_H if dcol == "dROC_M" else T.dPRC_H, f"{dl} vs {'dauROC_H' if dcol=='dROC_M' else 'dauPRC_H'} (within-human gain)")
    sp(T[dcol], T.topH_vs_nearestM_pcc, f"{dl} vs motif sim. baseline human PWM vs nearest mouse PWM")
    sp(T[dcol], T.median_H_vs_nearestM_pcc, f"{dl} vs motif sim. median(human PWM -> nearest mouse PWM)")
    sp(T[dcol], T.n_paralogs_approx, f"{dl} vs # paralogs (approx.)")
stats_md = "\n".join(lines)

fb = T.below_baseline_M
mw = []
if fb.any():
    for col, lab in [("n_train_pos", "# human training positives"), ("n_train_pos_M", "# mouse training positives"),
                     ("n_exp_H", "# human experiments"), ("n_exp_M", "# mouse experiments"),
                     ("PWM_auROC_M", "best-PWM auROC on mouse"), ("PWM_auPRC_M", "best-PWM auPRC on mouse"), ("dROC_H", "dauROC_H"),
                     ("topH_vs_nearestM_pcc", "baseline-human vs nearest-mouse PWM PCC"), ("median_H_vs_nearestM_pcc", "median human->mouse nearest PCC"),
                     ("n_paralogs_approx", "# paralogs (approx.)")]:
        u, p = stats.mannwhitneyu(T.loc[fb, col], T.loc[~fb, col], alternative="two-sided")
        mw.append(f"| {lab} | {T.loc[fb, col].median():.3g} | {T.loc[~fb, col].median():.3g} | {p:.3f} |")
mw_md = (f"| Mann-Whitney ({int(fb.sum())} below vs {int((~fb).sum())} others) | median (below) | median (others) | p |\n|---|---|---|---|\n"
         + "\n".join(mw)) if mw else "no TF below the baseline on the mouse test set"

print("below the baseline on the mouse test set (chr1/8/19):", below_M)
print(stats_md); print(mw_md)
print(T[["dROC_M", "dPRC_M", "n_train_pos_M", "topH_vs_nearestM_pcc", "median_H_vs_nearestM_pcc"]].round(3).to_string())
print("written:", os.path.join(OUT, "cross_species_table.csv"))
