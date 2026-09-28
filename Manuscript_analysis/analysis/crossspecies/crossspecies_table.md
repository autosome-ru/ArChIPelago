# Cross-species (human -> mouse) per-TF diagnostic table (mouse test set = chr1/8/19)

Performance columns: human values from the notebook 4 results table (RandomForestClassifier, PWM=mono+di rows); mouse values from hm_mm_results.csv (RF seed 0; best human monoPWM selected on the human training set, by auROC for the auROC columns and by auPRC for the auPRC columns). `below` marks the TFs below the best monoPWM on the mouse test set on >= 1 metric.

## Main table

| TF | below | family | #PWM H (mono+di) | #train pos H/M | #exp H/M | cell-type overlap | PWM auROC M | PWM auPRC M | ArCh auROC M | ArCh auPRC M | dauROC M | dauPRC M | dauROC H | dauPRC H | motif sim H->M |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| STA5A | **below** | STAT | 9+5 | 2550/10000 | 4/11 | yes | 0.876 | 0.385 | 0.870 | 0.401 | -0.006 | +0.016 | +0.028 | +0.036 | 0.952 |
| ERG | **below** | ETS | 37+17 | 6048/7106 | 6/4 | yes | 0.800 | 0.073 | 0.814 | 0.065 | +0.014 | -0.008 | +0.009 | +0.012 | 0.946 |
| P53 | **below** | p53 | 48+26 | 662/2144 | 11/4 | yes | 0.807 | 0.539 | 0.866 | 0.519 | +0.060 | -0.020 | +0.044 | +0.255 | 0.858 |
| TF65 |  | NF-kB/Rel | 72+36 | 2008/10000 | 22/10 | no | 0.803 | 0.114 | 0.803 | 0.115 | +0.000 | +0.001 | +0.037 | +0.056 | 0.962 |
| CTCF |  | C2H2 ZF CTCF | 198+99 | 10000/10000 | 57/51 | yes | 0.978 | 0.746 | 0.980 | 0.749 | +0.002 | +0.002 | +0.008 | +0.027 | 0.991 |
| USF2 |  | bHLH-ZIP USF | 12+6 | 6522/4260 | 5/3 | yes | 0.979 | 0.341 | 0.981 | 0.353 | +0.002 | +0.012 | +0.004 | +0.047 | 0.970 |
| COE1 |  | EBF | 10+5 | 3924/10000 | 3/5 | yes | 0.952 | 0.440 | 0.956 | 0.449 | +0.004 | +0.009 | +0.003 | +0.025 | 0.977 |
| CEBPB |  | bZIP C/EBP | 20+10 | 10000/10000 | 9/17 | partial | 0.922 | 0.450 | 0.929 | 0.495 | +0.007 | +0.045 | +0.010 | +0.092 | 0.990 |
| MAX |  | bHLH-ZIP MYC/MAX | 42+21 | 10000/5542 | 19/4 | yes | 0.861 | 0.147 | 0.869 | 0.162 | +0.008 | +0.015 | +0.025 | +0.088 | 0.936 |
| FLI1 |  | ETS | 18+8 | 1346/10000 | 1/3 | yes | 0.822 | 0.138 | 0.831 | 0.147 | +0.009 | +0.009 | +0.014 | +0.008 | 0.989 |
| PPARG |  | NR1C PPAR | 17+11 | 5435/10000 | 5/10 | yes | 0.868 | 0.421 | 0.878 | 0.434 | +0.010 | +0.012 | +0.014 | +0.021 | 0.983 |
| JUND |  | bZIP AP-1 (Jun) | 33+17 | 10000/8648 | 10/5 | partial? | 0.848 | 0.117 | 0.860 | 0.138 | +0.012 | +0.021 | +0.015 | +0.039 | 0.970 |
| HNF4A |  | NR2A HNF4 | 26+13 | 10000/10000 | 8/4 | yes | 0.927 | 0.313 | 0.940 | 0.381 | +0.012 | +0.069 | +0.012 | +0.055 | 0.984 |
| STAT3 |  | STAT | 5+5 | 3000/10000 | 3/4 | no | 0.794 | 0.062 | 0.807 | 0.073 | +0.013 | +0.011 | +0.015 | +0.001 | 0.914 |
| GATA1 |  | GATA | 28+14 | 10000/10000 | 13/18 | yes | 0.898 | 0.343 | 0.911 | 0.368 | +0.013 | +0.026 | +0.014 | +0.019 | 0.991 |
| TFE2 |  | bHLH E2A | 5+3 | 1939/10000 | 2/3 | yes | 0.882 | 0.185 | 0.896 | 0.242 | +0.014 | +0.057 | +0.045 | +0.026 | 0.970 |
| MYC |  | bHLH-ZIP MYC/MAX | 110+55 | 10000/10000 | 28/6 | yes | 0.875 | 0.167 | 0.888 | 0.214 | +0.014 | +0.047 | +0.022 | +0.090 | 0.949 |
| RUNX1 |  | RUNX | 39+20 | 10000/10000 | 9/8 | yes | 0.794 | 0.096 | 0.808 | 0.117 | +0.014 | +0.020 | +0.023 | +0.035 | 0.982 |
| IRF1 |  | IRF | 11+5 | 3280/5055 | 5/4 | partial | 0.789 | 0.223 | 0.803 | 0.241 | +0.015 | +0.018 | +0.035 | +0.062 | 0.973 |
| REST |  | C2H2 ZF REST | 70+36 | 10000/10000 | 33/5 | yes | 0.956 | 0.742 | 0.971 | 0.819 | +0.015 | +0.077 | +0.016 | +0.051 | 0.975 |
| PRGR |  | NR3C steroid receptor | 19+9 | 1260/3316 | 2/3 | no? | 0.862 | 0.113 | 0.881 | 0.119 | +0.019 | +0.007 | +0.036 | +0.009 | 0.920 |
| AP2A |  | AP-2 | 6+3 | 2406/5646 | 2/3 | no | 0.863 | 0.095 | 0.883 | 0.109 | +0.020 | +0.014 | +0.032 | +0.027 | 0.967 |
| SPI1 |  | ETS (SPI) | 42+21 | 10000/10000 | 11/18 | yes | 0.938 | 0.562 | 0.961 | 0.642 | +0.023 | +0.080 | +0.024 | +0.120 | 0.984 |
| ESR1 |  | NR3A estrogen receptor | 190+100 | 10000/10000 | 36/9 | yes | 0.869 | 0.371 | 0.894 | 0.388 | +0.024 | +0.017 | +0.020 | +0.042 | 0.956 |
| GCR |  | NR3C steroid receptor | 47+26 | 7758/10000 | 14/5 | no | 0.841 | 0.420 | 0.865 | 0.435 | +0.024 | +0.014 | +0.018 | +0.037 | 0.943 |
| TAL1 |  | bHLH TAL | 16+14 | 7020/10000 | 10/9 | yes | 0.892 | 0.161 | 0.918 | 0.195 | +0.026 | +0.034 | +0.047 | +0.032 | 0.966 |
| ANDR |  | NR3C steroid receptor | 196+108 | 10000/10000 | 31/13 | yes | 0.920 | 0.546 | 0.948 | 0.579 | +0.029 | +0.033 | +0.060 | +0.179 | 0.926 |
| MAFK |  | bZIP small MAF | 6+3 | 2458/4190 | 2/3 | yes | 0.909 | 0.533 | 0.941 | 0.658 | +0.032 | +0.125 | +0.048 | +0.180 | 0.946 |
| GATA3 |  | GATA | 44+19 | 8900/10000 | 9/4 | yes | 0.823 | 0.167 | 0.856 | 0.260 | +0.033 | +0.093 | +0.039 | +0.040 | 0.981 |
| SOX2 |  | SOX (B1) | 30+16 | 6374/10000 | 5/12 | yes | 0.808 | 0.130 | 0.848 | 0.185 | +0.040 | +0.054 | +0.024 | +0.064 | 0.981 |
| GATA2 |  | GATA | 37+18 | 9734/8367 | 9/11 | yes | 0.809 | 0.062 | 0.852 | 0.085 | +0.043 | +0.023 | +0.032 | +0.032 | 0.985 |
| E2F4 |  | E2F | 7+4 | 2686/5636 | 4/3 | partial? | 0.828 | 0.156 | 0.879 | 0.198 | +0.051 | +0.042 | +0.091 | +0.058 | 0.961 |
| SRF |  | MADS SRF | 14+11 | 8305/2084 | 11/3 | no | 0.777 | 0.103 | 0.830 | 0.167 | +0.053 | +0.064 | +0.098 | +0.112 | 0.939 |
| RXRA |  | NR2B RXR | 9+4 | 3918/10000 | 4/9 | yes | 0.764 | 0.182 | 0.818 | 0.254 | +0.054 | +0.073 | +0.079 | +0.007 | 0.963 |
| IRF4 |  | IRF | 12+6 | 6960/10000 | 6/3 | yes | 0.814 | 0.271 | 0.871 | 0.332 | +0.057 | +0.061 | +0.094 | +0.099 | 0.977 |
| STAT1 |  | STAT | 11+6 | 4588/10000 | 5/3 | no | 0.724 | 0.073 | 0.805 | 0.102 | +0.081 | +0.029 | +0.026 | +0.049 | 0.937 |

## Cell types (Sup. Table 1 experiments, cell line from metadata file)

| TF | below | human cell types | mouse cell types | overlap note |
|---|---|---|---|---|
| STA5A | **below** | Treg; K562; Tconv | mammary; T cells; n/a(FH); n/a(MH); Th1; 3T3-L1; MEF | T cells in both |
| ERG | **below** | VCaP; RWPE-1; CD34+ | prostate(ERG-tg); prostate(Pten-/-); spleen; GATA1s megakaryoblasts | prostate (VCaP/RWPE-1 vs ERG-tg prostate); CD34+ vs hematopoietic |
| P53 | **below** | U2OS; keratinocytes; H9-hESC; IMR90; LCL; BJ-hTERT | ESC; MEF | H9-hESC vs ESC; BJ/IMR90 fibroblasts vs MEF |
| TF65 |  | LCL; GM12878; fetal lung fibr.; HUVEC | macrophages; BMDC | LCL/fibroblast/HUVEC vs macrophages/BMDC |
| CTCF |  | LCL; hASC; fibroblast; HUVEC; MCF-7; CD36+ erythroid; GP5d; panc. islets; HMEC; HCT-116; CD14+ monocytes; lung; A549; HFF; CD20+ B; NHLF; HSMMtube; pancreas; GM12878; NHEK; MM1.S; fetal lung fibr.; WI-38; NH-A; HL-60; SK-N-SH; HVMF; HeLa-S3; T-47D; osteoblast; K562; AG09309; HCFaa; HMF; NB4; ECC-1; BJAB; CD34+; BL41; THP-1 | ESC; BMDC; bone marrow; pro-B; cortex; testis; cerebellum; lung; spleen; liver; heart; CH12; E14.5 limb; plasmablasts; DP thymocytes; pre-pro-B; diff. ESC; G1E-ER4; NIH-3T3; placenta; intestine; E14.5 heart; G1E; E14.5 brain; splenic B; DN thymocytes; limb; MEF; sm. intestine; brain; MEL; kidney; C2C12 | many shared lineages (ESC, erythroid, B, fibroblast, liver, lung) |
| USF2 |  | HepG2; GM12878; H1-hESC; K562; HeLa-S3 | MEL; CH12 | K562 vs MEL; GM12878 vs CH12 |
| COE1 |  | GM12878; ASC-adipocyte | pre-B; splenic B; bone marrow; pro-B; 3T3-L1 | GM12878 (B-LCL) vs pre-B/B cells; ASC vs 3T3-L1 |
| CEBPB |  | ECC-1; A549; HCT-116; HepG2; LoVo; monocytes | 3T3-L1; C2C12; n/a; tissue(C57BL6); MC3T3-E1; diff. hematopoietic; liver | HepG2 vs liver; monocytes vs hematopoietic; rest unmatched (adipocyte/myoblast) |
| MAX |  | P493-6; MCF-7; HepG2; HCT-116; H128; HeLa-S3; H1-hESC; NB4; HUVEC; K562; ECC-1; plasma cell; H2171 | C2C12; MEL; ESC | K562 vs MEL; H1-hESC vs ESC |
| FLI1 |  | CB megakaryocytes | HPC-7; megakaryocytes; diff. hematopoietic | megakaryocytes in both |
| PPARG |  | hASC; SGBS adipocytes; THP-1; adipocyte | 3T3-L1; 3T3-L1 adipocytes; WAT adipocytes; BAT adipocytes | adipocytes in both |
| JUND |  | HCT-116; SK-N-SH; GM12878; HepG2; MCF-7; T-47D; LoVo; GP5d | STHdh striatal; CH12; Th17; 3T3-L1 | GM12878 vs CH12 (B); SK-N-SH vs striatal; majority unmatched |
| HNF4A |  | HepG2; Caco-2; GP5d; LoVo | kidney; intestinal villus; intestinal villus(Cdx2-/-) | colon lines (Caco-2/LoVo/GP5d) vs intestinal villus |
| STAT3 |  | HeLa-S3; MCF10A-Src; U-2932 | AtT-20; ESC; CD4 T | HeLa/MCF10A/U-2932 vs AtT-20/ESC/CD4 T |
| GATA1 |  | K562; proerythroblasts; CD36+ erythroid; erythroblasts; CB megakaryocytes | G1E-ER4; G1ME; Ter119+ erythroid; erythroblasts; megakaryocytes; ES-EP; bone marrow; FDCP-mix; MEL | erythroid/megakaryocytic in both |
| TFE2 |  | GM12878; RPMI-8402(T-ALL) | thymocytes; C2C12 | RPMI-8402 (T-ALL) vs thymocytes; GM12878 unmatched |
| MYC |  | K562; P493-6; HeLa-S3; HepG2; GM12878; FB8470; BL41; Ramos; CA46; Blue1; MCF-7; H2171; plasma cell; HeLa; U2OS | ESC; MEL; CH12; pro-B; tissue(p53 mut.) | K562 vs MEL; LCL vs CH12/pro-B; ESC |
| RUNX1 |  | Jurkat; K562; Kasumi-1; CMK; CCRF-CEM | BMDC; diff. hematopoietic; LSK HPC; L8057(megakaryoblast); megakaryocytes; BM iFLT3 | hematopoietic/megakaryocytic leukaemia lines vs hematopoietic progenitors/megakaryocytes |
| IRF1 |  | CD14+ monocytes; K562 | BMDC | CD14+ monocytes vs BMDC (both LPS-stimulated myeloid) |
| REST |  | PANC-1; GM12878; PFSK-1; HepG2; GP5d; HCT-116; ECC-1; K562; U87; H295R; SK-N-SH; HeLa-S3; MCF-7; H1-hESC; HL-60; LoVo; APL blasts | ESC; C2C12 | H1-hESC vs ESC |
| PRGR |  | AB32; T-47D | uterus | T-47D/AB32 vs uterus (both PR-responsive reproductive tissue, not the same tissue) |
| AP2A |  | HeLa-S3; LoVo | epididymis; prostate | HeLa/LoVo vs epididymis/prostate |
| SPI1 |  | OCI-LY3; macrophages; CB CD133+; THP-1; K562; GM12878; OCI-LY10; H929 | macrophages; RAW264.7; B cells; 38B9(pre-B); MEL; HPC-7; DN2b thymocytes; 3T3-L1; FDCP-mix | macrophages/THP-1/B-lymphoma vs macrophages/B cells |
| ESR1 |  | breast tumor; ECC-1; MCF-7; T-47D; MDA-MB-134; Ishikawa; H3396; ZR-75-1 | liver; uterus; mammary 3134; mammary 7438 | breast/endometrium vs mammary/uterus |
| GCR |  | A549; LNCaP-1F5; ECC-1; HEK293 | mammary 3134; 3T3-L1 | A549/LNCaP/ECC-1 vs mammary 3134/3T3-L1 |
| TAL1 |  | proerythroblasts; CD36+ erythroid; T-ALL; CB megakaryocytes; K562; Jurkat; CCRF-CEM | G1E; Lin- HPC; Ter119+ erythroid; FL erythroblasts; G1E-ER4; erythroblasts; MEL; bone marrow | erythroid in both |
| ANDR |  | LNCaP; VCaP; LNCaP-1F5; DU145; LNCaP-AR; C4-2B; prostate ca. | epididymis; prostate; prostate(ERG-tg); prostate(Pten-/-); kidney | prostate (LNCaP/VCaP vs mouse prostate) |
| MAFK |  | K562; LoVo | MEL; CH12 | K562 vs MEL (erythroleukemia); LoVo/CH12 unmatched |
| GATA3 |  | MCF-7; T-47D; Jurkat; Th1; Th2; CCRF-CEM | CD4 T; DN1 thymocytes; DN2b thymocytes; DP thymocytes | T cells in both (Th1/Th2/Jurkat vs CD4 T/thymocytes) |
| SOX2 |  | ReN-VM(NPC); KYSE70; TT; HCC95; H9-hESC | ESC; diff. ESC; TS cells; ES-NPC; neural EB | H9-hESC/ReN-VM vs ESC/NPC |
| GATA2 |  | K562; HUVEC; CB megakaryocytes; CD34+ | FDCP-mix; G1E; Lin- HPC; HPC-7; G1ME | hematopoietic progenitors/K562 vs FDCP-mix/G1E |
| E2F4 |  | LCL; MSC-adipocyte; HeLa-S3; K562 | MEF(OSK reprogr.); CH12; C2C12 | LCL vs CH12 (B lineage); others unmatched (HeLa/K562 vs MEF/C2C12) |
| SRF |  | GM12878; HCT-116; MCF-7; HepG2; H1-hESC; K562; ECC-1 | HL-1; NIH-3T3; C2C12 | LCL/epithelial lines vs cardiomyocyte/fibroblast/myoblast |
| RXRA |  | LS180; HepG2 | F9 EC; liver; RMWT_227; IDG-SW3(osteocyte) | HepG2 vs liver |
| IRF4 |  | OCI-LY3; GM12878; OCI-LY10; H929 | B cells; pro-B | B-cell lymphoma/LCL/myeloma vs B/pro-B cells |
| STAT1 |  | K562; HeLa-S3 | BMDM | K562/HeLa vs BMDM |

## Statistics

| Spearman (n=36) | rho | p |
|---|---|---|
| dauROC_M vs # human training positives | -0.091 | 0.596 |
| dauROC_M vs # mouse training positives | -0.193 | 0.258 |
| dauROC_M vs # mouse experiments (Sup. Table 1) | -0.250 | 0.141 |
| dauROC_M vs # human experiments (Sup. Table 1) | -0.008 | 0.963 |
| dauROC_M vs # human monoPWMs | -0.060 | 0.728 |
| dauROC_M vs best-PWM auROC on mouse (ceiling) | -0.446 | 0.00639 |
| dauROC_M vs dauROC_H (within-human gain) | +0.689 | 3.37e-06 |
| dauROC_M vs motif sim. baseline human PWM vs nearest mouse PWM | -0.286 | 0.0907 |
| dauROC_M vs motif sim. median(human PWM -> nearest mouse PWM) | -0.276 | 0.103 |
| dauROC_M vs # paralogs (approx.) | -0.092 | 0.595 |
| dauPRC_M vs # human training positives | +0.354 | 0.0341 |
| dauPRC_M vs # mouse training positives | +0.188 | 0.273 |
| dauPRC_M vs # mouse experiments (Sup. Table 1) | +0.001 | 0.996 |
| dauPRC_M vs # human experiments (Sup. Table 1) | +0.046 | 0.791 |
| dauPRC_M vs # human monoPWMs | -0.129 | 0.454 |
| dauPRC_M vs best-PWM auPRC on mouse (ceiling) | +0.163 | 0.341 |
| dauPRC_M vs dauPRC_H (within-human gain) | +0.416 | 0.0116 |
| dauPRC_M vs motif sim. baseline human PWM vs nearest mouse PWM | +0.207 | 0.226 |
| dauPRC_M vs motif sim. median(human PWM -> nearest mouse PWM) | +0.125 | 0.468 |
| dauPRC_M vs # paralogs (approx.) | -0.243 | 0.153 |

| Mann-Whitney (3 below vs 33 others) | median (below) | median (others) | p |
|---|---|---|---|
| # human training positives | 2.55e+03 | 7.02e+03 | 0.080 |
| # mouse training positives | 7.11e+03 | 1e+04 | 0.207 |
| # human experiments | 6 | 9 | 0.841 |
| # mouse experiments | 4 | 5 | 0.931 |
| best-PWM auROC on mouse | 0.807 | 0.862 | 0.409 |
| best-PWM auPRC on mouse | 0.385 | 0.182 | 0.746 |
| dauROC_H | 0.0278 | 0.0239 | 0.914 |
| baseline-human vs nearest-mouse PWM PCC | 0.946 | 0.97 | 0.067 |
| median human->mouse nearest PCC | 0.963 | 0.964 | 0.864 |
| # paralogs (approx.) | 7 | 4 | 0.282 |
