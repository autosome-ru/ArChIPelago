### Median delta (RF - best single mono PWM, human-train-selected) vs number of PWMs k

Budget convention: TFs with P < k enter at k_eff = P; all 36 TFs at every k. Values: median over TFs of the per-TF mean over replicates (random design; top-k design in brackets).

| k | auROC, human test | auPRC, human test | auROC, mouse test | auPRC, mouse test |
|---|---|---|---|---|
| 1 | -0.0189 [+0.0058] | -0.0391 [+0.0026] | -0.0222 [+0.0011] | -0.0357 [-0.0068] |
| 2 | -0.0014 [+0.0101] | -0.0097 [+0.0142] | -0.0018 [+0.0078] | -0.0073 [+0.0118] |
| 4 | +0.0081 [+0.0145] | +0.0120 [+0.0146] | +0.0073 [+0.0103] | +0.0037 [+0.0126] |
| 8 | +0.0137 [+0.0181] | +0.0224 [+0.0232] | +0.0117 [+0.0121] | +0.0145 [+0.0151] |
| 16 | +0.0211 [+0.0205] | +0.0284 [+0.0284] | +0.0147 [+0.0136] | +0.0195 [+0.0208] |
| 32 | +0.0235 [+0.0219] | +0.0304 [+0.0314] | +0.0162 [+0.0155] | +0.0229 [+0.0214] |
| 64 | +0.0244 [+0.0236] | +0.0356 [+0.0366] | +0.0155 [+0.0153] | +0.0231 [+0.0210] |
| 128 | +0.0244 [+0.0244] | +0.0371 [+0.0362] | +0.0162 [+0.0147] | +0.0231 [+0.0210] |
| all (P) | +0.0244 [+0.0244] | +0.0371 [+0.0371] | +0.0162 [+0.0162] | +0.0210 [+0.0210] |

### k at which the median delta reaches 50 % / 90 % of the median full (k = P) delta

| metric | median full delta | 50 % (random) | 90 % (random) | 50 % (top-k) | 90 % (top-k) |
|---|---|---|---|---|---|
| auROC, human test | +0.0244 | 8 | 32 | 4 | 64 |
| auPRC, human test | +0.0371 | 8 | 64 | 8 | 64 |
| auROC, mouse test | +0.0162 | 8 | 16 | 4 | 32 |
| auPRC, mouse test | +0.0210 | 8 | 16 | 2 | 16 |

### Per-TF saturation k (first grid k at which the per-TF mean delta >= 90 % of that TF's full delta; TFs with full delta <= 0 excluded)

| metric | n TFs with full delta > 0 | median sat. k (random) | IQR | share saturated by k = 8 | share by k = 16 | median sat. k (top-k) | TFs where more PWMs hurt (> 0.005) |
|---|---|---|---|---|---|---|---|
| auROC, human test | 36 | 16 | 16-32 | 11 % | 50 % | 16 | 1/36 (P53) |
| auPRC, human test | 35 | 32 | 16-64 | 14 % | 49 % | 32 | 1/36 (PRGR) |
| auROC, mouse test | 35 | 16 | 16-32 | 20 % | 57 % | 16 | 3/36 (ERG, P53, TF65) |
| auPRC, mouse test | 34 | 16 | 8-32 | 26 % | 53 % | 16 | 6/36 (CTCF, ERG, ESR1, JUND, P53, USF2) |

(sat. k = 512 stands for 'only at k = P')

### Spearman correlation between P (number of PWMs) and the full-model gain over the best single mono PWM

| metric | rho | p |
|---|---|---|
| auROC, human test | -0.223 | 0.191 |
| auPRC, human test | 0.267 | 0.116 |
| auROC, mouse test | -0.084 | 0.626 |
| auPRC, mouse test | -0.103 | 0.551 |

### Top-k (train-auROC-ranked) minus random subsets, median over TFs of the delta difference

| k | auROC, human test | auPRC, human test | auROC, mouse test | auPRC, mouse test |
|---|---|---|---|---|
| 1 | +0.0287 (n=36) | +0.0315 (n=36) | +0.0229 (n=36) | +0.0262 (n=36) |
| 2 | +0.0135 (n=36) | +0.0267 (n=36) | +0.0094 (n=36) | +0.0139 (n=36) |
| 4 | +0.0067 (n=36) | +0.0109 (n=36) | +0.0028 (n=36) | +0.0068 (n=36) |
| 8 | +0.0020 (n=35) | +0.0046 (n=35) | +0.0010 (n=35) | -0.0006 (n=35) |
| 16 | +0.0004 (n=27) | +0.0022 (n=27) | -0.0009 (n=27) | +0.0003 (n=27) |
| 32 | -0.0003 (n=18) | -0.0000 (n=18) | -0.0007 (n=18) | -0.0011 (n=18) |
| 64 | +0.0005 (n=8) | -0.0007 (n=8) | +0.0002 (n=8) | +0.0039 (n=8) |
| 128 | -0.0034 (n=4) | -0.0043 (n=4) | -0.0020 (n=4) | -0.0028 (n=4) |

### RF seed-to-seed spread at k = P (|seed0 - seed1|, median / max over TFs)

| metric | median | max |
|---|---|---|
| auROC, human test | 0.0007 | 0.0091 |
| auPRC, human test | 0.0023 | 0.0140 |
| auROC, mouse test | 0.0005 | 0.0097 |
| auPRC, mouse test | 0.0016 | 0.0142 |

### At k = 1 and k = 2, RF on the subset vs the subset's best single PWM (median over TFs of RF - single)

| k | auROC, human test | auPRC, human test | auROC, mouse test | auPRC, mouse test |
|---|---|---|---|---|
| 1 | -0.0035 | -0.0084 | -0.0042 | -0.0133 |
| 2 | +0.0065 | +0.0008 | +0.0019 | -0.0025 |
| 4 | +0.0117 | +0.0097 | +0.0065 | +0.0060 |

Fits: 1272 rows (1000 random-subset, 200 top-k, 72 full); total RF fit+predict time 4.2 CPU-worker hours.
