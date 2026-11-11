"""assemble_fig4_s4.py - Figure_4_source_data.csv (human test set, SP=H) or Figure_S4_source_data.csv (mouse test set,
SP=M) from per-TF refit rows plus an hm_mm_results.csv.
Mapping, as in the package assembler analysis/slim_dichipmunk/assemble_Figure_4_S4.py (mouse half):
  base_auROC = HM_broc_test_M_auroc (best monoPWM selected by TRAIN auROC; = results_table roc_auc_test_M_PWM)
  base_auPRC = HM_bprc_test_M_auprc (best monoPWM selected by TRAIN auPRC; = results_table pr_auc_test_M_PWM)
  slim0/slim1/lslim5/munk = single_<m>_test_M of the rows; RF2f family = <m>.test_M of the rows
  RF_all = HM_rf_test_M_{auroc,auprc}_s0 (RF mono+di, random_state 0; = results_table RandomForest mono+di test_M)
The script also asserts rows best_mono_test_M auROC == HM_broc_test_M_auroc (same check as the package assembler).
Single models: single_<m>_test_<SP> of the rows (scan-based; for human = the notebook 2 table values).
usage: python assemble_fig4_s4.py ROWS_DIR HM_MM_RESULTS OUT_CSV N_TF SP    (refuses to overwrite OUT_CSV)"""
import glob, json, os, sys
import pandas as pd

ROWS, HM, OUT, N_TF, SP = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
assert SP in ("H", "M")
SET = {"H": "human", "M": "mouse"}[SP]
if os.path.exists(OUT):
    raise SystemExit("exists: " + OUT)
SINGLE = ("slim0", "slim1", "lslim5", "munk")
RF2F = ("RF2f", "RF2f_munk", "RF2f_slim1", "RF2f_lslim5", "RF2f_all5")
MODELS = list(SINGLE) + list(RF2F) + ["RF_all"]
LABEL = {"slim0": "Slim m=0", "slim1": "Slim m=1", "lslim5": "LSlim m=-5", "munk": "diChIPMunk",
         "RF2f": "ArChIPelago RF2f", "RF2f_munk": "ArChIPelago RF2f + diChIPMunk", "RF2f_slim1": "ArChIPelago RF2f + Slim m=1",
         "RF2f_lslim5": "ArChIPelago RF2f + LSlim m=-5",
         "RF2f_all5": "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk", "RF_all": "ArChIPelago RF on all PWMs"}
rows = {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(os.path.join(ROWS, "*.json")))}
hm = pd.read_csv(HM).set_index("TF")
tfs = sorted(rows)
if len(tfs) != N_TF or sorted(hm.index) != tfs:
    raise SystemExit("TF sets differ: rows %d, hm_mm %d, expected %d" % (len(tfs), len(hm), N_TF))
long = []
for tf in tfs:
    r, h = rows[tf], hm.loc[tf]
    base_roc, base_prc = float(h["HM_broc_test_%s_auroc" % SP]), float(h["HM_bprc_test_%s_auprc" % SP])
    assert abs(r["best_mono_test_%s" % SP][0] - base_roc) < 1e-9, (tf, r["best_mono_test_%s" % SP][0], base_roc)
    assert r["best_mono_col"] == h.HM_broc_col, (tf, r["best_mono_col"], h.HM_broc_col)
    vals = {m: tuple(r["single_%s_test_%s" % (m, SP)]) for m in SINGLE}
    vals.update({m: tuple(r[m]["test_%s" % SP]) for m in RF2F})
    vals["RF_all"] = (float(h["HM_rf_test_%s_auroc_s0" % SP]), float(h["HM_rf_test_%s_auprc_s0" % SP]))
    for m in MODELS:
        long.append(dict(TF=tf, test_set=SET, model=LABEL[m], model_key=m, auROC=vals[m][0], auPRC=vals[m][1],
                         base_auROC=base_roc, base_auPRC=base_prc, d_auROC=vals[m][0] - base_roc, d_auPRC=vals[m][1] - base_prc))
long = pd.DataFrame(long)
assert long[["auROC", "auPRC"]].notna().all().all() and len(long) == N_TF * len(MODELS)
long.to_csv(OUT, index=False)
print("wrote", OUT, len(long), "rows")
