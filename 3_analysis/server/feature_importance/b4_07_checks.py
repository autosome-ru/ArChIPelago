"""b4_07_checks.py - feature importance, step 7: checks that were made to fail once. Writes checks/B4_checks.csv (observed values printed, not 'OK').
 c1  column order of the loaded model: test auROC of the model vs the independent retrain scores in the score npz <TF>.npz
     (rf_H_s0, rf_H_s1); a wrong column order drops the auROC towards 0.5. Also test sequence names identical to the npz names.
 c2  Gini importances sum to 1 and number of features equals model.n_features_in_ (from the gini csv)
 c3  toy example for the strong-hit rule (b4_03.pvals): one value that must be strong, one that must not
 c4  wrong-order control: the same model with columns reversed for SRF (must give a much lower auROC, proves c1 can fail)
"""
import os, sys, joblib
from b4_common import *

roc_auc_score = ad.prroc_metrics.auroc  # full auROC by PRROC (y, s)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b4_03_usage import pvals

rows = []
for tf in TFS:
    C = np.load(os.path.join(CACHE, "%s_test.npz" % tf), allow_pickle=True)
    d = np.load(AE.work("H5_alltrain", "human_sym", "scores", "%s.npz" % tf), allow_pickle=True)  # same model -> identical auROC
    a_model = roc_auc_score(C["y"], C["rf"])
    a0, a1 = roc_auc_score(d["y_test_H"], d["rf_H_s0"]), roc_auc_score(d["y_test_H"], d["rf_H_s1"])
    G = pd.read_csv(os.path.join(OUT, "gini", "%s_gini.csv" % tf))
    rows.append(dict(TF=tf, same_names=bool((d["names_test_H"] == C["names"]).all()), auroc_model=a_model, auroc_retrain_s0=a0,
                     auroc_retrain_s1=a1, diff_vs_mean_retrain=a_model - (a0 + a1) / 2, gini_sum=G.gini.sum(), n_features=len(G),
                     n_cache_features=len(C["feats"])))
R = pd.DataFrame(rows)
# toy
Xt = np.r_[np.arange(100.0), [99.5, 50.0]][:, None]; yt = np.r_[np.zeros(100, int), [1, 1]]
P, _ = pvals(Xt.astype(np.float32), yt)
toy = dict(p_value_of_99p5=float(P[100, 0]), p_value_of_50=float(P[101, 0]), strong_99p5=bool(P[100, 0] <= 0.01), strong_50=bool(P[101, 0] <= 0.01))
# wrong-order control (SRF)
m = joblib.load(model_path("SRF"))
C = np.load(os.path.join(CACHE, "SRF_test.npz"), allow_pickle=True)
Xs = scale_test(C["X"], "SRF", list(C["feats"])); wrong = roc_auc_score(C["y"], m.predict_proba(Xs[:, ::-1])[:, 1])
os.makedirs(os.path.join(OUT, "checks"), exist_ok=True)
R.to_csv(os.path.join(OUT, "checks", "B4_checks_per_TF.csv"), index=False)
pd.DataFrame([dict(**toy, SRF_auroc_correct_order=float(R.set_index("TF").loc["SRF", "auroc_model"]), SRF_auroc_reversed_columns=wrong)]
             ).to_csv(os.path.join(OUT, "checks", "B4_checks_toy_and_wrong_order.csv"), index=False)
print(R.round(4).to_string()); print(toy, "SRF reversed-column auROC", wrong)
