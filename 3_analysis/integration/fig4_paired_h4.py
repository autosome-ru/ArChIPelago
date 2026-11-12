"""Paired per-TF comparison of Fig. 4 (human test set) and Fig. S4 (mouse test set): ArChIPelago RF on all PWMs vs each Slim /
diChIPMunk model. Inputs: the Fig. 4 / S4 source data (h3_inputs fig4_src, figS4_src =
results/server_h4_test_sets/human/fig4_refit_trainscale). Per test set and comparator: number of TFs, TFs with RF above the comparator
(strictly), ties, median paired difference (RF minus comparator) and the two-sided Wilcoxon signed-rank P (scipy default, as in
NUMBERS), for auROC and auPRC.
Output: 4_visualisation/main_figures/Figure_4/Figure_4_S4_paired_RF_vs_Slim.csv (Sup. Table 3 sheet 'Fig4_S4_paired_counts').
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/fig4_paired_h4.py"""
import os
import sys

import pandas as pd
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h3_inputs as H   # noqa: E402
import w2   # noqa: E402

RF = "ArChIPelago RF on all PWMs"
COMPARATORS = ("Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk")
OUT = os.path.join(os.path.dirname(HERE), "../4_visualisation/main_figures/Figure_4", "Figure_4_S4_paired_RF_vs_Slim.csv")


def main():
    src = {"human": H.IN["fig4_src"], "mouse": H.IN["figS4_src"]}
    rows = []
    for ts, path in src.items():
        d = pd.read_csv(path)
        if set(d.test_set) != {ts}:
            raise SystemExit("CHECK FAILED %s: test_set column %s" % (path, sorted(set(d.test_set))))
        rf = d[d.model == RF].set_index("TF")
        for cmp_ in COMPARATORS:
            c = d[d.model == cmp_].set_index("TF")
            if sorted(c.index) != sorted(rf.index) or len(rf) != 35:
                raise SystemExit("CHECK FAILED %s: %s and RF cover different TFs" % (path, cmp_))
            row = dict(test_set=ts, figure="Fig. 4" if ts == "human" else "Fig. S4", model=RF, comparator=cmp_, n_TF=len(rf))
            for m in ("auROC", "auPRC"):
                dd = rf[m] - c.loc[rf.index, m]
                row.update({"n_TF_RF_above_%s" % m: int((dd > 0).sum()), "n_TF_tied_%s" % m: int((dd == 0).sum()),
                            "median_difference_%s" % m: float(dd.median()), "wilcoxon_P_two_sided_%s" % m: float(wilcoxon(dd).pvalue)})
            rows.append(row)
    out = pd.DataFrame(rows)
    w2.backup(OUT)
    out.to_csv(OUT, index=False)
    w2.provenance(os.path.dirname(OUT), "3_analysis/integration/fig4_paired_h4.py", [OUT], list(src.values()),
                  params="{difference: RF on all PWMs minus comparator; above: strictly greater; P: scipy wilcoxon default, two-sided}")
    print(out.to_string())


if __name__ == "__main__":
    main()
