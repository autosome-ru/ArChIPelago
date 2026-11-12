"""Values for the [PENDING_H3B:<name>] tokens of the manuscript text producer (Abstract, Results 3.4), computed from the Fig. 4 / S4
stats tables that the figure producer writes (4_visualisation/main_figures/Figure_4, supplementary_figures/Figure_S4) and, for the abstract claim, from the
Fig. 4 source data (paired per TF). Token names: results/datasets/text_tokens/manuscript_text_tokens.txt. One row per token:
the text value, the numbers behind it, the first-submission text value it replaces, and whether the statement of the sentence holds.
Output: 4_visualisation/main_figures/Figure_4/text_tokens_h3b.csv, provenance block in the same folder (_provenance_h3.yaml).
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/fig4_tokens_h3b.py
"""
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402

R2, B13 = w2.R2, w2.B13
D = os.path.join(R2, "4_visualisation", "main_figures", "Figure_4")
IN = {"H": os.path.join(D, "Figure_4_stats.csv"),
      "M": os.path.join(R2, "4_visualisation", "supplementary_figures", "Figure_S4", "Figure_S4_stats.csv"),
      "src": os.path.join(B13, "../package", "Figures", "source_data", "Figure_4_source_data.csv"),
      # token lines of the text producer before the values were filled (the filled producer changes with the fill, so it cannot be an
      # input; tests/test_fig4_h3b.py checks that the filled producer has no token left and quotes each value)
      "text": os.path.join(R2, "results/datasets/text_tokens", "manuscript_text_tokens.txt")}
OUT = os.path.join(D, "text_tokens_h3b.csv")
ALPHA = 0.05
F5, ALL, RF2 = "ArChIPelago RF2f + Slim m=1, LSlim m=-5 + diChIPMunk", "ArChIPelago RF on all PWMs", "ArChIPelago RF2f"
SINGLE = ["Slim m=0", "Slim m=1", "LSlim m=-5", "diChIPMunk"]
SLIM = ["Slim m=0", "Slim m=1", "LSlim m=-5"]
w2.backup(OUT)
H = pd.read_csv(IN["H"]).set_index(["metric", "model"]); M = pd.read_csv(IN["M"]).set_index(["metric", "model"])
for t in (H, M):
    if (t.n != 35).any():
        raise SystemExit("CHECK FAILED Fig. 4 / S4 stats: n != 35")
sgn = lambda x: ("%+.3f" if abs(x) >= 0.0005 else "%+.4f") % x   # a median that rounds to 0.000 keeps its sign digit
med = lambda t, m, k="auROC": t.loc[(k, m), "med"]
archi = [m for m in H.index.get_level_values(1).unique() if m.startswith("ArChIPelago")]
src = pd.read_csv(IN["src"])
wide = {k: src.pivot(index="TF", columns="model", values=k) for k in ("auROC", "auPRC")}
above = {(k, s): int((wide[k][ALL] > wide[k][s]).sum()) for k in ("auROC", "auPRC") for s in SLIM}
rows = []


def add(tok, value, numbers, old, holds, note=""):
    rows.append(dict(token=tok, value=value, numbers=";".join("%.6g" % x for x in numbers), round1_text=old, statement_holds=holds, note=note))


add("Abstract_Slim_claim", "ArChIPelago RF (all PWMs) above %s in %s of 35 TFs by auROC and %s by auPRC (paired, human test set)" % (
        ", ".join(SLIM), ", ".join(str(above[("auROC", s)]) for s in SLIM), ", ".join(str(above[("auPRC", s)]) for s in SLIM)),
    [above[("auROC", s)] for s in SLIM] + [above[("auPRC", s)] for s in SLIM], "(no number)",
    min(above.values()) > 35 / 2, "claim: ArChIPelago outperforms Slim / LSlim models")
for tok, m, k, old in (("Fig4_RF2f_dauROC", RF2, "auROC", "+0.009"), ("Fig4_RF2f_dauPRC", RF2, "auPRC", "+0.004"),
                       ("Fig4_RF2f_diChIPMunk_dauROC", RF2 + " + diChIPMunk", "auROC", "+0.015"),
                       ("Fig4_RF2f_Slim1_dauROC", RF2 + " + Slim m=1", "auROC", "+0.020"),
                       ("Fig4_RF2f_LSlim_dauROC", RF2 + " + LSlim m=-5", "auROC", "+0.018"),
                       ("Fig4_5f_dauROC", F5, "auROC", "+0.023"), ("Fig4_5f_dauPRC", F5, "auPRC", "+0.030"),
                       ("Fig4_RFall_dauROC", ALL, "auROC", "+0.024"), ("Fig4_RFall_dauPRC", ALL, "auPRC", "+0.041")):
    add(tok, sgn(med(H, m, k)), [med(H, m, k)], old, med(H, m, k) > 0, "median delta, human test set (Fig. 4)")
for tok, k, old in (("Fig4_5f_P_auROC", "auROC", "1×10⁻¹¹"), ("Fig4_5f_P_auPRC", "auPRC", "2×10⁻⁹")):
    pv = H.loc[(k, F5), "p_two_sided"]
    add(tok, "%.1e" % pv, [pv], old, pv < ALPHA, "two-sided Wilcoxon signed-rank P, as printed in Fig. 4 (one-sided P %.1e)" % H.loc[(k, F5), "p"])
for tok, m, old in (("FigS4_5f_dauROC", F5, "+0.021"), ("FigS4_RFall_dauROC", ALL, "+0.015")):
    add(tok, sgn(med(M, m)), [med(M, m)], old, med(M, m) > 0, "median delta, mouse test set (Fig. S4)")
lead = "five-feature" if (med(M, F5) > med(M, ALL)) and (med(M, F5, "auPRC") > med(M, ALL, "auPRC")) else (
    "RF (all PWMs)" if (med(M, F5) < med(M, ALL)) and (med(M, F5, "auPRC") < med(M, ALL, "auPRC")) else "mixed")
add("FigS4_order_5f_vs_RFall", lead, [med(M, F5), med(M, ALL), med(M, F5, "auPRC"), med(M, ALL, "auPRC")], "five-feature", lead != "mixed",
    "medians auROC and auPRC, mouse test set: five-feature model vs RF (all PWMs)")
r_ = [med(H, m) for m in SINGLE]; p_ = [med(H, m, "auPRC") for m in SINGLE]
add("Fig4_single_models_sign", "(median ΔauROC %s to %s, ΔauPRC %s to %s)" % (sgn(min(r_)), sgn(max(r_)), sgn(min(p_)), sgn(max(p_))),
    [min(r_), max(r_), min(p_), max(p_)], "(no number)", max(r_) < 0 and sum(x < 0 for x in p_) >= 3,
    "Slim m=0, Slim m=1, LSlim m=-5, diChIPMunk, human test set; 'predominantly negative': all auROC medians and at least three auPRC medians below 0")
pmax = max(H.loc[(k, m), "p_two_sided"] for k in ("auROC", "auPRC") for m in archi)
def _ceil1(x):   # one significant digit, rounded up (an upper bound); 3 decimals at or above 0.001, else x.10^-k
    if x >= 1e-3:
        return "%.3f" % (float("%.1g" % x) if float("%.1g" % x) >= x else x)
    m, k = ("%.0e" % x).split("e")
    m = int(m) + (1 if float("%se%s" % (m, k)) < x else 0)
    m, k = (1, int(k) + 1) if m == 10 else (m, int(k))
    return "%d×10%s" % (m, str(k).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")))
add("Fig4_all_RF_models_P", "≤ " + _ceil1(pmax), [pmax], "< 0.05", pmax < ALPHA,
    "largest two-sided Wilcoxon P over the %d ArChIPelago models and both metrics, human test set" % len(archi))
sl = ["%s %s / %s" % (m, sgn(med(H, m)), sgn(med(H, m, "auPRC"))) for m in SLIM]
add("Fig4_Slim_sign", "(median ΔauROC / ΔauPRC: %s)" % "; ".join(sl), [x for m in SLIM for x in (med(H, m), med(H, m, "auPRC"))],
    "(no number)", all(med(H, m) <= 0.005 and med(H, m, "auPRC") <= 0.01 for m in SLIM),
    "'negative or near-zero': every Slim median at most +0.005 auROC and +0.01 auPRC")
add("Fig4_diChIPMunk_sign", "(median ΔauROC %s, ΔauPRC %s; %d and %d of 35 TFs above 0)" % (
        sgn(med(H, "diChIPMunk")), sgn(med(H, "diChIPMunk", "auPRC")), H.loc[("auROC", "diChIPMunk"), "n_pos"], H.loc[("auPRC", "diChIPMunk"), "n_pos"]),
    [med(H, "diChIPMunk"), med(H, "diChIPMunk", "auPRC"), H.loc[("auROC", "diChIPMunk"), "n_pos"], H.loc[("auPRC", "diChIPMunk"), "n_pos"]],
    "(no number)", med(H, "diChIPMunk") <= 0 and med(H, "diChIPMunk", "auPRC") <= 0, "failed to outperform the best monoPWM")
t = pd.DataFrame(rows)
toks = sorted(set(re.findall(r"\[PENDING_H3B:([A-Za-z0-9_]+)\]", open(IN["text"]).read())))
if sorted(t.token) != toks:
    raise SystemExit("CHECK FAILED tokens in the text %s differ from the computed ones %s" % (toks, sorted(t.token)))
t.to_csv(OUT, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(t[["token", "value", "round1_text", "statement_holds"]].to_string(index=False))
w2.provenance(D, "3_analysis/integration/fig4_tokens_h3b.py", [OUT], list(IN.values()), params="{alpha: %.2f, tokens: %d}" % (ALPHA, len(t)))
