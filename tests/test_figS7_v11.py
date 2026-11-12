"""Fig. S7 on HOCOMOCO v11: every tile label of the figure equals the nearest v11 model in
4_visualisation/supplementary_figures/Figure_S7/hocomoco11_check/subtype_vs_hocomoco11.csv (own TF on every tile; any TF on the co-factor tiles), the source
data carries the same v11 columns, and the figure draws STAT1, TAL1 (v11 subtypes) and IRF4, TFE2 (own motif plus co-factor).
Red: T1_MUTATE=figS7_v11 swaps the nearest own-TF model of STAT1 panel 3 (STAT1_HUMAN.H11MO.1.A) for subtype 0 in the csv.
"""
import os
import re
import subprocess
from collections import Counter

import numpy as np
import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
S7 = os.path.join(T.R2, "4_visualisation", "supplementary_figures", "Figure_S7")
FIGS = [os.path.join(S7, "Figure_S7.pdf")]
LABEL = re.compile(r"(v11|any TF:) ([A-Z0-9]+_HUMAN\.H11MO\.\d\.[A-D]), r (\d\.\d{3})")
KEY = ["TF", "panel", "pwm_name"]


def v11():
    d = pd.read_csv(os.path.join(S7, "hocomoco11_check", "subtype_vs_hocomoco11.csv"))
    if T.mutate("figS7_v11"):
        d.loc[(d.TF == "STAT1") & (d.panel == 3), "best_own_v11"] = "STAT1_HUMAN.H11MO.0.A"
    return d


def text(pdf):
    return subprocess.run([os.environ.get("PDFTOTEXT", "pdftotext"), "-layout", pdf, "-"], capture_output=True, text=True, check=True).stdout


def expected(gap=0.2):
    v = v11()
    v = v[v.TF.isin(["STAT1", "TAL1", "IRF4", "TFE2"])]
    lab = [("v11", r.best_own_v11, "%.3f" % r.r_best_own) for _, r in v.iterrows()]
    cof = v[(v.best_any_v11.str.split("_HUMAN").str[0] != v.TF) & (v.r_best_any - v.r_best_own >= gap)]
    return lab + [("any TF:", r.best_any_v11, "%.3f" % r.r_best_any) for _, r in cof.iterrows()], len(v), len(cof)


def test_label_regex_toy():
    assert LABEL.findall("v11 TAL1_HUMAN.H11MO.1.A, r 0.995  any TF: BATF_HUMAN.H11MO.0.A, r 0.975") == [
        ("v11", "TAL1_HUMAN.H11MO.1.A", "0.995"), ("any TF:", "BATF_HUMAN.H11MO.0.A", "0.975")]


@pytest.mark.parametrize("pdf", FIGS, ids=["Figure_S7"])
def test_every_tile_label_equals_the_nearest_v11_model(pdf):
    want, n_tiles, n_cof = expected()
    # tile and co-factor counts as the NUMBERS legend row states them (4 TFs x (4 mono + 3 di) = 28 tiles, 3 co-factor tiles)
    n = pd.read_csv(os.path.join(T.R2, "outputs/numbers/NUMBERS_R2_35TF.csv")).set_index("item")["new value"]
    leg = n["Fig. S7 legend counts: monoPWMs / diPWMs per TF; candidate TFs; TFs drawn; tiles drawn; co-factor tiles"].split("; ")
    assert (n_tiles, n_cof) == (int(leg[3]), int(leg[4])), (n_tiles, n_cof, leg)
    got = LABEL.findall(text(pdf))
    assert Counter(got) == Counter(want), (sorted(set(got) - set(want)), sorted(set(want) - set(got)))


def test_source_data_v11_columns_equal_the_csv():
    s = pd.read_csv(os.path.join(S7, "Figure_S7_source_data.csv"))
    m = s.merge(v11(), on=KEY, suffixes=("", "_csv"), validate="1:1")
    assert len(m) == len(s)
    bad = [(r.TF, r.panel, c) for _, r in m.iterrows() for c in ("best_own_v11", "best_any_v11") if r[c] != r[c + "_csv"]]
    bad += [(r.TF, r.panel, c) for _, r in m.iterrows() for c in ("r_best_own", "r_best_any") if abs(r[c] - r[c + "_csv"]) > 5e-4]
    assert not bad, bad


@pytest.mark.parametrize("pdf", FIGS, ids=["Figure_S7"])
def test_drawn_tfs_are_the_v11_examples(pdf):
    drawn = re.findall(r"^\s*[A-H]\s*\n\s*([A-Z0-9]+)\b", text(pdf), re.M)   # the TF name follows the panel letter
    v = pd.read_csv(os.path.join(S7, "subtype_verdicts_v11_h3.csv"))
    assert sorted(v.TF[v.drawn]) == ["IRF4", "STAT1", "TAL1", "TFE2"]
    assert drawn == ["STAT1", "TAL1", "IRF4", "TFE2"], drawn


def test_cluster_vs_v11_uses_the_ablation_cuts_and_its_numbers_recompute():
    """The clustering of cluster_vs_v11 is the ablation subtype (cut per TF from B4_ablation_summary.csv: 0.8, MAFK and
    SPI1 0.9), Rand index and ARI per TF recompute from the per-PWM table, informative = both partitions with >= 2 groups, and the
    summary rows are the median Rand over the 11 TFs and the median ARI over the informative TFs.
    Red: T1_MUTATE=cluster_cut reads every TF at the 0.8 cut."""
    from itertools import combinations
    from sklearn.metrics import adjusted_rand_score
    H = os.path.join(S7, "hocomoco11_check")
    per = pd.read_csv(os.path.join(H, "cluster_vs_v11_per_PWM.csv"))
    tab = pd.read_csv(os.path.join(H, "cluster_vs_v11.csv"))
    tab = tab.assign(key=tab.TF.fillna(tab.summary)).set_index("key")
    cut = pd.read_csv(os.path.join(T.R2, "3_analysis/feature_importance", "tables", "B4_ablation_summary.csv")).groupby("TF").sim_threshold.first()
    if T.mutate("cluster_cut"):
        cut[:] = 0.8
    sub = pd.read_csv(os.path.join(T.R2, "3_analysis/feature_importance", "tables", "B4_feature_subtypes.csv")).set_index(["TF", "feature"])
    assert len(cut) == 11 and sorted(cut.index[cut == 0.9]) == ["MAFK", "SPI1"]
    want = [sub.loc[(r.TF, r.feature), "subtype_%d" % round(100 * cut[r.TF])] for r in per.itertuples()]
    assert (per.cluster.values == np.array(want)).all(), "per-PWM clusters are not the ablation subtypes"
    rows = {}
    for tf, d in per.groupby("TF"):
        a, b = d.cluster.values, d.v11_subtype.values
        pr = list(combinations(range(len(a)), 2))
        rows[tf] = (sum((a[i] == a[j]) == (b[i] == b[j]) for i, j in pr) / len(pr), adjusted_rand_score(a, b),
                    len(set(a)) >= 2 and len(set(b)) >= 2)
    r = pd.DataFrame(rows, index=["rand", "ari", "inf"]).T
    assert len(r) == 11 and len(per) == 265
    assert np.allclose(tab.loc[r.index, "rand_index"].astype(float), r.rand.astype(float))
    assert np.allclose(tab.loc[r.index, "adjusted_rand"].astype(float), r.ari.astype(float))
    assert (tab.loc[r.index, "informative"].astype(str) == r.inf.astype(str)).all()
    inf = r[r.inf.astype(bool)]
    assert abs(tab.rand_index["median, all TFs"] - r.rand.astype(float).median()) < 1e-12
    assert abs(tab.adjusted_rand["median, informative TFs (%d)" % len(inf)] - inf.ari.astype(float).median()) < 1e-12


def test_named_reference_tfs_recompute_and_their_numbers_rows():
    """subtype_vs_named_v11.csv: for every Fig. S7 matrix, r to the best model of ETV2, ETV7, BATF and FLI1
    recomputes with pcc_max from the HOCOMOCO v11 pcms, quality letters equal the v11 annotation, gap = r_ref - r_own with r_own
    from subtype_vs_hocomoco11.csv, and the NUMBERS rows for TFE2 mono 1 / di 1 quote the csv (ETV2 0.903 / 0.939).
    Red: T1_MUTATE=named_v11 reads ETV2_HUMAN.H11MO.0.B with quality A (and r of TFE2 mono 1 + 0.01)."""
    import sys
    H = os.path.join(S7, "hocomoco11_check")
    sys.path.insert(0, H)
    import subtype_vs_hocomoco11 as V
    d = pd.read_csv(os.path.join(H, "subtype_vs_named_v11.csv"))
    ann = pd.read_csv(os.path.join(T.R2, "results/datasets/hocomoco11", "HOCOMOCOv11_full_annotation_HUMAN_mono.tsv"), sep="\t").set_index("Model")["Quality"]
    if T.mutate("named_v11"):
        ann = ann.copy()
        ann["ETV2_HUMAN.H11MO.0.B"] = "A"
        d.loc[(d.TF == "TFE2") & (d.panel == 1), "r_ETV2"] += 0.01
    own = v11().set_index(KEY).r_best_own
    src = pd.read_csv(os.path.join(S7, "Figure_S7_source_data.csv"))
    assert len(d) == len(src) and set(zip(d.TF, d.panel)) == set(zip(src.TF, src.panel))
    pcm = V.read_pcms(V.PCM)
    bad = []
    for ref in ("ETV2", "ETV7", "BATF", "FLI1"):
        z = {k: V.zrows(m) for k, m in pcm.items() if k.startswith(ref + "_HUMAN.")}
        for r, s in zip(d.itertuples(), d.merge(src, on=KEY, how="left").matrix_file):
            za = V.zrows(V.read_matrix(os.path.join(S7, s)))
            rr = {k: V.pcc_max(za, zz, V.MIN_OVERLAP)[0] for k, zz in z.items()}
            b = max(rr, key=rr.get)
            if getattr(r, "best_" + ref) != b or abs(getattr(r, "r_" + ref) - rr[b]) > 5e-4:
                bad.append((r.TF, r.panel, ref, "r"))
            if getattr(r, "quality_" + ref) != ann[b]:
                bad.append((r.TF, r.panel, ref, "quality"))
            if abs(getattr(r, "gap_" + ref) - (getattr(r, "r_" + ref) - own[(r.TF, r.panel, r.pwm_name)])) > 1.5e-3:
                bad.append((r.TF, r.panel, ref, "gap"))
    assert not bad, bad
    t = d[d.TF == "TFE2"].set_index("panel")
    assert (t.loc[1, "quality_ETV2"], t.loc[1, "quality_ETV7"]) == ("B", "D")
    n = pd.read_csv(os.path.join(T.R2, "outputs/numbers/NUMBERS_R2_35TF.csv")).set_index("item")["new value"]
    assert n["Fig. S7 TFE2: ETS-like PWMs (mono 1, di 1; RPMI-8402), r with ETV2_HUMAN.H11MO.0.B (quality B)"] == "%.3f and %.3f" % (
        t.loc[1, "r_ETV2"], t.loc[5, "r_ETV2"])
    assert n["Fig. S7 TFE2 mono 1: co-factor gap, r with ETV2 model minus r with TFE2 model / ETV7 model minus TFE2 model"] == "%.3f / %.3f" % (
        t.loc[1, "gap_ETV2"], t.loc[1, "gap_ETV7"])
