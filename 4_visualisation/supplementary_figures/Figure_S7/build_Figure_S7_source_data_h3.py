"""Fig. S7 source data on the final data. Same method as the package builder
(3_analysis/package/analysis/motif_subtypes/build_Figure_S7_source_data.py, imported for to_prob / consensus / pcc_max / ...),
with the final inputs:
  ranking      3_analysis/package/results_table.csv: RF mono+di seed 0 vs *_PWM (best monoPWM, asserted), human test set
  importances  h3_inputs b4_dir/gini/<TF>_gini.csv (Gini of the seed-0 mono+di model)
  PWM files    results/datasets/zenodo_v1.0.0/PWMs_{mono,di}_HUMAN (identical to the Zenodo v1.0.0 tars)
  counts       outputs/supplementary_tables Sup_Table_2 per_TF (Human monoPWMs / diPWMs) must equal the features of each model
Candidates = union of the top N_TOP by dauROC and by dauPRC (rule of the package builder); every candidate gets its panels.
HOCOMOCO v11 (subtype labels come from the PWM collection, not from logo reading): every panel matrix is matched
to the 769 HOCOMOCO v11 HUMAN mono models with pcc_max (hocomoco11_check/subtype_vs_hocomoco11.py, imported): nearest own-TF
model (its subtype index), nearest model of any TF, and the co-factor flag (any-TF r exceeds own-TF r by >= COFACTOR_GAP).
Figure TFs (FIG_TFS): subtype examples STAT1 and TAL1, own motif plus co-factor motif IRF4 and TFE2; the role
of each is checked against the v11 columns (>= 2 own v11 subtypes nearest without a tie / >= 1 co-factor panel).
Outputs (this folder): tf_ranking_h3.csv, Figure_S7_source_data.csv (with the v11 columns), subtype_verdicts_v11_h3.csv,
matrices/<TF>_<feature>.txt, _provenance_h3.yaml.
usage: ARCHI_STAGE=H5 ARCHI_ZENODO_DIR=results/datasets/zenodo_v1.0.0 python3 4_visualisation/supplementary_figures/Figure_S7/build_Figure_S7_source_data_h3.py
"""
import glob
import importlib.util
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

os.environ.setdefault("ARCHI_ZENODO_DIR", os.path.join(R2, "results/datasets", "zenodo_v1.0.0"))
_spec = importlib.util.spec_from_file_location("s7pkg", os.path.join(R2, "3_analysis", "package", "analysis", "motif_subtypes",
                                                                     "build_Figure_S7_source_data.py"))
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)
_v = importlib.util.spec_from_file_location("v11", os.path.join(HERE, "hocomoco11_check", "subtype_vs_hocomoco11.py"))
V11 = importlib.util.module_from_spec(_v)
_v.loader.exec_module(V11)
FIG_TFS = {"STAT1": "subtypes", "TAL1": "subtypes", "IRF4": "cofactor", "TFE2": "cofactor"}   # row order of the figure
COFACTOR_GAP = 0.2   # co-factor PWM: nearest v11 model of any TF beats the nearest own-TF model by at least this r
TIE = 0.01           # two own-TF v11 subtypes closer than this in r are a tie
PWM_ROOT = P.PWM_ROOT
N_TOP, N_MONO = P.N_TOP, P.N_MONO
# The top three diPWMs of every TF (one uniform rule). In the final data the AP-1 diPWM of IRF4 is
# its third diPWM (Gini 0.14918 against 0.14935 for the second), so two diPWMs would drop IRF4's co-factor panel.
# Override for a check run: ARCHI_S7_N_DI=<n>. The R figure reads the per-TF counts from the source data.
N_DI = int(os.environ.get("ARCHI_S7_N_DI", "3"))
IN = {"rt": os.path.join(R2, "3_analysis", "package", "results_table.csv"),
      "gini": os.path.join(H.IN["b4_dir"], "gini"),
      "t2": os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx"),
      "meta": P.META, "v11": V11.PCM}


def v11_match(tf, p, z11):
    """Nearest HOCOMOCO v11 models of one probability matrix (L x 4), as in subtype_vs_hocomoco11.py."""
    za = P.zrows(p)
    own = {k: V11.pcc_max(za, z, V11.MIN_OVERLAP)[0] for k, z in z11.items() if k.startswith(tf + "_HUMAN.")}
    if not own:
        raise SystemExit("CHECK FAILED no HOCOMOCO v11 model for %s" % tf)
    allm = {k: V11.pcc_max(za, z, V11.MIN_OVERLAP)[0] for k, z in z11.items()}
    r_own = sorted(own.values(), reverse=True)
    bo, ba = max(own, key=own.get), max(allm, key=allm.get)
    return dict(n_v11_subtypes=len(own), best_own_v11=bo, best_own_subtype=int(bo.split(".")[-2]), r_best_own=round(own[bo], 3),
                r_second_own=round(r_own[1], 3) if len(r_own) > 1 else np.nan,
                own_subtype_tie=len(r_own) > 1 and r_own[0] - r_own[1] < TIE,
                best_any_v11=ba, best_any_tf=ba.split("_HUMAN")[0], r_best_any=round(allm[ba], 3),
                cofactor=ba.split("_HUMAN")[0] != tf and allm[ba] - own[bo] >= COFACTOR_GAP)


def verdicts(src, rk):
    """One row per candidate TF from the v11 columns; FIG_TFS roles must hold in the data."""
    rows = []
    for tf, x in src.groupby("TF", sort=False):
        clean = x[~x.own_subtype_tie]
        sub = sorted(clean.best_own_subtype.unique())
        cof = x[x.cofactor]
        role = FIG_TFS.get(tf, "")
        rows.append(dict(TF=tf, rank_dauROC=int(rk.set_index("TF").rank_dauROC[tf]), n_v11_subtypes=int(x.n_v11_subtypes.iloc[0]),
                         own_subtypes_nearest=";".join("%d: panels %s" % (k, ",".join(map(str, x.panel[x.best_own_subtype == k])))
                                                       for k in sorted(x.best_own_subtype.unique())),
                         n_own_subtypes_without_tie=len(sub), tie_panels=",".join(map(str, x.panel[x.own_subtype_tie])),
                         cofactor_panels=",".join("%d %s r %.3f (own %.3f)" % (r.panel, r.best_any_v11, r.r_best_any, r.r_best_own)
                                                  for _, r in cof.iterrows()),
                         drawn=tf in FIG_TFS, role={"subtypes": "v11 subtype example", "cofactor": "own motif plus co-factor motif"}.get(role, ""),
                         figure_row=list(FIG_TFS).index(tf) + 1 if tf in FIG_TFS else 0))
        if role == "subtypes" and len(sub) < 2:
            raise SystemExit("CHECK FAILED %s: subtype example with %d own v11 subtype(s) nearest" % (tf, len(sub)))
        if role == "cofactor" and cof.empty:
            raise SystemExit("CHECK FAILED %s: co-factor example without a co-factor panel" % tf)
    v = pd.DataFrame(rows)
    out = os.path.join(HERE, "subtype_verdicts_v11_h3.csv")
    w2.backup(out)
    v.to_csv(out, index=False)
    return v


def ranking():
    t = pd.read_csv(IN["rt"], sep="\t")
    rf = t[(t.Model == "RandomForestClassifier") & (t.PWM == "mono+di")].set_index("TF_name").sort_index()
    mono = t[(t.Model == "Single best mono PWM") & (t.PWM == "mono")].drop_duplicates("TF_name").set_index("TF_name").loc[rf.index]
    for m in ("roc_auc", "pr_auc"):
        if not np.array_equal(rf[m + "_test_H_PWM"].values, mono[m + "_test_H_PWM_mono"].values):
            raise SystemExit("CHECK FAILED %s_test_H_PWM is not the best monoPWM (R2-F2)" % m)
    w2.PD.check_tfs(rf.index, "results_table")
    d = pd.DataFrame({"TF": rf.index, "delta_auROC": (rf.roc_auc_test_H - rf.roc_auc_test_H_PWM).values,
                      "delta_auPRC": (rf.pr_auc_test_H - rf.pr_auc_test_H_PWM).values,
                      "auROC_best_PWM": rf.roc_auc_test_H_PWM.values, "auROC_ArChIPelago": rf.roc_auc_test_H.values,
                      "auPRC_best_PWM": rf.pr_auc_test_H_PWM.values, "auPRC_ArChIPelago": rf.pr_auc_test_H.values})
    d["rank_dauROC"] = d.delta_auROC.rank(ascending=False, method="first").astype(int)
    d["rank_dauPRC"] = d.delta_auPRC.rank(ascending=False, method="first").astype(int)
    for c in ("delta_auROC", "delta_auPRC"):   # a tie at rank N_TOP would make the top-N set depend on the tie order
        s = d[c].sort_values(ascending=False).values
        if s[N_TOP - 1] == s[N_TOP]:
            raise SystemExit("CHECK FAILED tie at rank %d of %s" % (N_TOP, c))
    d["top8_dauROC"] = d.rank_dauROC <= N_TOP
    d["candidate"] = (d.rank_dauROC <= N_TOP) | (d.rank_dauPRC <= N_TOP)
    d = d.sort_values("rank_dauROC")
    out = os.path.join(HERE, "tf_ranking_h3.csv")
    w2.backup(out)
    d.to_csv(out, index=False)
    return d


def importances(tf):
    g = pd.read_csv(os.path.join(IN["gini"], "%s_gini.csv" % tf))
    if not (g.file_matches_header.all() and g.header_tf_is_target.all()):
        raise SystemExit("CHECK FAILED %s: gini file has unmatched PWM headers" % tf)
    if abs(g.gini.sum() - 1) > 1e-6:
        raise SystemExit("CHECK FAILED %s: Gini importances sum to %.6f" % (tf, g.gini.sum()))
    g = g.sort_values("gini", ascending=False).reset_index(drop=True)
    g["rank"] = np.arange(1, len(g) + 1)
    return g.rename(columns={"gini": "importance"})


def main():
    rk = ranking()
    cand = list(rk.TF[rk.candidate])
    if not set(FIG_TFS) <= set(cand):
        raise SystemExit("CHECK FAILED figure TFs not among the candidates: %s" % sorted(set(FIG_TFS) - set(cand)))
    z11 = {k: P.zrows(m) for k, m in V11.read_pcms(IN["v11"]).items()}
    print("top %d by dauROC: %s" % (N_TOP, " ".join(rk.TF[rk.top8_dauROC])))
    print("candidates (%d): %s" % (len(cand), " ".join(cand)))
    t2 = pd.read_excel(IN["t2"], sheet_name="per_TF").set_index("TF")
    meta = P.meta_lookup()
    mdir = os.path.join(HERE, "matrices")
    if os.path.isdir(mdir):
        bdir = mdir + "_backup_20261107"
        if not os.path.exists(bdir):
            shutil.copytree(mdir, bdir)
        shutil.rmtree(mdir)
    os.makedirs(mdir)
    src, check = [], []
    for tf in cand:
        g = importances(tf)
        n_mono = int((g.kind == "mono").sum()); n_di = len(g) - n_mono
        n_files = (len(glob.glob(os.path.join(PWM_ROOT, "PWMs_mono_HUMAN", tf, "*.pwm"))),
                   len(glob.glob(os.path.join(PWM_ROOT, "PWMs_di_HUMAN", tf, "*.dpwm"))))
        if (n_mono, n_di) != n_files or (n_mono, n_di) != (t2.loc[tf, "Human monoPWMs"], t2.loc[tf, "Human diPWMs"]):
            raise SystemExit("CHECK FAILED %s: %d/%d features vs %s PWM files vs Sup. Table 2 %d/%d" % (
                tf, n_mono, n_di, n_files, t2.loc[tf, "Human monoPWMs"], t2.loc[tf, "Human diPWMs"]))
        bad = sum(P.read_pwm(P.pwm_path(tf, f))[0] != h for f, h in zip(g.feature, g.pwm_header))
        if bad:
            raise SystemExit("CHECK FAILED %s: %d PWM headers differ from the gini table" % (tf, bad))
        check.append("%s: %d mono + %d di features = PWM files = Sup. Table 2, 0 header mismatches" % (tf, n_mono, n_di))
        sel = pd.concat([g[g.kind == "mono"].head(N_MONO), g[g.kind == "di"].head(N_DI)])
        row = rk[rk.TF == tf].iloc[0]
        probs = [(P.read_pwm(P.pwm_path(tf, f))) for f in sel.feature]
        probs = [(h, w, P.to_prob(f.split("_")[0], w)) for (h, w), f in zip(probs, sel.feature)]
        ref = P.zrows(probs[0][2])
        for j, ((_, r), (header, w, p)) in enumerate(zip(sel.iterrows(), probs), start=1):
            r_max, strand, _, _ = P.pcc_max(ref, P.zrows(p))
            if strand == "-":
                p = p[::-1, ::-1]
            mfile = os.path.join("matrices", "%s_%s.txt" % (tf, r.feature))
            pd.DataFrame(p.T, index=list(P.ACGT), columns=["pos%d" % (i + 1) for i in range(len(p))]).to_csv(os.path.join(HERE, mfile), sep="\t")
            peaks = re.search(r"PEAKS\d+", header).group(0)
            md = meta.get(peaks, {})
            src.append(dict(TF=tf, panel=j, pwm_class=r.kind, rank_overall=int(r["rank"]), feature=r.feature, pwm_name=header,
                            importance=float(r.importance), n_features=len(g), motif_length=len(p), log_odds_rows=len(w),
                            dataset=peaks, cell_line=md.get("cell_line", ""), treatment=md.get("treatment", ""),
                            length_class=header.split(".")[-1], strand_for_display=strand, r_max_vs_panel1=round(r_max, 3),
                            consensus=P.consensus(p), total_IC_bits=round(float(P.ic_bits(p).sum()), 2),
                            delta_auROC=float(row.delta_auROC), delta_auPRC=float(row.delta_auPRC),
                            rank_dauROC=int(row.rank_dauROC), rank_dauPRC=int(row.rank_dauPRC), matrix_file=mfile,
                            source_matrix=os.path.relpath(P.pwm_path(tf, r.feature), PWM_ROOT), **v11_match(tf, p, z11)))
    src = pd.DataFrame(src)
    src["rank_within_class"] = src.groupby(["TF", "pwm_class"]).cumcount() + 1
    cols = ["TF", "panel", "pwm_class", "rank_within_class", "rank_overall", "feature", "pwm_name", "importance", "n_features",
            "motif_length", "log_odds_rows", "length_class", "dataset", "cell_line", "treatment", "strand_for_display",
            "r_max_vs_panel1", "consensus", "total_IC_bits", "delta_auROC", "delta_auPRC", "rank_dauROC", "rank_dauPRC",
            "matrix_file", "source_matrix", "n_v11_subtypes", "best_own_v11", "best_own_subtype", "r_best_own", "r_second_own",
            "own_subtype_tie", "best_any_v11", "best_any_tf", "r_best_any", "cofactor"]
    out = os.path.join(HERE, "Figure_S7_source_data.csv")
    w2.backup(out)
    src[cols].to_csv(out, index=False)
    pd.set_option("display.width", 250)
    print(src[["TF", "panel", "feature", "importance", "motif_length", "strand_for_display", "r_max_vs_panel1", "consensus"]].to_string(index=False))
    print("\n".join(check))
    v = verdicts(src, rk)
    print(v.to_string(index=False))
    w2.provenance(HERE, "4_visualisation/supplementary_figures/Figure_S7/build_Figure_S7_source_data_h3.py",
                  [out, os.path.join(HERE, "tf_ranking_h3.csv"), os.path.join(HERE, "subtype_verdicts_v11_h3.csv")],
                  [IN["rt"], IN["t2"], IN["meta"], IN["v11"]] + sorted(glob.glob(os.path.join(IN["gini"], "*_gini.csv"))),
                  params="{N_TOP: %d, N_MONO: %d, N_DI: %d, candidates: union of top N_TOP by dauROC and dauPRC, figure: %s, "
                         "cofactor_gap: %.2f, tie: %.2f, v11_min_overlap: %d}" % (N_TOP, N_MONO, N_DI,
                         " ".join("%s %s" % kv for kv in FIG_TFS.items()), COFACTOR_GAP, TIE, V11.MIN_OVERLAP))


if __name__ == "__main__":
    main()
