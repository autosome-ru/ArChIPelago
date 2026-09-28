#!/usr/bin/env python3
"""Motif-subtype supplementary figure (Fig. S7): source data.

1. Ranks the 36 TFs by the gain of the mono+di ArChIPelago RF over the best monoPWM on the human
   test set (Figure_3_source_data.csv, rows PWM_set == "mono+di", test_set == "human"); candidates =
   union of the top 8 by dauROC and by dauPRC.
2. RF feature importances of the all-feature models
   `outputdir/<TF>/finalized_model_<TF>_HUMAN_RandomForestClassifier_mono_di_full_MODEL_all_features.sav`
   (`feature_importances_` in the column order of `code_H_base_1_RandomForestClassifier_mono_di.sh`), saved
   as importances/sav_model_importances_raw.tsv.
3. Feature -> PWM file: feature `mono_<k>` / `di_<k>` -> PWMs_{mono,di}_HUMAN/<TF>/<k>.{pwm,dpwm}; the
   header name of the file (dots) must equal the dataset name in the feature file name (underscores).
4. Log-odds -> probabilities per position: p_i(x) = 0.25*exp(w_i(x)) / sum (mono; uniform background);
   diPWM: P_i(xy) = exp(w_i(xy))/16 / sum over the 16 dinucleotides, then the mononucleotide marginal:
   position i (i = 1..L) = first-nucleotide marginal of row i, position L+1 = second-nucleotide marginal
   of row L (an L-row diPWM scores L+1 nt) -- the same projection HOCOMOCO uses to draw the mono-style
   logo of a dinucleotide PCM (summing the dinucleotide counts), here applied to the uniform-background
   probabilities. Consensus per position: the letter if p >= 0.6, the two-letter IUPAC code if the top two
   letters sum to >= 0.8, else N (lower case = 0.6 > p >= 0.4 single letter).

Outputs (all in this folder): tf_ranking.csv, importances/<TF>_HUMAN_forest_importances_all_features.csv,
Figure_S7_source_data.csv, matrices/<TF>_<feature>.txt (probability matrices, 4 x L, for R).
r_max_vs_panel1 of the source data = similarity of each panel to panel 1 = the top monoPWM (max over offsets/strands of the mean per-column Pearson r of the
per-position z-scored probability matrices, min overlap 5, as in ../crossspecies/motif_similarity.py; panels on the '-'
strand are reverse-complemented for display so that all logos of a TF are on the strand of panel 1).  Never writes outside this folder.
"""
import csv
import glob
import math
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
from archi_paths import RES, INPUTS, ZENODO as PWM_ROOT  # noqa: E402
META = os.path.join(INPUTS, "metadata_new_23_02_17_with_control.csv")
SRC = os.path.join(RES, "Figures", "source_data", "Figure_3_source_data.csv")
SUP5 = os.path.join(RES, "Sup_Tables", "Sup_Table_5_cross_species_and_mouse_trained.csv")
RAW = os.path.join(HERE, "importances", "sav_model_importances_raw.tsv")
N_TOP = 8
N_MONO, N_DI = 4, 2
IUPAC = {frozenset("AC"): "M", frozenset("AG"): "R", frozenset("AT"): "W", frozenset("CG"): "S",
         frozenset("CT"): "Y", frozenset("GT"): "K"}
ACGT = "ACGT"
DINUC = [a + b for a in ACGT for b in ACGT]   # row order of HOCOMOCO dinucleotide matrices: AA AC AG AT CA ...


# ---------------------------------------------------------------- 1. ranking
def ranking():
    d = pd.read_csv(SRC)
    d = d[(d.PWM_set == "mono+di") & (d.test_set == "human")].copy()
    if len(d) != 36:
        raise SystemExit("expected 36 TFs, got %d" % len(d))
    d["rank_dauROC"] = d.delta_auROC.rank(ascending=False).astype(int)
    d["rank_dauPRC"] = d.delta_auPRC.rank(ascending=False).astype(int)
    d["candidate"] = (d.rank_dauROC <= N_TOP) | (d.rank_dauPRC <= N_TOP)
    d = d.sort_values("rank_dauROC")
    d[["TF", "delta_auROC", "rank_dauROC", "delta_auPRC", "rank_dauPRC", "auROC_best_PWM", "auROC_ArChIPelago",
       "auPRC_best_PWM", "auPRC_ArChIPelago", "n_PWMs_total_mono_plus_di", "candidate"]].to_csv(
        os.path.join(HERE, "tf_ranking.csv"), index=False)
    return d


# ---------------------------------------------------------------- 2. importances
def pwm_path(tf, feat):
    kind, k = feat.split("_")
    return os.path.join(PWM_ROOT, "PWMs_%s_HUMAN" % kind, tf, "%s.%s" % (k, "pwm" if kind == "mono" else "dpwm"))


def read_pwm(path):
    with open(path) as fh:
        header = fh.readline().strip().lstrip(">")
        w = np.array([[float(v) for v in line.split()] for line in fh if line.strip()])
    return header, w


def all_feature_importances():
    raw = pd.read_csv(RAW, sep="\t", names=["TF", "feature", "feature_file", "importance"])
    out = {}
    for tf, g in raw.groupby("TF"):
        g = g.sort_values("importance", ascending=False).reset_index(drop=True)
        g["rank"] = np.arange(1, len(g) + 1)
        g["dataset_name"] = g.feature_file.str.replace("_feature_compare_table_control_cut.tab", "", regex=False)
        g["dataset_name"] = g.dataset_name.str.split("_").str[1:].str.join("_")   # drop the leading index
        g[["rank", "feature", "dataset_name", "importance"]].to_csv(
            os.path.join(HERE, "importances", "%s_HUMAN_forest_importances_all_features.csv" % tf), index=False)
        out[tf] = g
    return out


# ---------------------------------------------------------------- 3./4. matrices
def to_prob(kind, w):
    if kind == "mono":
        if w.shape[1] != 4:
            raise SystemExit("mono matrix with %d columns" % w.shape[1])
        p = 0.25 * np.exp(w)
        return p / p.sum(axis=1, keepdims=True)
    if w.shape[1] != 16:
        raise SystemExit("di matrix with %d columns" % w.shape[1])
    P = np.exp(w) / 16.0
    P = P / P.sum(axis=1, keepdims=True)
    P = P.reshape(len(w), 4, 4)                      # [row, first nt, second nt]
    first = P.sum(axis=2)                            # rows 1..L -> position 1..L
    last = P[-1].sum(axis=0)                         # second nt of the last row -> position L+1
    return np.vstack([first, last[None, :]])


def consensus(p):
    s = []
    for row in p:
        order = np.argsort(row)[::-1]
        if row[order[0]] >= 0.6:
            s.append(ACGT[order[0]])
        elif row[order[0]] + row[order[1]] >= 0.8:
            s.append(IUPAC[frozenset(ACGT[order[0]] + ACGT[order[1]])])
        elif row[order[0]] >= 0.4:
            s.append(ACGT[order[0]].lower())
        else:
            s.append("N")
    return "".join(s)


def ic_bits(p):
    return (2 + (p * np.log2(np.clip(p, 1e-12, None))).sum(axis=1))


def zrows(m):
    m = m - m.mean(axis=1, keepdims=True)
    sd = m.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return m / sd


def pcc_max(a, b, min_overlap=5):
    """As ../crossspecies/motif_similarity.py: max over offsets and both strands of the mean
    per-column Pearson r of two per-position z-scored 4-column matrices. Returns (r, strand, offset, overlap):
    strand "-" = b reverse-complemented; offset k = a[i] is aligned with b[i-k]."""
    best = (-1.0, "+", 0, 0)
    for strand, bb in (("+", b), ("-", b[::-1, ::-1])):
        R = a @ bb.T / 4.0
        la, lb = R.shape
        for k in range(-(lb - min_overlap), la - min_overlap + 1):
            d = np.diagonal(R, offset=-k)
            if len(d) >= min_overlap and d.mean() > best[0]:
                best = (float(d.mean()), strand, k, len(d))
    return best


def meta_lookup():
    m = {}
    with open(META, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) > 6 and row[6].startswith("PEAKS"):
                m[row[6]] = dict(cell_line=row[5], treatment=row[4], antibody=row[3], exp=row[0])
    return m


def main():
    rk = ranking()
    cand = list(rk.TF[rk.candidate])
    print("candidates (%d): %s" % (len(cand), " ".join(cand)))
    imp = all_feature_importances()
    missing = [tf for tf in cand if tf not in imp]
    if missing:
        raise SystemExit("no .sav importances for %s" % missing)
    meta = meta_lookup()
    sup5 = pd.read_csv(SUP5).set_index("TF")
    os.makedirs(os.path.join(HERE, "matrices"), exist_ok=True)

    src, check = [], []
    for tf in cand:
        g = imp[tf]
        # mapping check for EVERY feature of the TF, not only the plotted ones
        n_mono = (g.feature.str.startswith("mono")).sum(); n_di = len(g) - n_mono
        n_files_mono = len(glob.glob(os.path.join(PWM_ROOT, "PWMs_mono_HUMAN", tf, "*.pwm")))
        n_files_di = len(glob.glob(os.path.join(PWM_ROOT, "PWMs_di_HUMAN", tf, "*.dpwm")))
        if (n_mono, n_di) != (n_files_mono, n_files_di):
            raise SystemExit("%s: %d/%d features vs %d/%d PWM files" % (tf, n_mono, n_di, n_files_mono, n_files_di))
        if (n_mono, n_di) != (sup5.loc[tf, "Human monoPWMs"], sup5.loc[tf, "Human diPWMs"]):
            raise SystemExit("%s: feature counts differ from Sup. Table 5" % tf)
        bad = 0
        for feat, ds in zip(g.feature, g.dataset_name):
            header, _ = read_pwm(pwm_path(tf, feat))
            if header.replace(".", "_") != ds:
                bad += 1
                check.append("%s %s: header %s != feature %s" % (tf, feat, header, ds))
        check.append("%s: %d mono + %d di features, %d PWM files each class, %d header mismatches"
                     % (tf, n_mono, n_di, n_files_mono + n_files_di, bad))
        if bad:
            raise SystemExit("name mapping failed for %s" % tf)

        sel = pd.concat([g[g.feature.str.startswith("mono")].head(N_MONO), g[g.feature.str.startswith("di")].head(N_DI)])
        row = rk[rk.TF == tf].iloc[0]
        probs = []
        for _, r in sel.iterrows():
            header, w = read_pwm(pwm_path(tf, r.feature))
            probs.append((header, w, to_prob(r.feature.split("_")[0], w)))
        ref = zrows(probs[0][2])
        for j, (_, r) in enumerate(sel.iterrows(), start=1):
            kind = r.feature.split("_")[0]
            header, w, p = probs[j - 1]
            r_max, strand, offset, overlap = pcc_max(ref, zrows(p))
            if strand == "-":
                p = p[::-1, ::-1]                    # reverse complement for display (same strand as panel 1)
            mfile = os.path.join("matrices", "%s_%s.txt" % (tf, r.feature))
            pd.DataFrame(p.T, index=list(ACGT), columns=["pos%d" % (i + 1) for i in range(len(p))]).to_csv(
                os.path.join(HERE, mfile), sep="\t")
            peaks = re.search(r"PEAKS\d+", header).group(0)
            md = meta.get(peaks, {})
            src.append(dict(TF=tf, panel=j, pwm_class=kind, rank_overall=int(r["rank"]), feature=r.feature, pwm_name=header,
                            importance=float(r.importance), n_features=len(g),
                            motif_length=len(p), log_odds_rows=len(w),
                            dataset=peaks, cell_line=md.get("cell_line", ""), treatment=md.get("treatment", ""),
                            length_class=header.split(".")[-1],
                            strand_for_display=strand, r_max_vs_panel1=round(r_max, 3),
                            consensus=consensus(p), total_IC_bits=round(float(ic_bits(p).sum()), 2),
                            delta_auROC=float(row.delta_auROC), delta_auPRC=float(row.delta_auPRC),
                            rank_dauROC=int(row.rank_dauROC), rank_dauPRC=int(row.rank_dauPRC),
                            matrix_file=mfile,
                            source_matrix=os.path.relpath(pwm_path(tf, r.feature), os.path.expanduser("~"))))
    src = pd.DataFrame(src)
    src["rank_within_class"] = src.groupby(["TF", "pwm_class"]).cumcount() + 1
    src = src[["TF", "panel", "pwm_class", "rank_within_class", "rank_overall", "feature", "pwm_name", "importance",
               "n_features", "motif_length", "log_odds_rows", "length_class", "dataset", "cell_line", "treatment",
               "strand_for_display", "r_max_vs_panel1", "consensus", "total_IC_bits", "delta_auROC", "delta_auPRC", "rank_dauROC", "rank_dauPRC",
               "matrix_file", "source_matrix"]]
    src.to_csv(os.path.join(HERE, "Figure_S7_source_data.csv"), index=False)
    pd.set_option("display.width", 250)
    print(src[["TF", "panel", "feature", "pwm_name", "importance", "motif_length", "strand_for_display", "r_max_vs_panel1", "consensus"]]
          .to_string(index=False))
    print("\n".join(check))


if __name__ == "__main__":
    main()
