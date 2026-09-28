#!/usr/bin/env python3
"""Similarity of the baseline human monoPWM of every TF to the mouse monoPWMs of the same TF (Fig. S6C, Sup. Table 5).

Baseline human monoPWM = the best monoPWM by training auROC (column HM_broc_col of ../mouse_transfer/mouse_transfer_results.csv),
mapped to its HOCOMOCO name through ../inputs/single_PWM_features.csv (feature index -> PWM file header).
Similarity of two PWMs = maximum over strands and offsets (overlap >= 5 columns) of the mean per-column Pearson
correlation of the log-odds matrices. For every TF:
  topH_PWM                  name of the baseline human monoPWM
  topH_vs_nearestM_pcc      its similarity to the most similar mouse monoPWM of the TF
  median_H_vs_nearestM_pcc  median over the human monoPWMs of the similarity to their most similar mouse monoPWM
  median_H_vs_nearestH_pcc  the same within the human monoPWMs (self excluded)
  n_monoPWM_H_found, n_monoPWM_M_found   number of PWM files read

PWM files: $ARCHI_ZENODO_DIR/PWMs_mono_HUMAN/<TF>/*.pwm and PWMs_mono_MOUSE.tar.gz of the Zenodo archive.
Output: motif_similarity.csv in this folder (read by cross_species_table.py).
"""
import glob
import os
import sys
import tarfile
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
from archi_paths import BASIS, INPUTS, ZENODO  # noqa: E402


def read_pwm(path):
    return np.loadtxt(path, skiprows=1).reshape(-1, 4)


def zrows(m):
    """z-score each PWM column (position) across the 4 nucleotides"""
    m = m - m.mean(axis=1, keepdims=True)
    sd = m.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return m / sd


def pcc_max(a, b, min_overlap=5):
    """max over offsets and strands of the mean per-column Pearson r of two z-scored log-odds PWMs"""
    best = -1.0
    for bb in (b, b[::-1, ::-1]):
        R = a @ bb.T / 4.0
        la, lb = R.shape
        for k in range(-(lb - min_overlap), la - min_overlap + 1):
            d = np.diagonal(R, offset=-k)
            if len(d) >= min_overlap:
                best = max(best, d.mean())
    return best


def load_dir(d):
    files = glob.glob(os.path.join(d, "*.pwm"))
    if not files:
        raise SystemExit("no PWM files in %s" % d)
    return {open(f).readline().strip().lstrip(">"): zrows(read_pwm(f)) for f in files}


hm = pd.read_csv(os.path.join(BASIS, "mouse_transfer", "mouse_transfer_results.csv")).set_index("TF").sort_index()
names = pd.read_csv(os.path.join(INPUTS, "single_PWM_features.csv")).set_index(["TF", "feature"]).pwm_name
TFS = list(hm.index)
assert len(TFS) == 36

with tempfile.TemporaryDirectory() as tmp:
    with tarfile.open(os.path.join(ZENODO, "PWMs_mono_MOUSE.tar.gz")) as tf_:
        tf_.extractall(tmp)
    mouse_root = os.path.join(tmp, "PWMs_MOUSE_mono")
    E = pd.DataFrame(index=pd.Index(TFS, name="TF"))
    for tf in TFS:
        H = load_dir(os.path.join(ZENODO, "PWMs_mono_HUMAN", tf))
        M = load_dir(os.path.join(mouse_root, tf))
        hn, mn = list(H), list(M)
        S = np.array([[pcc_max(H[h], M[m]) for m in mn] for h in hn])
        nn_HM = S.max(axis=1)
        if len(hn) > 1:
            SH = np.array([[pcc_max(H[a], H[b]) if a != b else -1 for b in hn] for a in hn])
            nn_HH = SH.max(axis=1)
        else:
            nn_HH = np.array([np.nan])
        base = names.loc[(tf, hm.loc[tf, "HM_broc_col"])]
        if base not in hn:
            raise SystemExit("%s: baseline monoPWM %s not among the PWM files" % (tf, base))
        E.loc[tf, "n_monoPWM_H_found"] = len(hn)
        E.loc[tf, "n_monoPWM_M_found"] = len(mn)
        E.loc[tf, "topH_PWM"] = base
        E.loc[tf, "topH_vs_nearestM_pcc"] = nn_HM[hn.index(base)]
        E.loc[tf, "median_H_vs_nearestM_pcc"] = np.median(nn_HM)
        E.loc[tf, "median_H_vs_nearestH_pcc"] = np.nanmedian(nn_HH)
        print(tf, len(hn), len(mn), round(E.loc[tf, "topH_vs_nearestM_pcc"], 4), flush=True)

for c in ("n_monoPWM_H_found", "n_monoPWM_M_found"):
    E[c] = E[c].astype(int)
E.to_csv(os.path.join(HERE, "motif_similarity.csv"), float_format="%.4f")
print("written:", os.path.join(HERE, "motif_similarity.csv"))
