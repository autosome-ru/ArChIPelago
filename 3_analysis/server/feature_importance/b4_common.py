"""Shared helpers for the feature-importance analysis (b4_01 to b4_07).
Needs py3.8, sklearn 1.3.0 (2_training/VERSIONS.txt). Reads $ARCHI_RELEASE only.
"""
import os, sys, re, hashlib, subprocess
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..")); import archi_env as AE
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "..", "1_data_preparation"))  # human data build
import archi_data_v5 as ad  # loader with negatives v3

OUT = AE.work("H5_alltrain", "human_sym", "B4_importance")
MODELS = AE.work("H5_alltrain", "human_sym", "models")
CACHE = os.path.join(OUT, "cache")
TFS = ["ANDR", "AP2A", "CEBPB", "COE1", "CTCF", "E2F4", "ERG", "ESR1", "FLI1", "GATA1", "GATA2", "GATA3", "GCR",
       "HNF4A", "IRF1", "IRF4", "JUND", "MAFK", "MAX", "MYC", "PPARG", "PRGR", "REST", "RUNX1", "RXRA",
       "SOX2", "SPI1", "SRF", "STA5A", "STAT1", "STAT3", "TAL1", "TF65", "TFE2", "USF2"]
NEG_CAP = 100000   # negatives per test set; never triggered (largest human test set has 89,333 rows)


def model_path(tf):
    return os.path.join(MODELS, "%s_HUMAN_RandomForestClassifier_mono_di_seed0_human_sym_negv3.sav" % tf)


def feats_in_model_order(tf, model):
    import json  # the feature order is saved with the model (run_hm_mm_h5.py)
    f = json.load(open(model_path(tf)[:-4] + "_features.json"))["feature_names"]
    if len(f) != model.n_features_in_:
        raise ValueError("%s: .sh lists %d features, model has %d" % (tf, len(f), model.n_features_in_))
    return f


def pwm_dir(tf, kind):
    return os.path.join(ad.OUT, tf, "HUMAN_seq_HUMAN_pwm_" + kind)


def read_pwm_file(tf, feat):
    """(header, matrix) of the PWM file of a feature 'mono_k'/'di_k'. mono: Lx4 log-odds; di: Lx16."""
    kind, k = feat.split("_")
    p = os.path.join(pwm_dir(tf, kind), "%s.%s" % (k, "pwm" if kind == "mono" else "dpwm"))
    with open(p) as fh:
        header = fh.readline().strip().lstrip(">")
    m = np.loadtxt(p, skiprows=1)
    return header, m.reshape(-1, 4 if kind == "mono" else 16)


def feature_file_name(tf, feat):
    kind, k = feat.split("_")
    for _, p in ad.pwm_files(tf, "HUMAN", "HUMAN", kind, "control"):
        if int(os.path.basename(p).split("_")[0]) == int(k):
            return os.path.basename(p)
    raise FileNotFoundError(feat)


def load_test(tf, feats):
    """Raw (unscaled) human test matrix X (float32), y, names. Test = chr1/8/21 as in the paper."""
    X, y, meta = ad.load_matrix(tf, "HUMAN", "HUMAN", "control", feats)
    return X, y, meta


_TRAIN_SCALERS = {}


def train_scaler(tf, feats):
    """StandardScaler fitted on the human TRAINING matrix (train chromosomes, loader row selection) in the
        feature order `feats`; the same scaler as archi_data_v5.load_tf(..., scale="train")["scalers"]["train"]."""
    from sklearn.preprocessing import StandardScaler
    key = (tf, tuple(feats))
    if key not in _TRAIN_SCALERS:
        Xtr, _, _ = ad.load_matrix(tf, "HUMAN", "HUMAN", "train", list(feats))
        _TRAIN_SCALERS[key] = StandardScaler().fit(Xtr)
    return _TRAIN_SCALERS[key]


def scale_test(X, tf, feats):
    """the test matrix is standardised with the TRAINING scaler."""
    return train_scaler(tf, feats).transform(X).astype(np.float32)


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()
