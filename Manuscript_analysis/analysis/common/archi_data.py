"""
archi_data.py -- reusable loader for the ArChIPelago per-TF feature matrices.

Reads ONLY from the pipeline output directory $ARCHI_RELEASE_DIR (default ~/Release/TF-ML: the files written by
notebooks 0-2); never writes there.
Python 3.8 / sklearn 1.3.0 / pandas.  Re-implements the matrix assembly of notebook 2
(functions group_selection_GC / Scale_transform / model_building / rocauc_plotting) without the shell
`paste` step.

Data layout (per TF, e.g. $ARCHI_RELEASE_DIR/outputdir/SRF/):
  <SEQ>_seq_<PWMSP>_pwm_<mono|di>/<k>_<TF>_<SP>~C[MD]~..._feature_compare_table_<train|control>_cut.tab
      one column: best-hit log-odds of PWM k for every sequence of the GLOBAL table
      $ARCHI_RELEASE_DIR/out_tab_<SEQ>_10000_<train|control>.tab (same row order / row count).
  Global table columns: 0 running index, 1 sequence name ("<id>$>chrom: chrN; center: ...;"),
      2 TF-family code, 3 "<TF>_<SP>".  "train" = train chromosomes, "control" = test chromosomes
      (HUMAN: train = chr2-7, 9, 10, 13-20, control = chr1, chr8, chr21.
       MOUSE: train = chr2-7, 9, 10, 13-18, control = chr1, chr8, chr19.)
      Every caller uses split="train"/"control" in this logical sense.  In the pipeline output directory
      of the manuscript run the MOUSE file named "train" holds chr1/8/19 (170,571 rows) and the file named
      "control" holds chr2-7, 9, 10, 13-18 (345,500 rows); set ARCHI_MOUSE_FILES_SWAPPED=1 for such a
      directory and _file_split() maps the logical split to the other file name.  global_table() checks the
      chromosomes of every table it reads and stops if the split does not match TEST_CHROMS.
  Negative ids (BiasAway GC-matched, other families) per split:
      <mode>_neg_id_out_<SEQ>_result_<H|M>_<train|control>_full.csv  (mode = train|test)

Assembly (exactly as in the notebook):
  positives  = rows with column 3 == "<TF>_<SEQ>", .sample(n=10000, random_state=0) if more;
  negatives  = rows whose id is in the neg-id file,  .sample(n=1_000_000, random_state=0) if more;
  X = PWM columns, y = 1/0;  StandardScaler().fit_transform() applied to EACH matrix separately
  (train, test_H, test_M are each standardised with their own mean/sd -- this is what the
  release code does; use scale="train" for the conventional train-fitted scaler, or None).
  mono_di feature order in the release run: mono PWMs (random.shuffle'd, unseeded), then di PWMs
  sorted by index; here mono and di are both sorted by index (order is irrelevant for the
  RF except through the unseeded RNG, which is not reproducible anyway).
"""
import os
import glob
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, average_precision_score

REL = os.environ.get("ARCHI_RELEASE_DIR", os.path.expanduser("~/Release/TF-ML"))
MOUSE_FILES_SWAPPED = os.environ.get("ARCHI_MOUSE_FILES_SWAPPED", "0") == "1"
OUT = os.path.join(REL, "outputdir")
POS_CLASS = 10000
NEG_CLASS = POS_CLASS * 100
# Test chromosomes in the logical sense used by this loader.
TEST_CHROMS = {"HUMAN": ("chr1", "chr8", "chr21"), "MOUSE": ("chr1", "chr8", "chr19")}


def _file_split(seq_species, split):
    """Logical split ("train" = training chromosomes, "control" = test chromosomes) -> file name
    suffix (inverted for MOUSE when ARCHI_MOUSE_FILES_SWAPPED=1, see module docstring)."""
    if split not in ("train", "control"):
        raise ValueError(split)
    if seq_species == "MOUSE" and MOUSE_FILES_SWAPPED:
        return "control" if split == "train" else "train"
    return split
RF_PARAMS = {"max_depth": 6, "max_samples": 0.8, "n_estimators": 100}   # release notebook


# ----------------------------------------------------------------------------- raw readers
def global_table(seq_species, split):
    """out_tab_<SP>_10000_<train|control>.tab -> DataFrame(idx, name, family, tf); name w/o '$...'."""
    path = os.path.join(REL, "out_tab_%s_10000_%s.tab" % (seq_species, _file_split(seq_species, split)))
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: int, 1: str, 2: str, 3: str})
    df = df.dropna(subset=[1])
    df["chrom"] = df[1].str.extract(r"chrom: (chr\w+);")[0]
    df[1] = df[1].str.split("$").str[0]
    df.columns = ["idx", "name", "family", "tf", "chrom"]
    chroms = set(df.chrom.dropna())
    test = set(TEST_CHROMS[seq_species])
    if (split == "control" and not chroms <= test) or (split == "train" and chroms & test):
        raise ValueError("%s: chromosomes %s do not match the %s split of %s (test = %s); check ARCHI_MOUSE_FILES_SWAPPED"
                         % (path, sorted(chroms), split, seq_species, sorted(test)))
    return df


def pwm_files(tf, seq_species, pwm_species, kind, split):
    """Sorted (by leading integer) list of (k, path) score files for one feature dir."""
    d = os.path.join(OUT, tf, "%s_seq_%s_pwm_%s" % (seq_species, pwm_species, kind))
    if not os.path.isdir(d):
        raise FileNotFoundError(d)
    files = glob.glob(os.path.join(d, "*_feature_compare_table_%s_cut.tab" % _file_split(seq_species, split)))
    out = sorted((int(os.path.basename(f).split("_")[0]), f) for f in files)
    if not out:
        raise FileNotFoundError("no %s score files in %s" % (split, d))
    return out


def feature_list(tf, pwm_species, pwm_set):
    """Ordered feature names 'mono_<k>' / 'di_<k>' available for BOTH human and mouse sequences
    (the release code truncates to the shorter list; for all 36 TFs the lists are equal)."""
    kinds = {"mono": ["mono"], "di": ["di"], "mono_di": ["mono", "di"]}[pwm_set]
    feats = []
    for kind in kinds:
        h = [k for k, _ in pwm_files(tf, "HUMAN", pwm_species, kind, "train")]
        m = [k for k, _ in pwm_files(tf, "MOUSE", pwm_species, kind, "control")]
        if h != m:
            raise ValueError("%s %s %s: human/mouse PWM index lists differ (%d vs %d)"
                             % (tf, pwm_species, kind, len(h), len(m)))
        feats += ["%s_%d" % (kind, k) for k in h]
    return feats


def read_scores(tf, seq_species, pwm_species, split, feature_names):
    """float32 matrix (n_seq_in_global_table x n_features) in the given feature order."""
    paths = {}
    for kind in ("mono", "di"):
        if any(f.startswith(kind + "_") for f in feature_names):
            paths.update({"%s_%d" % (kind, k): p
                          for k, p in pwm_files(tf, seq_species, pwm_species, kind, split)})
    cols = []
    for f in feature_names:
        v = pd.read_csv(paths[f], header=None, sep="\t", usecols=[0], dtype=np.float32)[0].values
        cols.append(v)
    n = {len(c) for c in cols}
    if len(n) != 1:
        raise ValueError("score files of %s/%s/%s have different lengths: %s" % (tf, seq_species, split, n))
    return np.column_stack(cols)


def neg_ids(tf, seq_species, split):
    fsplit = _file_split(seq_species, split)
    mode = "train" if fsplit == "train" else "test"
    flag = "result_%s_%s_full" % (seq_species[0], fsplit)
    path = os.path.join(OUT, tf, "%s_neg_id_out_%s_%s.csv" % (mode, seq_species, flag))
    return pd.read_csv(path, header=None)[0].tolist()


def select_rows(gt, tf, seq_species, split, pos_class=POS_CLASS, neg_class=NEG_CLASS):
    """Replicates group_selection_GC: returns (row positions into the global table, y)."""
    pos = gt[gt["tf"] == "%s_%s" % (tf, seq_species)]
    if pos.shape[0] > pos_class:
        pos = pos.sample(n=pos_class, replace=False, random_state=0)
    neg = gt[gt["name"].isin(set(neg_ids(tf, seq_species, split)))]
    if neg.shape[0] > neg_class:
        neg = neg.sample(n=neg_class, replace=False, random_state=0)
    rows = np.concatenate([pos.index.values, neg.index.values])
    y = np.concatenate([np.ones(len(pos), dtype=np.int8), np.zeros(len(neg), dtype=np.int8)])
    return rows, y


def _scale(X, how, scaler=None):
    if how is None:
        return X, None
    if how == "per_matrix":
        sc = StandardScaler().fit(X)
        return sc.transform(X), sc
    if how == "train":
        if scaler is None:
            scaler = StandardScaler().fit(X)
        return scaler.transform(X), scaler
    raise ValueError(how)


# ----------------------------------------------------------------------------- public API
def load_matrix(tf, seq_species, pwm_species, split, feature_names):
    """One (X, y, meta) block: sequences of `seq_species` on `split` chromosomes scanned with
    `pwm_species` PWMs.  X is unscaled float32.  meta = global-table rows (name, chrom, tf)."""
    gt = global_table(seq_species, split)
    rows, y = select_rows(gt, tf, seq_species, split)
    S = read_scores(tf, seq_species, pwm_species, split, feature_names)
    if S.shape[0] != gt.shape[0]:
        raise ValueError("score rows %d != global table rows %d (%s %s)" % (S.shape[0], gt.shape[0], seq_species, split))
    return S[rows], y, gt.loc[rows, ["name", "chrom", "tf"]].reset_index(drop=True)  # labels == file positions


def load_tf(tf, species_train="HUMAN", pwm_species="HUMAN", pwm_set="mono_di", scale="per_matrix"):
    """
    Returns dict with X_train/y_train (species_train train chromosomes), X_test_H/y_test_H
    (human test chromosomes), X_test_M/y_test_M (mouse test chromosomes), all scanned with
    `pwm_species` PWMs; feature_names, is_mono mask, scalers (dict), meta_* frames.
    scale: "per_matrix" (release behaviour: each matrix standardised on itself),
           "train" (fit on train, applied to tests), or None (raw log-odds).
    species_train="MOUSE", pwm_species="MOUSE" gives the mouse->mouse control.
    """
    feats = feature_list(tf, pwm_species, pwm_set)
    is_mono = np.array([f.startswith("mono_") for f in feats])
    Xtr, ytr, mtr = load_matrix(tf, species_train, pwm_species, "train", feats)
    XH, yH, mH = load_matrix(tf, "HUMAN", pwm_species, "control", feats)
    XM, yM, mM = load_matrix(tf, "MOUSE", pwm_species, "control", feats)
    Xtr_s, sc_tr = _scale(Xtr, scale)
    XH_s, sc_H = _scale(XH, scale, sc_tr)
    XM_s, sc_M = _scale(XM, scale, sc_tr)
    return dict(X_train=Xtr_s, y_train=ytr, X_test_H=XH_s, y_test_H=yH, X_test_M=XM_s, y_test_M=yM,
                feature_names=feats, is_mono=is_mono,
                scalers={"train": sc_tr, "test_H": sc_H, "test_M": sc_M},
                meta_train=mtr, meta_test_H=mH, meta_test_M=mM, tf=tf,
                species_train=species_train, pwm_species=pwm_species, pwm_set=pwm_set, scale=scale)


def auprc_integral(y, s):
    """auPRC with Davis-Goadrich continuous interpolation between consecutive distinct-score
    points == PRROC::pr.curve()$auc.integral, the estimator behind the published tables
    (verified to all printed digits on the SRF single-PWM baselines; see REPORT.md)."""
    y = np.asarray(y, float); s = np.asarray(s, float)
    o = np.argsort(-s, kind="mergesort"); y = y[o]; s = s[o]
    last = np.r_[s[1:] != s[:-1], True]                      # end of each tie block
    tp = np.r_[0.0, np.cumsum(y)[last]]; fp = np.r_[0.0, np.cumsum(1 - y)[last]]
    P = tp[-1]
    a, b = tp[:-1], tp[:-1] + fp[:-1]                        # segment start TP, TP+FP
    c, e = np.diff(tp), np.diff(tp) + np.diff(fp)            # dTP, d(TP+FP)  (e > 0 always)
    with np.errstate(divide="ignore", invalid="ignore"):
        seg = np.where(b == 0, c / e, c / e + (a - b * c / e) / e * np.log((b + e) / np.where(b == 0, 1, b)))
    seg = np.where(c == 0, 0.0, seg)
    return float(np.sum(seg * c) / P)


def auroc_auprc(y, s, prc="integral"):
    """(auROC, auPRC).  prc="integral" (default): PRROC-style integral, PRROC::pr.curve()$auc.integral,
    the estimator of all tables of the manuscript; prc="trapz": sklearn auc(recall, precision). The two
    differ by <1e-3 for continuous RF scores, but for tied single-PWM scores the trapezoid can be far off
    (ERG mono_25: train auPRC 0.163 trapz vs 0.048 integral) and selects different auPRC baseline PWMs."""
    if prc == "trapz":
        p, r, _ = precision_recall_curve(y, s)
        pr = float(auc(r, p))
    elif prc == "integral":
        pr = auprc_integral(y, s)
    else:
        raise ValueError(prc)
    return float(roc_auc_score(y, s)), pr


def average_precision(y, s):
    return float(average_precision_score(y, s))


def best_pwm_baseline(X_train, y_train, X_tests, mask=None):
    """
    Single-PWM baseline as in the release notebook: the column with the highest TRAIN auROC
    (first max wins, as in the notebook loop) is chosen on the human train set and evaluated
    on every test matrix; independently the column with the highest TRAIN auPRC is chosen.
    X_tests: dict name -> (X, y).  mask: boolean column mask (e.g. is_mono) or None = all columns.
    Returns dict(by_roc=dict(col, train_auroc, tests={name: (auroc, auprc)}), by_prc=...).
    Column choice is invariant to the (monotone) standardisation.
    """
    cols = np.arange(X_train.shape[1]) if mask is None else np.flatnonzero(mask)
    tr = np.array([auroc_auprc(y_train, X_train[:, j]) for j in cols])
    out = {}
    for key, metric in (("by_roc", 0), ("by_prc", 1)):
        j = int(cols[int(np.argmax(tr[:, metric]))])
        out[key] = dict(col=j, train_auroc=float(tr[cols == j, 0][0]), train_auprc=float(tr[cols == j, 1][0]),
                        tests={n: auroc_auprc(y, X[:, j]) for n, (X, y) in X_tests.items()})
    return out


def fit_paper_rf(X, y, n_jobs=4, random_state=None):
    """RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100) -- the release
    notebook sets no random_state, so its numbers are not bit-reproducible."""
    m = RandomForestClassifier(n_jobs=n_jobs, random_state=random_state, **RF_PARAMS)
    m.fit(X, y)
    return m


def evaluate(model, d):
    """auROC/auPRC of a fitted model on the three matrices of a load_tf() dict."""
    res = {}
    for k in ("train", "test_H", "test_M"):
        s = model.predict_proba(d["X_" + k])[:, 1]
        res[k] = auroc_auprc(d["y_" + k], s)
    return res


def load_saved_model(tf, pwm_set="mono_di"):
    """finalized_model_<TF>_HUMAN_RandomForestClassifier_<set>_full_MODEL_all_features.sav
    (joblib; 25 features for SRF; column order of the release run = code_H_base_1_*.sh order)."""
    import joblib  # own model files from the published pipeline (trusted source)
    p = os.path.join(OUT, tf, "finalized_model_%s_HUMAN_RandomForestClassifier_%s_full_MODEL_all_features.sav" % (tf, pwm_set))
    return joblib.load(p)


def sh_feature_order(tf, pwm_set="mono_di"):
    """Feature order used when the saved .sav model was fitted, parsed from
    <TF>/code_H_base_1_RandomForestClassifier_<set>.sh (the paste command)."""
    p = os.path.join(OUT, tf, "code_H_base_1_RandomForestClassifier_%s.sh" % pwm_set)
    txt = open(p).read()
    feats = []
    for tok in txt.split("<( cut -f 1 ./")[1:]:
        d, f = tok.split(" ")[0].split("/")
        kind = "mono" if d.endswith("_mono") else "di"
        feats.append("%s_%s" % (kind, f.split("_")[0]))
    return feats


def reorder(d, feature_names):
    """Return copies of the X matrices of a load_tf() dict in a new feature order."""
    idx = [d["feature_names"].index(f) for f in feature_names]
    out = dict(d)
    for k in ("X_train", "X_test_H", "X_test_M"):
        out[k] = d[k][:, idx]
    out["feature_names"] = list(feature_names)
    out["is_mono"] = d["is_mono"][idx]
    return out
