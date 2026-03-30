"""SARUS PWM scanning and feature matrix construction.

Wraps the SARUS Java tool (command line unchanged) and builds feature
matrices from its output for use in model training and prediction.
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# SARUS scanning
# ---------------------------------------------------------------------------

def run_sarus(fasta_path, pwm_path, sarus_jar, output_path, java_bin="java",
              pwm_type="mono", xmx="2G"):
    """Run SARUS PWM scanner on a FASTA file.

    The shell command is kept exactly as used in the original notebooks.
    For mono PWMs uses `ru.autosome.SARUS`; for di-nucleotide PWMs uses
    `ru.autosome.di.SARUS`.

    Parameters
    ----------
    fasta_path : str or Path
        Input FASTA file.
    pwm_path : str or Path
        PWM file (.pwm for mono, .dpwm for di).
    sarus_jar : str or Path
        Path to sarus.jar.
    output_path : str or Path
        Where to write the SARUS score output.
    java_bin : str
        Java executable (default: 'java').
    pwm_type : str
        'mono' or 'di'.
    xmx : str
        Java heap size (default '2G').

    Returns
    -------
    Path
        Path to the SARUS output file.
    """
    fasta_path = Path(fasta_path)
    pwm_path = Path(pwm_path)
    sarus_jar = Path(sarus_jar)
    output_path = Path(output_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    if not pwm_path.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_path}")
    if not sarus_jar.exists():
        raise FileNotFoundError(
            f"SARUS jar not found: {sarus_jar}\n"
            "Download SARUS from https://github.com/autosome-ru/sarus"
        )

    sarus_class = "ru.autosome.SARUS" if pwm_type == "mono" else "ru.autosome.di.SARUS"

    # Command kept identical to notebooks to preserve reproducibility
    cmd = (
        f"{java_bin} -Xmx{xmx} -cp {sarus_jar} {sarus_class} "
        f"{fasta_path} {pwm_path} besthit --output-scoring-mode score"
    )

    with open(output_path, "w") as out_f:
        result = subprocess.run(
            cmd,
            shell=True,
            stdout=out_f,
            stderr=subprocess.PIPE,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"SARUS failed (exit {result.returncode}) for {pwm_path.name}:\n"
            f"{result.stderr.decode(errors='replace')}"
        )

    return output_path


def load_sarus_scores(score_file):
    """Load a SARUS output file as a pandas Series of scores.

    Parameters
    ----------
    score_file : str or Path

    Returns
    -------
    pd.Series
    """
    path = Path(score_file)
    if not path.exists() or path.stat().st_size == 0:
        return pd.Series(dtype=float)
    scores = pd.read_csv(path, header=None, sep="\t", comment=">")[0]
    return scores.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Feature matrix construction
# ---------------------------------------------------------------------------

def build_feature_matrix(scan_dir, mode="mono_di"):
    """Build a feature matrix from a directory of SARUS score files.

    Expects per-PWM score files named ``<pwm_name>.txt`` (or ``.tsv``) inside
    ``<scan_dir>/mono/`` and/or ``<scan_dir>/di/`` subdirectories depending on
    ``mode``.

    Parameters
    ----------
    scan_dir : str or Path
        Directory produced by scanning.  Must contain ``mono/`` and/or ``di/``
        subdirectories.
    mode : str
        One of 'mono', 'di', 'mono_di'.

    Returns
    -------
    pd.DataFrame
        Rows = sequences, columns = PWM feature names.
    """
    scan_dir = Path(scan_dir)
    dfs = []

    if mode in ("mono", "mono_di"):
        mono_dir = scan_dir / "mono"
        if not mono_dir.exists():
            raise FileNotFoundError(f"mono scan directory not found: {mono_dir}")
        for f in sorted(mono_dir.glob("*.txt")):
            scores = load_sarus_scores(f)
            dfs.append(scores.rename(f"mono_{f.stem}"))

    if mode in ("di", "mono_di"):
        di_dir = scan_dir / "di"
        if not di_dir.exists():
            raise FileNotFoundError(f"di scan directory not found: {di_dir}")
        for f in sorted(di_dir.glob("*.txt")):
            scores = load_sarus_scores(f)
            dfs.append(scores.rename(f"di_{f.stem}"))

    if not dfs:
        raise ValueError(
            f"No score files found in {scan_dir} for mode='{mode}'"
        )

    return pd.concat(dfs, axis=1).fillna(0.0)


def select_top_features(X, y, n=1000):
    """Select the top N features by Random Forest feature importance.

    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix.
    y : array-like
        Binary labels (0/1).
    n : int
        Number of features to select.

    Returns
    -------
    list[str]
        Selected column names (subset of X.columns).
    """
    from sklearn.ensemble import RandomForestClassifier

    n = min(n, X.shape[1])
    selector = RandomForestClassifier(
        n_estimators=100, max_depth=6, n_jobs=4, random_state=42
    )
    selector.fit(X.fillna(0), y)
    importances = pd.Series(selector.feature_importances_, index=X.columns)
    return importances.nlargest(n).index.tolist()
