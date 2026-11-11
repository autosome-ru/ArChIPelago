"""Reads the old (first-submission, 36-TF) values from the frozen copies results/previous_values/<git ref>/<file name>, which hold
the files as they were at that git ref before any rebuild. Live tables, which a rebuild rewrites, are never used. Nothing is
written to disk.
usage: import h3_old as O; O.csv("3_analysis/low_fpr_metrics/B2_summary.csv"); O.src(path) -> label for 'source file'
"""
import io
import os
import subprocess

import pandas as pd

TAG = "pre-mouse-sym-20261106"
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "previous_values"))


def _blob(path):
    f = os.path.join(ROOT, TAG, os.path.basename(path))
    if not os.path.isfile(f):
        raise SystemExit("CHECK FAILED old value: %s (git tag %s) has no frozen copy %s" % (path, TAG, f))
    with open(f, "rb") as fh:
        return fh.read()


def csv(path, **kw):
    return pd.read_csv(io.BytesIO(_blob(path)), **kw)


def excel(path, **kw):
    return pd.read_excel(io.BytesIO(_blob(path)), **kw)


def src(path):
    """'source file' label of an old value: 'old git <tag>:<path>' (the test tracer skips it; it holds no new value)."""
    return "old git %s:%s" % (TAG, path)
