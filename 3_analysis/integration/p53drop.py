"""Shared helpers for the P53 drop (human analysis on 35 TFs).

EXCLUDE_TFS comes from the environment variable ARCHI_EXCLUDE_TFS (comma separated, default "P53"); setting it to ""
keeps P53 and makes check_tfs() fail, which is the fail-once test of every patched producer.
The reference TF set is the TF column of results/datasets/peak_ids/mouse_sym_peak_ids.csv (35 TFs, made by the dataset rule),
so the check does not depend on EXCLUDE_TFS itself.
"""
import datetime
import hashlib
import os
import shutil
import sys

import pandas as pd

R2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B12 = os.path.join(R2, "3_analysis/integration")
HUMAN_IDS = os.path.join(R2, "results/datasets/peak_ids", "human_used_peak_ids.csv")
MOUSE_IDS = os.path.join(R2, "results/datasets/peak_ids", "mouse_sym_peak_ids.csv")
BACKUP_TAG = "backup_20261107"   # tag of the backups made by backup()
OLD_TAG = "backup_20261106"      # tag of the earlier backups, from which old values are read
EXCLUDE_TFS = tuple(t for t in os.environ.get("ARCHI_EXCLUDE_TFS", "P53").split(",") if t)


def reference_tfs():
    return sorted(set(pd.read_csv(MOUSE_IDS).TF))


def human_ids():
    """The 361 human peak sets used (human_used_peak_ids.csv without EXCLUDE_TFS)."""
    h = pd.read_csv(HUMAN_IDS)
    return h[~h.TF.isin(EXCLUDE_TFS)].reset_index(drop=True)


def mouse_ids():
    return pd.read_csv(MOUSE_IDS)


def drop(df, col="TF"):
    return df[~df[col].isin(EXCLUDE_TFS)].copy()


def check_tfs(tfs, label):
    """Stop unless the TF set equals the 35 reference TFs (P53 absent)."""
    got = sorted(set(tfs))
    ref = reference_tfs()
    if got != ref:
        extra, miss = sorted(set(got) - set(ref)), sorted(set(ref) - set(got))
        raise SystemExit("CHECK FAILED %s: %d TFs, extra %s, missing %s" % (label, len(got), extra, miss))
    print("check %s: %d TFs, P53 absent" % (label, len(got)))
    return len(got)


def backup(path, tier=""):
    """Copy path to <stem>_<BACKUP_TAG><tier><ext> once: an existing backup is never overwritten, so the first copy holds the
    state before the first run. A later run can pass its own tier to keep the state before that run."""
    if not os.path.exists(path):
        return None
    stem, ext = os.path.splitext(path)
    dst = "%s_%s%s%s" % (stem, BACKUP_TAG, tier, ext)
    if not os.path.exists(dst):
        shutil.copy2(path, dst)
        print("backup", os.path.relpath(dst, R2))
    return dst


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_provenance(folder, script, outputs, inputs, params=None, note=""):
    """Append one block to <folder>/_provenance_p53drop.yaml (outputs and inputs with md5)."""
    p = os.path.join(folder, "_provenance_p53drop.yaml")
    L = ["- script: %s" % script,
         "  run: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
         "  task: W1-L P53 drop (human 35 TFs)",
         "  python: %s" % sys.version.split()[0],
         "  exclude_tfs: [%s]" % ", ".join(EXCLUDE_TFS),
         "  reference_tf_list: {file: %s, md5: %s}" % (os.path.relpath(MOUSE_IDS, R2), md5(MOUSE_IDS)),
         "  human_used_ids: {file: %s, md5: %s}" % (os.path.relpath(HUMAN_IDS, R2), md5(HUMAN_IDS))]
    if params:
        L.append("  parameters: %s" % params)
    if note:
        L.append("  note: \"%s\"" % note.replace('"', "'"))
    L.append("  inputs:")
    L += ["    - {file: %s, md5: %s}" % (os.path.relpath(f, R2), md5(f)) for f in inputs]
    L.append("  outputs:")
    L += ["    - {file: %s, md5: %s}" % (os.path.relpath(f, R2), md5(f)) for f in outputs]
    with open(p, "a") as fh:
        fh.write("\n".join(L) + "\n")
    print("provenance ->", p)
