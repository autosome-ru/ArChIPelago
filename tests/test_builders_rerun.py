"""Local builders rerun in a scratch copy of the repository: each must reproduce the shipped table from the current inputs.

The copy (tmp dir) holds every folder the builders read; results/round1_package/Resubmission/Sup_Tables (the first-submission
tables) is copied to the matching place beside it, and results/round1_package/Manuscript_analysis is only read. The
builders write into the copy. A snapshot of the real tree (size and mtime of every file) taken before and after
proves that nothing in it was written.
A shipped output that the rerun does not reproduce is stale (an input changed after the table was written) or was edited
by hand.
"""
import hashlib
import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
SKIP_DIRS = {"tests", "C_manuscript", "E_figures", "00_reviews", "B1_literature"}
ROUND1 = os.path.join("results", "round1_package", "Resubmission", "Sup_Tables")
# (script relative to the repository root, args, outputs relative to the root), in dependency order
BUILDERS = [
    ("3_analysis/integration/top1000_check_h3.py", [], ["3_analysis/integration/top1000_check_h3.csv"]),
    ("3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py", ["1", "2", "3", "4", "6"], ["outputs/supplementary_tables/Sup_Table_1_ChIP-Seq_experiments.xlsx",
     "outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx", "outputs/supplementary_tables/Sup_Table_3_ArChIPelago_performance.xlsx",
     "outputs/supplementary_tables/Sup_Table_4_operational_metrics.xlsx", "outputs/supplementary_tables/Sup_Table_6_runtime_memory.xlsx"]),
    ("3_analysis/supplementary_tables/make_sup_table_5_R2.py", [], ["outputs/supplementary_tables/Sup_Table_5_cross_species_and_mouse_trained.xlsx"]),
    ("3_analysis/supplementary_tables/make_sup_tables_7_10.py", ["7", "8", "9", "10"], ["outputs/supplementary_tables/Sup_Table_7_low_FPR_metrics.xlsx",
     "outputs/supplementary_tables/Sup_Table_8_top_peaks_and_peak_widths.xlsx", "outputs/supplementary_tables/Sup_Table_9_feature_importance.xlsx",
     "outputs/supplementary_tables/Sup_Table_10_transfer_drivers_and_headroom.xlsx"]),
    ("3_analysis/integration/figure1_values.py", [], ["3_analysis/integration/figure1_values.csv"]),
    ("3_analysis/integration/table1_inputs.py", [], ["3_analysis/integration/table1_rows.csv"]),
    ("3_analysis/integration/pwm_origin_h3.py", [], ["3_analysis/integration/pwm_origin_h3.csv"]),
    ("3_analysis/integration/make_numbers.py", [], ["outputs/numbers/NUMBERS_human_35TF.csv"]),
    ("3_analysis/integration/make_numbers_w2.py", [], ["outputs/numbers/NUMBERS_R2_35TF.csv"]),
]


def _snapshot(root):
    snap = {}
    for d, dirs, files in os.walk(root):
        rel = os.path.relpath(d, root)
        if rel.split(os.sep)[0] == "tests":
            dirs[:] = []
            continue
        for f in files:
            p = os.path.join(d, f)
            try:
                st = os.lstat(p)
            except FileNotFoundError:
                continue
            snap[os.path.relpath(p, root)] = (st.st_size, st.st_mtime_ns)
    return snap


@pytest.fixture(scope="module")
def rerun(tmp_path_factory):
    before = _snapshot(T.R2)
    base = tmp_path_factory.mktemp("mirror")
    mirror = str(base / "manuscript")
    os.makedirs(mirror)
    for e in os.listdir(T.R2):
        src = os.path.join(T.R2, e)
        if e in SKIP_DIRS:
            continue
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(mirror, e), symlinks=True)
        else:
            shutil.copy2(src, os.path.join(mirror, e))
    shutil.copytree(os.path.join(T.R2, ROUND1), os.path.join(str(base), ROUND1), symlinks=True)
    mut = os.environ.get("T1_MUTATE", "")
    if "hm_mm_changed" in mut.split(","):             # red run: the rebuilt mouse results change after the tables were built
        f = os.path.join(mirror, os.path.relpath(T.h3in().IN["hm_mm"], T.R2))   # the final-run results (h3_inputs hm_mm)
        d = pd.read_csv(f)
        for c in [c for c in d.columns if c.startswith(("MM_rf_", "HM_rf_")) and c.endswith(("auroc_s0", "auprc_mean", "auroc_mean", "auprc_s0"))]:
            d[c] = d[c] * 0.9
        d.to_csv(f, index=False)
    if "pkg_t3_changed" in mut.split(","):            # red run: the package Sup. Table 3 summary changes (Fig. 1 values)
        f = os.path.join(mirror, "3_analysis", "package", "Sup_Tables", "Sup_Table_3_ArChIPelago_performance.xlsx")
        t3 = pd.read_excel(f, sheet_name=None)
        t3["summary"].loc[0, "ArChIPelago auROC (median)"] += 0.01
        with pd.ExcelWriter(f) as xw:
            for sh, d in t3.items():
                d.to_excel(xw, sheet_name=sh, index=False)
    if "builder_input_changed" in mut.split(","):     # red run: an input of Sup. Table 7 changes after the table was built
        f = os.path.join(mirror, "3_analysis", "low_fpr_metrics", "B2_metrics_wide.csv")
        d = pd.read_csv(f)
        d.iloc[0, d.columns.get_loc(d.select_dtypes("number").columns[0])] += 1
        d.to_csv(f, index=False)
    logs = {}
    env = dict(os.environ, OMP_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg")
    for script, args, _ in BUILDERS:
        path = os.path.join(mirror, script)
        p = subprocess.run([sys.executable, path] + args, capture_output=True, text=True, env=env, cwd=os.path.dirname(path))
        logs[script] = p
    if "write_real" in mut.split(","):               # red run of test_real_tree_not_written: a probe file in the real tree
        open(os.path.join(T.R2, "_t1_mutation_probe.txt"), "w").write("probe\n")
    after = _snapshot(T.R2)
    return dict(mirror=mirror, logs=logs, before=before, after=after)


def test_real_tree_not_written(rerun):
    changed = sorted(k for k in set(rerun["before"]) | set(rerun["after"]) if rerun["before"].get(k) != rerun["after"].get(k))
    assert not changed, changed[:20]


def _frames(path):
    if path.endswith(".xlsx"):
        return pd.read_excel(path, sheet_name=None)
    return {"csv": pd.read_csv(path)}


@pytest.mark.parametrize("script,args,outputs", BUILDERS, ids=[os.path.basename(b[0]) for b in BUILDERS])
def test_builder_reproduces_shipped_output(rerun, script, args, outputs):
    p = rerun["logs"][script]
    assert p.returncode == 0, "%s failed:\n%s\n%s" % (script, p.stdout[-2000:], p.stderr[-2000:])
    diffs = []
    for rel in outputs:
        new, old = _frames(os.path.join(rerun["mirror"], rel)), _frames(os.path.join(T.R2, rel))
        if list(new) != list(old):
            diffs.append("%s: sheets %s vs %s" % (rel, list(new), list(old)))
            continue
        for sh in new:
            try:
                pd.testing.assert_frame_equal(new[sh], old[sh], check_dtype=False)
            except AssertionError as e:
                diffs.append("%s/%s: %s" % (rel, sh, str(e).splitlines()[0:6]))
    assert not diffs, diffs
