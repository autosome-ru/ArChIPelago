"""Tests of Manuscript_analysis/: the headline numbers recomputed from results_table.csv, and the table builders
re-run on a copy of the folder, which must reproduce every tracked csv and json file byte for byte.

Run from the repository root: pytest Manuscript_analysis/tests
"""
import difflib
import filecmp
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.stats import wilcoxon

RES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# builders of rebuild_tables_and_figures.sh that need neither the Zenodo PWM files nor R
# (motif_similarity.py and build_Figure_S7_source_data.py read PWM files; their outputs are tracked inputs here)
BUILDERS = [
    "scripts/make_results_table.py",
    "analysis/saturation/analyze_saturation.py",
    "scripts/make_sup_tables_1_2.py",
    "scripts/make_sup_tables_3_to_6.py",
    "analysis/slim_dichipmunk/assemble_Figure_4_S4.py",
    "scripts/make_figure_source_data.py",
    "scripts/Figure_S2_RF_vs_best_monoPWM.py",
]


@pytest.fixture(scope="module")
def results():
    return pd.read_csv(os.path.join(RES, "results_table.csv"), sep="\t")


@pytest.fixture(scope="module")
def numbers():
    with open(os.path.join(RES, "numbers_in_text.json")) as fh:
        return json.load(fh)


def rf_monodi(results):
    rf = results[(results.Model == "RandomForestClassifier") & (results.PWM == "mono+di")]
    assert len(rf) == 36 and rf.TF_name.is_unique
    return rf.set_index("TF_name").sort_index()


def test_results_table_shape(results):
    assert results.shape == (900, 47)
    assert results.TF_name.nunique() == 36


@pytest.mark.parametrize("species", ["H", "M"])
def test_headline_numbers(results, numbers, species):
    rf = rf_monodi(results)
    n = numbers["RF_" + species]
    auroc, auprc = rf["roc_auc_test_" + species], rf["pr_auc_test_" + species]
    base_auroc, base_auprc = rf["roc_auc_test_%s_PWM" % species], rf["pr_auc_test_%s_PWM" % species]
    assert round(auroc.median(), 4) == n["rf_auroc_median"]
    assert round(auprc.median(), 4) == n["rf_auprc_median"]
    assert round(base_auroc.median(), 4) == n["base_auroc_median"]
    assert round(base_auprc.median(), 4) == n["base_auprc_median"]
    assert round(np.median(auroc - base_auroc), 4) == n["median_per_tf_auroc"]
    assert round(np.median(auprc - base_auprc), 4) == n["median_per_tf_auprc"]
    assert int((auroc > base_auroc).sum()) == n["n_improved_auroc"]
    assert int((auprc > base_auprc).sum()) == n["n_improved_auprc"]
    p = [wilcoxon(auroc, base_auroc).pvalue, wilcoxon(auprc, base_auprc).pvalue]
    assert np.allclose(p, n["wilcoxon_p"], rtol=1e-6)


def test_rebuild_reproduces_tracked_outputs(tmp_path):
    copy = str(tmp_path / "Manuscript_analysis")
    shutil.copytree(RES, copy, ignore=shutil.ignore_patterns("__pycache__", ".pytest_cache", "tests"))
    for script in BUILDERS:
        run = subprocess.run([sys.executable, os.path.join(copy, script)], capture_output=True, text=True,
                             env=dict(os.environ, MPLBACKEND="Agg"))
        assert run.returncode == 0, "%s failed:\n%s" % (script, run.stderr[-3000:])
    changed, compared = [], 0
    for root, _, files in os.walk(RES):
        if "tests" in root.split(os.sep) or "__pycache__" in root:
            continue
        for f in files:
            if f.endswith((".csv", ".json")):
                rel = os.path.relpath(os.path.join(root, f), RES)
                compared += 1
                if not filecmp.cmp(os.path.join(RES, rel), os.path.join(copy, rel), shallow=False):
                    changed.append(rel)
    assert compared > 100
    detail = []
    for rel in changed:
        with open(os.path.join(RES, rel)) as a, open(os.path.join(copy, rel)) as b:
            diff = difflib.unified_diff(a.read().splitlines(), b.read().splitlines(), rel, "rebuilt", n=0, lineterm="")
            detail.append("\n".join(list(diff)[:12]))
    assert not changed, "outputs differ from the tracked files: %s\n%s" % (changed, "\n".join(detail))
