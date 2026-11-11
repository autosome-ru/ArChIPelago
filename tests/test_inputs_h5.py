"""Stage switch of the input registry (ARCHI_STAGE=H4 or H5), the pre-flight, the split labels and the flip report. The final
run (ARCHI_STAGE=H5) adds chr11, chr12 (and human chr22) to training; its test sets are byte-identical to those of the
previous run (ARCHI_STAGE=H4).
Red: T1_MUTATE=h5_stage_leak points one staged final-run entry (hm_mm) back at results/server_h4_test_sets; h5_preflight_md5
changes one manifest md5 seen by the pre-flight; h5_split_old gives the final run the chromosome lists of the previous run;
flip_blind runs the flip report with alpha 0 (no significance flip)."""
import importlib.util
import os
import subprocess
import sys

import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
CODE = os.path.join(T.R2, "3_analysis", "integration")


def load(name, stage):
    """Module CODE/<name>.py imported fresh with ARCHI_STAGE=stage (h3_inputs is imported fresh for it too)."""
    old = os.environ.get("ARCHI_STAGE")
    os.environ["ARCHI_STAGE"] = stage
    try:
        sys.modules.pop("h3_inputs", None)
        sys.path.insert(0, CODE)
        spec = importlib.util.spec_from_file_location("%s_%s_t1" % (name, stage), os.path.join(CODE, name + ".py"))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m
    finally:
        sys.path.remove(CODE)
        sys.modules.pop("h3_inputs", None)
        if old is None:
            os.environ.pop("ARCHI_STAGE", None)
        else:
            os.environ["ARCHI_STAGE"] = old


def rel(p):
    return os.path.relpath(p, T.R2).replace(os.sep, "/")


def test_h5_resolves_every_staged_input_under_b17_with_the_b16_relative_path():
    h4, h5 = load("h3_inputs", "H4"), load("h3_inputs", "H5")
    IN5 = dict(h5.IN)
    if T.mutate("h5_stage_leak"):
        IN5["hm_mm"] = h4.IN["hm_mm"]
    assert (h4.STAGE, h4.ROOT, h5.STAGE, h5.ROOT) == ("H4", "results/server_h4_test_sets", "H5", "results/server")
    assert set(h4.IN) == set(IN5)
    bad, n = [], {"staged": 0, "data": 0, "not_in_H5": 0, "fixed": 0}
    for k in sorted(h4.IN):
        a, b = rel(h4.IN[k]), rel(IN5[k])
        if k in h5.FIXED:
            n["fixed"] += 1; ok = a == b and not a.startswith(("results/server/", "results/server_h4_test_sets/"))
        elif k in h5.NOT_IN_H5:
            n["not_in_H5"] += 1; ok = a == b and a.startswith("results/server_h4_test_sets/")
        elif k in h5.DATA_H5:
            n["data"] += 1; ok = b == h5.DATA_H5[k] and a == h4.H4_ONLY_DATA[k]
        else:
            n["staged"] += 1; ok = a.startswith("results/server_h4_test_sets/") and b == "results/server/" + a.split("/", 2)[2]
        if not ok:
            bad.append((k, a, b))
    print("H5 registry: %s" % n)
    assert not bad, bad
    assert n["staged"] >= 20 and all(h5.NOT_IN_H5.values()) and all(h5.FIXED.values())
    assert not set(h5.NOT_IN_H5) & set(h5.DATA_H5) and not set(h5.FIXED) & set(h5.DATA_H5)


def test_unknown_stage_stops():
    p = subprocess.run([sys.executable, "-c", "import h3_inputs"], cwd=CODE, env=dict(os.environ, ARCHI_STAGE="H6"),
                       capture_output=True, text=True)
    assert p.returncode != 0 and "ARCHI_STAGE='H6'" in (p.stdout + p.stderr)


def test_preflight_passes_on_the_h5_inputs():
    pf5 = load("preflight_inputs", "H5")
    if T.mutate("h5_preflight_md5"):
        real = pf5.manifests

        def fake():
            m = real()
            k = "human/hm_mm_results.csv"
            assert k in m["results/server"]
            m["results/server"][k] = "0" * 32
            return m
        pf5.manifests = fake
    r5 = {r["key"]: r for r in pf5.rows()}
    assert all(r["ok"] for r in r5.values()), [r for r in r5.values() if not r["ok"]]
    assert r5["b3_join_H"]["cls"] == "not_in_H5" and r5["b3_join_H"]["ok"]
    p = subprocess.run([sys.executable, os.path.join(CODE, "preflight_inputs.py")], env=dict(os.environ, ARCHI_STAGE="H5"),
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-500:]


def test_split_labels_h4_unchanged_and_h5_trains_on_every_autosome():
    h4, h5 = load("h3_inputs", "H4"), load("h3_inputs", "H5")
    lists5 = h5.SPLIT_LISTS["H4" if T.mutate("h5_split_old") else "H5"]
    t2 = pd.read_excel(os.path.join(T.R2, "outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx"), sheet_name="per_TF")
    cur, other = (h5, h4) if T.h3in().STAGE == "H5" else (h4, h5)
    if T.mutate("t2_stale_col"):
        t2[other.COL["HTR"]] = 0
    assert set(cur.COL.values()) <= set(t2.columns), sorted(set(cur.COL.values()) - set(t2.columns))
    # no count column of the other stage survives next to the current one
    stale = sorted(set(other.COL.values()) - set(cur.COL.values()) & set(t2.columns))
    assert not stale, stale
    assert h4.COL["HTR"] == "Human training positives (chr2-7, 9, 10, 13-20)"
    assert (h5.SPLIT["H_train"], h5.SPLIT["M_train"], h5.SPLIT["H_test"], h5.SPLIT["M_test"]) == (
        "chr2-7, 9-20, 22", "chr2-7, 9-18", "chr1, 8, 21", "chr1, 8, 19")
    ns = {}
    exec(open(h5.IN["build_positives"]).read(), ns)   # 1_data_preparation/h5_splits.py, the one split definition
    for sp, SP in (("H", "HUMAN"), ("M", "MOUSE")):
        tr, te = ["chr%d" % c for c in lists5[sp][0]], ["chr%d" % c for c in lists5[sp][1]]
        assert sorted(tr) == sorted(ns["LOGICAL"][SP]["train"]) and sorted(te) == sorted(ns["LOGICAL"][SP]["test"]), SP
        assert ns["partition_problems"](SP, tr, te) == [], SP


def test_flip_report_finds_the_known_h3a_to_h4_flips():
    fr = load("flip_report", "H4")
    d = pd.read_csv(os.path.join(T.R2, "outputs/numbers/NUMBERS_R2_35TF.csv"), dtype=str, keep_default_na=False)
    if T.h3in().STAGE == "H5":   # the previous run's values sit in "h4 value"; the flips are judged on them
        d = d.assign(**{"new value": d["h4 value"]})
    rep = fr.report(d, "h3a value", 0.0 if T.mutate("flip_blind") else 0.05)
    got = set(zip(rep.item, rep.kind))
    for want in (("Human TFs improved, auPRC (n of N)", "count"), ("Headroom Spearman P, dauPRC vs baseline auPRC", "signif"),
                 ("drivers: base_auROC_M vs dauROC H>M (rho, nominal P)", "signif"), ("Fig. S5 random k = 2, median delta_auroc_H", "sign")):
        assert want in got, want
    # GC control P 0.34 -> 0.33 is no flip; a row equal to its reference is never reported
    assert not rep.item.str.startswith("GC control: Spearman").any()
    if "h4 value" in d.columns:
        same = d[d["h4 value"] == d["new value"]]
        assert not set(fr.report(same, "h4 value", 0.05).item)


def test_flip_report_finds_the_known_h4_to_h5_flips():
    """The flips between the previous and the final run (results/server/human/checks/flips_h5_vs_h4_v2.csv) appear in the flip report of the final run."""
    if T.h3in().STAGE != "H5":
        pytest.skip("H5 stage only")
    fr = load("flip_report", "H5")
    d = pd.read_csv(os.path.join(T.R2, "outputs/numbers/NUMBERS_R2_35TF.csv"), dtype=str, keep_default_na=False)
    rep = fr.report(d, "h4 value", 0.0 if T.mutate("flip_blind") else 0.05)
    got = set(zip(rep.item, rep.kind))
    for want in (("Table 1, Mouse, mouse-trained, auROC", "count"),
                 ("USF2 mouse-trained dauROC on the mouse test set: mean of two seeds / seed 0 (N12)", "sign"),
                 ("ERG: mouse-trained gain on the mouse test set, auROC / auPRC (mean of two seeds); human-trained dauPRC H>M (mean of two seeds)", "sign"),
                 ("Fig. 4 auROC, LSlim m=-5: two-sided Wilcoxon P", "signif")):
        assert want in got, want
    below = d.set_index("item").loc["H>M TFs below the baseline on >= 1 metric"]
    assert below["h4 value"].startswith("ERG") and below["new value"].startswith("PRGR"), dict(below)


def test_numbers_items_name_the_stage_chromosomes():
    """With ARCHI_STAGE=H5 no NUMBERS item names the training chromosomes of the previous run (labels follow h3_inputs COL / SPLIT).
    Red: T1_MUTATE=numbers_old_label gives one item the label of the previous run."""
    if T.h3in().STAGE != "H5":
        pytest.skip("H5 stage only")
    d = pd.read_csv(os.path.join(T.R2, "outputs/numbers/NUMBERS_R2_35TF.csv"), dtype=str, keep_default_na=False)
    items = list(d.item)
    if T.mutate("numbers_old_label"):
        items[0] += " (chr2-7, 9, 10, 13-20)"
    old = ("10, 13-20", "9, 10, 13-18", "10,13-20", "9,10,13-18")
    bad = [it for it in items if any(o in it for o in old)]
    assert not bad, bad
    h5 = load("h3_inputs", "H5")
    assert any(h5.SPLIT["H_train"] in it for it in items), "no item names the H5 human training chromosomes"
