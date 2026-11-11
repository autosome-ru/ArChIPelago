"""3_analysis/integration/pwm_origin_h3.csv is recomputed here independently of its producer: every PWM file of the
four Zenodo v1.0.0 PWM archives is read, its source peak set is the PEAKS id in the first line, and it counts as 'from a listed
set' when that id is a GTRD peak set ID of the same TF and species in Sup. Table 1. A toy check pins the name parsing.
Red: T1_MUTATE=pwm_origin_drop_listed removes one listed set of SRF (human) before counting.
"""
import os
import re
import tarfile

import pandas as pd
import pytest

import t1_helpers as T

pytestmark = pytest.mark.r2
ZEN = os.path.join(T.R2, "results/datasets/zenodo_v1.0.0")
PEAKS = re.compile(r"PEAKS\d+")


def source_set(first_line):
    ids = PEAKS.findall(first_line)
    assert len(ids) == 1, first_line
    return ids[0]


def test_source_set_parsing_toy():
    assert source_set(">SRF_MOUSE~CM~SRF_MOUSE.PEAKS035686.macs.M") == "PEAKS035686"
    assert source_set("SRF_MOUSE~CM~SRF_MOUSE_PEAKS035686_macs_M") == "PEAKS035686"
    with pytest.raises(AssertionError):
        source_set(">SRF_HUMAN~CM~SRF_HUMAN.macs.M")


@pytest.fixture(scope="module")
def recomputed():
    if not os.path.isdir(ZEN):
        pytest.skip("Zenodo archives not present")
    t1 = pd.read_excel(os.path.join(T.R2, "outputs/supplementary_tables", "Sup_Table_1_ChIP-Seq_experiments.xlsx"), sheet_name=None)
    rows = {}
    for sp in ("human", "mouse"):
        listed = {tf: set(g["GTRD peak set ID"]) for tf, g in t1[sp].groupby("TF")}
        if T.mutate("pwm_origin_drop_listed") and sp == "human":
            listed["SRF"] = set(sorted(listed["SRF"])[1:])
        for k in ("mono", "di"):
            with tarfile.open(os.path.join(ZEN, "PWMs_%s_%s.tar.gz" % (k, sp.upper()))) as tar:
                for m in tar.getmembers():
                    if not (m.isfile() and m.name.endswith(("pwm", "dpwm"))):
                        continue
                    tf = m.name.split("/")[-2]
                    if tf not in listed:
                        continue
                    s = source_set(tar.extractfile(m).readline().decode())
                    r = rows.setdefault((tf, sp), [0, 0])
                    r[0] += 1
                    r[1] += s in listed[tf]
    return rows


def test_pwm_origin_table_matches_the_archives(recomputed):
    d = pd.read_csv(os.path.join(T.R2, "3_analysis/integration", "pwm_origin_h3.csv")).set_index(["TF", "species"])
    assert set(d.index) == set(recomputed), sorted(set(d.index) ^ set(recomputed))
    bad = [(k, v, (d.loc[k, "n_PWMs"], d.loc[k, "n_from_listed_sets"])) for k, v in recomputed.items()
           if (d.loc[k, "n_PWMs"], d.loc[k, "n_from_listed_sets"]) != tuple(v)]
    assert not bad, bad[:5]
    assert set(d.reset_index().TF) == set(T.TF35)
