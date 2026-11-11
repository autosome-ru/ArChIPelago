"""Copy the final server tables into the package inputs (3_analysis/package). The package producers are committed in
3_analysis/package, already patched for 35 TFs; this script only swaps their input tables for the registry inputs of
h3_inputs.py.

Input swaps (file in the package <- registry key; every swapped table is checked to hold exactly the 35 reference TFs; old
package files are backed up first):
  analysis/inputs/notebook4_results_table.csv        <- eval_HM (human-trained evaluation)
  analysis/model_evaluation/evaluation_table.csv     <- eval_HM
  analysis/model_evaluation/best_single_PWMs.csv     <- best_HM
  analysis/mouse_transfer/mouse_transfer_results.csv <- hm_mm
  analysis/operational/operational_metrics.csv       <- operational
  analysis/saturation/saturation_results.csv         <- saturation_H (human-trained)
  analysis/inputs/training_positives_per_TF.csv      <- pos_H (train_file_pos), checked against the Seq_count of the evaluation
  analysis/cross_species/motif_similarity.csv        -> removed; recomputed by motif_similarity.py (baselines changed)
  Figures/source_data/Figure_4_source_data.csv, Figure_S4_source_data.csv <- fig4_src, figS4_src
  analysis/cross_species/cross_species_table.csv     -> backed up and removed (rebuilt by the producer)
usage: python3 3_analysis/integration/make_pkg_copy.py   (then run_h3.py runs the producers)
"""
import difflib
import os
import shutil
import sys

import pandas as pd

B13 = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(B13))
sys.path.insert(0, os.path.join(R2, "3_analysis/integration"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p53drop as PD  # noqa: E402
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402  (final-data input registry)

P = os.path.join(R2, "3_analysis/package")
LOG = os.path.join(B13, "logs")
os.makedirs(LOG, exist_ok=True)

ins, outs = [], []


def swap(rel, src, sep, col):
    dst = os.path.join(P, rel)
    shutil.copy2(src, dst)
    PD.check_tfs(pd.read_csv(dst, sep=sep)[col], rel)
    ins.append(src); outs.append(dst)


H.check()
for _rel in ("analysis/inputs/notebook4_results_table.csv", "analysis/model_evaluation/evaluation_table.csv", "analysis/model_evaluation/best_single_PWMs.csv",
             "analysis/mouse_transfer/mouse_transfer_results.csv", "analysis/operational/operational_metrics.csv", "analysis/saturation/saturation_results.csv",
             "analysis/inputs/training_positives_per_TF.csv", "analysis/cross_species/motif_similarity.csv"):
    w2.backup(os.path.join(P, _rel))
# The human values come from the retrained evaluation (human-trained models, negatives v3). The first-submission notebook 4
# table is no longer an input: both make_results_table sources are the human-trained evaluation table.
swap("analysis/inputs/notebook4_results_table.csv", H.IN["eval_HM"], "\t", "TF_name")
swap("analysis/model_evaluation/evaluation_table.csv", H.IN["eval_HM"], "\t", "TF_name")
swap("analysis/model_evaluation/best_single_PWMs.csv", H.IN["best_HM"], "\t", "TF_name")
swap("analysis/mouse_transfer/mouse_transfer_results.csv", H.IN["hm_mm"], ",", "TF")
swap("analysis/operational/operational_metrics.csv", H.IN["operational"], ",", "tf")
swap("analysis/saturation/saturation_results.csv", H.IN["saturation_H"], ",", "tf")
swap("Figures/source_data/Figure_4_source_data.csv", H.IN["fig4_src"], ",", "TF")       # server mirror
swap("Figures/source_data/Figure_S4_source_data.csv", H.IN["figS4_src"], ",", "TF")     # server mirror
# training positives per TF (Fig. 2 filter Seq_count > 100): human training positives (file "train")
_tp = pd.read_csv(H.IN["pos_H"]).set_index("TF").sort_index()
PD.check_tfs(_tp.index, "B14 per_tf_positive_counts")
_tpd = os.path.join(P, "analysis/inputs/training_positives_per_TF.csv")
pd.DataFrame({"Names": _tp.index, "Seq_count": _tp.train_file_pos.values}).to_csv(_tpd, sep="\t")
ins.append(H.IN["pos_H"]); outs.append(_tpd)
_ev = pd.read_csv(H.IN["eval_HM"], sep="\t")
_sc = _ev.drop_duplicates("TF_name").set_index("TF_name").Seq_count.sort_index()
if not (_sc == _tp.train_file_pos).all():
    raise SystemExit("CHECK FAILED H3a Seq_count vs B14 training positives: %s" % list(_sc.index[_sc != _tp.train_file_pos]))
_ms = os.path.join(P, "analysis/cross_species/motif_similarity.csv")
if os.path.exists(_ms):
    os.remove(_ms)                               # backed up above; rewritten by the patched motif_similarity.py (run_h3.py)
cs = os.path.join(P, "analysis/cross_species/cross_species_table.csv")
R1CS = os.path.join(P, "analysis/cross_species/cross_species_table_round1.csv")
if os.path.exists(cs):
    if not os.path.exists(R1CS):
        os.rename(cs, R1CS)                      # first run: the released table is the first-submission table
    else:
        w2.backup(cs); os.remove(cs)             # later runs: keep the table as a backup, the producer refuses to overwrite

w2.provenance(P, "3_analysis/integration/make_pkg_copy.py", outs, ins,
                    note="H3: PKG copy with the final inputs (h3_inputs.py); patched producers (asserted replacements)")
