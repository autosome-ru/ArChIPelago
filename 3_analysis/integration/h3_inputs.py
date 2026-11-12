"""The one registry of final-data inputs: no producer reads a first-submission human table for a current value. Every
producer of the integration chain takes its inputs from here and records them in its provenance block.

Final data, selected by ARCHI_STAGE (any value other than H4 / H5 stops; H5 is the final data):
  ARCHI_STAGE=H5: chr11, chr12 (and human chr22) are in TRAINING; the test sets are byte-identical to the previous training run. Every
      model and evaluation input is read from results/server (human, mouse); the data-build inputs and the inputs that were
      mirrored only for the previous run come from DATA_H5 (results/server/human_data, mouse_data); exceptions with their
      reason in NOT_IN_H5.
  ARCHI_STAGE=H4: the previous training run (models retrained, negatives v3, test sets scaled with the training scaler), read from
      results/server_h4_test_sets. Its test-set tables are in the repo; its data-build inputs are not (see H4_ONLY_DATA).
OLD maps a quantity to the file that holds its first-submission or pre-rebuild value (previous values, round 1 package and
Resubmission tables).
usage: import h3_inputs as H; H.IN["hm_mm"]; H.check() raises when a registry path is missing; H.STAGE, H.SPLIT, H.COL.
Pre-flight (presence, md5 against the mirror manifest, split labels): python3 3_analysis/integration/preflight_inputs.py
"""
import os

R2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PKG = os.path.join(R2, "results/round1_package/Manuscript_analysis")
R1 = os.path.normpath(os.path.join(R2, "results", "round1_package", "Resubmission"))
STAGES = {"H4": "results/server_h4_test_sets", "H5": "results/server"}
STAGE = os.environ.get("ARCHI_STAGE", "H4")
if STAGE not in STAGES:
    raise SystemExit("CHECK FAILED h3_inputs: ARCHI_STAGE=%r, expected one of %s" % (STAGE, sorted(STAGES)))
ROOT = STAGES[STAGE]
MANIFEST = ROOT + "/MIRROR_MANIFEST_md5.txt"
# Every test-set evaluation comes from the run with the training-fitted scaler. Models, training metrics, PWM baselines, Gini
# importances and feature-importance subtypes are identical in the previous run.
H3H, H3M = ROOT + "/human", ROOT + "/mouse"
# Inputs that the previous run did not rerun or mirror because they involve no test-set scaling; they stay at the older mirror
# (listed so that tests/test_inputs_h4.py can tell an allowed path from a stale one)
H3A_ONLY = "B15_h3a/human_sym"
NOT_IN_H4 = {"gc_v3": "negatives v3 GC check (no model)", "neg_v3_stats": "negatives v3 statistics (no model)",
             "neg_v3_checks": "negatives v3 checks (no model)", "runtime": "Sup. Table 6 runtime benchmark (timing only, training scaler)",
             "b4_sim": "B4 PWM similarity matrices in Gini-rank order (PWMs only; .npy not mirrored; valid while the Gini files are "
                       "byte-identical, checked by b4_05_aggregate.py)"}
# Stage-specific inputs outside the model-run layout. ARCHI_STAGE=H4 reads the data builds and the older mirror; these inputs
# are not in the repo, so that branch cannot run from a checkout (it stays in the code). The final run rebuilds them, and its
# mirror provides them at the DATA_H5 paths (the pre-flight lists each one).
H4_ONLY_DATA = {
    "counts_H": "B14_human_sym/counts_per_tf.csv", "pos_H": "B14_human_sym/per_tf_positive_counts.csv",
    "counts_M": "B11_mouse_sym/counts_per_tf.csv", "pos_M": "B11_mouse_sym/per_tf_positive_counts.csv",
    "checkC_H": "B14_human_sym/logs/check_C.log", "checkC_M": "B11_mouse_sym/logs/check_C.log",
    "build_positives": "B14_human_sym/code/build_positives_sym.py",
    "gc_v3": H3A_ONLY + "/neg_v3/checks_neg_v3/gc_check_v3.csv", "neg_v3_stats": H3A_ONLY + "/neg_v3/checks_neg_v3/neg_v3_stats_all.csv",
    "neg_v3_checks": H3A_ONLY + "/neg_v3/checks_neg_v3/neg_v3_checks.csv", "runtime": H3A_ONLY + "/runtime",
    "b4_sim": H3A_ONLY + "/B4_importance/cache",
}
DATA_H5 = {
    "counts_H": "results/server/human_data/counts_per_tf.csv", "pos_H": "results/server/human_data/per_tf_positive_counts.csv",
    "counts_M": "results/server/mouse_data/counts_per_tf.csv", "pos_M": "results/server/mouse_data/per_tf_positive_counts.csv",
    "checkC_H": "results/server/human_data/logs/check_C.log", "checkC_M": "results/server/mouse_data/logs/check_C.log",
    "build_positives": "1_data_preparation/h5_splits.py",                  # the one split definition (replaces the literal chromosome lists)
    "gc_v3": "results/server/human/neg_v3/checks_neg_v3/gc_check_v3.csv",
    "neg_v3_stats": "results/server/human/neg_v3/checks_neg_v3/neg_v3_stats_all.csv",
    "neg_v3_checks": "results/server/human/neg_v3/checks_neg_v3/neg_v3_checks.csv",
    "runtime": "results/server/human/runtime",                     # make_h5_runtime.sh (training matrices grow)
    "b4_sim": "results/server/human/B4_importance/cache",          # .npy in the Gini order of the final run (the previous run did not mirror them)
}
if set(DATA_H5) != set(H4_ONLY_DATA):
    raise SystemExit("CHECK FAILED h3_inputs: DATA_H5 and H4_ONLY_DATA name different keys")
_D = DATA_H5 if STAGE == "H5" else H4_ONLY_DATA
# Inputs the final run does not reproduce, read from the mirror of the previous run on purpose (reason each). Every other staged key is from results/server.
NOT_IN_H5 = {
    "b3_widths_H": "peak widths of the test datasets: test sets byte-identical to H4 (PLAN_H5 step 8, not rerun)",
    "b3_widths_M": "as b3_widths_H, mouse",
    "b3_join_H": "MACS join of the test windows: test sets byte-identical to H4 (PLAN_H5 step 8, not rerun)",
    "b3_join_M": "as b3_join_H, mouse",
}
_B3H, _B3M = (STAGES["H4"] + "/human", STAGES["H4"] + "/mouse") if STAGE == "H5" else (H3H, H3M)
# inputs that do not depend on the stage (reason each); test_inputs_h5 asserts that every key is staged, DATA or FIXED
FIXED = {
    "ids_H": "dataset lists (D16 changes chromosomes only)", "ids_M": "dataset lists (D16 changes chromosomes only)",
    "b8_summary": "random split of the top 1,000 peaks per dataset: independent of the chromosome lists", "b8_top": "as b8_summary",
    "b8_per_tf": "as b8_summary", "mask_readme": "repeat mask description", "figS4_num": "round-1 Fig. S4 numbers (old text values)",
    "b9_summary": "B13 copy that Figure_S12.py writes from b9_h3b (staged)", "b9_per_fit": "as b9_summary",
}

# chromosome split: the previous run's lists as its data build made them (chr11, chr12 in neither set); the final run = h5_splits.LOGICAL
SPLIT_LISTS = {
    "H4": {"H": ([2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20], [1, 8, 21]),
           "M": ([2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18], [1, 8, 19])},
    "H5": {"H": ([i for i in range(1, 23) if i not in (1, 8, 21)], [1, 8, 21]),
           "M": ([i for i in range(1, 20) if i not in (1, 8, 19)], [1, 8, 19])},
}


def chrom_label(nums, sep=", "):
    """[2, 3, 4, 5, 6, 7, 9, 10, 13, ..., 20] -> 'chr2-7, 9, 10, 13-20' (runs of three or more written as ranges)."""
    nums, out, i = sorted(nums), [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out += ["%d-%d" % (nums[i], nums[j])] if j - i >= 2 else ["%d" % n for n in nums[i:j + 1]]
        i = j + 1
    return "chr" + sep.join(out)


SPLIT = {"%s_%s" % (sp, part): chrom_label(SPLIT_LISTS[STAGE][sp][k]) for sp in "HM" for k, part in ((0, "train"), (1, "test"))}
# labels of the 'chromosomes' column of gc_check_v3.csv (previous run's file: 'chr2-7,9,10,13-20', 'chr1/8/21'; the pre-flight checks them)
SPLIT_GC = {}
for _sp in "HM":
    SPLIT_GC[_sp + "_train"] = chrom_label(SPLIT_LISTS[STAGE][_sp][0], ",")
    SPLIT_GC[_sp + "_test"] = "chr" + "/".join(str(c) for c in SPLIT_LISTS[STAGE][_sp][1])
# pool names of gc_check_v3.csv (the final data build labels its pools *_sym_h5_v3)
GC_POOL = {"H4": {"H": "human_sym_v3", "M": "mouse_sym_v3"}, "H5": {"H": "human_sym_h5_v3", "M": "mouse_sym_h5_v3"}}[STAGE]
# Sup. Table 2 column names (read by Sup. Tables 2 and 5, Fig. 1 values, NUMBERS, checks)
COL = {"HTR": "Human training positives (%s)" % SPLIT["H_train"], "HTE": "Human test positives (%s)" % SPLIT["H_test"],
       "MTR": "Mouse training positives (%s)" % SPLIT["M_train"], "MTE": "Mouse test positives (%s)" % SPLIT["M_test"]}

IN = {k: os.path.join(R2, v) for k, v in {
    # models and evaluations
    "eval_HM": H3H + "/eval_models/HM/evaluation_table.csv",            # human-trained: test_H (H->H), test_M (H->M); tab
    "best_HM": H3H + "/eval_models/HM/best_single_PWMs.csv",
    "eval_MM": H3M + "/eval_models/MM/evaluation_table.csv",            # mouse-trained: test_M (M->M), test_H (M->H); tab
    "best_MM": H3M + "/eval_models/MM/best_single_PWMs.csv",
    "hm_mm": H3H + "/hm_mm_results.csv",                                # RF seeds 0/1, baselines, n per split (both species)
    "operational": H3H + "/operational/operational_metrics.csv",
    "saturation_H": H3H + "/saturation_HUMAN/saturation_results.csv",   # human-trained, test H and test M
    "saturation_M": H3M + "/saturation_MOUSE/saturation_results.csv",   # mouse-trained (M->M, M->H)
    "b2_wide": H3H + "/b2_v2/B2_metrics_wide.csv", "b2_long": H3H + "/b2_v2/B2_metrics_per_TF_long.csv",   # b2_v2: pAUC fix
    "b2_chrom": H3H + "/b2_v2/B2_per_chromosome.csv", "b2_dir": H3H + "/b2_v2",
    "b3_topN": H3H + "/b3/b3_topN_results.csv", "b3_signal": H3H + "/b3/b3_topN_signal_summary.csv",
    "b3_widths_H": _B3H + "/b3/b3_peak_widths_per_dataset_human.csv", "b3_widths_M": _B3M + "/b3/b3_peak_widths_per_dataset_mouse.csv",
    "b3_join_H": _B3H + "/b3/b3_join_test_per_tf_human.csv", "b3_join_M": _B3M + "/b3/b3_join_test_per_tf_mouse.csv",
    "b4_dir": H3H + "/B4_importance", "b4_sim": _D["b4_sim"],
    # data builds (positives, PWM counts) and negatives v3
    "counts_H": _D["counts_H"], "pos_H": _D["pos_H"],
    "counts_M": _D["counts_M"], "pos_M": _D["pos_M"], "checkC_H": _D["checkC_H"], "checkC_M": _D["checkC_M"],
    "gc_v3": _D["gc_v3"], "neg_v3_stats": _D["neg_v3_stats"],
    "neg_v3_checks": _D["neg_v3_checks"], "runtime": _D["runtime"],
    "build_positives": _D["build_positives"],    # chromosome lists of the split; final run: h5_splits.py
    # dataset lists
    "ids_H": "results/datasets/peak_ids/human_used_peak_ids.csv", "ids_M": "results/datasets/peak_ids/mouse_sym_peak_ids.csv",
    # unchanged by the retraining (peak-level analyses)
    "b8_summary": "results/datasets/split_leakage/b8w1g2_summary.csv", "b8_top": "results/datasets/split_leakage/b8w1g2_top1000_check.csv",
    "b8_per_tf": "results/datasets/split_leakage/b8w1g2_leakage_per_TF_masked.csv", "mask_readme": "results/datasets/repeat_mask/repeat_mask_counts.csv",
    # Fig. 4 / S4: RF2f refits with the training scaler; md5 of every mirrored file in MIRROR_MANIFEST_md5.txt
    "fig4_src": H3H + "/fig4_refit_trainscale/Figure_4_source_data.csv", "figS4_src": H3H + "/fig4_refit_trainscale/Figure_S4_source_data.csv",
    "fig4_num": H3H + "/fig4_refit_trainscale/Figure_4_and_S4_numbers.csv", "fig4_prov": MANIFEST,
    # first-submission Fig. S4 numbers (36 TFs; old text values of the Fig. S4 NUMBERS rows only)
    "figS4_num": "results/datasets/round1_fig_s4/Figure_S4_numbers.csv",
    # random split (Fig. S12, Sup. Table 7): retraining on the v3 negatives, copied and checked by 4_visualisation/supplementary_figures/Figure_S12/Figure_S12.py
    "b9_summary": "3_analysis/random_split/b9_summary.csv", "b9_per_fit": "3_analysis/random_split/b9_random_split_per_fit.csv",
    "b9_h3b": H3H + "/H3b/b9/b9_summary.csv", "b9_dir": H3H + "/H3b/b9",
}.items()}
PENDING = set()


def check():
    miss = [k for k, v in IN.items() if not os.path.exists(v)]
    if miss:
        raise SystemExit("CHECK FAILED h3_inputs: missing %s" % miss)
    return len(IN)


if __name__ == "__main__":
    print("h3_inputs: %d paths, all present" % check())
