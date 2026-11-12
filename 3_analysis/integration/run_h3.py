"""Integration runner: unpacks the Zenodo v1.0.0 PWM archives, then runs the producers of each stage in dependency order,
backing up every output file first (p53drop.backup), one process at a time under nice 15 with one BLAS thread.
Stops at the first failing command; logs go to 3_analysis/integration/logs/h3_<stage>_<n>.log.
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/run_h3.py all
       (or one or more stage names: see STAGES). The optional stage "gate" is not part of "all": it needs ARCHI_H5_SCORES, the
       unpacked scores/ folder of the Zenodo v2.0.0 archive scores_H5_human.tar.gz
"""
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402
import h3_inputs as H  # noqa: E402

R2, B13 = w2.R2, w2.B13
P = os.path.join(R2, "3_analysis/package")
PY = sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           ARCHI_ZENODO_DIR=os.path.join(R2, "results/datasets/zenodo_v1.0.0"), MPLBACKEND="Agg")


def j(*p):
    return os.path.join(R2, *p)


def backup_globs(*patterns):
    for pat in patterns:
        for f in sorted(glob.glob(j(pat))):
            if "_backup_" not in os.path.basename(f) and os.path.isfile(f):
                w2.backup(f)


STAGES = {
    # pkg producers (results table, saturation, motif similarity, cross-species, Sup. Tables 3-6 of the package, source data)
    "pkg": dict(backup=["3_analysis/package/results_table.csv", "3_analysis/package/Sup_Tables/*", "3_analysis/package/Figures/source_data/*",
                        "3_analysis/package/analysis/saturation/saturation_*.csv", "3_analysis/package/analysis/saturation/*.md",
                        "3_analysis/package/analysis/cross_species/*.md"],
                cmds=[([PY, "./make_pkg_copy.py"], B13),
                      ([PY, "scripts/make_results_table.py"], P),
                      ([PY, "analysis/saturation/analyze_saturation.py"], P),
                      ([PY, "analysis/cross_species/motif_similarity.py"], P),
                      ([PY, "analysis/cross_species/cross_species_table.py"], P),
                      ([PY, "scripts/make_sup_tables_3_to_6.py"], P),
                      ([PY, "scripts/make_figure_source_data.py"], P)]),
    # figures from the package tables: Fig. 2, 3 (human), S1, S3 (mouse), S2, S5, S6 (4_visualisation/figures/output)
    "figs": dict(backup=["3_analysis/integration/inputs/results_table.csv", "4_visualisation/figures/Fig2_S1/Figure_2*", "4_visualisation/figures/Fig3_S3/Figure_3_*",
                         "4_visualisation/figures/Fig3_S3/Figure_3D_*", "outputs/figures/Figure_3.p*", "4_visualisation/figures/output/S[1-6]/*", "4_visualisation/figures/output/Fig4/*"],
                 cmds=[([PY, "make_inputs.py"], j("3_analysis/integration")),
                       (["Rscript", "Fig2_S1/Figure_2_and_S1_restyled.R"], j("4_visualisation/figures"), dict(ARCHI_PKG=j("3_analysis/integration", "inputs"), ARCHI_DRAW="H")),
                       (["Rscript", "Fig2_S1/Figure_2_and_S1_restyled.R", j("4_visualisation/figures/output", "S1")], j("4_visualisation/figures"), dict(ARCHI_PKG=P, ARCHI_DRAW="M")),
                       (["Rscript", "Fig3_S3/Figure_3_and_S3_restyled.R"], j("4_visualisation/figures"), dict(ARCHI_PKG=j("3_analysis/integration", "inputs"), ARCHI_DRAW="H")),
                       (["Rscript", "Fig3_S3/Figure_3_and_S3_restyled.R", j("4_visualisation/figures/output", "S3")], j("4_visualisation/figures"), dict(ARCHI_PKG=P, ARCHI_DRAW="M")),
                       ([PY, "FigS2/Figure_S2_restyled.py", j("4_visualisation/figures/output", "S2")], j("4_visualisation/figures"), dict(ARCHI_PKG=P)),
                       (["Rscript", "FigS5/Figure_S5_restyled.R", os.path.join(P, "Figures", "source_data"), j("4_visualisation/figures/output", "S5")],
                        j("4_visualisation/figures"), dict(ARCHI_PKG=P)),
                       (["Rscript", "FigS6/Figure_S6_restyled.R", os.path.join(P, "Figures", "source_data"), j("4_visualisation/figures/output", "S6")],
                        j("4_visualisation/figures"), dict(ARCHI_PKG=P)),
                       # Fig. 4 / S4 from the package copy of the server mirror (make_pkg_copy swaps fig4_src / figS4_src)
                       (["Rscript", "Fig4_S4/Figure_4_and_S4_restyled.R", j("4_visualisation/figures/output/Fig4")], j("4_visualisation/figures"), dict(ARCHI_PKG=P)),
                       ([PY, "./fig4_tokens_h3b.py"], B13),
                       # paired RF on all PWMs vs Slim / diChIPMunk per test set (Sup. Table 3 sheet)
                       ([PY, "./fig4_paired_h4.py"], B13)]),
    # low-FPR metrics (Sup. Table 7), top peaks (Sup. Table 8, Fig. S9), feature importance (Sup. Table 9, Fig. S10) from the server mirrors
    "b234": dict(backup=["3_analysis/low_fpr_metrics/B2_*summary.csv", "3_analysis/low_fpr_metrics/B2_low_FPR_metrics.p*", "3_analysis/top_peaks/b3_*.csv", "3_analysis/top_peaks/B3_peaks_figure.p*",
                         "3_analysis/feature_importance/tables/B4_*.csv", "3_analysis/feature_importance/B4_importance_figure.p*",
                         "4_visualisation/figures/FigS7/hocomoco11_check/cluster_vs_v11*.csv"],
                 cmds=[([PY, "./make_b2_h3.py"], B13),
                       ([PY, "../low_fpr_metrics/b2_summary_w2.py"], B13),
                       ([PY, "../low_fpr_metrics/b2_figure_h3.py"], B13),
                       ([PY, "./make_b3_h3.py"], B13),
                       ([PY, "../top_peaks/b3_gain_vs_npos_h3.py", os.path.join(B13, "../top_peaks")], B13),
                       ([PY, "b3_summary_figure_w2.py"], os.path.join(B13, "../top_peaks")),
                       ([PY, "scripts/b4_05_aggregate.py"], j("3_analysis/feature_importance"), dict(B4_IN=H.IN["b4_dir"], B4_SIM=H.IN["b4_sim"])),
                       ([PY, "scripts/b4_06_figure.py"], j("3_analysis/feature_importance")),
                       # ablation subtypes vs HOCOMOCO v11; before "sup" (Sup. Table 9 sheets read it)
                       ([PY, "cluster_vs_v11.py"], j("4_visualisation/figures", "FigS7", "hocomoco11_check"))]),
    # Sup. Tables 1-5 and 7-10 (outputs/supplementary_tables)
    "sup": dict(backup=["outputs/supplementary_tables/Sup_Table_*.xlsx"],
                cmds=[([PY, "./top1000_check_h3.py"], B13),
                      ([PY, "3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py", "1", "2", "3", "4", "6"], R2),
                      ([PY, "3_analysis/supplementary_tables/make_sup_table_5_R2.py"], R2),
                      ([PY, "3_analysis/supplementary_tables/make_sup_tables_7_10.py", "7", "8", "9", "10"], R2)]),
    # Fig. 1 values, Table 1 rows, Fig. 1 (svg -> png 300 dpi, pdf); "copy" then collects all figures in outputs/figures
    "fig1": dict(backup=["4_visualisation/figure1/Figure_1_R2.*"],
                 cmds=[([PY, "./figure1_values.py"], B13),
                       ([PY, "./table1_inputs.py"], B13),
                       ([PY, "make_figure1_R2.py"], j("4_visualisation/figure1")),
                       (["rsvg-convert", "-d", "300", "-p", "300", "Figure_1_R2.svg", "-o", "Figure_1_R2.png"], j("4_visualisation/figure1")),
                       (["sips", "-s", "dpiWidth", "300", "-s", "dpiHeight", "300", "Figure_1_R2.png"], j("4_visualisation/figure1")),
                       (["rsvg-convert", "-d", "72", "-p", "72", "-f", "pdf", "Figure_1_R2.svg", "-o", "Figure_1_R2.pdf"], j("4_visualisation/figure1"))]),
    # Fig. S7: gain ranking from the package results_table, Gini of the models; HOCOMOCO v11 match per tile, TFs and
    # roles from subtype_verdicts_v11_h3.csv; then the v11 cross-check csv and the feature-importance clustering vs v11 agreement
    "s7": dict(backup=["4_visualisation/figures/FigS7/Figure_S7_motif_subtypes.p*", "4_visualisation/figures/FigS7/Figure_S7_source_data.csv", "4_visualisation/figures/FigS7/tf_ranking_h3.csv",
                       "4_visualisation/figures/FigS7/subtype_verdicts_v11_h3.csv", "4_visualisation/figures/FigS7/hocomoco11_check/*.csv",
                       "4_visualisation/figures/Fig3_S3/Figure_3D_*", "outputs/figures/Figure_3.p*"],
               cmds=[([PY, "build_Figure_S7_source_data_h3.py"], j("4_visualisation/figures", "FigS7")),
                     ([PY, "subtype_vs_hocomoco11.py"], j("4_visualisation/figures", "FigS7", "hocomoco11_check")),
                     # r and quality of the reference TFs the text names (ETV2, ETV7, BATF, FLI1), co-factor gaps
                     ([PY, "subtype_vs_named_v11.py"], j("4_visualisation/figures", "FigS7", "hocomoco11_check")),
                     (["Rscript", "Figure_S7_motif_subtypes.R"], j("4_visualisation/figures", "FigS7")),
                     # Fig. 3D from the Gini (same builder), then the vector composite of Fig. 3. Reads Sup. Table 2,
                     # so the stage runs after "sup"
                     ([PY, "Fig3_S3/build_Figure_3D_source_data_h3.py"], j("4_visualisation/figures")),
                     (["Rscript", "Fig3_S3/Figure_3D_logos.R"], j("4_visualisation/figures")),
                     ([PY, "Fig3_S3/assemble_figure3.py"], j("4_visualisation/figures"))]),
    # step-0 gate of the low-FPR metrics on b2_v2 and the best-mono baseline. With ARCHI_STAGE=H4 (inputs not in the repo) the 35 score npz of the
    # training-scaler run are not in the repository ($ARCHI_WORK/human_sym/H4_trainscale/scores); set
    # ARCHI_H4_SCORES to that folder (ARCHI_H3A_SCORES is still read, for the old name).
    "gate": dict(backup=["3_analysis/low_fpr_metrics/B2_gate_h3_*.csv"],
                 # ARCHI_STAGE=H5: ARCHI_H5_SCORES, the folder $ARCHI_WORK/H5_alltrain/human_sym/scores (unpacked scores/ folder of the Zenodo v2.0.0 archive scores_H5_human.tar.gz; no fallback to the previous run)
                 cmds=[([PY, "../low_fpr_metrics/b2_gate_h3.py", "--scores", (os.environ.get("ARCHI_H4_SCORES") or os.environ.get("ARCHI_H3A_SCORES")
                         or "ARCHI_H4_SCORES-not-set") if H.STAGE == "H4" else (os.environ.get("ARCHI_H5_SCORES") or "ARCHI_H5_SCORES-not-set"),
                         "--scores-origin", {"H4": "$ARCHI_WORK/human_sym/H4_trainscale/scores",
                                             "H5": "$ARCHI_WORK/H5_alltrain/human_sym/scores"}[H.STAGE]], B13)]),
    # random split (Fig. S12, Sup. Table 7 sheets); then rebuild Sup. Table 7 with the sup stage
    "b9": dict(backup=[], cmds=[([PY, "b9_figure_h3.py"], os.path.join(B13, "../random_split"))]),
    # PWM origin and the first-submission comparison tables (read Sup. Tables 1 and 2: run after "sup")
    "porig": dict(backup=[], cmds=[([PY, "./pwm_origin_h3.py"], B13)]),
    "v1": dict(backup=[], cmds=[([PY, "./v1_tables_h3.py"], B13)]),
    # NUMBERS: human 35-TF rows, then the full NUMBERS_R2_35TF.csv (reads every table above)
    "numbers": dict(backup=["outputs/numbers/NUMBERS_human_35TF*.csv", "outputs/numbers/NUMBERS_R2_35TF.csv"],
                    cmds=[([PY, "make_numbers.py"], j("3_analysis/integration")), ([PY, "./make_numbers_w2.py"], B13),
                          # rows whose TF count, sign or significance changed against the previous training run (flip_report_<stage>_vs_h4.csv)
                          ([PY, "./flip_report.py", "--ref", "h4"], B13)]),
    "copy": dict(backup=[], cmds=[([PY, "./copy_figures.py"], B13)]),
    # transfer drivers and headroom (Sup. Table 10, Fig. S11, Fig. 3C groups)
    "b5": dict(backup=["3_analysis/transfer_drivers/B5_*.csv", "3_analysis/transfer_drivers/B5_*.txt", "3_analysis/transfer_drivers/_provenance.yaml", "3_analysis/transfer_drivers/B5_figure.p*"],
               cmds=[([PY, "B5_01_drivers_and_ceiling.py", "--results-table", os.path.join(P, "results_table.csv"), "--mouse-transfer", H.IN["hm_mm"],
                       "--crossspecies", os.path.join(P, "analysis", "cross_species", "cross_species_table.csv")], j("3_analysis/transfer_drivers")),
                     ([PY, "B5_02_figure.py"], j("3_analysis/transfer_drivers"))]),
}


def run(stage):
    st = STAGES[stage]
    if stage == "gate":
        var = "ARCHI_H4_SCORES" if H.STAGE == "H4" else "ARCHI_H5_SCORES"
        d = os.environ.get(var) or (os.environ.get("ARCHI_H3A_SCORES") if H.STAGE == "H4" else None)
        if not d or not os.path.isdir(d) or not glob.glob(os.path.join(d, "*.npz")):
            raise SystemExit("CHECK FAILED stage gate: %s must name a folder with the score .npz files (for H5: the unpacked scores/ "
                             "folder of the Zenodo v2.0.0 archive scores_H5_human.tar.gz); got %r" % (var, d))
    H.check()
    backup_globs(*st["backup"])
    for i, c in enumerate(st["cmds"]):
        cmd, cwd, extra = c[0], c[1], (c[2] if len(c) > 2 else {})
        log = os.path.join(B13, "logs", "h3_%s_%02d.log" % (stage, i))
        print("+ (%s) %s  > %s" % (os.path.relpath(cwd, R2), " ".join(os.path.relpath(c, R2) if c.startswith(R2) else c for c in cmd),
                                   os.path.relpath(log, R2)), flush=True)
        with open(log, "w") as fh:
            p = subprocess.run(["nice", "-n", "15"] + cmd, cwd=cwd, env=dict(ENV, **extra), stdout=fh, stderr=subprocess.STDOUT)
        if p.returncode:
            print(open(log).read()[-3000:])
            raise SystemExit("CHECK FAILED stage %s step %d (rc %d), log %s" % (stage, i, p.returncode, log))


# "gate" is an optional extra stage (it needs the score files of the Zenodo v2.0.0 deposit) and is not in ALL.
# One order in which every consumer runs after its producers (xlsx before figure1_values, pwm_origin, v1, NUMBERS), so the
# input md5 recorded in each provenance block is the md5 of the file on disk at the end of the run
ALL = ["pkg", "figs", "b234", "b5", "b9", "sup", "s7", "porig", "v1", "fig1", "numbers", "copy"]

if __name__ == "__main__":
    stages = ALL if sys.argv[1:] == ["all"] else sys.argv[1:]
    # pre-flight: every registry input present, mirror files equal to their manifest md5, split labels in the GC file
    import preflight_inputs
    import zenodo_v1
    zenodo_v1.unpack(ENV["ARCHI_ZENODO_DIR"])
    if preflight_inputs.main([]) != 0:
        raise SystemExit("CHECK FAILED pre-flight (stage %s): see the table above" % H.STAGE)
    print("stage %s (inputs from %s)" % (H.STAGE, H.ROOT))
    for s in stages:
        run(s)
    print("done:", ", ".join(stages))
