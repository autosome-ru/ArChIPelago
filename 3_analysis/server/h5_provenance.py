"""h5_provenance.py - writes _provenance.yaml: chromosome lists, md5 of scripts, inputs (large sets as aggregate
md5 = md5 of the sorted 'md5  path' lines, in _provenance_filelists/<name>.md5.txt) and outputs. Refuses to overwrite.
usage: ARCHI_WORK=<root> ARCHI_MASK_HUMAN=<track_out_HUMAN.bed> ARCHI_GIT_HASH=<hash> python h5_provenance.py"""
import glob, hashlib, os, sys, time, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")); import archi_env as AE
RR = AE.work(); CODE = AE.root()
G, GM = f"{RR}/H5_alltrain/human_sym", f"{RR}/H5_alltrain/mouse_sym"
DH, DM = f"{RR}/human_sym_h5", f"{RR}/mouse_sym_h5"
sys.path.insert(0, AE.loader_dir())
import h5_splits as S
out = f"{G}/_provenance.yaml"
if os.path.exists(out):
    raise SystemExit("exists: " + out)
os.makedirs(f"{G}/_provenance_filelists", exist_ok=False)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def one(p):
    return "%s: {md5: %s, bytes: %d}" % (p.replace(RR, "$ARCHI_WORK"), md5(p), os.path.getsize(p))


def agg(name, files):
    files = sorted(files)
    txt = "".join("%s  %s\n" % (md5(f), f.replace(RR, "$ARCHI_WORK")) for f in files)
    open(f"{G}/_provenance_filelists/{name}.md5.txt", "x").write(txt)
    return "%s: {n_files: %d, aggregate_md5: %s, list: _provenance_filelists/%s.md5.txt}" % (name, len(files), hashlib.md5(txt.encode()).hexdigest(), name)


ver = subprocess.run([sys.executable, "-c", "import sklearn, xgboost, numpy, pandas, scipy; print(sklearn.__version__, xgboost.__version__, numpy.__version__, pandas.__version__, scipy.__version__)"],
                     capture_output=True, text=True).stdout.strip()
L = ["task: H5 all autosomes in training or test (the revision round 2, D16/D18; training scaler D10/D10a), written %s" % time.strftime("%Y-%m-%d %H:%M"),
     "python: %s (sklearn xgboost numpy pandas scipy = %s); <= 32 threads (RF stage), <= 96 threads (Slim queue); nice 10" % (sys.version.split()[0], ver),
     "git_hash_repo: %s" % os.environ.get("ARCHI_GIT_HASH", "not given"),
     "chromosomes (logical sets; mouse on-disk file 'train' = test set, D11):"]
for sp in ("HUMAN", "MOUSE"):
    L.append("  %s: {train: [%s], test: [%s], excluded: [chrX, chrY, other contigs]}" % (sp, ", ".join(S.LOGICAL[sp]["train"]), ", ".join(S.LOGICAL[sp]["test"])))
L += ["parameters: {scale: train, cap: '10,000 per split, pandas sample random_state=1 (build); loader positives <= 10,000 random_state=0', "
      "RF: 'max_depth 6, max_samples 0.8, 100 trees, seeds 0/1 (hm_mm)', auPRC: PRROC integral, TFs: 35, Slim: 'm=0,1,-5 Starts=20 threads=4', "
      "diChIPMunk: 'run_dichiphorde8.rb 15:15 filter y 1.0 200 20 1 4 random auto flat (stochastic)'}", "scripts:"]
scripts = [p for d in ("1_data_preparation", "2_training", "3_analysis/server") for p in glob.glob(os.path.join(CODE, d, "**", "*"), recursive=True)
           if os.path.isfile(p) and p.endswith((".py", ".sh", ".rb", ".txt")) and "__pycache__" not in p]
L += ["  " + one(p) for p in sorted(set(scripts))]
L += ["inputs:"]
L += ["  " + one(p) for p in (f"{RR}/human_sym/inputs/human_sym_peak_ids.csv", f"{RR}/mouse_sym/inputs/mouse_sym_peak_ids.csv",
                              AE.var("ARCHI_MASK_HUMAN", "track_out_HUMAN.bed"), f"{RR}/mouse_sym/ref/track_out_MOUSE_mm10.bed",
                              f"{DH}/TF-ML/out_tab_HUMAN_10000_train.tab", f"{DH}/TF-ML/out_tab_HUMAN_10000_control.tab",
                              f"{DM}/TF-ML/out_tab_MOUSE_10000_train.tab", f"{DM}/TF-ML/out_tab_MOUSE_10000_control.tab",
                              f"{DH}/TF-ML/per_tf_positive_counts.csv", f"{DM}/TF-ML/per_tf_positive_counts.csv")]
L += ["  " + agg("neg_ids_v3_human_mouse_h5", glob.glob(f"{DH}/TF-ML/outputdir/*/*_neg_id_out_*_v3.csv") + glob.glob(f"{DM}/TF-ML/outputdir/*/*_neg_id_out_*_v3.csv")),
      "  " + agg("feature_tables_human_sym_h5", glob.glob(f"{DH}/TF-ML/outputdir/*/*_seq_*_pwm_*/*_cut.tab")),
      "  " + agg("feature_tables_mouse_sym_h5", glob.glob(f"{DM}/TF-ML/outputdir/*/*_seq_*_pwm_*/*_cut.tab")),
      "  " + agg("slim_inputs", glob.glob(f"{G}/slim/inputs/*/*")),
      "  " + agg("slim_dichipmunk_models", [f for f in glob.glob(f"{G}/slim/models_server/*/**/*", recursive=True) if os.path.isfile(f)]),
      "  " + agg("slim_munk_scans", glob.glob(f"{G}/slim/scans/*/*/model1_predictions.txt") + glob.glob(f"{G}/slim/scans/*/munk_*.tab"))]
L += ["outputs:"]
outs = [p for p in glob.glob(f"{G}/**/*.csv", recursive=True) + glob.glob(f"{GM}/**/*.csv", recursive=True) if "/_provenance_filelists/" not in p]
L += ["  " + agg("output_tables_csv", outs),
      "  " + agg("hm_mm_rows_json", glob.glob(f"{G}/rows/*.json")), "  " + agg("hm_mm_scores_npz", glob.glob(f"{G}/scores/*.npz")),
      "  " + agg("models_sav_json", glob.glob(f"{G}/models/*")),
      "  " + agg("eval_rows_json", glob.glob(f"{G}/eval_models/rows_HM/*.json") + glob.glob(f"{GM}/eval_models/rows_MM/*.json")),
      "  " + agg("fig4_rows_json", glob.glob(f"{G}/fig4_refit_trainscale/rows/*.json"))]
for p in (f"{G}/hm_mm_results.csv", f"{G}/eval_models/HM/evaluation_table.csv", f"{GM}/eval_models/MM/evaluation_table.csv",
          f"{G}/operational/operational_metrics.csv", f"{G}/b2_v2/B2_metrics_wide.csv", f"{G}/b2_v2/B2_per_chromosome.csv",
          f"{G}/b3/b3_topN_results.csv", f"{G}/saturation_HUMAN/saturation_results.csv", f"{GM}/saturation_MOUSE/saturation_results.csv",
          f"{G}/b9/b9_summary.csv", f"{G}/fig4_refit_trainscale/Figure_4_source_data.csv", f"{G}/fig4_refit_trainscale/Figure_S4_source_data.csv",
          f"{G}/fig4_refit_trainscale/Figure_4_and_S4_numbers.csv", f"{G}/checks/check_h5_post.csv", f"{G}/checks/headline_h5_vs_h4.csv",
          f"{G}/checks/flips_h5_vs_h4.csv"):
    L.append("  " + (one(p) if os.path.exists(p) else "%s: MISSING" % p.replace(RR, "$ARCHI_WORK")))
open(out, "x").write("\n".join(L) + "\n")
print(open(out).read()[:3000])
