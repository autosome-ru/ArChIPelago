# results

Tracked inputs of the local analysis (`3_analysis/`, `4_visualisation/`). Nothing here is written by the local chain,
except the unpacked PWM folders in `datasets/zenodo_v1.0.0/` (not tracked).

| folder | contents | source |
|---|---|---|
| `server/` | result tables of the compute stages that the local stages read: `human/` (human-trained models: evaluation tables, `hm_mm_results.csv` with RF seeds 0 and 1 and the PWM baselines, operational metrics, saturation, low FPR metrics `b2_v2/`, top-N peaks `b3/`, feature importance `B4_importance/`, random split `H3b/b9/`, Fig. 4 / S4 refits `fig4_refit_trainscale/`, runtime benchmark `runtime/`, negatives checks `neg_v3/`, the 35 score files `scores/*.npz`), `mouse/` (mouse-trained evaluation and saturation), `human_data/`, `mouse_data/` (positive counts per TF and the data checks) | written by `1_data_preparation/`, `2_training/` and `3_analysis/server/`; `MIRROR_MANIFEST_md5.txt` holds the md5 of every file (folder names as in this repository) and is checked by `3_analysis/integration/preflight_inputs.py` |
| `server_h4_test_sets/` | per-dataset peak widths and the join of the test windows with the peak files (`human/b3/`, `mouse/b3/`), and `human/eval_models/HM/evaluation_table.csv` of the previous training run | the test sets of the final run are byte-identical to those of the previous run, so these test-set tables were not recomputed; `MIRROR_MANIFEST_md5.txt` as above |
| `datasets/peak_ids/` | `human_used_peak_ids.csv` (362 human peak sets, 361 without P53), `mouse_sym_peak_ids.csv` (276 mouse peak sets of the 35 TFs) | dataset selection rule of the paper |
| `datasets/zenodo_v1.0.0/` | the four PWM archives of the Zenodo deposit v1.0.0 (`PWMs_{mono,di}_{HUMAN,MOUSE}.tar.gz`) and its `md5sums.txt` | doi 10.5281/zenodo.23036475; unpacked at run time by `3_analysis/integration/zenodo_v1.py` after the md5 check |
| `datasets/hocomoco11/` | HOCOMOCO v11 full HUMAN mono PCMs (769 models) and annotation | HOCOMOCO v11 (Kulakovskiy et al. 2018) |
| `datasets/split_leakage/` | overlap of the top 1,000 peaks between datasets under a random split (masked and summary), and the top-1,000 check | random split check of the first revision |
| `datasets/repeat_mask/repeat_mask_counts.csv` | windows removed by the RepeatMasker mask (human, mouse) | parsed from the description of that check |
| `datasets/mouse_rebuild/mouse_rebuild_cpu_time.csv` | CPU time of the mouse data rebuild | parsed from the run record |
| `datasets/top_peaks/b3_missing_datasets.csv` | 19 human peak sets without a peak file | |
| `datasets/round1_fig_s4/`, `datasets/round1_b2_gate/` | Fig. S4 numbers and the low FPR gate of the first submission (old values) | |
| `datasets/text_tokens/manuscript_text_tokens.txt` | the lines of the manuscript text producer that hold the `[PENDING_H3B:<name>]` tokens; `3_analysis/integration/fig4_tokens_h3b.py` computes one value per token | manuscript text producer before the values were filled |
| `round1_package/` | inputs of the first submission that the old-value columns read: `Manuscript_analysis/` (released analysis package: results table, Sup. Table 6, Fig. 4 and S5 source data, motif similarity), `Resubmission/Sup_Tables/` (Sup. Tables 1, 2, 4 of the first submission), `analysis_mouse_swap/hm_mm_results.csv`, `Figure_1_original.png` | first submission |
| `previous_values/` | frozen copies of the files that the old-value columns of `outputs/numbers/NUMBERS_R2_35TF.csv` read from the git history of the manuscript repository; one folder per git commit or tag, file names as in that commit; `SOURCES.csv` lists the original path and md5 of each | exact bytes of `git show <ref>:<path>` |

The analysis package keeps its own inputs inside `3_analysis/package/analysis/inputs/` (GTRD experiment table and metadata,
TF table, single-PWM features); the final server tables are copied there by `3_analysis/integration/make_pkg_copy.py`.
