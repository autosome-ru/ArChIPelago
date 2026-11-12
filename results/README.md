# results

Tracked inputs of the local stages (`../3_analysis/`, `../4_visualisation/`). The local chain writes nothing here except
the unpacked PWM folders in `datasets/zenodo_v1.0.0/` (not tracked).

| folder | contents | written by or source |
|---|---|---|
| `server/` | result tables of the compute stages: `human/` (human-trained models: evaluation tables, `hm_mm_results.csv` with RF seeds 0 and 1 and the PWM baselines, operational metrics, saturation, low FPR metrics `b2_v2/`, top-N peaks `b3/`, feature importance `B4_importance/`, random split `human/H3b/b9/`, Fig. 4 and S4 refits `fig4_refit_trainscale/`, runtime benchmark `runtime/`, negatives checks `neg_v3/`), `mouse/` (mouse-trained evaluation and saturation), `human_data/`, `mouse_data/` (positive counts per TF, data checks) | `../1_data_preparation/`, `../2_training/`, `../3_analysis/server/`; `server/MIRROR_MANIFEST_md5.txt` holds the md5 of every file and is checked by `../3_analysis/integration/preflight_inputs.py` |
| `server_h4_test_sets/` | per-dataset peak widths and the join of the test windows with the peak files (`human/b3/`, `mouse/b3/`); the test sets of this reference run are identical to those of the final run | `server_h4_test_sets/MIRROR_MANIFEST_md5.txt` as above |
| `datasets/peak_ids/` | `human_used_peak_ids.csv` (362 human peak sets, 361 without P53), `mouse_sym_peak_ids.csv` (276 mouse peak sets of the 35 TFs) | dataset selection rule of the paper |
| `datasets/zenodo_v1.0.0/` | the four PWM archives `PWMs_{mono,di}_{HUMAN,MOUSE}.tar.gz` and `md5sums.txt` | Zenodo record 10.5281/zenodo.23036475; unpacked by `../3_analysis/integration/zenodo_v1.py` after the md5 check |
| `datasets/hocomoco11/` | HOCOMOCO v11 full HUMAN mono PCMs (769 models) and their annotation | HOCOMOCO v11 (Kulakovskiy et al. 2018) |
| `datasets/split_leakage/` | overlap of the top 1,000 peaks between datasets under a random split (masked and summary), and the top-1,000 check | read by `../3_analysis/supplementary_tables/make_sup_tables_7_10.py` and `../3_analysis/integration/make_numbers_w2.py` |
| `datasets/repeat_mask/repeat_mask_counts.csv` | windows removed by the RepeatMasker mask (human, mouse) | |
| `datasets/mouse_rebuild/mouse_rebuild_cpu_time.csv` | CPU time of the mouse data build | |
| `datasets/top_peaks/b3_missing_datasets.csv` | 19 human peak sets without a peak file | |
| `datasets/round1_fig_s4/`, `datasets/round1_b2_gate/` | earlier Fig. S4 numbers and low FPR gate summary, read for the comparison columns | |
| `datasets/text_tokens/manuscript_text_tokens.txt` | text lines with `[PENDING_H3B:<name>]` tokens; `../3_analysis/integration/fig4_tokens_h3b.py` computes one value per token | |
| `round1_package/` | earlier released tables read for the comparison columns of the numbers and tables: `Manuscript_analysis/` (results table, Sup. Table 6, Fig. 4 and S5 source data, motif similarity), `Resubmission/Sup_Tables/` (Sup. Tables 1, 2 and 4), `analysis_mouse_swap/hm_mm_results.csv` | |
| `previous_values/` | frozen copies of earlier versions of tables that the comparison columns of `../outputs/numbers/NUMBERS_R2_35TF.csv` read (through `../3_analysis/integration/h3_old.py`); one folder per git reference; `previous_values/SOURCES.csv` lists the md5 of each file | |

The analysis package keeps its own inputs in `../3_analysis/package/analysis/inputs/` (GTRD experiment table and
metadata, TF table, single-PWM features); `../3_analysis/integration/make_pkg_copy.py` copies the server tables there.
