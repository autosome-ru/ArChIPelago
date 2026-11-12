# outputs

Final deliverables of the local chain (`python3 3_analysis/integration/run_h3.py all`, see `3_analysis/README.md`).
Every file here is rewritten by a rerun; the tracked copies are those of the published run.

| folder | contents | producer |
|---|---|---|
| `figures/` | `Figure_2` to `Figure_4`, `Figure_S1` to `Figure_S12` (`.pdf`, `.png`) | panels by the producers listed in `4_visualisation/README.md`; copied by `3_analysis/integration/copy_figures.py`, Fig. 3 written by `4_visualisation/figures/Fig3_S3/assemble_figure3.py` |
| `supplementary_tables/` | `Sup_Table_1_ChIP-Seq_experiments.xlsx` to `Sup_Table_10_transfer_drivers_and_headroom.xlsx` | `3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py` (Sup. Tables 1-4 and 6), `make_sup_table_5_R2.py` (Sup. Table 5), `make_sup_tables_7_10.py` (Sup. Tables 7-10) |
| `numbers/` | `NUMBERS_R2_35TF.csv`: every number of the text with its source table, the value of the first submission, and the values of the two previous training runs (`h4 value`, `h3a value`); `NUMBERS_human_35TF.csv` and `NUMBERS_human_35TF_unchanged.csv`: the human rows; `flip_report_h5_vs_h4.csv`: rows whose TF count, sign or significance changed against the previous run | `3_analysis/integration/make_numbers.py`, `make_numbers_w2.py`, `flip_report.py --ref h4` |

Each folder holds the provenance log of its producers (`_provenance_h3.yaml`, `_provenance_p53drop.yaml`): inputs and
outputs with md5, git commit of the repository, run time.
