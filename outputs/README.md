# outputs

Supplementary tables and numbers written by the local chain (`python3 3_analysis/integration/run_h3.py all`, see
`../3_analysis/README.md`). A rerun rewrites every file here. The figures are in `../4_visualisation/`, next to their
producers.

| folder | contents | producer |
|---|---|---|
| `supplementary_tables/` | `Sup_Table_1_ChIP-Seq_experiments.xlsx` to `Sup_Table_10_transfer_drivers_and_headroom.xlsx` | `../3_analysis/supplementary_tables/make_sup_tables_1_3_R2.py` (Sup. Tables 1-4 and 6), `../3_analysis/supplementary_tables/make_sup_table_5_R2.py` (Sup. Table 5), `../3_analysis/supplementary_tables/make_sup_tables_7_10.py` (Sup. Tables 7-10) |
| `numbers/` | `NUMBERS_R2_35TF.csv`: every number of the text with its source table; `NUMBERS_human_35TF.csv` and `NUMBERS_human_35TF_unchanged.csv`: the human rows; `flip_report_h5_vs_h4.csv`: rows whose TF count, sign or significance differ from the reference run | `../3_analysis/integration/make_numbers_w2.py`, `../3_analysis/integration/make_numbers.py`, `../3_analysis/integration/flip_report.py` (`--ref h4`) |

`supplementary_tables/` holds the provenance logs of its producers (`_provenance_h3.yaml`, `_provenance_p53drop.yaml`):
inputs and outputs with md5, git commit of the repository, run time.
