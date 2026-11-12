# tests

Checks of the local chain: each test recomputes a value of the tables, figures or
`outputs/numbers/NUMBERS_R2_35TF.csv` from its source table and compares. 54 tests in 10 files.

## Run

From the repository root, after `python3 3_analysis/integration/run_h3.py all` or on the tracked outputs:

```
ARCHI_STAGE=H5 python3 -m pytest -q tests
```

Tests marked `r2` are skipped when `outputs/supplementary_tables/` is absent, tests marked `bedtools` when bedtools is
not on `PATH`. `pdftotext` is taken from `PATH` or from the environment variable `PDFTOTEXT`. `test_prroc_metrics.py`
needs R with PRROC and rpy2 and errors without them. GitHub Actions runs the suite on every push
(`../.github/workflows/tests.yml`).

| file | what it checks |
|---|---|
| `test_numbers_trace.py` | every value of `outputs/numbers/NUMBERS_R2_35TF.csv` whose source is a table occurs in that table |
| `test_provenance.py` | every `_provenance_h3.yaml` parses, and the md5 values it records equal the files on disk |
| `test_builders_rerun.py` | the supplementary table builders, rerun in a temporary copy, reproduce the tracked tables |
| `test_r2_sup_fixes.py`, `test_round3_tables.py` | cells of the supplementary tables and the group thresholds of Fig. 3C and Fig. S11B, recomputed from the producers' source tables |
| `test_figS7_v11.py` | Fig. S7 labels against the HOCOMOCO v11 comparison tables in `4_visualisation/supplementary_figures/Figure_S7/hocomoco11_check/` |
| `test_pwm_origin.py` | `3_analysis/integration/pwm_origin_h3.csv`, recomputed from the PWM archives in `results/datasets/zenodo_v1.0.0/` |
| `test_inputs_h5.py` | the input registry (`3_analysis/integration/h3_inputs.py`), the pre-flight and the flip report |
| `test_h5_splits.py` | the chromosome sets of training and test |
| `test_prroc_metrics.py` | auROC and auPRC of `1_data_preparation/prroc_metrics.py` (PRROC through rpy2) against a NumPy reference |

Shared helpers are in `conftest.py` and `t1_helpers.py`.

## Red runs

Eight of the ten test files can be made to fail on purpose (all but `test_builders_rerun.py` and
`test_prroc_metrics.py`): the environment variable `T1_MUTATE=<name>` injects one fault into the data a test reads (the
names are in each test's docstring). For example
`ARCHI_STAGE=H5 T1_MUTATE=numbers_value_changed python3 -m pytest -q tests/test_numbers_trace.py` must fail.
