# tests

Checks of the local chain: each test recomputes a value of the tables, figures or `outputs/numbers/NUMBERS_R2_35TF.csv`
from its source table and compares. `server_checks/` holds the checks that run on the
compute server.

## Run

From the repository root, after `3_analysis/integration/run_h3.py all` (or on the tracked outputs):

```
ARCHI_STAGE=H5 python3 -m pytest -q tests --ignore=tests/server_checks
```

Tests marked `r2` are skipped when `outputs/supplementary_tables/` is absent. `pdftotext` is taken from `PATH` or from the environment
variable `PDFTOTEXT`.

| file | what it checks |
|---|---|
| `test_numbers_trace.py` | every value of `NUMBERS_R2_35TF.csv` whose source is a table occurs in that table |
| `test_provenance.py` | every `_provenance_h3.yaml` parses, and the md5 values it records equal the files on disk |
| `test_builders_rerun.py` | the supplementary table builders, rerun in a temporary copy, reproduce the tracked tables |
| `test_r2_sup_fixes.py`, `test_round3_tables.py` | changed cells of the supplementary tables and the group thresholds of Fig. 3C and Fig. S11B, recomputed from the producers' source tables |
| `test_figS7_v11.py` | Fig. S7 labels against the HOCOMOCO v11 comparison tables in `4_visualisation/figures/FigS7/hocomoco11_check/` |
| `test_pwm_origin.py` | `3_analysis/integration/pwm_origin_h3.csv`, recomputed from the PWM archives in `results/datasets/zenodo_v1.0.0/` |
| `test_inputs_h5.py` | the input registry (`h3_inputs.py`), the pre-flight and the flip report |
| `test_h5_splits.py` | the chromosome sets of training and test |

## Red runs

Every test can be made to fail on purpose: the environment variable `T1_MUTATE=<name>` injects one fault into the data a
test reads (the names are in each test's docstring). For example
`ARCHI_STAGE=H5 T1_MUTATE=numbers_value_changed python3 -m pytest -q tests/test_numbers_trace.py` must fail.

## Result

52 tests, all pass on the tracked outputs. `test_inputs_h5.py` checks the pre-flight on the inputs of this repository only;
the inputs of the previous training run are not shipped.
