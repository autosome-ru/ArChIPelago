# Figure S10

Feature importance of the Random Forest models: importance concentration, PWM subtypes and ablation.

Producer: `Figure_S10.py`.

Inputs (paths from the repository root): `3_analysis/feature_importance/tables/B4_*.csv`.

Rebuilt by stage `b234` of `3_analysis/integration/run_h3.py`.

Panel C shows the ablation of the secondary subtype for the 29 TFs that have one (cut 0.8 for 20 TFs, 0.9 for 9; ten seeds per condition): change in auROC (upper strip) and auPRC (lower strip) after retraining without the secondary subtype and without the same number of random primary PWMs, against the full model of the same seed. The values drawn are in `Figure_S10C_plotdata.csv`.

Files: `Figure_S10.pdf`, `Figure_S10.png`, `Figure_S10.py`, `Figure_S10C_plotdata.csv`.
