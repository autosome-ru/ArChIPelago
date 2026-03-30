#!/usr/bin/env python3
"""
End-to-end ArChIPelago pipeline simulation.

Uses real CTCF FASTA data from the legacy directory and synthetic PWMs
to test every step of the pipeline without Zenodo downloads.

Usage:
    python run_simulation.py
"""
import os, sys, tempfile, shutil
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
os.chdir(Path(__file__).parent)

from archipielago import io, scanning, training

def create_synthetic_pwms(work_dir, rng):
    """Create 3 mono + 2 di PWMs for testing."""
    PWM_MONO = work_dir / "pwm_mono"; PWM_MONO.mkdir()
    PWM_DI = work_dir / "pwm_di"; PWM_DI.mkdir()

    # CTCF consensus-like: CCGCGNGGNGGCAG
    ctcf_pwm = np.array([
        [0.05,0.80,0.10,0.05],[0.05,0.80,0.10,0.05],[0.05,0.05,0.80,0.10],
        [0.05,0.80,0.05,0.10],[0.05,0.05,0.80,0.10],[0.25,0.25,0.25,0.25],
        [0.05,0.05,0.80,0.10],[0.05,0.05,0.80,0.10],[0.25,0.25,0.25,0.25],
        [0.05,0.05,0.80,0.10],[0.05,0.05,0.80,0.10],[0.05,0.80,0.05,0.10],
        [0.80,0.05,0.10,0.05],[0.05,0.05,0.80,0.10],
    ])

    for name, matrix in [("CTCF_motif", ctcf_pwm),
        ("random_1", np.clip(np.full((10,4),0.25)+rng.normal(0,0.02,(10,4)),0.01,0.99)),
        ("random_2", np.clip(np.full((8,4),0.25)+rng.normal(0,0.02,(8,4)),0.01,0.99))]:
        matrix = matrix / matrix.sum(axis=1, keepdims=True)
        lo = np.log2(matrix / 0.25)
        with open(PWM_MONO / f"{name}.pwm", "w") as f:
            f.write(f">{name}\n")
            for row in lo:
                f.write("\t".join(f"{v:.4f}" for v in row) + "\n")

    for i in range(2):
        m = np.clip(np.full((8,16),1/16)+rng.normal(0,0.005,(8,16)),0.001,0.999)
        m = m / m.sum(axis=1, keepdims=True)
        lo = np.log2(m / (1/16))
        with open(PWM_DI / f"di_pwm_{i}.dpwm", "w") as f:
            f.write(f">di_pwm_{i}\n")
            for row in lo:
                f.write("\t".join(f"{v:.4f}" for v in row) + "\n")

    return PWM_MONO, PWM_DI


def main():
    SARUS_JAR = Path("../sarus/releases/sarus-2.2.3.jar")
    CTCF_DIR = Path("../ArChIPelago-TFBS-finder/ArChIPelago-TFBS-finder legacy code/CTCF_data")

    if not SARUS_JAR.exists():
        sys.exit(f"SARUS jar not found: {SARUS_JAR}")
    if not CTCF_DIR.exists():
        sys.exit(f"CTCF data not found: {CTCF_DIR}")

    rng = np.random.default_rng(42)
    N, N_test = 200, 100

    # --- Step 1: Load FASTA ---
    print("STEP 1: Loading CTCF FASTA sequences")
    train_pos = io.load_fasta(CTCF_DIR / "train_pos_id_HUMAN.fasta")
    test_pos  = io.load_fasta(CTCF_DIR / "test_pos_id_HUMAN.fasta")
    test_neg  = io.load_fasta(CTCF_DIR / "test_neg_id_HUMAN.fasta")
    with open(CTCF_DIR / "train_H_0_SLIM.fasta") as f:
        train_neg_all = list(io.fasta_iter(f))

    train_pos_sub = [train_pos[i] for i in rng.choice(len(train_pos), N, replace=False)]
    train_neg_sub = [train_neg_all[i] for i in rng.choice(len(train_neg_all), N, replace=False)]
    test_pos_sub = [test_pos[i] for i in rng.choice(len(test_pos), N_test, replace=False)]
    test_neg_sub = [test_neg[i] for i in rng.choice(len(test_neg), N_test, replace=False)]
    print(f"  Train: {N} pos + {N} neg, Test: {N_test} pos + {N_test} neg")

    # --- Step 2: Create PWMs ---
    print("\nSTEP 2: Creating synthetic PWMs")
    WORK_DIR = Path(tempfile.mkdtemp(prefix="archip_sim_"))
    PWM_MONO, PWM_DI = create_synthetic_pwms(WORK_DIR, rng)
    print(f"  3 mono + 2 di PWMs")

    try:
        # --- Step 3: SARUS scanning ---
        print("\nSTEP 3: SARUS scanning")
        train_records = train_pos_sub + train_neg_sub
        train_labels = np.array([1]*N + [0]*N)
        test_records = test_pos_sub + test_neg_sub
        test_labels = np.array([1]*N_test + [0]*N_test)

        io.save_fasta(train_records, WORK_DIR / "train.fasta")
        io.save_fasta(test_records, WORK_DIR / "test.fasta")

        SCAN_TRAIN = WORK_DIR / "scan_train"
        SCAN_TEST = WORK_DIR / "scan_test"
        for d in [SCAN_TRAIN/"mono", SCAN_TRAIN/"di", SCAN_TEST/"mono", SCAN_TEST/"di"]:
            d.mkdir(parents=True)

        for pwm_file in sorted(PWM_MONO.glob("*.pwm")):
            for fasta, scan_dir in [(WORK_DIR/"train.fasta", SCAN_TRAIN), (WORK_DIR/"test.fasta", SCAN_TEST)]:
                scanning.run_sarus(fasta, pwm_file, SARUS_JAR, scan_dir/"mono"/f"{pwm_file.stem}.txt", pwm_type="mono")
            s = scanning.load_sarus_scores(SCAN_TRAIN/"mono"/f"{pwm_file.stem}.txt")
            print(f"  mono {pwm_file.stem}: {len(s)} scores, mean={s.mean():.3f}")

        for pwm_file in sorted(PWM_DI.glob("*.dpwm")):
            for fasta, scan_dir in [(WORK_DIR/"train.fasta", SCAN_TRAIN), (WORK_DIR/"test.fasta", SCAN_TEST)]:
                scanning.run_sarus(fasta, pwm_file, SARUS_JAR, scan_dir/"di"/f"{pwm_file.stem}.txt", pwm_type="di")
            s = scanning.load_sarus_scores(SCAN_TRAIN/"di"/f"{pwm_file.stem}.txt")
            print(f"  di   {pwm_file.stem}: {len(s)} scores, mean={s.mean():.3f}")

        # --- Step 4: Feature matrix ---
        print("\nSTEP 4: Building feature matrix")
        X_train = scanning.build_feature_matrix(SCAN_TRAIN, mode="mono_di")
        X_test = scanning.build_feature_matrix(SCAN_TEST, mode="mono_di")
        print(f"  Train: {X_train.shape}, Test: {X_test.shape}")

        # --- Step 5: Train RF ---
        print("\nSTEP 5: Training RandomForest")
        model = training.train_rf(X_train, train_labels, n_estimators=50, max_depth=6, n_jobs=4, random_state=42)
        imp = pd.Series(model.feature_importances_, index=X_train.columns).nlargest(3)
        print(f"  {model.n_estimators} trees, top feature: {imp.index[0]} ({imp.iloc[0]:.4f})")

        # --- Step 6: Evaluate ---
        print("\nSTEP 6: Evaluation")
        results = training.evaluate_model(model, X_test, test_labels)
        print(f"  ROC-AUC: {results['roc_auc']:.4f}")
        print(f"  PR-AUC:  {results['pr_auc']:.4f}")

        # --- Step 7: Predictions ---
        print("\nSTEP 7: Predictions")
        y_proba = model.predict_proba(X_test)[:, 1]
        n_correct = ((y_proba > 0.5) == test_labels).sum()
        print(f"  Accuracy: {n_correct}/{len(test_labels)} ({100*n_correct/len(test_labels):.1f}%)")

        # --- Step 8: Cross-validation ---
        print("\nSTEP 8: 3-fold cross-validation")
        from sklearn.ensemble import RandomForestClassifier
        cv = training.cross_validate_model(
            RandomForestClassifier(n_estimators=50, max_depth=6, random_state=42),
            X_train, train_labels, n_splits=3)
        print(f"  ROC-AUC: {cv['roc_auc_mean']:.4f} +/- {cv['roc_auc_std']:.4f}")
        print(f"  PR-AUC:  {cv['pr_auc_mean']:.4f} +/- {cv['pr_auc_std']:.4f}")

        print("\n" + "=" * 50)
        print("ALL 8 STEPS PASSED SUCCESSFULLY")
        print("=" * 50)

    finally:
        shutil.rmtree(WORK_DIR, ignore_errors=True)


if __name__ == "__main__":
    main()
