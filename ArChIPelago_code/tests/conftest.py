"""Shared pytest fixtures for archipielago tests."""

import io
import textwrap
import numpy as np
import pandas as pd
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# FASTA fixtures
# ---------------------------------------------------------------------------

SIMPLE_FASTA = textwrap.dedent("""\
    >seq1
    ACGTACGTACGT
    >seq2
    TTTTAAAACCCC
    >seq3
    GGGGCCCCTTTT
""")

MULTILINE_FASTA = textwrap.dedent("""\
    >seq1
    ACGT
    ACGT
    >seq2
    TTTT
    AAAA
""")

EMPTY_FASTA = ""


@pytest.fixture
def simple_fasta_handle():
    return io.StringIO(SIMPLE_FASTA)


@pytest.fixture
def multiline_fasta_handle():
    return io.StringIO(MULTILINE_FASTA)


@pytest.fixture
def tmp_fasta(tmp_path):
    """Write SIMPLE_FASTA to a temp file and return its Path."""
    p = tmp_path / "test.fasta"
    p.write_text(SIMPLE_FASTA)
    return p


# ---------------------------------------------------------------------------
# Feature matrix fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def small_feature_matrix():
    """3 sequences × 5 features."""
    np.random.seed(0)
    X = pd.DataFrame(
        np.random.rand(6, 10),
        columns=[f"mono_pwm{i}" for i in range(5)] + [f"di_pwm{i}" for i in range(5)],
    )
    y = np.array([1, 1, 1, 0, 0, 0])
    return X, y


@pytest.fixture
def mock_sarus_scan_dir(tmp_path):
    """Create a fake scan directory with 3 mono and 2 di score files."""
    mono_dir = tmp_path / "mono"
    di_dir = tmp_path / "di"
    mono_dir.mkdir()
    di_dir.mkdir()

    # 4 sequences × 3 mono PWMs
    for i in range(3):
        scores = "\n".join(str(v) for v in [0.1 * (i + 1), 0.2, 0.3, 0.4])
        (mono_dir / f"pwm_{i}.txt").write_text(scores + "\n")

    # 4 sequences × 2 di PWMs
    for i in range(2):
        scores = "\n".join(str(v) for v in [0.5, 0.6, 0.7 * (i + 1), 0.8])
        (di_dir / f"dpwm_{i}.txt").write_text(scores + "\n")

    return tmp_path


# ---------------------------------------------------------------------------
# Peaks / BED fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def peaks_df():
    """Small peaks DataFrame with 'chrom' column."""
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr8", "chr8", "chr21"],
        "start": [100, 200, 300, 400, 500, 600, 700, 800],
        "end":   [400, 500, 600, 700, 800, 900, 1000, 1100],
        "name":  list(range(8)),
    })
