"""Tests for archipielago.scoring module."""

import pytest
import numpy as np
from archipielago.scoring import (
    ConstantScorer,
    SklearnROCAUC,
    SklearnPRAUC,
    PRROC_ROCAUC,
    PRROC_PRAUC,
    ScorerInfo,
)


# ---------------------------------------------------------------------------
# ConstantScorer
# ---------------------------------------------------------------------------

def test_constant_scorer_returns_constant():
    scorer = ConstantScorer("const", 0.75)
    assert scorer.score() == pytest.approx(0.75)
    # Should ignore any positional/keyword arguments
    assert scorer.score([1, 2, 3], [0, 1, 0]) == pytest.approx(0.75)


def test_constant_scorer_name():
    scorer = ConstantScorer("my_const", 0.5)
    assert scorer.name == "my_const"
    assert scorer.const == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# SklearnROCAUC
# ---------------------------------------------------------------------------

def test_sklearn_rocauc_perfect():
    scorer = SklearnROCAUC("roc")
    y_score = [0.9, 0.8, 0.7, 0.2, 0.1, 0.0]
    y_real  = [1,   1,   1,   0,   0,   0  ]
    assert scorer.score(y_score, y_real) == pytest.approx(1.0)


def test_sklearn_rocauc_random():
    """Random predictions should give ~0.5 AUC (may vary; just check range)."""
    rng = np.random.default_rng(42)
    y_score = rng.random(200).tolist()
    y_real  = [1] * 100 + [0] * 100
    scorer = SklearnROCAUC("roc")
    result = scorer.score(y_score, y_real)
    assert 0.0 <= result <= 1.0


def test_sklearn_rocauc_worst():
    """Inverted predictions give AUC == 0.0."""
    scorer = SklearnROCAUC("roc")
    y_score = [0.0, 0.1, 0.2, 0.8, 0.9, 1.0]
    y_real  = [1,   1,   1,   0,   0,   0  ]
    assert scorer.score(y_score, y_real) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# SklearnPRAUC
# ---------------------------------------------------------------------------

def test_sklearn_prauc_perfect():
    scorer = SklearnPRAUC("pr")
    y_score = [0.9, 0.8, 0.1, 0.0]
    y_real  = [1,   1,   0,   0  ]
    assert scorer.score(y_score, y_real) == pytest.approx(1.0)


def test_sklearn_prauc_in_range():
    rng = np.random.default_rng(0)
    y_score = rng.random(100).tolist()
    y_real  = [1] * 50 + [0] * 50
    scorer = SklearnPRAUC("pr")
    result = scorer.score(y_score, y_real)
    assert 0.0 <= result <= 1.0


# ---------------------------------------------------------------------------
# PRROC_ROCAUC — dataclass instantiation only (rpy2 not required)
# ---------------------------------------------------------------------------

def test_prroc_rocauc_is_dataclass_instantiable():
    """BUGFIX H10: PRROC_ROCAUC was missing @dataclass — raised TypeError on instantiation."""
    scorer = PRROC_ROCAUC("rocauc")
    assert scorer.name == "rocauc"


def test_prroc_prauc_instantiable():
    scorer = PRROC_PRAUC("prauc", "integral")
    assert scorer.name == "prauc"
    assert scorer.type == "integral"


# ---------------------------------------------------------------------------
# ScorerInfo factory
# ---------------------------------------------------------------------------

def test_scorer_info_alias_defaults_to_name():
    """BUGFIX H9: __post_init__ was __attrs_post_init__ (never called for @dataclass)."""
    info = ScorerInfo(name="scikit_rocauc")
    assert info.alias == "scikit_rocauc"


def test_scorer_info_explicit_alias():
    info = ScorerInfo(name="scikit_prauc", alias="my_pr")
    assert info.alias == "my_pr"


def test_scorer_info_make_sklearn_roc():
    scorer = ScorerInfo(name="scikit_rocauc", alias="roc").make()
    assert isinstance(scorer, SklearnROCAUC)
    assert scorer.name == "roc"


def test_scorer_info_make_sklearn_pr():
    scorer = ScorerInfo(name="scikit_prauc", alias="pr").make()
    assert isinstance(scorer, SklearnPRAUC)


def test_scorer_info_make_constant():
    scorer = ScorerInfo(name="constant_scorer", alias="c", params={"cons": 0.42}).make()
    assert isinstance(scorer, ConstantScorer)
    assert scorer.score() == pytest.approx(0.42)


def test_scorer_info_to_dict_roundtrip():
    info = ScorerInfo(name="scikit_rocauc", alias="roc", params={"k": "v"})
    d = info.to_dict()
    info2 = ScorerInfo.from_dict(d)
    assert info2.name == info.name
    assert info2.alias == info.alias
    assert info2.params == info.params


# ---------------------------------------------------------------------------
# PRROC tests — skipped if rpy2 unavailable
# ---------------------------------------------------------------------------

rpy2 = pytest.importorskip("rpy2", reason="rpy2 not installed — skipping PRROC tests")


def test_prroc_rocauc_score():
    scorer = PRROC_ROCAUC("prroc_roc")
    y_score = [0.9, 0.8, 0.7, 0.2, 0.1, 0.0]
    y_real  = [1,   1,   1,   0,   0,   0  ]
    result = scorer.score(y_score, y_real)
    assert pytest.approx(result, abs=0.01) == 1.0


def test_prroc_prauc_integral_score():
    scorer = PRROC_PRAUC("prroc_pr", "integral")
    y_score = [0.9, 0.8, 0.1, 0.0]
    y_real  = [1,   1,   0,   0  ]
    result = scorer.score(y_score, y_real)
    assert 0.0 <= result <= 1.0
