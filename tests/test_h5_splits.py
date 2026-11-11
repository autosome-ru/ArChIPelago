"""Every autosome (human 1-22, mouse 1-19) is in exactly one of the two chromosome sets, and chrX, chrY and anything else
in neither, for both species: in the one definition (1_data_preparation/h5_splits.py), in the lists the build code
actually uses (build_positives_v5.SPECIES, on-disk order incl. the mouse swap of train and test files) and in the test
chromosomes of the loader (archi_data_v5.py).
Red: T1_MUTATE=h5_old_lists takes the lists from the earlier build code build_positives_sym.py
("11" "12" typo, human chr22 unlisted); expected red for both species.
"""
import importlib.util
import os
import sys
import types

import pytest

import t1_helpers as T

CODE5 = os.path.join(T.R2, "1_data_preparation")
sys.path.insert(0, CODE5)
import h5_splits as S  # noqa: E402

EXPECTED = {
    "HUMAN": {"train": [2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22], "test": [1, 8, 21]},
    "MOUSE": {"train": [2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], "test": [1, 8, 19]},
}


def _load(path, name):
    if "pybedtools" not in sys.modules:          # only the module-level lists are needed, not bedtools
        sys.modules["pybedtools"] = types.ModuleType("pybedtools")
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _build_lists(sp):
    """(logical train, logical test) as the build code uses them; file 'train' of MOUSE is the test set."""
    if T.mutate("h5_old_lists"):
        b = _load(os.path.join(T.R2, "B14_human_sym", "code", "build_positives_sym.py"), "bp_sym_old")
    else:
        b = _load(os.path.join(CODE5, "build_positives_v5.py"), "bp_v5")
    ftr, fte = b.SPECIES[sp]["TRAIN_CHR"], b.SPECIES[sp]["TEST_CHR"]
    return (fte, ftr) if sp == "MOUSE" else (ftr, fte)


@pytest.mark.parametrize("sp", ["HUMAN", "MOUSE"])
def test_definition_is_a_partition_of_the_autosomes(sp):
    tr, te = S.LOGICAL[sp]["train"], S.LOGICAL[sp]["test"]
    if T.mutate("h5_old_lists"):
        tr, te = _build_lists(sp)
    bad = S.partition_problems(sp, tr, te)
    print(sp, "train", tr, "test", te, "problems", bad)
    assert not bad, bad
    assert sorted(int(c[3:]) for c in tr) == EXPECTED[sp]["train"]
    assert sorted(int(c[3:]) for c in te) == EXPECTED[sp]["test"]
    assert not set(S.SEX) & (set(tr) | set(te))


@pytest.mark.parametrize("sp", ["HUMAN", "MOUSE"])
def test_build_code_uses_the_partition(sp):
    tr, te = _build_lists(sp)
    bad = S.partition_problems(sp, tr, te)
    print(sp, "build code: train", tr, "test", te, "problems", bad)
    assert not bad, bad
    assert sorted(tr) == sorted(S.LOGICAL[sp]["train"]) and sorted(te) == sorted(S.LOGICAL[sp]["test"])


def test_loader_test_chromosomes():
    ad = _load(os.path.join(CODE5, "archi_data_v5.py"), "ad_v5")
    for sp in ("HUMAN", "MOUSE"):
        assert sorted(ad.TEST_CHROMS[sp]) == sorted(S.LOGICAL[sp]["test"])
    assert ad.HUMAN_REL.endswith("/archi_work/human_sym_h5/TF-ML") and ad.MOUSE_REL.endswith("/archi_work/mouse_sym_h5/TF-ML")
    assert ad._file_split("MOUSE", "train") == "control" and ad._file_split("HUMAN", "train") == "train"   # mouse swap kept


def test_checker_flags_each_kind_of_error():
    good_tr, good_te = S.LOGICAL["HUMAN"]["train"], S.LOGICAL["HUMAN"]["test"]
    assert S.partition_problems("HUMAN", good_tr, good_te) == []
    assert "in neither set ['chr22']" in " ".join(S.partition_problems("HUMAN", [c for c in good_tr if c != "chr22"], good_te))
    assert "not an autosome ['chrX']" in " ".join(S.partition_problems("HUMAN", good_tr + ["chrX"], good_te))
    assert "in both sets ['chr21']" in " ".join(S.partition_problems("HUMAN", good_tr + ["chr21"], good_te))
    assert "duplicated" in " ".join(S.partition_problems("MOUSE", S.LOGICAL["MOUSE"]["train"] * 2, S.LOGICAL["MOUSE"]["test"]))
    assert "not an autosome ['chr20']" in " ".join(S.partition_problems("MOUSE", S.LOGICAL["MOUSE"]["train"] + ["chr20"], S.LOGICAL["MOUSE"]["test"]))
