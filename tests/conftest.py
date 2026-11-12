import os
import shutil
import sys

sys.dont_write_bytecode = True      # never leave __pycache__ in the production code folders
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import pytest

import t1_helpers as T


def pytest_configure(config):
    config.addinivalue_line("markers", "realdata: needs the real data under $ARCHI_WORK/{human,mouse}_sym")
    config.addinivalue_line("markers", "r2: needs the local manuscript tables")
    config.addinivalue_line("markers", "bedtools: needs bedtools on PATH")


def pytest_collection_modifyitems(config, items):
    skip_g = pytest.mark.skip(reason="real data not present (set ARCHI_WORK)")
    skip_r = pytest.mark.skip(reason="manuscript tables not present (local tables only)")
    skip_b = pytest.mark.skip(reason="bedtools not on PATH")
    has_bt = shutil.which("bedtools") is not None
    for it in items:
        if "realdata" in it.keywords and not T.REAL_DATA:
            it.add_marker(skip_g)
        if "r2" in it.keywords and not (T.R2 and os.path.isdir(os.path.join(T.R2, "outputs/supplementary_tables"))):
            it.add_marker(skip_r)
        if "bedtools" in it.keywords and not has_bt:
            it.add_marker(skip_b)
