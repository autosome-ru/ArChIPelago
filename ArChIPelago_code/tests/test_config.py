"""Tests for archipielago.config module."""

import pytest
from pathlib import Path
from unittest.mock import patch
from archipielago.config import load_config, get_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _yaml_available():
    try:
        import yaml  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------

def test_load_config_valid(tmp_path):
    pytest.importorskip("yaml", reason="PyYAML not installed")
    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text(
        "paths:\n  data: /data/myproject\n"
        "tools:\n  sarus: /tools/sarus.jar\n"
    )
    cfg = load_config(cfg_file)
    assert "paths" in cfg
    assert "tools" in cfg
    assert cfg["paths"]["data"] == "/data/myproject"
    assert cfg["tools"]["sarus"] == "/tools/sarus.jar"


def test_load_config_missing_file(tmp_path):
    pytest.importorskip("yaml", reason="PyYAML not installed")
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "nonexistent.yml")


def test_load_config_missing_required_key(tmp_path):
    pytest.importorskip("yaml", reason="PyYAML not installed")
    cfg_file = tmp_path / "config.yml"
    # Only 'paths' present — 'tools' is missing
    cfg_file.write_text("paths:\n  data: /data/myproject\n")
    with pytest.raises(KeyError):
        load_config(cfg_file)


def test_load_config_missing_both_sections(tmp_path):
    pytest.importorskip("yaml", reason="PyYAML not installed")
    cfg_file = tmp_path / "config.yml"
    cfg_file.write_text("other_key: value\n")
    with pytest.raises(KeyError):
        load_config(cfg_file)


# ---------------------------------------------------------------------------
# get_path
# ---------------------------------------------------------------------------

def test_get_path_valid(tmp_path):
    cfg = {
        "paths": {"genome": "/data/hg38.fa"},
        "tools": {},
    }
    result = get_path(cfg, "genome")
    assert result == Path("/data/hg38.fa")


def test_get_path_placeholder_raises(tmp_path):
    cfg = {
        "paths": {"genome": "/path/to/hg38.fa"},
        "tools": {},
    }
    with pytest.raises(ValueError):
        get_path(cfg, "genome")


def test_get_path_missing_key_raises():
    cfg = {
        "paths": {},
        "tools": {},
    }
    with pytest.raises(KeyError):
        get_path(cfg, "genome")
