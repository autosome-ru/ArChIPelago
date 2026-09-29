"""Configuration loading for ArChIPelago pipelines."""

import os
from pathlib import Path

try:
    import yaml
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False


def load_config(path=None):
    """Load ArChIPelago config.yml.

    Parameters
    ----------
    path : str or Path, optional
        Path to config.yml. Defaults to config.yml in the same directory as
        this package (ArChIPelago_code/config.yml).

    Returns
    -------
    dict
        Configuration with the sections 'paths' and 'tools'.

    Raises
    ------
    ImportError
        If PyYAML is not installed.
    FileNotFoundError
        If the config file does not exist.
    KeyError
        If a required top-level key is missing.
    """
    if not _HAS_YAML:
        raise ImportError(
            "PyYAML is required to load config.yml. "
            "Install it with: conda install pyyaml"
        )

    if path is None:
        path = Path(__file__).parent.parent / "config.yml"

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Config file not found: {path}\n"
            "Copy config.yml.template to config.yml and fill in your paths."
        )

    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    for required in ("paths", "tools"):
        if required not in cfg:
            raise KeyError(
                f"config.yml is missing required section '{required}'. "
                f"Check {path} against config.yml.template."
            )

    return cfg


def get_path(cfg, key):
    """Return a path from the config; raise an error if it is missing or a placeholder.

    Parameters
    ----------
    cfg : dict
        Config dict returned by load_config().
    key : str
        Key in cfg['paths'].

    Returns
    -------
    Path
    """
    value = cfg["paths"].get(key)
    if value is None:
        raise KeyError(
            f"config.yml is missing paths.{key}. "
            "Set this path in config.yml."
        )
    if "/path/to/" in str(value):
        raise ValueError(
            f"config.yml paths.{key} is still a placeholder: {value!r}. "
            "Set the path on your machine."
        )
    return Path(value)
