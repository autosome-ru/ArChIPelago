"""ArChIPelago Python package.

Provides reusable functions for the ArChIPelago transcription factor binding
site prediction pipeline. Import individual submodules as needed:

    from archipielago import io, scanning, training, scoring
    from archipielago.config import load_config

Submodules
----------
config   : Load and validate config.yml
io       : FASTA/BED reading, writing, and sequence extraction
scanning : SARUS wrapper and feature matrix construction
training : RandomForest training and evaluation helpers
scoring  : Scorer classes (ROC-AUC, PR-AUC via sklearn and PRROC)
"""

from . import config, io, scanning, training, scoring
from .config import load_config

__all__ = ["config", "io", "scanning", "training", "scoring", "load_config"]
