"""Backward-compatibility shim — all scorer classes now live in archipielago.scoring."""
from archipielago.scoring import (  # noqa: F401
    Scorer, ConstantScorer, BinaryScorer, SklearnScorer,
    SklearnROCAUC, SklearnPRAUC, PRROCScorer,
    PRROC_PRAUC, PRROC_ROCAUC, ScorerInfo, import_PRROC,
)
    
    def to_dict(self) -> dict:
        dt = {}
        dt['name'] = self.name
        dt['alias'] = self.alias
        dt['params'] = self.params
        return dt