"""auROC and auPRC computed with the R package PRROC through rpy2.

Reproduces the scorers of the round-1 wrapper: PRROC::pr.curve(scores.class0 = scores, weights.class0 = labels,
dg.compute = FALSE)$auc.integral (auPRC, integral interpolation) and PRROC::roc.curve(scores.class0 = scores,
weights.class0 = labels)$auc (auROC). The result elements are looked up by name.
PRROC: Grau J, Grosse I, Keilwagen J. PRROC: computing and visualizing precision-recall and receiver operating
characteristic curves in R. Bioinformatics 31:2595 (2015).

Requirements: R with PRROC 1.4 (install.packages("PRROC")) and rpy2 3.5.x. rpy2 must load the R whose library holds
PRROC (set R_HOME). rpy2 is imported on the first metric call, so importing this module needs no R.
"""
import numpy as np

_PKG = None

_INSTALL = ('PRROC through rpy2 is required for auROC/auPRC: install R, run install.packages("PRROC") in it '
            '(PRROC 1.4) and pip install rpy2==3.5.17. rpy2 must load the R whose library holds PRROC; '
            'set R_HOME to that R.')


def _prroc():
    global _PKG
    if _PKG is None:
        try:
            from rpy2.robjects.packages import importr
            _PKG = importr("PRROC")
        except ImportError as e:
            raise ImportError(_INSTALL + " (%s)" % e) from e
    return _PKG


def _vectors(y, s):
    from rpy2.robjects.vectors import FloatVector
    y = np.ascontiguousarray(np.asarray(y, dtype=np.float64))
    s = np.ascontiguousarray(np.asarray(s, dtype=np.float64))
    if y.shape != s.shape or y.ndim != 1:
        raise ValueError("y and s must be 1-D arrays of equal length, got %s and %s" % (y.shape, s.shape))
    return FloatVector(s), FloatVector(y)


def _element(res, name):
    names = list(res.names)
    if name not in names:
        raise KeyError("PRROC result has no element %r (elements: %s)" % (name, names))
    return float(res[names.index(name)][0])


def auroc(y, s):
    pkg = _prroc()
    from rpy2.rinterface_lib import openrlib
    with openrlib.rlock:
        scores, labels = _vectors(y, s)
        return _element(pkg.roc_curve(scores, weights_class0=labels), "auc")


def auprc(y, s):
    pkg = _prroc()
    from rpy2.rinterface_lib import openrlib
    with openrlib.rlock:
        scores, labels = _vectors(y, s)
        return _element(pkg.pr_curve(scores, weights_class0=labels, dg_compute=False), "auc.integral")


def auroc_auprc(y, s):
    """(auROC, auPRC) as Python floats."""
    return auroc(y, s), auprc(y, s)
