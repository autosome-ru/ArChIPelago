"""archi_style.py - matplotlib twin of archi_style.R (same hex values; check_palette.py compares them).
Series: best single monoPWM red (triangle), best single diPWM orange (triangle), ArChIPelago models circles in blue /
green / light blue / pink / indigo. Opaque marks, black text, Arial, no MathText (Unicode symbols only)."""
import matplotlib

MODEL_LABELS = ["Single best monoPWM", "Single best diPWM", "Logistic Regression", "Random Forest (RF)",
                "Gradient Boosting (XGBoost)", "Bagging (Gradient Boosting)", "Bagging (Logistic Regression)"]
PALETTE_ORDER = ["monoPWM", "diPWM", "LogReg", "RF", "XGB", "BagGB", "BagLR"]
PALETTE = dict(zip(PALETTE_ORDER, ["#D62728", "#FF7F0E", "#0072B2", "#009E73", "#56B4E9", "#CC79A7", "#332288"]))
MODEL_COL = dict(zip(MODEL_LABELS, PALETTE.values()))
MODEL_MARKER = dict(zip(MODEL_LABELS, ["^", "^", "o", "o", "o", "o", "o"]))
RF_COL = PALETTE["RF"]
# PWM feature sets of one model (Random Forest): the colours of best monoPWM (red), best diPWM (orange) and RF (green), black outline
PWMSET_LABELS = {"mono": "monoPWMs", "di": "diPWMs", "mono+di": "monoPWMs+diPWMs"}
PWMSET_COL = {"mono": PALETTE["monoPWM"], "di": PALETTE["diPWM"], "mono+di": PALETTE["RF"]}
GREY_PTS = "#BDBDBD"


def apply_rc():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 11,
        "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
        "axes.edgecolor": "black", "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0, "lines.linewidth": 1.6, "axes.titlecolor": "black",
        "legend.frameon": False, "axes.unicode_minus": False})


def save_both(fig, path_noext):
    fig.savefig(path_noext + ".pdf")
    fig.savefig(path_noext + ".png", dpi=600)


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmt_p(p, prefix="P = "):
    """P value: 'P = 0.68' for P >= 0.01, else 'P = 1.3x10^-8' with Unicode superscripts (matches R figures)."""
    if p >= 0.01:
        return prefix + ("%.2f" % p)
    m, e = ("%.1e" % p).split("e")
    return prefix + m + "×10" + str(int(e)).translate(_SUP)
