"""Flip report of outputs/numbers/NUMBERS_R2_35TF.csv: every row whose 'new value' differs from the reference column
('h4 value' by default, or 'h3a value') in a way the text has to follow:
  count  - a TF count 'n/N' or 'n of N' changes (crossing = one side is N of N, e.g. 35/35 -> 34/35);
  sign   - a number changes sign (direction of a gain or a correlation), P and q values excluded;
  signif - a P or q value crosses ALPHA (default 0.05): numbers labelled 'P' / 'q' in the values, or, when the values carry no
           label, every number in [0, 1] of a row whose item names a P or q;
  shape  - the two values hold a different number of numbers (a list or a statement changed; read the row).
Rows that are new, or equal to the reference, are not reported. Output: outputs/numbers/flip_report_<new>_vs_<ref>.csv and a
count per type on stdout.
usage: python3 3_analysis/integration/flip_report.py [--ref h4|h3a] [--alpha 0.05] [--numbers FILE] [--out FILE]"""
import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h3_inputs as H  # noqa: E402

NUM = re.compile(r"(?<![\w.])[-+−]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:e[-+]?\d+)?(?![\w])")
COUNT = re.compile(r"(?<![\d.,])(\d+)\s*(?:/|of)\s*(\d+)(?![\d.,])")


def counts(v):
    """TF counts 'n/N' or 'n of N' (n <= N, whole numbers only); 'a / b of N' and 'a, b, c of N' are expanded to one count each."""
    s = str(v)
    s = re.sub(r"(?<![\d.,])((?:\d+\s*(?:,|/)\s*)+\d+)\s+of\s+(\d+)(?![\d.,])",
               lambda m: "; ".join("%s of %s" % (x, m.group(2)) for x in re.split(r"\s*[,/]\s*", m.group(1))), s)
    return [(a, b) for a, b in COUNT.findall(s) if int(a) <= int(b)]
SUP = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
PWORD = re.compile(r"(?<![A-Za-z])[Pq](?![A-Za-z])")
PLAB = re.compile(r"(?<![A-Za-z])[Pq]\s*(?:=|<|≤)?\s*([0-9][0-9.e+-]*(?:\s*×\s*10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)?)")


def nums(v):
    s = str(v).replace("−", "-")
    s = re.sub(r"(\d(?:\.\d+)?)\s*×\s*10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)", lambda m: "%se%s" % (m.group(1), m.group(2).translate(SUP)), s)
    s = re.sub(r"(?<![\d.,])((?:\d+\s*(?:,|/)\s*)*\d+)\s+of\s+\d+(?![\d.,])", " ", s)   # counts are compared separately
    s = re.sub(r"(?<![\d.,])\d+\s*/\s*\d+(?![\d.,])", lambda m: " " if len(re.findall(r"\d+", m.group(0))) == 2 and
               int(re.findall(r"\d+", m.group(0))[0]) <= int(re.findall(r"\d+", m.group(0))[1]) else m.group(0).replace("/", " "), s)
    return [float(x.replace(",", "")) for x in NUM.findall(s)]


def flips(item, note, ref, new, alpha):
    if pd.isna(ref) or str(ref) == "" or str(ref) == str(new):
        return []
    out = []
    cr, cn = counts(ref), counts(new)
    for (a, n1), (b, n2) in zip(cr, cn):
        if (a, n1) != (b, n2):
            out.append(("count", "%s/%s -> %s/%s%s" % (a, n1, b, n2, " (crossing all)" if (a == n1) != (b == n2) else "")))
    if len(cr) != len(cn):
        out.append(("shape", "%d counts -> %d counts" % (len(cr), len(cn))))
    xr, xn = nums(ref), nums(new)
    if len(xr) != len(xn):
        out.append(("shape", "%d numbers -> %d numbers" % (len(xr), len(xn))))
        return out
    pr, pn = [nums(x)[0] for x in PLAB.findall(str(ref).replace("−", "-"))], [nums(x)[0] for x in PLAB.findall(str(new).replace("−", "-"))]
    if pr and len(pr) == len(pn):     # labelled P / q values in the cells
        ps = list(zip(pr, pn))
        rest = [(a, b) for a, b in zip(xr, xn) if not any(a == x and b == y for x, y in ps)]
    elif PWORD.search(str(item)):     # unlabelled: the item says the row is a P / q; numbers in [0, 1] are P
        ps = [(a, b) for a, b in zip(xr, xn) if 0 <= a <= 1 and 0 <= b <= 1]
        rest = [(a, b) for a, b in zip(xr, xn) if not (0 <= a <= 1 and 0 <= b <= 1)]
    else:
        ps, rest = [], list(zip(xr, xn))
    out += [("signif", "%g -> %g (alpha %g)" % (a, b, alpha)) for a, b in ps if (a < alpha) != (b < alpha)]
    out += [("sign", "%g -> %g" % (a, b)) for a, b in rest if a * b < 0]
    return out


def report(d, ref_col, alpha):
    rows = []
    for _, r in d.iterrows():
        for kind, det in flips(r["item"], r.get("statement_note", ""), r[ref_col], r["new value"], alpha):
            rows.append(dict(item=r["item"], kind=kind, detail=det, ref_value=r[ref_col], new_value=r["new value"],
                             source=r["source file"], inventory=r["inventory row"]))
    return pd.DataFrame(rows, columns=["item", "kind", "detail", "ref_value", "new_value", "source", "inventory"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="h4", choices=["h4", "h3a"])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--numbers", default=os.path.join(H.R2, "outputs/numbers/NUMBERS_R2_35TF.csv"))
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    d = pd.read_csv(a.numbers, dtype=str, keep_default_na=False)
    ref_col = a.ref + " value"
    if ref_col not in d.columns:
        raise SystemExit("CHECK FAILED %s has no column %r" % (a.numbers, ref_col))
    rep = report(d, ref_col, a.alpha)
    out = a.out or os.path.join(H.R2, "outputs/numbers", "flip_report_%s_vs_%s.csv" % (H.STAGE.lower(), a.ref))
    rep.to_csv(out, index=False)
    print("flip report %s vs %s: %d rows, %d flips (%s) -> %s" % (H.STAGE, a.ref, rep.item.nunique(), len(rep),
          ", ".join("%s %d" % (k, n) for k, n in rep.kind.value_counts().items()) or "none", os.path.relpath(out, H.R2)))
    return rep


if __name__ == "__main__":
    main()
