"""PWM origin table. Per TF and species: the same-species PWMs (mono + di) of the model, and how many of
them were derived from a peak set that is one of that TF's listed positive sets in Sup. Table 1 (the source peak set is the
PEAKS id in the PWM name, e.g. SRF_MOUSE~CM~SRF_MOUSE.PEAKS035686.macs.M).
Inputs: the four PWM archives of the Zenodo deposit v1.0.0 (read with tarfile, md5 checked against md5sums.txt there),
outputs/supplementary_tables/Sup_Table_1_ChIP-Seq_experiments.xlsx (human, mouse), outputs/supplementary_tables/Sup_Table_2_PWMs_and_datasets.xlsx (PWM counts
per TF, must equal the archive counts).
Output: 3_analysis/integration/pwm_origin_h3.csv (TF, species, n_monoPWMs, n_diPWMs, n_PWMs, n_from_listed_sets, frac_from_listed_sets,
n_source_sets, n_source_sets_listed, n_header_tf_mismatch).
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/pwm_origin_h3.py
"""
import os
import re
import sys
import tarfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w2  # noqa: E402

PD, R2, B13 = w2.PD, w2.R2, w2.B13
ZEN = os.environ.get("ARCHI_ZENODO_V1", os.path.join(R2, "results/datasets/zenodo_v1.0.0"))
TARS = {(sp, k): os.path.join(ZEN, "PWMs_%s_%s.tar.gz" % (k, sp)) for sp in ("HUMAN", "MOUSE") for k in ("mono", "di")}
T1 = os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_1_ChIP-Seq_experiments.xlsx")
T2 = os.path.join(R2, "outputs/supplementary_tables", "Sup_Table_2_PWMs_and_datasets.xlsx")
md5s = dict(l.split()[::-1] for l in open(os.path.join(ZEN, "md5sums.txt")) if l.strip())
for f in TARS.values():
    if md5s.get(os.path.basename(f)) != w2.md5(f):
        raise SystemExit("CHECK FAILED %s: md5 differs from the Zenodo md5sums.txt" % f)
t1 = pd.read_excel(T1, sheet_name=None)
t2 = pd.read_excel(T2, sheet_name="per_TF").set_index("TF")
REF = PD.reference_tfs()
NAME = re.compile(r"^>?(\w+?)_(HUMAN|MOUSE)~C[MD]~\w+?_(?:HUMAN|MOUSE)[._](PEAKS\d+)[._]")


def headers(tar):
    """{TF: [header, ...]} of every PWM file <dir>/<TF>/<k>.(d)pwm in the archive."""
    out = {}
    with tarfile.open(tar) as tf_:
        for m in tf_.getmembers():
            parts = m.name.strip("./").split("/")
            if m.isfile() and len(parts) >= 2 and re.search(r"\.d?pwm$", parts[-1]):
                out.setdefault(parts[-2], []).append(tf_.extractfile(m).readline().decode().strip())
    return out


rows = []
for sp, lab in (("HUMAN", "human"), ("MOUSE", "mouse")):
    listed = t1[lab].groupby("TF")["GTRD peak set ID"].apply(set)
    PD.check_tfs(listed.index, "Sup. Table 1 %s" % lab)
    h = {k: headers(TARS[(sp, k)]) for k in ("mono", "di")}
    for tf in REF:
        n = {k: len(h[k].get(tf, [])) for k in ("mono", "di")}
        want = (t2.loc[tf, "%s monoPWMs" % lab.capitalize()], t2.loc[tf, "%s diPWMs" % lab.capitalize()])
        if (n["mono"], n["di"]) != want:
            raise SystemExit("CHECK FAILED %s %s: archive %d/%d PWMs, Sup. Table 2 %d/%d" % (tf, sp, n["mono"], n["di"], *want))
        parsed = [NAME.match(x) for k in ("mono", "di") for x in h[k][tf]]
        if not all(parsed):
            raise SystemExit("CHECK FAILED %s %s: unparsed PWM names %s" % (tf, sp, [x for k in ("mono", "di") for x, p in zip(h[k][tf], parsed) if not p][:3]))
        src = [p.group(3) for p in parsed]
        bad_sp = sum(p.group(2) != sp for p in parsed)
        if bad_sp:
            raise SystemExit("CHECK FAILED %s %s: %d PWMs of the other species" % (tf, sp, bad_sp))
        nl = sum(s in listed[tf] for s in src)
        rows.append(dict(TF=tf, species=lab, n_monoPWMs=n["mono"], n_diPWMs=n["di"], n_PWMs=len(src), n_from_listed_sets=nl,
                         frac_from_listed_sets=nl / len(src), n_source_sets=len(set(src)), n_source_sets_listed=len(set(src) & listed[tf]),
                         n_header_tf_mismatch=sum(p.group(1) != tf for p in parsed)))
out = pd.DataFrame(rows)
dst = os.path.join(B13, "pwm_origin_h3.csv")
w2.backup(dst)
out.to_csv(dst, index=False)
print(out.groupby("species")[["n_PWMs", "n_from_listed_sets", "n_source_sets", "n_source_sets_listed", "n_header_tf_mismatch"]].sum())
w2.provenance(B13, "3_analysis/integration/pwm_origin_h3.py", [dst], list(TARS.values()) + [T1, T2],
              params="{source peak set: PEAKS id of the PWM name; listed: GTRD peak set ID of the TF in Sup. Table 1 R2, same species}")
