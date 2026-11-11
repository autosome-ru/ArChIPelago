"""check_E.py - check: no duplicate windows from repeated reads of a peak set.
(1) every peak set of the TF was read once (per_tf_positive_counts: peak_set_reads == peak_sets);
(2) identical windows (chrom, start, end) among the positives of a TF (train + test bed): count, and for each such
    window the number of distinct listed peak sets whose top 1,000 (tags >= 10) holds its summit (must be >= copies,
    i.e. the copies come from different sets);
(3) the same duplicate count on the published human beds (~/archi_release/TF-ML/{train,test}) as a control that must be red.
usage: python check_E.py ROOT_TFML PEAKS_CSV SPECIES OUTCSV   (refuses to overwrite)"""
import os, sys
from collections import Counter
import numpy as np, pandas as pd
root, peaks, sp, out = sys.argv[1:5]
if os.path.exists(out): raise SystemExit("refusing to overwrite " + out)
H = os.path.expanduser("~"); OLD = f"{H}/Release/TF-ML"
pk = pd.read_csv(peaks); cnt = pd.read_csv(f"{root}/per_tf_positive_counts.csv").set_index("TF")


def windows(r, tf):
    b = pd.concat([pd.read_csv(f"{r}/{s}/{tf}_{s}_set_{sp}.bed", sep="\t", index_col=0) for s in ("train", "test")], ignore_index=True)
    return b


rows = []
for tf in sorted(pk.TF.unique()):
    b = windows(root, tf)
    c = Counter(zip(b.chrom, b.start, b.end))
    dup = {k: v for k, v in c.items() if v > 1}
    src = {}
    for pid in sorted(set(pk.peak_id[pk.TF == tf])):
        d = pd.read_csv(f"{H}/macs/{pid}.interval", sep="\t").query("tags >= 10").sort_values(by=["-10*log10(pvalue)"], ascending=False).iloc[:1000]
        for k in set(zip("chr" + d["#CHROM"].astype(str), (d.START + d.summit - 500).clip(lower=0), d.START + d.summit + 500)):
            src.setdefault(k, set()).add(pid)
    bad = sum(1 for k, v in dup.items() if len(src.get(k, ())) < v)
    o = windows(OLD, tf) if sp == "HUMAN" else None
    oc = Counter(zip(o.chrom, o.start, o.end)) if o is not None else Counter()
    rows.append(dict(TF=tf, peak_sets=int(cnt.loc[tf, "peak_sets"]), peak_set_reads=int(cnt.loc[tf, "peak_set_reads"]),
                     positives=len(b), unique_windows=len(c), windows_with_copies=len(dup),
                     extra_copies=int(sum(v - 1 for v in dup.values())), dup_windows_from_fewer_sets_than_copies=bad,
                     orig_positives=len(o) if o is not None else np.nan, orig_unique_windows=len(oc) if o is not None else np.nan,
                     orig_extra_copies=int(sum(v - 1 for v in oc.values() if v > 1)) if o is not None else np.nan))
    print(rows[-1], flush=True)
df = pd.DataFrame(rows); df.to_csv(out, index=False)
print("TOTAL", df.drop(columns="TF").sum().to_dict())
