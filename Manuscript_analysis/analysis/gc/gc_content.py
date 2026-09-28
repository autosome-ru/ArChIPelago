"""Mean GC content (%) of the positive and negative sequences of every set used in the manuscript
(human training = chr2-7, 9, 10, 13-20; human test = chr1, 8, 21; mouse training = chr2-7, 9, 10, 13-18;
mouse test = chr1, 8, 19), for the 36 TFs, on exactly the rows that archi_data.select_rows() returns
(up to 10,000 positives and 1,000,000 negatives per TF and set, random_state=0).

Two sequence versions are read (same row order as out_tab_<SP>_10000_<train|control>.tab; the loader maps the logical split to the
file name, see ../common/archi_data.py):
  genomic     $ARCHI_RELEASE_DIR/all_mfa_file_<SP>_10000_<train|control>.fasta  -- 300-bp genome extracts,
              soft-masked (lowercase = RepeatMasker repeats), a few N; GC = (G+C)/(A+C+G+T), case-insensitive,
              N excluded.  These are the values reported in Sup. Table 2.
  modelinput  $ARCHI_RELEASE_DIR/all_mfa_file_<SP>_10000_<train|control>_no_NF_no_N.fasta -- the file scanned by
              the PWMs: identical to the genomic file except that every lowercase base and every N is 'A'
              (HUMAN: 10 % of all bases; MOUSE genome is not soft-masked, only N -> A).  GC of these sequences
              is reported for the record (columns *_modelinput_GC_pct), not in the table.
Output: gc_content.csv (one row per TF).  Read-only on $ARCHI_RELEASE_DIR; one thread, ~2 min."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import archi_data as ad
HERE = os.path.dirname(os.path.abspath(__file__))
TFS = sorted(pd.read_csv(os.path.join(HERE, "table1_ids.csv")).TF.unique())
assert len(TFS) == 36

def gc_array(path, allow_n):
    names, gc = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip())
            else:
                s = line.strip().upper()
                acgt = sum(s.count(c) for c in "ACGT")
                assert len(s) == 300 and (acgt == 300 or (allow_n and acgt + s.count("N") == 300)), (path, len(names), len(s), acgt)
                gc.append((s.count("G") + s.count("C")) / acgt)
    return np.array(names), np.array(gc)

res = {tf: {"TF": tf} for tf in TFS}
for sp, lab in (("HUMAN", "human"), ("MOUSE", "mouse")):
    for split, slab in (("train", "train"), ("control", "test")):
        fsplit = ad._file_split(sp, split)
        gt = ad.global_table(sp, split)
        assert set(gt["chrom"].unique()) == (set(ad.TEST_CHROMS[sp]) if split == "control" else set(gt["chrom"].unique()) - set(ad.TEST_CHROMS[sp]))
        sel = {tf: ad.select_rows(gt, tf, sp, split) for tf in TFS}
        for version, suffix, allow_n in (("genomic", "", True), ("modelinput", "_no_NF_no_N", False)):
            fasta = os.path.join(ad.REL, "all_mfa_file_%s_10000_%s%s.fasta" % (sp, fsplit, suffix))
            names, gc = gc_array(fasta, allow_n)
            assert len(gt) == len(gc) and (gt["name"].values == names).all(), fasta
            print(sp, slab, version, fasta, len(gc), "chroms:", sorted(gt["chrom"].unique(), key=lambda c: int(c[3:])), flush=True)
            for tf in TFS:
                rows, y = sel[tf]
                g = gc[rows]; r = res[tf]
                r["%s_%s_n_pos" % (lab, slab)] = int((y == 1).sum()); r["%s_%s_n_neg" % (lab, slab)] = int((y == 0).sum())
                tag = "" if version == "genomic" else "_modelinput"
                r["%s_%s_pos%s_GC_pct" % (lab, slab, tag)] = 100 * g[y == 1].mean()
                r["%s_%s_neg%s_GC_pct" % (lab, slab, tag)] = 100 * g[y == 0].mean()
out = pd.DataFrame([res[tf] for tf in TFS])
out.to_csv(os.path.join(HERE, "gc_content.csv"), index=False, float_format="%.4f")
print(out.round(3).to_string())
