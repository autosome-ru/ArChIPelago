#!/usr/bin/env python3
"""Sup. Table 1 (ChIP-Seq experiments) and Sup. Table 2 (PWMs, datasets, peaks, set sizes, GC content) of the
manuscript, built from verified sources and written in the style of Sup. Tables 3-6
(pandas to_excel, bold header row, descriptive column names, one sheet per logical block, no unlabelled rows).

Basis: mouse test set = mouse chromosomes 1, 8 and 19; mouse training set = chromosomes 2-7, 9, 10 and 13-18;
human test set = chromosomes 1, 8 and 21; human training set = chromosomes 2-7, 9, 10 and 13-20.

Sources
  Sup. Table 1
    inputs/Table_1_GTRD_experiments.xlsx          the list of experiments (TF, GTRD peak-set id) per species
    inputs/metadata_new_23_02_17_with_control.csv
                                                  the GTRD metadata used by the pipeline (notebook 0, cell 4): experiment id,
                                                  species, TFClass id, control experiment id, UniProt entry name and
                                                  accession, GEM peak count.  Parsed positionally (unquoted commas inside
                                                  fields shift the columns; every row has 32 fields).
  Sup. Table 2 (per TF)
    fig4/rows/<TF>.json (n_mono, n_di)      human monoPWMs / diPWMs = feature matrix of the models
    hm_mm/hm_mm_results.csv                       P_mono_H/P_di_H (cross-check), P_mono_M/P_di_M = mouse PWMs;
                                                  n_train_pos_H, n_test_pos_H, n_train_pos_M, n_test_pos_M = positives
    gc/peaks_raw_check.csv, gc/peaks_raw_check_mouse_allmodels.csv
                                                  pooled peak sets re-derived from the MACS peak files with the pipeline's
                                                  rules (peaks_after_repeats); must equal the counts of inputs/Table_2_TF_metadata.xlsx
    gc/gc_content.csv                             mean GC (%) of the positives and negatives of every set, genomic
                                                  sequences (300 bp), exactly the rows of archi_data.select_rows()
    inputs/Table_2_TF_metadata.xlsx               the TF table of the pipeline, used only for the verification report

Outputs (Manuscript_analysis/Sup_Tables/): Sup_Table_1_ChIP-Seq_experiments.xlsx (+ _human.csv, _mouse.csv),
Sup_Table_2_PWMs_and_datasets.xlsx (+ .csv, _summary.csv).  A verification report is printed to stdout.
"""
import csv
import json
import os
import re
import numpy as np
import pandas as pd

from archi_paths import INPUTS, BASIS, RES

OUT = os.path.join(RES, "Sup_Tables")
os.makedirs(OUT, exist_ok=True)
ORIG_T1 = os.path.join(INPUTS, "Table_1_GTRD_experiments.xlsx")
ORIG_T2 = os.path.join(INPUTS, "Table_2_TF_metadata.xlsx")
META = os.path.join(INPUTS, "metadata_new_23_02_17_with_control.csv")
GC = os.path.join(BASIS, "gc")
HUMAN_TEST, HUMAN_TRAIN = "chr1, 8, 21", "chr2-7, 9, 10, 13-20"
MOUSE_TEST, MOUSE_TRAIN = "chr1, 8, 19", "chr2-7, 9, 10, 13-18"

report = []


def note(msg):
    report.append(msg)
    print(msg)


# ------------------------------------------------------------------------------------------ GTRD metadata
def read_metadata():
    """PEAKS id -> fields.  Row layout (32 fields, commas inside fields shift the later ones): exp, species, TFClass, ...,
    PEAKS id, TF name, [control EXP id, control antibody], external refs (GEO:/PUBMED:), UniProt accession, UniProt entry
    name, control refs (GEO:/PUBMED:, may be empty), GEM, SISSRs, MACS, MACS-filtered, PICS peak counts, padding."""
    meta = {}
    with open(META, newline="") as fh:
        for row in csv.reader(fh):
            ip = [i for i, v in enumerate(row) if re.fullmatch(r"PEAKS\d{6}", v)]
            if not ip:
                continue
            ip = ip[0]
            ctrl = [v for v in row[ip + 1:] if re.fullmatch(r"EXP\d{6}", v)]
            iu = [i for i, v in enumerate(row) if re.fullmatch(r"[A-Z0-9]+_(HUMAN|MOUSE)", v)]
            if not iu:
                raise ValueError("no UniProt entry name in metadata row of %s" % row[ip])
            iu = iu[0]
            rest = row[iu + 1:]
            j = 0
            while j < len(rest) and ":" in rest[j]:       # control refs
                j += 1
            if j == 0:                                     # no control refs: the empty control_refs field comes first
                if rest[0] != "":
                    raise ValueError("unexpected token after UniProt entry name in row of %s: %r" % (row[ip], rest[:3]))
                j = 1
            counts = rest[j:j + 5]
            if not all(v == "" or v.isdigit() for v in counts):
                raise ValueError("peak counts not numeric in row of %s: %r" % (row[ip], counts))
            meta[row[ip]] = dict(exp=row[0], species=row[1], tfclass=row[2], control=ctrl[0] if ctrl else "",
                                 uniprot_ac=row[iu - 1], uniprot_id=row[iu], gem=int(counts[0]) if counts[0] else np.nan)
    return meta


meta = read_metadata()
note("GTRD metadata: %d peak sets parsed from %s" % (len(meta), os.path.basename(META)))

# ------------------------------------------------------------------------------------------ Sup. Table 1
T1_COLS = ["TF", "GTRD peak set ID", "GTRD experiment ID", "Species", "TFClass ID", "GTRD control experiment ID",
           "UniProt entry name", "UniProt accession", "GEM peaks (GTRD peak calls)", "UniProt short name"]
t1 = {}
for sp, sheet, ent in (("human", "Sheet1", "HUMAN"), ("mouse", "Sheet2", "MOUSE")):
    o = pd.read_excel(ORIG_T1, sheet_name=sheet)
    if o.peak_id.duplicated().any():
        raise ValueError("duplicated peak ids in Table_1.xlsx %s" % sheet)
    rows, diffs = [], []
    for r in o.itertuples(index=False):
        m = meta[r.peak_id]
        if m["uniprot_id"] != "%s_%s" % (r.TF, ent):
            raise ValueError("%s: metadata UniProt entry name %s != %s_%s" % (r.peak_id, m["uniprot_id"], r.TF, ent))
        rows.append([r.TF, r.peak_id, m["exp"], m["species"], m["tfclass"], m["control"], m["uniprot_id"], m["uniprot_ac"],
                     m["gem"], r.Uniprot_name])
        for col, key in (("exp", "exp"), ("specie", "species"), ("TFclass", "tfclass"), ("control_id", "control"),
                         ("Uniprot_AC", "uniprot_id")):
            if str(getattr(r, col)) != m[key]:
                diffs.append((r.TF, r.peak_id, col, getattr(r, col), m[key]))
        old_gem = getattr(r, "gem_peaks")
        if not (str(old_gem) == "abrakadabra" and pd.isna(m["gem"])) and not (str(old_gem).isdigit() and int(old_gem) == m["gem"]):
            diffs.append((r.TF, r.peak_id, "gem_peaks", old_gem, m["gem"]))
    t1[sp] = pd.DataFrame(rows, columns=T1_COLS)
    t1[sp]["GEM peaks (GTRD peak calls)"] = t1[sp]["GEM peaks (GTRD peak calls)"].astype("Int64")
    n_missing = int(t1[sp]["GEM peaks (GTRD peak calls)"].isna().sum())
    note("Sup. Table 1 %s: %d experiments, %d TFs; fields identical to the GTRD metadata for every experiment: %s; "
         "GEM count absent in the metadata for %d experiments (%s)"
         % (sp, len(t1[sp]), t1[sp].TF.nunique(), "yes" if not diffs else "NO: %s" % diffs, n_missing,
            ", ".join(t1[sp].loc[t1[sp]["GEM peaks (GTRD peak calls)"].isna(), "GTRD peak set ID"])))
    if diffs:
        raise SystemExit("Sup. Table 1: the published table disagrees with the GTRD metadata: %s" % diffs)
if not (len(t1["human"]) == 414 and len(t1["mouse"]) == 290 and set(t1["human"].TF) == set(t1["mouse"].TF) and t1["human"].TF.nunique() == 36):
    raise SystemExit("Sup. Table 1: expected 414 human + 290 mouse experiments of the same 36 TFs")
TFS = sorted(t1["human"].TF.unique())

with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_1_ChIP-Seq_experiments.xlsx")) as xw:
    t1["human"].to_excel(xw, sheet_name="human", index=False)
    t1["mouse"].to_excel(xw, sheet_name="mouse", index=False)
t1["human"].to_csv(os.path.join(OUT, "Sup_Table_1_ChIP-Seq_experiments_human.csv"), index=False)
t1["mouse"].to_csv(os.path.join(OUT, "Sup_Table_1_ChIP-Seq_experiments_mouse.csv"), index=False)

# ------------------------------------------------------------------------------------------ Sup. Table 2 sources
hm = pd.read_csv(os.path.join(BASIS, "hm_mm", "hm_mm_results.csv")).set_index("TF")
if sorted(hm.index) != TFS:
    raise SystemExit("hm_mm_results.csv does not hold the 36 TFs of Sup. Table 1")
gc = pd.read_csv(os.path.join(GC, "gc_content.csv")).set_index("TF").loc[TFS]
prH = pd.read_csv(os.path.join(GC, "peaks_raw_check.csv"))
prH = prH[prH.species == "HUMAN"].set_index("TF").loc[TFS]
prM = pd.read_csv(os.path.join(GC, "peaks_raw_check_mouse_allmodels.csv")).set_index("TF").loc[TFS]
old2 = pd.read_excel(ORIG_T2)
old2 = old2[old2.TF_name.notna()].set_index("TF_name").loc[TFS]

t2 = pd.DataFrame(index=TFS)
t2.index.name = "TF"
for tf in TFS:
    j = json.load(open(os.path.join(BASIS, "fig4", "rows", tf + ".json")))
    if (j["n_mono"], j["n_di"]) != (int(hm.loc[tf, "P_mono_H"]), int(hm.loc[tf, "P_di_H"])):
        raise SystemExit("%s: human PWM counts differ between fig4 rows and hm_mm_results" % tf)
    t2.loc[tf, "Human monoPWMs"] = j["n_mono"]
    t2.loc[tf, "Human diPWMs"] = j["n_di"]
t2["Human PWMs (mono + di)"] = t2["Human monoPWMs"] + t2["Human diPWMs"]
t2["Mouse monoPWMs"] = hm.loc[TFS, "P_mono_M"].values
t2["Mouse diPWMs"] = hm.loc[TFS, "P_di_M"].values
t2["Mouse PWMs (mono + di)"] = t2["Mouse monoPWMs"] + t2["Mouse diPWMs"]
t2["Human ChIP-Seq experiments (Sup. Table 1)"] = t1["human"].groupby("TF").size().loc[TFS].values
t2["Mouse ChIP-Seq experiments (Sup. Table 1)"] = t1["mouse"].groupby("TF").size().loc[TFS].values
t2["Human pooled peak set (peaks, all chromosomes)"] = prH["peaks_after_repeats"].values
t2["Mouse pooled peak set (peaks, all chromosomes)"] = prM["peaks_after_repeats"].values
t2["Human training positives (%s)" % HUMAN_TRAIN] = hm.loc[TFS, "n_train_pos_H"].values
t2["Human test positives (%s)" % HUMAN_TEST] = hm.loc[TFS, "n_test_pos_H"].values
t2["Mouse training positives (%s)" % MOUSE_TRAIN] = hm.loc[TFS, "n_train_pos_M"].values
t2["Mouse test positives (%s)" % MOUSE_TEST] = hm.loc[TFS, "n_test_pos_M"].values
GC_COLS = []
for lab, sp, tr, te in (("Human", "human", HUMAN_TRAIN, HUMAN_TEST), ("Mouse", "mouse", MOUSE_TRAIN, MOUSE_TEST)):
    for slab, chroms in (("training", tr), ("test", te)):
        for cls in ("positives", "negatives"):
            col = "%s %s %s, mean GC content (%%; %s)" % (lab, slab, cls, chroms)
            t2[col] = gc["%s_%s_%s_GC_pct" % (sp, "train" if slab == "training" else "test", cls[:3])].round(2).values
            GC_COLS.append(col)
count_cols = [c for c in t2.columns if c not in GC_COLS]
t2[count_cols] = t2[count_cols].astype(int)

# cross-checks between the sources
for col, key in (("Human training positives (%s)" % HUMAN_TRAIN, "human_train_n_pos"), ("Human test positives (%s)" % HUMAN_TEST, "human_test_n_pos"),
                 ("Mouse training positives (%s)" % MOUSE_TRAIN, "mouse_train_n_pos"), ("Mouse test positives (%s)" % MOUSE_TEST, "mouse_test_n_pos")):
    if not (t2[col].values == gc[key].values).all():
        raise SystemExit("positives of %s differ between hm_mm_results and gc_content" % col)
note("Sup. Table 2: positives of hm_mm_results.csv == rows counted by gc_content.py for all four sets: yes")
for sp, pr, lab in (("HUMAN", prH, "human"), ("MOUSE", prM, "mouse")):
    same = (pr["peaks_after_repeats"].values == old2["%s_peaks_raw" % sp].astype(int).values)
    note("Sup. Table 2: %s pooled peak set re-derived from the MACS peak files equals the published %s_peaks_raw for %d/36 TFs"
         % (lab, sp, int(same.sum())))
    if not same.all():
        raise SystemExit("%s peaks_raw not reproduced for %s" % (lab, list(np.array(TFS)[~same])))

# ------------------------------------------------------------------------------------------ summary sheet
summ = []
for stat, fn in (("Total over 36 TFs", lambda s: int(s.sum())), ("Median per TF", lambda s: float(s.median())), ("Mean per TF", lambda s: float(s.mean()))):
    row = {"Statistic": stat}
    for c in t2.columns:
        if c in GC_COLS and stat.startswith("Total"):
            row[c] = np.nan
        else:
            v = fn(t2[c])
            row[c] = round(v, 2) if isinstance(v, float) else v
    summ.append(row)
sum2 = pd.DataFrame(summ)

t2r = t2.reset_index()
with pd.ExcelWriter(os.path.join(OUT, "Sup_Table_2_PWMs_and_datasets.xlsx")) as xw:
    sum2.to_excel(xw, sheet_name="summary", index=False)
    t2r.to_excel(xw, sheet_name="per_TF", index=False)
t2r.to_csv(os.path.join(OUT, "Sup_Table_2_PWMs_and_datasets.csv"), index=False)
sum2.to_csv(os.path.join(OUT, "Sup_Table_2_PWMs_and_datasets_summary.csv"), index=False)

# ------------------------------------------------------------------------------------------ verification report: old vs new
note("")
note("Sup. Table 2, values that differ from the published Table_2.xlsx (TF: column, published -> new):")
pairs = [("HUMAN_monoPWMs", "Human monoPWMs"), ("HUMAN_diPWMs", "Human diPWMs"), ("MOUSE_monoPWMs", "Mouse monoPWMs"), ("MOUSE_diPWMs", "Mouse diPWMs"),
         ("HUMAN_datasets", "Human ChIP-Seq experiments (Sup. Table 1)"), ("MOUSE_datasets", "Mouse ChIP-Seq experiments (Sup. Table 1)"),
         ("HUMAN_peaks_raw", "Human pooled peak set (peaks, all chromosomes)"), ("MOUSE_peaks_raw", "Mouse pooled peak set (peaks, all chromosomes)"),
         ("HUMAN_train_set", "Human training positives (%s)" % HUMAN_TRAIN), ("HUMAN_test_set", "Human test positives (%s)" % HUMAN_TEST)]
n_diff = 0
for old, new in pairs:
    d = t2.index[old2[old].astype(int).values != t2[new].values]
    for tf in d:
        note("  %s: %s, %d -> %d" % (tf, new, int(old2.loc[tf, old]), int(t2.loc[tf, new])))
        n_diff += 1
note("  (%d differing count values)" % n_diff)
note("  MOUSE_test_set of the published table (positives of the file then used as the mouse test set) is replaced by the "
     "positives of the mouse test set (%s) and of the mouse training set (%s):" % (MOUSE_TEST, MOUSE_TRAIN))
for tf in TFS:
    note("  %s: %d -> test %d, training %d" % (tf, int(old2.loc[tf, "MOUSE_test_set"]),
                                              t2.loc[tf, "Mouse test positives (%s)" % MOUSE_TEST], t2.loc[tf, "Mouse training positives (%s)" % MOUSE_TRAIN]))
note("")
note("GC content, published -> new (mean over TFs of the difference new - published, and its range):")
for old, new in (("HUMAN_train_pos_GC%", GC_COLS[0]), ("HUMAN_train_neg_GC%", GC_COLS[1]), ("HUMAN_test_pos_GC%", GC_COLS[2]), ("HUMAN_test_neg_GC%", GC_COLS[3]),
                 ("MOUSE_test_pos_GC%", GC_COLS[4]), ("MOUSE_test_neg_GC%", GC_COLS[5])):
    d = t2[new].values - old2[old].astype(float).values
    note("  %s -> %s: mean %+.2f, range %+.2f to %+.2f, |diff| > 0.1 for %d/36 TFs" % (old, new, d.mean(), d.min(), d.max(), int((np.abs(d) > 0.1).sum())))
note("  (published MOUSE_test_* GC columns describe the %s positives, i.e. the current mouse training set; the mouse test set %s is new)"
     % (MOUSE_TRAIN, MOUSE_TEST))
note("")
note("Summary sheet: totals over 36 TFs: human %d mono + %d di = %d PWMs, mouse %d mono + %d di = %d PWMs; experiments %d human, %d mouse; "
     "medians per TF: human monoPWMs %.1f, diPWMs %.1f, mono + di %.1f; mouse %.1f, %.1f, %.1f"
     % (t2["Human monoPWMs"].sum(), t2["Human diPWMs"].sum(), t2["Human PWMs (mono + di)"].sum(),
        t2["Mouse monoPWMs"].sum(), t2["Mouse diPWMs"].sum(), t2["Mouse PWMs (mono + di)"].sum(),
        t2["Human ChIP-Seq experiments (Sup. Table 1)"].sum(), t2["Mouse ChIP-Seq experiments (Sup. Table 1)"].sum(),
        t2["Human monoPWMs"].median(), t2["Human diPWMs"].median(), t2["Human PWMs (mono + di)"].median(),
        t2["Mouse monoPWMs"].median(), t2["Mouse diPWMs"].median(), t2["Mouse PWMs (mono + di)"].median()))
note("written: %s" % ", ".join(sorted(f for f in os.listdir(OUT) if f.startswith("Sup_Table_1") or f.startswith("Sup_Table_2"))))
