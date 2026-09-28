"""Pooled human peak set of every TF (Sup. Table 2, "Human pooled peak set"), computed with the rules of notebook 0
(0_Data_preparation, cells 13/18) from the Sup. Table 1 experiments (experiment_ids.csv); also records which
datasets and how many model files (.M/.S) contributed peaks.  Read-only on ~/macs, ~/hocomoco11, ~/Cleaning_data;
writes pooled_peaks_human.csv next to this script."""
import os, sys
import numpy as np, pandas as pd
import pybedtools as pbt
HERE = os.path.dirname(os.path.abspath(__file__))
H = os.path.expanduser("~")
t1 = pd.read_csv(os.path.join(HERE, "experiment_ids.csv"))
sure = set(pd.read_csv(f"{H}/hocomoco11/curation/slices4bench_mono/out_sure.csv", header=None)[0])
sure_macs = [s for s in sure if "macs" in s]
rows = []
for organism, sp_name, rep in (("HUMAN", "Homo sapiens", "track_out.bed"),):
    ids_sp = t1[t1.specie == sp_name]
    files = sorted(os.listdir(f"{H}/hocomoco11/auc/mono/{organism}_datasets/"))
    repeats = pbt.BedTool(f"{H}/Cleaning_data/{rep}")
    for TF in sorted(ids_sp.TF.unique()):
        ids = set(ids_sp[ids_sp.TF == TF].peak_id)
        sel = [x for x in files if x.split("_")[0] == TF and organism in x and x.split("~")[-1].split(".txt")[0] in sure_macs
               and x.split(".")[1] in ids]
        tmp = []
        for j in sel:
            f = pd.read_csv(f"{H}/hocomoco11/auc/mono/{organism}_datasets/{j}", sep="\t", header=None)
            if np.percentile(f[1], 10) > 0.5:
                tmp.append(j)
        tmp = sorted(set(tmp))
        peaks = pd.DataFrame(); n_missing = 0
        for i in tmp:
            p = f"{H}/macs/{i.split('.')[1]}.interval"
            if not os.path.exists(p):
                n_missing += 1; continue
            df = pd.read_csv(p, sep="\t")
            df = df[df["tags"] >= 10].sort_values(by=["-10*log10(pvalue)"], ascending=False).iloc[:1000]
            peaks = pd.concat([peaks, df], ignore_index=True)
        n_after_repeats = 0; n_uniq = 0
        if peaks.shape[0] > 0:
            bed = pd.DataFrame({"chrom": peaks["#CHROM"].apply(lambda t: f"chr{t}"),
                                "start": peaks["START"] + peaks["summit"] - 500, "end": peaks["START"] + peaks["summit"] + 500,
                                "score": peaks["-10*log10(pvalue)"], "summit": peaks["START"] + peaks["summit"]})
            bed = bed[(bed.chrom != "chrnan") & (bed.chrom != "chrM") & (bed.chrom.str.len() < 6)]
            bt = pbt.BedTool.from_dataframe(bed)
            out = bt.subtract(repeats, f=0.7, N=True).to_dataframe()
            n_after_repeats = out.shape[0]
            n_uniq = out.drop_duplicates(subset=["chrom", "start", "end"]).shape[0]
            n_chromfilt = bed.shape[0]
        else:
            n_chromfilt = 0
        rows.append(dict(TF=TF, species=organism, datasets_table1=len(ids), model_files_selected=len(sel), model_files_p10=len(tmp),
                         datasets_p10=len({i.split(".")[1] for i in tmp}), interval_missing=n_missing,
                         peaks_top1000=peaks.shape[0], peaks_chromfilter=n_chromfilt, peaks_after_repeats=n_after_repeats,
                         peaks_unique=n_uniq))
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(os.path.join(HERE, "pooled_peaks_human.csv"), index=False)
