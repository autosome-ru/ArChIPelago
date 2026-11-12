"""
build_positives_v5.py - build the positive windows of the TF-ML data for one species (--species HUMAN|MOUSE).

Re-implements 1:1 the mouse and human branches of the release notebook 0 (cell numbers are given in the comments).
The chromosome lists are taken from h5_splits.py: human training chr2-7, 9-20, 22 / test chr1, 8, 21; mouse training
chr2-7, 9-18 / test chr1, 8, 19. On disk, mouse files named 'train' hold the test chromosomes 1, 8 and 19. Both
species use one code path; for HUMAN a "factor" column (= TF) is appended to the per-TF bed frames after the name join.

Two stages:
  per_tf : peak sets -> top 1,000 (tags >= 10, by -10log10 P) -> +-500 bp windows (find_peaks) -> repeat subtract
           (f=0.7, N=True) -> split (file "train" / file "test"/"control" as listed by h5_splits.on_disk()) -> cap
           10,000 per split (random_state=1) -> bed + per-TF mfa (bedtools getfasta, name=True).
  pooled : per-TF mfa -> all_mfa_file_<SPECIES>_10000_{train,control}.fasta (seq[350:650]), SLIM variants (1,000 bp,
           '$>chrom: ...' header), headerStr csvs, out_tab (family code from core_annotation_HUMAN_mono.csv),
           awk 'NF' -> _no_NF, sed non-ACGT -> A -> _no_NF_no_N.

usage:
  python build_positives_v5.py per_tf --species HUMAN|MOUSE --peaks CSV --mask BED --out TFMLDIR [--dup] [--tfs A B]
  python build_positives_v5.py pooled --species HUMAN|MOUSE --mfa-dir DIR_WITH_train_test --out TFMLDIR --tfs-from-mfa
--dup: read each peak set n_motif_files times (column of the peak-id csv), as the original notebook did.
Never overwrites: refuses if an output exists.
"""
import os, sys, csv, argparse, subprocess, glob
from itertools import groupby
from collections import defaultdict
import numpy as np, pandas as pd
import pybedtools as pbt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h5_splits  # the one definition of the split

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")); import archi_env as AE
MACS = AE.macs() + "/"
TEMPLATE = AE.macs("PEAKS041820.interval")
FAMILY = AE.family("core_annotation_HUMAN_mono.csv")
TOP, TAGS_MIN, WINDOW, CAP, SIZE_OF_SET = 1000, 10, 1000, 10000, 10000
SPECIES = {
    "MOUSE": dict(GENOME=AE.genomes("GRCm38.p6.genome.fa"),                               # notebook cell 61
                  TRAIN_CHR=h5_splits.on_disk("MOUSE")[0],  # on-disk "train" = mouse test set
                  TEST_CHR=h5_splits.on_disk("MOUSE")[1],  # on-disk "test" = mouse training set
                  FACTOR=False),
    "HUMAN": dict(GENOME=AE.genomes("hg38_genome.fasta"),                                 # notebook cell 34
                  TRAIN_CHR=h5_splits.on_disk("HUMAN")[0],  # chr2-7, 9-20, 22
                  TEST_CHR=h5_splits.on_disk("HUMAN")[1],  # chr1, 8, 21 (cell 28)
                  FACTOR=True),
}
GENOME = TRAIN_CHR = TEST_CHR = organism = FACTOR = None     # set by set_species()


def set_species(sp):
    global GENOME, TRAIN_CHR, TEST_CHR, organism, FACTOR
    c = SPECIES[sp]
    bad = h5_splits.partition_problems(sp, c["TRAIN_CHR"], c["TEST_CHR"])  # refuse an invalid split
    if bad:
        raise SystemExit("invalid chromosome split: %s" % bad)
    GENOME, TRAIN_CHR, TEST_CHR, FACTOR = c["GENOME"], c["TRAIN_CHR"], c["TEST_CHR"], c["FACTOR"]
    organism = sp


def no_clobber(p):
    if os.path.exists(p):
        raise SystemExit("refusing to overwrite " + p)
    return p


def find_peaks(df_in):                                                   # cell 21, verbatim
    df_tmp = pd.DataFrame(columns=list(pd.read_csv(TEMPLATE, index_col=None, sep="\t").columns) + ["ind"])
    window = WINDOW
    df_tmp["abs_summit"] = df_in["START"] + df_in["summit"]
    df_tmp["START"] = df_in["START"] + df_in["summit"] - window // 2
    df_tmp[df_tmp["START"] < 0] = 0
    df_tmp["END"] = df_in["START"] + df_in["summit"] + window // 2
    df_tmp = df_tmp.drop(columns=["summit", "tags", "FDR(%)", "fold_enrichment"])
    df_tmp['#CHROM'] = df_in['#CHROM'].apply(lambda t: f"chr{t}")
    df_tmp = df_tmp[(df_tmp["#CHROM"] != 'chrnan') & (df_tmp["#CHROM"] != 'chrM') & (df_tmp["#CHROM"].str.len() < 6)]
    df_tmp["-10*log10(pvalue)"] = df_in["-10*log10(pvalue)"]
    df_tmp = df_tmp[["#CHROM", "START", "END", "-10*log10(pvalue)", "abs_summit"]]
    return df_tmp


def repeats_cleaner(df_in, mask):                                        # cell 54, mask path as parameter
    replic = pbt.BedTool.from_dataframe(df_in)
    repeats = pbt.BedTool(mask)
    without_reps = replic.subtract(repeats, f=0.7, N=True)
    return without_reps.to_dataframe()


def fasta_iter(f):
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)


def per_tf(a):
    out = a.out
    for d in ("train", "test", "tmp_pbt"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    pbt.set_tempdir(os.path.join(out, "tmp_pbt"))
    peaks = pd.read_csv(a.peaks)
    tfs = sorted(peaks["TF"].unique()) if not a.tfs else a.tfs
    counts = []
    all_tr, all_te = pd.DataFrame(), pd.DataFrame()
    for TF in tfs:
        sub = peaks[peaks["TF"] == TF]
        if a.dup:
            ids = sorted(pid for pid, n in zip(sub["peak_id"], sub["n_motif_files"]) for _ in range(int(n)))
        else:
            ids = sorted(set(sub["peak_id"]))
        tmp = pd.DataFrame()
        missing = []
        for pid in ids:                                                    # cell 55
            fp = MACS + pid + ".interval"
            if not os.path.exists(fp):
                missing.append(pid)
                continue
            df = pd.read_csv(fp, sep="\t")
            df = df[df["tags"] >= TAGS_MIN]
            df = df.sort_values(by=["-10*log10(pvalue)"], ascending=False)
            df = df.iloc[:TOP]
            tmp = pd.concat([tmp, df], ignore_index=True)
        if missing:
            raise SystemExit("%s: peak files missing: %s" % (TF, missing))
        bed = find_peaks(tmp)
        n_pool = bed.shape[0]
        n_pool_unique = bed.drop_duplicates(subset=["#CHROM", "START", "END"]).shape[0]
        cleaned = repeats_cleaner(bed, a.mask)
        n_clean = cleaned.shape[0]
        tr = cleaned[[x in TRAIN_CHR for x in list(cleaned["chrom"])]]       # cell 58
        te = cleaned[[x in TEST_CHR for x in list(cleaned["chrom"])]]
        n_tr_pre, n_te_pre = tr.shape[0], te.shape[0]
        if tr.shape[0] > CAP:
            tr = tr.sample(n=CAP, random_state=1)
        if te.shape[0] > CAP:
            te = te.sample(n=CAP, random_state=1)
        tr = tr.copy(); te = te.copy()
        tr["name"] = tr[["name", "score"]].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
        te["name"] = te[["name", "score"]].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
        if FACTOR:                                                         # human cell 28
            tr["factor"] = TF
            te["factor"] = TF
        all_tr = pd.concat([all_tr, tr], ignore_index=True)                 # cell 61
        all_te = pd.concat([all_te, te], ignore_index=True)
        tr.to_csv(no_clobber(f"{out}/train/{TF}_train_set_{organism}.bed"), sep='\t')
        te.to_csv(no_clobber(f"{out}/test/{TF}_test_set_{organism}.bed"), sep='\t')
        pbt.BedTool.from_dataframe(tr).sequence(fi=GENOME, name=True).save_seqs(
            no_clobber(f"{out}/train/{TF}_{organism}.PEAKS000001.macs.train.mfa"))
        pbt.BedTool.from_dataframe(te).sequence(fi=GENOME, name=True).save_seqs(
            no_clobber(f"{out}/test/{TF}_{organism}.PEAKS000001.macs.test.mfa"))
        counts.append(dict(TF=TF, peak_sets=len(set(ids)), peak_set_reads=len(ids), pool_rows=n_pool,
                           pool_unique_windows=n_pool_unique, after_mask=n_clean,
                           removed_by_mask=n_pool - n_clean, train_chr1_8_19_before_cap=n_tr_pre,
                           control_before_cap=n_te_pre, train_file_pos=tr.shape[0], control_file_pos=te.shape[0]))
        print(counts[-1], flush=True)
        pbt.cleanup()
    all_tr.to_csv(no_clobber(f"{out}/all_bed_df_train_set_{organism}.bed"), sep='\t')
    all_te.to_csv(no_clobber(f"{out}/all_bed_df_test_set_{organism}.bed"), sep='\t')
    pd.DataFrame(counts).to_csv(no_clobber(f"{out}/per_tf_positive_counts.csv"), index=False)
    print("PER_TF DONE", flush=True)


def pooled(a):
    out, basicdir = a.out, a.out
    files_tr = sorted(glob.glob(os.path.join(a.mfa_dir, "train", f"*_{organism}.PEAKS000001.macs.train.mfa")))
    files_te = sorted(glob.glob(os.path.join(a.mfa_dir, "test", f"*_{organism}.PEAKS000001.macs.test.mfa")))
    selected = [os.path.basename(f).split(".")[0] for f in files_tr]
    fam = pd.read_csv(FAMILY, sep='\t')
    TF_family_dict = defaultdict(None)
    for i, j in zip(fam["Model"], fam["TF family"]):
        if type(j) == str:
            TF_family_dict[i.split('_')[0]] = j.split('{')[1][:3]
    names_all = {}
    for slim in (False, True):                                             # cell 41 (slim=False), cell 42 (True)
        tag = "SLIM_" if slim else ""
        for split, files in (("train", files_tr), ("control", files_te)):
            fa = no_clobber(f"{basicdir}/all_mfa_file_{organism}_{tag}{SIZE_OF_SET}_{split}.fasta")
            hcsv = no_clobber(f"{basicdir}/headerStr_names_{organism}_{tag}{SIZE_OF_SET}_{split}.csv")
            names = []
            index_global = 0
            with open(fa, "w") as fo:
                for i in sorted(selected):
                    tmp_file_list = []
                    for j in files:
                        if i in j:
                            headerStr_base = "_".join(j.split(".")[:-1]).split("/")[-1]
                            with open(j) as f:
                                index_local = 0
                                for headerStr, seq in fasta_iter(f):
                                    new = headerStr_base + "_globalInd_" + str(index_global) + "_localInd_" + str(index_local) + "_train"
                                    if slim:
                                        p = headerStr.split(":")
                                        new = new + "$" + f'>chrom: {p[3]}; center: {p[1]}; peak: 501; signal: {p[0]}'
                                        if headerStr_base.split("_")[1] == organism:
                                            tmp_file_list.append(">" + new + "\n" + seq + "\n")
                                        index_local += 1
                                    else:
                                        seq = seq[350:650]
                                        if len(seq) < 320 and len(seq) > 280:
                                            seq = seq[(len(seq) // 2) - 150:(len(seq) // 2) + 150]
                                            if headerStr_base.split("_")[1] == organism:
                                                tmp_file_list.append(">" + new + "\n" + seq + "\n")
                                            index_local += 1
                            index_global += 1
                    if len(tmp_file_list) == 0:
                        continue
                    if len(tmp_file_list) > SIZE_OF_SET:
                        raise SystemExit("more than %d sequences for %s (cap is applied upstream)" % (SIZE_OF_SET, i))
                    tmp_file_list = list(dict.fromkeys(tmp_file_list))
                    for k in tmp_file_list:
                        if len(k) != 0:
                            names.append(k[1:].split("\n")[0])
                            fo.write(k)
            with open(hcsv, "w") as output:
                w = csv.writer(output, lineterminator='\n')
                for v in names:
                    w.writerow([v])
            names_all[(tag, split)] = names
            print("pooled", tag or "300bp", split, len(names), flush=True)
    # cell 52 runs after cell 42 in the original kernel: out_tab is built from the SLIM names (with '$>chrom: ...')
    for split in ("train", "control"):
        n300 = [x for x in names_all[("", split)]]
        nsl = names_all[("SLIM_", split)]
        if [x.split("$")[0] for x in nsl] != n300:
            raise SystemExit("SLIM and 300 bp name lists differ for " + split)
        df = pd.DataFrame({'Names': nsl})
        fams = []
        for x in nsl:
            try:
                fams.append(TF_family_dict[x.split('_')[0]])
            except KeyError:
                fams.append(None)
        df["TF_family"] = np.array(fams)
        df["TF_name"] = np.array(["_".join(x.split("_")[:2]) for x in nsl])
        df.to_csv(no_clobber(f"{basicdir}/out_tab_{organism}_{SIZE_OF_SET}_{split}.tab"), sep='\t', header=False)
    for tag in ("", "SLIM_"):
        for split in ("train", "control"):
            src = f"{basicdir}/all_mfa_file_{organism}_{tag}{SIZE_OF_SET}_{split}.fasta"
            nf = no_clobber(src[:-6] + "_no_NF.fasta")
            nn = no_clobber(src[:-6] + "_no_NF_no_N.fasta")
            subprocess.run(f"awk 'NF' {src} > {nf}", shell=True, check=True)
            subprocess.run(f"sed '/^>/! s/[^ACTG]/A/g' {nf} > {nn}", shell=True, check=True)
    print("POOLED DONE", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["per_tf", "pooled"])
    ap.add_argument("--peaks"); ap.add_argument("--mask"); ap.add_argument("--out", required=True)
    ap.add_argument("--dup", action="store_true"); ap.add_argument("--tfs", nargs="*")
    ap.add_argument("--mfa-dir")
    ap.add_argument("--species", choices=["HUMAN", "MOUSE"], required=True)
    a = ap.parse_args()
    set_species(a.species)
    os.makedirs(a.out, exist_ok=True)
    per_tf(a) if a.stage == "per_tf" else pooled(a)
