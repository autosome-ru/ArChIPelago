"""
make_scan_jobs_sym.py - SARUS job list for the <SP>_seq_{HUMAN,MOUSE}_pwm_{mono,di} feature tables of one sequence
species (--species HUMAN|MOUSE).
The PWM lists, their k indices and the output file names are taken from the published feature folders
~/archi_release/TF-ML/outputdir/<TF>/<folder>/ (k.pwm / k.dpwm copies and <k>_<name>_feature_compare_table_<split>_cut.tab),
so the feature indices stay equal to the human tables. Each k.pwm is checked to be byte-identical to the HOCOMOCO
matrix it is named after (~/hocomoco11/models/pwm/{mono,di}/all/<TF_SP>/<name>.pwm|dpwm).
Command = release notebooks (MONO/DI cells, sarus-2.0.2.jar), plus single-CPU JVM flags (same output), piped through
the original awk truncation to two decimals.
usage: python make_scan_jobs_sym.py --tfml TFMLDIR --tfs A B ... --species HUMAN|MOUSE > jobs.txt
"""
import os, sys, glob, argparse, filecmp

H = os.path.expanduser("~")
ORIG = f"{H}/Release/TF-ML/outputdir"
HOCO = {"mono": f"{H}/hocomoco11/models/pwm/mono/all", "di": f"{H}/hocomoco11/models/pwm/di/all"}
JAR = f"{H}/sarus-2.0.2.jar"
CLS = {"mono": "ru.autosome.SARUS", "di": "ru.autosome.di.SARUS"}
EXT = {"mono": "pwm", "di": "dpwm"}

ap = argparse.ArgumentParser()
ap.add_argument("--tfml", required=True); ap.add_argument("--tfs", nargs="+", required=True)
ap.add_argument("--species", choices=["HUMAN", "MOUSE"], required=True)
a = ap.parse_args()
n = 0
summary = []
for tf in a.tfs:
    for pwm_sp in ("HUMAN", "MOUSE"):
        for kind in ("mono", "di"):
            folder = f"{a.species}_seq_{pwm_sp}_pwm_{kind}"
            od = os.path.join(ORIG, tf, folder)
            pre = {}
            for split in ("train", "control"):
                fs = glob.glob(os.path.join(od, f"*_feature_compare_table_{split}_cut.tab"))
                pre[split] = sorted(os.path.basename(f).split("_feature_compare_table_")[0] for f in fs)
            if not pre["train"] or pre["train"] != pre["control"]:
                raise SystemExit(f"{tf} {folder}: train/control PWM lists differ or empty")
            nd = os.path.join(a.tfml, "outputdir", tf, folder)
            os.makedirs(nd, exist_ok=True)
            for p in pre["train"]:
                k, name = p.split("_", 1)
                pwm = os.path.join(od, f"{k}.{EXT[kind]}")
                sp_dir = name.split("~")[0]                                    # e.g. SRF_HUMAN
                hoco_name = name.split("~")[0] + "~" + name.split("~")[1] + "~" + ".".join(name.split("~")[2].rsplit("_", 3)) + "." + EXT[kind]
                hoco = os.path.join(HOCO[kind], sp_dir, hoco_name)
                if not os.path.isfile(hoco) or not filecmp.cmp(pwm, hoco, shallow=False):
                    raise SystemExit(f"{tf} {folder} {k}: {pwm} is not identical to {hoco}")
                for split in ("train", "control"):
                    fa = os.path.join(a.tfml, f"all_mfa_file_{a.species}_10000_{split}_no_NF.fasta")
                    out = os.path.join(nd, f"{p}_feature_compare_table_{split}_cut.tab")
                    if os.path.exists(out):
                        raise SystemExit("refusing to overwrite " + out)
                    print(f"java -XX:ActiveProcessorCount=1 -XX:ParallelGCThreads=1 -Xmx2G -cp {JAR} {CLS[kind]} {fa} {pwm} "
                          f"--skipn --show-non-matching --output-scoring-mode score besthit | grep -v '>' | "
                          f"awk -F. '{{print $1\".\"substr($2,1,2)}}' > {out}")
                    n += 1
            summary.append(f"{tf}\t{folder}\t{len(pre['train'])}")
with open(os.path.join(a.tfml, "n_pwms_per_folder.tsv"), "x") as f:
    f.write("TF\tfolder\tn_pwm\n" + "\n".join(summary) + "\n")
print(f"{n} jobs", file=sys.stderr)
