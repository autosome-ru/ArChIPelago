"""Fig. S7 co-factor check against the HOCOMOCO v11 TFs the text names.
For every Fig. S7 matrix, r (pcc_max of the S7 builder, overlap >= 5) to the best HOCOMOCO v11 HUMAN mono model of each named
reference TF (ETV2, ETV7, BATF, FLI1 by default; --tfs to change), with the model's quality letter from
HOCOMOCOv11_full_annotation_HUMAN_mono.tsv, and the co-factor gap r_ref - r_own (r_own = r_best_own of subtype_vs_hocomoco11.csv).
Output: subtype_vs_named_v11.csv (one row per Fig. S7 matrix), _provenance_named.yaml.
usage: python3 4_visualisation/supplementary_figures/Figure_S7/hocomoco11_check/subtype_vs_named_v11.py [--tfs ETV2,ETV7,BATF,FLI1]"""
import argparse, os, datetime, subprocess
import pandas as pd

from subtype_vs_hocomoco11 import HERE, S7, R2, PCM, SRC, MIN_OVERLAP, read_pcms, read_matrix, md5, pcc_max, zrows

ANN = os.path.join(R2, "results/datasets/hocomoco11/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv")
OWN = os.path.join(HERE, "subtype_vs_hocomoco11.csv")
REF_TFS = "ETV2,ETV7,BATF,FLI1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", default=REF_TFS)
    tfs = ap.parse_args().tfs.split(",")
    v11 = read_pcms(PCM)
    qual = pd.read_csv(ANN, sep="\t").set_index("Model")["Quality"]
    if set(v11) != set(qual.index):
        raise SystemExit("pcm and annotation model sets differ")
    ref = {tf: {k: zrows(m) for k, m in v11.items() if k.startswith(tf + "_HUMAN.")} for tf in tfs}
    for tf, d in ref.items():
        if not d:
            raise SystemExit("no HOCOMOCO v11 model for %s" % tf)
    src = pd.read_csv(SRC)
    own = pd.read_csv(OWN).set_index(["TF", "panel", "pwm_name"])
    rows = []
    for _, r in src.iterrows():
        za = zrows(read_matrix(os.path.join(S7, r.matrix_file)))
        o = own.loc[(r.TF, r.panel, r.pwm_name)]
        row = dict(TF=r.TF, panel=r.panel, pwm_class=r.pwm_class, pwm_name=r.pwm_name, cell_line=r.cell_line,
                   best_own_v11=o.best_own_v11, quality_own=qual[o.best_own_v11], r_own=o.r_best_own)
        for tf, d in ref.items():
            rr = {k: pcc_max(za, z, MIN_OVERLAP)[0] for k, z in d.items()}
            b = max(rr, key=rr.get)
            row.update({"best_%s" % tf: b, "quality_%s" % tf: qual[b], "r_%s" % tf: round(rr[b], 3),
                        "gap_%s" % tf: round(round(rr[b], 3) - o.r_best_own, 3)})
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "subtype_vs_named_v11.csv"), index=False)
    git = subprocess.run(["git", "-C", R2, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    with open(os.path.join(HERE, "_provenance_named.yaml"), "w") as fh:
        fh.write("script: subtype_vs_named_v11.py\ndate: %s\ngit: %s\nmin_overlap: %d\nref_tfs: %s\n"
                 "gap: r_ref - r_own (both rounded to 3 decimals first)\ninputs:\n" % (
                     datetime.datetime.now().isoformat(timespec="seconds"), git, MIN_OVERLAP, ",".join(tfs)))
        for p in (PCM, ANN, SRC, OWN):
            fh.write("  - {path: %s, md5: %s}\n" % (os.path.relpath(p, R2), md5(p)))
    print(out[out.TF == "TFE2"].drop(columns=["pwm_name"]).to_string())


if __name__ == "__main__":
    main()
