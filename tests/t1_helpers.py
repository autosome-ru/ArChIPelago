"""Shared helpers of the test suite.

Paths (all overridable by environment variables, which is how the mutation runs point the suite at a mutated copy):
  T1_R2      repository root (local tables, NUMBERS file); default: the folder that holds tests/, if it has
             outputs/supplementary_tables/, otherwise unset (the tests that need it are skipped)
  T1_CODE_H  human_sym code (build_positives_sym.py, make_negatives_sym.py, archi_data.py, ...)
             default: $ARCHI_WORK/human_sym/code when the real data are present, else
             the earlier build code (not in this repository) under <T1_R2>/B14_human_sym/code
  T1_CODE_M  mouse_sym code (build_positives.py, make_negatives_v2.py, archi_data.py, ...), same rule with mouse_sym and
             the earlier build code (not in this repository) under <T1_R2>/B11_mouse_sym/code
  T1_MUTATE  name of an in-memory fault to inject into real data (comma separated; used only by the red runs of the real-data tests)
Python 3.8 compatible.
  ARCHI_WORK the working root; the real-data tests run only when it is set and $ARCHI_WORK/human_sym/TF-ML exists.
"""
import csv
import os
import random
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK_ROOT = os.environ.get("ARCHI_WORK")
REAL_DATA = bool(WORK_ROOT) and os.path.isdir(os.path.join(WORK_ROOT, "human_sym", "TF-ML"))


def _r2_default():
    cand = os.path.dirname(HERE)
    return cand if os.path.isdir(os.path.join(cand, "outputs/supplementary_tables")) else None


R2 = os.environ.get("T1_R2") or _r2_default()


def h3in():
    """The input registry module 3_analysis/integration/h3_inputs.py of the tree under test (its IN dict maps key to path)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("h3_inputs_t1", os.path.join(R2, "3_analysis", "integration", "h3_inputs.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


if REAL_DATA:
    _ch, _cm = os.path.join(WORK_ROOT, "human_sym", "code"), os.path.join(WORK_ROOT, "mouse_sym", "code")
else:
    _ch = os.path.join(R2, "B14_human_sym", "code") if R2 else None
    _cm = os.path.join(R2, "B11_mouse_sym", "code") if R2 else None
CODE_H = os.environ.get("T1_CODE_H") or _ch
CODE_M = os.environ.get("T1_CODE_M") or _cm

# the 35 modelled TFs (P53 dropped); the TFs of results/datasets/peak_ids/mouse_sym_peak_ids.csv
TF35 = ['ANDR', 'AP2A', 'CEBPB', 'COE1', 'CTCF', 'E2F4', 'ERG', 'ESR1', 'FLI1', 'GATA1', 'GATA2', 'GATA3', 'GCR', 'HNF4A',
        'IRF1', 'IRF4', 'JUND', 'MAFK', 'MAX', 'MYC', 'PPARG', 'PRGR', 'REST', 'RUNX1', 'RXRA', 'SOX2', 'SPI1', 'SRF', 'STA5A',
        'STAT1', 'STAT3', 'TAL1', 'TF65', 'TFE2', 'USF2']

CHROMS = {"HUMAN": ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]],
          "MOUSE": ["chr%s" % c for c in list(range(1, 20)) + ["X", "Y"]]}
CHROM_LEN = 60000
# expected split, from Methods; the on-disk file names per species
SPLIT = {"HUMAN": {"train_file": ["chr%d" % c for c in (2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20)],
                   "test_file": ["chr1", "chr8", "chr21"]},
         "MOUSE": {"train_file": ["chr1", "chr8", "chr19"],                     # the mouse test set lives in the *train* files
                   "test_file": ["chr%d" % c for c in (2, 3, 4, 5, 6, 7, 9, 10, 13, 14, 15, 16, 17, 18)]}}
INTERVAL_HEADER = ["#CHROM", "START", "END", "summit", "tags", "-10*log10(pvalue)", "fold_enrichment", "FDR(%)"]


def mutate(name):
    """True when the red run asks for fault `name` (T1_MUTATE, comma separated)."""
    return name in os.environ.get("T1_MUTATE", "").split(",")


# ------------------------------------------------------------------ synthetic input folders for the build scripts
def make_genome(path, species, seed=7):
    rng = random.Random(seed)
    seqs = {}
    with open(path, "w") as f:
        for c in CHROMS[species]:
            s = "".join(rng.choice("ACGT") for _ in range(CHROM_LEN))
            seqs[c] = s
            f.write(">%s\n" % c)
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    return seqs


def write_interval(path, rows):
    """rows: dicts with chrom (without 'chr'), START, summit, tags, p."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(INTERVAL_HEADER)
        for r in rows:
            w.writerow([r["chrom"], r["START"], r["START"] + 2 * r["summit"], r["summit"], r["tags"], r["p"], 5.0, 1.0])


FAMILY_ROWS = [  # Model, TF, TF family (TFClass code in braces; the code keeps the first 3 characters = class level)
    ("AAA_HUMAN.H11MO.0.A", "AAA", "Fam one{2.1.4}"),
    ("BBB_HUMAN.H11MO.0.A", "BBB", "Fam two, same class{2.1.3}"),
    ("CCC_HUMAN.H11MO.0.A", "CCC", "Other class{3.5.2}"),
]


def make_home(tmp, species="HUMAN"):
    """A fake folder tree with the inputs of the build scripts (ARCHI_MACS, ARCHI_GENOMES, ARCHI_FAMILY, ARCHI_SEQTK; see env_for)."""
    home = os.path.join(str(tmp), "home")
    for d in ("macs", "Genomic_data", "TF-ML_last_v_22", "seqtk"):
        os.makedirs(os.path.join(home, d), exist_ok=True)
    write_interval(os.path.join(home, "macs", "PEAKS041820.interval"), [])       # TEMPLATE (header only is read)
    with open(os.path.join(home, "TF-ML_last_v_22", "core_annotation_HUMAN_mono.csv"), "w") as f:
        f.write("Model\tTranscription factor\tTF family\n")
        for r in FAMILY_ROWS:
            f.write("\t".join(r) + "\n")
    gname = {"HUMAN": "hg38_genome.fasta", "MOUSE": "GRCm38.p6.genome.fa"}[species]
    seqs = make_genome(os.path.join(home, "Genomic_data", gname), species)
    return home, seqs


def env_for(home):
    e = dict(os.environ, HOME=home, ARCHI_MACS=os.path.join(home, "macs"), ARCHI_GENOMES=os.path.join(home, "Genomic_data"),
             ARCHI_FAMILY=os.path.join(home, "TF-ML_last_v_22"), ARCHI_SEQTK=os.path.join(home, "seqtk", "seqtk"), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
             PYTHONDONTWRITEBYTECODE="1")
    return e


def run(cmd, env=None, cwd=None, check=True):
    p = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=cwd)
    if check and p.returncode != 0:
        raise AssertionError("command failed rc=%d: %s\nSTDOUT:\n%s\nSTDERR:\n%s" % (p.returncode, " ".join(cmd), p.stdout[-3000:], p.stderr[-3000:]))
    return p


def build_positives(home, species, peaks_csv, mask, out, script=None, extra=()):
    """per_tf stage of build_positives_sym.py (or of the mouse-only build_positives.py when script is given)."""
    if script is None:
        cmd = [sys.executable, os.path.join(CODE_H, "build_positives_sym.py"), "per_tf", "--species", species]
    else:
        cmd = [sys.executable, script, "per_tf"]
    cmd += ["--peaks", peaks_csv, "--mask", mask, "--out", out] + list(extra)
    return run(cmd, env=env_for(home))


def pooled(home, species, out, script=None):
    if script is None:
        cmd = [sys.executable, os.path.join(CODE_H, "build_positives_sym.py"), "pooled", "--species", species]
    else:
        cmd = [sys.executable, script, "pooled"]
    return run(cmd + ["--mfa-dir", out, "--out", out], env=env_for(home))


def read_bed(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def read_fasta(path):
    out, name = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
            elif name is not None:
                out.append((name, line))
                name = None
    return out


def peaks_csv(path, pairs, n_motif=None):
    d = pd.DataFrame(pairs, columns=["TF", "peak_id"])
    if n_motif is not None:
        d["n_motif_files"] = n_motif
    d.to_csv(path, index=False)
    return path


def default_top(rows, top=1000, tags_min=10, kind=None):
    """Reference selection: tags >= 10, sort by -10log10 P descending with pandas' default sort (as the notebook), top 1,000."""
    d = pd.DataFrame(rows)
    d = d[d["tags"] >= tags_min]
    d = d.sort_values(by=["p"], ascending=False) if kind is None else d.sort_values(by=["p"], ascending=False, kind=kind)
    return d.iloc[:top]
