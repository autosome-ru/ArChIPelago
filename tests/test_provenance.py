"""Every _provenance_h3.yaml of the repository parses as YAML, every block carries a git commit, and the last
block that wrote each output records the md5 of its inputs and outputs as they are on disk now (a consumer that ran before its input was
rewritten is caught here). Files outside the repository that are absent (e.g. the md5-checked scratch copy of the server score npz)
are listed, not failed; a missing file inside the tree fails. Records written "zenodo:<archive>/<member>" (the score npz of the
Zenodo v2.0.0 deposit, not tracked) are listed as external and never looked up in the tree.
Red: T1_MUTATE=prov_md5 changes one recorded input md5; T1_MUTATE=prov_yaml appends an unquoted colon value to one file's text.
"""
import glob
import os

import pytest
import yaml

import t1_helpers as T

pytestmark = pytest.mark.r2


def prov_files():
    fs = sorted(f for f in glob.glob(os.path.join(T.R2, "**", "_provenance_h3.yaml"), recursive=True)
                if "/tests/" not in f and "_backup_" not in os.path.relpath(f, T.R2) and "/superseded" not in f)
    assert len(fs) >= 8, fs
    return fs


def load(f):
    txt = open(f).read()
    if T.mutate("prov_yaml") and f.endswith(os.path.join("3_analysis/integration", "_provenance_h3.yaml")):
        txt += "- script: x\n  note: a: b: c\n"
    return yaml.safe_load(txt)


def path(p):
    return p if os.path.isabs(p) else os.path.join(T.R2, p)


def test_every_provenance_file_is_valid_yaml_with_git_commit():
    bad = []
    for f in prov_files():
        try:
            d = load(f)
        except yaml.YAMLError as e:
            bad.append("%s: %s" % (os.path.relpath(f, T.R2), str(e).splitlines()[0]))
            continue
        if not isinstance(d, list) or not all(isinstance(b, dict) and len(str(b.get("git_commit", ""))) == 40 for b in d):
            bad.append("%s: not a list of blocks with a 40-character git_commit" % os.path.relpath(f, T.R2))
    assert not bad, bad


def test_recorded_md5_equal_the_files_on_disk():
    import hashlib

    def md5(p):
        h = hashlib.md5()
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()

    bad, external = [], []
    for f in prov_files():
        blocks = load(f)
        writer = {}   # output file -> index of the last block that wrote it; only those blocks describe the files on disk
        for n, b in enumerate(blocks):
            for e in b.get("outputs") or []:
                writer[e["file"]] = n
        current = sorted(set(writer.values()))
        for i, n in enumerate(current):
            b, script = blocks[n], blocks[n]["script"]
            entries = [("input", e) for e in b.get("inputs") or []] + [("output", e) for e in b.get("outputs") or [] if writer[e["file"]] == n]
            for k, (kind, e) in enumerate(entries):
                want = e["md5"]
                if T.mutate("prov_md5") and i == 0 and k == 0 and f == prov_files()[0]:
                    want = "0" * 32
                if e["file"].startswith("zenodo:"):
                    external.append("%s: %s %s (Zenodo archive)" % (script, kind, e["file"]))
                    continue
                p = path(e["file"])
                if not os.path.exists(p):
                    (external if not os.path.abspath(p).startswith(T.R2 + os.sep) else bad).append("%s: %s %s missing" % (script, kind, e["file"]))
                elif md5(p) != want:
                    bad.append("%s (%s): %s %s md5 differs" % (os.path.relpath(f, T.R2), script, kind, e["file"]))
    print("external files not present: %d" % len(external))
    assert not bad, bad[:20]


def test_h3b_fig4_mirror_equals_the_server_record():
    """The three Fig. 4 / S4 files mirrored from the server (fig4_refit_trainscale) have the md5 that the
    mirror manifest records for them. Red: T1_MUTATE=prov_h3b changes the recorded md5 of the Fig. 4 source data."""
    import hashlib
    b16 = os.path.join(T.R2, T.h3in().ROOT)   # mirror of the stage (results/server_h4_test_sets or results/server)
    d = os.path.join(b16, "human", "fig4_refit_trainscale")
    rec = {}
    for line in open(os.path.join(b16, "MIRROR_MANIFEST_md5.txt")):
        m, p = line.split(None, 1)
        if p.strip().startswith("human/fig4_refit_trainscale/"):
            rec[os.path.basename(p.strip())] = m
    if T.mutate("prov_h3b"):
        rec["Figure_4_source_data.csv"] = "0" * 32
    assert sorted(rec) == ["Figure_4_and_S4_numbers.csv", "Figure_4_source_data.csv", "Figure_S4_source_data.csv"]
    bad = [f for f, m in rec.items() if hashlib.md5(open(os.path.join(d, f), "rb").read()).hexdigest() != m]
    assert not bad, bad
