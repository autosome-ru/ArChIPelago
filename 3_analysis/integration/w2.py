"""Shared helpers of the integration scripts: backups, md5, provenance, TF checks (from p53drop.py)."""
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

B13 = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(B13))
sys.path.insert(0, os.path.join(R2, "3_analysis/integration"))
import p53drop as PD  # noqa: E402,F401

TAG = PD.BACKUP_TAG


def backup(path, tier=""):
    """p53drop.backup: <stem>_backup_<tag><tier><ext>, never overwritten."""
    return PD.backup(path, tier)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _git():
    """Commit of the git repository at the repo root and the number of changed tracked files; 'not a git repository' (with a printed warning) when there is none."""
    root = R2
    run = lambda *a: subprocess.run(["git", "-C", root] + list(a), capture_output=True, text=True)
    h = run("rev-parse", "HEAD")
    top = run("rev-parse", "--show-toplevel")
    if h.returncode or top.returncode or os.path.realpath(top.stdout.strip()) != os.path.realpath(root):
        print("provenance: NOT A GIT REPOSITORY at %s; git_commit is recorded as 'not a git repository'" % root)
        return "not a git repository", -1
    dirty = [l for l in run("status", "--porcelain", "--untracked-files=no").stdout.splitlines() if l.strip()]
    return h.stdout.strip(), len(dirty)


def provenance(folder, script, outputs, inputs, params=None, note="", external=None):
    """Append one block to <folder>/_provenance_h3.yaml. Every scalar is written as a JSON string (valid YAML), the block
    carries the git commit and the count of changed tracked files, and paths are relative to the repo root (absolute outside it).
    external: {local folder: label}, e.g. {scores_dir: "zenodo:scores_H5_human.tar.gz/scores"}, records files under that folder
    by label instead of their local path (the md5 is still taken from the local file)."""
    q = json.dumps
    def rel(f):
        a = os.path.abspath(f)
        for d, lab in (external or {}).items():
            if a.startswith(os.path.abspath(d) + os.sep):
                return lab + "/" + os.path.relpath(a, os.path.abspath(d))
        return os.path.relpath(f, R2) if a.startswith(R2 + os.sep) else a
    commit, dirty = _git()
    L = ["- script: %s" % q(script),
         "  run: %s" % q(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
         "  git_commit: %s" % q(commit),
         "  git_changed_tracked_files: %d" % dirty,
         "  task: %s" % q("H3 integration (final data: human B14 + H3a, mouse B11 + H3a, negatives v3)"),
         "  python: %s" % q(sys.version.split()[0]),
         "  exclude_tfs: [%s]" % ", ".join(q(t) for t in PD.EXCLUDE_TFS),
         "  human_used_ids: {file: %s, md5: %s}" % (q(rel(PD.HUMAN_IDS)), q(md5(PD.HUMAN_IDS))),
         "  mouse_used_ids: {file: %s, md5: %s}" % (q(rel(PD.MOUSE_IDS)), q(md5(PD.MOUSE_IDS)))]
    if params:
        L.append("  parameters: %s" % q(params))
    if note:
        L.append("  note: %s" % q(note))
    L.append("  inputs:" + ("" if inputs else " []"))
    L += ["    - {file: %s, md5: %s}" % (q(rel(f)), q(md5(f))) for f in inputs]
    L.append("  outputs:" + ("" if outputs else " []"))
    L += ["    - {file: %s, md5: %s}" % (q(rel(f)), q(md5(f))) for f in outputs]
    p = os.path.join(folder, "_provenance_h3.yaml")
    with open(p, "a") as fh:
        fh.write("\n".join(L) + "\n")
    print("provenance ->", os.path.relpath(p, R2))
