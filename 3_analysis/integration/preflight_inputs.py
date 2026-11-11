"""Pre-flight of the integration inputs. For every key of the registry h3_inputs.IN (stage from
ARCHI_STAGE): its class (staged = mirror of the stage, data = DATA_H5 / H4_ONLY_DATA, not_in_<stage> = read from the
mirror of the previous training run on purpose, fixed = stage independent), the path, present / missing, and for a file under a mirror folder its md5
against that mirror's MIRROR_MANIFEST_md5.txt (match / MISMATCH / not in manifest). A folder counts as present when it holds at
least one file; its files are md5-checked one by one. Also checks that the GC check file carries the split labels of the stage.
Fails (exit 1) on any missing path, any md5 mismatch, a staged file absent from its manifest, or a missing split label.
Output: the table on stdout and, with --out FILE, as csv.
usage: ARCHI_STAGE=H5 python3 3_analysis/integration/preflight_inputs.py [--out preflight_H5.csv]"""
import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h3_inputs as H  # noqa: E402


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def manifests():
    """{mirror folder (relative to R2): {path inside the mirror: md5}} for every stage root with a manifest."""
    out = {}
    for root in sorted(set(H.STAGES.values())):
        p = os.path.join(H.R2, root, "MIRROR_MANIFEST_md5.txt")
        if os.path.exists(p):
            out[root] = dict(reversed(l.split(None, 1)) for l in open(p).read().splitlines() if l.strip())
            out[root] = {k.strip(): v for k, v in out[root].items()}
    return out


def klass(k):
    if k in H.FIXED:
        return "fixed"
    if H.STAGE == "H5" and k in H.NOT_IN_H5:
        return "not_in_H5"
    if H.STAGE == "H4" and k in H.NOT_IN_H4:
        return "not_in_H4"
    if k in H.DATA_H5:
        return "data"
    return "staged"


def rows():
    man = manifests()
    out = []
    for k, path in sorted(H.IN.items()):
        rel = os.path.relpath(path, H.R2)
        root = next((r for r in H.STAGES.values() if rel.startswith(r + os.sep)), rel.split(os.sep, 1)[0])
        cl = klass(k)
        files = []
        if os.path.isdir(path):
            files = [os.path.join(d, f) for d, _, fs in os.walk(path) if "__pycache__" not in d for f in fs]   # local byte code
            present = "present (%d files)" % len(files) if files else "MISSING (empty folder)"
        elif os.path.exists(path):
            files, present = [path], "present"
        else:
            present = "MISSING"
        if files == [os.path.join(H.R2, root, "MIRROR_MANIFEST_md5.txt")]:
            mstat = "the manifest itself"
        elif root in man and files:
            inner = [os.path.relpath(f, os.path.join(H.R2, root)).replace(os.sep, "/") for f in files]
            got = [man[root].get(i) for i in inner]
            if any(g is None for g in got):
                mstat = "NOT IN MANIFEST (%d of %d files)" % (sum(g is None for g in got), len(got))
            else:
                bad = [i for i, f, g in zip(inner, files, got) if md5(f) != g]
                mstat = "MISMATCH %s" % bad[:3] if bad else "md5 match (%d)" % len(files)
        else:
            mstat = ("NO MANIFEST for %s" % root) if root in H.STAGES.values() and root not in man else (
                "outside the mirrors" if root not in man else "-")
        ok = present.startswith("present") and not mstat.startswith(("MISMATCH", "NOT IN MANIFEST", "NO MANIFEST"))
        out.append(dict(key=k, cls=cl, path=rel, present=present, manifest=mstat, ok=ok))
    # split labels of the stage in the GC check (Sup. Table 2 GC columns are keyed by them)
    gc = H.IN["gc_v3"]
    if os.path.exists(gc):
        import csv
        labs = {(r["pool"], r["chromosomes"]) for r in csv.DictReader(open(gc))}
        want = {(H.GC_POOL[k[0]], v) for k, v in H.SPLIT_GC.items()}   # (pool of the stage, split label)
        miss = sorted(want - labs)
        out.append(dict(key="split labels in gc_v3", cls="check", path=os.path.relpath(gc, H.R2),
                        present="labels %s" % sorted(H.SPLIT_GC.values()), manifest="MISSING %s" % miss if miss else "all present", ok=not miss))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    rs = rows()
    w = max(len(r["key"]) for r in rs)
    print("pre-flight, stage %s (root %s): %d inputs" % (H.STAGE, H.ROOT, len(rs)))
    for r in rs:
        print("%s  %-*s  %-10s %-22s %-34s %s" % ("ok  " if r["ok"] else "FAIL", w, r["key"], r["cls"], r["present"][:22], r["manifest"][:34], r["path"]))
    nbad = sum(not r["ok"] for r in rs)
    print("pre-flight %s: %d of %d inputs fail" % (H.STAGE, nbad, len(rs)))
    if a.out:
        import csv
        with open(a.out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rs[0]))
            wr.writeheader()
            wr.writerows(rs)
    return 1 if nbad else 0


if __name__ == "__main__":
    sys.exit(main())
