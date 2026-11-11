"""Unpack the human PWM archives of the Zenodo deposit v1.0.0 next to the archives.

The four PWM archives of the deposit and its md5sums.txt are tracked in results/datasets/zenodo_v1.0.0/. The local producers
read the PWM files from $ARCHI_ZENODO_DIR/PWMs_mono_HUMAN/<TF>/*.pwm and $ARCHI_ZENODO_DIR/PWMs_di_HUMAN/<TF>/*.dpwm, and
PWMs_mono_MOUSE.tar.gz directly. Every archive is checked against md5sums.txt first; a folder that is already unpacked is kept.
The unpacked folders are not tracked (results/datasets/zenodo_v1.0.0/.gitignore). run_h3.py calls unpack() before any stage.
usage: python3 3_analysis/integration/zenodo_v1.py [folder]      (default: results/datasets/zenodo_v1.0.0)
"""
import hashlib
import os
import sys
import tarfile

DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                       "results", "datasets", "zenodo_v1.0.0")
ARCHIVES = ("PWMs_mono_HUMAN", "PWMs_di_HUMAN", "PWMs_mono_MOUSE", "PWMs_di_MOUSE")
UNPACK = ("PWMs_mono_HUMAN", "PWMs_di_HUMAN")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def unpack(folder=DEFAULT):
    sums = dict(l.split()[::-1] for l in open(os.path.join(folder, "md5sums.txt")) if l.strip())
    for name in ARCHIVES:
        tar = os.path.join(folder, name + ".tar.gz")
        if sums.get(name + ".tar.gz") != md5(tar):
            raise SystemExit("CHECK FAILED %s: md5 differs from md5sums.txt of the Zenodo deposit" % tar)
    for name in UNPACK:
        dst = os.path.join(folder, name)
        if os.path.isdir(dst):
            continue
        with tarfile.open(os.path.join(folder, name + ".tar.gz")) as t:
            bad = [m.name for m in t.getmembers() if not (m.name == name or m.name.startswith(name + "/")) or ".." in m.name.split("/")
                   or not (m.isfile() or m.isdir())]
            if bad:
                raise SystemExit("CHECK FAILED %s.tar.gz: unexpected members %s" % (name, bad[:3]))
            t.extractall(folder)
        print("unpacked %s.tar.gz -> %s" % (name, dst))
    return folder


if __name__ == "__main__":
    unpack(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
