"""h5_splits.py - the one definition of the chromosome split of the human and mouse data builds.

Logical sets: "train" = training chromosomes, "test" = test chromosomes. Every autosome is in exactly one set; the sex
chromosomes (and anything else, e.g. chrM, unplaced contigs) are in neither.

On-disk file names: HUMAN file "train" = training set, file "test"/"control" = test set. MOUSE: the file named "train"
holds the mouse TEST chromosomes (1, 8, 19) and the file named "test"/"control" the mouse TRAINING chromosomes;
on_disk() returns the lists in file-name order so that build code keeps that convention.
"""

AUTOSOMES = {"HUMAN": 22, "MOUSE": 19}
SEX = ("chrX", "chrY")
LOGICAL = {
    "HUMAN": {"train": ["chr%d" % i for i in range(1, 23) if i not in (1, 8, 21)], "test": ["chr1", "chr8", "chr21"]},
    "MOUSE": {"train": ["chr%d" % i for i in range(1, 20) if i not in (1, 8, 19)], "test": ["chr1", "chr8", "chr19"]},
}


def on_disk(species):
    """(chromosomes of the file named 'train', chromosomes of the file named 'test'/'control')."""
    s = LOGICAL[species]
    if species == "MOUSE":
        return list(s["test"]), list(s["train"])  # mouse file "train" = test chromosomes
    return list(s["train"]), list(s["test"])


def partition_problems(species, train, test):
    """List of problems of a (train, test) pair of chromosome lists: every autosome of `species` must be in exactly one
    list, and nothing else (sex chromosomes, typo strings such as 'chr1112') in either. Empty list = valid."""
    auto = {"chr%d" % i for i in range(1, AUTOSOMES[species] + 1)}
    tr, te = list(train), list(test)
    out = []
    for name, lst in (("train", tr), ("test", te)):
        dup = sorted({c for c in lst if lst.count(c) > 1})
        if dup:
            out.append("%s %s: duplicated %s" % (species, name, dup))
        extra = sorted(set(lst) - auto)
        if extra:
            out.append("%s %s: not an autosome %s" % (species, name, extra))
    both = sorted(set(tr) & set(te))
    if both:
        out.append("%s: in both sets %s" % (species, both))
    missing = sorted(auto - set(tr) - set(te), key=lambda c: int(c[3:]))
    if missing:
        out.append("%s: in neither set %s" % (species, missing))
    return out
