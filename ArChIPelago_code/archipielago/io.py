"""FASTA/BED I/O helpers for ArChIPelago pipelines.

These functions are extracted from notebooks 0, 1, and 2 where fasta_iter
and related helpers were duplicated across multiple cells.
"""

from itertools import groupby
from pathlib import Path


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def fasta_iter(f):
    """Parse an open FASTA file handle into (header, sequence) tuples.

    Parameters
    ----------
    f : file-like
        Open text file handle positioned at the start of a FASTA file.

    Yields
    ------
    tuple[str, str]
        (header_without_gt, sequence_as_single_string)
    """
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
    for header in faiter:
        header_str = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header_str, seq


def load_fasta(path):
    """Read a FASTA file and return all (header, sequence) pairs as a list.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    list[tuple[str, str]]
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")
    with open(path, "r") as f:
        return list(fasta_iter(f))


def save_fasta(records, path):
    """Write (header, sequence) pairs to a FASTA file.

    Parameters
    ----------
    records : iterable[tuple[str, str]]
        (header, sequence) pairs.
    path : str or Path
    """
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n{seq}\n")


# ---------------------------------------------------------------------------
# BED / peak helpers
# ---------------------------------------------------------------------------

def make_train_test_beds(peaks_df, train_chr, test_chr, n_max=10000, random_state=1):
    """Split a peaks DataFrame into train and test by chromosome.

    Parameters
    ----------
    peaks_df : pd.DataFrame
        Must have a column named 'chrom' (or 'seqname') with chromosome names
        like 'chr1', 'chr2', etc.
    train_chr : list[str]
        Chromosomes to include in the training set.
    test_chr : list[str]
        Chromosomes to include in the test set.
    n_max : int
        Maximum number of peaks per split. Excess is randomly downsampled.
    random_state : int
        Random seed for reproducible downsampling.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        (train_set, test_set)
    """
    chrom_col = "chrom" if "chrom" in peaks_df.columns else "seqname"
    train_set = peaks_df[peaks_df[chrom_col].isin(train_chr)]
    test_set = peaks_df[peaks_df[chrom_col].isin(test_chr)]

    if len(train_set) > n_max:
        train_set = train_set.sample(n=n_max, random_state=random_state)
    if len(test_set) > n_max:
        test_set = test_set.sample(n=n_max, random_state=random_state)

    return train_set.reset_index(drop=True), test_set.reset_index(drop=True)


def extract_sequences(bed_df, genome_fasta, output_fasta):
    """Extract sequences from a genome FASTA using a BED DataFrame.

    Wraps pybedtools.BedTool.sequence() for use in the pipeline.

    Parameters
    ----------
    bed_df : pd.DataFrame
        BED-format DataFrame (chrom, start, end, ...).
    genome_fasta : str or Path
        Path to the indexed genome FASTA.
    output_fasta : str or Path
        Where to write the extracted sequences.

    Returns
    -------
    Path
        Path to the written output FASTA.
    """
    import pybedtools as pbt

    genome_fasta = str(genome_fasta)
    output_fasta = Path(output_fasta)

    peaks = pbt.BedTool.from_dataframe(bed_df)
    seqs = peaks.sequence(fi=genome_fasta, name=True)
    seqs.save_seqs(str(output_fasta))

    return output_fasta
