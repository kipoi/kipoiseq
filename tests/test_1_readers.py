from kipoi_dataloaders.utils import TsvReader, FastaReader
import pytest
from pybedtools import Interval
import numpy as np

@pytest.mark.parametrize("tsv_file", ["tests/data/sample_intervals.bed", "tests/data/sample_intervals_nochr.bed"])
@pytest.mark.parametrize("num_chr", [True, False])
@pytest.mark.parametrize("label_dtype", [str, np.int64])
def test_tsvreader(tsv_file, num_chr, label_dtype):
    reader = TsvReader(tsv_file, num_chr, label_dtype)
    interval, labels = reader[0]
    assert isinstance(interval, Interval)
    if not num_chr:
        assert interval.chrom.startswith("chr")
    assert isinstance(labels[0], label_dtype)
    assert interval.start == 2
    assert interval.end == 4


comp = {"A":"T", "C":"G", "G":"C", "T":"A"}
comp.update({k.lower(): v.lower() for k, v in comp.items()})

@pytest.mark.parametrize("use_strand", [True, False])
@pytest.mark.parametrize("force_upper", [True, False])
def test_fastareader(use_strand, force_upper):
    fp = "tests/data/sample.fasta"
    with open(fp, "r") as ifh:
        for i, s in enumerate(ifh):
            if i == 1:
                fasta_str = s.lstrip()
    fr = FastaReader(fp, use_strand, force_upper)
    test_intervals = [Interval("chr1", 0, 2, strand="-"), Interval("chr1", 3, 4)]
    seqs = fr(test_intervals)
    assert len(seqs) == len(test_intervals)
    for inter, seq in zip(test_intervals, seqs):
        ref_seq = fasta_str[inter.start:inter.end]
        if use_strand and inter.strand == "-":
            ref_seq = list(ref_seq)[::-1]
            ref_seq = "".join([comp[el] for el in ref_seq])
        if force_upper:
            assert seq == ref_seq.upper()
        else:
            assert seq == ref_seq




