from kipoiseq.extractors import FastaStringExtractor
import pytest
import numpy as np
from pybedtools import Interval


comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
comp.update({k.lower(): v.lower() for k, v in comp.items()})


@pytest.mark.parametrize("use_strand", [True, False])
@pytest.mark.parametrize("force_upper", [True, False])
def test_fastareader(use_strand, force_upper):
    fp = "tests/data/sample.fasta"
    with open(fp, "r") as ifh:
        for i, s in enumerate(ifh):
            if i == 1:
                fasta_str = s.lstrip()
    fr = FastaStringExtractor(fp, use_strand, force_upper)
    intervals = Interval("chr1", 0, 2, strand="-"), Interval("chr1", 3, 4)

    for interval in intervals:
        seq = fr.extract(interval)
        ref_seq = fasta_str[interval.start:interval.end]
        if use_strand and interval.strand == "-":
            ref_seq = list(ref_seq)[::-1]
            ref_seq = "".join([comp[el] for el in ref_seq])
        if force_upper:
            assert seq == ref_seq.upper()
        else:
            assert seq == ref_seq
