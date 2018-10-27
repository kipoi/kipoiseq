import pytest
import numpy as np
import copy
from kipoiseq.transforms.transforms import Compose, OneHot, SplitSplicingSeq, ReorderedOneHot
from kipoiseq.utils import DNA
from pybedtools import Interval


# --------------------------------------------
def test_compose():
    c = Compose([OneHot()])
    print(str(c))
    assert c("ACGT").shape == (4, 4)


def test_ReorderedOneHot():
    seqlen = 10
    seq = 'A' * seqlen

    test_pairs = [
        (('ACGT', 1, None), (10, 4)),
        (('ACGT', 2, 1), (10, 1, 4)),
        (('ACGTD', 2, 1), (10, 1, 5)),
        (('ACGT', 0, None), (4, 10)),
        (('ACGT', 0, 1), (4, 1, 10)),
        (('ACGT', 0, 2), (4, 10, 1)),
    ]

    for args, result in test_pairs:
        tr = ReorderedOneHot(alphabet=args[0], alphabet_axis=args[1], dummy_axis=args[2])
        out = tr(seq)
        assert out.shape == tr.get_output_shape(seqlen)
        assert out.shape == result

    with pytest.raises(ValueError):
        ReorderedOneHot(alphabet_axis=1, dummy_axis=1)

    with pytest.raises(ValueError):
        ReorderedOneHot(dummy_axis=1)


def test_SplitSplicingSeq():
    split = SplitSplicingSeq(exon_cut_l=0,
                             exon_cut_r=0,
                             intron5prime_cut=3,
                             intron3prime_cut=3,
                             acceptor_intron_len=2,
                             acceptor_exon_len=3,
                             donor_exon_len=3,
                             donor_intron_len=2
                             )

    # seq = 'TAAAG GTAGTAGA GTCCC'
    seq = 'TAAAGGTAGTAGAGTCCC'
    splited = split(seq, 5, 5)
    assert splited['intron5prime'] == 'TA'
    assert splited['acceptor'] == 'AGGTA'
    assert splited['exon'] == 'GTAGTAGA'
    assert splited['donor'] == 'AGAGT'
    assert splited['intron3prime'] == 'CC'


def test_ResizeInterval():
    pass
