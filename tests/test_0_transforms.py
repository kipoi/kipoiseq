import pytest
import numpy as np
import copy
from kipoiseq.transforms.functional import resize_interval
from kipoiseq.transforms.transforms import SplitSplicingSeq, ReorderedOneHot
from kipoiseq.utils import DNA
from pybedtools import Interval


# --------------------------------------------

@pytest.mark.parametrize("anchor", ['start', 'end', 'center'])
@pytest.mark.parametrize("ilen", [3, 4])
def test_resize_interval(anchor, ilen):
    import pybedtools
    dummy_start, dummy_end = 10, 20
    dummy_centre = int((dummy_start + dummy_end) / 2)

    dummy_inter = pybedtools.create_interval_from_list(['chr2', dummy_start, dummy_end, 'intname'])
    ret_inter = resize_interval(dummy_inter, ilen, anchor)

    # the original interval was left intact
    assert dummy_inter.chrom == 'chr2'
    assert dummy_inter.start == dummy_start
    assert dummy_inter.end == dummy_end
    assert dummy_inter.name == 'intname'

    # metadata kept
    assert ret_inter.chrom == dummy_inter.chrom
    assert ret_inter.name == 'intname'

    # desired output width
    assert ret_inter.length == ilen

    # correct anchor point
    if anchor == "start":
        assert ret_inter.start == dummy_start
    elif anchor == "end":
        assert ret_inter.end == dummy_end
    elif anchor == "centre":
        assert int((ret_inter.start + ret_inter.end) / 2) == dummy_centre


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

    with pytest.raises(Exception):
        ReorderedOneHot(alphabet_axis=1, dummy_axis=1)

    with pytest.raises(Exception):
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
