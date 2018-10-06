import pytest
import numpy as np
import copy
from kipoiseq.transforms.functional import resize_interval
from kipoiseq.transforms import TransformShape
from kipoiseq.utils import DNA
from pybedtools import Interval


# --------------------------------------------

@pytest.mark.parametrize("alphabet_axis", list(range(0, 5)))
@pytest.mark.parametrize("dummy_axis", [None] + list(range(0, 5)))
def test_reshape_seq(alphabet_axis, dummy_axis):
    n_samples, seq_len, alphabet_len = 20, 1000, 4
    in_array = np.zeros((n_samples, seq_len, alphabet_len))

    dummy_axis_int = dummy_axis
    if dummy_axis is None:
        dummy_axis_int = -1

    if (alphabet_axis == dummy_axis_int) or (alphabet_axis == 0) or (dummy_axis_int == 0) or \
            (alphabet_axis >= 4) or (dummy_axis_int >= 4) or ((alphabet_axis >= 3) and (dummy_axis is None)):
        with pytest.raises(Exception):
            reshaper = TransformShape(alphabet_axis, dummy_axis)
        return None

    reshaper = TransformShape(alphabet_axis, dummy_axis)
    num_axes = 3
    if dummy_axis is not None:
        num_axes += 1

    # Test batch reshaping
    reshaped = reshaper.reshape_batch(in_array)
    assert len(reshaped.shape) == num_axes
    for i in range(len(reshaped.shape)):
        if i == 0:
            assert reshaped.shape[i] == n_samples
        elif i == dummy_axis:
            assert reshaped.shape[i] == 1
        elif i == alphabet_axis:
            assert reshaped.shape[i] == alphabet_len
        else:
            assert reshaped.shape[i] == seq_len

    # test the single sample works
    reshaped = reshaper.reshape_single_sample(in_array[0, ...])
    for i in range(len(reshaped.shape)):
        i2 = i + 1
        if i2 == 0:
            assert reshaped.shape[i] == n_samples
        elif i2 == dummy_axis:
            assert reshaped.shape[i] == 1
        elif i2 == alphabet_axis:
            assert reshaped.shape[i] == alphabet_len
        else:
            assert reshaped.shape[i] == seq_len

    # Test if fails if input has wrong shape
    with pytest.raises(Exception):
        reshaper.reshape_batch(in_array[0, ...])


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
