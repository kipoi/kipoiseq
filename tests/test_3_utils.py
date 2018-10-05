import pytest
from kipoi_dataloaders.utils import resize_pybedtools_interval, get_onehot_shape
from pybedtools import Interval
from kipoi_dataloaders import utils


@pytest.mark.parametrize("how", ['start', 'end', 'center'])
@pytest.mark.parametrize("ilen", [3, 4])
def test_resize_pybedtools_interval(how, ilen):
    dummy_start, dummy_end = 10, 20
    dummy_centre = int((dummy_start + dummy_end) / 2)
    dummy_inter = Interval("chr1", dummy_start, dummy_end)
    ret_inter = resize_pybedtools_interval(dummy_inter, how, ilen)
    assert ret_inter.length == ilen
    if how == "start":
        assert ret_inter.start == dummy_start
    elif how == "end":
        assert ret_inter.end == dummy_end
    elif how == "centre":
        assert int((ret_inter.start + ret_inter.end) / 2) == dummy_centre


def test_get_alphabets():
    avail_alphabets = list(utils.alphabets.keys())
    assert all([getattr(utils, k) == utils.get_alphabet(k) for k in avail_alphabets])
    assert all([getattr(utils, k) == utils.get_alphabet(getattr(utils, k)) for k in avail_alphabets])


@pytest.mark.parametrize("alphabet_axis", list(range(1, 4)))
@pytest.mark.parametrize("dummy_axis", [None] + list(range(1, 4)))
def test_get_onehot_shape(alphabet_axis, dummy_axis):
    seq_len, alphabet_len = 1000, 4
    dummy_axis_int = dummy_axis
    if dummy_axis is None:
        dummy_axis_int = -1

    if (alphabet_axis == dummy_axis_int) or (alphabet_axis == 0) or (dummy_axis_int == 0) or \
            (alphabet_axis >= 4) or (dummy_axis_int >= 4) or ((alphabet_axis >= 3) and (dummy_axis is None)):
        return None

    num_axes = 3
    if dummy_axis is not None:
        num_axes += 1

    shape = get_onehot_shape(alphabet_axis, dummy_axis, seq_len, utils.DNA)
    for i in range(len(shape)):
        i2 = i + 1
        if i2 == dummy_axis:
            assert shape[i] == 1
        elif i2 == alphabet_axis:
            assert shape[i] == alphabet_len
        else:
            assert shape[i] == seq_len
