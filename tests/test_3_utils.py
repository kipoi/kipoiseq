import pytest
from kipoiseq.utils import get_onehot_shape
from pybedtools import Interval
from kipoiseq import utils


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
