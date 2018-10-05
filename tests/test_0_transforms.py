import pytest
import numpy as np
import copy
from kipoi_dataloaders.transforms import onehot, onehot_transform, TransformShape
from kipoi_dataloaders.utils import DNA


def test_onehot():
    seq = "ACGTN"
    dat = np.zeros((len(seq), 4))
    onehot(seq, dat, DNA)
    for i, letter in enumerate(seq):
        if letter in DNA:
            assert dat[i, :].sum() == 1
            assert dat[i, DNA.index(letter)].sum() == 1
        else:
            assert dat[i, :].sum() == 0


def test_onehot_trafo():
    seqs = ["ACGTN", "ACGTN"]
    dat = onehot_transform(seqs, DNA)
    assert len(dat.shape) == 3
    assert dat.shape[0] == 2
    assert dat.shape[1] == len(seqs[0])
    assert dat.shape[2] == len(DNA)


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
