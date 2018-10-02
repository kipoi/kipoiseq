import pytest
import numpy as np
import copy
from kipoi_dataloaders.utils import onehot, onehot_trafo, DNA, ReshapeSeq


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
    dat = onehot_trafo(seqs, DNA)
    assert len(dat.shape) == 3
    assert dat.shape[0] == 2
    assert dat.shape[1] == len(seqs[0])
    assert dat.shape[2] == len(DNA)


def test_reshape_seq():
    n_samples, seq_len, alphabet_len = 20, 1000, 4
    in_array = np.zeros((n_samples, seq_len, alphabet_len))
    num_axes, seq_axis, alphabet_axis = 10, 7, 4
    reshaper = ReshapeSeq(num_axes, seq_axis, alphabet_axis)

    # Test batch reshaping
    reshaped = reshaper.reshape_batch(in_array)
    assert len(reshaped.shape) == num_axes
    for i in range(len(reshaped.shape)):
        if i == 0:
            assert reshaped.shape[i] == n_samples
        elif i == seq_axis:
            assert reshaped.shape[i] == seq_len
        elif i == alphabet_axis:
            assert reshaped.shape[i] == alphabet_len
        else:
            assert reshaped.shape[i] == 1

    # test the single sample works
    reshaped = reshaper.reshape_single_sample(in_array[0,...])

    # Test if fails if input has wrong shape
    with pytest.raises(Exception):
        reshaper.reshape_batch(in_array[0,...])

    # Test if it raises an exception for invalid arguments
    args = [num_axes, seq_axis, alphabet_axis]
    for i in range(len(args)):
        args_new = copy.copy(args)
        args_new[i] = 0
        with pytest.raises(Exception):
            ReshapeSeq(*args_new)
        if i > 0:
            args_new[i] = num_axes
            with pytest.raises(Exception):
                ReshapeSeq(*args_new)
