import pytest
from kipoiseq.transforms.functional import tokenize, token2one_hot, one_hot, one_hot_dna, pad, trim, fixed_len
from kipoiseq.utils import DNA
import numpy as np


def test_tokenize():
    assert np.all(tokenize("ACGTTA", DNA, neutral_alphabet="N") == [0, 1, 2, 3, 3, 0])
    assert np.all(tokenize("ACGTGATGA", ["ACG", "TGA"], neutral_alphabet="NNN") == [0, 1, 1])
    assert np.all(tokenize("ACGTGATGA", ["ACG"], neutral_alphabet="TGA") == [0, -1, -1])
    with pytest.raises(Exception):
        tokenize("ACGTGATGA", ["ACG"], neutral_alphabet="NNN")


def test_token2one_hot():
    assert np.array_equal(token2one_hot(np.array([0, 1, -1]), 2), np.array([[1, 0],
                                                                            [0, 1],
                                                                            [0.25, .25]]))


def test_tokenize_one_hot():
    assert one_hot("ACG", DNA, "N").shape == (3, 4)

    et = tokenize("ACG", DNA, "N")
    assert et.shape == (3,)
    assert np.array_equal(et, np.array([0, 1, 2]))

    et = tokenize("TGTN", DNA, "N")
    assert np.array_equal(et, np.array([3, 2, 3, -1]))  # N mapped to -1


def test_one_hot():

    seq = "ACGTTTATNT"
    assert len(seq) == 10

    assert one_hot_dna(seq).shape == (10, 4)
    assert one_hot(seq).shape == (10, 4)

    assert np.all(one_hot_dna(seq) == one_hot(seq))

    assert one_hot(pad(seq, 20)).shape == (20, 4)

    assert one_hot(fixed_len(seq, 20)).shape == (20, 4)
    assert one_hot(fixed_len(seq, 5)).shape == (5, 4)
    assert trim(seq, 5) == 'TTTAT'
    assert trim(seq, 5, 'start') == 'ACGTT'
    assert trim(seq, 5, 'end') == 'TATNT'
    with pytest.raises(Exception):
        assert pad(seq, 5, 'end') == 'TATNT'

    assert np.all(one_hot(seq)[0] == np.array([1, 0, 0, 0]))
    assert np.all(one_hot(seq)[1] == np.array([0, 1, 0, 0]))
    assert np.all(one_hot(seq)[2] == np.array([0, 0, 1, 0]))
    assert np.all(one_hot(seq)[3] == np.array([0, 0, 0, 1]))
    assert np.all(one_hot(seq)[4] == np.array([0, 0, 0, 1]))
    assert np.all(one_hot(seq)[-1] == np.array([0, 0, 0, 1]))
    assert np.all(one_hot(seq)[-2] == np.array([0.25, 0.25, 0.25, 0.25]))


def test_fixed_len():
    seq = "ACGTTTATNT"
    assert len(fixed_len(seq, 20, value="N", anchor="end")) is 20
    assert fixed_len([1, 2, 3, 4], 6, value=[0], anchor="end") == [0, 0, 1, 2, 3, 4]
    assert fixed_len([1, 2, 3, 4], 2, value=[0], anchor="end") == [3, 4]

    # expect error
    assert fixed_len(seq, length=3, value="NNN", anchor="end") == "TNT"


def test_pad_sequences():
    seq = 'CTTACTCAGA'
    assert fixed_len(seq, 4, anchor="center", value="N") == 'ACTC'
    assert fixed_len(seq, 9, anchor="start", value="N") == 'CTTACTCAG'

    assert fixed_len(seq, 10, anchor="start", value="N") == seq
    assert fixed_len(seq, 10, anchor="end", value="N") == 'CTTACTCAGA'
