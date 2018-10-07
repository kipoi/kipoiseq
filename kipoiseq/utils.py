import numpy as np
from pybedtools import Interval


# alphabets:
DNA = ["A", "C", "G", "T"]
RNA = ["A", "C", "G", "U"]
AMINO_ACIDS = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H",
               "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

alphabets = {"DNA": DNA,
             "RNA": RNA,
             "AMINO_ACIDS": AMINO_ACIDS}


def get_alphabet(alphabet):
    if isinstance(alphabet, list):
        return alphabet
    else:
        return alphabets[alphabet]


def get_onehot_shape(alphabet_axis, dummy_axis, seq_len, alphabet):
    # these values are given with respect to batch, not the single sample:
    dim = 3 + int(dummy_axis is not None)
    non_seq_axis = [alphabet_axis] + [da for da in [dummy_axis] if da is not None]
    seq_axis = [i for i in range(1, dim) if i not in non_seq_axis][0]

    # now for the single sample assign the shape
    shape = [1] * (dim - 1)
    shape[alphabet_axis - 1] = len(alphabet)
    shape[seq_axis - 1] = seq_len
    return tuple(shape)


def to_scalar(obj):
    """Convert numpy scalar to native scalar
    """
    if isinstance(obj, np.generic):
        return np.asscalar(obj)
    else:
        return obj
