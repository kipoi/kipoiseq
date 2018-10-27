import numpy as np
from six import string_types


# alphabets:
DNA = ["A", "C", "G", "T"]
RNA = ["A", "C", "G", "U"]
AMINO_ACIDS = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H",
               "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

alphabets = {"DNA": DNA,
             "RNA": RNA,
             "AMINO_ACIDS": AMINO_ACIDS}


def to_scalar(obj):
    """Convert numpy scalar to native scalar
    """
    if isinstance(obj, np.generic):
        return np.asscalar(obj)
    else:
        return obj


def parse_alphabet(alphabet):
    if isinstance(alphabet, str):
        return list(alphabet)
    else:
        return alphabet


def parse_dtype(dtype):
    if isinstance(dtype, string_types):
        try:
            return eval(dtype)
        except Exception as e:
            raise ValueError("Unable to parse dtype: {}. \nException: {}".format(dtype, e))
    else:
        return dtype
