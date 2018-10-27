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


def parse_dtype(dtype):
    dtypes = {'int': int, 'string': str, 'float': float, 'bool': bool}
    if dtype is None:
        return None
    if dtype in list(dtypes.values()):
        return dtype
    if dtype not in dtypes:
        raise Exception("Datatype '{0}' not recognized. Allowed are: {1}".format(dtype, str(list(dtypes.keys()))))
    return dtypes[dtype]


def parse_alphabet(alphabet):
    if isinstance(alphabet, str):
        return list(alphabet)
    else:
        return alphabet


def parse_type(dtype):
    if isinstance(dtype, string_types):
        if dtype in dir(np):
            return getattr(np, dtype)
    else:
        return dtype
