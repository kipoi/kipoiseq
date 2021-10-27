from itertools import islice
from typing import Iterable, TypeVar

import numpy as np
from six import string_types

# alphabets:
from kipoiseq import Variant

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
        return obj.item()
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


T = TypeVar('T')


def batch_iter(items: Iterable[T], batch_size: int) -> Iterable[Iterable[T]]:
    # ensure this is an iterator
    item_iter = iter(items)
    while True:
        # create next `batch_size` number of items;
        batch = list(islice(item_iter, batch_size))
        if len(batch) == 0:
            break

        yield batch
