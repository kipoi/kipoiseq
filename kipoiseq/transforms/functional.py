from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

from kipoiseq.utils import DNA
from copy import deepcopy
import numpy as np
from six import string_types


try:
    # use the fast genomelake's one-hot-encode if it's installed
    from genomelake.util import one_hot_encode_sequence
except ImportError:
    one_hot_encode_sequence = None


# sequence -> array

def _get_alphabet_dict(alphabet):
    return {l: i for i, l in enumerate(alphabet)}


def _get_index_dict(alphabet):
    return {i: l for i, l in enumerate(alphabet)}


def one_hot2token(arr):
    return arr.argmax(axis=2)


def one_hot2string(arr, alphabet=DNA):
    """Convert a one-hot encoded array back to string
    """
    tokens = one_hot2token(arr)
    indexToLetter = _get_index_dict(alphabet)

    return [''.join([indexToLetter[x] for x in row]) for row in tokens]


def tokenize(seq, alphabet=DNA, neutral_alphabet=["N"]):
    """Convert sequence to integers

    # Arguments
       seq: Sequence to encode
       alphabet: Alphabet to use
       neutral_alphabet: Neutral alphabet -> assign those values to -1

    # Returns
       List of length `len(seq)` with integers from `-1` to `len(alphabet) - 1`
    """
    # Req: all alphabets have the same length
    if isinstance(neutral_alphabet, str):
        neutral_alphabet = [neutral_alphabet]

    nchar = len(alphabet[0])
    for l in alphabet + neutral_alphabet:
        assert len(l) == nchar
    assert len(seq) % nchar == 0  # since we are using striding

    alphabet_dict = _get_alphabet_dict(alphabet)
    for l in neutral_alphabet:
        alphabet_dict[l] = -1
    # current performance bottleneck
    return np.array([alphabet_dict[seq[(i * nchar):((i + 1) * nchar)]] for i in range(len(seq) // nchar)])


def token2one_hot(tokens, alphabet_size=4, neutral_value=.25, dtype=None):
    """
    Note: everything out of the alphabet is transformed into `np.zeros(alphabet_size)`
    """
    arr = np.zeros((len(tokens), alphabet_size), dtype=dtype)

    tokens_range = np.arange(len(tokens), dtype=int)
    arr[tokens_range[tokens >= 0], tokens[tokens >= 0]] = 1
    if neutral_value != 0:
        arr[tokens_range[tokens < 0], :] = neutral_value
    return arr


def one_hot(seq, alphabet=DNA, neutral_alphabet=['N'], neutral_value=.25, dtype=None):
    if not isinstance(seq, str):
        raise ValueError("seq needs to be a string")
    return token2one_hot(tokenize(seq, alphabet, neutral_alphabet), len(alphabet), neutral_value, dtype=dtype)


def one_hot_dna(seq, dtype=None):
    """One-hot encode DNA sequence
    """
    if not isinstance(seq, str):
        raise ValueError("seq needs to be a string")

    if one_hot_encode_sequence is not None:
        # genomelake's one_hot_encode_sequence could be imported
        out = np.zeros((len(seq), 4), dtype=np.float32)
        one_hot_encode_sequence(seq, out)
        return out.astype(dtype)
    else:
        return one_hot(seq, alphabet=DNA, neutral_alphabet=['N'], neutral_value=.25, dtype=dtype)

# sequence trimming


def pad(seq, length, value="N", anchor="center"):
    seq_len = len(seq)
    assert length >= seq_len
    if anchor is "end":
        n_left = length - seq_len
        n_right = 0
    elif anchor is "start":
        n_right = length - seq_len
        n_left = 0
    elif anchor is "center":
        n_left = (length - seq_len) // 2 + (length - seq_len) % 2
        n_right = (length - seq_len) // 2
    else:
        raise ValueError("anchor can be of: end, start or center")

    # normalize for the length
    n_left = n_left // len(value)
    n_right = n_right // len(value)

    return value * n_left + seq + value * n_right


def trim(seq, length, anchor="center"):
    seq_len = len(seq)

    assert length <= seq_len
    if anchor is "end":
        return seq[-length:]
    elif anchor is "start":
        return seq[0:length]
    elif anchor is "center":
        dl = seq_len - length
        n_left = dl // 2 + dl % 2
        n_right = seq_len - dl // 2
        return seq[n_left:n_right]
    else:
        raise ValueError("anchor can be of: end, start or center")


def fixed_len(seq, length, anchor="center", value="N"):
    """Pad and/or trim a list of sequences to have common length. Procedure:

        1. Pad the sequence with N's or any other string or list element (`value`)
        2. Subset the sequence

    # Note
        See also: https://keras.io/preprocessing/sequence/
        Aplicable also for lists of characters

    # Arguments
        sequence_vec: list of chars or lists
            List of sequences that can have various lengths
        value: Neutral element to pad the sequence with. Can be `str` or `list`.
        length: int or None; Final lenght of sequences.
             If None, length is set to the longest sequence length.
        anchor: character; 'start', 'end' or 'center'
            To which end to anchor the sequences when triming/padding. See examples bellow.

    # Returns
        List of sequences of the same class as sequence_vec

    # Example

        ```python
            >>> sequence = 'CTTACTCAGA'
            >>> pad_sequence(sequence, 10, anchor="start", value="N")
            'CTTACTCAGA'
            >>> pad_sequence(sequence, 10, anchor="end", value="N")
            'CTTACTCAGA'
            >>> pad_sequences(sequence, 4, anchor="center", value="N")
            'ACTC'

            >>> sequence = 'TCTTTA'
            >>> pad_sequence(sequence, 10, anchor="start", value="N")
            'TCTTTANNNN'
            >>> pad_sequence(sequence, 10, anchor="end", value="N")
            'NNNNTCTTTA'
            >>> pad_sequences(sequence, 4, anchor="center", value="N")
            'CTTT'
        ```
    """
    # neutral element type checking
    assert isinstance(value, list) or isinstance(value, str)
    assert isinstance(value, type(seq))
    assert isinstance(length, int)

    # pad and subset
    if len(seq) < length:
        return pad(seq, length, value=value, anchor=anchor)
    elif len(seq) > length:
        return trim(seq, length, anchor=anchor)
    else:
        return seq


def resize_interval(interval, width, anchor='center'):
    """Resize the Interval. Returns new Interval instance with correct length.

    Arguments:
        interval: pybedtools.Interval object or an object containing `start` and `end` attributes
        width: desired width of the output interval
        anchor (str): which part of the sequence should be anchored. Choices: 'start', 'center', or 'end'
    """
    interval = deepcopy(interval)

    if anchor == "start":
        interval.end = interval.start + width
    elif anchor == "end":
        interval.start = interval.end - width
    elif anchor == "center":
        center = int((interval.start + interval.end) / 2)
        half_len = int(width / 2)
        interval.start = center - half_len
        interval.end = center + half_len + width % 2
    else:
        raise Exception("Interval resizing anchor point can only be 'start', 'end' or 'center'")

    return interval
