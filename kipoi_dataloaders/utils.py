# alphabets:
DNA = ["A", "C", "G", "T"]
RNA = ["A", "C", "G", "U"]
AMINO_ACIDS = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H",
               "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

alphabets = {"DNA": DNA, "RNA": RNA, "AMINO_ACIDS": AMINO_ACIDS}

from pybedtools import Interval


def resize_pybedtools_interval(interval, how, new_len):
    """
    Resize the pybedtools interval. Returns new pybedtools Interval instance with correct length.  
    Arguments:
        interval: pybedtools interval object
        how: String indicating which part of the sequence should be maintained: 'start', 'center', or 'end'
        new_len: required length of the interal
    """
    i_props = {k: getattr(interval, k) for k in ["chrom", "start", "end", "name", "score", "strand"]
               if getattr(interval, k) is not None}
    if how == "start":
        i_props["end"] = interval.start + new_len
    elif how == "end":
        i_props["start"] = interval.end - new_len
    elif how == "center":
        center = int((interval.start + interval.end) / 2)
        half_len = int(new_len / 2)
        i_props["start"] = center - half_len
        i_props["end"] = center + half_len + new_len % 2
    else:
        raise Exception("Interval resizing anchor point can only be 'start', 'end' or 'center'")
    return Interval(**i_props)


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
