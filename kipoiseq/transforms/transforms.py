from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from kipoiseq.transforms import functional as F
from kipoiseq.utils import DNA, parse_alphabet, parse_dtype


class Compose(object):
    """Composes several transforms together.

    # Arguments

        transforms (list of ``Transform`` objects): list of transforms to compose.

    Example:
        >>> transforms.Compose([
        >>>     transforms.CenterCrop(10),
        >>>     transforms.ToTensor(),
        >>> ])
    """

    def __init__(self, transforms):
        self.transforms = transforms

    def __call__(self, img):
        for t in self.transforms:
            img = t(img)
        return img

    def __repr__(self):
        format_string = self.__class__.__name__ + '('
        for t in self.transforms:
            format_string += '\n'
            format_string += '    {0}'.format(t)
        format_string += '\n)'
        return format_string


# numpy wrappers

class DummyAxis(object):
    """np.expand_dims wrapper - Insert a dummy axis (calls np.expand_dims)
    """

    def __init__(self, axis=None):
        self.axis = axis

    def __call__(self, x):
        if self.axis is not None:
            return np.expand_dims(x, self.axis)
        else:
            return x


class SwapAxes(object):
    """np.swapaxes wrapper

    If any if the axis is None, do nothing.
    """

    def __init__(self, axis1=None, axis2=None):
        self.axis1 = axis1
        self.axis2 = axis2

    def __call__(self, x):
        if self.axis1 is None or self.axis2 is None:
            return x
        else:
            return np.swapaxes(x, self.axis1, self.axis2)


# Intervals


class ResizeInterval(object):
    """Resize the interval
    """

    def __init__(self, width, anchor='center'):
        self.width = width
        self.anchor = anchor

    def __call__(self, interval):
        return F.resize_interval(interval, self.width, self.anchor)


# Sequences

class OneHot(object):
    """One-hot encode the sequence

    # Arguments
      alphabet: alphabet to use for the one-hot encoding. This defines the order of the one-hot encoding.
          Can either be a list or a string: 'ACGT' or ['A, 'C', 'G', 'T']
      neutral_alphabet: which element to use
      neutral_value: value of the neutral element
      dtype: defines the numpy dtype of the returned array.
      alphabet_axis: axis along which the alphabet runs (e.g. A,C,G,T for DNA)
      dummy_axis: defines in which dimension a dummy axis should be added. None if no dummy axis is required.
    """

    def __init__(self, alphabet=DNA, neutral_alphabet='N', neutral_value=0.25, dtype=None):
        self.alphabet = alphabet
        if isinstance(neutral_alphabet, str):
            neutral_alphabet = [neutral_alphabet]
        self.neutral_alphabet = neutral_alphabet
        self.neutral_value = neutral_value
        self.dtype = dtype

    def __call__(self, seq):
        if self.alphabet == DNA and self.neutral_alphabet == ['N'] and self.neutral_value == 0.25:
            return F.one_hot_dna(seq, self.dtype)
        else:
            return F.one_hot(seq,
                             alphabet=self.alphabet,
                             neutral_alphabet=self.neutral_alphabet,
                             neutral_value=self.neutral_value,
                             dtype=self.dtype)


class ReorderedOneHot(object):
    """Flexible one-hot encoding class that can account for
    many different one-hot encoding formats.

    # Arguments
      alphabet: alphabet to use for the one-hot encoding. This defines the order of the one-hot encoding.
          Can either be a list or a string: 'ACGT' or ['A, 'C', 'G', 'T']
      neutral_alphabet: (single string character) neutral element representing
      neutral_value: value of the neutral element
      dtype: defines the numpy dtype of the returned array.
      alphabet_axis: axis along which the alphabet runs (e.g. A,C,G,T for DNA)
      dummy_axis: defines in which dimension a dummy axis should be added. None if no dummy axis is required.

    Examples (`None` = sequence axis):
      - `(None, 4)`: default
      - `(4, None)`: `alphabet_axis=0`
      - `(4, 1, None)`: `alphabet_axis=0, dummy_axis=1`
    """

    def __init__(self,
                 alphabet=DNA,
                 neutral_alphabet='N',
                 neutral_value=0.25,
                 dtype=None,
                 alphabet_axis=1,
                 dummy_axis=None):
        # make sure the alphabet axis and the dummy axis are valid:
        if dummy_axis is not None:
            if alphabet_axis == dummy_axis:
                raise ValueError("dummy_axis can't be the same as dummy_axis")
            if not (dummy_axis >= 0 and dummy_axis <= 2):
                raise ValueError("dummy_axis can be either 0,1 or 2")
        assert alphabet_axis >= 0 and (alphabet_axis < 2 or (
                alphabet_axis <= 2 and dummy_axis is not None))

        self.alphabet_axis = alphabet_axis
        self.dummy_axis = dummy_axis
        self.alphabet = parse_alphabet(alphabet)
        self.dtype = parse_dtype(dtype)
        self.neutral_alphabet = neutral_alphabet
        self.neutral_value = neutral_value

        # set the transform parameters correctly
        if dummy_axis is not None and dummy_axis < 2:
            # dummy axis is added somewhere in the middle, so the alphabet axis is at the end now
            existing_alphabet_axis = 2
        else:
            # alphabet axis stayed the same
            existing_alphabet_axis = 1

        # check if no swapping needed
        if existing_alphabet_axis == self.alphabet_axis:
            self.alphabet_axis = None

        # how to transform the input
        self.transform = Compose([
            OneHot(self.alphabet,
                   neutral_alphabet=self.neutral_alphabet,
                   neutral_value=self.neutral_value,
                   dtype=self.dtype),  # one-hot-encode
            DummyAxis(self.dummy_axis),  # optionally inject the dummy axis
            # put the alphabet axis elsewhere
            SwapAxes(existing_alphabet_axis, self.alphabet_axis),
        ])

    def __call__(self, seq):
        return self.transform(seq)

    def get_output_shape(self, seqlen=None):
        """Compute the output shape
        """
        if self.dummy_axis is not None and self.alphabet_axis == self.dummy_axis:
            raise ValueError("dummy_axis can't be the same as dummy_axis")

        # default
        output_shape = (seqlen, len(self.alphabet))
        alphabet_axis = self.alphabet_axis

        if self.dummy_axis is not None and self.dummy_axis < 2:
            # dummy axis is added somewhere in the middle, so the alphabet axis is at the end now
            existing_alphabet_axis = 2
        else:
            existing_alphabet_axis = 1

        if existing_alphabet_axis == alphabet_axis:
            alphabet_axis = None

        # inject the dummy axis
        if self.dummy_axis is not None:
            output_shape = output_shape[:self.dummy_axis] + \
                           (1,) + output_shape[self.dummy_axis:]

        # swap axes
        if alphabet_axis is not None:
            sh = list(output_shape)
            sh[alphabet_axis], sh[existing_alphabet_axis] = sh[existing_alphabet_axis], sh[alphabet_axis]
            output_shape = tuple(sh)

        return output_shape


# Splicing

class SplitSplicingSeq(object):
    """Split returned splice sequence (exon with flanking intron) to required format.
        It splits into ['intron5prime', 'acceptor', 'exon', 'donor', 'intron3prime'].
        'intron5prime' is the intron 5' of the exon, while 'intron3prime' is from the 3'.

    # Arguments
        exon_cut_l: when extract exon feature, how many base pair to cut out at the begining of an exon
        exon_cut_r: when extract exon feature, how many base pair to cut out at the end of an exon
           (cut out the part that is considered as acceptor site or donor site)
        intron5prime_cut: how many bp to cut out at the end of acceptor intron that consider as acceptor site
        intron3prime_cut: how many bp to cut out at the end of donor intron that consider as donor site
        acceptor_intron_len: what length in acceptor intron to consider for acceptor site model
        acceptor_exon_len: what length in acceptor exon to consider for acceptor site model
        donor_intron_len: what length in donor intron to consider for donor site model
        donor_exon_len: what length in donor exon to consider for donor site model
    """

    def __init__(self,
                 exon_cut_l=0,
                 exon_cut_r=0,
                 intron5prime_cut=6,
                 intron3prime_cut=6,
                 acceptor_intron_len=50,
                 acceptor_exon_len=3,
                 donor_exon_len=5,
                 donor_intron_len=13
                 ):

        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.intron5prime_cut = intron5prime_cut
        self.intron3prime_cut = intron3prime_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len

    def __call__(self,
                 x,
                 intron5prime_len,
                 intron3prime_len
                 ):
        """
        # Arguments
            x: a sequence to split
            intron5prime_len: 5' intronic sequence length to take.
            intron3prime_len: 5' intronic sequence length to take.
        """
        lackl = self.acceptor_intron_len - \
                intron5prime_len  # need to pad N if left seq not enough long
        if lackl >= 0:
            x = "N" * (lackl + 1) + x
            intron5prime_len += lackl + 1
        lackr = self.donor_intron_len - intron3prime_len
        if lackr >= 0:
            x = x + "N" * (lackr + 1)
            intron3prime_len += lackr + 1

        intron5prime = x[:intron5prime_len - self.intron5prime_cut]
        acceptor = x[(intron5prime_len - self.acceptor_intron_len)
                     :(intron5prime_len + self.acceptor_exon_len)]
        exon = x[(intron5prime_len + self.exon_cut_l)
                 :(-intron3prime_len - self.exon_cut_r)]
        donor = x[(-intron3prime_len - self.donor_exon_len)
                  :(-intron3prime_len + self.donor_intron_len)]
        intron3prime = x[-intron3prime_len + self.intron3prime_cut:]

        import warnings
        if donor[self.donor_exon_len:self.donor_exon_len + 2] != "GT":
            warnings.warn("None GT donor", UserWarning)
        if acceptor[self.acceptor_intron_len - 2:self.acceptor_intron_len] != "AG":
            warnings.warn("None AG donor", UserWarning)
        if len(exon) == 0:
            exon = 'N'

        return {
            "intron5prime": intron5prime,
            "acceptor": acceptor,
            "exon": exon,
            "donor": donor,
            "intron3prime": intron3prime
        }
