from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from kipoiseq.transforms import functional as F
from kipoiseq.utils import DNA


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

    def __init__(self, width, anchor='center'):
        self.width = width
        self.anchor = anchor

    def __call__(self, interval):
        return F.resize_interval(interval, self.width, self.anchor)


# Sequences

class OneHot(object):
    """One-hot encode the sequence
    """

    def __init__(self, alphabet=DNA, neutral_alphabet='N', neutral_value=0.25, dtype=None):
        self.alphabet = alphabet
        self.neutral_alphabet = neutral_alphabet
        self.neutral_value = neutral_value
        self.dtype = dtype

    def __call__(self, seq):
        if self.alphabet == DNA and self.neutral_alphabet == 'N' and self.neutral_value == 0.25:
            return F.one_hot_dna(seq, self.dtype)
        else:
            return F.one_hot(seq, self.alphabet, self.neutral_value, self.dtype)

class SplitSplicingSeq(object):
    """ Spilt returned splice sequence (exon with flanking intron) to required format.
        It splits into ['intron5prime', 'acceptor', 'exon', 'donor', 'intron3prime']. 
        'intron5prime' is the intron 5' of the exon, while 'intron3prime' is from the 3'.
        Args:
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
        ''' 
        Args:
            x: a sequence to split
            intron5prime_len: 5' intronic sequence length to take.
            intron3prime_len: 5' intronic sequence length to take.
        '''
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
        acceptor = x[(intron5prime_len - self.acceptor_intron_len):(intron5prime_len + self.acceptor_exon_len)]
        exon = x[(intron5prime_len + self.exon_cut_l):(-intron3prime_len - self.exon_cut_r)]
        donor = x[(-intron3prime_len - self.donor_exon_len):(-intron3prime_len + self.donor_intron_len)]
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