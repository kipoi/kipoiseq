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
