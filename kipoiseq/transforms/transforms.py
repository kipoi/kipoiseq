from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from . import functional as F


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


class ResizeInterval(object):

    def __init__(self, width, anchor='center'):
        self.width = width
        self.anchor = anchor

    def __call__(self, interval):
        return F.resize_interval(interval, self.width, self.anchor)


class OneHot(object):
    """One-hot encode the sequence
    """

    def __init__(self, alphabet, neutral_alphabet='N'):
        self.alphabet = alphabet
        self.neutral_alphabet = neutral_alphabet

    def __call__(self, seq):
        if self.alphabet == ["A", "C", "G", "T"] and self.neutral_alphabet == 'N':
            return F.one_hot_dna(seq)
        else:
            return F.one_hot(seq, self.alphabet, self.neutral_element)
