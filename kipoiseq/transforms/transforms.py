from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from . import functional as F
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


class ResizeInterval(object):

    def __init__(self, width, anchor='center'):
        self.width = width
        self.anchor = anchor

    def __call__(self, interval):
        return F.resize_interval(interval, self.width, self.anchor)


class OneHot(object):
    """One-hot encode the sequence
    """

    def __init__(self, alphabet=DNA, neutral_alphabet='N', neutral_value=0.25):
        self.alphabet = alphabet
        self.neutral_alphabet = neutral_alphabet
        self.neutral_value = neutral_value

    def __call__(self, seq):
        if self.alphabet == DNA and self.neutral_alphabet == 'N' and self.neutral_value == 0.25:
            return F.one_hot_dna(seq)
        else:
            return F.one_hot(seq, self.alphabet, self.neutral_value)


class TransformShape(object):
    def __init__(self, alphabet_axis, dummy_axis=None):
        """Reshape input array as defined by parameters

        # Arguments
            alphabet_axis: index of the alphabet axis. Cannot be 0.
            dummy_axis: index at which a dummy dimension should be added. If None no dummy axis is generated.

        # Returns
            Reshaped sequence array
        """
        num_axes = 3
        dummy_axis_int = dummy_axis
        if dummy_axis is not None:
            num_axes += 1
        else:
            dummy_axis_int = -1

        if alphabet_axis == 0 or dummy_axis_int == 0 or alphabet_axis >= num_axes or dummy_axis_int >= num_axes or \
                alphabet_axis == dummy_axis_int:
            raise Exception("alphabet_axis must be >0 and < {nax},"
                            " dummy_axis must either be None or >0 and < {nax}, and"
                            " alphabet_axis must be !=  dummy_axis".format(nax=num_axes))
        self.num_axes = num_axes
        self.dummy_axis = dummy_axis
        self.alphabet_axis = alphabet_axis
        # prepare the transformations to minimise calculation time on execution (since it may be applied to every
        # sample individually)
        self.transforms = []

        # global_axis_offset is necessary as we want to count from the right to left since the transformation can be
        # performed on batches or on single samples.
        global_axis_offset = - num_axes
        existing_alph_axis = 2
        if num_axes >= 3:
            self.transforms.append(lambda x: x.__getitem__(tuple([Ellipsis] + [None] * (self.num_axes - 3))))
            if dummy_axis is not None:
                self.transforms.append(lambda x: np.swapaxes(x, -1, dummy_axis + global_axis_offset))
                if dummy_axis == existing_alph_axis:
                    existing_alph_axis = num_axes - 1
        if alphabet_axis != existing_alph_axis:
            self.transforms.append(lambda x: np.swapaxes(x, existing_alph_axis + global_axis_offset,
                                                         self.alphabet_axis + global_axis_offset))

    def _apply_fns(self, input):
        out = input
        for fn in self.transforms:
            out = fn(out)
        return out

    def reshape_single_sample(self, input):
        """
        Reshape input array of single sample
        Arguments:
            input: Input numpy array. Requires dimension order: (<seq_axis>, <alphabet_axis>)

        returns:
            Reshaped sequence array
        """
        if len(input.shape) != 2:
            raise Exception("Input array must have dimensions: (<seq_axis>, <alphabet_axis>)")
        return self._apply_fns(input)

    def reshape_batch(self, input):
        """
        Reshape input array of a batch of samples
        Arguments:
            input: Input numpy array. Requires dimension order: (<samples>, <seq_axis>, <alphabet_axis>)

        returns:
            Reshaped sequence array
        """
        if len(input.shape) != 3:
            raise Exception("Input array must have dimensions: (<samples>, <seq_axis>, <alphabet_axis>)")
        return self._apply_fns(input)
