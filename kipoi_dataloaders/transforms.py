import numpy as np
from six import string_types


def get_string_transforms(trafo):
    if trafo is not None and isinstance(trafo, string_types):
        if trafo in TRANSFORMS:
            trafo = TRANSFORMS[trafo]
    return trafo


def onehot(seq, out, alphabet):
    for i, char in enumerate(seq):
        if char in alphabet:
            out[i, alphabet.index(char)] = 1


def onehot_transform(seqs, alphabet):
    """
    Transforms list of string sequences to one-hot encoded arrays
    
    Arguments:
        seqs (list): List of string sequences. All sequences have to be same length
        alphabet (list): characters present in the alphabet,
            e.g.: ["A", "C", "G", "T"] for DNA.
    """
    output_shape = (len(seqs), len(seqs[0]), len(alphabet))
    out = np.zeros(output_shape, dtype=np.float32)

    for i, seq in enumerate(seqs):
        onehot(seq, out[i, ...], alphabet)

    return out


class TransformShape(object):
    def __init__(self, alphabet_axis, dummy_axis=None):
        """
        Reshape input array as defined by parameters
        Arguments:
            alphabet_axis: index of the alphabet axis. Cannot be 0.
            dummy_axis: index at which a dummy dimension should be added. If None no dummy axis is generated.

        returns:
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
        self.trafos = []

        # global_axis_offset is necessary as we want to count from the right to left since the transformation can be
        # performed on batches or on single samples.
        global_axis_offset = - num_axes
        existing_alph_axis = 2
        if num_axes >= 3:
            self.trafos.append(lambda x: x.__getitem__(tuple([Ellipsis] + [None] * (self.num_axes - 3))))
            if dummy_axis is not None:
                self.trafos.append(lambda x: np.swapaxes(x, -1, dummy_axis + global_axis_offset))
                if dummy_axis == existing_alph_axis:
                    existing_alph_axis = num_axes - 1
        if alphabet_axis != existing_alph_axis:
            self.trafos.append(lambda x: np.swapaxes(x, existing_alph_axis + global_axis_offset,
                                                     self.alphabet_axis + global_axis_offset))

    def _apply_fns(self, input):
        out = input
        for fn in self.trafos:
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


TRANSFORMS = {"onehot_trafo": onehot_transform}
