import numpy as np
from six import string_types

def get_string_trafo(trafo):
    if trafo is not None and isinstance(trafo, string_types):
        if trafo in TRAFOS:
            trafo = TRAFOS[trafo]
    return trafo

def onehot(seq, out, alphabet):
    for i, char in enumerate(seq):
        if char in alphabet:
            out[i, alphabet.index(char)] = 1

def onehot_trafo(seqs, alphabet):
    """
    Transforms list of string sequences to one-hot encoded arrays
    
    Arguments:
        seqs (list): List of string sequences. All sequences have to be same length
        alphabet (list): characters present in the alphabet,
            e.g.: ["A", "C", "G", "T"] for DNA.
    """
    output_shape = (len(seqs), len(seqs[0]), len(alphabet))
    out = np.zeros(output_shape)

    for i, seq in enumerate(seqs):
        onehot(seq, out[i,...], alphabet)

    return out


class ReshapeSeq(object):
    def __init__(self, num_axes, seq_axis, alphabet_axis):
        """
        Reshape input array as defined by parameters
        Arguments:
            input: Input numpy array. Requires dimension order: (<samples>, <seq_axis>, <alphabet_axis>)
            num_axes: total dimensionality of output array including the sample axis. Must be >=3.
            seq_axis: index of the sequence axis. Cannot be 0.
            alphabet_axis: index of the alphabet axis. Cannot be 0.
        
        returns:
            Reshaped sequence array 
        """
        if num_axes <3 or seq_axis == 0 or alphabet_axis == 0 or seq_axis >= num_axes or alphabet_axis >= num_axes:
            raise Exception("'num_axes' must be >= 3, seq_axis must be >0 and < {nax},"
                            " alphabet_axis must be >0 and < {nax}".format(nax = num_axes))
        self.num_axes = num_axes
        self.seq_axis = seq_axis
        self.alphabet_axis = alphabet_axis
        # prepare the transformations to minimise calculation time on execution (since it may be applied to every
        # sample individually)
        self.trafos = []
        if num_axes >= 3:
            self.trafos.append(lambda x: x.__getitem__(tuple([Ellipsis] + [None]*(self.num_axes-3))))
        if seq_axis != 1:
            self.trafos.append(lambda x: np.swapaxes(x, 1, self.seq_axis))
        if alphabet_axis != 2:
            self.trafos.append(lambda x: np.swapaxes(x, 2, self.alphabet_axis))


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




TRAFOS = {"onehot_trafo": onehot_trafo}