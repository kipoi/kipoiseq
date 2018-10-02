from kipoi.specs import DataLoaderArgument
from .utils.readers import TsvReader, FastaReader
from kipoi.metadata import GenomicRanges
from kipoi.data import Dataset
from .utils import get_string_trafo, ReshapeSeq, DNA
import numpy as np
from collections import OrderedDict


# TODO - set the example files correctly
default_example_files = {
'intervals_file' : 'example_files/intervals_files_ENCSR000EMT_chr21_10000.tsv',
'fasta_file' : 'example_files/chr21.fa',
'num_chr_fasta' : True,
}

def parse_dtype(dtype):
    dtypes = {'int':int, 'string':str, 'float':float, 'bool':bool}
    if dtype is None:
        return None
    if dtype in list(dtypes.values()):
        return dtype
    if dtype not in dtypes:
        raise Exception("Datatype '{0}' not recognized. Allowed are: {1}".format(dtype, str(list(dtypes.keys()))))
    return dtypes[dtype]


class FastaBasedDataset(Dataset):
    """
    Args:
        intervals_file: bed3+<columns> file containing intervals+labels
        fasta_file: file path; Genome sequence
        num_chr_fasta: if True, the tsv-loader will make sure that the chromosomes
          don't start with chr
        label_dtype: label data type
        seq_len: required sequence length
        use_strand: reverse-complement fasta sequence if bed file defines negative strand
        force_upper: Force uppercase output of sequences
    """
    output_schema = None
    type = 'Dataset'
    defined_as = 'kipoi_dataloaders.FastaBasedDataset'
    info = None
    args = OrderedDict()
    args['intervals_file'] = DataLoaderArgument(
        doc="tsv file containing dna interval indices (chr, start, end) and optional additional labels",
        name='intervals_file', type='str', optional=False)
    args['fasta_file'] = DataLoaderArgument(doc="fasta file for dna intervals",
        name='fasta_file', type='str', optional=False)
    args['num_chr_fasta'] = DataLoaderArgument(
        doc="True, the tsv-loader will make sure that the chromosomes don't start with chr.",
        name='num_chr_fasta', type='bool', optional=True)
    args['label_dtype'] = DataLoaderArgument(
        doc="None, datatype of the task labels taken from the intervals_file. Allowed: "
            "'string', 'int', 'float', 'bool'",
        name='label_dtype', type='string', optional=True)
    args['seq_len'] = DataLoaderArgument(
        doc="None, required sequence length.",
        name='seq_len', type='int', optional=True)
    args['use_strand'] = DataLoaderArgument(
        doc="False, reverse-complement fasta sequence if bed file defines negative strand.",
        name='use_strand', type='bool', optional=True)
    args['force_upper'] = DataLoaderArgument(
        doc="True, Force uppercase output of sequences.",
        name='force_upper', type='bool', optional=True)
    dependencies = None
    postprocessing = OrderedDict()
    _yaml_path = None
    source = None
    source_dir = None
    # TODO - enable postprocessing
    # TODO - the examples have to be set in the `DataLoaderArgument` objects


    def __init__(self, intervals_file, fasta_file, num_chr_fasta=False, label_dtype = None,
                 seq_len=None, use_strand=False, force_upper=True):

        self.num_chr_fasta = num_chr_fasta
        self.intervals_file = intervals_file
        self.fasta_file = fasta_file
        self.seq_len = seq_len
        self.use_strand = use_strand
        self.force_upper = force_upper

        self.tsv = TsvReader(self.intervals_file,
                             num_chr=self.num_chr_fasta,
                             label_dtype=parse_dtype(label_dtype))
        self.fasta_reader = None


        # TODO - optionally
        # - once you know the schema, modify it (optionally)


    def __len__(self):
        return len(self.tsv)


    def __getitem__(self, idx):
        if self.fasta_reader is None:
            self.fasta_reader = FastaReader(self.fasta_file, use_strand=self.use_strand, force_upper=self.force_upper)

        interval, labels = self.tsv[idx]

        if self.seq_len is not None:
            # Intervals need to be 1000bp wide
            assert interval.stop - interval.start == self.seq_len

        # Run the fasta extractor and transform if necessary
        seq = self.fasta_reader([interval])[0]
        # TODO - check: is it allowed to just return a single string as dataloader?

        return {
            "inputs": seq,
            "targets": labels,
            "metadata": {
                "ranges": GenomicRanges(interval.chrom, interval.start, interval.stop, str(idx))
            }
        }


class SeqDataset(FastaBasedDataset):
    """
    Args:
        intervals_file: bed3+<columns> file containing intervals+labels
        fasta_file: file path; Genome sequence
        num_chr_fasta: if True, the tsv-loader will make sure that the chromosomes
          don't start with chr
        label_dtype: label data type
        seq_len: required sequence length
        use_strand: reverse-complement fasta sequence if bed file defines negative strand
        num_axes: dimensionality of returned array (this value includes the sample axis)
        seq_axis: axis along which the sequence is defined
        alphabet_axis: axis along which the alphabet runs (e.g. A,C,G,T for DNA)
    """
    defined_as = 'kipoi_dataloaders.SeqDataset'
    args = OrderedDict()
    args['intervals_file'] = DataLoaderArgument(
        doc="tsv file containing dna interval indices (chr, start, end) and optional additional labels",
        name='intervals_file', type='str', optional=False)
    args['fasta_file'] = DataLoaderArgument(doc="fasta file for dna intervals",
                                            name='fasta_file', type='str', optional=False)
    args['num_chr_fasta'] = DataLoaderArgument(
        doc="True, the tsv-loader will make sure that the chromosomes don't start with chr.",
        name='num_chr_fasta', type='bool', optional=True)
    args['label_dtype'] = DataLoaderArgument(
        doc="None, datatype of the task labels taken from the intervals_file. Allowed: "
            "'string', 'int', 'float', 'bool'",
        name='label_dtype', type='string', optional=True)
    args['seq_len'] = DataLoaderArgument(
        doc="None, required sequence length.",
        name='seq_len', type='int', optional=True)
    args['use_strand'] = DataLoaderArgument(
        doc="False, reverse-complement fasta sequence if bed file defines negative strand.",
        name='use_strand', type='bool', optional=True)
    args['num_axes'] = DataLoaderArgument(
        doc="3, dimensionality of returned array (this value includes the sample axis)",
        name='num_axes', type='int', optional=True)
    args['seq_axis'] = DataLoaderArgument(
        doc="1, axis along which the sequence is defined.",
        name='seq_axis', type='int', optional=True)
    args['alphabet_axis'] = DataLoaderArgument(
        doc="2, axis along which the alphabet runs (e.g. A,C,G,T for DNA).",
        name='alphabet_axis', type='int', optional=True)

    def __init__(self, intervals_file, fasta_file, num_chr_fasta=False, label_dtype = None, seq_len=None,
                 use_strand=False, num_axes=3, seq_axis=1, alphabet_axis=2):

        super(SeqDataset, self).__init__(intervals_file, fasta_file, num_chr_fasta=num_chr_fasta,
                                         label_dtype = label_dtype, seq_len=seq_len,
                                         use_strand=use_strand, force_upper=True)

        self.num_axes = num_axes
        self.seq_axis = seq_axis
        self.alphabet_axis = alphabet_axis
        self.reshaper = ReshapeSeq(num_axes, seq_axis, alphabet_axis)
        self.out_trafo = get_string_trafo("onehot_trafo")

    def __getitem__(self, idx):
        ret = super(SeqDataset, self).__getitem__(idx)
        ret["inputs"] = self.out_trafo([ret["inputs"]], alphabet=DNA)[0]
        ret["inputs"] = self.reshaper.reshape_single_sample(ret["inputs"])
        return ret