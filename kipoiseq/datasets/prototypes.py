from collections import OrderedDict
from kipoi.specs import DataLoaderArgument, DataLoaderSchema, ArraySchema, MetadataStruct, MetadataType
from kipoi.data import Dataset

default_example_files = {
    'intervals_file': 'example_files/intervals_files_ENCSR000EMT_chr21_10000.tsv',
    'fasta_file': 'example_files/chr21.fa',
}

args_prototype = OrderedDict()
args_prototype['intervals_file'] = DataLoaderArgument(
    doc="tsv file containing dna interval indices (chr, start, end) and optional additional labels",
    name='intervals_file', type='str', optional=False, example=default_example_files['intervals_file'])
args_prototype['fasta_file'] = DataLoaderArgument(doc="fasta file for dna intervals",
                                                  name='fasta_file', type='str', optional=False,
                                                  example=default_example_files['fasta_file'])
args_prototype['num_chr_fasta'] = DataLoaderArgument(
    doc="True, the tsv-loader will make sure that the chromosomes don't start with chr.",
    name='num_chr_fasta', type='bool', optional=True)
args_prototype['label_dtype'] = DataLoaderArgument(
    doc="None, datatype of the task labels taken from the intervals_file. Allowed: "
        "'string', 'int', 'float', 'bool'",
    name='label_dtype', type='string', optional=True)
args_prototype['required_seq_len'] = DataLoaderArgument(
    doc="None, required sequence length.",
    name='required_seq_len', type='int', optional=True)
args_prototype['max_seq_len'] = DataLoaderArgument(
    doc="None, maximum allowed sequence length.",
    name='max_seq_len', type='int', optional=True)
args_prototype['use_strand'] = DataLoaderArgument(
    doc="False, reverse-complement fasta sequence if bed file defines negative strand.",
    name='use_strand', type='bool', optional=True)
args_prototype['force_upper'] = DataLoaderArgument(
    doc="True, Force uppercase output of sequences.",
    name='force_upper', type='bool', optional=True)
args_prototype['auto_resize'] = DataLoaderArgument(
    doc="None, Automatically resize the given bed input to the required_seq_len. Allowed arguments: "
        "'start': keeps the start coordinate, 'end', 'center' accordingly.",
    name='auto_resize', type='string', optional=True)
args_prototype['alphabet_axis'] = DataLoaderArgument(
    doc="2, axis along which the alphabet runs (e.g. A,C,G,T for DNA).",
    name='alphabet_axis', type='int', optional=True)
args_prototype['dummy_axis'] = DataLoaderArgument(
    doc="None, in which dimension a dummy axis should be added.",
    name='dummy_axis', type='int', optional=True)
args_prototype['alphabet'] = DataLoaderArgument(
    doc="'DNA', Alphabet to use for the one-hot encoding. This defines the order of the one-hot encoding."
        " Can either be a list or a string: 'DNA', 'RNA', 'AMINO_ACIDS'.",
    name='alphabet', type='string', optional=True)
