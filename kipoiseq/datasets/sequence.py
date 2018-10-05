from collections import OrderedDict

import pandas as pd
import numpy as np
from kipoi.metadata import GenomicRanges
from kipoi.specs import DataLoaderArgument, ArraySpecialType
from kipoi.plugin import is_installed

import copy
from pybedtools import BedTool, Interval
from kipoiseq.extractors import TsvExtractor, FastaStringExtractor
from kipoiseq.transforms import TransformShape, onehot_transform
from kipoiseq.utils import resize_pybedtools_interval, get_alphabet, get_onehot_shape
from kipoiseq.datasets.prototypes import SeqDatasetPrototype, args_prototype, get_seq_dataset_output_schema


def parse_dtype(dtype):
    dtypes = {'int': int, 'string': str, 'float': float, 'bool': bool}
    if dtype is None:
        return None
    if dtype in list(dtypes.values()):
        return dtype
    if dtype not in dtypes:
        raise Exception("Datatype '{0}' not recognized. Allowed are: {1}".format(dtype, str(list(dtypes.keys()))))
    return dtypes[dtype]


# TODO - this should be the dataset
#   TODO - this should be called: Bed3Dataset
class TsvExtractor:
    """Reads a tsv file in the following format:
    chr  start  stop  task1  task2 ...
    Args:
      tsv_file: tsv file type
      label_dtype: data type of the labels
    """

    def __init__(self, tsv_file, num_chr=False, label_dtype=None):
        self.tsv_file = tsv_file
        self.num_chr = num_chr
        self.label_dtype = label_dtype
        df_peek = pd.read_table(self.tsv_file,
                                header=None,
                                nrows=1,
                                sep='\t')
        self.n_tasks = df_peek.shape[1] - 3
        assert self.n_tasks >= 0
        self.df = pd.read_table(self.tsv_file,
                                header=None,
                                dtype={i: d
                                       for i, d in enumerate([str, int, int] +
                                                             [self.label_dtype] * self.n_tasks)},
                                sep='\t')
        if self.num_chr and self.df.iloc[0][0].startswith("chr"):
            self.df[0] = self.df[0].str.replace("^chr", "")
        if not self.num_chr and not self.df.iloc[0][0].startswith("chr"):
            self.df[0] = "chr" + self.df[0]

    def __getitem__(self, idx):
        """Returns (pybedtools.Interval, labels)
        """
        row = self.df.iloc[idx]
        interval = Interval(row[0], row[1], row[2])

        if self.n_tasks == 0:
            labels = {}
        else:
            labels = row.iloc[3:].values.astype(self.label_dtype)
        return interval, labels

    def __len__(self):
        return len(self.df)


class SeqStringDataset(SeqDatasetPrototype):
    """
    Dataloader for a combination of fasta and tab-delimited input files such as bed files. The dataloader extracts
    regions from the fasta file as defined in the tab-delimited `intervals_file`. Returned sequences are of the type
    np.array([<str>]).
    Arguments:
        intervals_file: bed3+<columns> file containing intervals+labels
        fasta_file: file path; Genome sequence
        num_chr_fasta: if True, the tsv-loader will make sure that the chromosomes
          don't start with chr
        label_dtype: label data type
        required_seq_len: required sequence length
        max_seq_len: maximum allowed sequence length
        use_strand: reverse-complement fasta sequence if bed file defines negative strand
        force_upper: Force uppercase output of sequences
        auto_resize: Automatically resize the given bed input to the required_seq_len. Allowed arguments: 
            'start': keeps the start coordinate, 'end', 'center' accordingly.
    """
    defined_as = 'kipoiseq.SeqStringDataset'
    # todo: use is_installed
    if is_installed("kipoi_veff"):
        from kipoi_veff.specs import VarEffectDataLoaderArgs
        postprocessing = OrderedDict([('variant_effects', VarEffectDataLoaderArgs(bed_input=["intervals_file"]))])
    args = OrderedDict([(k, args_prototype[k]) for k in ["intervals_file", "fasta_file", "num_chr_fasta", "label_dtype",
                                                         "required_seq_len", "max_seq_len", "use_strand", "force_upper",
                                                         "auto_resize"]])

    output_schema_params = {"inputs_special_type": ArraySpecialType.DNAStringSeq}
    output_schema = get_seq_dataset_output_schema(**output_schema_params)

    def __init__(self, intervals_file, fasta_file, num_chr_fasta=False, label_dtype=None,
                 required_seq_len=None, max_seq_len=None, use_strand=False, force_upper=True,
                 auto_resize=None):

        self.num_chr_fasta = num_chr_fasta
        self.intervals_file = intervals_file
        self.fasta_file = fasta_file
        self.required_seq_len = required_seq_len
        self.use_strand = use_strand
        self.force_upper = force_upper
        self.max_seq_len = max_seq_len
        self.auto_resize = auto_resize

        self.tsv = TsvExtractor(self.intervals_file,
                                num_chr=self.num_chr_fasta,
                                label_dtype=parse_dtype(label_dtype))
        self.fasta_reader = None

        # correct the output schema
        self.output_schema_params = copy.copy(self.output_schema_params)

        self.output_schema_params['inputs_shape'] = (1,)
        if self.tsv.n_tasks != 0:
            self.output_schema_params['targets_shape'] = (self.tsv.n_tasks,)

        self.output_schema = get_seq_dataset_output_schema(**self.output_schema_params)

    def __len__(self):
        return len(self.tsv)

    def __getitem__(self, idx):
        if self.fasta_reader is None:
            self.fasta_reader = FastaStringExtractor(self.fasta_file, use_strand=self.use_strand,
                                                     force_upper=self.force_upper)

        interval, labels = self.tsv[idx]

        if self.required_seq_len is not None:
            if not interval.stop - interval.start == self.required_seq_len:
                if self.auto_resize is not None:
                    interval = resize_pybedtools_interval(interval, self.auto_resize, self.required_seq_len)
                else:
                    raise Exception("Sequence interval in intervals_file does not match required model sequence "
                                    "length. Update intervals_file or use the 'auto_resize' argument.")

        if self.max_seq_len is not None:
            assert interval.stop - interval.start <= self.max_seq_len

        # Run the fasta extractor and transform if necessary
        seq = self.fasta_reader([interval])

        return {
            "inputs": np.array(seq),
            "targets": labels,
            "metadata": {
                "ranges": GenomicRanges(interval.chrom, interval.start, interval.stop, str(idx))
            }
        }


# TODO - check lzamparo's dataloader:
# - https://github.com/kipoi/kipoiseq/issues/1#issuecomment-427412487
# - https://raw.githubusercontent.com/lzamparo/bindspace_revisions/master/deepbind/src/dataloader.py

class SeqDataset(SeqStringDataset):
    """
    Dataloader for a combination of fasta and tab-delimited input files such as bed files. The dataloader extracts
    regions from the fasta file as defined in the tab-delimited `intervals_file` and converts them into one-hot encoded
    format. Returned sequences are of the type np.array with the shape inferred from the arguments: `alphabet_axis`
    and `dummy_axis`.

    Arguments:
        intervals_file: bed3+<columns> file containing intervals+labels
        fasta_file: file path; Genome sequence
        num_chr_fasta: if True, the tsv-loader will make sure that the chromosomes
          don't start with chr
        label_dtype: label data type
        required_seq_len: required sequence length
        use_strand: reverse-complement fasta sequence if bed file defines negative strand
        auto_resize: Automatically resize the given bed input to the required_seq_len. Allowed arguments:
            'start': keeps the start coordinate, 'end', 'center' accordingly.
        alphabet_axis: axis along which the alphabet runs (e.g. A,C,G,T for DNA)
        dummy_axis: defines in which dimension a dummy axis should be added. None if no dummy axis is required.
        alphabet: alphabet to use for the one-hot encoding. This defines the order of the one-hot encoding.
            Can either be a list or a string: 'DNA', 'RNA', 'AMINO_ACIDS'.
    """
    defined_as = 'kipoiseq.SeqDataset'
    args = OrderedDict([(k, args_prototype[k]) for k in ["intervals_file", "fasta_file", "num_chr_fasta", "label_dtype",
                                                         "required_seq_len", "max_seq_len", "use_strand", "auto_resize",
                                                         "alphabet_axis", "dummy_axis", "alphabet"]])

    output_schema_params = {"inputs_special_type": ArraySpecialType.DNASeq}
    output_schema = get_seq_dataset_output_schema(**output_schema_params)

    def __init__(self, intervals_file, fasta_file, num_chr_fasta=False, label_dtype=None, required_seq_len=None,
                 max_seq_len=None, use_strand=False, auto_resize=None, alphabet_axis=2, dummy_axis=None,
                 alphabet="DNA"):
        super(SeqDataset, self).__init__(intervals_file, fasta_file, num_chr_fasta=num_chr_fasta,
                                         label_dtype=label_dtype, required_seq_len=required_seq_len,
                                         max_seq_len=max_seq_len, auto_resize=auto_resize,
                                         use_strand=use_strand, force_upper=True)

        self.alphabet_axis = alphabet_axis
        self.dummy_axis = dummy_axis
        self.reshaper = TransformShape(alphabet_axis, dummy_axis)
        self.alphabet = get_alphabet(alphabet)

        # correct the output schema
        self.output_schema_params = copy.copy(self.output_schema_params)

        self.output_schema_params['inputs_shape'] = get_onehot_shape(self.alphabet_axis, self.dummy_axis,
                                                                     self.required_seq_len, self.alphabet)
        if self.tsv.n_tasks != 0:
            self.output_schema_params['targets_shape'] = (self.tsv.n_tasks,)

        self.output_schema = get_seq_dataset_output_schema(**self.output_schema_params)

    def __getitem__(self, idx):
        ret = super(SeqDataset, self).__getitem__(idx)
        ret["inputs"] = onehot_transform([ret["inputs"][0]], alphabet=self.alphabet)[0]
        ret["inputs"] = self.reshaper.reshape_single_sample(ret["inputs"])
        return ret
