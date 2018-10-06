from collections import OrderedDict

import pandas as pd
import numpy as np
from kipoi.metadata import GenomicRanges
from kipoi.specs import DataLoaderArgument, ArraySpecialType
from kipoi.plugin import is_installed

import copy
import pybedtools
from pybedtools import BedTool, Interval
from kipoiseq.extractors import FastaStringExtractor
from kipoiseq.transforms import TransformShape
from kipoiseq.transforms.functional import one_hot, resize_interval
from kipoiseq.utils import get_alphabet, get_onehot_shape, to_scalar
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


class BedDataset:
    """Reads a tsv file in the following format:
    ```
    chr  start  stop  task1  task2 ...
    ```

    # Arguments
      tsv_file: tsv file type
      bed_columns: number of columns corresponding to the bed file. All the columns
        after that will be parsed as targets
      num_chr: if specified, 'chr' in the chromosome name will be dropped
      label_dtype: specific data type for labels
      ambiguous_mask: if specified, rows containing only ambiguous_mask values will be skipped
      incl_chromosomes: exclusive list of chromosome names to include in the final dataset.
        if not None, only these will be present in the dataset
      excl_chromosomes: list of chromosome names to omit from the dataset.
    """

    # bed types accorging to
    # https://www.ensembl.org/info/website/upload/bed.html
    bed_types = [str,  # chrom
                 int,  # chromStart
                 int,  # chromEnd
                 str,  # name
                 float,  # score
                 str,  # strand
                 int,  # thickStart
                 int,  # thickEnd
                 str,  # itemRbg
                 int,  # blockCount
                 int,  # blockSizes
                 int]  # blockStarts

    def __init__(self, tsv_file,
                 bed_columns=3,
                 label_dtype=None,
                 num_chr=False,
                 ambiguous_mask=None,
                 incl_chromosomes=None,
                 excl_chromosomes=None):
        self.tsv_file = tsv_file
        self.bed_columns = bed_columns
        self.num_chr = num_chr
        self.label_dtype = label_dtype
        self.ambiguous_mask = ambiguous_mask
        self.incl_chromosomes = incl_chromosomes
        self.excl_chromosomes = excl_chromosomes

        df_peek = pd.read_table(self.tsv_file,
                                header=None,
                                nrows=1,
                                sep='\t')
        self.n_tasks = df_peek.shape[1] - self.bed_columns
        assert self.n_tasks >= 0
        self.df = pd.read_table(self.tsv_file,
                                header=None,
                                dtype={i: d
                                       for i, d in enumerate(self.bed_types[:self.bed_columns] +
                                                             [self.label_dtype] * self.n_tasks)},
                                sep='\t')
        if self.num_chr and self.df.iloc[0][0].startswith("chr"):
            self.df[0] = self.df[0].str.replace("^chr", "")
        if not self.num_chr and not self.df.iloc[0][0].startswith("chr"):
            self.df[0] = "chr" + self.df[0]

        if ambiguous_mask is not None:
            # exclude regions where only ambigous labels are present
            self.df = self.df[~np.all(self.df.iloc[:, self.bed_columns:] == ambiguous_mask, axis=1)]

            # omit data outside chromosomes
        if incl_chromosomes is not None:
            self.df = self.df[self.df[0].isin(incl_chromosomes)]
        if excl_chromosomes is not None:
            self.df = self.df[~self.df[0].isin(excl_chromosomes)]

    def __getitem__(self, idx):
        """Returns (pybedtools.Interval, labels)
        """
        row = self.df.iloc[idx]
        interval = pybedtools.create_interval_from_list([to_scalar(x) for x in row.iloc[:self.bed_columns]])

        if self.n_tasks == 0:
            labels = {}
        else:
            labels = row.iloc[self.bed_columns:].values.astype(self.label_dtype)
        return interval, labels

    def __len__(self):
        return len(self.df)

    def get_targets(self):
        return self.df.iloc[:, self.bed_columns:].values.astype(self.label_dtype)


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

        self.tsv = BedDataset(self.intervals_file,
                              num_chr=self.num_chr_fasta,
                              label_dtype=parse_dtype(label_dtype))
        self.fasta_extractors = None

        # correct the output schema
        self.output_schema_params = copy.copy(self.output_schema_params)

        self.output_schema_params['inputs_shape'] = (1,)
        if self.tsv.n_tasks != 0:
            self.output_schema_params['targets_shape'] = (self.tsv.n_tasks,)

        self.output_schema = get_seq_dataset_output_schema(**self.output_schema_params)

    def __len__(self):
        return len(self.tsv)

    def __getitem__(self, idx):
        if self.fasta_extractors is None:
            self.fasta_extractors = FastaStringExtractor(self.fasta_file, use_strand=self.use_strand,
                                                         force_upper=self.force_upper)

        interval, labels = self.tsv[idx]

        if self.required_seq_len is not None:
            if not interval.stop - interval.start == self.required_seq_len:
                if self.auto_resize is not None:
                    interval = resize_interval(interval, self.auto_resize, self.required_seq_len)
                else:
                    raise Exception("Sequence interval in intervals_file does not match required model sequence "
                                    "length. Update intervals_file or use the 'auto_resize' argument.")

        if self.max_seq_len is not None:
            assert interval.stop - interval.start <= self.max_seq_len

        # Run the fasta extractor and transform if necessary
        seq = self.fasta_extractors.extract(interval)

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

        # TODO - re-use one-hot from above

        # correct the output schema
        self.output_schema_params = copy.copy(self.output_schema_params)

        self.output_schema_params['inputs_shape'] = get_onehot_shape(self.alphabet_axis, self.dummy_axis,
                                                                     self.required_seq_len, self.alphabet)
        if self.tsv.n_tasks != 0:
            self.output_schema_params['targets_shape'] = (self.tsv.n_tasks,)

        self.output_schema = get_seq_dataset_output_schema(**self.output_schema_params)

    def __getitem__(self, idx):
        ret = super(SeqDataset, self).__getitem__(idx)

        # TODO - abandon inheritence and use special transforms instead?
        ret["inputs"] = one_hot(ret["inputs"], alphabet=self.alphabet)
        ret["inputs"] = self.reshaper.reshape_single_sample(ret["inputs"])
        return ret
