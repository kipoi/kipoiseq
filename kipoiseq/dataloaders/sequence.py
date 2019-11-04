import pandas as pd
import numpy as np
from copy import deepcopy
from kipoi.metadata import GenomicRanges
from kipoi.data import Dataset, kipoi_dataloader
from kipoi_conda.dependencies import Dependencies
from kipoiseq.transforms import ReorderedOneHot
from kipoi.specs import Author
from kipoi_utils.utils import default_kwargs
from kipoiseq.extractors import FastaStringExtractor
from kipoiseq.transforms.functional import resize_interval
from kipoiseq.utils import to_scalar, parse_dtype


# general dependencies
# bioconda::genomelake', TODO - add genomelake again once it gets released with pyfaidx to bioconda
deps = Dependencies(conda=['bioconda::pybedtools', 'bioconda::pyfaidx', 'numpy', 'pandas'],
                    pip=['kipoiseq'])
package_authors = [Author(name='Ziga Avsec', github='avsecz'),
                   Author(name='Roman Kreuzhuber', github='krrome')]

# Object exported on import *
__all__ = ['SeqIntervalDl', 'StringSeqIntervalDl', 'BedDataset']


class BedDataset(object):
    """Reads a tsv file in the following format:
    ```
    chr  start  stop  task1  task2 ...
    ```

    # Arguments
      tsv_file: tsv file type
      bed_columns: number of columns corresponding to the bed file. All the columns
        after that will be parsed as targets
      num_chr: if specified, 'chr' in the chromosome name will be dropped
      label_dtype: specific data type for labels, Example: `float` or `np.float32`
      ambiguous_mask: if specified, rows containing only ambiguous_mask values will be skipped
      incl_chromosomes: exclusive list of chromosome names to include in the final dataset.
        if not None, only these will be present in the dataset
      excl_chromosomes: list of chromosome names to omit from the dataset.
      ignore_targets: if True, target variables are ignored
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
                 label_dtype=None,
                 bed_columns=3,
                 num_chr=False,
                 ambiguous_mask=None,
                 incl_chromosomes=None,
                 excl_chromosomes=None,
                 ignore_targets=False):
        # TODO - `chrom` column: use pd.Categorical for memory efficiency
        self.tsv_file = tsv_file
        self.bed_columns = bed_columns
        self.num_chr = num_chr
        self.label_dtype = label_dtype
        self.ambiguous_mask = ambiguous_mask
        self.incl_chromosomes = incl_chromosomes
        self.excl_chromosomes = excl_chromosomes
        self.ignore_targets = ignore_targets

        df_peek = pd.read_table(self.tsv_file,
                                header=None,
                                nrows=1,
                                sep='\t')
        found_columns = df_peek.shape[1]
        self.n_tasks = found_columns - self.bed_columns
        if self.n_tasks < 0:
            raise ValueError("BedDataset requires at least {} bed columns. Found only {} columns".
                             format(self.bed_columns, found_columns))

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
            self.df = self.df[~np.all(
                self.df.iloc[:, self.bed_columns:] == ambiguous_mask, axis=1)]

            # omit data outside chromosomes
        if incl_chromosomes is not None:
            self.df = self.df[self.df[0].isin(incl_chromosomes)]
        if excl_chromosomes is not None:
            self.df = self.df[~self.df[0].isin(excl_chromosomes)]

    def __getitem__(self, idx):
        """Returns (pybedtools.Interval, labels)
        """
        row = self.df.iloc[idx]

        # TODO: use kipoiseq.dataclasses.interval instead of pybedtools
        import pybedtools
        interval = pybedtools.create_interval_from_list(
            [to_scalar(x) for x in row.iloc[:self.bed_columns]])

        if self.ignore_targets or self.n_tasks == 0:
            labels = {}
        else:
            labels = row.iloc[self.bed_columns:].values.astype(
                self.label_dtype)
        return interval, labels

    def __len__(self):
        return len(self.df)

    def get_targets(self):
        return self.df.iloc[:, self.bed_columns:].values.astype(self.label_dtype)


@kipoi_dataloader(override={"dependencies": deps, 'info.authors': package_authors})
class StringSeqIntervalDl(Dataset):
    """
    info:
        doc: >
           Dataloader for a combination of fasta and tab-delimited input files such as bed files. The dataloader extracts
           regions from the fasta file as defined in the tab-delimited `intervals_file`. Returned sequences are of the type
           np.array([str]).
    args:
        intervals_file:
            doc: bed3+<columns> file path containing intervals + (optionally) labels
            example:
              url: https://raw.githubusercontent.com/kipoi/kipoiseq/master/tests/data/intervals_51bp.tsv
              md5: a76e47b3df87fd514860cf27fdc10eb4
        fasta_file:
            doc: Reference genome FASTA file path.
            example:
              url: https://raw.githubusercontent.com/kipoi/kipoiseq/master/tests/data/hg38_chr22_32000000_32300000.fa
              md5: 01320157a250a3d2eea63e89ecf79eba
        num_chr_fasta:
            doc: True, the the dataloader will make sure that the chromosomes don't start with chr.
        label_dtype:
            doc: None, datatype of the task labels taken from the intervals_file. Example - str, int, float, np.float32
        auto_resize_len:
            doc: None, required sequence length.
        # max_seq_len:
        #     doc: maximum allowed sequence length
        # use_strand:
        #     doc: reverse-complement fasta sequence if bed file defines negative strand
        force_upper:
            doc: Force uppercase output of sequences
        ignore_targets:
            doc: if True, don't return any target variables
    output_schema:
        inputs:
            name: seq
            shape: ()
            doc: DNA sequence as string
            special_type: DNAStringSeq
            associated_metadata: ranges
        targets:
            shape: (None,)
            doc: (optional) values following the bed-entry - chr  start  end  target1   target2 ....
        metadata:
            ranges:
                type: GenomicRanges
                doc: Ranges describing inputs.seq
    postprocessing:
        variant_effects:
          bed_input:
            - intervals_file
    """

    def __init__(self,
                 intervals_file,
                 fasta_file,
                 num_chr_fasta=False,
                 label_dtype=None,
                 auto_resize_len=None,
                 # max_seq_len=None,
                 # use_strand=False,
                 force_upper=True,
                 ignore_targets=False):

        self.num_chr_fasta = num_chr_fasta
        self.intervals_file = intervals_file
        self.fasta_file = fasta_file
        self.auto_resize_len = auto_resize_len
        # self.use_strand = use_strand
        self.force_upper = force_upper
        # self.max_seq_len = max_seq_len

        # if use_strand:
        #     # require a 6-column bed-file if strand is used
        #     bed_columns = 6
        # else:
        #     bed_columns = 3

        self.bed = BedDataset(self.intervals_file,
                              num_chr=self.num_chr_fasta,
                              bed_columns=3,
                              label_dtype=parse_dtype(label_dtype),
                              ignore_targets=ignore_targets)
        self.fasta_extractors = None

    def __len__(self):
        return len(self.bed)

    def __getitem__(self, idx):
        if self.fasta_extractors is None:
            self.fasta_extractors = FastaStringExtractor(self.fasta_file, use_strand=False,  # self.use_strand,
                                                         force_upper=self.force_upper)

        interval, labels = self.bed[idx]

        if self.auto_resize_len:
            # automatically resize the sequence to cerat
            interval = resize_interval(
                interval, self.auto_resize_len, anchor='center')

        # QUESTION: @kromme - why to we need max_seq_len?
        # if self.max_seq_len is not None:
        #     assert interval.stop - interval.start <= self.max_seq_len

        # Run the fasta extractor and transform if necessary
        seq = self.fasta_extractors.extract(interval)

        return {
            "inputs": np.array(seq),
            "targets": labels,
            "metadata": {
                "ranges": GenomicRanges(interval.chrom, interval.start, interval.stop, str(idx))
            }
        }

    @classmethod
    def get_output_schema(cls):
        output_schema = deepcopy(cls.output_schema)
        kwargs = default_kwargs(cls)
        ignore_targets = kwargs['ignore_targets']
        if ignore_targets:
            output_schema.targets = None
        return output_schema


# TODO - properly deal with samples outside of the genome


@kipoi_dataloader(override={"dependencies": deps, 'info.authors': package_authors})
class SeqIntervalDl(Dataset):
    """
    info:
        doc: >
            Dataloader for a combination of fasta and tab-delimited input files such as bed files. The dataloader extracts
            regions from the fasta file as defined in the tab-delimited `intervals_file` and converts them into one-hot encoded
            format. Returned sequences are of the type np.array with the shape inferred from the arguments: `alphabet_axis`
            and `dummy_axis`.
    args:
        intervals_file:
            doc: bed3+<columns> file path containing intervals + (optionally) labels
            example:
              url: https://raw.githubusercontent.com/kipoi/kipoiseq/master/tests/data/intervals_51bp.tsv
              md5: a76e47b3df87fd514860cf27fdc10eb4
        fasta_file:
            doc: Reference genome FASTA file path.
            example:
              url: https://raw.githubusercontent.com/kipoi/kipoiseq/master/tests/data/hg38_chr22_32000000_32300000.fa
              md5: 01320157a250a3d2eea63e89ecf79eba
        num_chr_fasta:
            doc: True, the the dataloader will make sure that the chromosomes don't start with chr.
        label_dtype:
            doc: 'None, datatype of the task labels taken from the intervals_file. Example: str, int, float, np.float32'
        auto_resize_len:
            doc: None, required sequence length.
        # use_strand:
        #     doc: reverse-complement fasta sequence if bed file defines negative strand
        alphabet_axis:
            doc: axis along which the alphabet runs (e.g. A,C,G,T for DNA)
        dummy_axis:
            doc: defines in which dimension a dummy axis should be added. None if no dummy axis is required.
        alphabet:
            doc: >
                alphabet to use for the one-hot encoding. This defines the order of the one-hot encoding.
                Can either be a list or a string: 'ACGT' or ['A, 'C', 'G', 'T']. Default: 'ACGT'
        dtype:
            doc: 'defines the numpy dtype of the returned array. Example: int, np.int32, np.float32, float'
        ignore_targets:
            doc: if True, don't return any target variables

    output_schema:
        inputs:
            name: seq
            shape: (None, 4)
            doc: One-hot encoded DNA sequence
            special_type: DNASeq
            associated_metadata: ranges
        targets:
            shape: (None,)
            doc: (optional) values following the bed-entry - chr  start  end  target1   target2 ....
        metadata:
            ranges:
                type: GenomicRanges
                doc: Ranges describing inputs.seq
    postprocessing:
        variant_effects:
          bed_input:
            - intervals_file
    """

    def __init__(self,
                 intervals_file,
                 fasta_file,
                 num_chr_fasta=False,
                 label_dtype=None,
                 auto_resize_len=None,
                 # max_seq_len=None,
                 # use_strand=False,
                 alphabet_axis=1,
                 dummy_axis=None,
                 alphabet="ACGT",
                 ignore_targets=False,
                 dtype=None):
        # core dataset, not using the one-hot encoding params
        self.seq_dl = StringSeqIntervalDl(intervals_file, fasta_file, num_chr_fasta=num_chr_fasta,
                                          label_dtype=label_dtype, auto_resize_len=auto_resize_len,
                                          # use_strand=use_strand,
                                          ignore_targets=ignore_targets)

        self.input_transform = ReorderedOneHot(alphabet=alphabet,
                                               dtype=dtype,
                                               alphabet_axis=alphabet_axis,
                                               dummy_axis=dummy_axis)

    def __len__(self):
        return len(self.seq_dl)

    def __getitem__(self, idx):
        ret = self.seq_dl[idx]
        ret['inputs'] = self.input_transform(str(ret["inputs"]))
        return ret

    @classmethod
    def get_output_schema(cls):
        """Get the output schema. Overrides the default `cls.output_schema`
        """
        output_schema = deepcopy(cls.output_schema)

        # get the default kwargs
        kwargs = default_kwargs(cls)

        # figure out the input shape
        mock_input_transform = ReorderedOneHot(alphabet=kwargs['alphabet'],
                                               dtype=kwargs['dtype'],
                                               alphabet_axis=kwargs['alphabet_axis'],
                                               dummy_axis=kwargs['dummy_axis'])
        input_shape = mock_input_transform.get_output_shape(
            kwargs['auto_resize_len'])

        # modify it
        output_schema.inputs.shape = input_shape

        # (optionally) get rid of the target shape
        if kwargs['ignore_targets']:
            output_schema.targets = None

        return output_schema
