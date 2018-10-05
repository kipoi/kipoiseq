import pandas as pd
import numpy as np
import pyfaidx
from pybedtools import BedTool, Interval


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


class FastaStringExtractor(object):
    dtype = np.float32

    def __init__(self, datafile, use_strand=False, force_upper=False):
        """Fasta file extractor

        NOTE: The extractor is not thread-save.
        If you with to use it with multiprocessing,
        create a new extractor object in each process.

        Args:
          datafile (str): path to the bigwig file
          use_strand (bool): if True, the extracted sequence
            is reverse complemented in case interval.strand == "-"
          force_upper (bool): Force uppercase output
        """
        self._datafile = datafile
        self.use_strand = use_strand
        self.fasta = pyfaidx.Fasta(self._datafile)
        self.force_upper = force_upper

    def _extract(self, intervals):
        seqs = []
        for index, interval in enumerate(intervals):
            chrom = interval.chrom
            # reverse-complement seq the negative strand
            rc = self.use_strand and interval.strand == "-"
            # pyfaidx wants a 1-based interval
            seq = str(self.fasta.get_seq(chrom, interval.start + 1, interval.stop, rc=rc).seq)
            if self.force_upper:
                seq = seq.upper()

            seqs.append(seq)
        return seqs

    def __call__(self, intervals):
        return self._extract(intervals)
