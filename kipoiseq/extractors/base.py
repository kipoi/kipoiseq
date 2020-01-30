import abc

from kipoiseq import Interval

__all__ = ["BaseExtractor", "FastaStringExtractor"]  # "BigWigExtractor"]


class BaseExtractor(object):
    __metaclass__ = abc.ABCMeta

    # main method
    @abc.abstractmethod
    def extract(self, interval: Interval, *args, **kwargs) -> str:
        raise NotImplementedError

    # closing files
    def __del__(self):
        return self.close()

    def close(self):
        # implemented by the subclass
        pass


class FastaStringExtractor(BaseExtractor):
    """Fasta file extractor

    NOTE: The extractor is not thread-save.
    If you with to use it with multiprocessing,
    create a new extractor object in each process.

    # Arguments
      fasta_file (str): path to the fasta_file
      use_strand (bool): if True, the extracted sequence
        is reverse complemented in case interval.strand == "-"
      force_upper (bool): Force uppercase output
    """

    def __init__(self, fasta_file, use_strand=False, force_upper=False):
        from pyfaidx import Fasta

        self.fasta_file = fasta_file
        self.use_strand = use_strand
        self.fasta = Fasta(self.fasta_file)
        self.force_upper = force_upper

    def extract(self, interval: Interval, **kwargs) -> str:
        """
        Returns the FASTA sequence in some given interval as string

        Args:
            interval: the interval to query
            **kwargs:

        Returns:
            sequence of requested interval

        """
        # reverse-complement seq the negative strand
        rc = self.use_strand and interval.strand == "-"

        # pyfaidx wants a 1-based interval
        seq = str(self.fasta.get_seq(
            interval.chrom,
            interval.start + 1,
            interval.stop,
            rc=rc
        ).seq)

        # optionally, force upper-case letters
        if self.force_upper:
            seq = seq.upper()
        return seq

    def close(self):
        return self.fasta.close()

# class BigWigExtractor(BaseExtractor):
#     """Big-wig file extractor

#     NOTE: The extractor is not thread-save.
#     If you with to use it with multiprocessing,
#     create a new extractor object in each process.

#     # Arguments
#       bigwig_file: path to the bigwig file
#     """

#     def __init__(self, bigwig_file):
#         from genomelake.extractors import BigwigExtractor

#         self.bigwig_file = bigwig_file
#         self.batch_extractor = BigwigExtractor(self.bigwig_file)

#     def extract(self, interval):
#         return self.batch_extractor([interval])[0]

#     def close(self):
#         return self.batch_extractor.close()
