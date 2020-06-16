import abc
from typing import List, Union

from kipoiseq import Interval

__all__ = [
    "BaseExtractor",
    "BaseMultiIntervalFetcher",
    "BaseMultiIntervalSeqExtractor",
    "GenericMultiIntervalSeqExtractor",
]

from kipoiseq.transforms.functional import rc_dna


class BaseExtractor(object):
    __metaclass__ = abc.ABCMeta

    _use_strand: bool

    # main method
    @abc.abstractmethod
    def extract(self, interval: Interval, *args, **kwargs) -> str:
        raise NotImplementedError

    @property
    def use_strand(self):
        return self._use_strand

    # closing files
    def __del__(self):
        return self.close()

    def close(self):
        # implemented by the subclass
        pass


class BaseMultiIntervalSeqExtractor:

    def __init__(self, extractor: BaseExtractor, use_strand=True):
        self.extractor = extractor
        self._use_strand = use_strand

    @property
    def use_strand(self):
        return self._use_strand

    # def __len__(self):
    #     return len(self.interval_fetcher)

    @classmethod
    def _prepare_seq(
            cls,
            seqs: List[str],
            # intervals: List[Interval],
            reverse_complement: Union[str, bool],
            *args,
            **kwargs
    ) -> str:
        """
        Prepare the dna sequence in the final variant, which should be translated in amino acid sequence

        :param seqs: current dna sequence
        :param intervals: the list of intervals corresponding to the sequence snippets
        :param reverse_complement: should the dna be reverse-complemented?
        """
        seq = "".join(seqs)
        if reverse_complement is True or reverse_complement == "-":
            # optionally reverse complement
            seq = rc_dna(seq)
        return seq

    def extract(self, intervals: List[Interval], *args, **kwargs):
        """
        Extract and concatenate the sequence for a list of intervals
        :param intervals: list of intervals
        :param kwargs: will be passed on to `_prepare_seq()`
        :return: concatenated dna sequence of the intervals
        """
        seqs = [self.extractor.extract(i) for i in intervals]

        reverse_strand = False
        if self.use_strand:
            reverse_strand = intervals[0].neg_strand
            if self.extractor.use_strand:
                # If the fasta extractor already does the reversion, we do not have to do it manually.
                # Instead, we need to reverse the list of sequences.
                seqs.reverse()
                reverse_strand = False

        return self._prepare_seq(seqs, reverse_strand, **kwargs)


class BaseMultiIntervalFetcher(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def keys(self) -> List[object]:
        raise NotImplementedError()

    @abc.abstractmethod
    def get_intervals(self, key) -> List[Interval]:
        """
        Returns a list of Intervals on request of some key
        :param key: The identifier for the requested intervals
        :return: list of intervals
        """
        raise NotImplementedError()

    def sel(self, *args, **kwargs):
        """
        alias for get_intervals
        """
        return self.get_intervals(*args, **kwargs)

    def isel(self, idx):
        """
        Get interval by index
        :param idx: index in the range (0, len(self))
        """
        return self.get_intervals(self.keys()[idx])

    def __len__(self):
        return len(self.keys())

    def __getitem__(self, idx):
        return self.isel(idx)

    def items(self):
        for i in self.keys():
            yield i, self.get_intervals(i)

    # def __iter__(self):
    #     for i in self.keys():
    #         yield self.get_intervals(i)


class GenericMultiIntervalSeqExtractor(BaseMultiIntervalSeqExtractor):

    def __init__(self, extractor: BaseExtractor, interval_fetcher: BaseMultiIntervalFetcher, use_strand=True):
        super().__init__(
            extractor=extractor,
            use_strand=use_strand,
        )
        self.interval_fetcher = interval_fetcher

    def keys(self) -> List[object]:
        return self.interval_fetcher.keys()

    def get_intervals(self, key) -> List[Interval]:
        """
        Returns a list of Intervals on request of some key
        :param key: The identifier for the requested intervals
        :return: list of intervals
        """
        return self.interval_fetcher.get_intervals(key)

    def sel(self, *args, **kwargs):
        """
        extract sequence by identifier
        """
        return self.extract(self.interval_fetcher.sel(*args, **kwargs))

    def isel(self, idx):
        """
        extract sequence by  index
        :param idx: index in the range (0, len(self))
        """
        return self.extract(self.interval_fetcher.isel(idx))

    def __len__(self):
        return len(self.interval_fetcher)

    def __getitem__(self, idx):
        return self.isel(idx)

    def items(self):
        for i in self.keys():
            yield i, self.sel(i)

    # def __iter__(self):
    #     for i in self.keys():
    #         yield self.sel(i)
