import abc
from typing import List, Union

from kipoiseq import Interval
from kipoiseq.transforms.functional import rc_dna

import logging

log = logging.getLogger(__name__)

__all__ = [
    "BaseExtractor",
    "BaseMultiIntervalFetcher",
    "BaseMultiIntervalSeqExtractor",
    "GenericMultiIntervalSeqExtractor",
    "BaseVCFSeqExtractor"
]

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kipoiseq.extractors.vcf_matching import BaseVariantMatcher
    from kipoiseq.extractors.vcf_query import VariantIntervalQueryable, MultiSampleVCF


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

    def get_seq(self, key: Union[object, List[object]]):
        """
        Get one or multiple sequences by key.
        """
        if isinstance(key, list):
            return [self.sel(k) for k in key]
        else:
            return self.sel(key)

    def iter_seq(self, key: Union[object, List[object]]):
        """
        iterate over (list of) keys
        """
        if isinstance(key, list):
            for k in key:
                yield self.sel(k)
        else:
            yield self.sel(key)

    def __len__(self):
        return len(self.interval_fetcher)

    def __getitem__(self, idx):
        return self.isel(idx)

    def items(self):
        for i in self.keys():
            yield i, self.sel(i)

    def extract_all(self):
        """
        Extract all sequences; alias for self.items()
        """
        yield from self.items()

    # def __iter__(self):
    #     for i in self.keys():
    #         yield self.sel(i)


class BaseVCFSeqExtractor(GenericMultiIntervalSeqExtractor, metaclass=abc.ABCMeta):

    def __init__(
            self,
            interval_fetcher: BaseMultiIntervalFetcher,
            reference_seq_extractor: BaseExtractor,
            variant_matcher: 'BaseVariantMatcher',
            multi_sample_VCF: 'MultiSampleVCF',
    ):
        self.variant_matcher = variant_matcher
        self.multi_sample_VCF = multi_sample_VCF

        # takes reference sequence, interval + variants to extract the alternative sequence
        from kipoiseq.extractors import VariantSeqExtractor
        self.variant_seq_extractor = VariantSeqExtractor(reference_sequence=reference_seq_extractor)

        super().__init__(interval_fetcher=interval_fetcher, extractor=reference_seq_extractor)

    @property
    def reference_seq_extractor(self):
        return self.extractor

    @classmethod
    def from_pyranges(
            cls,
            pyranges,
            concat_by: Union[str, List[str]] = "transcript_id",
            **kwargs
    ):
        from kipoiseq.extractors import GTFMultiIntervalFetcher
        interval_fetcher = GTFMultiIntervalFetcher(region_df=pyranges.df.set_index(concat_by))
        return cls(interval_fetcher=interval_fetcher, **kwargs)

    @staticmethod
    def _unstrand(intervals: List[Interval]):
        """
        Set strand of list of intervals to default - '.'
        """
        return [i.unstrand() for i in intervals]

    def extract(self, intervals: List[Interval], sample_id=None, *args, **kwargs):
        """
        Extract cds with variants in their dna sequence. It depends on the
        child class if a sequence have all variants inserted or only one variant
        is inserted per dna sequence
        :param intervals_df: GTF dataframe as given by PyRanges containing at least the following columns:
            - Chromosome
            - Start
            - End
            - Strand
        :param sample_id: optional sample id
        :return: alternative sequence with variant information
        """
        reverse_complement = False
        if self.use_strand:
            reverse_complement = intervals[0].neg_strand
        intervals = self._unstrand(intervals)

        variant_interval_queryable = self.multi_sample_VCF.query_variants(intervals, sample_id=sample_id)

        # extract_query can return one sequence per variant or one sequence per sample
        iter_seqs = self.extract_query(variant_interval_queryable, sample_id=sample_id)
        for seqs, variant_info in iter_seqs:
            # 1st seq, 2nd variant info
            yield self._prepare_seq(seqs, reverse_complement), variant_info

        return None

    @abc.abstractmethod
    def extract_query(self, variant_interval_queryable: 'VariantIntervalQueryable', sample_id=None):
        raise NotImplementedError()

    def _reference_sequence(self, variant_interval_queryable):
        """
        Extract reference sequence without variants
        """
        intervals = variant_interval_queryable.iter_intervals()
        return [self.reference_seq_extractor.extract(interval) for interval in intervals]

    def _filter_snv(self, variants):
        for variant in variants:
            if len(variant.ref) == len(variant.alt) == 1:  # only SOVs supported
                yield variant
            elif len(variant.ref) == len(variant.alt) > 1:
                log.warning('Current version of extractor works only for len(variant.ref)'
                            ' == len(variant.alt) == 1, but the len was: ' + str(len(variant.alt)))
            else:
                log.warning('Current version of extractor ignores indel'
                            ' to avoid shift in frame')
