import abc
from typing import Iterable, Union, Iterator

from kipoiseq import Interval, Variant
from kipoiseq.utils import batch_iter


class VariantFetcher(Iterable, metaclass=abc.ABCMeta):
    """
    Base class for all variant-returning data sources.
    """

    @abc.abstractmethod
    def fetch_variants(self, interval: Union[Interval, Iterable[Interval]]) -> Iterator[Variant]:
        """
        Fetch variants that intersect the provided interval(s).

        :param interval:  One or multiple Interval objects
        :return: Iterator of Variant objects that intersect `interval`
        """
        raise NotImplementedError

    def batch_iter(self, batch_size=10000) -> Iterator[Iterable[Variant]]:
        """
        Fetch variants in batches.

        :param batch_size: Number of variants per batch
        :return: Iterator that yields batches of Variant objects
        """
        yield from batch_iter(self, batch_size)

    @abc.abstractmethod
    def __iter__(self) -> Iterator[Variant]:
        """
        Fetch variants in batches.

        :return: Iterator of Variant objects
        """
        raise NotImplementedError
