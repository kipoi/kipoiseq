import os
from typing import List, Union, Iterable, Iterator
import pandas as pd
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.variant_source import VariantFetcher

try:
    from pyranges import PyRanges
except ImportError:
    from typing import Any

    PyRanges = Any

__all__ = [
    'variants_to_pyranges',
    'BaseVariantMatcher',
    'SingleVariantMatcher',
    'MultiVariantsMatcher',
]


def variants_to_pyranges(variants: List[Variant]) -> PyRanges:
    """
    Create pyrange object given list of variant objects.

    Args:
      variants: list of variant objects have CHROM, POS, REF, ALT properties.
    """
    import pyranges
    df = pd.DataFrame([
        (
            v.chrom,
            v.start,
            v.end,
            v
        )
        for v in variants
    ], columns=['Chromosome', 'Start', 'End', 'variant'])
    return pyranges.PyRanges(df)


def pyranges_to_intervals(pr: PyRanges, interval_attrs: List[str] = None):
    """
    Convert pyranges into list of intervals.


    Args:
      pr: pyranges.PyRanges
      interval_attrs: attribute of interval which should substituted from pr.

    Returns:
      List[Interval]: list of intervals
    """
    interval_attrs = interval_attrs or list()
    for chrom, df in pr:
        for _, row in df.iterrows():
            attrs = {i: row[i] for i in interval_attrs}

            yield Interval(row.Chromosome, row.Start, row.End,
                           strand=row.get('Strand', '.'), attrs=attrs)


def intervals_to_pyranges(intervals):
    """
    Convert list of intervals to pyranges

    Args:
      intervals List[Interval]: list of intervals

    Returns:
      pyranges.Pyranges: Pyranges object.
    """
    import pyranges
    chromosomes, starts, ends, strands = zip(*[
        (i.chrom, i.start, i.end, i.strand)
        for i in intervals
    ])
    return pyranges.PyRanges(
        chromosomes=chromosomes,
        strands=strands,
        starts=starts,
        ends=ends
    )


class PyrangesVariantFetcher(VariantFetcher):

    def __init__(self, variants: List[Variant]):
        self.variants = variants
        self._variants_pr = None

    @property
    def variants_pr(self):
        # convert to PyRanges on demand
        if self._variants_pr is None:
            self._variants_pr = variants_to_pyranges(self.variants)
        return self._variants_pr

    def fetch_variants(self, interval: Union[Interval, Iterable[Interval]]) -> Iterator[Variant]:
        if isinstance(interval, Interval):
            interval = [interval]
        # convert interval(s) to PyRanges object
        interval_pr: PyRanges = intervals_to_pyranges(interval)
        # join with variants
        pr_join = interval_pr.join(self.variants_pr, suffix='_variant')

        yield from pr_join.df["variant"]

    def __iter__(self) -> Iterator[Variant]:
        yield from self.variants


class BaseVariantMatcher:
    """
    Base variant intervals matcher
    """

    def __init__(
            self,
            vcf_file: str = None,
            variants: List[Variant] = None,
            variant_fetcher: VariantFetcher = None,
            gtf_path: str = None,
            bed_path: str = None,
            pranges: PyRanges = None,
            intervals: List[Interval] = None,
            interval_attrs: List[str] = None,
            vcf_lazy: bool = True,
            variant_batch_size: int = 10000
    ):
        """

        Args:
          vcf_file: (optional) path of vcf file
          variants: (optional) readily processed variants
          gtf_path: (optional) path of gtf file contains features
          bed_path: (optional) path of bed file
          pranges: (optional) pyranges object
          intervals: (optional) list of intervals
          interval_attrs: attr of intervals should read from files or
            pyranges object. This argument is not valid with intervals.
            Currently unused
        """
        self.variant_fetcher = self._read_variants(vcf_file, variants, variant_fetcher, vcf_lazy)
        self.interval_attrs = interval_attrs
        self.pr = self._read_intervals(gtf_path, bed_path, pranges,
                                       intervals, interval_attrs, duplicate_attr=True)
        self.variant_batch_size = variant_batch_size

    @staticmethod
    def _read_variants(
            vcf_file=None,
            variants=None,
            variant_fetcher=None,
            vcf_lazy: bool = True,
    ) -> VariantFetcher:
        if vcf_file is not None:
            from kipoiseq.extractors import MultiSampleVCF
            return MultiSampleVCF(vcf_file, lazy=vcf_lazy)
        elif variant_fetcher is not None:
            assert isinstance(variant_fetcher, VariantFetcher), \
                "Wrong type of variant fetcher: %s" % type(variant_fetcher)
            return variant_fetcher
        elif variants is not None:
            return PyrangesVariantFetcher(variants)
        else:
            raise ValueError("No source of variants was specified!")

    @staticmethod
    def _read_intervals(gtf_path=None, bed_path=None, pranges=None,
                        intervals=None, interval_attrs=None, duplicate_attr=False):
        alternatives = [bed_path, pranges, intervals, gtf_path]
        if sum(i is not None for i in alternatives) != 1:
            raise ValueError('only one of `gth_path`, `bed_path`, `pranges`,'
                             '`intervals` or should given as input.')
        if gtf_path:
            import pyranges
            pranges = pyranges.read_gtf(gtf_path, duplicate_attr=duplicate_attr)

        elif bed_path:
            import pyranges
            pranges = pyranges.read_bed(bed_path)

        elif intervals:
            if interval_attrs is not None:
                raise ValueError(
                    '`interval_attrs` is not valid with `intervals`')

            pranges = intervals_to_pyranges(intervals)
            pranges.intervals = intervals

        return pranges

    def __iter__(self):
        raise NotImplementedError()


class SingleVariantMatcher(BaseVariantMatcher):
    """
    Match and iterate variants with intervals.

    """

    def __init__(self, *args, **kwargs):
        """

        Args:
          vcf_file: path of vcf file
          gtf_path: (optional) path of gtf file contains features
          bed_path: (optional) path of bed file
          pranges: (optional) pyranges object
          intervals: (optional) list of intervals
          interval_attrs: attr of intervals should read from files or
            pyranges object. This argument is not valid with intervals.
        """
        super().__init__(*args, **kwargs)

    def _read_vcf_pyranges(self, batch_size=10000):
        """
        Reads vcf and returns batch of pyranges objects.

        Args:
          batch_size: size of each batch.
        """
        for batch in self.variant_fetcher.batch_iter(batch_size):
            yield variants_to_pyranges(batch)

    def iter_pyranges(self) -> PyRanges:
        """

        Iter matched variants with intervals as pyranges.

        Returns:

        """
        for pr_variants in self._read_vcf_pyranges():
            pr_join = self.pr.join(pr_variants, suffix='_variant')
            if not hasattr(pr_join, 'intervals'):
                pr_join.intervals = list(pyranges_to_intervals(
                    pr_join, interval_attrs=self.interval_attrs))
            yield pr_join

    def iter_rows(self):
        """
        Iter matched variants with intervals as pandas series.
        """
        for pr in self.iter_pyranges():
            for _, df in pr:
                df = df.sort_values(['Start_variant', 'Start'])
                for _, row in df.iterrows():
                    yield row

    def __iter__(self) -> (Interval, Variant):
        """
        Iterate interval and variant object.
        """
        for row in self.iter_rows():
            yield row['intervals'], row['variant']


class MultiVariantsMatcher(BaseVariantMatcher):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if hasattr(self.pr, 'intervals'):
            self.intervals = self.pr.intervals
        else:
            self.intervals = pyranges_to_intervals(self.pr)

    def __iter__(self):
        for i in self.intervals:
            yield i, self.variant_fetcher.fetch_variants(i)
