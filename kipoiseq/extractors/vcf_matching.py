from typing import List
import pandas as pd
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors import MultiSampleVCF

# pyranges is optional
try:
    import pyranges
except ImportError:
    pyranges = None

__all__ = [
    'variants_to_pyranges',
    'SingleVariantMatcher',
    'MultiVariantsMatcher'
]


def variants_to_pyranges(variants: List[Variant]) -> pyranges.PyRanges:
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
            v.start + max(len(v.ref), len(v.alt)),
            v
        )
        for v in variants
    ], columns=['Chromosome', 'Start', 'End', 'variant'])
    return pyranges.PyRanges(df)


def pyranges_to_intervals(pr: pyranges.PyRanges, interval_attrs: List[str] = None):
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


class BaseVariantMatcher:
    """
    Base variant intervals matcher
    """

    def __init__(
            self,
            vcf_file: str,
            gtf_path: str = None,
            bed_path: str = None,
            pranges: pyranges.PyRanges = None,
            intervals: List[Interval] = None,
            interval_attrs: List[str] = None,
            vcf_lazy: bool = True,
            variant_batch_size: int = 10000
    ):
        """

        Args:
          vcf_file: path of vcf file
          gtf_path: (optional) path of gtf file contains features
          bed_path: (optional) path of bed file
          pranges: (optional) pyranges object
          intervals: (optional) list of intervals
          interval_attrs: attr of intervals should read from files or
            pyranges object. This argument is not valid with intervals.
            Currently unused
        """
        self.vcf = MultiSampleVCF(vcf_file, lazy=vcf_lazy)
        self.interval_attrs = interval_attrs
        self.pr = self._read_intervals(gtf_path, bed_path, pranges,
                                       intervals, interval_attrs)
        self.variant_batch_size = variant_batch_size

    @staticmethod
    def _read_intervals(gtf_path=None, bed_path=None, pranges=None,
                        intervals=None, interval_attrs=None):
        alternatives = [bed_path, pranges, intervals, gtf_path]
        if sum(i is not None for i in alternatives) != 1:
            raise ValueError('only one of `gth_path`, `bed_path`, `pranges`,'
                             '`intervals` or should given as input.')
        if gtf_path:
            import pyranges
            pranges = pyranges.read_gtf(gtf_path)

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
        for batch in self.vcf.batch_iter(batch_size):
            yield variants_to_pyranges(batch)

    def iter_pyranges(self) -> pyranges.PyRanges:
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

    def __iter__(self):
        """
        Iterate interval and variant object.
        """
        for row in self.iter_rows():
            yield row['intervals'], row['variant']


class MultiVariantsMatcher(BaseVariantMatcher):

    def __init__(
            self,
            vcf_file,
            gtf_path=None,
            pranges=None,
            intervals=None,
            interval_attrs=None,
            vcf_lazy=True,
            variant_batch_size=10000
    ):
        super().__init__(
            vcf_file,
            gtf_path=gtf_path,
            pranges=pranges,
            intervals=intervals,
            interval_attrs=interval_attrs,
            vcf_lazy=vcf_lazy,
            variant_batch_size=variant_batch_size
        )
        if hasattr(self.pr, 'intervals'):
            self.intervals = self.pr.intervals
        else:
            self.intervals = pyranges_to_intervals(self.pr)

    def __iter__(self):
        for i in self.intervals:
            yield i, self.vcf.fetch_variants(i)
