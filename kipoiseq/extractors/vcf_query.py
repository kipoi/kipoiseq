import abc
from itertools import islice
from typing import Tuple, Iterable, List
from tqdm import tqdm
from kipoiseq.dataclasses import Variant, Interval

__all__ = [
    'BaseVariantQuery',
    'VariantQuery',
    'FilterVariantQuery',
    'BaseVariantIntervalQuery',
    'VariantIntervalQuery',
    'NumberVariantQuery',
    'VariantIntervalQueryable',
]


class BaseVariantQuery:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __call__(self, variant: Variant):
        raise NotImplementedError

    def __or__(self, other):
        return VariantQuery(lambda variant: self(variant) or other(variant))

    def __and__(self, other):
        return VariantQuery(lambda variant: self(variant) and other(variant))


class VariantQuery(BaseVariantQuery):

    def __init__(self, func):
        self.func = func

    def __call__(self, variant: Variant):
        return self.func(variant)

    def __or__(self, other):
        return VariantQuery(lambda variant: self(variant) or other(variant))

    def __and__(self, other):
        return VariantQuery(lambda variant: self(variant) and other(variant))


class FilterVariantQuery(BaseVariantQuery):

    def __init__(self, filter='PASS'):
        self.filter = filter

    def __call__(self, variant):
        return variant.filter == self.filter


class BaseVariantIntervalQuery:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __call__(self, variants: List[Variant], interval: Interval):
        raise NotImplementedError

    def __or__(self, other):
        return VariantIntervalQuery(
            lambda variants, interval: [
                i or j for i, j in zip(
                    self(variants, interval),
                    other(variants, interval)
                )
            ])

    def __and__(self, other):
        return VariantIntervalQuery(
            lambda variants, interval: [
                i and j for i, j in zip(
                    self(variants, interval),
                    other(variants, interval)
                )
            ])


class VariantIntervalQuery(BaseVariantIntervalQuery):

    def __init__(self, func):
        self.func = func

    def __call__(self, variants: List[Variant], interval: Interval):
        return self.func(variants, interval)


class NumberVariantQuery(BaseVariantIntervalQuery):
    """
    Closure for variant query. Filter variants for interval
      if number of variants in given limits.
    """

    def __init__(self, max_num=float('inf'), min_num=0):
        # TODO: sample speficity
        self.max_num = max_num
        self.min_num = min_num

    def __call__(self, variants, interval):
        if self.max_num >= len(variants) >= self.min_num:
            return [True] * len(variants)
        else:
            return [False] * len(variants)


_VariantIntervalType = List[Tuple[Iterable[Variant], Interval]]


class VariantIntervalQueryable:

    def __init__(self, vcf, variant_intervals: _VariantIntervalType,
                 progress=False):
        """
        Query object of variants.

        Args:
          vcf: cyvcf2.VCF objects.
          variants: iter of (variant, interval) tuples.
        """
        self.vcf = vcf
        self.variant_intervals = variant_intervals
        self.progress = progress

    def __iter__(self):
        variant_intervals = tqdm(self.variant_intervals) if self.progress \
            else self.variant_intervals

        for variants, interval in variant_intervals:
            yield from variants

    def iter_intervals(self):
        for _, interval in self.variant_intervals:
            yield interval

    def filter(self, query: VariantQuery):
        """
        Filters variant given conduction.

        Args:
          query: function which get a variant as input and filtered iter of
            variants.
        """
        self.variant_intervals = [
            (filter(query, variants), Interval)
            for variants, interval in self.variant_intervals
        ]
        return self

    def filter_range(self, query: VariantIntervalQuery):
        """
        Filters variant given conduction.

        Args:
          query: function which get variants and an interval as input
            and filtered iter of variants.
        """
        self.variant_intervals = list(self._filter_range(query))
        return self

    def _filter_range(self, query: VariantIntervalQuery):
        for variants, interval in self.variant_intervals:
            variants = list(variants)
            yield (
                v
                for v, cond in zip(variants, query(variants, interval))
                if cond
            ), interval

    def batch_iter(self, batch_size=10000) -> Iterable:
        """Iterates variatns and interval in query object.

        # Arguments
            batch_size: size of each batch.
        """
        variants = iter(self)
        batch = list()

        for v in variants:
            if len(batch) == 0:
                # 1-based variant to 0-based interval
                interval = Interval(v.chrom, v.pos - 1, v.pos)

            elif v.chrom != interval.chrom:
                # if new chrom return current batch and create new batch
                yield VariantIntervalQueryable(self.vcf, [(batch, interval)])
                batch = list()
                interval = Interval(v.chrom, v.pos - 1, v.pos)

            interval._end = v.pos
            batch.append(v)

            if len(batch) == batch_size:
                yield VariantIntervalQueryable(self.vcf, [(batch, interval)])
                batch = list()

        if len(batch) > 0:
            yield VariantIntervalQueryable(self.vcf, [(batch, interval)])

    def to_vcf(self, path, remove_samples=False, clean_info=False):
        """
        Parse query result as vcf file.

        Args:
          path: path of the file.
          remove_samples: remove sample columns from vcf file
          clean_info: clean info fields from vcf file
        """
        from cyvcf2 import Writer

        header = self.vcf.raw_header

        # remove sample columns  from header
        # 2th last column in the header is columns
        # all the columns after 8th is for samples
        if remove_samples:
            header = header.split('\n')
            columns = header[-2].strip().split('\t')
            header[-2] = '\t'.join(columns[:8])
            header = '\n'.join(header)

        writer = Writer.from_string(path, header)

        for v in self:
            variant = v.source

            if remove_samples or clean_info:
                variant = str(variant).strip().split()

                if remove_samples:
                    # all the columns after 8th is for samples so remove them
                    variant = variant[:8]

                if clean_info:
                    # 7th column contains info fields
                    # replace it with N/A
                    variant[7] = '.'

                variant = writer.variant_from_string('\t'.join(variant))

            writer.write_record(variant)
