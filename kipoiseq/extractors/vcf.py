from itertools import islice
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_query import VariantIntervalQueryable

try:
    from cyvcf2 import VCF
except ImportError:
    VCF = object


__all__ = [
    'MultiSampleVCF'
]


class MultiSampleVCF(VCF):

    def __init__(self, *args, **kwargs):
        from cyvcf2 import VCF
        super(MultiSampleVCF, self).__init__(*args, **kwargs)
        self.sample_mapping = dict(zip(self.samples, range(len(self.samples))))

    def fetch_variants(self, interval, sample_id=None):
        for v in self(self._region(interval)):
            v = Variant.from_cyvcf(v)
            if sample_id is None or self.has_variant(v, sample_id):
                yield v

    def _region(self, interval):
        return '%s:%d-%d' % (interval.chrom, interval.start, interval.end)

    def has_variant(self, variant, sample_id):
        gt_type = variant.source.gt_types[self.sample_mapping[sample_id]]
        return self._has_variant_gt(gt_type)

    def _has_variant_gt(self, gt_type):
        return gt_type != 0 and gt_type != 2

    def __next__(self):
        return Variant.from_cyvcf(super().__next__())

    def batch_iter(self, batch_size=10000):
        """Iterates variatns in vcf file.

        # Arguments
            vcf_file: path of vcf file.
            batch_size: size of each batch.
        """
        variants = iter(self)
        batch = list(islice(variants, batch_size))

        while batch:
            yield batch
            batch = list(islice(variants, batch_size))

    def query_variants(self, intervals, sample_id=None, progress=False):
        """Fetch variants for given multi-intervals from vcf file
        for sample if sample id is given.

        # Arguments
            intervals (List[pybedtools.Interval]): list of Interval objects
            sample_id (str, optional): sample id in vcf file.

        # Returns
          VCFQueryable: queryable object whihc allow you to query the
            fetched variatns.

        # Example
            To fetch variants if only single variant present in interval.
            ```
              >>> MultiSampleVCF(vcf_path) \
                    .query_variants(intervals) \
                    .filter(lambda variant: variant.qual > 10) \
                    .filter_range(NumberVariantQuery(max_num=1)) \
                    .to_vcf(output_path)
            ```
        """
        pairs = ((self.fetch_variants(i, sample_id=sample_id), i)
                 for i in intervals)
        return VariantIntervalQueryable(self, pairs, progress=progress)

    def get_variant(self, variant):
        """Returns variant from vcf file. Let you use vcf file as dict.

        # Arguments:
            vcf: cyvcf2.VCF file
            variant: variant object or variant id as string.

        # Returns
            Variant object.

        # Example
            ```python
              >>> MultiSampleVCF(vcf_path).get_variant("chr1:4:T:['C']")
            ```
        """
        if type(variant) == str:
            variant = Variant.from_str(variant)

        variants = self.fetch_variants(
            Interval(variant.chrom, variant.pos, variant.pos))
        for v in variants:
            if v.ref == variant.ref and v.alt == variant.alt:
                return v
        raise KeyError('Variant %s not found in vcf file.' % str(variant))

    def get_samples(self, variant):
        """Fetchs sample names which have given variants

        # Arguments
            variant: variant object.

        # Returns
            Dict[str, int]: Dict of sample which have variant and gt as value.
        """
        return dict(filter(lambda x: self._has_variant_gt(x[1]),
                           zip(self.samples, variant.source.gt_types)))
