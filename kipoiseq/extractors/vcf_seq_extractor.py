from cyvcf2 import VCF
from kipoiseq.extractors import FastaSeqExtractor


class IntervalStringBuilder:
    def add():
        pass

    def restore(interval, seq):
        pass

    def concat():
        pass


class VariantSeqExtractor:

    def __init__(self, fasta_file):
        """
        Args:
          fasta_file: path to the fasta file (can be gzipped)
        """
        self.fasta = FastaSeqExtractor(fasta_file)

    def extract(self, interval, variants, anchor, fixed_len=True):
        """

        Args:
          interval: pybedtools.Interval Region of interest from
            which to query the sequence. 0-based
          variants List[cyvcf2.Variant]: variants overlapping the `interval`.
            can also be indels. 1-based
          anchor: position w.r.t. the interval start. (0-based). E.g.
            for an interval of `chr1:10-20` the anchor of 0 denotes
            the point chr1:10 in the 0-based coordinate system. Similarly,
            `anchor=5` means the anchor point is right in the middle
            of the sequence e.g. first half of the sequence (5nt) will be
            upstream of the anchor and the second half (5nt) will be
            downstream of the anchor.
          fixed_len: if True, the return sequence will have the same length
            as the `interval` (e.g. `interval.end - interval.start`)

        Returns:
          A single sequence (`str`) with all the variants applied.
        """
        anchor_abs = anchor + interval.start

        # 1. split the variants hitting the anchor into two
        variants = self._split_overlapping(variants, anchor_abs)

        # 2. split the variants into upstream and downstream
        # 3. sort the variants in each interval
        upstream_variants = sorted(
            filter(variants, lambda x: x.pos < anchor_abs))

        downstream_variants = sorted(
            filter(variants, lambda x: x.pos < anchor_abs), reversed=True)

        # 4. Iterate from the anchor point outwards. At each
        # register the interval from which to take the reference sequence
        # as well as the interval for the variant
        # [
        # (start, end)
        # (seq)
        # (start, end)
        #  ]
        # Also, update the start and end of the interval
        up_sb = IntervalStringBuilder()
        down_sb = IntervalStringBuilder()

        # 5. fetch the sequence
        seq = self.fasta.extract(interval)

        # 6. loop through them, if you see an interval, cut it out
        # from the reference sequence. Store the sequences into a list
        up_sb.restore(interval, seq)
        down_sb.restore(interval, seq)

        # 7. Concate sequences from the upstream and downstream splits. Concat
        # upstream and downstream sequence
        return up_sb.concat() + down_sb.concat()

    def _split_overlapping(variant, anchor_abs):
        pass


class MultiSampleVCF:

    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.vcf = VCF(vcf_file)

        # self._single_sample

    def fetch_variants(self, interval, sample_id=None):
        return self._variant_in(interval, sample_id)

    def _variant_in(self, interval, sample_id=None):
        for v in self.vcf(interval):
            if sample_id is None or self._has_variant(v, sample_id):
                return v

    def _has_variant(self, variant, sample_id):
        # TODO
        return variant.gt_types[self.individual_mapping[sample_id]] != 0


class SingleVariantVCFSeqExtractor:

    def __init__(self, fasta_file, vcf_file):
        '''
        :params strategy: 'one' or 'all'
        '''
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.variant_extractor = VariantSeqExtractor(fasta_file)
        self.vcf = MultiSampleVCF(vcf_file)

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        for variant in self.vcf.fetch_variants(interval, sample_id):
            yield self.variant_extractor.extract(interval,
                                                 variants=[variant],
                                                 anchor=anchor,
                                                 fixed_len=fixed_len)


class SingleSeqVCFSeqExtractor:

    def __init__(self, fasta_file, vcf_file):
        '''
        :params strategy: 'one' or 'all'
        '''
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.variant_extractor = VariantSeqExtractor(fasta_file)
        self.vcf = MultiSampleVCF(vcf_file)

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        return self.variant_extractor.extract(interval,
                                              variants=self.fetch_variants(
                                                  interval, sample_id),
                                              anchor=anchor,
                                              fixed_len=fixed_len)
