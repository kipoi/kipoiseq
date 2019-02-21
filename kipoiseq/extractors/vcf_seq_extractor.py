from cyvcf2 import VCF
from pyfaidx import Sequence
from pybedtools import Interval
from kipoiseq.extractors import BaseExtractor, FastaStringExtractor


__all__ = [
    'VariantSeqExtractor',
    'MultiSampleVCF',
    'SingleVariantVCFSeqExtractor',
    'SingleSeqVCFSeqExtractor'
]


class MultiSampleVCF(VCF):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sample_mapping = dict(zip(self.samples, range(len(self.samples))))

    def fetch_variants(self, interval, sample_id=None):
        for v in self(self._region(interval)):
            if sample_id is None or self._has_variant(v, sample_id):
                yield v

    def _region(self, interval):
        return '%s:%d-%d' % (interval.chrom, interval.start, interval.end)

    def _has_variant(self, variant, sample_id):
        return variant.gt_types[self.sample_mapping[sample_id]] != 0


class IntervalSeqBuilder(list):
    """
    String builder for `pyfaidx.Sequence` and `Interval` objects.
    """

    def restore(self, sequence):
        """
        Args:
          seq: `pyfaidx.Sequence` which convert all interval inside
            to `Seqeunce` objects.
        """
        for i, interval in enumerate(self):
            if type(self[i]) == Interval:
                start = interval.start - sequence.start
                end = start + interval.length
                self[i] = sequence[start: end]

    def _concat(self):
        for sequence in self:
            if type(sequence) != Sequence:
                raise TypeError('Intervals should be restored with `restore`'
                                ' method before calling concat method!')
            yield sequence.seq

    def concat(self):
        """
        Build the string from sequence objects.

        Returns:
          str: the final sequence.
        """
        return ''.join(self._concat())


class VariantSeqExtractor(BaseExtractor):

    def __init__(self, fasta_file):
        """
        Args:
          fasta_file: path to the fasta file (can be gzipped)
        """
        self.fasta = FastaStringExtractor(fasta_file)

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
        pairs_seq = zip(*self._variant_to_sequence(variants))

        anchor_abs = anchor + interval.start
        variants = self._split_overlapping(pairs_seq, anchor_abs)

        # 2. split the variants into upstream and downstream
        # 3. sort the variants in each interval
        # upstream_variants = sorted(
        #     filter(variants, lambda x: x.pos < anchor_abs))

        # downstream_variants = sorted(
        #     filter(variants, lambda x: x.pos < anchor_abs), reversed=True)

        # 4. Iterate from the anchor point outwards. At each
        # register the interval from which to take the reference sequence
        # as well as the interval for the variant
        # [
        # (start, end)
        # (seq)
        # (start, end)
        #  ]
        # Also, update the start and end of the interval
        # up_sb = IntervalSeqBuilder()
        # down_sb = IntervalSeqBuilder()

        # 5. fetch the sequence
        # seq = self.fasta.extract(interval)

        # 6. loop through them, if you see an interval, cut it out
        # from the reference sequence. Store the sequences into a list
        # up_sb.restore(interval, seq)
        # down_sb.restore(interval, seq)

        # 7. Concate sequences from the upstream and downstream splits. Concat
        # upstream and downstream sequence
        # return up_sb.concat() + down_sb.concat()

    def _variant_to_sequence(self, variants):
        """
        Convert `cyvcf2.Variant` objects to `pyfaidx.Seqeunce` objects
        for reference and variants.
        """
        for v in variants:
            ref = Sequence(name=v.CHROM, seq=v.REF,
                           start=v.POS, end=v.POS + len(v.REF))
            # TO DO: consider alternative alleles.
            alt = Sequence(name=v.CHROM, seq=v.ALT[0],
                           start=v.POS, end=v.POS + len(v.ALT[0]))
            yield ref, alt

    def _split_overlapping(pairs_seq, anchor_abs):
        """
        Split the variants hitting the anchor into two
        """
        pass


class SingleVariantVCFSeqExtractor(BaseExtractor):

    def __init__(self, fasta_file, vcf_file):
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


class SingleSeqVCFSeqExtractor(BaseExtractor):

    def __init__(self, fasta_file, vcf_file):
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
