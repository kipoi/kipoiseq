from pybedtools import Interval
from pyfaidx import Sequence, complement
from kipoiseq.extractors import BaseExtractor, FastaStringExtractor
try:
    from cyvcf2 import VCF
except ImportError:
    VCF = object

__all__ = [
    'VariantSeqExtractor',
    'MultiSampleVCF',
    'SingleVariantVCFSeqExtractor',
    'SingleSeqVCFSeqExtractor'
]


class MultiSampleVCF(VCF):

    def __init__(self, *args, **kwargs):
        from cyvcf2 import VCF
        super(MultiSampleVCF, self).__init__(*args, **kwargs)
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
            # interval.end can be bigger than interval.start
            interval_len = max(0, interval.end - interval.start)

            if type(self[i]) == Interval:
                start = interval.start - sequence.start
                end = start + interval_len
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
        self.fasta = FastaStringExtractor(fasta_file, use_strand=True)

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
        # Preprocessing
        anchor = max(min(anchor, interval.end), interval.start)
        variant_pairs = self._variant_to_sequence(variants)

        # 1. Split variants overlapping with anchor
        variant_pairs = list(self._split_overlapping(variant_pairs, anchor))

        # 2. split the variants into upstream and downstream
        # and sort the variants in each interval
        upstream_variants = sorted(
            filter(lambda x: x[0].start >= anchor, variant_pairs),
            key=lambda x: x[0].start)

        downstream_variants = sorted(
            filter(lambda x: x[0].start < anchor, variant_pairs),
            key=lambda x: x[0].start, reverse=True)

        # 3. Extend start and end position for deletions
        if fixed_len:
            istart, iend = self._updated_interval(
                interval, upstream_variants, downstream_variants)
        else:
            istart, iend = interval.start, interval.end

        # 4. Iterate from the anchor point outwards. At each
        # register the interval from which to take the reference sequence
        # as well as the interval for the variant
        down_sb = self._downstream_builder(
            downstream_variants, interval, anchor, istart)

        up_sb = self._upstream_builder(
            upstream_variants, interval, anchor, iend)

        # 5. fetch the sequence and restore intervals in builder
        seq = self._fetch(interval, istart, iend)
        up_sb.restore(seq)
        down_sb.restore(seq)

        # 6. Concate sequences from the upstream and downstream splits. Concat
        # upstream and downstream sequence. Cut to fix the length.
        down_str = down_sb.concat()
        up_str = up_sb.concat()

        if fixed_len:
            down_str, up_str = self._cut_to_fix_len(
                down_str, up_str, interval, anchor)

        seq = down_str + up_str

        if interval.strand == '-':
            seq = complement(seq)[::-1]

        return seq

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

    def _split_overlapping(self, variant_pairs, anchor):
        """
        Split the variants hitting the anchor into two
        """
        for ref, alt in variant_pairs:
            if ref.start < anchor < ref.end or alt.start < anchor < alt.end:
                mid = anchor - ref.start
                yield ref[:mid], alt[:mid]
                yield ref[mid:], alt[mid:]
            else:
                yield ref, alt

    def _updated_interval(self, interval, up_variants, down_variants):
        istart = interval.start
        iend = interval.end

        for ref, alt in up_variants:
            diff_len = len(alt) - len(ref)
            if diff_len < 0:
                iend -= diff_len

        for ref, alt in down_variants:
            diff_len = len(alt) - len(ref)
            if diff_len < 0:
                istart += diff_len

        return istart, iend

    def _downstream_builder(self, down_variants, interval, anchor, istart):
        down_sb = IntervalSeqBuilder()

        prev = anchor
        for ref, alt in down_variants:
            if ref.end <= istart:
                break
            down_sb.append(Interval(interval.chrom, ref.end, prev))
            down_sb.append(alt)
            prev = ref.start
        down_sb.append(Interval(interval.chrom, istart, prev))
        down_sb.reverse()

        return down_sb

    def _upstream_builder(self, up_variants, interval, anchor, iend):
        up_sb = IntervalSeqBuilder()

        prev = anchor
        for ref, alt in up_variants:
            if ref.start > iend:
                break
            up_sb.append(Interval(interval.chrom, prev, ref.start))
            up_sb.append(alt)
            prev = ref.end
        up_sb.append(Interval(interval.chrom, prev, iend))

        return up_sb

    def _fetch(self, interval, istart, iend):
        seq = self.fasta.extract(Interval(interval.chrom, istart, iend))
        seq = Sequence(name=interval.chrom, seq=seq, start=istart, end=iend)
        return seq

    def _cut_to_fix_len(self,  down_str, up_str, interval, anchor):
        down_len = anchor - interval.start
        up_len = interval.end - anchor
        down_str = down_str[-down_len:] if down_len else ''
        up_str = up_str[:up_len] if up_len else ''
        return down_str, up_str


class BaseVCFSeqExtractor(BaseExtractor):
    def __init__(self, fasta_file, vcf_file):
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.variant_extractor = VariantSeqExtractor(fasta_file)
        self.vcf = MultiSampleVCF(vcf_file)


class SingleVariantVCFSeqExtractor(BaseVCFSeqExtractor):

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        for variant in self.vcf.fetch_variants(interval, sample_id):
            yield self.variant_extractor.extract(interval,
                                                 variants=[variant],
                                                 anchor=anchor,
                                                 fixed_len=fixed_len)


class SingleSeqVCFSeqExtractor(BaseVCFSeqExtractor):

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        return self.variant_extractor.extract(
            interval, variants=self.vcf.fetch_variants(interval, sample_id),
            anchor=anchor, fixed_len=fixed_len)
