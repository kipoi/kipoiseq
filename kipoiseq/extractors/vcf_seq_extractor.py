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
        anchor_abs = anchor + interval.start
        istart = interval.start
        iend = interval.end

        variant_pairs = self._variant_to_sequence(variants)

        # 1. Split variants overlapping with anchor
        variant_pairs = list(self._split_overlapping(
            variant_pairs, anchor_abs))

        # 2. split the variants into upstream and downstream
        # and sort the variants in each interval
        upstream_variants = sorted(
            filter(lambda x: x[0].start >= anchor_abs, variant_pairs),
            key=lambda x: x[0].start)

        downstream_variants = sorted(
            filter(lambda x: x[0].start < anchor_abs, variant_pairs),
            key=lambda x: x[0].start, reverse=True)

        # 3. Iterate from the anchor point outwards. At each
        # register the interval from which to take the reference sequence
        # as well as the interval for the variant
        # Also, update the start and end of the interval
        up_sb = IntervalSeqBuilder()
        down_sb = IntervalSeqBuilder()

        prev = anchor_abs
        for ref, alt in upstream_variants:
            if ref.end > iend:
                if ref.start > iend:
                    break
                else:
                    ref, alt = self._cut_half((ref, alt), iend)[0]
            if fixed_len:
                iend -= len(alt) - len(ref)
            up_sb.append(Interval(interval.chrom, prev,
                                  ref.start, interval.strand))
            up_sb.append(alt)
            prev = ref.end
        up_sb.append(Interval(interval.chrom, prev, iend, interval.strand))

        prev = anchor_abs
        for ref, alt in downstream_variants:
            if ref.start < istart:
                if ref.end < istart:
                    break
                else:
                    ref, alt = self._cut_half((ref, alt), istart)[1]
            if fixed_len:
                istart -= len(alt) - len(ref)
            down_sb.append(Interval(interval.chrom, ref.end,
                                    prev, interval.strand))
            down_sb.append(alt)
            prev = ref.start
        down_sb.append(Interval(interval.chrom, istart, prev, interval.strand))
        down_sb.reverse()

        # 4. fetch the sequence and restore intervals in builder
        t_interval = Interval(interval.chrom, istart, iend, interval.strand)
        seq = self.fasta.extract(t_interval)
        seq = Sequence(name=interval.chrom, seq=seq, start=istart, end=iend)
        up_sb.restore(seq)
        down_sb.restore(seq)

        # 5. Concate sequences from the upstream and downstream splits. Concat
        # upstream and downstream sequence
        return down_sb.concat() + up_sb.concat()

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

    def _split_overlapping(self, variant_pairs, anchor_abs):
        """
        Split the variants hitting the anchor into two
        """
        for ref, alt in variant_pairs:
            if ref.start < anchor_abs < ref.end \
               or alt.start < anchor_abs < alt.end:
                fhalf, shalf = self._cut_half((ref, alt), anchor_abs)
                yield fhalf
                yield shalf
            else:
                yield ref, alt

    def _cut_half(self, variant_pair, cutting_point):
        """
        Cuts variant pairs based on cutting point.
        """
        ref, alt = variant_pair
        mid = cutting_point - ref.start
        return (
            (ref[:mid], alt[:mid]),
            (ref[mid:], alt[mid:])
        )


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
