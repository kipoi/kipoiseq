import abc
from typing import Union

from pyfaidx import Sequence, complement
from cyvcf2 import VCF
from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.extractors import (
    BaseExtractor,
    FastaStringExtractor,
    MultiSampleVCF,
)

from kipoiseq import __version__
from deprecation import deprecated
import math

__all__ = [
    'VariantSeqExtractor',
    'SingleVariantVCFSeqExtractor',
    'SingleSeqVCFSeqExtractor',
    'SampleSeqExtractor'
]


class IntervalSeqBuilder(list):
    """
    String builder for `pyfaidx.Sequence` and `Interval` objects.
    """

    def restore(self, sequence: Sequence):
        """
        Args:
          sequence: `pyfaidx.Sequence` which convert all interval inside
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

    def __init__(self, fasta_file: str = None, reference_sequence: BaseExtractor = None, use_strand=True):
        """
        Sequence extractor which allows to obtain the alternative sequence,
        given some interval and variants inside this interval.

        Args:
            fasta_file: path to the fasta file (can be gzipped)
            reference_sequence: extractor returning the reference sequence given some interval
            use_strand (bool): if True, the extracted sequence
                is reverse complemented in case interval.strand == "-"
        """
        self._use_strand = use_strand

        if fasta_file is not None:
            if reference_sequence is not None:
                raise ValueError(
                    "either fasta_file or ref_seq_extractor have to be specified")
            self._ref_seq_extractor = FastaStringExtractor(
                fasta_file, use_strand=False)
        else:
            if reference_sequence is None:
                raise ValueError(
                    "either fasta_file or ref_seq_extractor have to be specified")
            self._ref_seq_extractor = reference_sequence

    @property
    @deprecated(deprecated_in="1.0",
                # removed_in="2.0",
                current_version=__version__,
                details="Use `ref_seq_extractor` instead")
    def fasta(self):
        return self._ref_seq_extractor

    @property
    def ref_seq_extractor(self) -> BaseExtractor:
        """

        Returns:
            The reference sequence extractor of this object
        """
        return self._ref_seq_extractor

    def extract(self, interval, variants, anchor, fixed_len=True, use_strand=None, chrom_len=math.inf,
                is_padding=False, **kwargs):
        """
        Args:
            interval: pybedtools.Interval Region of interest from
                which to query the sequence. 0-based
            variants: List[cyvcf2.Variant]: variants overlapping the `interval`.
                can also be indels. 1-based
            anchor: absolution position w.r.t. the interval start. (0-based).
                E.g. for an interval of `chr1:10-20` the anchor of 10 denotes
                the point chr1:10 in the 0-based coordinate system.
            fixed_len: if True, the return sequence will have the same length
                as the `interval` (e.g. `interval.end - interval.start`)
            use_strand (bool, optional): if True, the extracted sequence
                is reverse complemented in case interval.strand == "-".
                Overrides `self.use_strand`
            chrom_len: length of the chromosome. If chrom_len == math.inf, the length of the chromosome is not checked.
            is_padding: if True, the sequence is padded with 'N's if sequence can't extend to the fixed length,

        Returns:
            A single sequence (`str`) with all the variants applied.
        """
        # Preprocessing
        anchor = max(min(anchor, interval.end), interval.start)
        variant_pairs = self._variant_to_sequence(variants)

        # 1. Split variants overlapping with anchor
        # and interval start end if not fixed_len
        variant_pairs = self._split_overlapping(variant_pairs, anchor)

        if not fixed_len:
            variant_pairs = self._split_overlapping(
                variant_pairs, interval.start, which='right')
            variant_pairs = self._split_overlapping(
                variant_pairs, interval.end, which='left')

        variant_pairs = list(variant_pairs)

        # 2. split the variants into upstream and downstream
        # and sort the variants in each interval
        upstream_variants = sorted(
            filter(lambda x: x[0].start >= anchor, variant_pairs),
            key=lambda x: x[0].start
        )

        downstream_variants = sorted(
            filter(lambda x: x[0].start < anchor, variant_pairs),
            key=lambda x: x[0].start,
            reverse=True
        )

        # 3. Extend start and end position for deletions
        if fixed_len:
            istart, iend = self._updated_interval(
                interval, upstream_variants, downstream_variants)
        else:
            istart, iend = interval.start, interval.end

        istart = max(istart, 0)
        iend = min(iend, chrom_len - 1)

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
                down_str, up_str, interval, anchor, is_padding=is_padding)

        seq = down_str + up_str

        if use_strand is None:
            use_strand = self.use_strand
        if use_strand and interval.strand == '-':
            # reverse-complement
            seq = complement(seq)[::-1]

        return seq

    @staticmethod
    def _variant_to_sequence(variants):
        """
        Convert `cyvcf2.Variant` objects to `pyfaidx.Seqeunce` objects
        for reference and variants.
        """
        for v in variants:
            ref = Sequence(name=v.chrom, seq=v.ref,
                           start=v.start, end=v.start + len(v.ref))
            alt = Sequence(name=v.chrom, seq=v.alt,
                           start=v.start, end=v.start + len(v.alt))
            yield ref, alt

    @staticmethod
    def _split_overlapping(variant_pairs, anchor, which='both'):
        """
        Split the variants hitting the anchor into two
        """
        for ref, alt in variant_pairs:
            if ref.start < anchor < ref.end:
                mid = anchor - ref.start
                if which == 'left' or which == 'both':
                    yield ref[:mid], alt[:mid]
                if which == 'right' or which == 'both':
                    yield ref[mid:], alt[mid:]
            else:
                yield ref, alt

    @staticmethod
    def _updated_interval(interval, up_variants, down_variants):
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

    @staticmethod
    def _downstream_builder(down_variants, interval, anchor, istart):
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

    @staticmethod
    def _upstream_builder(up_variants, interval, anchor, iend):
        up_sb = IntervalSeqBuilder()

        prev = anchor
        for ref, alt in up_variants:
            if ref.start >= iend:
                break
            up_sb.append(Interval(interval.chrom, prev, ref.start))
            up_sb.append(alt)
            prev = ref.end
        up_sb.append(Interval(interval.chrom, prev, iend))

        return up_sb

    def _fetch(self, interval, istart, iend):
        # fetch interval, ignore strand
        seq = self.ref_seq_extractor.extract(
            Interval(interval.chrom, istart, iend))
        seq = Sequence(name=interval.chrom, seq=seq, start=istart, end=iend)
        return seq

    @staticmethod
    def _cut_to_fix_len(down_str, up_str, interval, anchor, is_padding=False):
        down_len = anchor - interval.start
        down_diff = len(down_str) - down_len
        if down_diff > 0:
            down_str = down_str[-down_len:]
        elif down_diff < 0:
            if is_padding:
                down_str = 'N' * abs(down_diff) + down_str
            else:
                raise ValueError(f"padding should be set to True, if the sequence can't extend to the fixed length")

        up_len = interval.end - anchor
        up_diff = len(up_str) - up_len
        if up_diff > 0:
            up_str = up_str[: up_len]
        elif up_diff < 0:
            if is_padding:
                up_str = up_str + 'N' * abs(up_diff)
            else:
                raise ValueError(f"padding should be set to True, if the sequence can't extend to the fixed length")

        return down_str, up_str


class _BaseVCFSeqExtractor(BaseExtractor, metaclass=abc.ABCMeta):
    """
    Base class to fetch sequence in which variants applied based
    on given vcf file.
    """

    def __init__(self, fasta_file, vcf_file):
        """
        Args:
          fasta_file: path to the fasta file (can be gzipped)
          vcf_file: path to the fasta file (need be bgzipped and indexed)
        """
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.variant_extractor = VariantSeqExtractor(fasta_file)
        self.vcf = MultiSampleVCF(vcf_file)

    @abc.abstractmethod
    def extract(self, interval: Interval, *args, **kwargs) -> str:
        raise NotImplementedError()


class SingleVariantVCFSeqExtractor(_BaseVCFSeqExtractor):
    """
    Fetch list of sequence in which each variant applied based
    on given vcf file.
    """

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        for variant in self.vcf.fetch_variants(interval, sample_id):
            yield self.variant_extractor.extract(
                interval,
                variants=[variant],
                anchor=anchor,
                fixed_len=fixed_len
            )


class SingleSeqVCFSeqExtractor(_BaseVCFSeqExtractor):
    """
    Fetch sequence in which all variant applied based on given vcf file.
    """

    def extract(self, interval, anchor=None, sample_id=None, fixed_len=True):
        return self.variant_extractor.extract(
            interval,
            variants=self.vcf.fetch_variants(interval, sample_id),
            anchor=anchor,
            fixed_len=fixed_len
        )


class SampleSeqExtractor(VariantSeqExtractor):
    def __init__(self, fasta_file, vcf_file):
        """Sequence extractor which can extract an alternate sequence for a
        given interval and the variants corresponding to a given
        sample and phase.

        Args:
          fasta_file: Path to the fasta file containing the reference
            sequence (can be gzipped)
          vcf_file: Path to the VCF file containing phased genotype information

        """
        self.vcf = VCF(vcf_file)
        self._sample_indices = dict(zip(self.vcf.samples,
                                        range(len(self.vcf.samples))))

        super().__init__(fasta_file)

    def extract(self, interval, sample, phase, anchor,
                fixed_len=True, **kwargs):
        """Extracts an alternate sequence for a given interval and the
        variants corresponding to a given sample.

        Args:
          interval: `kipoiseq.dataclasses.Interval`, Region of
            interest from which to query the sequence. 0-based.
          sample: `str`, Sample from the VCF file for which variants should be
            extracted.
          phase: `0` or `1`, Phase for which sequence should be extracted
          anchor: `int`, Absolution position w.r.t. the interval
            start. (0-based).  E.g. for an interval of `chr1:10-20`
            the anchor of 10 denotes the point chr1:10 in the 0-based
            coordinate system.
          fixed_len: `bool`, If True, the return sequence will have the
            same length as the `interval` (e.g. `interval.end -
            interval.start`)
          kwargs: Additional keyword arguments to pass to
            `SampleSeqExtractor.extract`

        Returns:
          A single sequence (`str`) with all the variants applied.
        """
        variants = []
        if sample is not None:
            if sample not in self.vcf.samples:
                raise ValueError(f'Sample \'{sample}\' '
                                 'not present in VCF file')

            if phase not in (0, 1):
                raise ValueError('phase argument must be in (0, 1)'
                                 ' if sample is not None')

            # Interval is  0-based, cyvcf2 positions are 1-based: need to add 1
            variants = self._get_sample_variants(
                self.vcf(
                    f'{interval.chrom}:'
                    + f'{interval.start + 1}-{interval.end + 1}'
                ),
                sample,
                phase
            )

        return super(SampleSeqExtractor, self).extract(
            interval, variants, anchor, fixed_len, **kwargs)

    def _get_sample_variants(self, variants, sample, phase):
        """Given a list of `cyvcf2.Variant`, returns all those present for a
        given sample and phase and converts them to
        `kipoiseq.dataclasses.Variant`

        Args:
          variants: List of `cyvcf2.Variant`, Variants of interest
          sample: `str`, Sample for which to filter genotypes
          phase: `0` or `1`, Phase for which to filter genotypes

        Returns:
          List of `kipoiseq.dataclasses.Variant`
        """
        sample_index = self._sample_indices[sample]
        return [
            Variant.from_cyvcf(v) for v in variants
            if v.genotypes[sample_index][phase]
        ]
