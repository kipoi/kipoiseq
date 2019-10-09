import pytest
from conftest import vcf_file, sample_5kb_fasta_file
from cyvcf2 import VCF
from pyfaidx import Sequence
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_seq import IntervalSeqBuilder, \
    VariantSeqExtractor, SingleSeqVCFSeqExtractor, SingleVariantVCFSeqExtractor

fasta_file = sample_5kb_fasta_file


@pytest.fixture
def interval_seq_builder():
    return IntervalSeqBuilder([
        Interval('chr1', 10, 13),
        Interval('chr1', 13, 14),
        Sequence(seq='TAGC', start=14, end=18),
        Interval('chr1', 18, 20)
    ])


def test_interval_seq_builder_restore(interval_seq_builder):
    sequence = Sequence(seq='CCCCATCGTT', start=10, end=20)
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[0].seq == 'CCC'
    assert interval_seq_builder[1].seq == 'C'
    assert interval_seq_builder[2].seq == 'TAGC'
    assert interval_seq_builder[3].seq == 'TT'

    interval_seq_builder.append(Interval('chr1', 5, 10))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[4].seq == ''

    interval_seq_builder.append(Interval('chr1', 20, 25))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[5].seq == ''

    interval_seq_builder.append(Interval('chr1', 10, 5))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[6].seq == ''

    interval_seq_builder.append(Interval('chr1', 25, 20))
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder[7].seq == ''


def test_interval_seq_builder_concat(interval_seq_builder):
    with pytest.raises(TypeError):
        interval_seq_builder.concat()

    sequence = Sequence(seq='CCCCATCGNN', start=10, end=20)
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder.concat() == 'CCCCTAGCNN'


@pytest.fixture
def variant_seq_extractor():
    return VariantSeqExtractor(fasta_file)


def test__split_overlapping(variant_seq_extractor):
    pair = (Sequence(seq='AAA', start=3, end=6),
            Sequence(seq='T', start=3, end=4))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 5))

    assert splited_pairs[0][0].seq == 'AA'
    assert splited_pairs[0][1].seq == 'T'
    assert splited_pairs[1][0].seq == 'A'
    assert splited_pairs[1][1].seq == ''

    pair = (Sequence(seq='TT', start=3, end=5),
            Sequence(seq='AAA', start=3, end=6))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 4))

    assert splited_pairs[0][0].seq == 'T'
    assert splited_pairs[0][1].seq == 'A'
    assert splited_pairs[1][0].seq == 'T'
    assert splited_pairs[1][1].seq == 'AA'


def test_extract(variant_seq_extractor):
    variants = [Variant.from_cyvcf(v) for v in VCF(vcf_file)]

    interval = Interval('chr1', 2, 9)

    seq = variant_seq_extractor.extract(interval, variants, anchor=5)
    assert len(seq) == interval.end - interval.start
    assert seq == 'CGAACGT'

    interval = Interval('chr1', 2, 9, strand='-')
    seq = variant_seq_extractor.extract(interval, variants, anchor=5)
    assert len(seq) == interval.end - interval.start
    assert seq == 'ACGTTCG'

    interval = Interval('chr1', 4, 14)
    seq = variant_seq_extractor.extract(interval, variants, anchor=7)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 4, 14)
    seq = variant_seq_extractor.extract(interval, variants, anchor=4)
    assert len(seq) == interval.end - interval.start
    assert seq == 'GAACGTAACG'

    interval = Interval('chr1', 2, 5)
    seq = variant_seq_extractor.extract(interval, variants, anchor=3)
    assert len(seq) == interval.end - interval.start
    assert seq == 'GCG'

    interval = Interval('chr1', 24, 34)
    seq = variant_seq_extractor.extract(interval, variants, anchor=27)
    assert len(seq) == interval.end - interval.start
    assert seq == 'TGATAACGTA'

    interval = Interval('chr1', 25, 35)
    seq = variant_seq_extractor.extract(interval, variants, anchor=34)
    assert len(seq) == interval.end - interval.start
    assert seq == 'TGATAACGTA'

    interval = Interval('chr1', 34, 44)
    seq = variant_seq_extractor.extract(interval, variants, anchor=37)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 34, 44)
    seq = variant_seq_extractor.extract(interval, variants, anchor=100)
    assert len(seq) == interval.end - interval.start
    assert seq == 'AACGTAACGT'

    interval = Interval('chr1', 5, 11, strand='+')
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=False)
    assert seq == 'ACGTAA'

    interval = Interval('chr1', 0, 3, strand='+')
    seq = variant_seq_extractor.extract(
        interval, variants, anchor=10, fixed_len=False)
    assert seq == 'ACG'


@pytest.fixture
def single_variant_vcf_seq_extractor():
    return SingleVariantVCFSeqExtractor(fasta_file, vcf_file)


def test_single_variant_vcf_seq_extract(single_variant_vcf_seq_extractor):
    interval = Interval('chr1', 2, 9)
    seqs = single_variant_vcf_seq_extractor.extract(interval, anchor=3)
    assert next(seqs) == 'GCAACGT'
    assert next(seqs) == 'GTGAACG'


@pytest.fixture
def single_seq_vcf_seq_extractor():
    return SingleSeqVCFSeqExtractor(fasta_file, vcf_file)


def test_single_seq_vcf_seq_extract(single_seq_vcf_seq_extractor):
    interval = Interval('chr1', 2, 9)
    seq = single_seq_vcf_seq_extractor.extract(interval, anchor=3)
    assert seq == 'GCGAACG'
