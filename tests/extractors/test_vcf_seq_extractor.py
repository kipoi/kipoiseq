import pytest
from cyvcf2 import VCF
from pyfaidx import Sequence
from pybedtools import Interval
from kipoiseq.extractors.vcf_seq import IntervalSeqBuilder
from kipoiseq.extractors import *

fasta_file = 'tests/data/sample.5kb.fa'
vcf_file = 'tests/data/test.vcf.gz'


@pytest.fixture
def multi_sample_vcf():
    return MultiSampleVCF(vcf_file)


def test_multi_sample_vcf_fetch_variant(multi_sample_vcf):
    interval = Interval('chr1', 3, 5)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 2
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 1

    interval = Interval('chr1', 7, 12)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 0
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 0


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

    pair = (Sequence(seq='T', start=3, end=4),
            Sequence(seq='AAA', start=3, end=6))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 5))

    assert splited_pairs[0][0].seq == 'T'
    assert splited_pairs[0][1].seq == 'AA'
    assert splited_pairs[1][0].seq == ''
    assert splited_pairs[1][1].seq == 'A'


def test_extract(variant_seq_extractor):
    interval = Interval('chr1', 2, 9)
    variants = list(VCF(vcf_file)())
    seq = variant_seq_extractor.extract(interval, variants, anchor=5)
    assert len(seq) == interval.end - interval.start
    assert seq == 'GCGAACG'

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
