import pytest
from cyvcf2 import VCF
from pyfaidx import Sequence
from pybedtools import Interval
from kipoiseq.extractors.vcf_seq_extractor import IntervalSeqBuilder
from kipoiseq.extractors import *

fasta_file = 'tests/data/sample.fasta'
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


def test_interval_seq_builder_concat(interval_seq_builder):
    sequence = Sequence(seq='CCCCATCGNN', start=10, end=20)
    interval_seq_builder.restore(sequence)
    assert interval_seq_builder.concat() == 'CCCCTAGCNN'


@pytest.fixture
def variant_seq_extractor():
    return VariantSeqExtractor(fasta_file)


# def test__split_overlapping(variant_seq_extractor):
#     interval = Interval('chr1', 2, 9)
#     variants = list(VCF(vcf_file)())
#     variant_seq_extractor.extract(interval, variants, anchor=5)

# def test_extract(variant_seq_extractor):
#     pass
