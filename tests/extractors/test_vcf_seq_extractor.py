import pytest
from cyvcf2 import VCF
from pyfaidx import Sequence
from pybedtools import Interval
from kipoiseq.extractors.vcf_seq import IntervalSeqBuilder, VariantQueryable
from kipoiseq.extractors import *

fasta_file = 'tests/data/sample.5kb.fa'
vcf_file = 'tests/data/test.vcf.gz'


intervals = [
    Interval('chr1', 4, 10),
    Interval('chr1', 5, 30),
    Interval('chr1', 20, 30)
]


@pytest.fixture
def multi_sample_vcf():
    return MultiSampleVCF(vcf_file)


def test_multi_sample_vcf_fetch_variant(multi_sample_vcf):
    interval = Interval('chr1', 3, 5)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 2
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 1
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00001'))) == 0

    interval = Interval('chr1', 7, 12)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 0
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 0


def test_multi_sample_vcf_fetch_samples_with_variants(multi_sample_vcf):
    intervals = [Interval('chr1', 3, 10)]
    d = multi_sample_vcf.fetch_samples_with_variants(intervals)
    assert len(d) == 1
    assert len(d['NA00003']) == 1

    intervals = [Interval('chr1', 3, 10), Interval('chr1', 4, 7)]
    d = multi_sample_vcf.fetch_samples_with_variants(intervals)
    assert len(d) == 1
    assert len(d['NA00003']) == 1


def test_multi_sample_query_samples(multi_sample_vcf):
    intervals = [Interval('chr1', 3, 10)]
    d = list(multi_sample_vcf.query_samples(intervals))
    assert len(d[0]) == 1
    assert len(d[0]['NA00003']) == 1


def test_query_variants(multi_sample_vcf):
    vq = multi_sample_vcf.query_variants(intervals)
    variants = list(vq)
    assert len(variants) == 5
    assert variants[0].end == 4
    assert variants[1].end == 5


def test_get_samples(multi_sample_vcf):
    variants = list(multi_sample_vcf)
    samples = multi_sample_vcf.get_samples(variants[0])
    assert samples == {'NA00003': 3}


def test_get_variant_by_id(multi_sample_vcf):
    variant = multi_sample_vcf.get_variant_by_id("chr1:4:T:['C']")
    assert variant.CHROM == 'chr1'
    assert variant.POS == 4
    assert variant.REF == 'T'
    assert variant.ALT[0] == 'C'


@pytest.fixture
def variant_queryable(multi_sample_vcf):
    variants = [(multi_sample_vcf.fetch_variants(i), i) for i in intervals]
    return VariantQueryable(multi_sample_vcf, variants)


def test_VariantQueryable__iter__(variant_queryable):
    variants = list(variant_queryable)
    assert len(variants) == 5
    assert variants[0].REF == 'T'
    assert variants[0].ALT[0] == 'C'


def test_VariantQueryable_filter_all(variant_queryable):
    assert 2 == len(list(variant_queryable.filter_all(
        lambda variants, interval: (v.REF == 'A' for v in variants))))


def test_VariantQueryable_filter_by_num_max(variant_queryable):
    assert 1 == len(list(variant_queryable.filter_by_num(max_num=1)))


def test_VariantQueryable_filter_by_num_min(variant_queryable):
    assert 4 == len(list(variant_queryable.filter_by_num(min_num=2)))


def test_VariantQueryable_to_vcf(tmpdir, variant_queryable):
    output_vcf_file = str(tmpdir / 'output.vcf')

    variant_queryable \
        .filter_by_num(max_num=1) \
        .to_vcf(output_vcf_file)

    vcf = MultiSampleVCF(output_vcf_file)
    variants = list(vcf)
    assert len(variants) == 1
    assert variants[0].REF == 'AACG'
    assert variants[0].ALT[0] == 'GA'


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

    pair = (Sequence(seq='TT', start=3, end=5),
            Sequence(seq='AAA', start=3, end=6))
    splited_pairs = list(variant_seq_extractor._split_overlapping([pair], 4))

    assert splited_pairs[0][0].seq == 'T'
    assert splited_pairs[0][1].seq == 'A'
    assert splited_pairs[1][0].seq == 'T'
    assert splited_pairs[1][1].seq == 'AA'


def test_extract(variant_seq_extractor):
    variants = list(VCF(vcf_file)())

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
