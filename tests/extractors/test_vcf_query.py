import pytest
from conftest import vcf_file
import pandas as pd
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from kipoiseq.extractors.vcf_query import *


@pytest.fixture
def query_true():
    return VariantQuery(lambda v: True)


@pytest.fixture
def query_false():
    return VariantQuery(lambda v: False)


def test_VariantQuery__and__(query_false, query_true):
    assert not (query_false & query_true)(None)


def test_VariantQuery__or__(query_false, query_true):
    assert (query_false | query_true)(None)


@pytest.fixture
def query_interval_true():
    return VariantIntervalQuery(lambda vs, i: [True, True, False, False])


@pytest.fixture
def query_interval_false():
    return VariantIntervalQuery(lambda vs, i: [True, False, True, False])


def test_VariantIntervalQuery__and__(query_interval_false,
                                     query_interval_true):
    assert (query_interval_false & query_interval_true)(None, None) == [
        True, False, False, False]


def test_VariantIntervalQuery__or__(query_interval_false,
                                    query_interval_true):
    assert (query_interval_false | query_interval_true)(None, None) == [
        True, True, True, False]


@pytest.fixture
def variant_queryable():
    vcf = MultiSampleVCF(vcf_file)
    return VariantIntervalQueryable(vcf, [
        (
            [
                Variant('chr1', 12, 'A', 'T'),
                Variant('chr1', 18, 'A', 'C', filter='q10'),
            ],
            Interval('chr1', 10, 20)
        ),
        (
            [
                Variant('chr2', 120, 'AT', 'AAAT'),
            ],
            Interval('chr2', 110, 200)
        )
    ])


def test_variant_queryable__iter__(variant_queryable):
    variants = list(variant_queryable)
    assert len(variants) == 3
    assert variants[0].ref == 'A'
    assert variants[0].alt == 'T'


def test_variant_queryable_filter_1(variant_queryable):
    assert len(list(variant_queryable.filter(lambda v: v.alt == 'T'))) == 1


def test_variant_queryable_filter_2(variant_queryable):
    assert len(list(variant_queryable.filter(lambda v: v.ref == 'A'))) == 2


def test_variant_filter_range(variant_queryable):
    assert 2 == len(list(variant_queryable.filter_range(
        lambda variants, interval: (v.ref == 'A' for v in variants))))


def test_VariantQueryable_filter_by_num_max(variant_queryable):
    assert 1 == len(list(variant_queryable.filter_range(
        NumberVariantQuery(max_num=1))))


def test_VariantQueryable_filter_by_num_min(variant_queryable):
    assert 2 == len(list(variant_queryable.filter_range(
        NumberVariantQuery(min_num=2))))


def test_VariantQueryable_filter_variant_query_2(variant_queryable):
    assert 2 == len(list(variant_queryable.filter(FilterVariantQuery())))


def test_VariantQueryable_filter_variant_query_3(variant_queryable):
    assert 3 == len(list(variant_queryable.filter(
        FilterVariantQuery() | FilterVariantQuery(filter='q10'))))


def test_VariantQueryable_batch_iter():
    vcf = MultiSampleVCF(vcf_file)

    variant_queryable = vcf.query_all()
    batches = list(variant_queryable.batch_iter(batch_size=1))
    assert len(batches) == 3

    vcf = MultiSampleVCF(vcf_file)
    variant_queryable = vcf.query_all()
    batches = list(variant_queryable.batch_iter(batch_size=2))

    assert len(batches) == 2
    assert len(batches[0].variant_intervals[0][0]) == 2
    assert len(batches[1].variant_intervals[0][0]) == 1

    vcf = MultiSampleVCF(vcf_file)
    variant_queryable = vcf.query_all()
    batches = list(variant_queryable.batch_iter(batch_size=10))
    assert len(batches) == 1

    variants, interval = batches[0].variant_intervals[0]
    assert interval == Interval('chr1', 3, 25)
    assert len(variants) == 3


def test_VariantQueryable_to_vcf(tmp_path):
    path = str(tmp_path / 'a.vcf')
    vcf = MultiSampleVCF(vcf_file)
    variant_queryable = vcf.query_all()
    variant_queryable.to_vcf(path, remove_samples=True, clean_info=True)

    vcf = MultiSampleVCF(path)
    assert len(vcf.samples) == 0


def test_VariantQueryable_to_sample_csv(tmp_path):
    vcf = MultiSampleVCF(vcf_file)

    variant_queryable = vcf.query_all()

    path = str(tmp_path / 'sample.csv')
    variant_queryable.to_sample_csv(path)

    df = pd.read_csv(path)
    df_expected = pd.DataFrame({
        'variant': ['chr1:4:T>C', 'chr1:25:AACG>GA'],
        'sample': ['NA00003', 'NA00002'],
        'genotype': [3, 3],
        'GT': ['1/1', '1/1'],
        'HQ': ['51,51', '10,10']
    })
    pd.testing.assert_frame_equal(df, df_expected)
