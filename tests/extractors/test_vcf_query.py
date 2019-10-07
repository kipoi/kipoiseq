import pytest
from conftest import vcf_file
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_seq import MultiSampleVCF
from kipoiseq.extractors.vcf_query import *


@pytest.fixture
def query_true():
    return VariantQuery(lambda v: True)


@pytest.fixture
def query_false():
    return VariantQuery(lambda v: False)


def test_base_query__and__(query_false, query_true):
    assert not (query_false & query_true)(None)


def test_base_query__or__(query_false, query_true):
    assert (query_false | query_true)(None)


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
