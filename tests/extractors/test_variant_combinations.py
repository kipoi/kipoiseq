import pytest
from conftest import example_intervals_bed, sample_5kb_fasta_file
import pyranges as pr
from kipoiseq import Interval
from kipoiseq.extractors import VariantCombinator


@pytest.fixture
def variant_combinator():
    return VariantCombinator(sample_5kb_fasta_file, example_intervals_bed)


def test_VariantCombinator_combination_variants(variant_combinator):
    interval = Interval('chr1', 20, 30)
    variants = list(variant_combinator.combination_variants(interval, 'snv'))
    assert len(variants) == 30

    interval = Interval('chr1', 20, 22)
    variants = list(variant_combinator.combination_variants(interval, 'snv'))
    assert variants[0].chrom == 'chr1'
    assert variants[0].ref == 'A'
    assert variants[0].alt == 'C'
    assert variants[1].alt == 'G'
    assert variants[2].alt == 'T'

    assert variants[3].ref == 'C'
    assert variants[3].alt == 'A'
    assert variants[4].alt == 'G'
    assert variants[5].alt == 'T'

    interval = Interval('chr1', 20, 22)
    variants = list(variant_combinator.combination_variants(interval, 'in'))
    len(variants) == 32
    assert variants[0].ref == 'A'
    assert variants[0].alt == 'AA'
    assert variants[15].alt == 'TT'

    assert variants[16].ref == 'C'
    assert variants[16].alt == 'AA'
    assert variants[31].alt == 'TT'

    interval = Interval('chr1', 20, 22)
    variants = list(variant_combinator.combination_variants(
        interval, 'del', del_length=2))
    assert len(variants) == 3
    assert variants[0].ref == 'A'
    assert variants[0].alt == ''
    assert variants[1].ref == 'AC'
    assert variants[1].alt == ''
    assert variants[2].ref == 'C'
    assert variants[2].alt == ''

    variants = list(variant_combinator.combination_variants(
        interval, 'all', in_length=2, del_length=2))
    assert len(variants) == 6 + 32 + 3


def test_VariantCombinator_iter(variant_combinator):
    variants = list(variant_combinator)
    df = pr.read_bed(example_intervals_bed).df
    num_snv = (df['End'] - df['Start']).sum() * 3
    assert len(variants) == num_snv
