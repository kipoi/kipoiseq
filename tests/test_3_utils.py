from conftest import vcf_file, sample_5kb_fasta_file, example_intervals_bed
import pytest
import numpy as np
import pyranges
from kipoiseq.extractors import FastaStringExtractor, MultiSampleVCF
from kipoiseq.utils import parse_alphabet, parse_dtype, \
    compare_chrom_annotation


def test_parse_alphabet():
    assert parse_alphabet(['A', 'C']) == ['A', 'C']
    assert parse_alphabet('AC') == ['A', 'C']


def test_parse_type():
    with pytest.raises(Exception):
        parse_dtype('string')
    with pytest.raises(Exception):
        parse_dtype('np.float322')

    assert parse_dtype('float') == float
    assert parse_dtype(float) == float
    assert parse_dtype("np.float32") == np.float32


def test_compare_chrom_annotation():
    sources = [
        MultiSampleVCF(vcf_file),
        FastaStringExtractor(sample_5kb_fasta_file),
        pyranges.read_bed(example_intervals_bed)
    ]

    with pytest.raises(ValueError):
        assert compare_chrom_annotation([])

    with pytest.raises(ValueError):
        assert compare_chrom_annotation([object()])

    assert compare_chrom_annotation(sources) == {'chr1'}
    assert compare_chrom_annotation(sources, strategy='all') == {'chr1'}

    with pytest.raises(AssertionError) as exception:
        sources[1].fasta = {'chr1': '', 'chr2': '', 'chr3': ''}
        compare_chrom_annotation(sources, strategy='all')

    assert str(exception.value) == 'chroms annotations are not all same.'

    assert compare_chrom_annotation(sources) == {'chr1'}

    with pytest.raises(AssertionError) as exception:
        sources[1].fasta = {'chr2': '', 'chr3': ''}
        compare_chrom_annotation(sources)
    assert str(exception.value) == 'there is not intersection between chromosomes.'
