import pytest
from conftest import vcf_file, gtf_file, example_intervals_bed
import pyranges
from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_matching import variants_to_pyranges, \
    pyranges_to_intervals, intervals_to_pyranges, BaseVariantMatcher, \
    SingleVariantMatcher, MultiVariantsMatcher

intervals = [
    Interval('chr1', 1, 10, strand='+'),
    Interval('chr1', 23, 30, strand='-')
]

variants = [
    Variant('chr1', 4, 'T', 'C'),
    Variant('chr1', 5, 'A', 'GA'),
    Variant('chr1', 25, 'AACG', 'GA')
]

pr = pyranges.PyRanges(
    chromosomes='chr1',
    starts=[1, 23, 5],
    ends=[10, 30, 50],
    strands=['+', '-', '.']
)


def test_variants_to_pyranges():
    vcf = MultiSampleVCF(vcf_file)
    variants = list(vcf)
    df = variants_to_pyranges(variants).df
    assert df.shape[0] == len(variants)

    v = df.iloc[0]
    assert v.Chromosome == 'chr1'
    assert v.Start == 3
    assert v.End == 4
    assert v.variant.ref == 'T'
    assert v.variant.alt == 'C'


def test_pyranges_to_intervals():
    pranges = pyranges.read_gtf(gtf_file)
    intervals = list(pyranges_to_intervals(pranges, interval_attrs=[
        'gene_id', 'gene_name', 'transcript_id', 'exon_id']))

    assert len(intervals) == 5
    assert intervals[4].attrs['gene_id'] == 'ENSG00000012048'
    assert intervals[4].attrs['gene_name'] == 'BRCA1'
    assert intervals[4].attrs['transcript_id'] == 'ENST00000357654'
    assert intervals[4].attrs['exon_id'] == 'ENSE00003510592'

    pranges = pyranges.read_bed(example_intervals_bed)
    intervals = list(pyranges_to_intervals(pranges))

    assert len(intervals) == 4
    assert intervals[0].start == 2

    assert pranges.Chromosome.tolist() == ['chr1'] * 4
    assert pranges.Start.tolist() == [2, 2, 2, 602]
    assert pranges.End.tolist() == [1000, 5000, 1002, 604]


def test_intervals_to_pyranges():
    pr = intervals_to_pyranges(intervals)

    assert pr.df.shape[0] == 2
    assert pr.df.Chromosome.tolist() == ['chr1', 'chr1']
    assert pr.df.Start.tolist() == [1, 23]
    assert pr.df.End.tolist() == [10, 30]
    assert pr.df.Strand.tolist() == ['+', '-']


def test_BaseVariantMatcher__read_intervals():
    pranges = pyranges.read_gtf(gtf_file)

    with pytest.raises(ValueError):
        pr = BaseVariantMatcher._read_intervals(
            pranges=pranges, gtf_path=gtf_file)

    with pytest.raises(ValueError):
        pr = BaseVariantMatcher._read_intervals(
            intervals=intervals, interval_attrs=['gene_id'])

    pr = BaseVariantMatcher._read_intervals(gtf_path=gtf_file)
    assert pr.Chromosome.tolist() == ['chr1'] * 5
    assert pr.Start.tolist() == [200, 200, 200, 1049, 3029]
    assert pr.End.tolist() == [4230, 4230, 402, 1340, 4230]
    # assert len(pr.intervals.tolist()) == 5

    pr = BaseVariantMatcher._read_intervals(bed_path=example_intervals_bed)
    assert pr.Chromosome.tolist() == ['chr1'] * 4
    assert pr.Start.tolist() == [2, 2, 2, 602]
    assert pr.End.tolist() == [1000, 5000, 1002, 604]
    # assert len(pr.intervals.tolist()) == 4

    pr = BaseVariantMatcher._read_intervals(pranges=pranges)
    assert pr.Chromosome.tolist() == ['chr1'] * 5
    assert pr.Start.tolist() == [200, 200, 200, 1049, 3029]
    assert pr.End.tolist() == [4230, 4230, 402, 1340, 4230]
    # assert len(pr.intervals.tolist()) == 5

    pr = BaseVariantMatcher._read_intervals(intervals=intervals)
    assert pr.df.Chromosome.tolist() == ['chr1', 'chr1']
    assert pr.df.Start.tolist() == [1, 23]
    assert pr.df.End.tolist() == [10, 30]
    assert pr.df.Strand.tolist() == ['+', '-']
    assert len(pr.intervals.tolist()) == 2


def test_SingleVariantMatcher__iter__():
    inters = intervals + [Interval('chr1', 5, 50)]

    matcher = SingleVariantMatcher(vcf_file, intervals=inters)
    pairs = list(matcher)

    assert (inters[0], variants[0]) in pairs
    assert (inters[0], variants[1]) in pairs
    assert (inters[1], variants[2]) in pairs
    assert (inters[2], variants[2]) in pairs

    assert len(pairs) == 4

    matcher = SingleVariantMatcher(vcf_file, pranges=pr)
    pairs = list(matcher)

    assert (inters[0], variants[0]) in pairs
    assert (inters[0], variants[1]) in pairs
    assert (inters[1], variants[2]) in pairs
    assert (inters[2], variants[2]) in pairs
    assert len(pairs) == 4


def test_MultiVariantMatcher__iter__():
    matcher = MultiVariantsMatcher(vcf_file, intervals=intervals)
    pairs = list(matcher)

    assert pairs[0][0] == intervals[0]
    assert list(pairs[0][1]) == [variants[0], variants[1]]
    assert pairs[1][0] == intervals[1]
    assert list(pairs[1][1]) == [variants[2]]

    matcher = MultiVariantsMatcher(vcf_file, pranges=pr)
    pairs = list(matcher)

    assert pairs[0][0] == intervals[0]
    assert list(pairs[0][1]) == [variants[0], variants[1]]
    assert pairs[1][0] == intervals[1]
    assert list(pairs[1][1]) == [variants[2]]
