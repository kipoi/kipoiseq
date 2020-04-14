import pytest
from conftest import vcf_file, sample_5kb_fasta_file, test_with_multiple_variants
from kipoiseq.dataclasses import Variant, Interval
from kipoiseq.extractors.vcf_query import NumberVariantQuery
from kipoiseq.extractors.vcf import MultiSampleVCF

fasta_file = sample_5kb_fasta_file

intervals = [
    Interval('chr1', 3, 10),
    Interval('chr1', 4, 30),
    Interval('chr1', 19, 30)
]


@pytest.fixture
def multi_sample_vcf():
    return MultiSampleVCF(vcf_file)


def test_MultiSampleVCF__next__(multi_sample_vcf):
    variant = next(multi_sample_vcf)
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'


def test_MultiSampleVCF__iter__(multi_sample_vcf):
    variant = list(multi_sample_vcf)[0]
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'


def test_MultiSampleVCF_fetch_variant(multi_sample_vcf):
    interval = Interval('chr1', 3, 5)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 2
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 1
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00001'))) == 0

    interval = Interval('chr1', 7, 12)
    assert len(list(multi_sample_vcf.fetch_variants(interval))) == 0
    assert len(list(multi_sample_vcf.fetch_variants(interval, 'NA00003'))) == 0


def test_MultiSampleVCF_query_variants(multi_sample_vcf):
    vq = multi_sample_vcf.query_variants(intervals)
    variants = list(vq)
    
    assert len(variants) == 5
    assert variants[0].pos == 4
    assert variants[1].pos == 5
    
    
    msvcf = MultiSampleVCF(test_with_multiple_variants)
    vq = msvcf.query_variants([Interval('chr1', 3, 10)])
    variants = list(vq)
    
    assert len(variants) == 5
    assert variants[0].ref == 'T'
    assert variants[1].ref == 'T'
    assert variants[2].ref == 'T'
    assert variants[0].alt == 'C'
    assert variants[1].alt == 'A'
    assert variants[2].alt == 'G'
    assert variants[4].alt == ''
    
    vq = msvcf.query_variants([Interval('chr1', 11, 14)])
    variants = list(vq)
    
    assert len(variants) == 1
    assert variants[0].ref == 'T'
    assert variants[0].alt == ''


def test_MultiSampleVCF_get_samples(multi_sample_vcf):
    variants = list(multi_sample_vcf)
    samples = multi_sample_vcf.get_samples(variants[0])
    assert samples == {'NA00003': 3}


def test_MultiSampleVCF_get_variant(multi_sample_vcf):

    variant = multi_sample_vcf.get_variant("chr1:4:T>C")
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'

    variant = multi_sample_vcf.get_variant(Variant('chr1', 4, 'T', 'C'))
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'

    with pytest.raises(KeyError):
        multi_sample_vcf.get_variant("chr1:4:A>C")


def test_MultiSampleVCF_get_variants(multi_sample_vcf):
    variants = multi_sample_vcf.get_variants(["chr1:4:T>C"], intervals)
    assert len(variants) == 1

    variant = variants[0]
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'

    variants = multi_sample_vcf.get_variants(["chr1:4:T>C", "chr1:25:AACG>GA"])
    assert len(variants) == 2

    variant = variants[0]
    assert variant.chrom == 'chr1'
    assert variant.pos == 4
    assert variant.ref == 'T'
    assert variant.alt == 'C'


def test_MultiSampleVCF__regions_from_variants(multi_sample_vcf):
    variants = [
        Variant('chr1', 4, 'T', 'C'),
        Variant('chr1', 25, 'AACG', 'GA'),
        Variant('chr1', 55525, 'AACG', 'GA'),
        Variant('chr10', 55525, 'AACG', 'GA')

    ]
    regions = multi_sample_vcf._regions_from_variants(variants)

    assert set(regions) == set([
        Interval('chr1', 3, 25),
        Interval('chr1', 55524, 55525),
        Interval('chr10', 55524, 55525)
    ])


def test_MultiSampleVCF_VariantQueryable_to_vcf(tmpdir, multi_sample_vcf):
    output_vcf_file = str(tmpdir / 'output.vcf')

    multi_sample_vcf \
        .query_variants(intervals) \
        .filter_range(NumberVariantQuery(max_num=1)) \
        .to_vcf(output_vcf_file)

    vcf = MultiSampleVCF(output_vcf_file)
    variants = list(vcf)
    assert len(variants) == 1
    assert variants[0].ref == 'AACG'
    assert variants[0].alt == 'GA'


def test_batch_iter_vcf(multi_sample_vcf):
    batchs = list(multi_sample_vcf.batch_iter(10))
    assert sum(len(i) for i in batchs) == 3
