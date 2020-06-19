import pytest
from kipoiseq.dataloaders import SingleVariantProteinDataLoader
from kipoiseq.dataloaders import SingleVariantUTRDataLoader

gtf_file = 'tests/data/sample_1_protein.gtf'
fasta_file = 'tests/data/demo_dna_seq.fa'
vcf_file = 'tests/data/singleVar_vcf_enst_test2.vcf.gz'


@pytest.fixture
def single_variant_protein_dataLoader():
    return SingleVariantProteinDataLoader(gtf_file, fasta_file, vcf_file)


def test_single_variant_protein_dataLoader(single_variant_protein_dataLoader):
    units = list(single_variant_protein_dataLoader)

    assert len(units) == 3

    assert units[2]['metadata']['transcript_id'] == 'enst_test2'

    txt_file = 'tests/data/Output_singleSeq_vcf_enst_test2.txt'
    expected_seq = open(txt_file).readline()
    assert units[0]['inputs']['ref_seq'][1:] == expected_seq[1:]

    txt_file = 'tests/data/Output_singleVar_vcf_enst_test2.txt'
    expected_seq = open(txt_file).read().splitlines()
    assert units[0]['inputs']['alt_seq'] == expected_seq[0]
    assert units[1]['inputs']['alt_seq'] == expected_seq[1]
    assert units[2]['inputs']['alt_seq'] == expected_seq[2]

    assert type(units[2]['metadata']) == dict
    assert len(units[2]) == 2
    assert len(units[2]['inputs']) == 2
    assert len(units[2]['metadata']) == 2  # number of columns

    assert units[0]['metadata']['variant']['ref'] == "C"
    assert units[0]['metadata']['variant']['alt'] == "G"
    assert units[0]['metadata']['variant']['pos'] == 596
    assert units[0]['metadata']['variant']['id'] == "9"

    assert units[1]['metadata']['variant']['ref'] == "A"
    assert units[1]['metadata']['variant']['alt'] == "G"
    assert units[1]['metadata']['variant']['pos'] == 597
    assert units[1]['metadata']['variant']['id'] == "10"

    assert units[2]['metadata']['variant']['ref'] == "T"
    assert units[2]['metadata']['variant']['alt'] == "G"
    assert units[2]['metadata']['variant']['pos'] == 598
    assert units[2]['metadata']['variant']['id'] == "11"


chr22_fasta_file = 'tests/data/chr22.fa.gz'
chr22_gtf_file = 'tests/data/chr22_ENST00000319363.gtf'
chr22_vcf_file = 'tests/data/chr22_ENST00000319363.vcf.gz'
chr22_5UTR_vcf_file = 'tests/data/chr22_ENST00000319363_5UTR.vcf.gz'
chr22_5UTR_ref_seq = 'tests/data/chr22_ENST00000319363_5UTR.ref_seq.txt'
chr22_5UTR_alt_seqs = 'tests/data/chr22_ENST00000319363_5UTR.alt_seqs.txt'
chr22_3UTR_vcf_file = 'tests/data/chr22_ENST00000319363_3UTR.vcf.gz'
chr22_3UTR_ref_seq = 'tests/data/chr22_ENST00000319363_3UTR.ref_seq.txt'
chr22_3UTR_alt_seqs = 'tests/data/chr22_ENST00000319363_3UTR.alt_seqs.txt'


def test_single_variant_5UTR_dataLoader():
    single_variant_5UTR_dataLoader = SingleVariantUTRDataLoader(chr22_gtf_file, chr22_fasta_file, chr22_5UTR_vcf_file,
                                                                feature_type="5UTR")

    units = list(single_variant_5UTR_dataLoader)

    assert len(units) == 4

    expected_ref_seq = open(chr22_5UTR_ref_seq).read().splitlines()[0]
    expected_alt_seqs = open(chr22_5UTR_alt_seqs).read().splitlines()
    for unit, expected_alt_seq in zip(units, expected_alt_seqs):
        assert unit['metadata']['transcript_id'] == 'ENST00000319363'
        assert unit['inputs']['ref_seq'] == expected_ref_seq
        assert unit['inputs']['alt_seq'] == expected_alt_seq

        assert type(unit['metadata']) == dict
        assert len(unit) == 2
        assert len(unit['inputs']) == 2
        assert len(unit['metadata']) == 2  # number of columns

    assert units[0]['metadata']['variant']['ref'] == "G"
    assert units[0]['metadata']['variant']['alt'] == "C"
    assert units[0]['metadata']['variant']['pos'] == 17565851

    assert units[1]['metadata']['variant']['ref'] == "G"
    assert units[1]['metadata']['variant']['alt'] == "A"
    assert units[1]['metadata']['variant']['pos'] == 17565853

    assert units[2]['metadata']['variant']['ref'] == "TGAA"
    assert units[2]['metadata']['variant']['alt'] == "C"
    assert units[2]['metadata']['variant']['pos'] == 17565854

    assert units[3]['metadata']['variant']['ref'] == "T"
    assert units[3]['metadata']['variant']['alt'] == "CAGG"
    assert units[3]['metadata']['variant']['pos'] == 17565854


def test_single_variant_3UTR_dataLoader():
    single_variant_3UTR_dataLoader = SingleVariantUTRDataLoader(chr22_gtf_file, chr22_fasta_file, chr22_3UTR_vcf_file,
                                                                feature_type="3UTR")

    units = list(single_variant_3UTR_dataLoader)

    assert len(units) == 4

    expected_ref_seq = open(chr22_3UTR_ref_seq).read().splitlines()[0]
    expected_alt_seqs = open(chr22_3UTR_alt_seqs).read().splitlines()
    for unit, expected_alt_seq in zip(units, expected_alt_seqs):
        assert unit['metadata']['transcript_id'] == 'ENST00000319363'
        assert unit['inputs']['ref_seq'] == expected_ref_seq
        assert unit['inputs']['alt_seq'] == expected_alt_seq

        assert type(unit['metadata']) == dict
        assert len(unit) == 2
        assert len(unit['inputs']) == 2
        assert len(unit['metadata']) == 2  # number of columns

    assert units[0]['metadata']['variant']['ref'] == "AAAATAAAT"
    assert units[0]['metadata']['variant']['alt'] == "A"
    assert units[0]['metadata']['variant']['pos'] == 17596132

    assert units[1]['metadata']['variant']['ref'] == "T"
    assert units[1]['metadata']['variant']['alt'] == "A"
    assert units[1]['metadata']['variant']['pos'] == 17596136

    assert units[2]['metadata']['variant']['ref'] == "T"
    assert units[2]['metadata']['variant']['alt'] == "C"
    assert units[2]['metadata']['variant']['pos'] == 17596167

    assert units[3]['metadata']['variant']['ref'] == "T"
    assert units[3]['metadata']['variant']['alt'] == "TATG"
    assert units[3]['metadata']['variant']['pos'] == 17596167
