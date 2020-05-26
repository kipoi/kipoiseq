import pytest
from kipoiseq.dataloaders import SingleVariantProteinDataLoader


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
    assert units[0]['input']['ref_seq'][1:] == expected_seq[1:]

    txt_file = 'tests/data/Output_singleVar_vcf_enst_test2.txt'
    expected_seq = open(txt_file).read().splitlines()
    assert units[0]['input']['alt_seq'] == expected_seq[0]
    assert units[1]['input']['alt_seq'] == expected_seq[1]
    assert units[2]['input']['alt_seq'] == expected_seq[2]

    assert type(units[2]['metadata']) == dict
    assert len(units[2]) == 2
    assert len(units[2]['input']) == 2
    assert len(units[2]['metadata']) == 18  # number of columns
    
    assert units[0]['metadata']['variants']['ref'] == "C"
    assert units[0]['metadata']['variants']['alt'] == "G"
    assert units[0]['metadata']['variants']['pos'] == 596
    assert units[0]['metadata']['variants']['id'] == "9"
    
    assert units[1]['metadata']['variants']['ref'] == "A"
    assert units[1]['metadata']['variants']['alt'] == "G"
    assert units[1]['metadata']['variants']['pos'] == 597
    assert units[1]['metadata']['variants']['id'] == "10"
    
    assert units[2]['metadata']['variants']['ref'] == "T"
    assert units[2]['metadata']['variants']['alt'] == "G"
    assert units[2]['metadata']['variants']['pos'] == 598
    assert units[2]['metadata']['variants']['id'] == "11"
