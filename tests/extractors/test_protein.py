import pytest
from pytest_mock import mocker
import pandas as pd
from kipoiseq.transforms.functional import translate, rc_dna
from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.extractors.protein import cut_transcript_seq, gtf_row2interval, \
    CDSFetcher, TranscriptSeqExtractor, ProteinSeqExtractor, \
    ProteinVCFSeqExtractor, SingleSeqProteinVCFSeqExtractor, \
    SingleVariantProteinVCFSeqExtractor

gtf_file = 'tests/data/sample_1_protein.gtf'
fasta_file = 'tests/data/demo_dna_seq.fa'
transcript_id = 'enst_test1'
vcf_file = 'tests/data/singleVar_vcf_enst_test2.vcf.gz'


intervals = [
    Interval('22', 580, 596, strand='+', attrs={'tag': 'cds_end_NF'}),
    Interval('22', 597, 610, strand='+', attrs={'tag': 'cds_end_NF'})
]


def test_cut_seq():
    seq = 'ATCGATG'
    seq = cut_transcript_seq(seq, 'cds_end_NF')
    assert len(seq) == 6


def test_gtf_row2interval():
    row = pd.Series({
        'Chromosome': '22',
        'Start': 10,
        'End': 20,
        'Strand': '-',
        'tag': 'cds_end_NF'
    })
    expected_interval = Interval(chrom='22', start=10,
                                 end=20, name='', strand='-', attrs={'tag': 'cds_end_NF'})

    assert gtf_row2interval(row) == expected_interval


def test_CDSFetcher__read_cds():
    cds = CDSFetcher._read_cds(gtf_file)
    assert cds.shape[0] == 7

    assert cds.iloc[0].Chromosome == '22'
    assert cds.iloc[0].Start == 598
    assert cds.iloc[0].End == 3050

    assert cds.iloc[3].Start == 3
    assert cds.iloc[3].End == 300


@pytest.fixture
def cds_fetcher():
    return CDSFetcher(gtf_file)


def test_CDSFetcher__len__(cds_fetcher):
    assert len(cds_fetcher) == 3


def test_CDSFetcher_get_cds(cds_fetcher):
    intervals = cds_fetcher.get_cds(transcript_id)
    intervals[0] == Interval(chrom='22', start=598,
                             end=3196, name='', strand='+')
    # TODO: Improve testcase with adding transcript with 2 cds


@pytest.fixture
def transcript_seq_extractor():
    return TranscriptSeqExtractor(gtf_file, fasta_file)


def test_get_protein_seq(transcript_seq_extractor):
    transcript_id = 'enst_test2'
    seq = transcript_seq_extractor.get_protein_seq(transcript_id)
    txt_file = 'tests/data/Output_singleSeq_vcf_enst_test2.txt'
    expected_seq = open(txt_file).readline()
    assert seq[1:] == expected_seq[1:]  # no expected mutation here


def test_TranscriptSeqExtractor_prepare_seq():
    seqs = ['ATCGATG']
    assert 'ATCGAT' == TranscriptSeqExtractor._prepare_seq(
        seqs, '+', 'cds_end_NF')
    assert 'CATCGA' == TranscriptSeqExtractor._prepare_seq(
        seqs, '-', 'cds_end_NF')


def test_TranscriptSeqExtractor_get_seq(transcript_seq_extractor):
    seq = transcript_seq_extractor.get_seq(transcript_id)
    assert len(seq) == 3196 - 598


def test_TranscriptSeqExtractor_get_item(transcript_seq_extractor):
    assert transcript_seq_extractor[0] == transcript_seq_extractor.get_seq(
        transcript_id)


@pytest.fixture
def protein_seq_extractor():
    return ProteinSeqExtractor(gtf_file, fasta_file)


def test_ProteinSeqExtractor_prepare_seq(protein_seq_extractor):
    seqs = ['ATCGATG']

    pro_seq = protein_seq_extractor._prepare_seq(seqs, '+', 'cds_end_NF')
    assert pro_seq == 'ID'

    pro_seq = protein_seq_extractor._prepare_seq(seqs, '-', 'cds_end_NF')
    assert pro_seq == 'HR'


def test_ProteinVCFSeqExtractor__unstrand():
    unstrand_intervals = ProteinVCFSeqExtractor._unstrand(intervals)
    assert all(i.strand == '.' for i in unstrand_intervals)

# TODO: write test for with sample_id


@pytest.fixture
def protein_vcf_seq(mocker):
    extractor = ProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
    extractor.extract_query = mocker.MagicMock(
        return_value=iter((['ATC', 'GATG'], ['CATC', 'GAT'])))
    return extractor


def test_ProteinVCFSeqExtractor_extract_cds(protein_vcf_seq):
    protein_seqs = list(protein_vcf_seq.extract_cds(intervals))

    assert protein_seqs[0] == 'ID'
    assert protein_seqs[1] == 'HR'

    query = list(protein_vcf_seq.extract_query
                 .call_args[0][0].variant_intervals)

    variants = list(query[0][0])
    assert len(variants) == 1
    assert variants[0].pos == 596
    interval = query[0][1]
    assert interval.start == 580

    variants = list(query[1][0])

    assert len(variants) == 1
    assert variants[0].pos == 598
    interval = query[1][1]
    assert interval.start == 597


def test_ProteinVCFSeqExtractor_extract(protein_vcf_seq):
    transcript_id = 'enst_test2'
    protein_seqs = list(protein_vcf_seq.extract(transcript_id))
    assert protein_seqs[0] == 'HR'
    assert protein_seqs[1] == 'ID'


@pytest.fixture
def single_seq_protein():
    vcf_file = 'tests/data/singleVar_vcf_enst_test2.vcf.gz'
    return SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)


def test_SingleSeqProteinVCFSeqExtractor_extract(single_seq_protein, transcript_seq_extractor):
    transcript_id = 'enst_test2'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/Output_singleSeq_vcf_enst_test2.txt'
    expected_seq = open(txt_file).readline()
    assert seq == expected_seq

    vcf_file = 'tests/data/singleVar_vcf_enst_test1_diff_type_of_variants.vcf.gz'
    transcript_id = 'enst_test1'
    single_seq_protein = SingleSeqProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    seq = single_seq_protein.extract(transcript_id)
    ref_seq = transcript_seq_extractor.get_protein_seq(transcript_id)

    assert len(seq) == len(ref_seq)
    count = diff_between_two_seq(seq, ref_seq)
    assert count == 1, 'Expected diff of 1 AA, but it was: '+str(count)

    """ no mutationn for enst_test1
    transcript_id = 'enst_test1'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/dna_seq_enst_test1.txt'
    expected_seq = translate(cut_transcript_seq(open(txt_file).readline(), 'cds_end_NF'))
    assert seq == expected_seq
    """

    """ this test doesnt make sense, no mutations in the vcf file
    for current test gtf file
    vcf_file = 'tests/data/singleSeq_vcf_enst_test2.vcf.gz'

    single_seq_protein = SingleSeqProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    transcript_id = 'enst_test2'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/dna_seq_enst_test2.txt'
    expected_seq = translate(cut_transcript_seq(
        rc_dna(open(txt_file).readline()), 'cds_end_NF'))
    assert seq == expected_seq
    """


@pytest.fixture
def single_variant_seq():
    vcf_file = 'tests/data/singleVar_vcf_enst_test2.vcf.gz'
    return SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)


def diff_between_two_seq(seq1, seq2):
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count


def test_SingleVariantProteinVCFSeqExtractor_extract(single_variant_seq, transcript_seq_extractor):
    transcript_id = 'enst_test2'
    seqs = list(single_variant_seq.extract(transcript_id))
    txt_file = 'tests/data/Output_singleVar_vcf_enst_test2.txt'
    expected_seq = open(txt_file).read().splitlines()
    assert seqs[0] == expected_seq[0]
    assert seqs[1] == expected_seq[1]
    assert seqs[2] == expected_seq[2]

    seqs = list(single_variant_seq.extract_all())
    counter = 0
    for t_id in seqs:
        t_id = list(t_id)
        counter += len(t_id)
        for i, seq in enumerate(t_id):
            assert seq == expected_seq[i]
    assert counter == 3, 'Number of variants in vcf 3, but # of seq was: ' + \
        str(counter)

    vcf_file = 'tests/data/singleVar_vcf_enst_test1_diff_type_of_variants.vcf.gz'
    transcript_id = 'enst_test1'
    single_seq_protein = SingleVariantProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    seqs = list(single_seq_protein.extract(transcript_id))
    ref_seq = transcript_seq_extractor.get_protein_seq(transcript_id)

    assert len(seqs) == 1
    for seq in seqs:
        assert len(seq) == len(ref_seq)
        count = diff_between_two_seq(seq, ref_seq)
        assert count == 1, 'Expected diff of 1 AA, but it was: '+str(count)

# TODO: add for all proteins.pep.all.fa
