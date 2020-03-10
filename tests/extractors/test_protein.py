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
transcript_id = 'ENST00000319363'
vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'


intervals = [
    Interval('22', 580, 596, strand='+'),
    Interval('22', 597, 610, strand='+')
]


def test_cut_seq():
    seq = 'ATCGATG'
    seq = cut_transcript_seq(seq, 'normal_tag')
    assert len(seq) == 6


def test_gtf_row2interval():
    row = pd.Series({
        'Chromosome': '22',
        'Start': 10,
        'End': 20,
        'Strand': '-'
    })
    expected_interval = Interval(chrom='22', start=10,
                                 end=20, name='', strand='-')

    assert gtf_row2interval(row) == expected_interval


def test_CDSFetcher__read_cds():
    cds = CDSFetcher._read_cds(gtf_file)
    assert cds.shape[0] == 2

    assert cds.iloc[0].Chromosome == '22'
    assert cds.iloc[0].Start == 598
    assert cds.iloc[0].End == 3196

    assert cds.iloc[1].Start == 3
    assert cds.iloc[1].End == 598


@pytest.fixture
def cds_fetcher():
    return CDSFetcher(gtf_file)


def test_CDSFetcher__len__(cds_fetcher):
    assert len(cds_fetcher) == 2


def test_CDSFetcher_get_cds(cds_fetcher):
    intervals = cds_fetcher.get_cds(transcript_id)
    intervals[0] == Interval(chrom='22', start=598,
                             end=3196, name='', strand='+')
    # TODO: Improve testcase with adding transcript with 2 cds


@pytest.fixture
def transcript_seq_extractor():
    return TranscriptSeqExtractor(gtf_file, fasta_file)


def test_TranscriptSeqExtractor_prepare_seq():
    seqs = ['ATCGATG']
    assert 'ATCGAT' == TranscriptSeqExtractor._prepare_seq(seqs, '+', 'normal_tag')
    assert 'CATCGA' == TranscriptSeqExtractor._prepare_seq(seqs, '-', 'normal_tag')


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

    pro_seq = protein_seq_extractor._prepare_seq(seqs, '+', 'normal_tag')
    assert pro_seq == 'ID'

    pro_seq = protein_seq_extractor._prepare_seq(seqs, '-', 'normal_tag')
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
    assert len(variants) == 2
    assert variants[0].pos == 597
    assert variants[1].pos == 598
    interval = query[1][1]
    assert interval.start == 597


def test_ProteinVCFSeqExtractor_extract(protein_vcf_seq):
    protein_seqs = list(protein_vcf_seq.extract(transcript_id))
    assert protein_seqs[0] == 'ID'
    assert protein_seqs[1] == 'HR'


@pytest.fixture
def single_seq_protein():
    vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'
    return SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)


def test_SingleSeqProteinVCFSeqExtractor_extract(single_seq_protein):
    transcript_id = 'ENST00000381176'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/Output_singleSeq_vcf_ENST000000381176.txt'
    expected_seq = open(txt_file).readline()
    assert seq == expected_seq
    
    #import pdb
    #pdb.set_trace()

    transcript_id = 'ENST00000319363'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/dna_seq_ENST00000319363.txt'
    expected_seq = translate(cut_transcript_seq(open(txt_file).readline()))
    assert seq == expected_seq

    vcf_file = 'tests/data/singleSeq_vcf_ENST000000381176.vcf.gz'
    single_seq_protein = SingleSeqProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    transcript_id = 'ENST00000381176'
    seq = single_seq_protein.extract(transcript_id)
    txt_file = 'tests/data/dna_seq_ENST00000381176.txt'
    expected_seq = translate(cut_transcript_seq(
        rc_dna(open(txt_file).readline())))
    assert seq == expected_seq


@pytest.fixture
def single_variant_seq():
    vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'
    return SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)


def test_SingleVariantProteinVCFSeqExtractor_extract(single_variant_seq):
    transcript_id = 'ENST00000381176'
    seqs = list(single_variant_seq.extract(transcript_id))
    txt_file = 'tests/data/Output_singleVar_vcf_ENST000000381176.txt'
    expected_seq = open(txt_file).read().splitlines()
    assert seqs[0] == expected_seq[0]
    assert seqs[1] == expected_seq[1]
    assert seqs[2] == expected_seq[2]


# TODO: add for all proteins.pep.all.fa
