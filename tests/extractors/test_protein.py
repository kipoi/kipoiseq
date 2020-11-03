import pytest
from pytest_mock import mocker
import pandas as pd
from kipoiseq.transforms.functional import translate, rc_dna
from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.extractors.protein import cut_transcript_seq, TranscriptSeqExtractor, ProteinSeqExtractor, \
    ProteinVCFSeqExtractor, SingleSeqProteinVCFSeqExtractor, \
    SingleVariantProteinVCFSeqExtractor
from kipoiseq.extractors import CDSFetcher, gtf_row2interval, UTRFetcher

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

    seq = 'ATCGATG'
    seq = cut_transcript_seq(seq, 'cds_end_NF,cds_start_NF')
    assert len(seq) == 3

    seq = 'ATCGATG'
    seq = cut_transcript_seq(seq, 'cds_start_NF')
    assert len(seq) == 9

    seq = 'ATCGATG'
    seq = cut_transcript_seq(seq, 'no_tag')
    assert len(seq) == 3


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
    cds = CDSFetcher._read_cds(gtf_file, duplicate_attr=True)
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
    intervals = cds_fetcher.get_intervals(transcript_id)
    intervals[0] == Interval(chrom='22', start=598, end=3196, name='', strand='+')
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
        seqs, intervals, '+')
    assert 'CATCGA' == TranscriptSeqExtractor._prepare_seq(
        seqs, intervals, '-')


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

    pro_seq = protein_seq_extractor._prepare_seq(seqs, intervals, '+')
    assert pro_seq == 'ID'

    pro_seq = protein_seq_extractor._prepare_seq(seqs, intervals, '-')
    assert pro_seq == 'HR'


def test_ProteinVCFSeqExtractor__unstrand():
    unstrand_intervals = ProteinVCFSeqExtractor._unstrand(intervals)
    assert all(i.strand == '.' for i in unstrand_intervals)


# TODO: write test for with sample_id

class TestExtractor(ProteinVCFSeqExtractor):
    def extract_query(
            self,
            variant_interval_queryable,
            ref_seqs,
            intervals,
            reverse_complement: bool,
            sample_id=None
    ):
        for seqs, variant_info in [
            (['ATC', 'GATG'], ['Var_Mutation_Mock']),
            (['CATC', 'GAT'], ['Var_Mutation_Mock']),
        ]:
            alt_seq = self._prepare_seq(
                seqs=seqs,
                intervals=intervals,
                reverse_complement=reverse_complement
            )
            yield alt_seq, variant_info


@pytest.fixture
def protein_vcf_seq(mocker):
    extractor = TestExtractor(gtf_file, fasta_file, vcf_file)
    return extractor


def test_ProteinVCFSeqExtractor_extract_cds(protein_vcf_seq):
    protein_ref_seq, protein_alt_seqs = protein_vcf_seq.extract(intervals)
    protein_alt_seqs = list(protein_alt_seqs)

    assert protein_alt_seqs[0][0] == 'ID'
    assert protein_alt_seqs[1][0] == 'HR'


def test_ProteinVCFSeqExtractor_correctquery(mocker):
    protein_vcf_seq = TestExtractor(gtf_file, fasta_file, vcf_file)
    protein_vcf_seq.extract_query = mocker.MagicMock(
        return_value=(
            'LATTGLWGP',
            iter([
                ('ID', ['Var_Mutation_Mock']),
                ('HR', ['Var_Mutation_Mock'])
            ])
        )
    )
    protein_ref_seq, protein_alt_seqs = protein_vcf_seq.extract(intervals)

    query = list(protein_vcf_seq.extract_query.call_args[1]["variant_interval_queryable"].variant_intervals)

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
    protein_ref_seq, protein_alt_seqs = protein_vcf_seq.get_seq(transcript_id)
    protein_alt_seqs = list(protein_alt_seqs)
    assert protein_alt_seqs[0][0] == 'HR'
    assert protein_alt_seqs[1][0] == 'ID'


@pytest.fixture
def single_seq_protein():
    vcf_file = 'tests/data/singleVar_vcf_enst_test2.vcf.gz'
    return SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)


def test_SingleSeqProteinVCFSeqExtractor_extract(single_seq_protein, transcript_seq_extractor):
    transcript_id = 'enst_test2'
    txt_file = 'tests/data/Output_singleSeq_vcf_enst_test2.txt'
    expected_seq = open(txt_file).readline()

    ref_seq, (alt_seq, variant_info) = single_seq_protein.get_seq(transcript_id)
    assert alt_seq == expected_seq

    vcf_file = 'tests/data/singleVar_vcf_enst_test1_diff_type_of_variants.vcf.gz'
    transcript_id = 'enst_test1'
    single_seq_protein = SingleSeqProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    ref_seq, (alt_seq, variant_info) = single_seq_protein.get_seq(transcript_id)
    assert ref_seq == transcript_seq_extractor.get_protein_seq(transcript_id)
    assert len(alt_seq) == len(ref_seq)
    count = diff_between_two_seq(alt_seq, ref_seq)
    assert count == 1, 'Expected diff of 1 AA, but it was: ' + str(count)

    vcf_file = 'tests/data/singleSeq_vcf_enst_test2.vcf.gz'
    single_seq_protein = SingleSeqProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)

    # transcripts without variants return the reference sequence
    alt_seqs = [alt_seq for t_id, (ref_seq, (alt_seq, variants)) in single_seq_protein.extract_all() if
                len(variants) > 0]
    assert len(alt_seqs) == 0


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
    txt_file = 'tests/data/Output_singleVar_vcf_enst_test2.txt'
    expected_seq = open(txt_file).read().splitlines()

    # test single
    transcript_id = 'enst_test2'
    ref_seq, alt_seqs = single_variant_seq.get_seq(transcript_id)
    alt_seqs = list(alt_seqs)
    assert alt_seqs[0][0] == expected_seq[0]
    assert alt_seqs[1][0] == expected_seq[1]
    assert alt_seqs[2][0] == expected_seq[2]

    # test multiple
    transcript_id = ['enst_test1', 'enst_test2']
    transcript_seqs = single_variant_seq.get_seq(transcript_id)
    assert isinstance(transcript_seqs, list)
    transcript_seqs = [list(alt_seqs) for ref_seq, alt_seqs in transcript_seqs]
    assert transcript_seqs[1][0][0] == expected_seq[0]
    assert transcript_seqs[1][1][0] == expected_seq[1]
    assert transcript_seqs[1][2][0] == expected_seq[2]

    transcript_seqs = single_variant_seq.iter_seq(transcript_id)
    assert not isinstance(transcript_seqs, list)
    transcript_seqs = [list(alt_seqs) for ref_seq, alt_seqs in transcript_seqs]
    assert transcript_seqs[1][0][0] == expected_seq[0]
    assert transcript_seqs[1][1][0] == expected_seq[1]
    assert transcript_seqs[1][2][0] == expected_seq[2]

    counter = 0
    for tr_id, (ref_seq, t_id_seqs) in single_variant_seq.extract_all():
        t_id_seqs = [seq for seq, info in list(t_id_seqs)]
        if len(t_id_seqs) == 0:
            continue
        counter += len(t_id_seqs)
        for i, seq in enumerate(t_id_seqs):
            assert seq == expected_seq[i]
        assert tr_id == 'enst_test2'
    assert counter == 3, 'Number of variants in vcf 3, but # of seq was: ' + \
                         str(counter)

    vcf_file = 'tests/data/singleVar_vcf_enst_test1_diff_type_of_variants.vcf.gz'
    transcript_id = 'enst_test1'
    single_var_protein = SingleVariantProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)
    ref_seq, alt_seqs = single_var_protein.get_seq(transcript_id)
    alt_seqs = list(alt_seqs)
    # alt_seqs = [seq for seq, info in list(single_var_protein.get_seq(transcript_id))]
    # ref_seq = transcript_seq_extractor.get_protein_seq(transcript_id)

    assert len(alt_seqs) == 1
    for alt_seq, variant_info in alt_seqs:
        assert len(alt_seq) == len(ref_seq)
        count = diff_between_two_seq(alt_seq, ref_seq)
        assert count == 1, 'Expected diff of 1 AA, but it was: ' + str(count)

    # this test should result in 0 sequences yielded
    vcf_file = 'tests/data/singleSeq_vcf_enst_test2.vcf.gz'
    single_var_protein = SingleVariantProteinVCFSeqExtractor(
        gtf_file, fasta_file, vcf_file)
    length = 0
    for tr_id, (ref_seq, t_id_seqs) in single_var_protein.extract_all():
        t_id_seqs = [seq for seq, info in list(t_id_seqs)]
        length += len(t_id_seqs)
    assert length == 0


# TODO: add for all proteins.pep.all.fa

# chr22_fasta_file = 'tests/data/chr22.fa.gz'
chr22_gtf_file = 'tests/data/chr22_ENST00000319363.gtf'


# chr22_5UTR_vcf_file = 'tests/data/chr22_ENST00000319363_5UTR.vcf.gz'


def test_5UTRFetcher__read_utr():
    utr5 = UTRFetcher._read_utr(chr22_gtf_file, feature_type="5UTR")

    assert utr5.shape == (1, 12)

    assert utr5.iloc[0].Chromosome == 'chr22'
    assert utr5.iloc[0].Start == 17565848
    assert utr5.iloc[0].End == 17565981
    assert utr5.iloc[0].Strand == "+"

    utr5_from_cds = UTRFetcher._read_utr(chr22_gtf_file, feature_type="5UTR", infer_from_cds=True)

    pd.testing.assert_frame_equal(left = utr5.drop(['exon_number', 'exon_id'], axis=1), right = utr5_from_cds.drop(['exon_number', 'exon_id'], axis=1), check_dtype=False)


def test_3UTRFetcher__read_utr():
    utr3 = UTRFetcher._read_utr(chr22_gtf_file, feature_type="3UTR")

    assert utr3.shape == (1, 12)

    assert utr3.iloc[0].Chromosome == 'chr22'
    assert utr3.iloc[0].Start == 17590710
    assert utr3.iloc[0].End == 17596583
    assert utr3.iloc[0].Strand == "+"

    utr3_from_cds = UTRFetcher._read_utr(chr22_gtf_file, feature_type="3UTR", infer_from_cds=True)

    pd.testing.assert_frame_equal(left=utr3.drop(['exon_number', 'exon_id'], axis=1),
                                  right=utr3_from_cds.drop(['exon_number', 'exon_id'], axis=1), check_dtype=False)
