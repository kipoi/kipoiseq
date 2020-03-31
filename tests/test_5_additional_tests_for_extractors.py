from tqdm import tqdm_notebook as tqdm
from kipoiseq.extractors.protein import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor, SingleSeqProteinVCFSeqExtractor
from kipoiseq.transforms.functional import translate
from pyfaidx import Fasta
from conftest import gtf_file_GRCh38, fasta_file_GRCh38, vcf_file_for_testing_synonymous_mutations, protein_file_GRCh38, uniprot_seq_ref
import pytest
import os


pytestmark_gtf = pytest.mark.skipif(not os.path.isfile(gtf_file_GRCh38),
                                    reason="File does not exist")
pytestmark_fasta = pytest.mark.skipif(not os.path.isfile(fasta_file_GRCh38),
                                      reason="File does not exist")
pytestmark_vcf = pytest.mark.skipif(not os.path.isfile(vcf_file_for_testing_synonymous_mutations),
                                    reason="File does not exist")
pytestmark_protein = pytest.mark.skipif(not os.path.isfile(protein_file_GRCh38),
                                        reason="File does not exist")


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
@pytest.fixture
def ssp():
    return SingleSeqProteinVCFSeqExtractor(gtf_file_GRCh38, fasta_file_GRCh38, vcf_file_for_testing_synonymous_mutations)


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
@pytest.fixture
def svp():
    return SingleVariantProteinVCFSeqExtractor(gtf_file_GRCh38, fasta_file_GRCh38, vcf_file_for_testing_synonymous_mutations)


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
@pytest.fixture
def tse():
    return TranscriptSeqExtractor(gtf_file_GRCh38, fasta_file_GRCh38)


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
def test_vcf_single_variant_synonymous_mutations(tse, svp):
    transcript_id = 'ENST00000356175'
    ref_seq = translate(tse.get_seq(transcript_id), True)
    single_var_seq = list(svp.extract(transcript_id))
    for seq in single_var_seq:
        assert seq == ref_seq, seq
    assert len(single_var_seq) == 337, 'Number of sequences != number of variants'

    count = 0
    single_var_seq = list(svp.extract_all())
    for t_id in single_var_seq:
        count += len(list(t_id))

    assert count == 825


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
def test_vcf_single_seq_variants(ssp):
    single_seqs = list(ssp.extract_all())
    assert len(single_seqs) == 5


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
def test_ensembl_uniprot_seq(tse):
    id_and_seq = {}
    with open(uniprot_seq_ref, 'r+') as f:
        key = ""
        for line in f:
            if '>' in line:
                key = (line.replace('>', '')).rstrip()
            else:
                id_and_seq[key] = line.rstrip()

    for transkript_id, ref_seq in tqdm(id_and_seq.items()):
        test_seq = translate(tse.get_seq(transkript_id), True)
        assert test_seq == ref_seq, test_seq
