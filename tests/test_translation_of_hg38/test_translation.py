from tqdm import tqdm_notebook as tqdm
from kipoiseq.extractors.protein import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor, SingleSeqProteinVCFSeqExtractor
from kipoiseq.transforms.functional import translate
from pyfaidx import Fasta
from conftest import gtf_file_GRCh38, fasta_file_GRCh38, vcf_file_for_testing_synonymous_mutations, protein_file_GRCh38
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
def tse():
    return TranscriptSeqExtractor(gtf_file_GRCh38, fasta_file_GRCh38)


def read_pep_fa(protein_file):
    import pandas as pd
    proteins = Fasta(str(protein_file))
    pl = []
    for v in proteins:
        names = v.long_name.split(" ", 8)
        d = {"protein_id": names[0], 'protein_type': names[1]}
        d = {**d, **dict([n.split(":", 1) for n in names[2:]])}
        d['seq'] = str(proteins[v.name])
        pl.append(d)
    return pd.DataFrame(pl)


@pytestmark_gtf
@pytestmark_fasta
@pytestmark_vcf
@pytestmark_protein
@pytest.mark.xfail
def test_hg38(tse):
    dfp = read_pep_fa(protein_file)
    dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
    #assert not dfp['transcript_id'].duplicated().any()
    dfp = dfp.set_index("transcript_id")
    #dfp = dfp[~dfp.chromosome.isnull()]
    assert len(tse) > 100
    assert tse.transcripts.isin(dfp.index).all()
    div3_error = 0
    seq_mismatch_err = 0
    err_transcripts = []
    for transcript_id in tqdm(tse.transcripts):
        # make sure all ids can be found in the proteome
        dna_seq = tse.get_seq(transcript_id)
        if dna_seq == "NNN":
            print(f"{transcript_id} has an ambiguous start and end.Skip!")
            continue
        # dna_seq = dna_seq[:(len(dna_seq) // 3) * 3]
        # if len(dna_seq) % 3 != 0:
         #   div3_error += 1
          #  print("len(dna_seq) % 3 != 0: {}".format(transcript_id))
           # err_transcripts.append({"transcript_id": transcript_id, "div3_err": True})
            # continue
        if len(dna_seq) % 3 != 0:
            print(transcript_id)
            continue
        prot_seq = translate(dna_seq, hg38=True)
        if dfp.loc[transcript_id].seq != prot_seq:
            seq_mismatch_err += 1
            print("seq.mismatch: {}".format(transcript_id))
            n_mismatch = 0
            for i in range(len(prot_seq)):
                a = dfp.loc[transcript_id].seq[i]
                b = prot_seq[i]
                if a != b:
                    n_mismatch += 1
                    print("{} {} {}/{}".format(a, b, i, len(prot_seq)))
            err_transcripts.append({"transcript_id": transcript_id, "div3_err": False,
                                    "n-seq-mismatch": n_mismatch})
