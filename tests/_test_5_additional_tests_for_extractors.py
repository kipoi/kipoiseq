from pathlib import Path
from tqdm import tqdm_notebook as tqdm
from kipoiseq.extractors.protein import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor, SingleSeqProteinVCFSeqExtractor
from kipoiseq.transforms.functional import translate
from pyfaidx import Fasta


ddir = Path('/s/genomes/human/hg38/ensembl_GRCh38')
gtf_file = ddir / 'Homo_sapiens.GRCh38.99.gtf'
#gtf_file = 'tests/data/38_sample.gtf.txt'
protein_file = ddir / 'Homo_sapiens.GRCh38.pep.all.fa'
fasta_file = ddir / 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
vcf_file = 'data/vcf_file_for_testing_synonymous_mutations.vcf.gz'  # NF1

ssp = SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
svp = SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
tse = TranscriptSeqExtractor(gtf_file, fasta_file)


def test_vcf_single_variant_synonymous_mutations():
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


def test_vcf_single_seq_variants():
    single_seqs = list(ssp.extract_all())
    assert len(single_seqs) == 5


def test_ensembl_uniprot_seq():
    ref_path = '/data/nasif12/home_if12/nonchev/code/kipoiseq/tests/data/uniprot_ref.txt'
    id_and_seq = {}
    with open(ref_path, 'r+') as f:
        key = ""
        for line in f:
            if '>' in line:
                key = (line.replace('>', '')).rstrip()
            else:
                id_and_seq[key] = line.rstrip()

    for transkript_id, ref_seq in tqdm(id_and_seq.items()):
        test_seq = translate(tse.get_seq(transkript_id), True)
        assert test_seq == ref_seq, test_seq


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


"""
def test_hg38():
    dfp = read_pep_fa(protein_file)
    dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
    #assert not dfp['transcript_id'].duplicated().any()
    dfp = dfp.set_index("transcript_id")
    #dfp = dfp[~dfp.chromosome.isnull()]
    assert len(tse) >100
    assert tse.transcripts.isin(dfp.index).all()
    div3_error = 0
    seq_mismatch_err = 0
    err_transcripts = []
    for transcript_id in tqdm(tse.transcripts):
        # make sure all ids can be found in the proteome
        dna_seq = tse.get_seq(transcript_id)
        if dna_seq == "XXX":
            print(f"{transcript_id} has an ambiguous start and end.Skip!")
            continue
        # dna_seq = dna_seq[:(len(dna_seq) // 3) * 3]
        #if len(dna_seq) % 3 != 0:
         #   div3_error += 1
          #  print("len(dna_seq) % 3 != 0: {}".format(transcript_id))
           # err_transcripts.append({"transcript_id": transcript_id, "div3_err": True})
            #continue
        if len(dna_seq)%3 != 0:
            print(transcript_id)
            continue
        prot_seq = translate(dna_seq,hg38=True)
        if dfp.loc[transcript_id].seq != prot_seq:
            seq_mismatch_err += 1
            print("seq.mismatch: {}".format(transcript_id))
            n_mismatch = 0
            for i in range(len(prot_seq)):
                a = dfp.loc[transcript_id].seq[i]
                b = prot_seq[i]
                if a != b:
                    n_mismatch += 1
                    print("{} {} {}/{}".format(a,b,i,len(prot_seq)))
            err_transcripts.append({"transcript_id": transcript_id, "div3_err": False,
                                "n-seq-mismatch": n_mismatch})
            # print("prot:", dfp.loc[transcript_id].seq)
            # print("seq: ", prot_seq)
    err_transcripts = pd.DataFrame(err_transcripts)
"""
