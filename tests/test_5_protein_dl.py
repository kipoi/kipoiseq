"""Test protein dataloader
"""
from pybedtools import Interval
import pytest
import pyranges as pr
from concise.utils.fasta import read_fasta
import pandas as pd
import numpy as np
from kipoiseq.transforms.functional import translate, rc_dna
from kipoiseq.extractors import FastaStringExtractor
from kipoiseq.extractors.vcf_seq import SingleSeqVCFSeqExtractor,SingleVariantVCFSeqExtractor
from tqdm import tqdm


from pathlib import Path
ddir = Path('/s/genomes/human/hg19/ensembl_GRCh37.p13_release75')

gtf_file = ddir / 'Homo_sapiens.GRCh37.75.chr22.gtf'
gtf_full_file = ddir / 'Homo_sapiens.GRCh37.75.gtf'
fasta_file = ddir / 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
protein_file = ddir / 'Homo_sapiens.GRCh37.75.pep.all.fa'


# gff_file = 'data/protein/Homo_sapiens.GRCh38.97.chromosome.22.gff3.gz'
# gtf_file = 'data/protein/Homo_sapiens.GRCh38.97.chr.chr22.gtf.gz'
# gtf_full_file = 'data/protein/Homo_sapiens.GRCh38.97.chr.gtf.gz'
# fasta_file = 'data/protein/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
# protein_file = 'data/protein/Homo_sapiens.GRCh38.pep.all.fa'

# gtf = pr.read_gtf(gtf_file, output_df=True)

# gtf_full = pr.read_gtf(gtf_full_file, output_df=True)

# gff = pr.read_gff(gff_file, full=True)

# df = gtf_full

# df = df[(df.transcript_biotype == 'protein_coding')]
# df[(df.transcript_biotype == 'protein_coding') & (~df.protein_id.isna()) & (df.exon_number == '1')]
# np.sum((df.transcript_biotype == 'protein_coding') & (df.Feature == 'transcript'))
# dict(df.Feature.value_counts())
# df[df.Feature == 'transcript'].transcript_biotype.value_counts()
# # 43469

# assert not df[df.Feature == 'transcript'].transcript_id.duplicated().any()

# len(df)


def read_pep_fa(protein_file):
    proteins = read_fasta(protein_file)
    pl = []
    for k, v in proteins.items():
        names = k.split(" ", 8)
        d = {"protein_id": names[0], 'protein_type': names[1]}
        d = {**d, **dict([n.split(":", 1) for n in names[2:]])}
        d['seq'] = v
        pl.append(d)
    return pd.DataFrame(pl)


# dfp = read_pep_fa(protein_file)
# dfp.transcript_biotype.value_counts()
# assert not dfp.transcript.duplicated().any()

# df[df.Feature == 'transcript'].transcript_id.duplicated().any()

# cds = df[(df.Feature == 'CDS')].set_index('transcript_id')

# transcript_id = 'ENST00000252835'
# transcript_id = 'ENST00000395590'


def gtf_row2interval(row):
    """Note: GTF is 1-based
    """
    return Interval(str(row.Chromosome),
                    int(row.Start) - 1,
                    int(row.End),
                    strand=str(row.Strand))


class GenomeCDSSeq:

    def __init__(self, gtf_file, fasta_file):
        """Protein sequences in the genome
        """
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)

        self.fae = FastaStringExtractor(self.fasta_file, use_strand=False)
        df = pr.read_gtf(self.gtf_file, output_df=True)
        biotype_str = 'transcript_biotype' if 'transcript_biotype' in df else 'gene_biotype'
        self.cds = (df
                    .query(f"{biotype_str} == 'protein_coding'")
                    .query("(Feature == 'CDS') | (Feature == 'CCDS')")
                    .query("(tag == 'basic') | (tag == 'CCDS')")
                    .set_index('transcript_id'))
        # filter valid transcripts
        if 'transcript_support_level' in self.cds:
            self.cds = self.cds[~self.cds.transcript_support_level.isnull()]
            self.cds = self.cds[self.cds.transcript_support_level != 'NA']
            self.cds.transcript_support_level = self.cds.transcript_support_level.astype(int)
            self.cds = self.cds[~self.cds.transcript_support_level.isna()]
        else:
            print("Transcript support level not in gtf. Skipping the associated filters.")
        self.cds_exons = pr.PyRanges(self.cds.reset_index())
        self.transcripts = self.cds.index.unique()

    def __len__(self):
        return len(self.transcripts)

    def get_cds_exons(self, transcript_id):
        cds_exons = self.cds.loc[transcript_id]

        # get cds intervals
        if isinstance(cds_exons, pd.Series):
            # single exon
            strand = cds_exons.Strand
            intervals = [gtf_row2interval(cds_exons)]
        else:
            # multiple exons
            strand = cds_exons.iloc[0].Strand
            assert np.all(strand == cds_exons.Strand)

            intervals = [gtf_row2interval(row)
                         for i, row in cds_exons.loc[transcript_id].sort_values("Start").iterrows()]
        return intervals, strand

    def get_seq(self, transcript_id):
        exons, strand = self.get_cds_exons(transcript_id)
        # merge the sequences

        seq = "".join([self.fae.extract(exon) for exon in exons])

        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        while(len(seq)%3!=0):
            seq=seq[:-1]
        return seq

    def get_seq_variants(self, transcript_id, variants):
        exons, strand = self.get_cds_exons(transcript_id)
        # merge the sequences

        # TODO - insert genetic variants here
        seq = "".join([self.fae.extract(exon) for exon in exons])

        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        return seq

    def __getitem__(self, idx):
        return self.get_seq(self.transcripts[idx])

    def overlaped_exons(self, variant):
        """Which exons are overlapped by a variant

        Overall strategy:
        1. given the variant, get all the affected transcripts
        2. Given the transcript and the variants,
          fetch the ref and alt sequences for the transcripts
        """
        # TODO - perform a join between variants and exons
        # https://github.com/gagneurlab/MMSplice/blob/master/mmsplice/vcf_dataloader.py#L136-L190
        # this will generate (variant, cds_exon) pairs
        # cds_exon will contain also the information about the order in the transcript
        return self.cds_exons.join(variant)
    
class VarianceSeq:
    
    def __init__(self,fasta_file,vcf_file,gtf_file):
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)
        self.vcf_file=str(vcf_file)
        self.gps = GenomeCDSSeq(self.gtf_file, self.fasta_file)
        self.sse=SingleSeqVCFSeqExtractor(self.fasta_file, self.vcf_file)
        self.sve=SingleVariantVCFSeqExtractor(self.fasta_file, self.vcf_file)

    def coding_seq_chr(self):
        for transcript_id in tqdm(self.gps.transcripts):
            intervals,strand = self.gps.get_cds_exons(transcript_id)
            seq=""
            for interval in intervals:
                interval = Interval(interval[0],int(interval[1]),int(interval[2]))
                seq="".join([seq,self.sse.extract(interval,anchor=1)])
            if strand == '-':
                seq = rc_dna(seq)
            # optionally reverse complement
            seq=self.cut_seq(seq)
            yield translate(seq)
            
    def cut_seq(self,seq):
        while(len(seq)%3!=0):
            seq=seq[:-1]
        return seq    
    def coding_seq_chr_var(self):
        import itertools

        for transcript_id in tqdm(self.gps.transcripts):
            seq_list=[""]
            intervals,strand = self.gps.get_cds_exons(transcript_id)
            for interval in intervals:
                interval = Interval(interval[0],int(interval[1]),int(interval[2]))
                gen_of_seq=self.sve.extract(interval,anchor=1)
                list_of_seq = []
                for gs in gen_of_seq: list_of_seq.append(gs)
                if len(list_of_seq)==0: list_of_seq.append(self.sse.extract(interval,anchor=1))
                result=itertools.product(seq_list,list_of_seq)
                list_of_seq.clear()
                for r in result: list_of_seq.append(r)
                seq_list.clear()
                for min_seq in list_of_seq: seq_list.append("".join([seq for seq in min_seq]))
            if strand == '-':
            # optionally reverse complement
               seq_list = [rc_dna(seq) for seq in seq_list]
            seq_list=[self.cut_seq(seq) for seq in seq_list]
            seq_list=[translate(seq) for seq in seq_list]
            yield seq_list
            
            
    #function to test seq with all variants
    def test_coding_seq_chr(self,transcript_id):
        intervals,strand = self.gps.get_cds_exons(transcript_id)
        seq=""
        for interval in intervals:
            interval = Interval(interval[0],int(interval[1]),int(interval[2]))
            seq="".join([seq,self.sse.extract(interval,anchor=1)])
        if strand == '-':
            seq = rc_dna(seq)
        # optionally reverse complement
        seq=self.cut_seq(seq)
        return translate(seq)
    
    #function to test single variant in seq
    def test_coding_seq_chr_var(self,transcript_id):
        import itertools
        intervals,strand = self.gps.get_cds_exons(transcript_id)
        seq_list=[""]
        for interval in intervals:
            interval = Interval(interval[0],int(interval[1]),int(interval[2]))
            gen_of_seq=self.sve.extract(interval,anchor=1)
            list_of_seq = []
            for gs in gen_of_seq:
                list_of_seq.append(gs)
            if len(list_of_seq)==0:
                list_of_seq.append(self.sse.extract(interval,anchor=1))
            result=itertools.product(seq_list,list_of_seq)
            list_of_seq.clear()
            for r in result:
                list_of_seq.append(r)
            seq_list.clear()
            for min_seq in list_of_seq: seq_list.append("".join([se for se in min_seq]))
        if strand == '-':
        # optionally reverse complement
            seq_list = [rc_dna(seq) for seq in seq_list]
        seq_list=[self.cut_seq(seq)for seq in seq_list]
        seq_list=[translate(seq) for seq in seq_list]
        return seq_list

<<<<<<< HEAD
def test_mutation_in_each_exon_all_variance():
    from pathlib import Path
    ddir = Path('/s/genomes/human/hg19/ensembl_GRCh37.p13_release75')
    gtf_file = ddir / 'Homo_sapiens.GRCh37.75.chr22.gtf'
    fasta_file = ddir / 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    vcf_file = 'tests/data/test_1.2.vcf.gz'
    protein_file = ddir / 'Homo_sapiens.GRCh37.75.pep.all.fa'
    vs = VarianceSeq(fasta_file,vcf_file,gtf_file)  
    transcript_id = 'ENST00000381176'
    s2 = vs.test_coding_seq_chr(transcript_id)
    dfp = read_pep_fa(protein_file)
    dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
    assert not dfp['transcript_id'].duplicated().any()
    dfp = dfp.set_index("transcript_id")
    dfp = dfp[~dfp.chromosome.isnull()]
    s1 = dfp.loc[transcript_id].seq
    for i in range(len(s1)):
        if s1[i]!=s2[i]:
            assert s2[i] in ['L','A','K'],'Unexepcted AA: '+s2[i]+'. Should be: ' + s1[i]

def test_mutation_single_variance():
    from pathlib import Path
    ddir = Path('/s/genomes/human/hg19/ensembl_GRCh37.p13_release75')
    gtf_file = ddir / 'Homo_sapiens.GRCh37.75.chr22.gtf'
    fasta_file = ddir / 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    vcf_file2 = 'tests/data/test_1.1.vcf.gz'
    protein_file = ddir / 'Homo_sapiens.GRCh37.75.pep.all.fa'
    vs2 = VarianceSeq(fasta_file,vcf_file2,gtf_file)  
    transcript_id = 'ENST00000381176'
    s2 = vs2.test_coding_seq_chr_var(transcript_id)
    dfp = read_pep_fa(protein_file)
    dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
    assert not dfp['transcript_id'].duplicated().any()
    dfp = dfp.set_index("transcript_id")
    dfp = dfp[~dfp.chromosome.isnull()]
    s1 = dfp.loc[transcript_id].seq
    for seq in s2:
        for a in range(len(seq)):
            if s1[a]!=seq[a]:
                assert seq[a] in ['I','T','L'],'Unexepcted AA: '+ seq[a]




"""
=======


>>>>>>> 86ee4a3eb63859d18054898baa6d583e5c6db0b1
dfp = read_pep_fa(protein_file)
dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
assert not dfp['transcript_id'].duplicated().any()
dfp = dfp.set_index("transcript_id")
dfp = dfp[~dfp.chromosome.isnull()]

gps = GenomeCDSSeq(gtf_full_file, fasta_file)
assert len(gps) > 100
assert gps.transcripts.isin(dfp.index).all()

transcript_id = 'ENST00000485079'
div3_error = 0
seq_mismatch_err = 0
err_transcripts = []
for transcript_id in tqdm(gps.transcripts):
    # make sure all ids can be found in the proteome
    dna_seq = gps.get_seq(transcript_id)
    # dna_seq = dna_seq[:(len(dna_seq) // 3) * 3]
    if len(dna_seq) % 3 != 0:
        div3_error += 1
        print(f"len(dna_seq) % 3 != 0: {transcript_id}")
        err_transcripts.append({"transcript_id": transcript_id, "div3_err": True})
        continue
    prot_seq = translate(dna_seq)
    if dfp.loc[transcript_id].seq != prot_seq:
        seq_mismatch_err += 1
        print(f"seq.mismatch: {transcript_id}")
        n_mismatch = 0
        for i in range(len(prot_seq)):
            a = dfp.loc[transcript_id].seq[i]
            b = prot_seq[i]
            if a != b:
                n_mismatch += 1
                print(f"{a} {b} {i}/{len(prot_seq)}")
        err_transcripts.append({"transcript_id": transcript_id, "div3_err": False,
                                "n-seq-mismatch": n_mismatch})
        # print("prot:", dfp.loc[transcript_id].seq)
        # print("seq: ", prot_seq)
err_transcripts = pd.DataFrame(err_transcripts)
#err_cds.to_csv("data/protein/err_cds.csv")
#err_transcripts.to_csv("data/protein/err_transcripts.csv")
# len(dna_seq) % 3 != 0: ENST00000625258
# len(dna_seq) % 3 != 0: ENST00000485079
# len(dna_seq) % 3 != 0: ENST00000638880
# len(dna_seq) % 3 != 0: ENST00000606910
# len(dna_seq) % 3 != 0: ENST00000412982
# seq.mismatch: ENST00000606505
# TODO - run this for all the sequences

# they use U for the stop

# M,L code at the begining


# 359 26 1802

# now: 0 115 1802
#print(div3_error, seq_mismatch_err, len(gps))
# TODO - fix all the errors

#def test_translate():
#    assert translate("TGAATGGAC") == '_MD'
#    assert translate("TTTATGGAC") == 'FMD'
#    with pytest.raises(ValueError):
#        translate("TGAATGGA")
<<<<<<< HEAD
"""
=======
>>>>>>> 86ee4a3eb63859d18054898baa6d583e5c6db0b1
