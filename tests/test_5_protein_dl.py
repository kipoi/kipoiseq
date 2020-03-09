from itertools import chain, islice
import numpy as np
from kipoiseq.dataclasses import Interval
from kipoiseq.transforms.functional import rc_dna, translate
from kipoiseq.extractors.base import FastaStringExtractor
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_seq import VariantSeqExtractor
import os.path
import pytest
from pyfaidx import Fasta
import pandas as pd
import numpy as np
from tqdm import tqdm


# TODO: convert print to logs
# TODO: documentation

def read_pep_fa(protein_file):
    proteins = Fasta(str(protein_file))
    pl = []
    for v in proteins:
        names = v.long_name.split(" ", 8)
        d = {"protein_id": names[0], 'protein_type': names[1]}
        d = {**d, **dict([n.split(":", 1) for n in names[2:]])}
        d['seq'] = str(proteins[v.name])
        pl.append(d)
    return pd.DataFrame(pl)

def cut_transcript_seq(seq):
    # if the dna seq is not %3==0, there are unnecessary bases at the end
    # should be called only after all exons are connected!
    while(len(seq) % 3 != 0):
        seq = seq[:-1]
    return seq

def gtf_row2interval(row):
    """
    Convert gtf row object into interval class.
    """
    return Interval(str(row.Chromosome),
                    int(row.Start) - 1,
                    int(row.End),
                    strand=str(row.Strand))


class CDSFetcher:

    def __init__(self, gtf_file):
        """Protein sequences in the genome
        """
        self.gtf_file = str(gtf_file)
        self.cds = self._read_cds(self.gtf_file)
        self.transcripts = self.cds.index.unique()

    @staticmethod
    def _read_cds(gtf_file):
        import pyranges
        df = pyranges.read_gtf(gtf_file, output_df=True)

        cds = CDSFetcher._get_cds_from_gtf(df)
        cds = CDSFetcher._filter_valid_transcripts(cds)
        return cds

    @staticmethod
    def _get_cds_from_gtf(df):
        biotype_str = CDSFetcher._get_biotype_str(df)
        return (df
                .query("{} == 'protein_coding'".format(biotype_str))
                .query("(Feature == 'CDS') | (Feature == 'CCDS')")
                .query("(tag == 'basic') | (tag == 'CCDS')")
                .set_index('transcript_id'))

    @staticmethod
    def _get_biotype_str(df):
        if 'transcript_biotype' in df:
            return 'transcript_biotype'
        elif 'gene_biotype' in df:
            return 'gene_biotype'
        else:
            raise ValueError('Cannot obtain `biotype_str` from gtf file')

    @staticmethod
    def _filter_valid_transcripts(cds):
        if 'transcript_support_level' in cds:
            cds = cds[~cds.transcript_support_level.isnull()]
            cds = cds[cds.transcript_support_level != 'NA']
            cds.transcript_support_level = cds.transcript_support_level.astype(
                int)
            cds = cds[~cds.transcript_support_level.isna()]
        else:
            print('Transcript support level not in gtf.'
                  'Skipping the associated filters.')
        return cds

    def __len__(self):
        return len(self.transcripts)

    def get_cds(self, transcript_id):
        cds = self.cds.loc[[transcript_id]]
        assert np.all(cds.iloc[0].Strand == cds.Strand)

        return [
            gtf_row2interval(row)
            for i, row in cds.sort_values("Start").iterrows()
        ]


class TranscriptSeqExtractor:

    def __init__(self, gtf_file, fasta_file):
        self.fasta_file = str(fasta_file)
        self.gtf_file = str(gtf_file)
        self.fasta = FastaStringExtractor(self.fasta_file, use_strand=False)
        self.cds_fetcher = CDSFetcher(self.gtf_file)
        self.transcripts = self.cds_fetcher.transcripts

    def __len__(self):
        return len(self.cds_fetcher)
    
    def get_seq1(self, transcript_id):
        cds = self.cds_fetcher.get_cds(transcript_id)
        seqs = [self.fasta.extract(i) for i in cds]
        return self._prepare_seq(seqs, cds[0].strand)

    @staticmethod
    def _prepare_seq(seqs, strand):
        seq = "".join(seqs)
        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        seq = cut_transcript_seq(seq)
        return seq
    
    def seq_translation_positive_strand(self, transcript_id, cds):
        interval_ref = cds[len(cds)-1]
        interval_with_end = Interval(interval_ref.chrom, interval_ref.start, interval_ref.end+3, interval_ref.name, interval_ref.score, interval_ref.strand, interval_ref.attrs)
        end_codon = self.fasta.extract(interval_with_end)[-3:]
        if end_codon not in ['TGA', 'TAA', 'TAG', 'TG.','TA.']:
            index = 1
            while True and index<=1:
                interval = Interval(interval_with_end.chrom, interval_with_end.start, interval_with_end.end+index, interval_with_end.name, interval_with_end.score, interval_with_end.strand, interval_with_end.attrs)
                end_codon = self.fasta.extract(interval)[-3:]
                if end_codon in ['TGA', 'TAA', 'TAG']:
                    cds[len(cds)-1] = Interval(interval_ref.chrom, interval_ref.start, interval_ref.end+index, interval_ref.name, interval_ref.score, interval_ref.strand, interval_ref.attrs)
                    break
                index +=1
        seqs = "".join([self.fasta.extract(i) for i in cds])
        flag = False
        if len(seqs) % (1 or 0) != 0 and end_codon in ['TGA', 'TAA', 'TAG', 'TG.','TA.']:
            flag = True
            while True:
                seqs = seqs[1:]
                if len(seqs) % 3 == 0:
                    break
        seq_end_codon_fixed = self._prepare_seq(seqs, cds[0].strand)    
        sequence_acid = translate(seq_end_codon_fixed)
        if '_' in sequence_acid and sequence_acid[0] not in ['M','L']:
            index=1
            while '_' in sequence_acid and index<=2:
                flag = False
                seq = ['X']
                seq.append(translate(seq_end_codon_fixed[index:-1-(index%2)]))
                sequence_acid = "".join(seq)
                index+=1
        if flag:
            sequence_acid = "".join(["X",sequence_acid])
        return sequence_acid
    
    def seq_translation_negative_strand(self, transcript_id, cds):
        interval_ref = cds[0]
        interval_with_end = Interval(interval_ref.chrom, interval_ref.start-3, interval_ref.end, interval_ref.name, interval_ref.score, interval_ref.strand, interval_ref.attrs)
        end_codon = rc_dna(self.fasta.extract(interval_with_end)[:3])
        if end_codon not in ['TGA', 'TAA', 'TAG', 'TG.','TA.']:
            index = 1
            while True and index <= 1:
                interval = Interval(interval_with_end.chrom, interval_with_end.start-index, interval_with_end.end, interval_with_end.name, interval_with_end.score, interval_with_end.strand, interval_with_end.attrs)
                end_codon = rc_dna(self.fasta.extract(interval)[:3])
                if end_codon in ['TGA', 'TAA', 'TAG']:
                    cds[0] = Interval(interval_ref.chrom, interval_ref.start-index, interval_ref.end, interval_ref.name, interval_ref.score, interval_ref.strand, interval_ref.attrs)
                    break
                index +=1
        seqs = "".join([self.fasta.extract(i) for i in cds])
        flag = False
        if len(seqs) % 3 != (1 or 0) and end_codon in ['TGA', 'TAA', 'TAG', 'TG.','TA.']:
            flag = True
            while True:
                seqs = seqs[:-1]
                if len(seqs) % 3 == 0:
                    break
        seq_end_codon_fixed = self._prepare_seq(seqs, cds[0].strand)    
        sequence_acid = translate(seq_end_codon_fixed)
        if '_' in sequence_acid and sequence_acid[0] not in ['M','L']:
            index=1
            while '_' in sequence_acid and index<=2:
                flag = False
                seq = ['X']
                seq.append(translate(seq_end_codon_fixed[index:-1-(index%2)]))
                sequence_acid = "".join(seq)
                index+=1
        if flag:
            sequence_acid = "".join(["X",sequence_acid])
        return sequence_acid

    def get_seq(self, transcript_id):
        cds = self.cds_fetcher.get_cds(transcript_id)
        seqs = [self.fasta.extract(i) for i in cds]
        seq = translate(self._prepare_seq(seqs, cds[0].strand))
        if seq[0] in ['M', 'L', 'V', 'K'] and '_' not in seq:
            return seq
        if cds[0].strand == '-':
            return self.seq_translation_negative_strand(transcript_id, cds)
        else:
            return self.seq_translation_positive_strand(transcript_id, cds)

    def __getitem__(self, idx):
        return self.get_seq(self.transcripts[idx])

    def overlaped_cds(self, variants):
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
        raise NotImplementedError()


class ProteinSeqExtractor(TranscriptSeqExtractor):

    @staticmethod
    def _prepare_seq(seqs, strand):
        return translate(TranscriptSeqExtractor._prepare_seq(seqs, strand))


class ProteinVCFSeqExtractor:

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)
        self.vcf_file = str(vcf_file)
        self.cds_fetcher = CDSFetcher(self.gtf_file)
        self.transcripts = self.cds_fetcher.transcripts
        self.fasta = FastaStringExtractor(self.fasta_file)
        self.multi_sample_VCF = MultiSampleVCF(self.vcf_file)
        self.variant_seq_extractor = VariantSeqExtractor(self.fasta_file)

    @staticmethod
    def _unstrand(intervals):
        return [i.unstrand() for i in intervals]

    def extract_cds(self, cds, sample_id=None):
        intervals = self._unstrand(cds)

        variant_interval_queryable = self.multi_sample_VCF.query_variants(
            intervals, sample_id=sample_id)

        iter_seqs = self.extract_query(variant_interval_queryable,
                                       sample_id=sample_id)

        for seqs in iter_seqs:
            yield ProteinSeqExtractor._prepare_seq(seqs, cds[0].strand)

    def extract(self, transcript_id, sample_id=None):
        return self.extract_cds(self.cds_fetcher.get_cds(transcript_id),
                                sample_id=sample_id)

    def _ref_cds_seq(self, variant_interval_queryable):
        intervals = variant_interval_queryable.iter_intervals()
        return [self.fasta.extract(interval) for interval in intervals]


class SingleSeqProteinVCFSeqExtractor(ProteinVCFSeqExtractor):

    def _extract_query(self, variant_interval_queryable, sample_id=None):
        for variants, interval in variant_interval_queryable.variant_intervals:
            for variant in variants:
                if len(variant.ref) == len(variant.alt):
                    yield self.variant_seq_extractor.extract(
                        interval, variants, anchor=0)
                else:
                    print('Current version of extractor ignores indel'
                          ' to avoid shift in frame')

    def extract_query(self, variant_interval_queryable, sample_id=None):
        cds_seqs = list(self._extract_query(variant_interval_queryable,
                                            sample_id=sample_id))
        return cds_seqs if cds_seqs \
            else self._ref_cds_seq(variant_interval_queryable)

    def extract_cds(self, cds, sample_id=None):
        return next(super().extract_cds(cds, sample_id=sample_id))


class SingleVariantProteinVCFSeqExtractor(ProteinVCFSeqExtractor):

    def extract_query(self, variant_interval_queryable, sample_id=None):
        ref_cds_seq = self._ref_cds_seq(variant_interval_queryable)

        for i, (variants, interval) in enumerate(
                variant_interval_queryable.variant_intervals):

            for variant in variants:
                if len(variant.ref) == len(variant.alt):
                    yield [
                        *ref_cds_seq[:i],
                        self.variant_seq_extractor.extract(
                            interval, [variant], anchor=0),
                        *ref_cds_seq[(i+1):],
                    ]
                else:
                    print('Current version of extractor ignores indel'
                          ' to avoid shift in frame')

def test_mutation_in_each_exon_all_variance():
    gtf_file = 'tests/data/sample_1_protein.gtf'
    fasta_file = 'tests/data/demo_dna_seq.fa'
    vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'
    protein_file = 'tests/data/demo_proteins.pep.all.fa'
    vs = SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
    vs = SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)

    transcript_id = 'ENST00000381176'
    seq = vs.extract(transcript_id)
    txt_file = 'tests/data/Output_singleSeq_vcf_ENST000000381176.txt'
    f = open(txt_file)
    control_seq = f.readline()
    assert seq == control_seq, 'Seq mismatch'


def test_mutation_single_variance():
    gtf_file = 'tests/data/sample_1_protein.gtf'
    fasta_file = 'tests/data/demo_dna_seq.fa'
    vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'
    protein_file = 'tests/data/demo_proteins.pep.all.fa'
    vs = SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)

    transcript_id = 'ENST00000381176'
    seq_list = []
    seq_list_all = vs.extract(transcript_id)
    for seq in seq_list_all:
        seq_list.append(seq)
    assert len(seq_list) == 3, 'Number of seq!=number of variances'
    txt_file = 'tests/data/Output_singleVar_vcf_ENST000000381176.txt'
    f = open(txt_file)
    control_seq_list = []
    for seq in f:
        control_seq_list.append(seq.replace('\n', ''))
    for seq in seq_list:
        assert seq in control_seq_list, 'Seq mismatch'
        control_seq_list.remove(seq)
    assert len(control_seq_list) == 0, '# of expected varSeq!= # of generated varSeq'
    
def test_cut_seq():
    seq = 'ATCGATG'
    seq = cut_seq(seq)
    assert len(seq) == 6, 'cut_seq does not work proper! Expected length 6, but was ' + str(len(seq))

def test_strand_positive():
    gtf_file = 'tests/data/sample_1_protein.gtf'
    fasta_file = 'tests/data/demo_dna_seq.fa'
    vcf_file = 'tests/data/singleVar_vcf_ENST000000381176.vcf.gz'
    protein_file = 'tests/data/demo_proteins.pep.all.fa'
    txt_file = 'tests/data/dna_seq_ENST00000319363.txt'


    'Interval.strand = "+"'
    transcript_id = 'ENST00000319363'
    f = open(txt_file)
    control_seq_list = []
    for seq in f:
        control_seq_list.append(seq.replace('\n', ''))
        
    vs = SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
    test_dna_seq = vs.extract(transcript_id)

    assert test_dna_seq == translate(cut_seq("".join(control_seq_list))), "Seq mismatch for Interval.strand = +"

def test_strand_negative():
    gtf_file = 'tests/data/sample_1_protein.gtf'
    fasta_file = 'tests/data/demo_dna_seq.fa'
    vcf_file = 'tests/data/singleSeq_vcf_ENST000000381176.vcf.gz'
    protein_file = 'tests/data/demo_proteins.pep.all.fa'
    txt_file = 'tests/data/dna_seq_ENST00000381176.txt'

    'Interval.strand = "-"'
    transcript_id = 'ENST00000381176'
    
    f = open(txt_file)
    control_seq_list = []
    for seq in f:
        control_seq_list.append(seq.replace('\n', ''))
        
    vs = SingleSeqProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
    test_dna_seq = vs.extract(transcript_id)

    ref_dna_seq = translate(cut_seq(rc_dna("".join(control_seq_list))))
    assert test_dna_seq == ref_dna_seq, "Seq mismatch for Interval.strand = -"

from pathlib import Path
ddir = Path('/s/genomes/human/hg38/ensembl_GRCh38')
gtf_file = ddir / 'Homo_sapiens.GRCh38.99.gtf'
#gtf_file = 'tests/data/38_sample.gtf.txt'
protein_file = ddir / 'Homo_sapiens.GRCh38.pep.all.fa'

fasta_file = ddir / 'Homo_sapiens.GRCh38.dna.alt.fa'

fasta_file = ddir / 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'

@pytest.mark.skipif(not os.path.isfile(protein_file),
                    reason="No vep result file to test")

def test_all_proteins_translation():

    dfp = read_pep_fa(protein_file)
    dfp['transcript_id'] = dfp.transcript.str.split(".", n=1, expand=True)[0]
    assert not dfp['transcript_id'].duplicated().any()
    dfp = dfp.set_index("transcript_id")
    dfp = dfp[~dfp.chromosome.isnull()]

    gps = TranscriptSeqExtractor(gtf_file, fasta_file)
    assert len(gps) >100
    #assert gps.transcripts.isin(dfp.index).all()

    transcript_id = 'ENST00000485079'
    div3_error = 0
    seq_mismatch_err = 0
    err_transcripts = []
    for transcript_id in tqdm(gps.transcripts):
        # make sure all ids can be found in the proteome
        dna_seq = gps.get_seq(transcript_id)
        # dna_seq = dna_seq[:(len(dna_seq) // 3) * 3]
        #if len(dna_seq) % 3 != 0:
         #   div3_error += 1
          #  print("len(dna_seq) % 3 != 0: {}".format(transcript_id))
           # err_transcripts.append({"transcript_id": transcript_id, "div3_err": True})
            #continue
        prot_seq = dna_seq
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
# err_cds.to_csv("data/protein/err_cds.csv")
# err_transcripts.to_csv("data/protein/err_transcripts.csv")
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