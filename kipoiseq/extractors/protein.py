from itertools import chain, islice
import numpy as np
from kipoiseq.dataclasses import Interval
from kipoiseq.transforms.functional import rc_dna, translate
from kipoiseq.extractors.base import FastaStringExtractor
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_seq import VariantSeqExtractor

# TODO: convert print to logs
# TODO: documentation


def cut_transcript_seq(seq, tag):
    # if the dna seq is not %3==0, there are unnecessary bases at the end
    # should be called only after all exons are connected!
    if "cds_end_NF" in tag and "cds_start_NF" not in tag:
        while(len(seq) % 3 != 0):
            seq = seq[:-1]
        if seq[-3:] in ["TAA", "TAG", "TGA"]:
            seq = seq[:-3]
    elif "cds_end_NF" not in tag and "cds_start_NF" in tag and len(seq) % 3 != 0:
        while(len(seq) % 3 != 0):
            seq = seq[1:]
        seq = "XXX"+seq
    elif "cds_end_NF" in tag and "cds_start_NF" in tag:
        print("Ambiguous start and end!")
        seq = "XXX"
    elif "cds_end_NF" not in tag and "cds_start_NF" not in tag and len(seq) % 3 != 0:
        print("No tags for ambiguous start and end, but len % 3 != 0")
        seq = "XXX"
        print(tag)
        
    return seq


def gtf_row2interval(row):
    """
    Convert gtf row object into interval class.
    """
    return Interval(str(row.Chromosome),
                    int(row.Start),
                    int(row.End),
                    strand=str(row.Strand),
                   attrs={'tag': str(row.tag)})


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
        df = pyranges.read_gtf(gtf_file, output_df=True, duplicate_attr=True)

        cds = CDSFetcher._get_cds_from_gtf(df)
        cds = CDSFetcher._filter_valid_transcripts(cds)
        return cds

    @staticmethod
    def _get_cds_from_gtf(df):
        biotype_str = CDSFetcher._get_biotype_str(df)
        df = (df
                .query("{} == 'protein_coding'".format(biotype_str))
                .query("(Feature == 'CDS') | (Feature == 'CCDS')")
                )
        return df[df["tag"].str.contains("basic|CCDS")].set_index('transcript_id')

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

    @staticmethod
    def _prepare_seq(seqs, strand, tag):
        seq = "".join(seqs)
        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        seq = cut_transcript_seq(seq, tag)
        return seq

    def get_seq(self, transcript_id):
        cds = self.cds_fetcher.get_cds(transcript_id)
        seqs = [self.fasta.extract(i) for i in cds]
        return self._prepare_seq(seqs, cds[0].strand, cds[0].attrs["tag"])

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
    def _prepare_seq(seqs, strand, tag):
        return translate(TranscriptSeqExtractor._prepare_seq(seqs, strand, tag))


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
            yield ProteinSeqExtractor._prepare_seq(seqs, cds[0].strand, cds[0].attrs['tag'])

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
