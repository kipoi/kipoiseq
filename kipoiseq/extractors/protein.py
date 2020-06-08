# from itertools import chain, islice
import numpy as np
from kipoiseq.dataclasses import Interval
from kipoiseq.transforms.functional import rc_dna, translate
from kipoiseq.extractors.base import FastaStringExtractor
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_seq import VariantSeqExtractor
from kipoiseq.extractors.vcf_matching import SingleVariantMatcher
from typing import List

import logging

log = logging.getLogger(__name__)


# TODO: documentation

def cut_transcript_seq(seq: str, tag: str):
    """
    Some of the sequences contain length % 3 != 0, because they have ambiguous
    start and/or end. If this is the case, they should be cut until length % 3 == 0
    There are sequences which have both ambiguous start and end => no solution yet
    :param seq: dna sequences of the current protein
    :param tag: tags, which contain information about ambiguous start and/or end
    :return: correct dna sequence with length % 3 == 0
    if ambiguous start and end or no tags provided, but the sequence has length % 3 != 0
    seq = 'NNN'
    """
    if "cds_end_NF" in tag and "cds_start_NF" not in tag:
        # remove suffix
        seq_modulo = len(seq) % 3
        if seq_modulo != 0:
            seq = seq[:-seq_modulo]

        if seq[-3:] in ["TAA", "TAG", "TGA"]:
            seq = seq[:-3]
    elif "cds_end_NF" not in tag and "cds_start_NF" in tag and len(seq) % 3 != 0:
        # remove prefix
        seq_modulo = len(seq) % 3
        if seq_modulo != 0:
            seq = seq[seq_modulo:]

        seq = "XXX" + seq
    elif "cds_end_NF" in tag and "cds_start_NF" in tag:
        log.warning("Ambiguous start and end! Skip seq!")
        seq = "NNN"  # NNN will be translated as empty string
    elif "cds_end_NF" not in tag and "cds_start_NF" not in tag and len(seq) % 3 != 0:
        log.warning("No tags for ambiguous start and end, but len % 3 != 0. Skip seq!")
        seq = "NNN"  # NNN will be translated as empty string

    return seq


def gtf_row2interval(row, interval_attrs: List[str] = None):
    """
    Convert gtf row object into interval class.
    """
    if interval_attrs:
        interval_attrs = {i: row[i] for i in interval_attrs if i in row}
    return Interval(str(row.Chromosome),
                    int(row.Start),
                    int(row.End),
                    strand=str(row.Strand),
                    attrs=interval_attrs)


def _filter_valid_transcripts(gtf_df):
    if 'transcript_support_level' in gtf_df:
        gtf_df = gtf_df[~gtf_df.transcript_support_level.isnull()]
        gtf_df = gtf_df[gtf_df.transcript_support_level != 'NA']
        gtf_df.transcript_support_level = gtf_df.transcript_support_level.astype(
            int)
        gtf_df = gtf_df[~gtf_df.transcript_support_level.isna()]
    else:
        raise ValueError('Transcript support level not in gtf. '
                         'Skipping the associated filters.')
    return gtf_df


def _get_biotype_str(df):
    if 'transcript_biotype' in df:
        return 'transcript_biotype'
    elif 'gene_biotype' in df:
        return 'gene_biotype'
    else:
        raise ValueError('Cannot obtain `biotype_str` from gtf file')


def _filter_biotype(gtf_df, value):
    """
    Gets the biotype column and checks whether it equals `value`.
    :param gtf_df: genome annotation GTF dataframe
    :param value: The value to check for
    :return: filtered dataframe
    """
    biotype_str = _get_biotype_str(gtf_df)
    return gtf_df.query(
        "{key} == '{value}'".format(key=biotype_str, value=value)
    )


def _filter_biotype_proteincoding(gtf_df):
    return _filter_biotype(gtf_df, value="protein_coding")


def _filter_tag(gtf_df, regex_contains):
    return gtf_df.query(
        "tag.notnull() & tag.str.contains('{}')".format(regex_contains)
    )


class CDSFetcher:

    def __init__(
            self,
            gtf_file,
            filter_valid_transcripts=True,
            filter_biotype=True,
            filter_tag=True,
            on_error_warn=True,
    ):
        """
        Protein sequences in the genome
        :param gtf_file:
        """
        self.gtf_file = str(gtf_file)
        self.cds = self._read_cds(
            self.gtf_file,
            filter_valid_transcripts=filter_valid_transcripts,
            filter_biotype=filter_biotype,
            filter_tag=filter_tag,
            on_error_warn=on_error_warn
        )
        self.transcripts = self.cds.index.unique()

    @staticmethod
    def _read_cds(
            gtf_file,
            filter_valid_transcripts=False,
            filter_biotype=False,
            filter_tag=False,
            duplicate_attr=None,
            on_error_warn=True,
    ):
        """
        Read, extract and filter valid cds from the given gtf_file
        :param gtf_file: path to the GTF file
        """
        import pyranges

        if duplicate_attr == None:
            # One row in the GTF file can have multiple tags;
            # therefore, to filter them we have to allow duplicate attrs.
            duplicate_attr = filter_tag

        df = pyranges.read_gtf(gtf_file, as_df=True, duplicate_attr=duplicate_attr)

        cds = CDSFetcher.get_cds_from_gtf(
            df,
            filter_valid_transcripts=filter_valid_transcripts,
            filter_biotype=filter_biotype,
            filter_tag=filter_tag,
            on_error_warn=on_error_warn
        )

        cds = cds.set_index("transcript_id")
        return cds

    @staticmethod
    def get_cds_from_gtf(
            df,
            filter_valid_transcripts=False,
            filter_biotype=False,
            filter_tag=False,
            on_error_warn=True,
    ):
        """
        Create DataFrame with valid cds
        :param df: the GTF dataframe
        :param filter_valid_transcripts: Filter for "transcript_support_level" column in GTF
        :param filter_biotype: Filter for "[gene|transcript]_biotype == 'protein_coding'" in GTF
        :param filter_tag: Filter for 'basic' or 'CCDS' tag
        """
        cds = df.query("(Feature == 'CDS') | (Feature == 'CCDS')")

        if filter_biotype:
            try:
                cds = _filter_biotype_proteincoding(cds)
            except ValueError as e:
                log.warning("Error during filtering biotype: %s", e)
        if filter_valid_transcripts:
            try:
                cds = _filter_valid_transcripts(cds)
            except ValueError as e:
                log.warning("Error during filtering for valid transcripts: %s", e)
        if filter_tag:
            try:
                cds = _filter_tag(cds, regex_contains='basic|CCDS')
            except ValueError as e:
                log.warning("Error during filtering for tag: %s", e)

        return cds

    def __len__(self):
        return len(self.transcripts)

    def get_cds(self, transcript_id: str):
        """
        Create Interval objects of the cds for given transcript_id
        """
        cds = self.cds.loc[[transcript_id]]
        assert np.all(cds.iloc[0].Strand == cds.Strand)

        return [
            gtf_row2interval(row, ['tag'])
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
    def _prepare_seq(seqs: List[str], strand: str, tag: str):
        """
        Prepare the dna sequence in the final variant, which should be
        translated in amino acid sequence
        :param seqs: current dna sequence
        :param strand: dna strand, where the gene is located
        :param tag: tags, which contain information about ambiguous start/end
        :return: prepared dna sequence ready for translation into amino acid
        sequence
        """
        seq = "".join(seqs)
        if strand == '-':
            # optionally reverse complement
            seq = rc_dna(seq)
        seq = cut_transcript_seq(seq, tag)
        return seq

    def get_seq(self, transcript_id: str):
        """
        Extract the dna sequence for given transcript_id
        and prepare it in its final shape
        :param transcript_id:
        :return: dna sequence for the given transcript_id
        """
        cds = self.cds_fetcher.get_cds(transcript_id)
        seqs = [self.fasta.extract(i) for i in cds]

        return self._prepare_seq(seqs, cds[0].strand, cds[0].attrs["tag"])

    def __getitem__(self, idx):
        return self.get_seq(self.transcripts[idx])

    def overlaped_cds(self, variants):
        """
        19.03.20
        This method is implemented in the ProteinVCFSeqExtractor as
        extract_all()

        The information below is old.

        """
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

    def get_protein_seq(self, transcript_id: str):
        """
        Extract amino acid sequence for given transcript_id
        :param transcript_id: 
        :return: amino acid sequence
        """
        return translate(self.get_seq(transcript_id), hg38=True)


class ProteinSeqExtractor(TranscriptSeqExtractor):

    @staticmethod
    def _prepare_seq(seqs: List[str], strand: str, tag: str):
        """
        Prepare the dna sequence and translate it into amino acid sequence
        :param seqs: current dna sequence
        :param strand: dna strand, where the gene is located
        :param tag: tags, which contain information about ambiguous start/end
        :return: amino acid sequence
        """
        return translate(TranscriptSeqExtractor._prepare_seq(seqs, strand, tag), True)


class ProteinVCFSeqExtractor:

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)
        self.vcf_file = str(vcf_file)
        self.cds_fetcher = CDSFetcher(self.gtf_file)
        # dataframe to pyranges
        import pyranges
        pr_cds = pyranges.PyRanges(self.cds_fetcher.cds.reset_index())
        # match variant with transcript_id
        self.single_variant_matcher = SingleVariantMatcher(
            self.vcf_file, pranges=pr_cds)

        self.fasta = FastaStringExtractor(self.fasta_file)
        self.multi_sample_VCF = MultiSampleVCF(self.vcf_file)
        self.variant_seq_extractor = VariantSeqExtractor(self.fasta_file)

    @staticmethod
    def _unstrand(intervals: List[Interval]):
        """
        Set strand of list of intervals to default - '.'
        """
        return [i.unstrand() for i in intervals]

    @staticmethod
    def _prepare_variants(variants: 'List of variants'):
        variants_dict = dict()
        # fill dict with variants (as dict)
        for index, v in enumerate(variants):
            variants_dict[index] = dict(
                (key.replace('_', ''), value) for key, value in v.__dict__.items())
        # if single varint, unpack dict
        if len(variants_dict) == 1:
            variants_dict = variants_dict[0]
        return variants_dict

    def extract_cds(self, cds: List[Interval], sample_id=None):
        """
        Extract cds with variants in their dna sequence. It depends on the
        child class if a sequence have all variants inserted or only one variant
        is inserted per dna sequence
        :param cds: list of Intervals
        :param sample_id:
        :return: sequence with variants
        """
        intervals = self._unstrand(cds)

        variant_interval_queryable = self.multi_sample_VCF.query_variants(
            intervals, sample_id=sample_id)

        iter_seqs = self.extract_query(variant_interval_queryable,
                                       sample_id=sample_id)
        for seqs in iter_seqs:
            # 1st seq, 2nd variant info
            yield ProteinSeqExtractor._prepare_seq(seqs[0], cds[0].strand, cds[0].attrs['tag']), seqs[1]

    def extract_all(self):
        """
        Extract all amino acid sequences for transcript_ids with variants
        given into the vcf_file
        """
        for pr in self.single_variant_matcher.iter_pyranges():
            # check if variants exist
            if len(pr) > 0:
                for transcript_id in pr.df.transcript_id.drop_duplicates():
                    yield transcript_id, self.extract(transcript_id)
            else:
                log.warning('No matched variants with transcript_ids.')

    def extract_list(self, list_with_transcript_id: List[str]):
        """
        Extract all amino acid sequences for transcript_id given in the list
        :param list_with_transcript_id: list which contains transcript_ids
        :return: sequences with variants
        """
        for transcript_id in list_with_transcript_id:
            yield transcript_id, self.extract(transcript_id)

    def extract(self, transcript_id, sample_id=None):
        """
        Extract all amino acid sequences for transcript_id
        """
        return self.extract_cds(self.cds_fetcher.get_cds(transcript_id),
                                sample_id=sample_id)

    def _ref_cds_seq(self, variant_interval_queryable):
        """
        Extract amino acid sequence without variants, which can be used as
        reference
        """
        intervals = variant_interval_queryable.iter_intervals()
        return [self.fasta.extract(interval) for interval in intervals]

    def _filter_snv(self, variants):
        for variant in variants:
            if len(variant.ref) == len(variant.alt) == 1:  # only SOVs supported
                yield variant
            elif len(variant.ref) == len(variant.alt) > 1:
                log.warning('Current version of extractor works only for len(variant.ref)'
                            ' == len(variant.alt) == 1, but the len was: ' + str(len(variant.alt)))
            else:
                log.warning('Current version of extractor ignores indel'
                            ' to avoid shift in frame')


class SingleSeqProteinVCFSeqExtractor(ProteinVCFSeqExtractor):

    def _extract_query(self, variant_interval_queryable, sample_id=None):
        """
        Iterate through all intervals and extract dna sequence with all
        variants inserted into it
        :param variant_interval_queryable: Object which contains information
        about the variants for current sequence
        :param sample_id:
        :return: dna sequence with all variants
        """
        seqs = []
        flag = True
        variants_info = list()
        for variants, interval in variant_interval_queryable.variant_intervals:
            variants = list(self._filter_snv(variants))
            if len(variants) > 0:
                flag = False
                variants_info.extend(variants)
            seqs.append(self.variant_seq_extractor.extract(
                interval, variants, anchor=0))
        if flag:
            seqs = []
            variants_info = []
        yield "".join(seqs), self._prepare_variants(variants_info)

    def extract_query(self, variant_interval_queryable, sample_id=None):
        """
        Extract dna sequence with all variants inserted
        """
        cds_seqs = list(self._extract_query(variant_interval_queryable,
                                            sample_id=sample_id))
        return cds_seqs

    def extract_cds(self, cds: 'list of Intervals', sample_id=None):
        """
        Call parent method which inserts all variants into the dna sequence
        """
        return next(super().extract_cds(cds, sample_id=sample_id))


class SingleVariantProteinVCFSeqExtractor(ProteinVCFSeqExtractor):

    def extract_query(self, variant_interval_queryable, sample_id=None):
        """
        Iterate through all variants and creates a sequence for
        each variant individually
        :param variant_interval_queryable: Object which contains information
        about the variants for current sequence
        :param sample_id:
        :return: for each variant a sequence with a single variant
        """
        ref_cds_seq = self._ref_cds_seq(variant_interval_queryable)
        for i, (variants, interval) in enumerate(
                variant_interval_queryable.variant_intervals):
            variants = self._filter_snv(variants)
            for variant in variants:
                yield [
                          *ref_cds_seq[:i],
                          self.variant_seq_extractor.extract(
                              interval, [variant], anchor=0),
                          *ref_cds_seq[(i + 1):],
                      ], self._prepare_variants([variant])
