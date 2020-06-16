# from itertools import chain, islice
import abc

import pandas as pd
from kipoiseq.extractors.gtf import CDSFetcher, gtf_row2interval
from kipoiseq.extractors.vcf_query import VariantIntervalQueryable

from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.transforms.functional import rc_dna, translate
from kipoiseq.extractors.base import BaseExtractor, BaseMultiIntervalSeqExtractor, GenericMultiIntervalSeqExtractor
from kipoiseq.extractors.fasta import FastaStringExtractor
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_seq import VariantSeqExtractor
from kipoiseq.extractors.vcf_matching import SingleVariantMatcher, BaseVariantMatcher
from typing import List, Tuple, Union

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
    # if not tag:
    #     return seq

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


class TranscriptSeqExtractor(GenericMultiIntervalSeqExtractor):

    def __init__(self, gtf_file, fasta_file):
        self.fasta_file = str(fasta_file)
        self.gtf_file = str(gtf_file)

        cds_fetcher = CDSFetcher(self.gtf_file)
        self.cds = cds_fetcher.df

        # self.transcripts = self.cds_fetcher.keys()
        extractor = FastaStringExtractor(self.fasta_file, use_strand=False)

        super().__init__(
            extractor=extractor,
            interval_fetcher=cds_fetcher
        )

    @classmethod
    def _prepare_seq(
            cls,
            seqs: List[str],
            reverse_complement: Union[str, bool],
            tag=None,
            **kwargs
    ) -> str:
        """
        Prepare the dna sequence in the final variant, which should be
        translated in amino acid sequence
        :param seqs: current dna sequence
        :param reverse_complement: should the dna be reverse-complemented?
        :param tag: string containing tag list
        :return: prepared dna sequence ready for translation into amino acid sequence
        """
        seq = super()._prepare_seq(
            seqs=seqs,
            # intervals=intervals,
            reverse_complement=reverse_complement,
        )
        # tag = intervals[0].attrs["tag"]
        seq = cut_transcript_seq(seq, tag)
        return seq

    def get_seq(self, transcript_id: str):
        """
        Extract the dna sequence for given transcript_id
        and prepare it in its final shape
        :param transcript_id:
        :return: dna sequence for the given transcript_id
        """
        return self.sel(transcript_id)
        # cds_intervals = self.cds_fetcher.get_intervals(transcript_id)
        # return self.extract(cds_intervals, tag=cds_intervals[0].attrs["tag"])

        # return self._prepare_seq(seqs, cds_intervals[0].strand, cds_intervals[0].attrs["tag"])

    def extract(self, intervals: List[Interval], **kwargs) -> str:
        """
        Extract and concatenate the sequence for a list of intervals
        :param intervals: list of intervals
        :param kwargs: will be passed on to `_prepare_seq()`
        :return: concatenated dna sequence of the intervals
        """
        seqs = [self.extractor.extract(i) for i in intervals]

        reverse_strand = False
        if self.use_strand:
            reverse_strand = intervals[0].neg_strand
            if self.extractor.use_strand:
                # If the fasta extractor already does the reversion, we do not have to do it manually.
                # Instead, we need to reverse the list of sequences.
                seqs.reverse()
                reverse_strand = False

        tag = intervals[0].attrs["tag"]

        return self._prepare_seq(seqs, reverse_strand, tag=tag, **kwargs)

    # def __getitem__(self, idx):
    #     return self.get_seq(self.transcripts[idx])

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

    @classmethod
    def _prepare_seq(cls, *args, **kwargs):
        """
        Prepare the dna sequence and translate it into amino acid sequence
        :param seqs: current dna sequence
        :param reverse_complement: should the dna be reverse-complemented?
        :param tag: string containing tag list
        :return: amino acid sequence
        """
        return translate(super()._prepare_seq(*args, **kwargs), True)


class VCFSeqExtractor(metaclass=abc.ABCMeta):

    def __init__(
            self,
            ranges_df: pd.DataFrame,
            reference_seq_extractor: BaseExtractor,
            variant_matcher: BaseVariantMatcher,
            multi_sample_VCF: MultiSampleVCF,
    ):
        self.ranges_df = ranges_df
        self.reference_seq = reference_seq_extractor
        self.variant_matcher = variant_matcher
        self.multi_sample_VCF = multi_sample_VCF

        # takes reference sequence, interval + variants to extract the alternative sequence
        self.variant_seq_extractor = VariantSeqExtractor(reference_sequence=self.reference_seq)

    @classmethod
    def from_pyranges(
            cls,
            pyranges,
            concat_by: Union[str, List[str]] = "transcript_id",
            **kwargs
    ):
        return cls(ranges_df=pyranges.cds.set_index(concat_by), **kwargs)

    def keys(self):
        return self.ranges_df.index.drop_duplicates()

    def extract_query(self, variant_interval_queryable: VariantIntervalQueryable, sample_id=None):
        raise NotImplementedError()

    @staticmethod
    def _unstrand(intervals: List[Interval]):
        """
        Set strand of list of intervals to default - '.'
        """
        return [i.unstrand() for i in intervals]

    @classmethod
    def _prepare_seq(
            cls,
            seqs: List[str],
            reverse_complement: Union[str, bool],
            **kwargs
    ) -> str:
        """
        Prepare the dna sequence in the final variant by concatenating and eventually reverse-complementing it

        :param seqs: current dna sequence
        :param reverse_complement: if True, the sequence will be reversed
        """
        seq = "".join(seqs)
        if reverse_complement:
            # optionally reverse complement
            seq = rc_dna(seq)
        return seq

    def extract_sequence(self, intervals_df: pd.DataFrame, sample_id=None) -> Tuple[str, object]:
        """
        Extract cds with variants in their dna sequence. It depends on the
        child class if a sequence have all variants inserted or only one variant
        is inserted per dna sequence
        :param intervals_df: GTF dataframe as given by PyRanges containing at least the following columns:
            - Chromosome
            - Start
            - End
            - Strand
        :param sample_id: optional sample id
        :return: alternative sequence with variant information
        """
        reverse_complement = intervals_df.iloc[0].Strand == "-"
        intervals = [
            gtf_row2interval(row) for i, row in intervals_df.sort_values("Start").iterrows()
        ]
        intervals = self._unstrand(intervals)

        variant_interval_queryable = self.multi_sample_VCF.query_variants(intervals, sample_id=sample_id)

        # extract_query can return one sequence per variant or one sequence per sample
        iter_seqs = self.extract_query(variant_interval_queryable, sample_id=sample_id)
        for seqs, variant_info in iter_seqs:
            # 1st seq, 2nd variant info
            yield self._prepare_seq(seqs, intervals, reverse_complement), variant_info

        return None

    def extract_all(self):
        """
        Extract all amino acid sequences for transcript_ids with variants
        given into the vcf_file
        """
        yield from self.extract_multiple(self.keys())

    def extract_multiple(self, keys: List[str]):
        """
        Extract all amino acid sequences for transcript_id given in the list
        :param keys: list which contains transcript_ids
        :return: sequences with variants
        """
        for key in keys:
            retval = self.extract(key)
            if retval is not None:
                yield key, retval

    def extract(self, key, sample_id=None):
        """
        Extract all amino acid sequences for transcript_id
        """
        return self.extract_sequence(self.ranges_df.loc[[key]], sample_id=sample_id)

    def _reference_sequence(self, variant_interval_queryable):
        """
        Extract reference sequence without variants
        """
        intervals = variant_interval_queryable.iter_intervals()
        return [self.reference_seq.extract(interval) for interval in intervals]

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


class ProteinVCFSeqExtractor(VCFSeqExtractor, metaclass=abc.ABCMeta):

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)
        self.vcf_file = str(vcf_file)

        reference_seq = FastaStringExtractor(self.fasta_file)
        multi_sample_VCF = MultiSampleVCF(self.vcf_file)

        ranges_df = CDSFetcher(self.gtf_file).df

        # dataframe to pyranges
        import pyranges
        # match variant with transcript_id
        variant_matcher = SingleVariantMatcher(
            self.vcf_file,
            pranges=pyranges.PyRanges(ranges_df.reset_index())
        )

        super().__init__(
            ranges_df=ranges_df,
            reference_seq_extractor=reference_seq,
            variant_matcher=variant_matcher,
            multi_sample_VCF=multi_sample_VCF,
        )

    def extract_query(self, variant_interval_queryable: VariantIntervalQueryable, sample_id=None):
        raise NotImplementedError()

    @staticmethod
    def _prepare_variants(variants: List[Variant]):
        variants_dict = dict()
        # fill dict with variants (as dict)
        for index, v in enumerate(variants):
            variants_dict[index] = dict(
                (key.replace('_', ''), value) for key, value in v.__dict__.items()
            )
        # if single varint, unpack dict
        if len(variants_dict) == 1:
            variants_dict = variants_dict[0]
        return variants_dict

    def extract_sequence(self, intervals_df: pd.DataFrame, sample_id=None):
        """
        Extract cds with variants in their dna sequence. It depends on the
        child class if a sequence have all variants inserted or only one variant
        is inserted per dna sequence
        :param cds: list of Intervals
        :param sample_id:
        :return: sequence with variants
        """
        # strand = intervals_df.iloc[0].Strand
        # tags = intervals_df.iloc[0].tag
        intervals = [
            gtf_row2interval(row, ["tag"]) for i, row in intervals_df.sort_values("Start").iterrows()
        ]

        yield from self.extract_sequence_by_intervals(intervals)

    def extract_sequence_by_intervals(
            self,
            intervals: List[Interval],
            sample_id: List[str] = None,
            **kwargs
    ):
        reverse_complement = intervals[0].neg_strand
        tag = intervals[0].attrs["tag"]
        # remove strand information
        intervals = [i.unstrand() for i in intervals]

        variant_interval_queryable = self.multi_sample_VCF.query_variants(intervals, sample_id=sample_id)

        iter_seqs = self.extract_query(variant_interval_queryable, sample_id=sample_id)
        for seqs, variant_info in iter_seqs:
            # 1st seq, 2nd variant info
            yield ProteinSeqExtractor._prepare_seq(seqs, reverse_complement, tag=tag), variant_info

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
        :return: dna sequence with all variants. If no variants match, will return the reference sequence.
        """
        seqs = []
        variants_info = []
        for variants, interval in variant_interval_queryable.variant_intervals:
            variants = list(self._filter_snv(variants))
            if len(variants) > 0:
                flag = False
                variants_info.extend(variants)
            seqs.append(
                self.variant_seq_extractor.extract(interval, variants, anchor=0)
            )
        yield "".join(seqs), self._prepare_variants(variants_info)

    def extract_query(self, variant_interval_queryable, sample_id=None):
        """
        Extract dna sequence with all variants inserted
        """
        cds_seqs = list(self._extract_query(variant_interval_queryable, sample_id=sample_id))
        return cds_seqs

    def extract_sequence(self, intervals_df: pd.DataFrame, sample_id=None):
        """
        Call parent method which inserts all variants into the dna sequence
        """
        return next(super().extract_sequence(intervals_df, sample_id=sample_id))


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
        ref_cds_seq = self._reference_sequence(variant_interval_queryable)
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
