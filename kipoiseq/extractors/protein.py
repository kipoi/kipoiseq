# from itertools import chain, islice
import abc

from kipoiseq.extractors import CDSFetcher, UTRFetcher

from kipoiseq.dataclasses import Interval, Variant
from kipoiseq.transforms.functional import translate
from kipoiseq.extractors.multi_interval import (
    GenericMultiIntervalSeqExtractor,
    BaseMultiIntervalVCFSeqExtractor,
    SingleVariantExtractorMixin,
    SingleSeqExtractorMixin,
)
from kipoiseq.extractors.fasta import FastaStringExtractor
from kipoiseq.extractors.vcf import MultiSampleVCF
from kipoiseq.extractors.vcf_matching import SingleVariantMatcher
from typing import List, Union

import logging

log = logging.getLogger(__name__)

__all__ = [
    "cut_transcript_seq",
    "UTRSeqExtractor",
    "TranscriptSeqExtractor",
    "ProteinSeqExtractor",
    "TranscriptVCFSeqExtractor",
    "ProteinVCFSeqExtractor",
    "SingleSeqProteinVCFSeqExtractor",
    "SingleVariantProteinVCFSeqExtractor",
]


class UTRSeqExtractor(GenericMultiIntervalSeqExtractor):
    def __init__(
            self,
            gtf_file,
            fasta_file,
            feature_type="5UTR",
            infer_from_cds=False,
            on_error_warn=True,
    ):
        """
        Reference sequence extractor for UTR's

        :param fasta_file: fasta file for reference sequence input
        :param gtf_file: path to the GTF file
        :param feature_type: type of the feature that will be filtered for. In general '5UTR' or '3UTR'.
        :param infer_from_cds: Substract the CDS from the exon regions to infer the UTR regions.
            Will use 'feature_type' to decide whether '5UTR' or '3UTR' should be returned.
        :param on_error_warn: print warning instead of throwing an error
        """
        self.fasta_file = str(fasta_file)
        self.gtf_file = str(gtf_file)

        utr_fetcher = UTRFetcher(self.gtf_file, feature_type=feature_type, infer_from_cds=infer_from_cds,
                                 on_error_warn=on_error_warn)

        extractor = FastaStringExtractor(self.fasta_file, use_strand=False)

        super().__init__(
            extractor=extractor,
            interval_fetcher=utr_fetcher
        )

    @property
    def df(self):
        return self.extractor.df


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


# TODO: documentation
class TranscriptSeqExtractor(GenericMultiIntervalSeqExtractor):

    def __init__(self, gtf_file, fasta_file):
        self.fasta_file = str(fasta_file)
        self.gtf_file = str(gtf_file)

        cds_fetcher = CDSFetcher(self.gtf_file)

        extractor = FastaStringExtractor(self.fasta_file, use_strand=False)

        super().__init__(
            extractor=extractor,
            interval_fetcher=cds_fetcher
        )

    @property
    def df(self):
        return self.extractor.df

    @property
    def cds(self):
        return self.extractor.df

    @classmethod
    def _prepare_seq(
            cls,
            seqs: List[str],
            intervals: List[Interval],
            reverse_complement: Union[str, bool],
            # **kwargs
    ) -> str:
        """
        Prepare the dna sequence in the final variant, which should be
        translated in amino acid sequence
        :param seqs: current dna sequence
        :param intervals: the list of intervals corresponding to the sequence snippets
        :param reverse_complement: should the dna be reverse-complemented?
        :return: prepared dna sequence ready for translation into amino acid sequence
        """
        seq = super()._prepare_seq(
            seqs=seqs,
            intervals=intervals,
            reverse_complement=reverse_complement,
        )
        tag = intervals[0].attrs["tag"]
        seq = cut_transcript_seq(seq, tag)
        return seq

    # def extract(self, intervals: List[Interval], **kwargs) -> str:
    #     """
    #     Extract and concatenate the sequence for a list of intervals
    #     :param intervals: list of intervals
    #     :param kwargs: will be passed on to `_prepare_seq()`
    #     :return: concatenated dna sequence of the intervals
    #     """
    #     seqs = [self.extractor.extract(i) for i in intervals]
    #
    #     reverse_strand = False
    #     if self.use_strand:
    #         reverse_strand = intervals[0].neg_strand
    #         if self.extractor.use_strand:
    #             # If the fasta extractor already does the reversion, we do not have to do it manually.
    #             # Instead, we need to reverse the list of sequences.
    #             seqs.reverse()
    #             reverse_strand = False
    #
    #     tag = intervals[0].attrs["tag"]
    #
    #     return self._prepare_seq(seqs, reverse_strand, tag=tag, **kwargs)

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
        :param intervals: the list of intervals corresponding to the sequence snippets
        :param reverse_complement: should the dna be reverse-complemented?
        :return: amino acid sequence
        """
        return translate(super()._prepare_seq(*args, **kwargs), hg38=True)


class TranscriptVCFSeqExtractor(BaseMultiIntervalVCFSeqExtractor, metaclass=abc.ABCMeta):
    @classmethod
    def _prepare_seq(
            cls,
            seqs: List[str],
            intervals: List[Interval],
            reverse_complement: Union[str, bool],
            # **kwargs
    ) -> str:
        """
        Prepare the dna sequence in the final variant, which should be
        translated in amino acid sequence
        :param seqs: current dna sequence
        :param intervals: the list of intervals corresponding to the sequence snippets
        :param reverse_complement: should the dna be reverse-complemented?
        :return: prepared dna sequence ready for translation into amino acid sequence
        """
        seq = super()._prepare_seq(
            seqs=seqs,
            intervals=intervals,
            reverse_complement=reverse_complement,
        )
        tag = intervals[0].attrs["tag"]
        seq = cut_transcript_seq(seq, tag)
        return seq


class ProteinVCFSeqExtractor(TranscriptVCFSeqExtractor, metaclass=abc.ABCMeta):

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.gtf_file = str(gtf_file)
        self.fasta_file = str(fasta_file)
        self.vcf_file = str(vcf_file)

        cds_fetcher = CDSFetcher(self.gtf_file)
        self.cds = cds_fetcher.df

        reference_seq = FastaStringExtractor(self.fasta_file)
        multi_sample_VCF = MultiSampleVCF(self.vcf_file)

        # dataframe to pyranges
        import pyranges
        # match variant with transcript_id
        variant_matcher = SingleVariantMatcher(
            self.vcf_file,
            pranges=pyranges.PyRanges(cds_fetcher.df.reset_index())
        )

        super().__init__(
            interval_fetcher=cds_fetcher,
            reference_seq_extractor=reference_seq,
            variant_matcher=variant_matcher,
            multi_sample_VCF=multi_sample_VCF,
        )

    @classmethod
    def _prepare_seq(cls, *args, **kwargs):
        """
        Prepare the dna sequence and translate it into amino acid sequence
        :param seqs: current dna sequence
        :param intervals: the list of intervals corresponding to the sequence snippets
        :param reverse_complement: should the dna be reverse-complemented?
        :return: amino acid sequence
        """
        return translate(super()._prepare_seq(*args, **kwargs), hg38=True)

    # def extract(
    #         self,
    #         intervals: List[Interval],
    #         sample_id: List[str] = None,
    #         **kwargs
    # ):
    #     reverse_complement = intervals[0].neg_strand
    #     tag = intervals[0].attrs["tag"]
    #     # remove strand information
    #     intervals = [i.unstrand() for i in intervals]
    #
    #     variant_interval_queryable = self.multi_sample_VCF.query_variants(intervals, sample_id=sample_id)
    #
    #     iter_seqs = self.extract_query(variant_interval_queryable, sample_id=sample_id)
    #     for seqs, variant_info in iter_seqs:
    #         # 1st seq, 2nd variant info
    #         yield ProteinSeqExtractor._prepare_seq(seqs, reverse_complement, tag=tag), variant_info

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


class SingleSeqProteinVCFSeqExtractor(SingleSeqExtractorMixin, ProteinVCFSeqExtractor):
    pass
    # def extract_query(
    #         self,
    #         variant_interval_queryable,
    #         ref_seqs,
    #         reverse_complement,
    #         sample_id=None
    # ):
    #     """
    #     Iterate through all intervals and extract dna sequence with all variants inserted into it
    #     :return: dna sequence with all variants. If no variants match, will return the reference sequence.
    #     """
    #     """
    #             Iterate through all intervals and extract dna sequence with all
    #             variants inserted into it
    #             :param variant_interval_queryable: Object which contains information
    #             about the variants for current sequence
    #             :param sample_id:
    #             :return: dna sequence with all variants. If no variants match, will return the reference sequence.
    #             """
    #     seqs = []
    #     variants_info = []
    #     for variants, interval in variant_interval_queryable.variant_intervals:
    #         variants = list(self._filter_snv(variants))
    #         if len(variants) > 0:
    #             flag = False
    #             variants_info.extend(variants)
    #         seqs.append(
    #             self.variant_seq_extractor.extract(interval, variants, anchor=0)
    #         )
    #
    #     alt_seq = "".join(seqs)
    #     variants_info = self._prepare_variants(variants_info)
    #     return (
    #         alt_seq,  # the final sequence
    #         variants_info,  # dictionary of variants
    #     )
    #
    #     # cds_seqs = list(self._extract_query(variant_interval_queryable, sample_id=sample_id))
    #     # return cds_seqs


class SingleVariantProteinVCFSeqExtractor(SingleVariantExtractorMixin, ProteinVCFSeqExtractor):
    pass
