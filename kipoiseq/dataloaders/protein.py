from kipoiseq.extractors import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor
from kipoiseq.extractors import GenericMultiIntervalSeqExtractor, GenericSingleVariantMultiIntervalVCFSeqExtractor, \
    FastaStringExtractor, UTRFetcher, MultiSampleVCF, SingleVariantMatcher
from kipoi.data import SampleIterator

__all__ = [
    'SingleVariantProteinDataLoader',
    'SingleVariantUTRDataLoader',
]


class SingleVariantProteinDataLoader(SampleIterator):

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.protein_vcf_extractor = SingleVariantProteinVCFSeqExtractor(
            gtf_file, fasta_file, vcf_file)
        cds = self.protein_vcf_extractor.cds
        # only needed metadata
        self.metadatas = (
            (
                cds.loc[~cds.index.duplicated(keep='first')]
            ).drop(columns=['Start', 'End'])
        )
        # generator for all sequences with variants
        self.sequences = self._extractor()

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.sequences)

    def _extractor(self):
        """
        return ref_seq, alt_seq, metadata for all
        transcript_ids with variants
        Returns: {
            'input': {
                'ref_seq': ref_seq,
                'alt_seq': alt_seq,
            },
            # dict
            'metadata': metadata
        }
        """
        for transcript_id, (ref_seq, alt_seqs) in self.protein_vcf_extractor.extract_all():
            for (alt_seq, variant) in alt_seqs:
                yield {
                    'input': {
                        'ref_seq': ref_seq,
                        'alt_seq': alt_seq,
                    },
                    'metadata': self.get_metadata(transcript_id, variant)
                }

    def get_metadata(self, transcript_id: str, variant: dict):
        """
        get metadata for given transcript_id
        """
        row = self.metadatas.loc[transcript_id]
        metadata = self.metadatas.loc[transcript_id].to_dict()
        metadata['transcript_id'] = row.name
        metadata['variants'] = variant
        return metadata


class SingleVariantUTRDataLoader(SampleIterator):

    def __init__(
            self,
            gtf_file,
            fasta_file,
            vcf_file,
            feature_type="5UTR",
            infer_from_cds=False,
            on_error_warn=True,
    ):
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.feature_type = feature_type
        self.infer_from_cds = infer_from_cds
        self.on_error_warn = on_error_warn

        self.interval_fetcher = UTRFetcher(
            gtf_file=gtf_file,
            feature_type=feature_type,
            infer_from_cds=infer_from_cds,
            on_error_warn=on_error_warn
        )
        self.multi_sample_VCF = MultiSampleVCF(vcf_file)
        self.reference_seq_extractor = FastaStringExtractor(fasta_file)

        df = self.interval_fetcher.df
        import pyranges
        # match variant with transcript_id
        self.variant_matcher = SingleVariantMatcher(
            self.vcf_file,
            pranges=pyranges.PyRanges(df.reset_index())
        )

        self.extractor = GenericSingleVariantMultiIntervalVCFSeqExtractor(
            interval_fetcher=self.interval_fetcher,
            reference_seq_extractor=self.reference_seq_extractor,
            variant_matcher=self.variant_matcher,
            multi_sample_VCF=self.multi_sample_VCF,
        )

        # only needed metadata
        self.metadatas = (
            (
                df.loc[~df.index.duplicated(keep='first')]
            ).drop(columns=['Start', 'End'])
        )
        # generator for all sequences with variants
        self.sequences = self._extractor()

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.sequences)

    def _extractor(self):
        """
        return ref_seq, alt_seq, metadata for all
        transcript_ids with variants
        Returns: {
            'input': {
                'ref_seq': ref_seq,
                'alt_seq': alt_seq,
            },
            # dict
            'metadata': metadata
        }
        """
        for transcript_id, (ref_seq, alt_seqs) in self.extractor.items():
            for (alt_seq, variant) in alt_seqs:
                yield {
                    'input': {
                        'ref_seq': ref_seq,
                        'alt_seq': alt_seq,
                    },
                    'metadata': self.get_metadata(transcript_id, variant)
                }

    def get_metadata(self, transcript_id: str, variant: dict):
        """
        get metadata for given transcript_id
        """
        row = self.metadatas.loc[transcript_id]
        metadata = self.metadatas.loc[transcript_id].to_dict()
        metadata['transcript_id'] = row.name
        metadata['variants'] = variant
        return metadata
