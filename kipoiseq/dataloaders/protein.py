from kipoiseq.extractors import SingleVariantProteinVCFSeqExtractor, \
    TranscriptSeqExtractor
from kipoi.data import SampleIterator

__all__ = [
    'SingleVariantProteinDataLoader'
]


class SingleVariantProteinDataLoader(SampleIterator):

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.protein_vcf_extractor = SingleVariantProteinVCFSeqExtractor(
            gtf_file, fasta_file, vcf_file)
        self.transcript_extractor = TranscriptSeqExtractor(
            gtf_file, fasta_file)
        cds = self.transcript_extractor.cds
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
        for transcript_id, seqs in self.protein_vcf_extractor.extract_all():
            # reference sequence
            ref_seq = self.transcript_extractor.get_protein_seq(transcript_id)
            for (alt_seq, variant) in seqs:
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
