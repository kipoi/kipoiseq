from kipoiseq.extractors import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor
from kipoiseq.extractors import GenericMultiIntervalSeqExtractor, GenericSingleVariantMultiIntervalVCFSeqExtractor, \
    FastaStringExtractor, UTRFetcher, MultiSampleVCF, SingleVariantMatcher

from kipoi.data import SampleIterator, kipoi_dataloader
from kipoi_conda.dependencies import Dependencies
from kipoi.specs import Author

__all__ = [
    'SingleVariantProteinDataLoader',
    'SingleVariantUTRDataLoader',
]
deps = Dependencies(
    conda=[
        'bioconda::pybedtools',
        'bioconda::pyfaidx',
        'bioconda::pyranges',
        'bioconda::biopython',
        'numpy',
        'pandas',
    ],
    pip=['kipoiseq']
)
package_authors = [
    Author(name='Florian R. Hölzlwimmer', github='hoeze'),
    Author(name='Kalin Nonchev', github='KalinNonchev')
]


class SingleVariantProteinDataLoader(SampleIterator):
    """
    info:
        doc: >
            Dataloader for protein sequence models. With inputs as gtf annotation file and fasta file,
            each output is a protein sequence with flanking intronic seuqences. Intronic sequnce
            lengths specified by the users. Returned sequences are of the type np.array([str])
    type: SampleIterator
    args:
        gtf_file:
            doc: file path; Genome annotation GTF file
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22_ENST00000319363.gtf?raw=true
                md5: 8a1f158e17379773fcab21628fc3910f
        fasta_file:
            doc: Reference Genome sequence in fasta format
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22.fa.gz?raw=true
                md5: 5ebe034256ecc5689989a96387c5a65e
        vcf_file:
            doc: Genomic variants to evaluate in VCF format
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22_ENST00000319363.vcf.gz?raw=true
                md5: c45e75fb75326c2be514d2dcea52e585
    output_schema:
        inputs:
            ref_seq:
                name: ref_seq
                shape: ()
                special_type: DNAStringSeq
                doc: reference sequence of UTR
                associated_metadata: ranges
            alt_seq:
                name: alt_seq
                doc: alternative sequence of 5' UTR
                shape: ()
                special_type: DNAStringSeq
                associated_metadata: ranges, variants
        metadata:
            transcript_id:
                type: str
                doc: transcript id
            variant:
                CHROM:
                    type: str
                    doc: chromsome of variant
                POS:
                    type: int
                    doc: variant position
                REF:
                    type: str
                    doc: variant reference
                ALT:
                    type: str
                    doc: variant alternative string
                STR:
                    type: str
                    doc: string representation of the variant
    """

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.protein_vcf_extractor = SingleVariantProteinVCFSeqExtractor(
            gtf_file, fasta_file, vcf_file)

        # # only needed metadata
        # cds = self.protein_vcf_extractor.cds
        # self.metadatas = (
        #     (
        #         cds.loc[~cds.index.duplicated(keep='first')]
        #     ).drop(columns=['Start', 'End'])
        # )

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
            'inputs': {
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
                    'inputs': {
                        'ref_seq': ref_seq,
                        'alt_seq': alt_seq,
                    },
                    'metadata': self.get_metadata(transcript_id, variant)
                }

    # def get_metadata(self, transcript_id: str, variant: dict):
    #     """
    #     get metadata for given transcript_id
    #     """
    #     row = self.metadatas.loc[transcript_id]
    #     metadata = self.metadatas.loc[transcript_id].to_dict()
    #     metadata['transcript_id'] = row.name
    #     metadata['variants'] = variant
    #     return metadata

    def get_metadata(self, transcript_id: str, variant: dict):
        """
        get metadata for given transcript_id
        """
        metadata = dict()
        metadata['transcript_id'] = transcript_id
        variant_str_repr = ":".join([
            variant["chrom"],
            str(variant["pos"]),
            variant["ref"],
            variant["alt"],
        ])
        metadata['variant'] = {
            "chrom": variant["chrom"],
            "pos": variant["pos"],
            "ref": variant["ref"],
            "alt": variant["alt"],
            "id": variant['id'] if "id" in variant else variant_str_repr,
            "str": variant_str_repr
        }
        return metadata


@kipoi_dataloader(override={"dependencies": deps, 'info.authors': package_authors})
class SingleVariantUTRDataLoader(SampleIterator):
    """
    info:
        doc: >
            Dataloader for splicing models. With inputs as gtf annotation file and fasta file,
            each output is an exon sequence with flanking intronic seuqences. Intronic sequnce
            lengths specified by the users. Returned sequences are of the type np.array([str])
    type: SampleIterator
    args:
        gtf_file:
            doc: file path; Genome annotation GTF file
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22_ENST00000319363.gtf?raw=true
                md5: 8a1f158e17379773fcab21628fc3910f
                name: gtf_file.gtf
        fasta_file:
            doc: Reference Genome sequence in fasta format
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22.fa.gz?raw=true
                md5: 5ebe034256ecc5689989a96387c5a65e
                name: fasta_file.fa.gz
        vcf_file:
            doc: Genomic variants to evaluate in VCF format
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22_ENST00000319363.vcf.gz?raw=true
                md5: c45e75fb75326c2be514d2dcea52e585
                name: vcf_file.vcf.gz
        vcf_file_tbi:
            doc: tabix index of vcf (just to make kipoi tests work - leave as None in normal usage)
            example:
                url: https://github.com/kipoi/kipoiseq/blob/ddeb4eefc15ebf8a9b88fca4ce99d9b315d54f34/tests/data/chr22_ENST00000319363.vcf.gz.tbi?raw=true
                md5: 9aebc88287a3d6b8517ace9e0fc427af
                name: vcf_file.vcf.gz.tbi
        feature_type:
            doc: Either 5UTR or 3UTR
            example: 5UTR
            type: str
        infer_from_cds:
            doc: infer UTR regions from coding sequence
            optional: True
            default: False
            example: False
            type: bool
        on_error_warn:
            doc: print warning instead of throwing an error on malformed input
            optional: True
            default: True
            example: True
            type: bool
    output_schema:
        inputs:
            ref_seq:
                name: ref_seq
                shape: ()
                special_type: DNAStringSeq
                doc: reference sequence of UTR
                associated_metadata: ranges
            alt_seq:
                name: alt_seq
                doc: alternative sequence of 5' UTR
                shape: ()
                special_type: DNAStringSeq
                associated_metadata: ranges, variants
        metadata:
            transcript_id:
                type: str
                doc: transcript id
            variant:
                chrom:
                    type: str
                    doc: chromsome of variant
                pos:
                    type: int
                    doc: variant position
                ref:
                    type: str
                    doc: variant reference
                alt:
                    type: str
                    doc: variant alternative string
                id:
                    type: str
                    doc: variant id
                str:
                    type: str
                    doc: string representation of the variant
    """

    def __init__(
            self,
            gtf_file,
            fasta_file,
            vcf_file,
            feature_type,
            infer_from_cds=False,
            on_error_warn=True,
            vcf_file_tbi=None,
            **kwargs
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

        # # only needed metadata
        # self.metadatas = (
        #     (
        #         df.loc[~df.index.duplicated(keep='first')]
        #     ).drop(columns=['Start', 'End'])
        # )
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
                    'inputs': {
                        'ref_seq': ref_seq,
                        'alt_seq': alt_seq,
                    },
                    'metadata': self.get_metadata(transcript_id, variant)
                }

    def get_metadata(self, transcript_id: str, variant: dict):
        """
        get metadata for given transcript_id
        """
        metadata = dict()
        metadata['transcript_id'] = transcript_id
        variant_str_repr = ":".join([
            variant["chrom"],
            str(variant["pos"]),
            variant["ref"],
            variant["alt"],
        ])
        metadata['variant'] = {
            "chrom": variant["chrom"],
            "pos": variant["pos"],
            "ref": variant["ref"],
            "alt": variant["alt"],
            "id": variant['id'] if "id" in variant else variant_str_repr,
            "str": variant_str_repr
        }
        return metadata
