from kipoiseq.extractors import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor
from kipoi.data import SampleIterator, kipoi_dataloader
from kipoi.specs import Author, Dependencies

package_authors = [Author(name='Jun Cheng', github='s6juncheng')]  # TODO check if this is correct
deps = Dependencies(conda=['bioconda::pyfaidx', 'numpy', 'pandas'],
                    pip=['kipoiseq', 'kipoi'])


#@kipoi_dataloader(override={"dependencies": deps, 'info.authors': package_authors}) # TODO check if this is correct
class SingleVariantProteinDataLoader(SampleIterator):

    def __init__(self, gtf_file, fasta_file, vcf_file):
        self.single_variant_protein_VCF_seq_extractor = SingleVariantProteinVCFSeqExtractor(gtf_file, fasta_file, vcf_file)
        self.transcript_seq_extractor = TranscriptSeqExtractor(gtf_file, fasta_file)
        cds = self.transcript_seq_extractor.cds_fetcher.cds
        # only needed metadata
        self.metadatas = (cds.loc[~cds.index.duplicated(keep='first')]).drop(columns=['Start', 'End'])
        self.sequences = self._extractor()

    def __iter__(self):
        return self
    
    def _extractor(self):
        for transcript_id, seqs in self.single_variant_protein_VCF_seq_extractor.extract_all():
            ref_seq = self.transcript_seq_extractor.get_protein_seq(transcript_id)
            metadata = self.metadatas.loc[transcript_id]
            for alt_seq in seqs:
                yield {
                    'input': {
                        'ref_seq': ref_seq,
                        'alt_seq': alt_seq,
                    },
                    # pandas.core.series.Series
                    'metadata': metadata
                }
                
    def __next__(self):
        for unit in self.sequences:
            # empty generator
            if type(unit) == None:
                break
            return unit
        raise StopIteration
