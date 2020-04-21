from kipoiseq.extractors import SingleVariantProteinVCFSeqExtractor, TranscriptSeqExtractor
from kipoi.data import SampleIterator, kipoi_dataloader
from kipoi.specs import Author, Dependencies
import pandas
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
        self.metadatas = ((cds.loc[~cds.index.duplicated(keep='first')]).drop(columns=['Start', 'End'])).reset_index()
        # generator for all sequences with variants
        self.sequences = self._extractor()

    def __iter__(self):
        return self


    def __next__(self):
        for unit in self.sequences:
            # empty generator
            if type(unit) == None:
                break
            return unit
        raise StopIteration


    def _extractor(self):
        for transcript_id, seqs in self.single_variant_protein_VCF_seq_extractor.extract_all():
            # reference sequence
            ref_seq = self.transcript_seq_extractor.get_protein_seq(transcript_id)
            # get information for this transcript_id
            metadata = self.get_metadata(transcript_id)
            for alt_seq in seqs:
                yield self.pack_information(ref_seq, alt_seq, metadata)
    
    def get_metadata(self, transcript_id: str):
        return self.metadatas.loc[self.metadatas.transcript_id == transcript_id]
    
    def pack_information(self, ref_seq: str, alt_seq: str, metadata: pandas.core.frame.DataFrame):
        return {
                'input': {
                    'ref_seq': ref_seq,
                    'alt_seq': alt_seq,
                    },
                # pandas.core.frame.DataFrame
                'metadata': metadata
                }
