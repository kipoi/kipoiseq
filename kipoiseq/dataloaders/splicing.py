from kipoi.data import Dataset, kipoi_dataloader
from kipoi.metadata import GenomicRanges
from kipoi.specs import Author, Dependencies
from kipoi.data import SampleIterator

import gffutils
from pyfaidx import Fasta
import pickle

# general dependencies
# bioconda::genomelake', TODO - add genomelake again once it gets released with pyfaidx to bioconda
deps = Dependencies(conda=['bioconda::pyfaidx', 'numpy', 'pandas'],
                    pip=['kipoiseq', 'kipoi'])
package_authors = [Author(name='Jun Cheng', github='s6juncheng')]

__all__ = ['ExonInterval', 'generate_exons', 'MMSpliceDl']

# python 2.7 compatibility

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

try:
    ModuleNotFoundError
except NameError:
    ModuleNotFoundError = ImportError
# ------------


class ExonInterval(gffutils.Feature):

    def __init__(self, order=-1, **kwargs):
        super(ExonInterval, self).__init__(**kwargs)
        self.order = order
        self.name = self.attributes["exon_id"][0]
        self.transcript_id = self.attributes["transcript_id"][0]
        self.gene_id = self.attributes["gene_id"][0]
        self.Exon_Start = self.start  # with overhang, will be set later
        self.Exon_End = self.end  # with overhang
        self.overhang = (0, 0)
        self._isLast = False

    @property
    def isLast(self):
        return self._isLast

    @isLast.setter
    def isLast(self, value):
        self._isLast = value

    @property
    def isFirst(self):
        return self.order == 1

    @property
    def grange(self):
        return GenomicRanges(self.chrom,
                             self.start,
                             self.end,
                             self.transcript_id,
                             self.strand)

    def __str__(self):
        return '{0}_{1}_{2}:{3}'.format(
            self.chrom, self.Exon_Start, self.Exon_End, self.strand)

    @property
    def to_dict(self):
        return {'isLast': self.isLast,
                'isFirst': self.isFirst,
                'order': self.order,
                'name': self.name,
                'gene_id': self.gene_id,
                'Exon_Start': self.Exon_Start,
                'Exon_End': self.Exon_End,
                'intronl_len': self.overhang[0],
                'intronr_len': self.overhang[1],
                'seqid': self.seqid,
                'strand': self.strand,
                'start': self.start,
                'end': self.end}

    @classmethod
    def from_feature(cls,
                     feature,
                     overhang=(0, 0)):
        # convert gffutils.Feature to ExonInterval
        iv = cls(seqid=feature.chrom,
                 source=feature.source,
                 start=feature.start,  # exon start
                 end=feature.end,  # exon end
                 strand=feature.strand,
                 frame=feature.frame,
                 attributes=feature.attributes,
                 order=int(feature.attributes['exon_number'][0]))
        if iv.strand == "+":
            iv.start = iv.Exon_Start - overhang[0]
            iv.end = iv.Exon_End + overhang[1]
        else:
            iv.start = iv.Exon_Start - overhang[1]
            iv.end = iv.Exon_End + overhang[0]
        iv.overhang = overhang
        return iv

    @classmethod
    def from_exonfile(cls, exon, attributes, overhang=(0, 0)):
        iv = cls(seqid=exon.CHROM,
                 source='',
                 start=exon.Exon_Start,  # exon start
                 end=exon.Exon_End,  # exon end
                 strand=exon.strand,
                 frame='',
                 attributes=attributes,
                 order=attributes['order'])
        if iv.strand == "+":
            iv.start = iv.Exon_Start - overhang[0]
            iv.end = iv.Exon_End + overhang[1]
        else:
            iv.start = iv.Exon_Start - overhang[1]
            iv.end = iv.Exon_End + overhang[0]
        iv.overhang = overhang
        return iv

    def get_seq(self, fasta, use_strand=True):
        seq = self.sequence(fasta, use_strand=use_strand)
        seq = seq.upper()
        return seq


def generate_exons(gtf_file,
                   overhang=(100, 100),  # overhang from the exon
                   gtf_db_path=":memory:",
                   out_file=None,
                   disable_infer_transcripts=True,
                   disable_infer_genes=True,
                   firstLastNoExtend=True,
                   source_filter=None):
    """
    Build IntervalTree object from gtf file for one feature unit (e.g. gene, exon). If give out_file, pickle it.
    Args:
        gtf_file: gtf format file or pickled Intervaltree object.
        overhang: flanking intron length to take along with exon. Corresponding to left (acceptor side) and right (donor side)
        gtf_db_path: (optional) gtf database path. Database for one gtf file only need to be created once
        out_file: (optional) file path to store the pickled Intervaltree obejct. Next time run it can be given to `gtf_file`
        disable_infer_transcripts: option to disable infering transcripts. Can be True if the gtf file has transcripts annotated.
        disable_infer_genes: option to disable infering genes. Can be True if the gtf file has genes annotated.
        firstLastNoExtend: if True, overhang is not taken for 5' of the first exon, or 3' of the last exon of a gene.
        source_filter: gene source filters, such as "protein_coding" filter for protein coding genes
    """
    try:
        gtf_db = gffutils.interface.FeatureDB(gtf_db_path)
    except ValueError:
        gtf_db = gffutils.create_db(
            gtf_file,
            gtf_db_path,
            disable_infer_transcripts=disable_infer_transcripts,
            disable_infer_genes=disable_infer_genes)

    genes = gtf_db.features_of_type('gene')
    default_overhang = overhang
    for gene in genes:
        if source_filter is not None:
            if gene.source != source_filter:
                continue
        for exon in gtf_db.children(gene, featuretype='exon'):
            isLast = False  # track whether is last exon
            if firstLastNoExtend:
                if (gene.strand == "+" and exon.end == gene.end) or (gene.strand == "-" and exon.start == gene.start):
                    overhang = (overhang[0], 0)
                    isLast = True
                elif (gene.strand == "+" and exon.start == gene.start) or (gene.strand == "-" and exon.end == gene.end):
                    overhang = (0, overhang[1])
            exon = ExonInterval.from_feature(exon, overhang)
            exon.isLast = isLast
            overhang = default_overhang
            yield exon


@kipoi_dataloader(override={"dependencies": deps, 'info.authors': package_authors})
class MMSpliceDl(SampleIterator):
    """
    info:
        doc: >
            Dataloader for splicing models. With inputs as gtf annotation file and fasta file,
            each output is an exon sequence with flanking intronic seuqences. Intronic sequnce
            lengths specified by the users. Returned sequences are of the type np.array([str])
    args:
        gtf_file:
            doc: file path; Genome annotation GTF file
            example:
                url: https://raw.githubusercontent.com/kipoi/models/master/MMSplice/tests/data/test.gtf
                md5: b20607afe91ec20d6ee79ed95ab0e85b
        fasta_file:
            doc: Reference Genome sequence in fasta format
            example:
                url: https://raw.githubusercontent.com/kipoi/models/master/MMSplice/tests/data/hg19.nochr.chr17.fa
                md5: e3f6630a8323c4306469fdfe8d8b9448
        intron5prime_len:
            doc: 5' intronic sequence length to take.
        intron3prime_len:
            doc: 3' intronic sequence length to take.
        transform:
            doc: >
                transformation operation applied to the returned sequence. It needs to take seq,
                intron5prime_len and intron3prime_len as arguments.
    output_schema:
        inputs:
            name: seq
            shape: ()
            special_type: DNAStringSeq
            doc: exon sequence with flanking intronic sequence
            associated_metadata: ranges
        metadata:
            geneID:
                type: str
                doc: gene ID
            transcriptID:
                type: str
                doc: transcript id
            ranges:
                type: GenomicRanges
                doc: ranges that the sequences were extracted
    postprocessing:
        variant_effects:
          bed_input:
            - gtf_file
    """

    def __init__(self,
                 gtf_file,
                 fasta_file,
                 intron5prime_len=100,
                 intron3prime_len=100,
                 transform=None,
                 **kwargs):

        try:
            with open(gtf_file, 'rb') as f:
                self.exons = pickle.load(f)
        except (FileNotFoundError, pickle.UnpicklingError, ModuleNotFoundError):
            self.exons = generate_exons(gtf_file=gtf_file,
                                        overhang=(intron5prime_len, intron3prime_len),
                                        **kwargs)
        import six
        if isinstance(fasta_file, six.string_types):
            fasta = Fasta(fasta_file, as_raw=False)
        self.fasta = fasta
        self.transform = transform

    def __iter__(self):
        return self

    def __next__(self):
        exon = next(self.exons)
        seq = exon.get_seq(self.fasta).upper()
        if self.transform:
            seq = self.transform(seq,
                                 exon.overhang[0],  # intron5prime_len
                                 exon.overhang[1]  # intron3prime_len
                                 )
        return{
            'inputs': {
                'seq': seq
            },
            'metadata': {
                'ranges': exon.grange,
                'geneID': exon.gene_id,
                'transcriptID': exon.transcript_id
            }
        }

    # python 2.7 compatibility
    next = __next__


# TODO - implement
# class SpliceDonorSeqDl(SampleIterator):
#     def __init__(self,
#                  gtf_file,
#                  fasta_file,
#                  n_upstream=10,
#                  n_downstream=10,
#                  seq_transform=None):
#         raise NotImplementedError


# class SpliceAcceptorSeqDl(SampleIterator):
#
#     def __init__(self,
#                  gtf_file,
#                  fasta_file,
#                  n_upstream=40,
#                  n_downstream=10,
#                  seq_transform=None):
#         raise NotImplementedError
#     # seq_transform: ?
#     #   name: sdasd
#     #   args:
#     #     a: 1
#     #     b: 3
#     #
