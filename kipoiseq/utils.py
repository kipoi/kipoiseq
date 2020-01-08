from six import string_types
import numpy as np
import pyranges
from kipoiseq.extractors import MultiSampleVCF, FastaStringExtractor


# alphabets:
DNA = ["A", "C", "G", "T"]
RNA = ["A", "C", "G", "U"]
AMINO_ACIDS = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H",
               "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

alphabets = {"DNA": DNA,
             "RNA": RNA,
             "AMINO_ACIDS": AMINO_ACIDS}


def to_scalar(obj):
    """Convert numpy scalar to native scalar
    """
    if isinstance(obj, np.generic):
        return np.asscalar(obj)
    else:
        return obj


def parse_alphabet(alphabet):
    if isinstance(alphabet, str):
        return list(alphabet)
    else:
        return alphabet


def parse_dtype(dtype):
    if isinstance(dtype, string_types):
        try:
            return eval(dtype)
        except Exception as e:
            raise ValueError(
                "Unable to parse dtype: {}. \nException: {}".format(dtype, e))
    else:
        return dtype


def _get_chrom_annotation(source):
    if type(source) == FastaStringExtractor:
        return set(source.fasta.keys())
    elif type(source) == MultiSampleVCF:
        return set(source.seqnames)
    elif type(source) == pyranges.PyRanges:
        return set(source.Chromosome)
    else:
        raise ValueError('source `%s` is not valid is not valid because '
                         ' source type `%s` is not supported.'
                         % (repr(source), type(source)))


def compare_chrom_annotation(sources, strategy='some', core_chroms=None):
    """Compares chromosome annotations from different sources.
        Throws exception iif annotations are not compatible.

    # Arguments:
        sources: list of different objects. vcf, fasta, pyranges are valid.
        strategy: comparison strategy. `some` means some intersection excepted
          or  `all` all chromosomes should be same.
        core_chroms: chromosomes must exist.

    # Returns:
        chroms common cross files.

    # Example:
        ```python
          >>> sources = [
                MultiSampleVCF(...),
                FastaStringExtractor(...),
                pyranges,
                pyranges,
                MultiSampleVCF(...)
              ]
          >>> compare_chrom_annotation(sources, strategy='all')
        ```
    """
    if not len(sources) > 1:
        raise ValueError(
            'At least two item should gived as sources to compare')

    chroms = list(map(_get_chrom_annotation, sources))

    if strategy == 'all':
        assert all(chroms[0] == i for i in chroms), \
            'chroms annotations are not all same.'
        return chroms[0]
    elif strategy == 'some':
        chrom_intersect = set.intersection(*chroms)
        assert len(chrom_intersect) > 0, \
            'there is not intersection between chromosomes.'
        return chrom_intersect
