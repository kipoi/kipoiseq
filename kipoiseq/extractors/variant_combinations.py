from typing import Iterable
from itertools import product
from kipoiseq import Interval, Variant
from kipoiseq.utils import alphabets
from kipoiseq.extractors import FastaStringExtractor


class VariantCombinator:

    def __init__(self, fasta_file: str, bed_file: str = None,
                 variant_type='snv', alphabet='DNA'):
        if variant_type not in {'all', 'snv', 'in', 'del'}:
            raise ValueError("variant_type should be one of "
                             "{'all', 'snv', 'in', 'del'}")

        self.bed_file = bed_file
        self.fasta = fasta_file
        self.fasta = FastaStringExtractor(fasta_file)
        self.variant_type = variant_type
        self.alphabet = alphabets[alphabet]

    def combination_variants_snv(self, interval: Interval) -> Iterable[Variant]:
        """Returns all the possible variants in the regions.

          interval: interval of variants
        """
        seq = self.fasta.extract(interval)
        for pos, ref in zip(range(interval.start, interval.end), seq):
            pos = pos + 1  # 0 to 1 base
            for alt in self.alphabet:
                if ref != alt:
                    yield Variant(interval.chrom, pos, ref, alt)

    def combination_variants_insertion(self, interval, length=2) -> Iterable[Variant]:
        """Returns all the possible variants in the regions.

          interval: interval of variants
          length: insertions up to length
        """
        if length < 2:
            raise ValueError('length argument should be larger than 1')

        seq = self.fasta.extract(interval)
        for pos, ref in zip(range(interval.start, interval.end), seq):
            pos = pos + 1  # 0 to 1 base
            for l in range(2, length + 1):
                for alt in product(self.alphabet, repeat=l):
                    yield Variant(interval.chrom, pos, ref, ''.join(alt))

    def combination_variants_deletion(self, interval, length=1) -> Iterable[Variant]:
        """Returns all the possible variants in the regions.

          interval: interval of variants
          length: deletions up to length
        """
        if length < 1 and length <= interval.width:
            raise ValueError('length argument should be larger than 0'
                             ' and smaller than interval witdh')

        seq = self.fasta.extract(interval)
        for i, pos in enumerate(range(interval.start, interval.end)):
            pos = pos + 1  # 0 to 1 base
            for j in range(1, length + 1):
                if i + j <= len(seq):
                    yield Variant(interval.chrom, pos, seq[i:i + j], '')

    def combination_variants(self, interval, variant_type='snv',
                             in_length=2, del_length=2) -> Iterable[Variant]:
        if variant_type in {'snv', 'all'}:
            yield from self.combination_variants_snv(interval)
        if variant_type in {'indel', 'in', 'all'}:
            yield from self.combination_variants_insertion(
                interval, length=in_length)
        if variant_type in {'indel', 'del', 'all'}:
            yield from self.combination_variants_deletion(
                interval, length=del_length)

    def __iter__(self) -> Iterable[Variant]:
        for line in open(self.bed_file):
            line = line.split('\t')
            interval = Interval(line[0], int(line[1]), int(line[2]))
            yield from self.combination_variants(interval, self.variant_type)
