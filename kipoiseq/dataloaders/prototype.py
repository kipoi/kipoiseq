import abc
from typing import List, Union, Iterable, Tuple, Dict, Optional, Type, Iterator

from kipoi.data import SampleIterator
from kipoiseq.transforms import OneHot
from kipoiseq.extractors import GenericMultiIntervalSeqExtractor, BaseMultiIntervalFetcher, \
    GTFMultiIntervalFetcher, BaseExtractor, FastaStringExtractor, VariantSeqExtractor, \
    MultiSampleVCF
from kipoiseq.dataclasses import Interval, Variant
from kipoi.metadata import GenomicRanges

import pandas as pd
import pyranges

from kipoi_conda.dependencies import Dependencies
from kipoi.specs import Author

__all__ = [
    'GenericSingleSeqDataloader'
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
    Author(name='Alexander Karollus', github='Karollus'),
    Author(name='Florian R. Hölzlwimmer', github='hoeze')
]

##### Utility Functions #####

def _get_interval_metadata(
    key, 
    intervals : List[Interval],
    interval_attrs
):
    start = intervals[0].start # start of the encompassing genomic range
    end = intervals[-1].end # end of the encompassing genomic range
    # assemble the metadata
    # as we get a list of intervals, we take metadata from first (maybe questionable)
    # one could define a metadata aggregation function
    metadata_dict = {k:intervals[0].attrs.get(k, '') for k in interval_attrs}
    metadata_dict["ranges"] = GenomicRanges(intervals[0].chrom,
                                            start,
                                            end,
                                            key,
                                            intervals[0].strand)
    return metadata_dict

class IdentityTransform:
    
    def __call__(self, x):
        return x

##### Generic Sequence Dataloader #####

class GenericSingleSeqDataloader(SampleIterator):
    
    def __init__(
        self,
        interval_fetcher : BaseMultiIntervalFetcher,
        reference_sequence_extractor : BaseExtractor,
        sequence_transformer, # need a base class for this
        interval_attrs = []
    ):
        # Generic dataloader ingredients:
        # A source of (lists of) intervals (called as iterator)
        self.interval_fetcher = interval_fetcher.items()
        # A sequence source (usually a fasta extractor)
        self.reference_sequence_extractor = reference_sequence_extractor
        # A sequence transformer object (could be identity transform)
        self.sequence_transformer = sequence_transformer
        # Metadata attributes
        self.interval_attrs = interval_attrs
    
    def __next__(self):
        # Generic single seq dataloader step:
        # Fetch the next list of intervals from the interval source
        key, intervals = next(self.interval_fetcher)
        # extract and assemble the corresponding sequence
        intervals.sort(key=lambda x: x.start)
        seq_parts = [self.reference_sequence_extractor.extract(interval) for interval in intervals]
        if self.reference_sequence_extractor.use_strand and intervals[0].strand == "-":
            seq_parts.reverse()
        seq = "".join(seq_parts)
        # transform the sequence
        transformed_seq = self.sequence_transformer(seq)
        # get metadata
        metadata_dict = _get_interval_metadata(key, intervals, self.interval_attrs)
        # return all as dictionary
        return {
            "inputs": {
                "seq": transformed_seq,
            },
            "metadata": metadata_dict
        }
    
    def __iter__(self):
        return self        

##### Variant Dataloader #####

### Variant Fetcher ###

class BaseVariantFetcher:
        
    def fetch_variants(
        self,
        interval : Interval
    ) -> Iterator[Variant]:
        raise NotImplementedError

# Like many existing variant loading classes, this is not super efficient
# One could probably make it more efficient using a smart prefetching/batching strategy
class VCFVariantFetcher(BaseVariantFetcher):
    
    def __init__(
        self,
        vcf_file : str
    ):
        self.vcf = MultiSampleVCF(vcf_file)
    
    def fetch_variants(
        self,
        interval : Interval
    ) -> Iterator[Variant]:
        return self.vcf.fetch_variants(interval)        

### Insertion Strategy ###    
    
class VariantStrategy:
    
    def generate_alt_seq_templates(
        self,
        intervals_with_variants : List[Tuple[Interval, Iterator[Variant]]]
    ) -> Iterator[List[Tuple[Interval, List[Variant]]]]:
        raise NotImplementedError()

class MultiVariantStrategy(VariantStrategy):
    """
    This class will yield a template to build one alt-seq
    which has all matching variants inserted
    If no variants intersect the intervals, nothing is returned.
    The class assumes that there are no "mutually exclusive" variants
    """
    
    def generate_alt_seq_templates(
        self,
        intervals_with_variants : List[Tuple[Interval, Iterator[Variant]]]
    ) -> Iterator[List[Tuple[Interval, List[Variant]]]]:
        alt_seq_intervals_with_variants = []
        has_variant = False
        for interval, var_iter in intervals_with_variants:
            variants = [var for var in var_iter]
            if len(variants) > 0:
                has_variant = True
            alt_seq_intervals_with_variants.append((interval, variants))
        if has_variant:
            yield alt_seq_intervals_with_variants

class SingleVariantStrategy(VariantStrategy):
    """
    This class will yield a template to build one alt-seq 
    for every variant that intersects the interval
    """
    
    def generate_alt_seq_templates(
        self,
        intervals_with_variants : List[Tuple[Interval, Iterator[Variant]]]
    ) -> Iterator[List[Tuple[Interval, List[Variant]]]]:
        for idx, (interval, var_iter) in enumerate(intervals_with_variants):
            for variant in var_iter:
                alt_seq_intervals_with_variant = [(interval, []) for 
                                                   interval, _ in intervals_with_variants]
                alt_seq_intervals_with_variant[idx] = (interval, [variant])
                yield alt_seq_intervals_with_variant
                
### Variant Sequence Extractor ###    

class VariantSequenceExtractor:
    
    def __init__(
        self,
        reference_sequence_extractor : BaseExtractor,
        anchor = 0,
        fixed_len = True
    ):
        self.anchor = anchor
        self.fixed_len = fixed_len
        self.extractor = VariantSeqExtractor(
            reference_sequence = reference_sequence_extractor
        )
        
    def extract(
        self,
        interval, 
        variants
    ):
        return self.extractor.extract(interval, variants, 
                anchor=self.anchor, fixed_len=self.fixed_len)

### The Generic Dataloader ###
        
class GenericVariantDataloader(SampleIterator):
    
    def __init__(
        self,
        interval_fetcher,
        variant_source,
        variant_insertion_strategy,
        reference_sequence_extractor,
        variant_sequence_extractor,
        sequence_transformer,
        interval_attrs,
        anchor = 0,
        fixed_len = False
    ):
        # Generic variant dataloader ingredients:
        # A interval source (called as iterator)
        self.interval_fetcher = interval_fetcher.items()
        # A variant source
        self.variant_source = variant_source
        #A variant insertion strategy
        self.variant_insertion_strategy = variant_insertion_strategy
        # A sequence source
        self.reference_sequence_extractor = reference_sequence_extractor
        # A variant sequence extractor
        self.variant_sequence_extractor = variant_sequence_extractor
        # A sequence transformer object (could be identity transform)
        self.sequence_transformer = sequence_transformer
        # Build iterator object
        self.generator = self._generator()
        # interval attributes
        self.interval_attrs = interval_attrs
    
    def _build_sequences(self, alt_seq_template):
        ref_seq_parts = [self.reference_sequence_extractor.extract(interval) 
                        for interval, variants in alt_seq_template]
        alt_seq_parts = []
        for idx, (interval, variants) in enumerate(alt_seq_template):
            if len(variants) > 0:
                alt_seq_parts.append(self.variant_sequence_extractor.extract(interval, variants))
            else:
                alt_seq_parts.append(ref_seq_parts[idx])
        if self.reference_sequence_extractor.use_strand and alt_seq_template[0][0].strand == "-":
            ref_seq_parts.reverse()
            alt_seq_parts.reverse()
        ref_seq = "".join(ref_seq_parts)
        alt_seq = "".join(alt_seq_parts)
        return ref_seq, alt_seq
            
    def _generator(self):
        # Iterate through the interval source
        for key, intervals in self.interval_fetcher:
            intervals.sort(key=lambda x: x.start)
            # Fetch variants corresponding to the intervals
            intervals_with_variants = [(interval, self.variant_source.fetch_variants(interval))
                                      for interval in intervals]
            # Get an alternative sequence template from the variant insertion strategy
            alt_seq_template_iterator = (self.variant_insertion_strategy
                                         .generate_alt_seq_templates(intervals_with_variants))
            # Get templates and build ref_seq, alt_seq pairs to supply to model
            for alt_seq_template in alt_seq_template_iterator:
                # extract and assemble the corresponding ref and alt sequences
                ref_seq, alt_seq = self._build_sequences(alt_seq_template)
                # transform the sequences
                transformed_ref_seq = self.sequence_transformer(ref_seq)
                transformed_alt_seq = self.sequence_transformer(alt_seq)
                # get metadata
                metadata_dict = _get_interval_metadata(key, intervals, self.interval_attrs)
                # yield all as dictionary
                yield {
                    "inputs": {
                        "ref_seq": transformed_ref_seq,
                        "alt_seq": transformed_alt_seq,
                    },
                    "metadata": metadata_dict
                }    
                
    def __next__(self):
        return next(self.generator)
    
    def __iter__(self):
        return self        
