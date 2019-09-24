"""Test kipoiseq.dataclasses

Tests to perform:

# Both
- make sure the immutable objects are really immutable
- str() and from_str
- =, hashable
-
# VCF
- from cyvcf2 

# Interval
- validation of the interval
- access to all the attributes
- from_pybedtools and to_pybedtools
- shift, swapt_strand, trim, etc
"""


from kipoiseq.dataclasses import Variant, Interval

def test_variant():
    v = Variant("chr1", 10, 'C', 'T')

    #
    pass
