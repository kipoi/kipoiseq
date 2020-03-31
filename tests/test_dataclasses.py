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
import pytest
import cyvcf2
import pybedtools


def test_variant():
    v = Variant("chr1", 10, 'C', 'T')

    assert v.start == 9
    assert v.chrom == 'chr1'
    assert v.pos == 10
    assert v.ref == 'C'
    assert v.alt == 'T'
    assert isinstance(v.info, dict)
    assert len(v.info) == 0
    assert v.qual == 0
    assert v.filter == 'PASS'
    v.info['test'] = 10
    assert v.info['test'] == 10
    assert isinstance(str(v), str)

    # make sure the original got unchangd
    v2 = v.copy()
    v.info['test'] = 20
    assert v2.info['test'] == 10
    v.__repr__()

    # __str__, from_str
    assert v == Variant.from_str(str(v))

    # hash test
    assert isinstance(hash(v), int)
    assert hash(v) == hash(Variant.from_str(str(v)))

    # fixed arguments
    with pytest.raises(AttributeError):
        v.chrom = 'asd'
    with pytest.raises(AttributeError):
        v.pos = 10
    with pytest.raises(AttributeError):
        v.ref = 'asd'
    with pytest.raises(AttributeError):
        v.alt = 'asd'

    # non-fixed arguments
    v.id = 'asd'
    v.qual = 10
    v.filter = 'asd'
    v.source = 2

    assert isinstance(Variant("chr1", '10', 'C', 'T').pos, int)

    # from cyvcf2
    vcf = cyvcf2.VCF('tests/data/test.vcf.gz')
    cv = list(vcf)[0]

    v2 = Variant.from_cyvcf(cv)
    assert isinstance(v2.source, cyvcf2.Variant)


def test_interval():
    interval = Interval("chr1", 10, 20, strand='-')
    interval.__repr__()

    assert interval.start == 10
    assert interval.end == 20
    assert interval.chrom == 'chr1'
    assert interval.name == ''
    assert isinstance(interval.attrs, dict)
    assert len(interval.attrs) == 0
    interval.attrs['test'] = 10
    assert interval.attrs['test'] == 10
    assert isinstance(str(interval), str)
    assert interval.neg_strand

    assert interval.width() == 10
    assert len(interval) == 10

    # __str__, from_str
    assert interval == Interval.from_str(str(interval))

    # make sure the original got unchangd
    i2 = interval.copy()
    interval.attrs['test'] = 20
    assert i2.attrs['test'] == 10

    # hash test
    assert isinstance(hash(interval), int)
    assert hash(interval) == hash(Interval.from_str(str(interval)))

    # fixed arguments
    with pytest.raises(AttributeError):
        interval.chrom = 'asd'
    with pytest.raises(AttributeError):
        interval.start = 10
    with pytest.raises(AttributeError):
        interval.end = 300
    with pytest.raises(AttributeError):
        interval.strand = '+'
    assert interval.strand == '-'

    # non-fixed arguments
    interval.name = 'asd'
    interval.score = 10

    assert interval.unstrand().strand == '.'

    assert interval == Interval.from_pybedtools(interval.to_pybedtools())
    assert isinstance(interval.to_pybedtools(), pybedtools.Interval)

    i2 = interval.shift(10, use_strand=False)

    # original unchanged
    assert interval.start == 10
    assert interval.end == 20

    assert i2.start == 20
    assert i2.end == 30

    i2 = interval.shift(10)  # use_strand = True by default
    assert i2.start == 0
    assert i2.end == 10

    assert not interval.shift(20, use_strand=True).is_valid()

    i2 = interval.shift(15, use_strand=True).truncate()
    assert i2.start == 0
    assert i2.end == 5

    assert interval.center() == 15

    # resize
    i2 = interval.resize(11)
    assert i2.start == 10 and i2.end == 21

    i2 = interval.resize(12)
    assert i2.start == 9 and i2.end == 21

    i2 = interval.resize(9)
    assert i2.start == 11 and i2.end == 20

    i2 = interval.swap_strand()
    assert interval.strand == "-"
    assert i2.strand == "+"
    assert i2.strand == '+'
    assert len(i2) == 10
    i2 = i2.resize(11)
    assert i2.start == 9 and i2.end == 20

    i2 = interval.copy()
    assert i2.center(use_strand=True) == 15
    assert i2.center(use_strand=False) == 15

    i2 = interval.swap_strand()
    assert i2.strand == "+"

    i3 = i2.resize(11).shift(1)
    assert i3.center(use_strand=True) == 16
    assert i3.center(use_strand=False) == 15

    i3 = i2.resize(11).shift(1)
    assert i3.center(use_strand=True) == 16
    assert i3.center(use_strand=False) == 15

    i2 = interval.trim(1, 10)
    assert i2.start == 10 and i2.end == 19

    i2 = interval.trim(1, 10, use_strand=False)
    assert i2.start == 11 and i2.end == 20
