"""Test BedDataset
"""
from kipoiseq.datasets.sequence import BedDataset
import numpy as np
import pytest
import pybedtools
from pybedtools import Interval


def write_tmp(string, tmpdir):
    p = tmpdir.mkdir("bed-files").join("file2.txt")
    p.write(string)
    return str(p)


def test_bed3(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\nchr1\t1\t3', tmpdir)
    bt = BedDataset(bed_file)
    assert bt.n_tasks == 0
    assert len(bt) == 2
    assert np.all(bt.df[0] == 'chr1')
    assert bt[0] == (Interval("chr1", 1, 2), {})
    assert bt[1] == (Interval("chr1", 1, 3), {})


def test_bed3_labels(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr1\t1\t3\t0\t1', tmpdir)
    bt = BedDataset(bed_file)
    assert np.all(bt.get_targets() == np.array([[1, 0],
                                                [0, 1]]))
    assert len(bt) == 2
    assert bt.n_tasks == 2
    assert np.all(bt.df[0] == 'chr1')
    assert bt[0][0] == Interval("chr1", 1, 2)
    assert np.all(bt[0][1] == np.array([1, 0]))

    assert bt[1][0] == Interval("chr1", 1, 3)
    assert np.all(bt[1][1] == np.array([0, 1]))
    assert len(bt) == 2


def test_incl_excl_chromosomes(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr2\t1\t3\t0\t1\nchr3\t1\t3\t0\t1', tmpdir)
    bt = BedDataset(bed_file)
    assert len(bt) == 3

    bt = BedDataset(bed_file, incl_chromosomes=['chr1'])
    assert len(bt) == 1
    assert bt[0][0] == Interval("chr1", 1, 2)

    bt = BedDataset(bed_file, excl_chromosomes=['chr1'])
    assert len(bt) == 2
    assert bt[0][0] == Interval("chr2", 1, 3)


def test_ambiguous_mask(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr2\t1\t3\t0\t1\nchr3\t1\t3\t0\t-1', tmpdir)
    bt = BedDataset(bed_file)
    assert len(bt) == 3
    assert np.all(bt[2][1] == np.array([0, -1]))

    # same as before
    bt = BedDataset(bed_file, ambiguous_mask=-1)
    assert len(bt) == 3
    assert np.all(bt[2][1] == np.array([0, -1]))
    assert np.all(bt.get_targets().max(axis=1) >= 0)


def test_ambiguous_mask2(tmpdir):
    # only ambigous regions are present
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr2\t1\t3\t0\t1\nchr3\t1\t3\t-1\t-1', tmpdir)
    bt = BedDataset(bed_file, ambiguous_mask=-1)
    assert len(bt) == 2
    assert np.all(bt.get_targets().max(axis=1) >= 0)


def test_num_chr(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr2\t1\t3\t0\t1\nchr3\t1\t3\t0\t-1', tmpdir)
    bt = BedDataset(bed_file, num_chr=True)
    assert len(bt) == 3
    assert bt[0][0].chrom == '1'


def test_label_dtype(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\t1\t0\nchr2\t1\t3\t0\t1', tmpdir)
    bt = BedDataset(bed_file, label_dtype=bool)
    assert len(bt) == 2
    assert bt[0][1].dtype == bool
    assert bt.get_targets().dtype == bool


def test_more_columns(tmpdir):
    bed_file = write_tmp('chr1\t1\t2\tinterval1\t1\t0\nchr2\t1\t3\tinterval2\t0\t1', tmpdir)
    with pytest.raises(Exception):
        bt = BedDataset(bed_file, label_dtype=bool)
    bt = BedDataset(bed_file, bed_columns=4, label_dtype=bool)
    assert bt[0][0].name == 'interval1'
    assert bt[1][0].name == 'interval2'

    with pytest.raises(Exception):
        bt = BedDataset(bed_file)
