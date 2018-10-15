import os
import numpy as np
import pytest
from pybedtools import Interval
from kipoiseq.transforms.functional import one_hot_dna
from kipoiseq.datasets.sequence import SeqStringDataset, SeqDataset, parse_dtype, BedDataset


@pytest.fixture
def fasta_file():
    return "tests/data/sample.fasta"


@pytest.fixture
def intervals_file():
    return "tests/data/sample_intervals.bed"


# fasta_file = fasta_file()
# intervals_file = intervals_file()


def test_min_props():
    # minimal set of properties that need to be specified on the object
    min_set_props = ["output_schema", "type", "defined_as", "info", "args", "dependencies", "postprocessing",
                     "source", "source_dir"]

    for Dl in [SeqStringDataset, SeqDataset]:
        props = dir(Dl)
        assert all([el in props for el in min_set_props])


def test_parse_dtype():
    dtypes = {'int': int, 'string': str, 'float': float, 'bool': bool}
    assert all([parse_dtype(dt) == dtypes[dt] for dt in dtypes.keys()])
    assert all([parse_dtype(dt) == dt for dt in dtypes.values()])
    with pytest.raises(Exception):
        parse_dtype("int8")
    assert parse_dtype(None) is None


def test_fasta_based_dataset(intervals_file, fasta_file):
    # just test the functionality
    dl = SeqStringDataset(intervals_file, fasta_file)
    ret_val = dl[0]
    assert isinstance(ret_val["inputs"], np.ndarray)
    assert ret_val["inputs"].shape == ()
    # # test with set wrong seqlen:
    # dl = SeqStringDataset(intervals_file, fasta_file, required_seq_len=3)
    # with pytest.raises(Exception):
    #     dl[0]

    dl = SeqStringDataset(intervals_file, fasta_file, label_dtype="string")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.str_)
    dl = SeqStringDataset(intervals_file, fasta_file, label_dtype="int")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.int_)
    dl = SeqStringDataset(intervals_file, fasta_file, label_dtype="bool")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.bool_)
    vals = dl.load_all()
    assert vals['inputs'][0] == 'GT'


def test_seq_dataset(intervals_file, fasta_file):
    dl = SeqDataset(intervals_file, fasta_file)
    ret_val = dl[0]

    assert np.all(ret_val['inputs'] == one_hot_dna("GT"))
    assert isinstance(ret_val["inputs"], np.ndarray)
    assert ret_val["inputs"].shape == (2, 4)


# download example files
@pytest.mark.parametrize("cls", [SeqStringDataset, SeqDataset])
def test_examples_exist(cls):
    ex = cls.init_example()
    example_kwargs = cls.example_kwargs
    bed_fh = open(example_kwargs['intervals_file'], "r")
    bed_entries = len([el for el in bed_fh])
    bed_fh.close()

    dl_entries = 0
    for out in ex:
        dl_entries += 1
    assert dl_entries == len(ex)
    assert len(ex) == bed_entries
