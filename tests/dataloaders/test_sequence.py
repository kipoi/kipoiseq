import numpy as np
import pytest
from copy import deepcopy
from kipoi.utils import override_default_kwargs
from kipoiseq.transforms.functional import one_hot_dna
from kipoiseq.dataloaders.sequence import StringSeqIntervalDl, SeqIntervalDl


@pytest.fixture
def fasta_file():
    return "tests/data/sample.fasta"


@pytest.fixture
def intervals_file():
    return "tests/data/sample_intervals.bed"

@pytest.fixture
def intervals_file_strand():
    return "tests/data/sample_interval_strand.bed"

# fasta_file = fasta_file()
# intervals_file = intervals_file()


def test_min_props():
    # minimal set of properties that need to be specified on the object
    min_set_props = ["output_schema", "type", "defined_as", "info", "args", "dependencies", "postprocessing",
                     "source", "source_dir"]

    for Dl in [StringSeqIntervalDl, SeqIntervalDl]:
        props = dir(Dl)
        assert all([el in props for el in min_set_props])


def test_fasta_based_dataset(intervals_file, fasta_file):
    # just test the functionality
    dl = StringSeqIntervalDl(intervals_file, fasta_file)
    ret_val = dl[0]
    assert isinstance(ret_val["inputs"], np.ndarray)
    assert ret_val["inputs"].shape == ()
    # # test with set wrong seqlen:
    # dl = StringSeqIntervalDl(intervals_file, fasta_file, required_seq_len=3)
    # with pytest.raises(Exception):
    #     dl[0]

    dl = StringSeqIntervalDl(intervals_file, fasta_file, label_dtype="str")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.str_)
    dl = StringSeqIntervalDl(intervals_file, fasta_file, label_dtype="int")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.int_)
    dl = StringSeqIntervalDl(intervals_file, fasta_file, label_dtype="bool")
    ret_val = dl[0]
    assert isinstance(ret_val['targets'][0], np.bool_)
    vals = dl.load_all()
    assert vals['inputs'][0] == 'GT'

def test_use_strand(intervals_file_strand, fasta_file):
    dl = StringSeqIntervalDl(intervals_file_strand, fasta_file, use_strand=True)
    vals = dl.load_all()
    assert vals['inputs'][0] == 'AC'

def test_seq_dataset(intervals_file, fasta_file):
    dl = SeqIntervalDl(intervals_file, fasta_file)
    ret_val = dl[0]

    assert np.all(ret_val['inputs'] == one_hot_dna("GT"))
    assert isinstance(ret_val["inputs"], np.ndarray)
    assert ret_val["inputs"].shape == (2, 4)


@pytest.fixture
def example_kwargs():
    return SeqIntervalDl.example_kwargs


@pytest.mark.parametrize("alphabet_axis", list(range(0, 4)))
@pytest.mark.parametrize("dummy_axis", [None] + list(range(0, 4)))
def test_seq_dataset_reshape(alphabet_axis, dummy_axis, example_kwargs):
    seq_len, alphabet_len = 3, 4

    kwargs = example_kwargs
    kwargs['auto_resize_len'] = seq_len
    kwargs['alphabet_axis'] = alphabet_axis
    kwargs['dummy_axis'] = dummy_axis

    dummy_axis_int = dummy_axis
    if dummy_axis is None:
        dummy_axis_int = -2

    if (alphabet_axis == dummy_axis_int) or (alphabet_axis == -1) or (dummy_axis_int == -1) or \
            (alphabet_axis >= 3) or (dummy_axis_int >= 3) or ((alphabet_axis >= 2) and (dummy_axis is None)):
        with pytest.raises(Exception):
            seq_dataset = SeqIntervalDl(**kwargs)
        return None

    seq_dataset = SeqIntervalDl(**kwargs)

    # test the single sample works
    reshaped = seq_dataset[0]['inputs']
    for i in range(len(reshaped.shape)):
        if i == dummy_axis:
            assert reshaped.shape[i] == 1
        elif i == alphabet_axis:
            assert reshaped.shape[i] == alphabet_len
        else:
            assert reshaped.shape[i] == seq_len


# download example files
@pytest.mark.parametrize("cls", [StringSeqIntervalDl, SeqIntervalDl])
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


def test_output_schape():
    Dl = deepcopy(SeqIntervalDl)
    assert Dl.get_output_schema().inputs.shape == (None, 4)
    Dlc = override_default_kwargs(Dl, {"auto_resize_len": 100})
    assert Dlc.get_output_schema().inputs.shape == (100, 4)

    # original left intact
    assert Dl.get_output_schema().inputs.shape == (None, 4)

    Dlc = override_default_kwargs(
        Dl, {"auto_resize_len": 100, "dummy_axis": 1, "alphabet_axis": 2})
    assert Dlc.get_output_schema().inputs.shape == (100, 1, 4)
    Dlc = override_default_kwargs(
        Dl, {"auto_resize_len": 100, "dummy_axis": 2})
    assert Dlc.get_output_schema().inputs.shape == (100, 4, 1)
    # original left intact
    assert Dl.get_output_schema().inputs.shape == (None, 4)

    Dlc = override_default_kwargs(
        Dl, {"auto_resize_len": 100, "alphabet": "ACGTD"})
    assert Dlc.get_output_schema().inputs.shape == (100, 5)

    Dlc = override_default_kwargs(
        Dl, {"auto_resize_len": 160, "dummy_axis": 2, "alphabet_axis": 0})
    assert Dlc.get_output_schema().inputs.shape == (4, 160, 1)

    Dlc = override_default_kwargs(
        Dl, {"auto_resize_len": 160, "dummy_axis": 2, "alphabet_axis": 1})
    assert Dlc.get_output_schema().inputs.shape == (160, 4, 1)
    targets = Dlc.get_output_schema().targets
    assert targets.shape == (None,)

    Dlc = override_default_kwargs(Dl, {"ignore_targets": True})
    assert Dlc.get_output_schema().targets is None
    # reset back

    # original left intact
    assert Dl.get_output_schema().inputs.shape == (None, 4)
    assert Dl.get_output_schema().targets.shape == (None, )
