import pytest


@pytest.mark.parametrize("alphabet_axis", list(range(0, 4)))
@pytest.mark.parametrize("dummy_axis", [None] + list(range(0, 4)))
def test_seq_dataset_reshape(alphabet_axis, dummy_axis):
    from kipoiseq.datasets.sequence import SeqDataset
    seq_len, alphabet_len = 3, 4

    kwargs = SeqDataset.example_kwargs
    kwargs['auto_resize_len'] = seq_len
    kwargs['alphabet_axis'] = alphabet_axis
    kwargs['dummy_axis'] = dummy_axis

    dummy_axis_int = dummy_axis
    if dummy_axis is None:
        dummy_axis_int = -2

    if (alphabet_axis == dummy_axis_int) or (alphabet_axis == -1) or (dummy_axis_int == -1) or \
            (alphabet_axis >= 3) or (dummy_axis_int >= 3) or ((alphabet_axis >= 2) and (dummy_axis is None)):
        with pytest.raises(Exception):
            seq_dataset = SeqDataset(**kwargs)
        return None

    seq_dataset = SeqDataset(**kwargs)

    # test the single sample works
    reshaped = seq_dataset[0]['inputs']
    for i in range(len(reshaped.shape)):
        if i == dummy_axis:
            assert reshaped.shape[i] == 1
        elif i == alphabet_axis:
            assert reshaped.shape[i] == alphabet_len
        else:
            assert reshaped.shape[i] == seq_len

