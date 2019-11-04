import pytest
import numpy as np
from kipoiseq.utils import parse_alphabet, parse_dtype


def test_parse_alphabet():
    assert parse_alphabet(['A', 'C']) == ['A', 'C']
    assert parse_alphabet('AC') == ['A', 'C']


def test_parse_type():
    with pytest.raises(Exception):
        parse_dtype('string')
    with pytest.raises(Exception):
        parse_dtype('np.float322')

    assert parse_dtype('float') == float
    assert parse_dtype(float) == float
    assert parse_dtype("np.float32") == np.float32
