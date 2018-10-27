"""Test MMSpliceDl
"""
from kipoiseq.dataloaders.splicing import MMSpliceDl
from kipoiseq.transforms.transforms import SplitSplicingSeq
import pytest


@pytest.fixture
def fasta_file():
    return "tests/data/sample.5kb.fa"


@pytest.fixture
def gtf_file():
    return "tests/data/sample.gtf"


def test_splicedataset(gtf_file,
                       fasta_file,
                       transform=SplitSplicingSeq()):
    dl = MMSpliceDl(gtf_file,
                    fasta_file,
                    intron5prime_len=100,
                    intron3prime_len=100,
                    transform=transform)
    dt = next(dl)
    assert dt['inputs']['seq']['intron5prime'] == 'GTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTA'
    assert dt['inputs']['seq']['donor'] == 'TACGTNNNNNNNNNNNNN'
    assert dt['inputs']['seq']['intron3prime'] == 'NNNNNNNN'
    assert dt['inputs']['seq']['acceptor'] == 'GTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTT'
