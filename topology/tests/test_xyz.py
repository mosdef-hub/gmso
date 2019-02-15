import pytest

from topology.formats.xyz import read_xyz
from topology.utils.io import get_fn


def test_read_xyz():
    top = read_xyz(get_fn('ethane.xyz'))
    assert top.n_sites == 8

    top = read_xyz(get_fn('cu_block.xyz'))
    assert top.n_sites == 108

def test_wrong_n_atoms():
    with pytest.raises(ValueError):
        read_xyz(get_fn('too_few_atoms.xyz'))
    with pytest.raises(ValueError):
        read_xyz(get_fn('too_many_atoms.xyz'))
