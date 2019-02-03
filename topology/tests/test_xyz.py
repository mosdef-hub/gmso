from topology.formats.xyz import read_xyz
from topology.utils.io import get_fn


def test_read_xyz():
    top = read_xyz(get_fn('ethane.xyz'))

    assert top.n_sites == 8
