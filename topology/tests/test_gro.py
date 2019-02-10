import numpy as np

from topology.formats.gro import read_gro
from topology.utils.io import get_fn


def test_read_gro():
    top = read_gro(get_fn('acn.gro'))

    assert top.name == 'ACN'
    assert top.n_sites == 6
    assert np.allclose(top.box.lengths, 4*np.ones(3))

    top = read_gro(get_fn('350-waters.gro'))

    assert top.name == 'Generic title'
    assert top.n_sites == 1050
    assert np.allclose(top.box.lengths, 2.20866*np.ones(3))
