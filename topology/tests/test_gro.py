import numpy as np

from topology.formats.gro import read_gro
from topology.utils.io import get_fn


def test_read_gro():
    top = read_gro(get_fn('acn.gro'))
    import pdb; pdb.set_trace()
    assert top.name == 'ACN'
    assert top.n_sites == 6
    assert np.allclose(top.box.lengths, 4*np.ones(3))
