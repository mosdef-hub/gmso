import numpy as np
import unyt as u
import pytest

from topology.formats.gro import read_gro
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn
from topology.testing.utils import allclose


class TestGro(BaseTest):
    def test_read_gro(self):
        top = read_gro(get_fn('acn.gro'))

        assert top.name == 'ACN'
        assert top.n_sites == 6
        assert allclose(top.box.lengths, 4*np.ones(3)*u.nm)

        top = read_gro(get_fn('350-waters.gro'))

        assert top.name == 'Generic title'
        assert top.n_sites == 1050
        assert allclose(top.box.lengths, 2.20866*np.ones(3)*u.nm)

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            read_gro(get_fn('too_few_atoms.gro'))
        with pytest.raises(ValueError):
            read_gro(get_fn('too_many_atoms.gro'))
