import numpy as np
import unyt as u
import parmed as pmd
import pytest

from topology.formats.gro import read_gro, write_gro
from topology.external.convert_parmed import from_parmed
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

    def test_write_gro(self):
        top = from_parmed(pmd.load_file(get_fn('ethane.gro'), structure=True))

        write_gro(top, 'out.gro')

    def test_write_gro_non_orthogonal(self):
        top = from_parmed(pmd.load_file(get_fn('ethane.gro'), structure=True))
        top.box.angles = u.degree * [90, 90, 120]

        write_gro(top, 'out.gro')
