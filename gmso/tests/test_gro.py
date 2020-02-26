import numpy as np
import unyt as u
import pytest

from gmso.formats.gro import read_gro, write_gro
from gmso.external.convert_parmed import from_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, import_, has_parmed
from gmso.utils.testing import allclose


if has_parmed:
    pmd = import_('parmed')

@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
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
