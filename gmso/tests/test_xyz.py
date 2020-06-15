import unyt as u
import pytest

from gmso.formats.xyz import read_xyz, write_xyz
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn
from unyt.testing import assert_allclose_units

class TestXYZ(BaseTest):
    def test_read_xyz(self):
        top = read_xyz(get_fn('ethane.xyz'))
        assert top.n_sites == 8
        assert top.n_connections == 0
        assert set([type(site.position) for site in top.sites]) == {u.unyt_array}
        assert set([site.position.units for site in top.sites]) == {u.nm}

        top = read_xyz(get_fn('cu_block.xyz'))
        assert top.n_sites == 108

        assert top.n_connections == 0
        assert set([type(site.position) for site in top.sites]) == {u.unyt_array}
        assert set([site.position.units for site in top.sites]) == {u.nm}

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            read_xyz(get_fn('too_few_atoms.xyz'))
        with pytest.raises(ValueError):
            read_xyz(get_fn('too_many_atoms.xyz'))

    def test_write_xyz(self):
        top = read_xyz(get_fn('ethane.xyz'))
        write_xyz(top, 'tmp.xyz')

    def test_full_io(self):
        original_top = read_xyz(get_fn('ethane.xyz'))

        write_xyz(original_top, 'full_conversion.xyz')
        new_top = read_xyz('full_conversion.xyz')

        assert original_top.n_sites == new_top.n_sites
        assert original_top.n_connections == new_top.n_connections
        assert_allclose_units(original_top.positions, new_top.positions, rtol=1e-5, atol=1e-8)
