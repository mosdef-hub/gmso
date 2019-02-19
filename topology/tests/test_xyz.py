import unyt as u
import pytest

from topology.formats.xyz import read_xyz
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn


class TestXYZ(BaseTest):
    def test_read_xyz(self):
        top = read_xyz(get_fn('ethane.xyz'))
        assert top.n_sites == 8
        assert top.n_connections == 0
        assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
        assert set([site.position.units for site in top.site_list]) == {u.nm}

        top = read_xyz(get_fn('cu_block.xyz'))
        assert top.n_sites == 108

        assert top.n_connections == 0
        assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
        assert set([site.position.units for site in top.site_list]) == {u.nm}

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            read_xyz(get_fn('too_few_atoms.xyz'))
        with pytest.raises(ValueError):
            read_xyz(get_fn('too_many_atoms.xyz'))
