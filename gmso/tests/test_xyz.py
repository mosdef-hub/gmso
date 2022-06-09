import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn


class TestXYZ(BaseTest):
    def test_read_xyz(self):
        top = Topology.load(get_fn("ethane.xyz"))
        assert top.n_sites == 8
        assert top.n_connections == 0
        assert set([type(site.position) for site in top.sites]) == {
            u.unyt_array
        }
        assert set([site.position.units for site in top.sites]) == {u.nm}

        top = Topology.load(get_fn("cu_block.xyz"))
        assert top.n_sites == 108

        assert top.n_connections == 0
        assert set([type(site.position) for site in top.sites]) == {
            u.unyt_array
        }
        assert set([site.position.units for site in top.sites]) == {u.nm}

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            Topology.load(get_fn("too_few_atoms.xyz"))
        with pytest.raises(ValueError):
            Topology.load(get_fn("too_many_atoms.xyz"))

    def test_write_xyz(self):
        top = Topology.load(get_fn("ethane.xyz"))
        top.save("tmp.xyz")

    def test_full_io(self):
        original_top = Topology.load(get_fn("ethane.xyz"))

        original_top.save("full_conversion.xyz")
        new_top = Topology.load("full_conversion.xyz")

        assert original_top.n_sites == new_top.n_sites
        assert original_top.n_connections == new_top.n_connections
        assert_allclose_units(
            original_top.positions, new_top.positions, rtol=1e-5, atol=1e-8
        )
