import numpy as np
import pytest
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.box import Box
from topology.tests.base_test import BaseTest
from topology.testing.utils import allclose


class TestTopology(BaseTest):
    def test_new_topology(self):
        top = Topology(name='mytop')
        assert top.name == 'mytop'

    def test_add_site(self):
        top = Topology()
        site = Site(name='site')

        assert top.n_sites == 0
        top.add_site(site)
        assert top.n_sites == 1

    def test_add_connection(self):
        top = Topology()
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        connect = Connection(site1=site1, site2=site2)

        top.add_site(site1)
        top.add_site(site2)

        top.update_connection_list()

        assert len(top.connection_list) == 1

    def test_add_box(self):
        top = Topology()
        box = Box(2*u.nm*np.ones(3))

        assert top.box is None
        top.box = box
        assert top.box is not None
        assert allclose(top.box.lengths, u.nm*2*np.ones(3))

    def test_positions_dtype(self):
        top = Topology()
        site1 = Site(name='site1')
        top.add_site(site1)

        assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
        assert set([site.position.units for site in top.site_list]) == {u.nm}

        assert top.positions().dtype == float
        assert top.positions().units == u.nm
        assert isinstance(top.positions(), u.unyt_array)
