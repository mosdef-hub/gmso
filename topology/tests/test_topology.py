import numpy as np
import pytest
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.atom_type import AtomType
from topology.core.connection import Connection
from topology.core.connection_type import ConnectionType
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

    def test_top_update(self):
        top = Topology()
        top.update_top()
        assert top.n_sites == 0
        assert len(top.atom_types) == 0
        assert len(top.atom_type_functionals) == 0
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_functionals) == 0

        atomtype = AtomType()
        site1 = Site(name='site1', atom_type=atomtype)
        top.add_site(site1)
        site2 = Site(name='site2', atom_type=atomtype)
        top.add_site(site2)
        assert top.n_sites == 2
        assert len(top.atom_types) == 0
        assert len(top.atom_type_functionals) == 0
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_functionals) == 0
        top.update_atom_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_functionals) == 1
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_functionals) == 0


        ctype = ConnectionType()
        connection_12 = Connection(site1=site1, site2=site2, 
                connection_type=ctype)
        top.add_connection(connection_12)
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_functionals) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 0
        assert len(top.connection_type_functionals) == 0
        top.update_connection_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_functionals) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_functionals) == 1

        site1.atom_type = AtomType(expression='sigma*epsilon')
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_functionals) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_functionals) == 1
        top.update_atom_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 2
        assert len(top.atom_type_functionals) == 2
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_functionals) == 1

