import numpy as np
import pytest
import unyt as u

from topology import *

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
        connect = Bond(connection_members=[site1,site2])

        top.add_connection(connect)
        top.add_site(site1)
        top.add_site(site2)


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
        assert len(top.atom_type_expressions) == 0
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_expressions) == 0

        atomtype = AtomType()
        site1 = Site(name='site1', atom_type=atomtype)
        top.add_site(site1)
        site2 = Site(name='site2', atom_type=atomtype)
        top.add_site(site2)
        assert top.n_sites == 2
        #assert len(top.atom_types) == 0
        #assert len(top.atom_type_expressions) == 0
        #assert top.n_connections == 0
        #assert len(top.connection_types) == 0
        #assert len(top.connection_type_expressions) == 0
        #top.update_atom_types()
        #assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_expressions) == 0


        ctype = BondType()
        connection_12 = Bond(connection_members=[site1, site2], 
                connection_type=ctype)
        top.add_connection(connection_12)
        #assert top.n_sites == 2
        #assert len(top.atom_types) == 1
        #assert len(top.atom_type_expressions) == 1
        #assert top.n_connections == 1
        #assert len(top.connection_types) == 0
        #assert len(top.connection_type_expressions) == 0
        #top.update_connection_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1

        site1.atom_type = AtomType(expression='sigma*epsilon')
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1
        top.update_atom_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 2
        assert len(top.atom_type_expressions) == 2
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1

    def test_atomtype_update(self):
        top = Topology()

        assert top.n_sites == 0
        assert top.n_bonds == 0
        assert top.n_connections == 0

        atype1 = AtomType(expression='sigma + epsilon')
        atype2 = AtomType(expression='sigma * epsilon')
        site1 = Site('a', atom_type=atype1)
        site2 = Site('b', atom_type=atype2)
        top.add_site(site1)
        top.add_site(site2)
        #assert top.n_sites == 2
        #assert len(top.atom_types) == 0
        #assert len(top.atom_type_expressions) == 0

        #top.update_atom_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 2
        assert len(top.atom_type_expressions) == 2

    def test_bond_bondtype_update(self):
        top = Topology()

        atype1 = AtomType(expression='sigma + epsilon')
        atype2 = AtomType(expression='sigma * epsilon')
        site1 = Site('a', atom_type=atype1)
        site2 = Site('b', atom_type=atype2)
        btype = BondType()
        bond = Bond(connection_members=[site1, site2], connection_type=btype)
        top.add_site(site1)
        top.add_site(site2)
        top.add_connection(bond) 

        #assert top.n_connections == 1
        #assert top.n_bonds == 0
        #assert len(top.bond_types) == 0 
        #assert len(top.bond_type_expressions) == 0 

        #top.update_bond_list()
        #assert top.n_bonds == 1
        #assert len(top.bond_types) == 0 
        #assert len(top.bond_type_expressions) == 0 

        #top.update_bond_types()
        assert top.n_bonds == 1
        assert len(top.bond_types) == 1 
        assert len(top.bond_type_expressions) == 1

    def test_angle_angletype_update(self):
        top = Topology()

        atype1 = AtomType(expression='sigma + epsilon')
        atype2 = AtomType(expression='sigma * epsilon')
        site1 = Site('a', atom_type=atype1)
        site2 = Site('b', atom_type=atype2)
        site3 = Site('c', atom_type=atype2)
        atype = AngleType() 
        angle = Angle(connection_members=[site1, site2, site3], connection_type=atype)
        top.add_site(site1)
        top.add_site(site2)
        top.add_site(site3)
        top.add_connection(angle) 

        #assert top.n_connections == 1
        #assert top.n_angles == 0
        #assert len(top.angle_types) == 0 
        #assert len(top.angle_type_expressions) == 0 

        #top.update_angle_list()
        #assert top.n_angles == 1
        #assert len(top.angle_types) == 0 
        #assert len(top.angle_type_expressions) == 0 

        #top.update_angle_types()
        assert top.n_angles == 1
        assert len(top.angle_types) == 1 
        assert len(top.angle_type_expressions) == 1

