import pytest

from gmso.core.topology import Topology
from gmso.core.bond import Bond
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest


class TestConnectivity(BaseTest):

    def test_methane_connectivity(self, methane):
        assert methane.n_bonds == 4
        assert methane.n_angles == 0
        assert methane.n_dihedrals == 0
        assert methane.n_impropers == 0

        methane.identify_connections()

        assert methane.n_bonds == 4
        assert methane.n_angles == 6
        assert methane.n_dihedrals == 0
        assert methane.n_impropers == 4

    def test_ethane_connectivity(self, ethane):
        assert ethane.n_bonds == 7
        assert ethane.n_angles == 0
        assert ethane.n_dihedrals == 0
        assert ethane.n_impropers == 0

        ethane.identify_connections()

        assert ethane.n_bonds == 7
        assert ethane.n_angles == 12
        assert ethane.n_dihedrals == 9
        assert ethane.n_impropers == 8

    def test_square(self):
        mytop = Topology()
        s1 = Atom(name='1')
        s2 = Atom(name='2')
        s3 = Atom(name='3')
        s4 = Atom(name='4')
        c12 = Bond(connection_members=[s1, s2])
        c23 = Bond(connection_members=[s2, s3])
        c34 = Bond(connection_members=[s3, s4])
        c41 = Bond(connection_members=[s4, s1])

        for site in [s1, s2, s3, s4]:
            mytop.add_site(site, update_types=False)

        for conn in [c12, c23, c34, c41]:
            mytop.add_connection(conn, update_types=False)

        assert mytop.n_bonds == 4
        assert mytop.n_angles == 0
        assert mytop.n_dihedrals == 0
        assert mytop.n_impropers == 0

        mytop.identify_connections()

        assert mytop.n_bonds == 4
        assert mytop.n_angles == 4
        assert mytop.n_dihedrals == 4
        assert mytop.n_impropers == 0

    def test_square_with_bridge(self):
        mytop = Topology()
        s1 = Atom(name='1')
        s2 = Atom(name='2')
        s3 = Atom(name='3')
        s4 = Atom(name='4')
        c12 = Bond(connection_members=[s1, s2])
        c23 = Bond(connection_members=[s2, s3])
        c34 = Bond(connection_members=[s3, s4])
        c41 = Bond(connection_members=[s4, s1])
        c24 = Bond(connection_members=[s2, s4])

        mytop.add_site(s1, update_types=False)
        mytop.add_site(s2, update_types=False)
        mytop.add_site(s3, update_types=False)
        mytop.add_site(s4, update_types=False)

        mytop.add_connection(c12, update_types=False)
        mytop.add_connection(c23, update_types=False)
        mytop.add_connection(c34, update_types=False)
        mytop.add_connection(c41, update_types=False)
        mytop.add_connection(c24, update_types=False)

        assert mytop.n_bonds == 5
        assert mytop.n_angles == 0
        assert mytop.n_dihedrals == 0
        assert mytop.n_impropers == 0

        mytop.identify_connections()

        assert mytop.n_bonds == 5
        assert mytop.n_angles == 8
        assert mytop.n_dihedrals == 6
        assert mytop.n_impropers == 2
