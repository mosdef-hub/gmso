import mbuild as mb
import pytest

from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.topology import Topology
from gmso.external import from_mbuild
from gmso.tests.base_test import BaseTest
from gmso.utils.connectivity import generate_pairs_lists, identify_connections


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
        s1 = Atom(name="1")
        s2 = Atom(name="2")
        s3 = Atom(name="3")
        s4 = Atom(name="4")
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
        s1 = Atom(name="1")
        s2 = Atom(name="2")
        s3 = Atom(name="3")
        s4 = Atom(name="4")
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

    def test_index_only(self):
        atom1 = Atom(name="A")
        atom2 = Atom(name="B")
        atom3 = Atom(name="C")
        atom4 = Atom(name="D")
        atom5 = Atom(name="E")
        atom6 = Atom(name="F")

        bond1 = Bond(connection_members=[atom1, atom2])

        bond2 = Bond(connection_members=[atom2, atom3])

        bond3 = Bond(connection_members=[atom3, atom4])

        bond4 = Bond(connection_members=[atom2, atom5])

        bond5 = Bond(connection_members=[atom2, atom6])

        top = Topology()
        for site in [atom1, atom2, atom3, atom4, atom5, atom6]:
            top.add_site(site, update_types=False)

        for conn in [bond1, bond2, bond3, bond4, bond5]:
            top.add_connection(conn, update_types=False)

        top.update_topology()

        indices = identify_connections(top, index_only=True)
        assert len(indices["angles"]) == 7
        angle_indices = [
            (0, 1, 2),
            (1, 2, 3),
            (0, 1, 4),
            (5, 1, 2),
            (2, 1, 4),
            (4, 1, 5),
            (0, 1, 5),
        ]

        for idx_tuple in angle_indices:
            assert (
                idx_tuple in indices["angles"]
                or (idx_tuple[-1], idx_tuple[-2], idx_tuple[-3])
                in indices["angles"]
            )

        assert len(indices["dihedrals"]) == 3
        dihedral_indices = [(0, 1, 2, 3), (3, 2, 1, 5), (3, 2, 1, 4)]

        for idx_tuple in dihedral_indices:
            assert (
                idx_tuple in indices["dihedrals"]
                or (idx_tuple[-1], idx_tuple[-2], idx_tuple[-3], idx_tuple[-4])
                in indices["dihedrals"]
            )

        assert len(indices["impropers"]) == 4
        improper_indices = [
            (1, 0, 4, 5),
            (1, 0, 5, 2),
            (1, 0, 4, 2),
            (1, 4, 5, 2),
        ]

        for idx_tuple in improper_indices:
            assert all(
                idx_tuple[0] == members_tuple[0]
                for members_tuple in indices["impropers"]
            )

    def test_generate_pairs_list(self):
        # Methane with no 1-4 pair
        methane = mb.load("C", smiles=True)
        methane_top = from_mbuild(methane)
        methane_top.identify_connections()
        methane_pairs = generate_pairs_lists(
            methane_top, refer_from_scaling_factor=False
        )
        assert len(methane_pairs["pairs14"]) == len(methane_top.dihedrals) == 0

        # Ethane with 9 1-4 pairs
        ethane = mb.load("CC", smiles=True)
        ethane_top = from_mbuild(ethane)
        ethane_top.identify_connections()
        ethane_pairs = generate_pairs_lists(
            ethane_top, refer_from_scaling_factor=False
        )
        assert len(ethane_pairs["pairs14"]) == len(ethane_top.dihedrals) == 9

        # Cyclobutadiene with 16 dihedrals and 8 pairs (due to cyclic structure)
        cyclobutadiene = mb.load("C1=CC=C1", smiles=True)
        cyclobutadiene_top = from_mbuild(cyclobutadiene)
        cyclobutadiene_top.identify_connections()
        cyclobutadiene_top_pairs = generate_pairs_lists(
            cyclobutadiene_top, refer_from_scaling_factor=False
        )
        assert len(cyclobutadiene_top.dihedrals) == 16
        assert len(cyclobutadiene_top_pairs["pairs14"]) == 8

        # Cyclopentane with 45 dihedrals and 40 pairs (due to cyclic structure)
        cyclopentane = mb.load("C1CCCC1", smiles=True)
        cyclopentane_top = from_mbuild(cyclopentane)
        cyclopentane_top.identify_connections()
        cyclopentane_top_pairs = generate_pairs_lists(
            cyclopentane_top, refer_from_scaling_factor=False
        )

        assert len(cyclopentane_top.dihedrals) == 45
        assert len(cyclopentane_top_pairs["pairs14"]) == 40
