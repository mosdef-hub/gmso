import mbuild as mb
import pytest

from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization.subtopology_utils import (
    _members_in_subtop,
    assert_no_boundary_bonds,
    subtop_angles,
    subtop_bonds,
    subtop_dihedrals,
    subtop_impropers,
)
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)
from gmso.utils.connectivity import identify_connections


class TestSubTopologyUtils(ParameterizationBaseTest):
    @pytest.fixture(scope="session")
    def ethane_box_gmso(self):
        ethane_box = mb.fill_box(
            mb.lib.molecules.Ethane(), n_compounds=20, density=2
        )
        ethane_box_gmso = from_mbuild(ethane_box)
        identify_connections(ethane_box_gmso)
        return ethane_box_gmso

    def test_no_boundary_bonds_ethane(self, ethane):
        with pytest.raises(AssertionError):
            assert_no_boundary_bonds(ethane.subtops[0])

    def test_no_boundary_bonds_ethane_box(self, ethane_box_gmso):
        for subtop in ethane_box_gmso.subtops:
            assert_no_boundary_bonds(subtop)

    def test_subtopology_bonds(self, ethane_box_gmso):
        for subtop in ethane_box_gmso.subtops:
            bonds = list(subtop_bonds(subtop))
            assert len(bonds) == 7
            for bond in bonds:
                assert _members_in_subtop(bond, subtop)

            bond_members = map(
                lambda b: tuple(map(lambda s: s.name, b.connection_members)),
                bonds,
            )
            expected_members = {("C", "H"), ("C", "C"), ("H", "C")}
            assert all(
                b_member in expected_members for b_member in bond_members
            )

    def test_subtopology_angles(self, ethane_box_gmso):
        for subtop in ethane_box_gmso.subtops:
            angles = list(subtop_angles(subtop))
            assert len(list(angles)) == 12
            for angle in angles:
                assert _members_in_subtop(angle, subtop)

            angle_members = map(
                lambda a: tuple(map(lambda s: s.name, a.connection_members)),
                angles,
            )
            expected_members = {
                ("H", "C", "H"),
                ("H", "C", "C"),
                ("C", "C", "H"),
            }
            assert all(
                a_member in expected_members for a_member in angle_members
            )

    def test_subtopology_dihedrals(self, ethane_box_gmso):
        for subtop in ethane_box_gmso.subtops:
            dihedrals = list(subtop_dihedrals(subtop))
            assert len(dihedrals) == 9
            for dihedral in dihedrals:
                assert _members_in_subtop(dihedral, subtop)

            dihedral_members = map(
                lambda d: tuple(map(lambda s: s.name, d.connection_members)),
                dihedrals,
            )
            expected_members = {("H", "C", "C", "H")}
            assert all(
                a_member in expected_members for a_member in dihedral_members
            )

    def test_subtopology_impropers(self, ethane_box_gmso):
        for subtop in ethane_box_gmso.subtops:
            impropers = list(subtop_impropers(subtop))
            assert len(impropers) == 8
            for improper in impropers:
                assert _members_in_subtop(improper, subtop)

            improper_members = list(
                map(
                    lambda i: tuple(
                        map(lambda s: s.name, i.connection_members)
                    ),
                    impropers,
                )
            )
            expected_members = {
                ("C", "C", "H", "H"),
                ("C", "H", "H", "C"),
                ("C", "H", "H", "H"),
                ("C", "H", "C", "H"),
            }
            assert all(
                a_member in expected_members for a_member in improper_members
            )
