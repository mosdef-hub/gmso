import mbuild as mb
import pytest

from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization.molecule_utils import (
    _conn_in_molecule,
    assert_no_boundary_bonds,
    molecule_angles,
    molecule_bonds,
    molecule_dihedrals,
    molecule_impropers,
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

    def test_no_boundary_bonds(self, benzene_ua_box):
        benzene_ua_box.sites[0].molecule = benzene_ua_box.sites[6].molecule
        with pytest.raises(AssertionError):
            assert_no_boundary_bonds(benzene_ua_box)

    def test_no_boundary_bonds_ethane_box(self, ethane_box_gmso):
        assert_no_boundary_bonds(ethane_box_gmso)

    def test_molecule_bonds(self, ethane_box_gmso):
        for molecule in ethane_box_gmso.unique_site_labels("molecule"):
            bonds = list(molecule_bonds(ethane_box_gmso, molecule))
            assert len(bonds) == 7
            for bond in bonds:
                assert _conn_in_molecule(bond, molecule)

            bond_members = map(
                lambda b: tuple(map(lambda s: s.name, b.connection_members)),
                bonds,
            )
            expected_members = {("C", "H"), ("C", "C"), ("H", "C")}
            assert all(
                b_member in expected_members for b_member in bond_members
            )

    def test_molecule_angles(self, ethane_box_gmso):
        for molecule in ethane_box_gmso.unique_site_labels("molecule"):
            angles = list(molecule_angles(ethane_box_gmso, molecule))
            assert len(list(angles)) == 12
            for angle in angles:
                assert _conn_in_molecule(angle, molecule)

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

    def test_molecule_dihedrals(self, ethane_box_gmso):
        for molecule in ethane_box_gmso.unique_site_labels("molecule"):
            dihedrals = list(molecule_dihedrals(ethane_box_gmso, molecule))
            assert len(dihedrals) == 9
            for dihedral in dihedrals:
                assert _conn_in_molecule(dihedral, molecule)

            dihedral_members = map(
                lambda d: tuple(map(lambda s: s.name, d.connection_members)),
                dihedrals,
            )
            expected_members = {("H", "C", "C", "H")}
            assert all(
                a_member in expected_members for a_member in dihedral_members
            )

    def test_molecule_impropers(self, ethane_box_gmso):
        for molecule in ethane_box_gmso.unique_site_labels("molecule"):
            impropers = list(molecule_impropers(ethane_box_gmso, molecule))
            assert len(impropers) == 8
            for improper in impropers:
                assert _conn_in_molecule(improper, molecule)

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
