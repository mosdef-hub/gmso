import foyer
import mbuild as mb
import pytest

from gmso.external.convert_mbuild import from_mbuild
from gmso.external.convert_parmed import from_parmed
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
    def usable_compound(self):
        return mb.fill_box(mb.lib.molecules.Ethane(), n_compounds=20, density=2)

    @pytest.fixture(scope="session")
    def top_from_mbuild(self, usable_compound):
        top = from_mbuild(usable_compound)
        identify_connections(top)
        return top

    @pytest.fixture(scope="session")
    def top_from_parmed(self, usable_compound):
        pmd_structure = usable_compound.to_parmed()
        top = from_parmed(pmd_structure)
        identify_connections(top)
        return top

    @pytest.fixture(scope="session")
    def top_from_typed_parmed(self, usable_compound):
        ff = foyer.Forcefield(name="oplsaa")
        pmd_structure = ff.apply(usable_compound)
        top = from_parmed(pmd_structure)
        return top

    def test_no_boundary_bonds(self, benzene_ua_box):
        benzene_ua_box.sites[0].molecule = benzene_ua_box.sites[6].molecule
        with pytest.raises(AssertionError):
            assert_no_boundary_bonds(benzene_ua_box)

    def test_no_boundary_bonds_ethane_box(self, top_from_mbuild):
        assert_no_boundary_bonds(top_from_mbuild)

    def test_molecule_bonds(self, top_from_mbuild):
        for molecule in top_from_mbuild.unique_site_labels("molecule"):
            bonds = list(molecule_bonds(top_from_mbuild, molecule))
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

    def test_molecule_angles(self, top_from_mbuild):
        for molecule in top_from_mbuild.unique_site_labels("molecule"):
            angles = list(molecule_angles(top_from_mbuild, molecule))
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

    def test_molecule_dihedrals(self, top_from_mbuild):
        for molecule in top_from_mbuild.unique_site_labels("molecule"):
            dihedrals = list(molecule_dihedrals(top_from_mbuild, molecule))
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

    def test_molecule_impropers(self, top_from_mbuild):
        for molecule in top_from_mbuild.unique_site_labels("molecule"):
            impropers = list(molecule_impropers(top_from_mbuild, molecule))
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

    def test_molecule_numbers(
        self, top_from_mbuild, top_from_parmed, top_from_typed_parmed
    ):
        assert top_from_mbuild.sites[0].molecule.number == 0
        assert top_from_parmed.sites[0].molecule.number == 0
        assert top_from_typed_parmed.sites[0].molecule.number == 0
