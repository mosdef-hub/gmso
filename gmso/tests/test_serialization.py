from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.angle import Angle
from gmso.core.dihedral import Dihedral
from gmso.core.improper import Improper
from gmso.core.atom_type import AtomType
from gmso.core.bond_type import BondType
from gmso.core.angle_type import AngleType
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper_type import ImproperType
from gmso.tests.base_test import BaseTest


class TestSerialization(BaseTest):

    def test_atom_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        atoms_to_test = typed_ethane.sites
        for atom in atoms_to_test:
            atom_json = atom.json()
            atom_copy = Atom.parse_raw(atom_json)
            assert are_equivalent_atoms(atom, atom_copy)

    def test_atom_types_to_json_loop(self, typed_ethane):
        atom_types_to_test = typed_ethane.atom_types
        for atom_type in atom_types_to_test:
            atom_type_json = atom_type.json()
            atom_type_copy = AtomType.parse_raw(atom_type_json)
            atom_type_copy.topology = atom_type.topology
            assert atom_type_copy == atom_type

    def test_bond_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for bond in typed_ethane.bonds:
            bond_json = bond.json()
            bond_copy = Bond.parse_raw(bond_json)
            assert bond_copy.name == bond.name
            for member1, member2 in zip(bond.connection_members, bond_copy.connection_members):
                assert are_equivalent_atoms(member1, member2)
            assert bond_copy.bond_type == bond.bond_type

    def test_bond_type_to_json_loop(self, typed_ethane):
        bond_types_to_test = typed_ethane.bond_types
        for bond_type in bond_types_to_test:
            bond_type_json = bond_type.json()
            bond_type_copy = BondType.parse_raw(bond_type_json)
            bond_type_copy.topology = bond_type.topology
            assert bond_type_copy == bond_type

    def test_angle_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for angle in typed_ethane.angles:
            angle_json = angle.json()
            angle_copy = Angle.parse_raw(angle_json)
            for member1, member2 in zip(angle.connection_members, angle_copy.connection_members):
                assert are_equivalent_atoms(member1, member2)
            assert angle.angle_type == angle_copy.angle_type

    def test_angle_type_to_json_loop(self, typed_ethane):
        angle_types_to_test = typed_ethane.angle_types
        for angle_type in angle_types_to_test:
            angle_type_json = angle_type.json()
            angle_type_copy = AngleType.parse_raw(angle_type_json)
            angle_type_copy.topology = angle_type.topology
            assert angle_type_copy == angle_type

    def test_dihedral_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for dihedral in typed_ethane.dihedrals:
            dihedral_json = dihedral.json()
            dihedral_copy = Dihedral.parse_raw(dihedral_json)
            for member1, member2 in zip(dihedral.connection_members, dihedral_copy.connection_members):
                assert are_equivalent_atoms(member1, member2)
            assert dihedral.dihedral_type == dihedral_copy.dihedral_type

    def test_dihedral_types_to_json_loop(self, typed_ethane):
        dihedral_types_to_test = typed_ethane.dihedral_types
        for dihedral_type in dihedral_types_to_test:
            dihedral_type_json = dihedral_type.json()
            dihedral_type_copy = DihedralType.parse_raw(dihedral_type_json)
            dihedral_type_copy.topology = dihedral_type.topology
            assert dihedral_type_copy == dihedral_type

    def test_improper_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for improper in typed_ethane.impropers:
            improper_json = improper.json()
            improper_copy = Improper.parse_raw(improper_json)
            for member1, member2 in zip(improper_copy.connection_members, improper.connection_members):
                assert are_equivalent_atoms(member1, member2)
            assert improper_copy.improper_type == improper.improper_type

    def test_improper_types_to_json_loop(self, typed_ethane):
        improper_types_to_test = typed_ethane.improper_types
        for improper_type in improper_types_to_test:
            improper_type_json = improper_type.json()
            improper_type_copy = ImproperType.parse_raw(improper_type_json)
            improper_type_copy.topology = improper_type.topology
            assert improper_type_copy == improper_type
