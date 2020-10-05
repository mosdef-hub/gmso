from gmso.core.atom_type import AtomType
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest


class TestSerialization(BaseTest):

    def test_atom_to_json_loop(self, typed_ethane):
        atoms_to_test = typed_ethane.sites
        for atom in atoms_to_test:
            atom_json = atom.json()
            print(atom_json)
            atom_copy = Atom.parse_raw(atom_json)
            assert atom_copy.name == atom.name
            assert atom_copy.label == atom.label
            assert atom_copy.charge == atom.charge
            assert atom_copy.element == atom.element
            assert atom_copy.atom_type == atom.atom_type

    def test_atom_types_to_json_loop(self, typed_ethane):
        atom_types_to_test = typed_ethane.atom_types
        for atom_type in atom_types_to_test:
            atom_type_json = atom_type.json()
            atom_type_copy = AtomType.parse_raw(atom_type_json)
            atom_type_copy.topology = atom_type.topology
            assert atom_type_copy == atom_type