import json

import pytest
import unyt as u

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.element import element_by_symbol
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.topology import Topology
from gmso.formats.formats_registry import UnsupportedFileFormatError
from gmso.tests.base_test import BaseTest


class TestSerialization(BaseTest):
    @pytest.fixture(scope="module")
    def full_atom_type(self):
        return AtomType(
            name="test_atom_type",
            expression="a*b+c*d+e**2",
            independent_variables={"a"},
            parameters={
                "b": 2.0 * u.amu,
                "c": 3.0 * u.nm / u.kg**2,
                "d": 5.0 * u.kJ / u.mol,
                "e": 1.0 * u.C,
            },
            mass=1.0 * u.amu,
            charge=1.0 * u.elementary_charge,
            atomclass="test_atom_class",
            doi="https://dx.doi.org/110.200.300",
            overrides={"A", "B", "C"},
            definition="CX_6",
            description="A test AtomType object",
            tags={"tag1": 10, "tag2": 10 * u.nm},
        )

    def test_atom_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        atoms_to_test = typed_ethane.sites
        for atom in atoms_to_test:
            atom_json = atom.model_dump_json()
            atom_copy = Atom.model_validate(json.loads(atom_json))
            assert are_equivalent_atoms(atom, atom_copy)

    def test_atom_types_to_json_loop(self, typed_ethane):
        atom_types_to_test = typed_ethane.atom_types
        for atom_type in atom_types_to_test:
            atom_type_json = atom_type.model_dump_json()
            atom_type_copy = AtomType.model_validate(json.loads(atom_type_json))
            assert atom_type_copy == atom_type

    def test_bond_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for bond in typed_ethane.bonds:
            bond_json = bond.model_dump_json()
            bond_copy = Bond.model_validate(json.loads(bond_json))
            assert bond_copy.name == bond.name
            for member1, member2 in zip(
                bond.connection_members, bond_copy.connection_members
            ):
                assert are_equivalent_atoms(member1, member2)
            assert bond_copy.bond_type == bond.bond_type

    def test_bond_type_to_json_loop(self, typed_ethane):
        bond_types_to_test = typed_ethane.bond_types
        for bond_type in bond_types_to_test:
            bond_type_json = bond_type.model_dump_json()
            bond_type_copy = BondType.model_validate(json.loads(bond_type_json))
            assert bond_type_copy == bond_type

    def test_angle_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for angle in typed_ethane.angles:
            angle_json = angle.model_dump_json()
            angle_copy = Angle.model_validate(json.loads(angle_json))
            for member1, member2 in zip(
                angle.connection_members, angle_copy.connection_members
            ):
                assert are_equivalent_atoms(member1, member2)
            assert angle.angle_type == angle_copy.angle_type

    def test_angle_type_to_json_loop(self, typed_ethane):
        angle_types_to_test = typed_ethane.angle_types
        for angle_type in angle_types_to_test:
            angle_type_json = angle_type.model_dump_json()
            angle_type_copy = AngleType.model_validate(
                json.loads(angle_type_json)
            )
            assert angle_type_copy == angle_type

    def test_dihedral_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for dihedral in typed_ethane.dihedrals:
            dihedral_json = dihedral.model_dump_json()
            dihedral_copy = Dihedral.model_validate(json.loads(dihedral_json))
            for member1, member2 in zip(
                dihedral.connection_members, dihedral_copy.connection_members
            ):
                assert are_equivalent_atoms(member1, member2)
            assert dihedral.dihedral_type == dihedral_copy.dihedral_type

    def test_dihedral_types_to_json_loop(self, typed_ethane):
        dihedral_types_to_test = typed_ethane.dihedral_types
        for dihedral_type in dihedral_types_to_test:
            dihedral_type_json = dihedral_type.model_dump_json()
            dihedral_type_copy = DihedralType.model_validate(
                json.loads(dihedral_type_json)
            )
            assert dihedral_type_copy == dihedral_type

    def test_improper_to_json_loop(self, typed_ethane, are_equivalent_atoms):
        for improper in typed_ethane.impropers:
            improper_json = improper.model_dump_json()
            improper_copy = Improper.model_validate(json.loads(improper_json))
            for member1, member2 in zip(
                improper_copy.connection_members, improper.connection_members
            ):
                assert are_equivalent_atoms(member1, member2)
            assert improper_copy.improper_type == improper.improper_type

    def test_improper_types_to_json_loop(self, typed_ethane):
        improper_types_to_test = typed_ethane.improper_types
        for improper_type in improper_types_to_test:
            improper_type_json = improper_type.model_dump_json()
            improper_type_copy = ImproperType.model_validate(
                json.loads(improper_type_json)
            )
            improper_type_copy.topology = improper_type.topology
            assert improper_type_copy == improper_type

    def test_atom_every_field_set(self, full_atom_type, are_equivalent_atoms):
        atom = Atom(
            name="test_atom",
            label="test_label",
            position=[0.0, 0.0, 0.0],
            charge=1.5,
            mass=2.0,
            element=element_by_symbol("C"),
            atom_type=full_atom_type,
        )

        atom_copy = Atom.model_validate(atom.json_dict())
        assert are_equivalent_atoms(atom, atom_copy)

    def test_bond_every_field_set(self, full_atom_type, are_equivalent_atoms):
        atom1 = Atom(
            name="test_atom1",
            label="test_label1",
            position=[0.0, 0.0, 0.0],
            charge=1.5,
            mass=2.0,
            element=element_by_symbol("C"),
            atom_type=full_atom_type,
        )

        atom2 = Atom(
            name="test_atom1",
            label="test_label1",
            position=[0.1, 0.4, 0.5],
            charge=5,
            mass=2.6,
            element=element_by_symbol("H"),
            atom_type=full_atom_type,
        )

        bond = Bond(name="test_bond1", connection_members=(atom1, atom2))

        bond_type = BondType(
            name="test_bond_type",
            expression="a*b+c**2",
            parameters={"a": 10 * u.nm, "b": 20 * u.angstrom},
            independent_variables={"c"},
        )

        bond_copy = Bond.model_validate(bond.json_dict())
        assert bond_copy.name == bond.name
        for member1, member2 in zip(
            bond.connection_members, bond_copy.connection_members
        ):
            assert are_equivalent_atoms(member1, member2)
        assert bond_copy.bond_type == bond.bond_type

    def test_include_and_exclude(self):
        atom = Atom(mass=2.0 * u.g / u.mol, charge=30.0 * u.C, name="TestAtom")
        atom_json = atom.model_dump_json(exclude={"mass"})
        assert "mass" not in atom_json
        atom_json = atom.model_dump_json(exclude={"mass_"})
        assert "mass" not in atom_json
        atom_json = atom.model_dump_json(include={"mass"})
        assert "name" not in atom_json
        atom_json = atom.model_dump_json(include={"mass_"})
        assert "name" not in atom_json

    def test_full_serialization(
        self,
        typed_ethane,
        are_equivalent_topologies,
    ):
        typed_ethane.save("eth.json", types=True)
        typed_ethane_copy = Topology.load("eth.json")
        assert are_equivalent_topologies(typed_ethane_copy, typed_ethane)

    def test_serialization_with_box(
        self, n_typed_xe_mie, are_equivalent_topologies
    ):
        top = n_typed_xe_mie(n_sites=20)
        top.save("n_typed_xe_mie_20.json")
        top_copy = Topology.load("n_typed_xe_mie_20.json")
        assert are_equivalent_topologies(top, top_copy)

    def test_serialization_with_pairpotential_types(
        self, pairpotentialtype_top, are_equivalent_topologies
    ):
        pairpotentialtype_top.save("pptype.json", types=True)
        pptop_copy = Topology.load("pptype.json")
        assert are_equivalent_topologies(pptop_copy, pairpotentialtype_top)

    def test_serialization_unsupported_file_format(self, ethane_from_scratch):
        with pytest.raises(UnsupportedFileFormatError):
            ethane_from_scratch.save("ethane_from_scratch.zip")

    def test_serialization_untyped_with_types_info(self, ethane_from_scratch):
        with pytest.warns(UserWarning):
            ethane_from_scratch.save("ethane_from_scratch.json", types=True)

    def test_serialization_overwrite(self, ethane_from_scratch):
        ethane_from_scratch.save("ethane_from_scratch.json", overwrite=False)
        with pytest.raises(FileExistsError):
            ethane_from_scratch.save(
                "ethane_from_scratch.json", overwrite=False
            )
        ethane_from_scratch.save("ethane_from_scratch.json", overwrite=True)
