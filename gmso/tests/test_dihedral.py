import pytest
from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestDihedral(BaseTest):
    def test_dihedral_nonparametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom1")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4])

        assert connect.connection_type is None

    def test_dihedral_parametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        dihedral_type = DihedralType()

        connect = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
            dihedral_type=dihedral_type,
            name="dihedral_name",
        )

        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.dihedral_type is not None
        assert connect.name == "dihedral_name"

    def test_dihedral_fake(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(TypeError):
            Dihedral(connection_members=["fakeatom1", "fakeatom2", 4.2])

    def test_dihedral_fake_dihedraltype(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(ValidationError):
            Dihedral(
                connection_members=[atom1, atom2, atom3, atom4],
                dihedral_type="Fake dihedraltype",
            )

    def test_dihedral_constituent_types(self):
        atom1 = Atom(
            name="atom1", position=[0, 0, 0], atom_type=AtomType(name="A")
        )
        atom2 = Atom(
            name="atom2", position=[1, 0, 0], atom_type=AtomType(name="B")
        )
        atom3 = Atom(
            name="atom3", position=[1, 1, 0], atom_type=AtomType(name="C")
        )
        atom4 = Atom(
            name="atom4", position=[1, 1, 4], atom_type=AtomType(name="D")
        )
        dihtype = DihedralType(
            member_types=[
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
                atom4.atom_type.name,
            ]
        )
        dih = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )
        dih.dihedral_type = dihtype
        assert "A" in dih.connection_type.member_types
        assert "B" in dih.connection_type.member_types
        assert "C" in dih.connection_type.member_types
        assert "D" in dih.connection_type.member_types

    def test_dihedral_eq(self):
        atom1 = Atom(name="atom1", position=[0, 0, 0])
        atom2 = Atom(name="atom2", position=[1, 0, 0])
        atom3 = Atom(name="atom3", position=[1, 1, 0])
        atom4 = Atom(name="atom4", position=[1, 1, 1])

        ref_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        same_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        diff_dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        assert ref_dihedral != same_dihedral
        assert ref_dihedral != diff_dihedral

    def test_add_equivalent_connections(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        dihedral = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        dihedral_eq = Dihedral(connection_members=[atom4, atom3, atom2, atom1])
        dihedral_not_eq = Dihedral(
            connection_members=[atom4, atom2, atom3, atom1]
        )

        top = Topology()
        top.add_connection(dihedral)
        top.add_connection(dihedral_eq)
        assert top.n_dihedrals == 1

        top.add_connection(dihedral_not_eq)
        assert top.n_dihedrals == 2

    def test_equivalent_members_set(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        dihedral = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        dihedral_eq = Dihedral(connection_members=[atom4, atom3, atom2, atom1])
        dihedral_not_eq = Dihedral(
            connection_members=[atom4, atom2, atom3, atom1]
        )

        assert (
            tuple(dihedral_eq.connection_members)
            in dihedral.equivalent_members()
        )
        assert (
            tuple(dihedral.connection_members)
            in dihedral_eq.equivalent_members()
        )
        assert not (
            tuple(dihedral.connection_members)
            in dihedral_not_eq.equivalent_members()
        )

    def test_sort_dihedral_types(self):
        from gmso.utils.sorting import sort_by_classes, sort_by_types

        atom1 = Atom(
            name="atom1", position=[0, 0, 0], atom_type=AtomType(name="A")
        )
        atom2 = Atom(
            name="atom2", position=[1, 0, 0], atom_type=AtomType(name="B")
        )
        atom3 = Atom(
            name="atom3", position=[1, 1, 0], atom_type=AtomType(name="C")
        )
        atom4 = Atom(
            name="atom4", position=[1, 1, 4], atom_type=AtomType(name="D")
        )

        consituentList = [
            atom2.atom_type.name,
            atom4.atom_type.name,
            atom3.atom_type.name,
            atom1.atom_type.name,
        ]
        dihtype = DihedralType(
            member_types=consituentList, member_classes=consituentList
        )

        expected_sortingList = tuple(
            [atom.atom_type.name for atom in [atom1, atom3, atom4, atom2]]
        )
        assert sort_by_classes(dihtype) == expected_sortingList
        assert sort_by_types(dihtype) == expected_sortingList
