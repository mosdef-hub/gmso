import pytest
from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestImproper(BaseTest):
    def test_improper_nonparametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        connect = Improper(connection_members=[atom1, atom2, atom3, atom4])

        assert connect.connection_type is None

    def test_improper_parametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")

        improper_type = ImproperType()

        connect = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
            improper_type=improper_type,
            name="improper_name",
        )

        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.improper_type is not None
        assert connect.name == "improper_name"

    def test_improper_fake(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(TypeError):
            Improper(connection_members=["fakeatom1", "fakeatom2", 4.2])

    def test_improper_fake_impropertype(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        atom4 = Atom(name="atom4")
        with pytest.raises(ValidationError):
            Improper(
                connection_members=[atom1, atom2, atom3, atom4],
                connection_type="Fake impropertype",
            )

    def test_improper_constituent_types(self):
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
        imptype = ImproperType(
            member_types=[
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
                atom4.atom_type.name,
            ]
        )
        imp = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )
        imp.improper_type = imp.connection_type = imptype
        assert "A" in imp.connection_type.member_types
        assert "B" in imp.connection_type.member_types
        assert "C" in imp.connection_type.member_types
        assert "D" in imp.connection_type.member_types

    def test_improper_eq(self):
        atom1 = Atom(name="atom1", position=[0, 0, 0])
        atom2 = Atom(name="atom2", position=[1, 0, 0])
        atom3 = Atom(name="atom3", position=[1, 1, 0])
        atom4 = Atom(name="atom4", position=[1, 1, 1])

        ref_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        same_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        diff_improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4],
        )

        assert ref_improper != same_improper
        assert ref_improper != diff_improper

    def test_add_equivalent_connections(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        improper = Improper(connection_members=[atom1, atom2, atom3, atom4])
        improper_eq = Improper(connection_members=[atom1, atom3, atom2, atom4])
        improper_not_eq = Improper(
            connection_members=[atom2, atom3, atom1, atom4]
        )

        top = Topology()
        top.add_connection(improper)
        top.add_connection(improper_eq)
        assert top.n_impropers == 1
        top.add_connection(improper_not_eq)
        assert top.n_impropers == 2

    def test_equivalent_members_set(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")
        atom4 = Atom(name="AtomD")

        improper = Improper(connection_members=[atom1, atom2, atom3, atom4])
        improper_eq = Improper(connection_members=[atom1, atom3, atom2, atom4])
        improper_not_eq = Improper(
            connection_members=[atom2, atom3, atom1, atom4]
        )

        assert (
            tuple(improper_eq.connection_members)
            in improper.equivalent_members()
        )
        assert (
            tuple(improper.connection_members)
            in improper_eq.equivalent_members()
        )
        assert not (
            tuple(improper.connection_members)
            in improper_not_eq.equivalent_members()
        )

    def test_sort_improper_types(self):
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
        imptype = ImproperType(
            member_types=consituentList, member_classes=consituentList
        )

        expected_sortingList = tuple(
            [atom.atom_type.name for atom in [atom2, atom3, atom4, atom1]]
        )
        assert sort_by_classes(imptype) == expected_sortingList
        assert sort_by_types(imptype) == expected_sortingList
