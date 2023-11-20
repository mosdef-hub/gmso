import pytest
from pydantic import ValidationError

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.topology import Topology
from gmso.tests.base_test import BaseTest


class TestAngle(BaseTest):
    def test_angle_nonparametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")

        connect = Angle(connection_members=[atom1, atom2, atom3])
        assert connect.angle_type is None

    def test_angle_parametrized(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")

        angle_type = AngleType()

        connect = Angle(
            connection_members=[atom1, atom2, atom3],
            angle_type=angle_type,
            name="angle_name",
        )

        assert len(connect.connection_members) == 3
        assert connect.angle_type is not None
        assert connect.name == "angle_name"

    def test_angle_fake(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        with pytest.raises(TypeError):
            Angle(connection_members=["fakesite1", "fakesite2", 4.2])

    def test_angle_fake_angletype(self):
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        atom3 = Atom(name="atom3")
        with pytest.raises(ValidationError):
            Angle(
                connection_members=[atom1, atom2, atom3],
                angle_type="Fake angletype",
            )

    def test_angle_constituent_types(self):
        atom1 = Atom(
            name="atom1", position=[0, 0, 0], atom_type=AtomType(name="A")
        )
        atom2 = Atom(
            name="atom2", position=[1, 0, 0], atom_type=AtomType(name="B")
        )
        atom3 = Atom(
            name="atom3", position=[1, 1, 0], atom_type=AtomType(name="C")
        )
        angtype = AngleType(
            member_types=[
                atom1.atom_type.name,
                atom2.atom_type.name,
                atom3.atom_type.name,
            ]
        )
        ang = Angle(
            connection_members=[atom1, atom2, atom3], angle_type=angtype
        )
        assert "A" in ang.angle_type.member_types
        assert "B" in ang.angle_type.member_types
        assert "C" in ang.angle_type.member_types

    def test_angle_eq(self):
        atom1 = Atom(name="atom1", position=[0, 0, 0])
        atom2 = Atom(name="atom2", position=[1, 1, 1])
        atom3 = Atom(name="atom3", position=[1, 1, 1])

        ref_angle = Angle(
            connection_members=[atom1, atom2, atom3],
        )

        same_angle = Angle(
            connection_members=[atom1, atom2, atom3],
        )

        diff_angle = Angle(
            connection_members=[atom3, atom2, atom1],
        )

        assert ref_angle != same_angle
        assert ref_angle != diff_angle

    def test_add_equivalent_connections(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")

        angle = Angle(connection_members=[atom1, atom2, atom3])
        angle_eq = Angle(connection_members=[atom3, atom2, atom1])
        angle_not_eq = Angle(connection_members=[atom1, atom3, atom2])

        top = Topology()
        top.add_connection(angle)
        top.add_connection(angle_eq)
        assert top.n_angles == 1

        top.add_connection(angle_not_eq)
        assert top.n_angles == 2

    def test_equivalent_members_set(self):
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        atom3 = Atom(name="AtomC")

        angle = Angle(connection_members=[atom1, atom2, atom3])
        angle_eq = Angle(connection_members=[atom3, atom2, atom1])
        angle_not_eq = Angle(connection_members=[atom1, atom3, atom2])

        assert tuple(angle_eq.connection_members) in angle.equivalent_members()
        assert tuple(angle.connection_members) in angle_eq.equivalent_members()
        assert not (
            tuple(angle.connection_members) in angle_not_eq.equivalent_members()
        )
