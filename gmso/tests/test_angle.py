import pytest

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom_type import AtomType
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestAngle(BaseTest):
    def test_angle_nonparametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')

        connect = Angle(connection_members=[atom1, atom2, atom3])
        assert atom1 in connect.connection_members
        assert atom2 in connect.connection_members
        assert atom3 in connect.connection_members
        assert connect.connection_type is None

    def test_angle_parametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')

        angle_type = AngleType()

        connect = Angle(connection_members=[atom1, atom2, atom3],
                        connection_type=angle_type,
                        name='angle_name')

        assert len(connect.connection_members) == 3
        assert connect.connection_type is not None
        assert connect.name == 'angle_name'

    def test_angle_fake(self):
        with pytest.raises(GMSOError):
            Angle(connection_members=['fakeatom1', 'fakeatom2', 4.2])

    def test_angle_fake_angletype(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        with pytest.raises(GMSOError):
            Angle(connection_members=[atom1, atom2, atom3],
                  connection_type='Fake angletype')

    def test_angle_constituent_types(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0], atom_type=AtomType(name='A'))
        atom2 = Atom(name='atom2', position=[1, 0, 0], atom_type=AtomType(name='B'))
        atom3 = Atom(name='atom3', position=[1, 1, 0], atom_type=AtomType(name='C'))
        angtype = AngleType(member_types=[atom1.atom_type.name, atom2.atom_type.name,
            atom3.atom_type.name])
        ang = Angle(connection_members=[atom1, atom2, atom3], 
                connection_type=angtype)
        assert 'A' in ang.connection_type.member_types
        assert 'B' in ang.connection_type.member_types
        assert 'C' in ang.connection_type.member_types

    def test_angle_eq(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0])
        atom2 = Atom(name='atom2', position=[1, 1, 1])
        atom3 = Atom(name='atom3', position=[1, 1, 1])

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
