import pytest

from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.atom_type import AtomType
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestImproper(BaseTest):
    def test_improper_nonparametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')

        connect = Improper(connection_members=[atom1, atom2, atom3, atom4])

        assert atom1 in connect.connection_members
        assert atom2 in connect.connection_members
        assert atom3 in connect.connection_members
        assert atom4 in connect.connection_members

        assert connect.connection_type is None

    def test_improper_parametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')

        improper_type = ImproperType()

        connect = Improper(connection_members=[atom1, atom2, atom3, atom4],
                        connection_type=improper_type,
                        name='improper_name')

        assert len(connect.connection_members) == 4
        assert atom1 in connect.connection_members
        assert atom2 in connect.connection_members
        assert atom3 in connect.connection_members
        assert atom4 in connect.connection_members
        assert connect.connection_type is not None
        assert connect.name == 'improper_name'

    def test_improper_fake(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')
        with pytest.raises(GMSOError):
            Improper(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_improper_fake_impropertype(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')
        with pytest.raises(GMSOError):
            Improper(connection_members=[atom1, atom2, atom3, atom4],
                  connection_type='Fake impropertype')

    def test_improper_constituent_types(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0], atom_type=AtomType(name='A'))
        atom2 = Atom(name='atom2', position=[1, 0, 0], atom_type=AtomType(name='B'))
        atom3 = Atom(name='atom3', position=[1, 1, 0], atom_type=AtomType(name='C'))
        atom4 = Atom(name='atom4', position=[1, 1, 4], atom_type=AtomType(name='D'))
        imptype = ImproperType(member_types=[atom1.atom_type.name,
                                             atom2.atom_type.name,
                                             atom3.atom_type.name,
                                             atom4.atom_type.name])
        imp = Improper(connection_members=[atom1, atom2, atom3, atom4],
                connection_type=imptype)
        assert 'A' in imp.connection_type.member_types
        assert 'B' in imp.connection_type.member_types
        assert 'C' in imp.connection_type.member_types
        assert 'D' in imp.connection_type.member_types

    def test_improper_eq(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0])
        atom2 = Atom(name='atom2', position=[1, 0, 0])
        atom3 = Atom(name='atom3', position=[1, 1, 0])
        atom4 = Atom(name='atom4', position=[1, 1, 1])

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
