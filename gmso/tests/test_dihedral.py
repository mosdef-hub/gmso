import pytest

from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.atom_type import AtomType
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestDihedral(BaseTest):
    def test_dihedral_nonparametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')

        connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        assert atom1 in connect.connection_members
        assert atom2 in connect.connection_members
        assert atom3 in connect.connection_members
        assert atom4 in connect.connection_members
        assert connect.connection_type is None

    def test_dihedral_parametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')

        dihedral_type = DihedralType()

        connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4],
                           connection_type=dihedral_type,
                           name='dihedral_name')

        assert len(connect.connection_members) == 4
        assert atom1 in connect.connection_members
        assert atom2 in connect.connection_members
        assert atom3 in connect.connection_members
        assert atom4 in connect.connection_members
        assert connect.connection_type is not None
        assert connect.name == 'dihedral_name'

    def test_dihedral_fake(self):
        with pytest.raises(GMSOError):
            Dihedral(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_dihedral_fake_dihedraltype(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        atom3 = Atom(name='atom3')
        atom4 = Atom(name='atom4')
        with pytest.raises(GMSOError):
            Dihedral(connection_members=[atom1, atom2, atom3, atom4],
                  connection_type='Fake dihedraltype')

    def test_dihedral_constituent_types(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0], atom_type=AtomType(name='A'))
        atom2 = Atom(name='atom2', position=[1, 0, 0], atom_type=AtomType(name='B'))
        atom3 = Atom(name='atom3', position=[1, 1, 0], atom_type=AtomType(name='C'))
        atom4 = Atom(name='atom4', position=[1, 1, 4], atom_type=AtomType(name='D'))
        dihtype = DihedralType(member_types=[atom1.atom_type.name,
                                             atom2.atom_type.name,
                                             atom3.atom_type.name,
                                             atom4.atom_type.name])
        dih = Dihedral(connection_members=[atom1, atom2, atom3, atom4],
                connection_type=dihtype)
        assert 'A' in dih.connection_type.member_types
        assert 'B' in dih.connection_type.member_types
        assert 'C' in dih.connection_type.member_types
        assert 'D' in dih.connection_type.member_types

    def test_dihedral_eq(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0])
        atom2 = Atom(name='atom2', position=[1, 0, 0])
        atom3 = Atom(name='atom3', position=[1, 1, 0])
        atom4 = Atom(name='atom4', position=[1, 1, 1])

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
