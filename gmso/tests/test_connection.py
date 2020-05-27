import pytest

from gmso.core.connection import Connection
from gmso.core.potential import Potential
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestConnection(BaseTest):
    def test_connection_nonparametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')

        connect = Connection(connection_members=[atom1, atom2])

        assert connect.connection_type is None

    def test_connection_parametrized(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')

        c_type = Potential()
        name = 'name'

        connect = Connection(connection_members=[atom1, atom2],
                             connection_type=c_type,
                             name=name)

        assert len(connect.connection_members) == 2
        assert connect.connection_type is not None
        assert connect.name == name

    def test_connection_fake(self):
        with pytest.raises(GMSOError):
            Connection(connection_members=['fakesite1', 'fakesite2'])

    def test_bond_fake_ctype(self):
        atom1 = Atom(name='atom1')
        atom2 = Atom(name='atom2')
        with pytest.raises(GMSOError):
            Connection(connection_members=[atom1, atom2],
                       connection_type='Fake ctype',
                       name='fake')

    def test_connection_eq(self):
        atom1 = Atom(name='atom1', position=[0, 0, 0])
        atom2 = Atom(name='atom2', position=[1, 1, 1])

        ref_connection = Connection(
            connection_members=[atom1, atom2],
        )

        same_connection = Connection(
            connection_members=[atom1, atom2],
        )
        # Two connections are never equal
        assert ref_connection != same_connection

    def test_add_connection_same_sites(self):
        atom1 = Atom()
        atom2 = Atom()
        with pytest.raises(GMSOError):
            bond1 = Connection([atom1, atom1])
            angle1 = Connection([atom1, atom2, atom1])
