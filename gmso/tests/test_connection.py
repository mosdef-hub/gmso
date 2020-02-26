import pytest

from gmso.core.connection import Connection
from gmso.core.potential import Potential
from gmso.core.site import Site
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestConnection(BaseTest):
    def test_connection_nonparametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')

        assert site1.n_connections == 0
        assert site2.n_connections == 0

        connect = Connection(connection_members=[site1, site2])

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert connect.connection_type is None

    def test_connection_parametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        c_type = Potential()
        name = 'name'

        connect = Connection(connection_members=[site1, site2],
                             connection_type=c_type,
                             name=name)

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert len(connect.connection_members) == 2
        assert connect.connection_type is not None
        assert connect.name == name

    def test_connection_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(GMSOError):
            Connection(connection_members=['fakesite1', 'fakesite2'])

    def test_bond_fake_ctype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(GMSOError):
            Connection(connection_members=[site1, site2],
                       connection_type='Fake ctype',
                       name='fake')

    def test_connection_eq(self):
        site1 = Site(name='site1', position=[0, 0, 0])
        site2 = Site(name='site2', position=[1, 1, 1])

        ref_connection = Connection(
            connection_members=[site1, site2],
        )

        same_connection = Connection(
            connection_members=[site1, site2],
        )
        # Two connections are never equal
        assert ref_connection != same_connection

    def test_add_connection_same_sites(self):
        site1 = Site()
        site2 = Site()
        with pytest.raises(GMSOError):
            bond1 = Connection([site1, site1])
            angle1 = Connection([site1, site2, site1])
