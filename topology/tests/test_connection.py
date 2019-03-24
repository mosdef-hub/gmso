import pytest

from topology.core.connection import Connection
from topology.core.potential import Potential
from topology.core.site import Site
from topology.tests.base_test import BaseTest
from topology.exceptions import TopologyError


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

        connect = Connection(connection_members=[site1, site2], 
                connection_type=c_type)

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert len(connect.connection_members) == 2
        assert connect.connection_type is not None

    def test_connection_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(TopologyError):
            Connection(connection_members=['fakesite1', 'fakesite2'])

    def test_bond_fake_ctype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(TopologyError):
            Connection(connection_members=[site1,site2],
                       connection_type='Fake ctype')

