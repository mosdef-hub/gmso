from topology.core.connection import Connection
from topology.core.site import Site
from topology.tests.base_test import BaseTest


class TestConnection(BaseTest):
    def test_new_connect(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')

        assert site1.n_connections == 0
        assert site2.n_connections == 0

        connect = Connection(site1=site1, site2=site2)

        assert site1.n_connections == 1
        assert site2.n_connections == 1
