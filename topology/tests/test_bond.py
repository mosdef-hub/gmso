import pytest

from topology.core.bond import Bond
from topology.core.bond_type import BondType
from topology.core.site import Site
from topology.tests.base_test import BaseTest
from topology.exceptions import TopologyError


class TestBond(BaseTest):
    def test_bond_nonparametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')

        assert site1.n_connections == 0
        assert site2.n_connections == 0

        connect = Bond(connection_members=[site1, site2])

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert connect.connection_type is None

    def test_bond_parametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        bond_type = BondType()

        connect = Bond(connection_members=[site1, site2],
                connection_type=bond_type)

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert len(connect.connection_members) == 2
        assert connect.connection_type is not None

    def test_bond_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(TopologyError):
            Bond(connection_members=['fakesite1', 'fakesite2'])

    def test_bond_fake_bondtype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        with pytest.raises(TopologyError):
            Bond(connection_members=[site1, site2],
                 connection_type='Fake bondtype')

