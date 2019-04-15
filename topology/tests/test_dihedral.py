import pytest

from topology.core.dihedral import Dihedral
from topology.core.dihedral_type import DihedralType
from topology.core.site import Site
from topology.tests.base_test import BaseTest
from topology.exceptions import TopologyError


class TestDihedral(BaseTest):
    def test_dihedral_nonparametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0
        assert site4.n_connections == 0

        connect = Dihedral(connection_members=[site1, site2, site3, site4])

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert site4.n_connections == 1
        assert connect.connection_type is None

    def test_dihedral_parametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0
        assert site4.n_connections == 0
        dihedral_type = DihedralType()

        connect = Dihedral(connection_members=[site1, site2, site3, site4],
                        connection_type=dihedral_type)

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert site4.n_connections == 1
        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None

    def test_dihedral_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(TopologyError):
            Dihedral(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_dihedral_fake_dihedraltype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(TopologyError):
            Dihedral(connection_members=[site1, site2, site3, site4],
                  connection_type='Fake dihedraltype')

