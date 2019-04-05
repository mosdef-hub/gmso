import pytest

from topology.core.angle import Angle
from topology.core.angle_type import AngleType
from topology.core.site import Site
from topology.tests.base_test import BaseTest
from topology.exceptions import TopologyError


class TestAngle(BaseTest):
    def test_angle_nonparametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0

        connect = Angle(connection_members=[site1, site2, site3])

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert connect.connection_type is None

    def test_angle_parametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0
        angle_type = AngleType()

        connect = Angle(connection_members=[site1, site2, site3],
                        connection_type=angle_type,
                        connection_name='angle_name')

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert len(connect.connection_members) == 3
        assert connect.connection_type is not None
        assert connect.connection_name is not None

    def test_angle_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        with pytest.raises(TopologyError):
            Angle(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_angle_fake_angletype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        with pytest.raises(TopologyError):
            Angle(connection_members=[site1, site2, site3],
                  connection_type='Fake angletype')

