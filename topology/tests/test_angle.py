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
                        name='angle_name')

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert len(connect.connection_members) == 3
        assert connect.connection_type is not None
        assert connect.name == "angle_name"

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

    def test_angle_eq(self):
        site1 = Site(name='site1', position=[0, 0, 0])
        site2 = Site(name='site2', position=[1, 1, 1])
        site3 = Site(name='site3', position=[1, 1, 1])

        ref_angle = Angle(
            connection_members=[site1, site2, site3],
        )

        same_angle = Angle(
            connection_members=[site1, site2, site3],
        )

        diff_angle = Angle(
            connection_members=[site2, site2, site1],
        )

        assert ref_angle == same_angle
        assert ref_angle != diff_angle
