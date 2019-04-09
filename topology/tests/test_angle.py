import pytest

from topology.core.angle import Angle
from topology.core.angle_type import AngleType
from topology.core.atom_type import AtomType
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

    def test_angle_constituent_types(self):
        site1 = Site(name='site1', position=[0,0,0], atom_type=AtomType(name='A'))
        site2 = Site(name='site2', position=[1,0,0], atom_type=AtomType(name='B'))
        site3 = Site(name='site3', position=[1,1,0], atom_type=AtomType(name='C'))
        angtype = AngleType(types=[site1.atom_type.name, site2.atom_type.name,
            site3.atom_type.name])
        ang = Angle(connection_members=[site1, site2,site3], 
                connection_type=angtype)
        assert 'A' in ang.connection_type.types
        assert 'B' in ang.connection_type.types
        assert 'C' in ang.connection_type.types

