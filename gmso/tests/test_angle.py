import pytest

from gmso.core.topology import Topology
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom_type import AtomType
from gmso.core.site import Site
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


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
        assert connect.name == 'angle_name'

    def test_angle_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        with pytest.raises(GMSOError):
            Angle(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_angle_fake_angletype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        with pytest.raises(GMSOError):
            Angle(connection_members=[site1, site2, site3],
                  connection_type='Fake angletype')

    def test_angle_constituent_types(self):
        site1 = Site(name='site1', position=[0,0,0], atom_type=AtomType(name='A'))
        site2 = Site(name='site2', position=[1,0,0], atom_type=AtomType(name='B'))
        site3 = Site(name='site3', position=[1,1,0], atom_type=AtomType(name='C'))
        angtype = AngleType(member_types=[site1.atom_type.name, site2.atom_type.name,
            site3.atom_type.name])
        ang = Angle(connection_members=[site1, site2, site3],
                connection_type=angtype)
        assert 'A' in ang.connection_type.member_types
        assert 'B' in ang.connection_type.member_types
        assert 'C' in ang.connection_type.member_types

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
            connection_members=[site3, site2, site1],
        )

        assert ref_angle != same_angle
        assert ref_angle != diff_angle

    def test_add_equivalent_connections(self):
        site1 = Site(name="SiteA")
        site2 = Site(name="SiteB")
        site3 = Site(name="SiteC")

        angle = Angle([site1, site2, site3])
        angle_eq = Angle([site3, site2, site1])
        angle_not_eq = Angle([site1, site3, site2])

        top = Topology()
        top.add_connection(angle)
        top.add_connection(angle_eq)
        assert top.n_angles == 1

        top.add_connection(angle_not_eq)
        assert top.n_angles == 2

    def test_equivalent_members_set(self):
        site1 = Site(name="SiteA")
        site2 = Site(name="SiteB")
        site3 = Site(name="SiteC")

        angle = Angle([site1, site2, site3])
        angle_eq = Angle([site3, site2, site1])
        angle_not_eq = Angle([site1, site3, site2])

        assert (tuple(angle_eq.connection_members)
                in angle.equivalent_members())
        assert (tuple(angle.connection_members)
                in angle_eq.equivalent_members())
        assert not (tuple(angle.connection_members)
                in angle_not_eq.equivalent_members())
