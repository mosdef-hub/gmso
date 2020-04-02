import pytest

from gmso.core.topology import Topology
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.atom_type import AtomType
from gmso.core.site import Site
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestImproper(BaseTest):
    def test_improper_nonparametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0
        assert site4.n_connections == 0

        connect = Improper(connection_members=[site1, site2, site3, site4])

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert site4.n_connections == 1
        assert connect.connection_type is None

    def test_improper_parametrized(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')

        assert site1.n_connections == 0
        assert site2.n_connections == 0
        assert site3.n_connections == 0
        assert site4.n_connections == 0
        improper_type = ImproperType()

        connect = Improper(connection_members=[site1, site2, site3, site4],
                        connection_type=improper_type,
                        name='improper_name')

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert site4.n_connections == 1
        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.name == 'improper_name'

    def test_improper_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(GMSOError):
            Improper(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_improper_fake_impropertype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(GMSOError):
            Improper(connection_members=[site1, site2, site3, site4],
                  connection_type='Fake impropertype')

    def test_improper_constituent_types(self):
        site1 = Site(name='site1', position=[0,0,0], atom_type=AtomType(name='A'))
        site2 = Site(name='site2', position=[1,0,0], atom_type=AtomType(name='B'))
        site3 = Site(name='site3', position=[1,1,0], atom_type=AtomType(name='C'))
        site4 = Site(name='site4', position=[1,1,4], atom_type=AtomType(name='D'))
        imptype = ImproperType(member_types=[site1.atom_type.name,
                                             site2.atom_type.name,
                                             site3.atom_type.name,
                                             site4.atom_type.name])
        imp = Improper(connection_members=[site1, site2, site3, site4],
                connection_type=imptype)
        assert 'A' in imp.connection_type.member_types
        assert 'B' in imp.connection_type.member_types
        assert 'C' in imp.connection_type.member_types
        assert 'D' in imp.connection_type.member_types

    def test_improper_eq(self):
        site1 = Site(name='site1', position=[0, 0, 0])
        site2 = Site(name='site2', position=[1, 0, 0])
        site3 = Site(name='site3', position=[1, 1, 0])
        site4 = Site(name='site4', position=[1, 1, 1])

        ref_improper = Improper(
            connection_members=[site1, site2, site3, site4],
        )

        same_improper = Improper(
            connection_members=[site1, site2, site3, site4],
        )

        diff_improper = Improper(
            connection_members=[site1, site2, site3, site4],
        )

        assert ref_improper != same_improper
        assert ref_improper != diff_improper

    def test_add_equivalent_connections(self):
        site1 = Site(name="SiteA")
        site2 = Site(name="SiteB")
        site3 = Site(name="SiteC")
        site4 = Site(name="SiteD")

        improper = Improper([site1, site2, site3, site4])
        improper_eq = Improper([site1, site3, site2, site4])
        improper_not_eq = Improper([site2, site3, site1, site4])

        top = Topology()
        top.add_connection(improper)
        top.add_connection(improper_eq)
        assert top.n_impropers == 1
        top.add_connection(improper_not_eq)
        assert top.n_impropers == 2
