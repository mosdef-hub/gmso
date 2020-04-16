import pytest

from gmso.core.topology import Topology
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.atom_type import AtomType
from gmso.core.site import Site
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


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
                        connection_type=dihedral_type,
                        name='dihedral_name')

        assert site1.n_connections == 1
        assert site2.n_connections == 1
        assert site3.n_connections == 1
        assert site4.n_connections == 1
        assert len(connect.connection_members) == 4
        assert connect.connection_type is not None
        assert connect.name == 'dihedral_name'

    def test_dihedral_fake(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(GMSOError):
            Dihedral(connection_members=['fakesite1', 'fakesite2', 4.2])

    def test_dihedral_fake_dihedraltype(self):
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        site3 = Site(name='site3')
        site4 = Site(name='site4')
        with pytest.raises(GMSOError):
            Dihedral(connection_members=[site1, site2, site3, site4],
                  connection_type='Fake dihedraltype')

    def test_dihedral_constituent_types(self):
        site1 = Site(name='site1', position=[0,0,0], atom_type=AtomType(name='A'))
        site2 = Site(name='site2', position=[1,0,0], atom_type=AtomType(name='B'))
        site3 = Site(name='site3', position=[1,1,0], atom_type=AtomType(name='C'))
        site4 = Site(name='site4', position=[1,1,4], atom_type=AtomType(name='D'))
        dihtype = DihedralType(member_types=[site1.atom_type.name,
                                             site2.atom_type.name,
                                             site3.atom_type.name,
                                             site4.atom_type.name])
        dih = Dihedral(connection_members=[site1, site2, site3, site4],
                connection_type=dihtype)
        assert 'A' in dih.connection_type.member_types
        assert 'B' in dih.connection_type.member_types
        assert 'C' in dih.connection_type.member_types
        assert 'D' in dih.connection_type.member_types

    def test_dihedral_eq(self):
        site1 = Site(name='site1', position=[0, 0, 0])
        site2 = Site(name='site2', position=[1, 0, 0])
        site3 = Site(name='site3', position=[1, 1, 0])
        site4 = Site(name='site4', position=[1, 1, 1])

        ref_dihedral = Dihedral(
            connection_members=[site1, site2, site3, site4],
        )

        same_dihedral = Dihedral(
            connection_members=[site1, site2, site3, site4],
        )

        diff_dihedral = Dihedral(
            connection_members=[site1, site2, site3, site4],
        )

        assert ref_dihedral != same_dihedral
        assert ref_dihedral != diff_dihedral

    def test_add_equivalent_connections(self):
        site1 = Site(name="SiteA")
        site2 = Site(name="SiteB")
        site3 = Site(name="SiteC")
        site4 = Site(name="SiteD")

        dihedral = Dihedral([site1, site2, site3, site4])
        dihedral_eq = Dihedral([site4, site3, site2, site1])
        dihedral_not_eq = Dihedral([site4, site2, site3, site1])

        top = Topology()
        top.add_connection(dihedral)
        top.add_connection(dihedral_eq)
        assert top.n_dihedrals == 1

        top.add_connection(dihedral_not_eq)
        assert top.n_dihedrals == 2

    def test_equivalent_members_set(self):
        site1 = Site(name="SiteA")
        site2 = Site(name="SiteB")
        site3 = Site(name="SiteC")
        site4 = Site(name="SiteD")

        dihedral = Dihedral([site1, site2, site3, site4])
        dihedral_eq = Dihedral([site4, site3, site2, site1])
        dihedral_not_eq = Dihedral([site4, site2, site3, site1])


        assert (tuple(dihedral_eq.connection_members)
                in dihedral.equivalent_members())
        assert (tuple(dihedral.connection_members)
                in dihedral_eq.equivalent_members())
        assert not (tuple(dihedral.connection_members)
                in dihedral_not_eq.equivalent_members())

