import pytest

from topology.core.dihedral import Dihedral
from topology.core.dihedral_type import DihedralType
from topology.core.atom_type import AtomType
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