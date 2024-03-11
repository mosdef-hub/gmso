import foyer
import mbuild as mb
import pytest
import unyt as u
from forcefield_utilities.xml_loader import FoyerFFs
from mbuild.lib.molecules import Ethane, Methane

from gmso.external.convert_mbuild import from_mbuild
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class ParameterizationBaseTest(BaseTest):
    @pytest.fixture(scope="session")
    def xml_loader(self):
        return FoyerFFs()

    @pytest.fixture(scope="session")
    def oplsaa_gmso(self, xml_loader):
        return xml_loader.load("oplsaa").to_gmso_ff()

    @pytest.fixture(scope="session")
    def trappe_ua_gmso(self, xml_loader):
        return xml_loader.load("trappe-ua").to_gmso_ff()

    @pytest.fixture(scope="session")
    def fake_improper_ff_gmso(self, xml_loader):
        return xml_loader.load(
            get_path("fake_ethane_impropers.xml")
        ).to_gmso_ff()

    @pytest.fixture(scope="session")
    def benzene_alkane_aa_ff_gmso(self, xml_loader):
        return xml_loader.load(
            get_path("benzene_and_alkane_branched_benzene_aa.xml")
        ).to_gmso_ff()

    @pytest.fixture(scope="session")
    def oplsaa_foyer(self):
        return foyer.forcefields.load_OPLSAA()

    @pytest.fixture(scope="session")
    def trappe_ua_foyer(self):
        return foyer.forcefields.load_TRAPPE_UA()

    @pytest.fixture(scope="session")
    def assert_same_connection_params(self):
        def _assert_same_connection_params(top1, top2, connection_type="bonds"):
            """Match connection parameters between two gmso topologies."""
            connection_types_top1 = {}
            for connection in getattr(top1, connection_type):
                eq_connsList = connection.equivalent_members()
                indexList = [
                    tuple(map(lambda x: top1.get_index(x), conn))
                    for conn in eq_connsList
                ]
                atom_indicesList = sorted(indexList)[0]
                connection_types_top1[atom_indicesList] = connection
            connection_types_top2 = {}
            for connection in getattr(top2, connection_type):
                eq_connsList = connection.equivalent_members()
                indexList = [
                    tuple(map(lambda x: top2.get_index(x), conn))
                    for conn in eq_connsList
                ]
                atom_indicesList = sorted(indexList)[0]
                connection_types_top2[atom_indicesList] = connection

            # for connection in getattr(top2, connection_type):
            #    connection_types_mirror[
            #        tuple(
            #            top2.get_index(member)
            #            for member in sort_connection_members(connection, "atom_type")
            #        )
            #    ] = connection

            # for connection in getattr(top1, connection_type):
            #    connection_types_original[
            #        tuple(
            #            top1.get_index(member)
            #            for member in sort_connection_members(connection, "atom_type")
            #        )
            #    ] = connection

            for key in connection_types_top1:
                conn1 = connection_types_top1[key]
                conn2 = connection_types_top2[key]
                conn_type_attr = connection_type[:-1] + "_type"
                conn_type1 = getattr(conn1, conn_type_attr)
                conn_type2 = getattr(conn2, conn_type_attr)
                for param in conn_type1.parameters:
                    assert u.allclose_units(
                        conn_type2.parameters[param],
                        conn_type1.parameters[param],
                    )

        return _assert_same_connection_params

    @pytest.fixture(scope="session")
    def assert_same_atom_params(self):
        def _assert_same_atom_params(top1, top2):
            """Match atom parameters between two gmso topologies.

            Notes
            -----
            This is specific
            """
            for atom, mirror in zip(top1.sites, top2.sites):
                assert atom.name == mirror.name
                assert u.allclose_units(atom.mass, mirror.mass, 1e-3)

                atom_params = atom.atom_type.get_parameters()
                mirror_params = mirror.atom_type.get_parameters()

                for k in atom_params:
                    assert u.allclose_units(atom_params[k], mirror_params[k])

        return _assert_same_atom_params

    @pytest.fixture
    def ethane_methane_top(self):
        cmpd = mb.Compound()
        cmpd.add(Ethane())
        cmpd.add(Methane())
        gmso_top = from_mbuild(cmpd)
        gmso_top.identify_connections()
        return gmso_top

    @pytest.fixture
    def ethane_box_with_methane(self):
        cmpd_box = mb.fill_box([Ethane(), Methane()], [50, 50], density=1.0)
        return from_mbuild(cmpd_box)
