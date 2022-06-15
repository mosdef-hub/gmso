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
        return xml_loader.load("oplsaa", rel_to_module=True).to_gmso_ff()

    @pytest.fixture(scope="session")
    def trappe_ua_gmso(self, xml_loader):
        return xml_loader.load("trappe_ua", rel_to_module=True).to_gmso_ff()

    @pytest.fixture(scope="session")
    def fake_improper_ff_gmso(self, xml_loader):
        return xml_loader.load(
            get_path("fake_ethane_impropers.xml"), rel_to_module=True
        ).to_gmso_ff()

    @pytest.fixture(scope="session")
    def benzene_alkane_aa_ff_gmso(self, xml_loader):
        return xml_loader.load(
            get_path("benzene_and_alkane_branched_benzene_aa.xml"),
            rel_to_module=True,
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
            connection_types_original = {}
            connection_types_mirror = {}
            for connection in getattr(top2, connection_type):
                connection_types_mirror[
                    tuple(
                        top2.get_index(member)
                        for member in connection.connection_members
                    )
                ] = connection

            for connection in getattr(top1, connection_type):
                connection_types_original[
                    tuple(
                        top1.get_index(member)
                        for member in connection.connection_members
                    )
                ] = connection

            for key in connection_types_original:
                conn = connection_types_original[key]
                conn_mirror = connection_types_mirror[key]
                conn_type_attr = connection_type[:-1] + "_type"
                conn_type_mirror = getattr(conn_mirror, conn_type_attr)
                conn_type = getattr(conn, conn_type_attr)
                for param in conn_type.parameters:
                    assert u.allclose_units(
                        conn_type_mirror.parameters[param],
                        conn_type.parameters[param],
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

    @pytest.fixture(scope="session")
    def ethane_methane_top(self):
        cmpd = mb.Compound()
        cmpd.add(Ethane())
        cmpd.add(Methane())
        gmso_top = from_mbuild(cmpd)
        gmso_top.identify_connections()
        return gmso_top

    @pytest.fixture(scope="session")
    def ethane_box_with_methane(self):
        cmpd_box = mb.fill_box([Ethane(), Methane()], [50, 50], density=1.0)
        return from_mbuild(cmpd_box)
