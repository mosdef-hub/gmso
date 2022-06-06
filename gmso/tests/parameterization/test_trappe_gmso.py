import glob
from pathlib import Path

import foyer
import pytest
import unyt as u
from forcefield_utilities.xml_loader import FoyerFFs
from pkg_resources import resource_filename

import gmso
from gmso.external.convert_parmed import from_parmed, to_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.base_test import BaseTest


def get_foyer_trappe_test_dirs():
    all_dirs = glob.glob(resource_filename("foyer", "trappe_validation") + "/*")
    with open(
        resource_filename("foyer", "tests/implemented_trappe_tests.txt")
    ) as impl_file:
        correctly_implemented = set(impl_file.read().strip().split("\n"))

    parent_dirs = map(Path, all_dirs)
    parent_dirs = list(
        filter(
            lambda p: p.name in correctly_implemented
            and (p / f"{p.name}.mol2").exists(),
            parent_dirs,
        )
    )
    return parent_dirs


def _match_connection_parameters(
    from_gmso_top, from_parmed_top, connection_type="bonds"
):
    connection_types_original = {}
    connection_types_mirror = {}
    for connection in getattr(from_parmed_top, connection_type):
        connection_types_mirror[
            tuple(
                from_parmed_top.get_index(member)
                for member in connection.connection_members
            )
        ] = connection

    for connection in getattr(from_gmso_top, connection_type):
        connection_types_original[
            tuple(
                from_gmso_top.get_index(member)
                for member in connection.connection_members
            )
        ] = connection

    for key in connection_types_original:
        conn = connection_types_original[key]
        print(key)
        conn_mirror = connection_types_mirror[key]
        for param in conn_mirror.bond_type.parameters:
            assert u.allclose_units(
                conn_mirror.bond_type.parameters[param],
                conn.bond_type.parameters[param],
            )


class TestTrappeGMSO(BaseTest):
    @pytest.fixture(scope="session")
    def trappe_ua_openmm_foyer(self):
        return foyer.forcefields.load_TRAPPE_UA()

    @pytest.mark.parametrize(
        "system_dir", get_foyer_trappe_test_dirs(), ids=lambda p: p.name
    )
    def test_foyer_trappe_files(
        self,
        system_dir,
        trappe_ua_openmm_foyer,
        trappe_ua_gmso,
        are_equivalent_connections,
    ):
        mol2_file = system_dir / f"{system_dir.name}.mol2"
        gmso_top = gmso.Topology.load(mol2_file)
        struct_pmd = trappe_ua_openmm_foyer.apply(to_parmed(gmso_top))
        apply(gmso_top, trappe_ua_gmso, identify_connected_components=False)
        gmso_top_from_parmeterized_pmd = from_parmed(struct_pmd)
        for atom, mirror in zip(
            gmso_top.sites, gmso_top_from_parmeterized_pmd.sites
        ):
            assert atom.name == mirror.name
            assert u.allclose_units(atom.mass, mirror.mass, 1e-4)

            assert u.allclose_units(atom.atom_type.charge, mirror.charge, 1e-4)

            atom_params = atom.atom_type.get_parameters()
            mirror_params = mirror.atom_type.get_parameters()
            for k in atom_params:
                assert u.allclose_units(atom_params[k], mirror_params[k])

        _match_connection_parameters(gmso_top, gmso_top_from_parmeterized_pmd)
        _match_connection_parameters(
            gmso_top, gmso_top_from_parmeterized_pmd, "angles"
        )
        _match_connection_parameters(
            gmso_top, gmso_top_from_parmeterized_pmd, "dihedrals"
        )
