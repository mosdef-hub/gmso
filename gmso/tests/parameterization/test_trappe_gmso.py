import glob
from pathlib import Path

import foyer
import parmed as pmd
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


class TestTrappeGMSO(BaseTest):
    @pytest.fixture(scope="session")
    def xml_loader(self):
        return FoyerFFs()

    @pytest.fixture(scope="session")
    def trappe_ua_gmso(self, xml_loader):
        return xml_loader.load("trappe_ua", rel_to_module=True).to_gmso_ff()

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
        top_gmso = gmso.Topology.load(mol2_file)
        struct_pmd = trappe_ua_openmm_foyer.apply(to_parmed(top_gmso))
        apply(top_gmso, trappe_ua_gmso, identify_connected_components=False)
        gmso_top_from_parmeterized_pmd = from_parmed(struct_pmd)
        for atom, mirror in zip(
            top_gmso.sites, gmso_top_from_parmeterized_pmd.sites
        ):
            assert atom.name == mirror.name
            assert u.allclose_units(atom.mass, mirror.mass, 1e-4)

            assert u.allclose_units(atom.atom_type.charge, mirror.charge, 1e-4)

            atom_params = atom.atom_type.get_parameters()
            mirror_params = mirror.atom_type.get_parameters()
            for k in atom_params:
                assert u.allclose_units(atom_params[k], mirror_params[k])
