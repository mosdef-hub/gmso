import glob
from pathlib import Path

import parmed as pmd
import pytest
from forcefield_utilities.xml_loader import FoyerFFs
from pkg_resources import resource_filename

from gmso.external.convert_parmed import from_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.base_test import BaseTest


def get_foyer_opls_test_dirs():
    all_dirs = glob.glob(resource_filename("foyer", "opls_validation") + "/*")
    with open(
        resource_filename("foyer", "tests/implemented_opls_tests.txt")
    ) as impl_file:
        correctly_implemented = set(impl_file.read().strip().split("\n"))

    parent_dirs = map(Path, all_dirs)
    parent_dirs = list(
        filter(
            lambda p: p.name in correctly_implemented
            and (p / f"{p.name}.top").exists(),
            parent_dirs,
        )
    )
    return parent_dirs


class TestOPLSGMSO(BaseTest):
    @pytest.fixture(scope="session")
    def xml_loader(self):
        return FoyerFFs()

    @pytest.fixture(scope="session")
    def oplsaa(self, xml_loader):
        return xml_loader.load("oplsaa", rel_to_module=True).to_gmso_ff()

    @pytest.mark.parametrize(
        "system_dir", get_foyer_opls_test_dirs(), ids=lambda p: p.name
    )
    def test_foyer_oplsaa_files(
        self, system_dir, oplsaa, are_equivalent_topologies
    ):
        top_file = str(system_dir / f"{system_dir.name}.top")
        gro_file = str(system_dir / f"{system_dir.name}.gro")
        struct = pmd.load_file(top_file, xyz=gro_file)
        gmso_top_unparam = from_parmed(struct, refer_type=False)
        apply(gmso_top_unparam, oplsaa, identify_connected_components=False)
        assert are_equivalent_topologies(
            gmso_top_unparam, from_parmed(struct, refer_type=True)
        ), f"{system_dir} Failed"
