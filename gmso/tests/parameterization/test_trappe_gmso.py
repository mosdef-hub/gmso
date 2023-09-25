import glob
from pathlib import Path

import pytest
from pkg_resources import resource_filename

from gmso.core.topology import Topology
from gmso.external.convert_parmed import from_parmed, to_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)


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


class TestTrappeGMSO(ParameterizationBaseTest):
    @pytest.mark.parametrize(
        "system_dir", get_foyer_trappe_test_dirs(), ids=lambda p: p.name
    )
    def test_foyer_trappe_files(
        self,
        system_dir,
        trappe_ua_foyer,
        trappe_ua_gmso,
        assert_same_connection_params,
        assert_same_atom_params,
    ):
        mol2_file = system_dir / f"{system_dir.name}.mol2"
        gmso_top = Topology.load(mol2_file)
        struct_pmd = trappe_ua_foyer.apply(to_parmed(gmso_top))
        apply(
            gmso_top,
            trappe_ua_gmso,
            speedup_by_molgraph=False,
            identify_connections=True,
        )
        gmso_top_from_parmeterized_pmd = from_parmed(struct_pmd)

        assert_same_atom_params(gmso_top_from_parmeterized_pmd, gmso_top)
        assert_same_connection_params(gmso_top, gmso_top_from_parmeterized_pmd)
        assert_same_connection_params(
            gmso_top, gmso_top_from_parmeterized_pmd, "angles"
        )
        assert_same_connection_params(
            gmso_top, gmso_top_from_parmeterized_pmd, "dihedrals"
        )
