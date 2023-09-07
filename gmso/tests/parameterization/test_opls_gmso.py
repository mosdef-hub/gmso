import glob
from pathlib import Path

import parmed as pmd
import pytest
from pkg_resources import resource_filename

from gmso.external.convert_parmed import from_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)


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


class TestOPLSGMSO(ParameterizationBaseTest):
    @pytest.mark.parametrize(
        "system_dir", get_foyer_opls_test_dirs(), ids=lambda p: p.name
    )
    def test_foyer_oplsaa_files(
        self,
        system_dir,
        oplsaa_gmso,
        oplsaa_foyer,
        assert_same_connection_params,
        assert_same_atom_params,
    ):
        top_file = str(system_dir / f"{system_dir.name}.top")
        gro_file = str(system_dir / f"{system_dir.name}.gro")
        struct = oplsaa_foyer.apply(pmd.load_file(top_file, xyz=gro_file))

        gmso_top_from_pmd = from_parmed(struct, refer_type=True)
        gmso_top = from_parmed(struct, refer_type=False)
        apply(gmso_top, oplsaa_gmso, speedup_by_molgraph=False)

        assert_same_atom_params(gmso_top, gmso_top_from_pmd)
        assert_same_connection_params(gmso_top, gmso_top_from_pmd)
        assert_same_connection_params(gmso_top, gmso_top_from_pmd, "angles")
        assert_same_connection_params(gmso_top, gmso_top_from_pmd, "dihedrals")
