from pathlib import Path

import importlib_resources
import parmed as pmd
import pytest

from gmso.external.convert_parmed import from_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)


def get_foyer_opls_test_dirs():
    fn = importlib_resources.files("foyer") / "opls_validation"
    all_dirs = fn.glob("*")
    tests_fn = importlib_resources.files("foyer") / "tests/implemented_opls_tests.txt"
    with importlib_resources.as_file(tests_fn) as tempPath:
        with open(tempPath) as impl_file:
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
        struct = oplsaa_foyer.apply(
            pmd.load_file(top_file, xyz=gro_file, parametrize=False)
        )

        gmso_top_from_pmd = from_parmed(struct, refer_type=True)
        gmso_top = from_parmed(struct, refer_type=False)
        apply(gmso_top, oplsaa_gmso, speedup_by_molgraph=False)

        assert_same_atom_params(gmso_top, gmso_top_from_pmd)
        assert_same_connection_params(gmso_top, gmso_top_from_pmd)
        assert_same_connection_params(gmso_top, gmso_top_from_pmd, "angles")
        assert_same_connection_params(gmso_top, gmso_top_from_pmd, "dihedrals")
