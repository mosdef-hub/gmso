from importlib.resources import as_file, files
from pathlib import Path

import parmed as pmd
import pytest

from gmso import ForceField
from gmso.external.convert_parmed import from_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)
from gmso.tests.utils import get_path


def get_foyer_opls_test_dirs():
    fn = files("foyer") / "opls_validation"
    all_dirs = fn.glob("*")
    tests_fn = files("foyer") / "tests/implemented_opls_tests.txt"
    with as_file(tests_fn) as tempPath:
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


class TestGeneralParameterizations(ParameterizationBaseTest):
    def test_wildcards(self, ethane_methane_top):
        from gmso.core.views import PotentialFilters

        ff = ForceField(get_path("alkanes_wildcards.xml"))
        ptop = apply(ethane_methane_top, ff, identify_connections=True)
        assert ptop.is_fully_typed()
        assert len(ptop.bond_types) == 11
        assert len(ptop.bond_types(PotentialFilters.UNIQUE_NAME_CLASS)) == 1
        assert ptop.bonds[0].bond_type.member_types == ("*", "*")

        assert len(ptop.angle_types) == 12 + 6  # ethane + methane
        assert len(ptop.angle_types(PotentialFilters.UNIQUE_NAME_CLASS)) == 1
        assert ptop.angles[0].angle_type.member_types == ("*", "*", "*")

        assert len(ptop.dihedral_types) == 9
        assert len(ptop.dihedral_types(PotentialFilters.UNIQUE_NAME_CLASS)) == 1
        assert ptop.dihedrals[0].dihedral_type.member_types == ("*", "*", "*", "*")
