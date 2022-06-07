import glob
from pathlib import Path

import parmed as pmd
import pytest
import unyt as u
from pkg_resources import resource_filename

from gmso.external.convert_parmed import from_parmed
from gmso.parameterization.parameterize import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import match_connection_parameters


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
    return parent_dirs[0:5]


class TestOPLSGMSO(BaseTest):
    @pytest.mark.parametrize(
        "system_dir", get_foyer_opls_test_dirs(), ids=lambda p: p.name
    )
    def test_foyer_oplsaa_files(
        self, system_dir, oplsaa_gmso, oplsaa_foyer, are_equivalent_topologies
    ):
        top_file = str(system_dir / f"{system_dir.name}.top")
        gro_file = str(system_dir / f"{system_dir.name}.gro")
        struct = pmd.load_file(top_file, xyz=gro_file)
        gmso_top = from_parmed(struct, refer_type=False)
        apply(gmso_top, oplsaa_gmso, identify_connected_components=False)

        pmd_struct_param = oplsaa_foyer.apply(struct)
        gmso_top_from_pmd = from_parmed(pmd_struct_param, refer_type=True)

        for atom, mirror in zip(gmso_top.sites, gmso_top_from_pmd.sites):
            assert atom.name == mirror.name
            assert u.allclose_units(atom.mass, mirror.mass, 1e-3)

            assert u.allclose_units(atom.atom_type.charge, mirror.charge, 1e-3)

            atom_params = atom.atom_type.get_parameters()
            mirror_params = mirror.atom_type.get_parameters()

            for k in atom_params:
                assert u.allclose_units(atom_params[k], mirror_params[k])

        match_connection_parameters(gmso_top, gmso_top_from_pmd)
        match_connection_parameters(gmso_top, gmso_top_from_pmd, "angles")
        match_connection_parameters(gmso_top, gmso_top_from_pmd, "dihedrals")
