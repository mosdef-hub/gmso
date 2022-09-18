import forcefield_utilities as ffutils
import parmed as pmd
import pytest
import unyt as u

import gmso
from gmso.exceptions import EngineIncompatibilityError
from gmso.formats.top import write_top
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_mbuild


class TestTop(BaseTest):
    def test_write_top(self, typed_ar_system):
        top = typed_ar_system
        top.save("ar.top")

    @pytest.mark.parametrize(
        "top", ["typed_ar_system", "typed_water_system", "typed_ethane"]
    )
    def test_pmd_loop(self, top, request):
        top = request.getfixturevalue(top)
        top.save("system.top")
        pmd.load_file("system.top")

    def test_modified_potentials(self, ar_system):
        top = ar_system

        ff = gmso.ForceField(get_fn("ar.xml"))

        for site in top.sites:
            site.atom_type = ff.atom_types["Ar"]

        top.update_topology()

        list(top.atom_types)[0].set_expression("sigma + epsilon*r")

        with pytest.raises(EngineIncompatibilityError):
            top.save("out.top")

        alternate_lj = "4*epsilon*sigma**12/r**12 - 4*epsilon*sigma**6/r**6"
        list(top.atom_types)[0].set_expression(alternate_lj)

        top.save("ar.top")

    def test_water_top(self, water_system):
        top = water_system

        ff = gmso.ForceField(get_path("tip3p.xml"))

        for site in top.sites:
            site.atom_type = ff.atom_types[site.name]

        top.update_atom_types()

        for bond in top.bonds:
            bond.bond_type = bond.connection_type = ff.bond_types[
                "opls_111~opls_112"
            ]

        top.update_bond_types()

        for molecule in top.unique_site_labels("molecule"):
            angle = gmso.core.angle.Angle(
                connection_members=[
                    site for site in top.iter_sites("molecule", molecule)
                ],
                name="opls_112~opls_111~opls_112",
                angle_type=ff.angle_types["opls_112~opls_111~opls_112"],
            )
            top.add_connection(angle)

        top.update_angle_types()

        top.save("water.top")

    def test_ethane_periodic(self, typed_ethane):
        from gmso.core.dihedral_type import DihedralType
        from gmso.lib.potential_templates import PotentialTemplateLibrary

        per_torsion = PotentialTemplateLibrary()["PeriodicTorsionPotential"]
        params = {
            "k": 10 * u.Unit("kJ / mol"),
            "phi_eq": 15 * u.Unit("degree"),
            "n": 3 * u.Unit("dimensionless"),
        }
        periodic_dihedral_type = DihedralType.from_template(
            potential_template=per_torsion, parameters=params
        )
        for dihedral in typed_ethane.dihedrals:
            dihedral.connection_type = periodic_dihedral_type

        typed_ethane.update_connection_types()

        typed_ethane.save("system.top")
        struct = pmd.load_file("system.top")
        assert len(struct.dihedrals) == 9

    def test_custom_defaults(self, typed_ethane):
        typed_ethane.save(
            "system.top",
            top_vars={"gen-pairs": "yes", "fudgeLJ": 0.5, "fudgeQQ": 0.5},
        )
        struct = pmd.load_file("system.top")
        assert struct.defaults.gen_pairs == "yes"
        assert struct.defaults.fudgeLJ == 0.5
        assert struct.defaults.fudgeQQ == 0.5

    @pytest.mark.skipif(not has_mbuild, reason="mBuild not installed.")
    def test_benzene_top(self):
        import mbuild as mb
        from mbuild.packing import fill_box

        from gmso.external import from_mbuild

        benzene = mb.load(get_fn("benzene.mol2"))
        benzene.children[0].name = "Benzene"
        box_of_benzene = fill_box(compound=benzene, n_compounds=5, density=1)
        top = from_mbuild(box_of_benzene)
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top=top, forcefields=oplsaa)
        top.save("benzene.top")

        with open("benzene.top") as f:
            f_cont = f.readlines()

        with open(get_path("benzene.top")) as ref:
            ref_cont = ref.readlines()

        assert len(f_cont) == len(ref_cont)
