import forcefield_utilities as ffutils
import parmed as pmd
import pytest
import unyt as u

import gmso
from gmso.exceptions import EngineIncompatibilityError
from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_mbuild, import_

if has_mbuild:
    mb = import_("mbuild")


class TestTop(BaseTest):
    def test_write_top(self, typed_ar_system):
        top = typed_ar_system
        top.save("ar.top")

    def test_top(self, spce_water):
        spce_water.save("spce_restraint.top", overwrite=True)
        with open("spce_restraint.top", "r") as f:
            contents = f.readlines()
        with open(get_path("spce_restraint.top"), "r") as f:
            ref_contents = f.readlines()

        assert len(contents) == len(ref_contents)
        for line, ref_line in zip(contents[1:], ref_contents[1:]):
            assert line == ref_line

    @pytest.mark.parametrize(
        "top",
        [
            "typed_ar_system",
            "typed_water_system",
            "typed_ethane",
            "typed_benzene_aa_system",
        ],
    )
    def test_pmd_loop(self, top, request):
        fname = f"{top}.top"
        top = request.getfixturevalue(top)
        top.save(fname, overwrite=True)
        pmd.load_file(fname)

    @pytest.mark.parametrize(
        "top",
        [
            "typed_ar_system",
            "typed_water_system",
            "typed_ethane",
            "typed_benzene_aa_system",
        ],
    )
    def test_against_ref(self, top, request):
        fname = top
        top = request.getfixturevalue(top)
        top.save(f"{fname}.top", overwrite=True)
        with open(f"{fname}.top") as f:
            conts = f.readlines()
        with open(get_path(f"{fname}_ref.top")) as f:
            ref_conts = f.readlines()

        assert len(conts) == len(ref_conts)

        for cont, ref_cont in zip(conts[1:], ref_conts[1:]):
            assert cont == ref_cont

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
        top = apply(top, ff)

        for site in top.sites:
            site.atom_type = ff.atom_types[site.atom_type.name]

        for bond in top.bonds:
            bond.bond_type = bond.connection_type = ff.bond_types[
                "opls_111~opls_112"
            ]

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

        for i in range(typed_ethane.n_impropers - 1, -1, -1):
            if not typed_ethane.impropers[i].improper_type:
                typed_ethane._impropers.pop(i)

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

    def test_benzene_restraints(self, typed_benzene_ua_system):
        top = typed_benzene_ua_system

        for bond in top.bonds:
            bond.restraint = {
                "r_eq": bond.bond_type.parameters["r_eq"],
                "k": 1000 * u.kJ / (u.mol * u.nm**2),
            }
        for angle in top.angles:
            # Apply restraint for angle
            angle.restraint = {
                "theta_eq": angle.angle_type.parameters["theta_eq"],
                "k": 1000 * u.kJ / u.mol,
                "n": 1,
            }

        for dihedral in top.dihedrals:
            # Apply restraint fot dihedral
            dihedral.restraint = {
                "phi_eq": 0 * u.degree,
                "delta_phi": 0 * u.degree,
                "k": 1000 * u.kJ / (u.mol * u.rad**2),
            }
        top.save("restrained_benzene_ua.top")

        with open("restrained_benzene_ua.top") as f:
            f_cont = f.readlines()

        with open(get_path("restrained_benzene_ua.top")) as ref:
            ref_cont = ref.readlines()

        assert len(f_cont) == len(ref_cont)

        ref_sections = dict()
        sections = dict()
        current_section = None
        for line, ref in zip(f_cont[1:], ref_cont[1:]):
            if line.startswith("["):
                assert line == ref
                current_section = line
                sections[current_section] = set()
                ref_sections[current_section] = set()
            elif line.startswith("#"):
                assert line == ref
            elif current_section is not None:
                sections[current_section].add(line)
                ref_sections[current_section].add(ref)

        for section, ref_section in zip(sections, ref_sections):
            assert section == ref_section
            if "dihedral" in section:
                # Need to deal with these separatelt due to member's order issue
                # Each dict will have the keys be members and values be their parameters
                members = dict()
                ref_members = dict()
                for line, ref in zip(
                    sections[section], ref_sections[ref_section]
                ):
                    line = line.split()
                    ref = ref.split()
                    members["-".join(line[:4])] = line[4:]
                    members["-".join(reversed(line[:4]))] = line[4:]
                    ref_members["-".join(ref[:4])] = ref[4:]
                    ref_members["-".join(reversed(ref[:4]))] = ref[4:]

                assert members == ref_members
                for member in members:
                    assert members[member] == ref_members[member]

            else:
                assert sections[section] == ref_sections[ref_section]
