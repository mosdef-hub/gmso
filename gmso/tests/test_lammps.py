import copy
import os

import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

import gmso
from gmso import Topology
from gmso.core.box import Box
from gmso.core.views import PotentialFilters

pfilter = PotentialFilters.UNIQUE_SORTED_NAMES
from gmso.exceptions import EngineIncompatibilityError
from gmso.external import from_parmed, to_parmed
from gmso.formats.formats_registry import UnsupportedFileFormatError
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


def compare_lammps_files(fn1, fn2, skip_linesList=[], offsets=None):
    """Check for line by line equality between lammps files, by any values.

    offsets = [file1: [(start, step)], file2: [(start, step)]
    """
    with open(fn1, "r") as f:
        line1 = f.readlines()
    with open(fn2, "r") as f:
        line2 = f.readlines()
    length1 = len(line1)
    length2 = len(line2)
    line_counter1 = 0
    line_counter2 = 0
    while True:
        # check for offsets
        if offsets is not None:
            if len(offsets[0]) == 0:
                pass
            elif offsets[0][0][0] == line_counter1:
                line_counter1 += offsets[0][0][1]
                offsets[0].pop(0)
            if len(offsets[1]) == 0:
                pass
            elif offsets[1][0][0] == line_counter2:
                line_counter2 += offsets[1][0][1]
                offsets[1].pop(0)
        if line_counter1 in skip_linesList and line_counter1 == line_counter2:
            line_counter1 += 1
            line_counter2 += 1
            continue
        l1 = line1[line_counter1]
        l2 = line2[line_counter2]
        print(f"############### ({line_counter1}) ({line_counter2})")
        print(l1, l2)

        for arg1, arg2 in zip(l1.split(), l2.split()):
            try:
                comp1 = float(arg1)
                comp2 = float(arg2)
            except ValueError:
                comp1 = str(arg1)
                comp2 = str(arg2)
            if isinstance(comp1, float):
                assert np.isclose(
                    comp1, comp2, 1e-3
                ), f"The following two lines have not been found to have equality {l1} and {l2}"
        line_counter1 += 1
        line_counter2 += 1
        if line_counter1 >= length1 or line_counter2 >= length2:
            break
    return True


class TestLammpsWriter(BaseTest):
    @pytest.mark.parametrize(
        "fname", ["data.lammps", "data.data", "data.lammpsdata"]
    )
    def test_write_lammps(self, fname, typed_ar_system):
        typed_ar_system.save(fname)

    def test_ethane_lammps_conversion(
        self, typed_ethane, are_equivalent_topologies
    ):
        typed_ethane.save("ethane.lammps")
        read_top = Topology.load("ethane.lammps")
        assert are_equivalent_topologies(read_top, typed_ethane)

    def test_opls_lammps(self, typed_ethane_opls, are_equivalent_topologies):
        typed_ethane_opls.save("ethane.lammps")
        read_top = Topology.load("ethane.lammps")
        assert are_equivalent_topologies(read_top, typed_ethane_opls)

    def test_water_lammps(self, typed_water_system, are_equivalent_topologies):
        typed_water_system.save("water.lammps")
        read_top = Topology.load("water.lammps")
        assert are_equivalent_topologies(read_top, typed_water_system)

    def test_read_lammps(self, filename=get_path("data.lammps")):
        gmso.Topology.load(filename)

    def test_read_box(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, typed_ar_system):
        typed_ar_system.save("ar.lammps")
        read = gmso.Topology.load("ar.lammps")
        assert read.n_sites == 100

    def test_read_mass(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        masses = [i.mass for i in read.atom_types]

        assert_allclose_units(
            masses, u.unyt_array(1.0079, u.g / u.mol), rtol=1e-5, atol=1e-8
        )

    def test_read_charge(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        charge = [i.charge for i in read.atom_types]

        assert_allclose_units(
            charge, u.unyt_array(0, u.C), rtol=1e-5, atol=1e-8
        )

    def test_read_sigma(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["sigma"], u.unyt_array(3, u.angstrom), rtol=1e-5, atol=1e-8
        )

    def test_read_epsilon(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["epsilon"],
            u.unyt_array(0.0717, (u.kcal / u.mol)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_water(self, typed_water_system):
        typed_water_system.save("water.lammps")
        water = gmso.Topology.load("water.lammps")

        assert_allclose_units(
            water.sites[0].charge,
            u.unyt_array(-0.834, u.elementary_charge),
            rtol=1e-5,
            atol=1e-8,
        )
        assert water.n_sites == 3
        assert water.n_connections == 3

    def test_read_lammps_triclinic(self, typed_ar_system):
        typed_ar_system.box = Box(lengths=[1, 1, 1], angles=[60, 90, 120])
        typed_ar_system.save("triclinic.lammps")

        read = gmso.Topology.load("triclinic.lammps")
        assert_allclose_units(
            read.box.lengths,
            u.unyt_array([1, 1, 1], u.nm),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            read.box.angles,
            u.unyt_array([60, 90, 120], u.degree),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_n_bonds(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_bonds == 7

    def test_read_n_angles(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_angles == 12

    def test_read_bond_params(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")
        bond_params = [i.parameters for i in read.bond_types(filter_by=pfilter)]

        assert_allclose_units(
            bond_params[0]["k"],
            u.unyt_array(680, (u.kcal / u.mol / u.angstrom / u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[0]["r_eq"],
            u.unyt_array(1.09, (u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[1]["k"],
            u.unyt_array(536, (u.kcal / u.mol / u.angstrom / u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[1]["r_eq"],
            u.unyt_array(1.529, (u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_angle_params(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")
        angle_params = [
            i.parameters for i in read.angle_types(filter_by=pfilter)
        ]

        assert_allclose_units(
            angle_params[0]["k"],
            u.unyt_array(75, (u.kcal / u.mol / u.radian / u.radian)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[0]["theta_eq"],
            u.unyt_array(110.7, (u.degree)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[1]["k"],
            u.unyt_array(66, (u.kcal / u.mol / u.radian / u.radian)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[1]["theta_eq"],
            u.unyt_array(107.8, (u.degree)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_n_diherals(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_dihedrals == 9
        assert len(read.dihedral_types(filter_by=pfilter)) == 1

    # TODO: would be good to create a library of molecules and styles to test
    # Test potential styles that are directly comparable to ParmEd writers.
    @pytest.mark.parametrize(
        "top",
        [
            "typed_ethane",
            "typed_methylnitroaniline",
            "typed_methaneUA",
        ],
    )
    def test_lammps_vs_parmed_by_mol(self, top, request):
        """Parmed LAMMPSDATA Compare outputs.

        atom_style = 'full', 'atomic', 'charge', 'molecular'
        unit_style = 'real', 'lj'
        improper_style = "cvff"
        dihedral_style = 'CHARMM', 'OPLS',
        angle_style = 'harmonic', 'urey_bradleys'
        bond_style = 'harmonic
        pair_style = 'lj
        """
        top = request.getfixturevalue(top)
        pmd_top = to_parmed(top)
        top.save("gmso.lammps")
        pmd_top.impropers = []
        from mbuild.formats.lammpsdata import (
            write_lammpsdata as mb_write_lammps,
        )

        mb_write_lammps(
            structure=pmd_top,
            filename="pmd.lammps",
            detect_forcefield_style=True,
            use_dihedrals=False,
            use_rb_torsions=True,
            mins=[0, 0, 0],
            maxs=top.box.lengths.convert_to_units(u.nm),
        )
        assert compare_lammps_files(
            "gmso.lammps",
            "pmd.lammps",
            skip_linesList=[0],
            offsets=[[[1, 1], ["none", "none"]], [["none", "none"]]],
        )

    @pytest.mark.parametrize(
        "atom_style", ["atomic", "charge", "molecular", "full"]
    )
    def test_lammps_vs_parmed_by_styles(
        self, atom_style, typed_ethane, parmed_ethane
    ):
        """Test all support styles in lammps writer.
        _______References_______
        See https://docs.lammps.org/atom_style.html for more info.
        """
        typed_ethane.save("gmso.lammps", atom_style=atom_style)
        from mbuild.formats.lammpsdata import (
            write_lammpsdata as mb_write_lammps,
        )

        mb_write_lammps(
            structure=parmed_ethane,
            filename="pmd.lammps",
            atom_style=atom_style,
            detect_forcefield_style=True,
            use_dihedrals=False,
            use_rb_torsions=True,
            mins=[0, 0, 0],
            maxs=typed_ethane.box.lengths.convert_to_units(u.nm),
        )
        assert compare_lammps_files(
            "gmso.lammps",
            "pmd.lammps",
            skip_linesList=[0],
            offsets=[[[0, 1]], [["none", "none"]]],
        )

    def test_lammps_default_conversions(
        self, typed_ethane, harmonic_parmed_types_charmm
    ):
        """Test for parameter intraconversions with potential styles.

        These include:
        bonds: factor of 2 harmonic k
        angles: factor of 2 harmonic k
        dihedrals: RB torsions to OPLS
        impropers: factor of 2 harmonic k
        pairs:
        additional: All styles to zero and none
        """
        typed_ethane.save("opls.lammps")
        with open("opls.lammps", "r") as f:
            lines = f.readlines()
        assert lines[39:42] == [
            "Dihedral Coeffs #OPLSTorsionPotential\n",
            "#\tk1 (kcal/mol)\tk2 (kcal/mol)\tk3 (kcal/mol)\tk4 (kcal/mol)\n",
            "1\t0.000000\t-0.000000\t0.300000\t-0.000000\t# opls_140\topls_135\topls_135\topls_140\n",
        ]

        struc = harmonic_parmed_types_charmm
        from mbuild.formats.lammpsdata import (
            write_lammpsdata as mb_write_lammps,
        )

        mb_write_lammps(struc, "pmd.lammps")
        top = from_parmed(struc)
        top.save("gmso.lammps")
        assert compare_lammps_files(
            "gmso.lammps",
            "pmd.lammps",
            skip_linesList=[0],
            offsets=[[[0, 1], [17, 1]], []],
        )
        out_lammps = open("gmso.lammps", "r").readlines()
        found_impropers = False
        for i, line in enumerate(out_lammps):
            if "Improper Coeffs" in line:
                assert "HarmonicTorsionPotential" in line
                assert "k" in out_lammps[i + 1]
                assert "phi_eq" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 3
                assert out_lammps[i + 2].split("#")[0].split()[0] == "1"
                found_impropers = True
        assert found_impropers

    def test_lammps_strict_true(self, typed_ethane):
        with pytest.raises(EngineIncompatibilityError):
            typed_ethane.save("error.lammps", strict_potentials=True)
        typed_ethane = typed_ethane.convert_potential_styles(
            {
                "dihedrals": "OPLSTorsionPotential",
                "bonds": "LAMMPSHarmonicBondPotential",
                "angles": "LAMMPSHarmonicAnglePotential",
            }
        )
        typed_ethane.save("test2.lammps", strict_potentials=True)

    # TODO: Need to add a list of forcefield templates to check with
    def test_lammps_potential_styles(self, typed_ethane):
        """Test for parameter handling of potential styles.

        ______Styles______
        The GMSO topology potentials must be written in one of these formats in order to be compatable in lammps
        bond_styles: ["none", "zero", "fene", "gromos", "morse", "harmonic"]
        angle_styles: ["none", "zero", "amoeba", "charmm", "cosine", "fourier", "harmonic"]
        dihedral_styles: ["none", "zero", "charmm", "fourier", "harmonic", "opls", "quadratic"]
        improper_styles: ["none", "zero", "amoeba", "fourier", "harmonic", "umbrella"]
        pair_styles: ["none", "zero", "amoeba", "buck", "coul", "dpd", "gauss", "harmonic", "lj", "mie", "morse", "yukawa"]
        Additional Styles: ['special_bonds']

        _______References_______
         See https://docs.lammps.org/Commands_category.html#force-fields for more info.
        """
        # TODO: Create a library of molecules that use the above styles
        typed_ethane.save("test.lammps")
        # TODO: Read and check test.lammps for correct writing

    @pytest.mark.parametrize(
        "unit_style",
        ["real", "metal", "si", "cgs", "electron", "micro", "nano", "lj"],
    )
    def test_lammps_units(self, typed_ethane, unit_style):
        """Generate topoogy with different units and check the output.
        Supporte styles are: ["real", "lj", "metal", "si", "cgs", "electron", "micro", "nano"]
        _______References_______
        https://docs.lammps.org/units.html
        """
        # check the initial set of units
        from gmso.formats.lammpsdata import get_units
        from gmso.utils.units import LAMMPS_UnitSystems

        base_unyts = LAMMPS_UnitSystems(unit_style)

        # real units should be in: [g/mol, angstroms, fs, kcal/mol, kelvin, electon charge, ...]
        mass_multiplierDict = {
            "si": 1,
            "cgs": 1 / 1000,
            "micro": 1e-15,
            "nano": 1e-21,
        }
        length_multiplierDict = {"si": 1e10, "cgs": 1e8}
        if unit_style in ["si", "cgs", "micro", "nano"]:
            typed_ethane.box.lengths *= length_multiplierDict.get(unit_style, 1)
            for atype in typed_ethane.atom_types:
                atype.mass = 12 * mass_multiplierDict[unit_style] * u.kg
        typed_ethane.save("ethane.lammps", unit_style=unit_style)
        real_top = Topology().load("ethane.lammps", unit_style=unit_style)
        energy_unit = get_units(base_unyts, "energy")
        angle_unit = get_units(base_unyts, "angle_eq")
        length_unit = get_units(base_unyts, "length")
        charge_unit = get_units(base_unyts, "charge")
        assert (
            real_top.dihedrals[0].dihedral_type.parameters["k1"].units
            == energy_unit
        )
        assert (
            real_top.angles[0].angle_type.parameters["theta_eq"].units
            == angle_unit
        )
        assert (
            real_top.bonds[0].bond_type.parameters["r_eq"].units == length_unit
        )
        assert real_top.sites[0].charge.units == charge_unit
        if unit_style == "lj":
            largest_eps = max(
                list(
                    map(
                        lambda x: x.parameters["epsilon"],
                        typed_ethane.atom_types,
                    )
                )
            )
            largest_sig = max(
                list(
                    map(
                        lambda x: x.parameters["sigma"], typed_ethane.atom_types
                    )
                )
            )
            assert_allclose_units(
                real_top.dihedrals[0].dihedral_type.parameters["k1"],
                (
                    typed_ethane.dihedrals[0].dihedral_type.parameters["k1"]
                    / largest_eps
                ),
                rtol=1e-5,
                atol=1e-8,
            )
            assert (
                real_top.dihedrals[0].dihedral_type.parameters["k1"]
                == typed_ethane.dihedrals[0].dihedral_type.parameters["k1"]
                / largest_eps
            )
            assert_allclose_units(
                real_top.bonds[0].bond_type.parameters["r_eq"],
                (
                    typed_ethane.bonds[0].bond_type.parameters["r_eq"]
                    / largest_sig
                ),
                rtol=1e-5,
                atol=1e-8,
            )

    def test_lammps_errors(self, typed_ethane):
        with pytest.raises(UnsupportedFileFormatError):
            typed_ethane.save("e.lammmps")
        missing_bonds_top = copy.deepcopy(typed_ethane)
        for bond in missing_bonds_top.bonds:
            bond.bond_type = None
        with pytest.raises(AttributeError):
            missing_bonds_top.save("e.lammps")
        with pytest.raises(ValueError):
            typed_ethane.save(
                "e.lammps", unit_style="lj", lj_cfactorsDict={"bonds": 1}
            )
        with pytest.raises(ValueError):
            typed_ethane.save("e.lammps", lj_cfactorsDict={"energy": "kJ/mol"})

        with pytest.raises(ValueError):
            typed_ethane.save("error.lammps", atom_style="None")

        with pytest.raises(ValueError):
            typed_ethane.save("error.lammps", unit_style="None")

    def test_lammps_validate_units(self, typed_methylnitroaniline):
        from gmso.formats.lammpsdata import _validate_unit_compatibility
        from gmso.utils.units import LAMMPS_UnitSystems

        base_unyts = LAMMPS_UnitSystems("si")
        with pytest.raises(AssertionError):
            _validate_unit_compatibility(typed_methylnitroaniline, base_unyts)

        usys = LAMMPS_UnitSystems("real")
        _validate_unit_compatibility(typed_methylnitroaniline, usys)

    def test_units_in_headers(self, typed_ethane):
        """Make sure units are written out properly."""
        typed_ethane.save("ethane.lammps")
        with open("ethane.lammps", "r") as f:
            lines = f.readlines()

        unitsDict = {
            "Pair": {"epsilon": "kcal/mol", "sigma": "Å"},
            "Bond": {"k": "kcal/(mol*Å**2)", "r_eq": "Å"},
            "Angle": {"k": "kcal/(mol*rad**2)", "theta_eq": "degrees"},
            "Dihedral": {"k1": "kcal/mol"},
            "Improper": {"k": "kcal/(mol*rad**2)", "phi_eq": "degrees"},
        }
        for i, line in enumerate(lines):
            if "Coeffs" in line:
                units = lines[i + 1].split(" \n")
                for j in range(len(units[1:-1:2])):
                    assert units[j * 2 + 1] == unitsDict[units[j * 2 + 2]]

    def test_atom_style_printing(self, typed_ethane):
        """Check writers for correctly printing potential eqn."""
        typed_ethane.save("ethane.lammps")
        with open("ethane.lammps", "r") as f:
            lines = f.readlines()

        stylesDict = {
            "Pair": "4*epsilon*(-sigma**6/r**6+sigma**12/r**12)",
            "Bond": "#LAMMPSHarmonicBondPotential",
            "Angle": "#LAMMPSHarmonicAnglePotential",
            "Dihedral": "#OPLSTorsionPotential",
            "Improper": "#HarmonicImproperPotential",
        }
        for i, line in enumerate(lines):
            if "Coeffs" in line:
                styleLine = lines[i].split()
                if styleLine[0] == "Pair":
                    assert "".join(styleLine[-3:]) == stylesDict[styleLine[0]]
                else:
                    assert styleLine[-1] == stylesDict[styleLine[0]]

    def test_lj_passed_units(self, typed_ethane):
        largest_eps = max(
            list(
                map(
                    lambda x: x.parameters["epsilon"],
                    typed_ethane.atom_types,
                )
            )
        )
        typed_ethane.save(
            "ethane.lammps",
            unit_style="lj",
            lj_cfactorsDict={"energy": largest_eps * 2},
        )
        with open("ethane.lammps", "r") as f:
            lines = f.readlines()
        start = 0
        end = 1
        for i in range(len(lines)):
            if "Pair Coeffs" in lines[i]:
                start = i
            if start > 0 and lines[i] == "\n":
                end = i
                break
        largest_eps_written = max(
            [
                obj
                for obj in map(
                    lambda x: float(x.split()[1]), lines[start + 2 : end]
                )
            ]
        )
        assert largest_eps_written == 0.5

    def test_unit_style_factor(self):
        from gmso.utils.units import LAMMPS_UnitSystems

        for styleStr in [
            "real",
            "metal",
            "si",
            "cgs",
            "electron",
            "micro",
            "nano",
        ]:
            assert (
                LAMMPS_UnitSystems(styleStr).usystem.name
                == "lammps_" + styleStr
            )
        from gmso.exceptions import NotYetImplementedWarning

        with pytest.raises(NotYetImplementedWarning):
            LAMMPS_UnitSystems("None")

    def test_charmm_style(self, harmonic_parmed_types_charmm):
        top = from_parmed(harmonic_parmed_types_charmm)
        top.save("test.lammps")
        assert os.path.exists("test.lammps")

    def test_charmm_improper_ff(self):
        import mbuild as mb

        ff = gmso.ForceField(get_path("tfa_charmm.xml"))

        cpd = mb.load("C(=O)(C(F)(F)F)N", smiles=True)  # TFA molecule
        top = cpd.to_gmso()
        from gmso.parameterization import apply

        ptop = apply(top, ff, identify_connections=True)
        ptop.save("test.lammps", unit_style="real", overwrite=True)
        assert compare_lammps_files("test.lammps", get_path("charmm.lammps"))
