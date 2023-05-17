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
from gmso.external import to_parmed
from gmso.formats.formats_registry import UnsupportedFileFormatError
from gmso.formats.lammpsdata import read_lammpsdata, write_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


def compare_lammps_files(fn1, fn2, skip_linesList=[]):
    """Check for line by line equality between lammps files, by any values."""
    with open(fn1, "r") as f:
        line1 = f.readlines()
    with open(fn2, "r") as f:
        line2 = f.readlines()
    for lnum, (l1, l2) in enumerate(zip(line1, line2)):
        print(l1, l2)
        print("###############")
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
    return True


class TestLammpsWriter(BaseTest):
    @pytest.mark.parametrize(
        "fname", ["data.lammps", "data.data", "data.lammpsdata"]
    )
    def test_write_lammps(self, fname, typed_ar_system):
        typed_ar_system.save(fname)

    def test_write_lammps_triclinic(self, typed_ar_system):
        typed_ar_system.box = Box(lengths=[1, 1, 1], angles=[60, 90, 120])
        typed_ar_system.save("triclinic.lammps")

    def test_ethane_lammps(self, typed_ethane):
        typed_ethane.save("ethane.lammps")

    def test_opls_lammps(self, typed_ethane_opls):
        typed_ethane_opls.save("ethane.lammps")

    def test_water_lammps(self, typed_water_system):
        typed_water_system.save("data.lammps")

    def test_read_lammps(self, filename=get_path("data.lammps")):
        top = gmso.Topology.load(filename)

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
            "typed_water_system",
        ],
    )
    def test_lammps_vs_parmed_by_mol(self, top, request):
        """Parmed LAMMPSDATA Compare outputs.

        atom_style = 'full', 'atomic', 'charge', 'molecular'
        unit_style = 'real', 'lj'
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
        )

    # TODO: Test parameters that have intraconversions between them
    def test_lammps_default_conversions(self, typed_ethane):
        """Test for parameter intraconversions with potential styles.

        These include:
        bonds: factor of 2 harmonic k
        angles: factor of 2 harmonic k
        dihedrals: RB torsions to OPLS
        impropers:
        pairs:
        additional: All styles to zero and none
        """
        typed_ethane.save("opls.lammps")
        with open("opls.lammps", "r") as f:
            lines = f.readlines()
        assert lines[38:41] == [
            "Dihedral Coeffs #OPLSTorsionPotential\n",
            "#\tk1 (kcal/mol)\tk2 (kcal/mol)\tk3 (kcal/mol)\tk4 (kcal/mol)\n",
            "1\t 0.00000\t-0.00000\t 0.30000\t-0.00000\t# opls_140\topls_135\topls_135\topls_140\n",
        ]

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

    # TODO: Test potential styles that are not supported by parmed
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
        Supporte styles are: ["real", "lj"]
        TODO: ["metal", "si", "cgs", "electron", "micro", "nano"]
        _______References_______
        https://docs.lammps.org/units.html
        """
        # check the initial set of units
        from gmso.formats.lammpsdata import get_units

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
        energy_unit = get_units(unit_style, "energy")
        angle_unit = get_units(unit_style, "angle_eq")
        length_unit = get_units(unit_style, "length")
        charge_unit = get_units(unit_style, "charge")
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

    # TODO: Test for warning handling
    def test_lammps_warnings(self, typed_ethane):
        with pytest.warns(
            UserWarning, match="Call to function write_lammpsdata is WIP."
        ):
            """check for warning about WIP"""
            typed_ethane.save("warning.lammps")

    # TODO: Test for error handling
    from gmso.exceptions import EngineIncompatibilityError

    def test_lammps_errors(self, typed_ethane):
        with pytest.raises(UnsupportedFileFormatError):
            typed_ethane.save("error.lammmps")
