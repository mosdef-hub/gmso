import numpy as np
import pytest
import sympy
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.forcefield import ForceField
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class TestForceFieldFromXML(BaseTest):
    def test_carbon_force_field(self):
        carbon = ForceField(get_path("carbon.xml"))

        assert len(carbon.atom_types) == 1
        assert len(carbon.bond_types) == 1
        assert len(carbon.angle_types) == 1
        assert len(carbon.dihedral_types) == 1

        # Store expected expressions in list
        ref_exprs = [
            sympy.sympify(expr)
            for expr in [
                "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
                "0.5 * k * (r-r_eq)**2",
                "0.5 * k * (theta-theta_eq)**2",
                "k * (1 + cos(n * theta - theta_0))",
            ]
        ]

        assert carbon.atom_types["C"].charge.value == 0

        assert (
            sympy.simplify(carbon.atom_types["C"].expression - ref_exprs[0])
            == 0
        )
        assert (
            sympy.simplify(carbon.bond_types["C~C"].expression - ref_exprs[1])
            == 0
        )
        assert (
            sympy.simplify(
                carbon.angle_types["C~C~C"].expression - ref_exprs[2]
            )
            == 0
        )
        assert (
            sympy.simplify(
                carbon.dihedral_types["*~C~C~*"].expression - ref_exprs[3]
            )
            == 0
        )

        assert_allclose_units(
            carbon.atom_types["C"].parameters["sigma"],
            0.339966950842 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            carbon.atom_types["C"].parameters["epsilon"],
            0.359824 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            carbon.bond_types["C~C"].parameters["r_eq"],
            0.1324 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            carbon.bond_types["C~C"].parameters["k"],
            493460.96 * u.Unit("kJ/(mol*nm**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            carbon.angle_types["C~C~C"].parameters["theta_eq"],
            2.12598556185 * u.radian,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            carbon.angle_types["C~C~C"].parameters["k"],
            584.42112 * u.Unit("kJ/(mol*rad**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            carbon.dihedral_types["*~C~C~*"].parameters["k"],
            27.8236 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            carbon.dihedral_types["*~C~C~*"].parameters["n"],
            2 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            carbon.dihedral_types["*~C~C~*"].parameters["theta_0"],
            np.pi * u.radian,
            rtol=1e-5,
            atol=1e-8,
        )

    def test_tip3p_force_field(self):
        water = ForceField(get_path("tip3p.xml"))
        assert len(water.atom_types) == 2
        assert len(water.bond_types) == 1
        assert len(water.angle_types) == 1
        assert len(water.dihedral_types) == 0

        # Store expected expressions in list
        ref_exprs = [
            sympy.sympify(expr)
            for expr in [
                "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
                "0.5 * k * (r-r_eq)**2",
                "0.5 * k * (theta-theta_eq)**2",
            ]
        ]

        assert water.atom_types["opls_111"].charge.value == -0.834
        assert water.atom_types["opls_112"].charge.value == 0.417

        assert (
            sympy.simplify(
                water.atom_types["opls_111"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                water.atom_types["opls_112"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                water.bond_types["opls_111~opls_112"].expression - ref_exprs[1]
            )
            == 0
        )
        assert (
            sympy.simplify(
                water.angle_types["opls_112~opls_111~opls_112"].expression
                - ref_exprs[2]
            )
            == 0
        )

        assert_allclose_units(
            water.atom_types["opls_111"].parameters["sigma"],
            0.315061 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            water.atom_types["opls_111"].parameters["epsilon"],
            0.636386 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            water.atom_types["opls_112"].parameters["sigma"],
            1.0 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            water.atom_types["opls_112"].parameters["epsilon"],
            0.0 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            water.bond_types["opls_111~opls_112"].parameters["r_eq"],
            0.09572 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            water.bond_types["opls_111~opls_112"].parameters["k"],
            502416.0 * u.Unit("kJ/(mol*nm**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            water.angle_types["opls_112~opls_111~opls_112"].parameters[
                "theta_eq"
            ],
            1.824218134 * u.radian,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            water.angle_types["opls_112~opls_111~opls_112"].parameters["k"],
            682.02 * u.Unit("kJ/(mol*rad**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_spce_xml(self):
        spce = ForceField(get_path("spce.xml"))

        assert len(spce.atom_types) == 2
        assert len(spce.bond_types) == 1
        assert len(spce.angle_types) == 1
        assert len(spce.dihedral_types) == 0

        # Store expected expressions in list
        ref_exprs = [
            sympy.sympify(expr)
            for expr in [
                "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
                "0.5 * k * (r-r_eq)**2",
                "0.5 * k * (theta-theta_eq)**2",
            ]
        ]

        assert_allclose_units(
            spce.atom_types["opls_116"].charge,
            -0.8476 * u.elementary_charge,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.atom_types["opls_117"].charge,
            0.4238 * u.elementary_charge,
            rtol=1e-5,
            atol=1e-8,
        )

        assert (
            sympy.simplify(
                spce.atom_types["opls_116"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                spce.atom_types["opls_117"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                spce.bond_types["opls_116~opls_117"].expression - ref_exprs[1]
            )
            == 0
        )
        assert (
            sympy.simplify(
                spce.angle_types["opls_117~opls_116~opls_117"].expression
                - ref_exprs[2]
            )
            == 0
        )

        assert_allclose_units(
            spce.atom_types["opls_116"].parameters["sigma"],
            0.316557 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.atom_types["opls_116"].parameters["epsilon"],
            0.650194 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.atom_types["opls_117"].parameters["sigma"],
            0.1 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.atom_types["opls_117"].parameters["epsilon"],
            0.0 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            spce.bond_types["opls_116~opls_117"].parameters["r_eq"],
            0.1 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.bond_types["opls_116~opls_117"].parameters["k"],
            345000.0 * u.Unit("kJ/(mol*nm**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            spce.angle_types["opls_117~opls_116~opls_117"].parameters[
                "theta_eq"
            ],
            109.47 * u.degree,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            spce.angle_types["opls_117~opls_116~opls_117"].parameters["k"],
            383.0 * u.Unit("kJ/mol/rad**2"),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_noble_mie_xml(self):
        ff = ForceField(get_path("noble_mie.xml"))
        templates = PotentialTemplateLibrary()
        ref_expr = templates["MiePotential"].expression

        assert len(ff.atom_types) == 4
        assert len(ff.bond_types) == 0
        assert len(ff.angle_types) == 0
        assert len(ff.dihedral_types) == 0

        for name, atom_type in ff.atom_types.items():
            assert sympy.simplify(atom_type.expression - ref_expr) == 0

        assert_allclose_units(
            ff.atom_types["Ne"].parameters["epsilon"],
            0.26855713 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ne"].parameters["sigma"],
            2.794 * u.Angstrom,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ne"].parameters["n"],
            11 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ne"].parameters["m"],
            6 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert ff.atom_types["Ne"].charge.value == 0

        assert_allclose_units(
            ff.atom_types["Ar"].parameters["epsilon"],
            1.01519583 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ar"].parameters["sigma"],
            3.405 * u.Angstrom,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ar"].parameters["n"],
            13 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Ar"].parameters["m"],
            6 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert ff.atom_types["Ar"].charge.value == 0

        assert_allclose_units(
            ff.atom_types["Kr"].parameters["epsilon"],
            1.46417678 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Kr"].parameters["sigma"],
            3.645 * u.Angstrom,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Kr"].parameters["n"],
            14 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Kr"].parameters["m"],
            6 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert ff.atom_types["Kr"].charge.value == 0

        assert_allclose_units(
            ff.atom_types["Xe"].parameters["epsilon"],
            2.02706587 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Xe"].parameters["sigma"],
            3.964 * u.Angstrom,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Xe"].parameters["n"],
            14 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ff.atom_types["Xe"].parameters["m"],
            6 * u.dimensionless,
            rtol=1e-5,
            atol=1e-8,
        )
        assert ff.atom_types["Xe"].charge.value == 0

    def test_ethylene_forcefield(self):
        ethylene = ForceField(get_path("ethylene.xml"))

        assert len(ethylene.atom_types) == 2
        assert len(ethylene.bond_types) == 2
        assert len(ethylene.angle_types) == 2
        assert len(ethylene.dihedral_types) == 1

        # Store expected expressions in list
        ref_exprs = [
            sympy.sympify(expr)
            for expr in [
                "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
                "0.5 * k * (r-r_eq)**2",
                "0.5 * k * (theta-theta_eq)**2",
                "c0 + c1 * cos(phi) + c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + c5 * cos(phi)**5",
            ]
        ]

        assert_allclose_units(
            ethylene.atom_types["opls_143"].charge,
            -0.23 * u.elementary_charge,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.atom_types["opls_144"].charge,
            0.115 * u.elementary_charge,
            rtol=1e-5,
            atol=1e-8,
        )

        assert (
            sympy.simplify(
                ethylene.atom_types["opls_143"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                ethylene.atom_types["opls_144"].expression - ref_exprs[0]
            )
            == 0
        )
        assert (
            sympy.simplify(
                ethylene.bond_types["opls_143~opls_144"].expression
                - ref_exprs[1]
            )
            == 0
        )
        assert (
            sympy.simplify(
                ethylene.angle_types["opls_144~opls_143~opls_144"].expression
                - ref_exprs[2]
            )
            == 0
        )
        assert (
            sympy.simplify(
                ethylene.dihedral_types[
                    "opls_144~opls_143~opls_143~opls_144"
                ].expression
                - ref_exprs[3]
            )
            == 0
        )
        assert_allclose_units(
            ethylene.atom_types["opls_143"].parameters["sigma"],
            0.317984 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.atom_types["opls_143"].parameters["epsilon"],
            0.355 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.atom_types["opls_144"].parameters["sigma"],
            0.12552 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.atom_types["opls_144"].parameters["epsilon"],
            0.242 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            ethylene.bond_types["opls_143~opls_143"].parameters["r_eq"],
            0.134 * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.bond_types["opls_143~opls_144"].parameters["k"],
            284512.0 * u.Unit("kJ/(mol*nm**2)"),
            rtol=1e-5,
            atol=1e-8,
        )

        assert_allclose_units(
            ethylene.angle_types["opls_143~opls_143~opls_144"].parameters[
                "theta_eq"
            ],
            2.0943951 * u.rad,
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.angle_types["opls_144~opls_143~opls_144"].parameters["k"],
            292.88 * u.Unit("kJ/mol/rad**2"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.dihedral_types[
                "opls_144~opls_143~opls_143~opls_144"
            ].parameters["c0"],
            58.576 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.dihedral_types[
                "opls_144~opls_143~opls_143~opls_144"
            ].parameters["c2"],
            -58.576 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            ethylene.dihedral_types[
                "opls_144~opls_143~opls_143~opls_144"
            ].parameters["c5"],
            0.0 * u.Unit("kJ/mol"),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_error_duplicated_types(self):
        # Temporarily opt out, pending new forcefield-utilities release
        # with pytest.raises(ValueError) as e:
        with pytest.raises(Exception) as e:
            ForceField(get_path("ff-nonunique-dihedral.xml"))
            # assert (
            #     e
            #     == "Duplicate identifier found for DihedralTypes: ('CT', 'CT', 'CT', 'HC')"
            # )
