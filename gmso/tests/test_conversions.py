from copy import deepcopy

import numpy as np
import pytest
import sympy
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.tests.base_test import BaseTest
from gmso.utils.conversions import (
    convert_kelvin_to_energy_units,
    convert_params_units,
)


def _convert_potential_types(top, connStr, expected_units_dim, base_units):
    potentials = getattr(top, connStr + "_types")
    ref_values = {"energy": "kJ/mol", "length": "nm", "angle": "radians"}
    convert_params_units(potentials, expected_units_dim, base_units, ref_values)
    return potentials


class TestKelvinToEnergy(BaseTest):
    def test_convert_potential_styles(self, typed_ethane):
        from sympy import sympify

        rb_expr = sympify(
            "c0 * cos(phi)**0 + c1 * cos(phi)**1 + c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + c5 * cos(phi)**5"
        )
        assert typed_ethane.dihedrals[0].dihedral_type.expression == rb_expr
        for dihedral in typed_ethane.dihedrals:
            dihedral.dihedral_type.name = "RyckaertBellemansTorsionPotential"
        typed_ethane.convert_potential_styles(
            {"dihedrals": "OPLSTorsionPotential"}
        )
        opls_expr = sympify(
            "0.5 * k1 * (1 + cos(phi)) + 0.5 * k2 * (1 - cos(2*phi)) + \
            0.5 * k3 * (1 + cos(3*phi)) + 0.5 * k4 * (1 - cos(4*phi))"
        )
        assert typed_ethane.dihedrals[0].dihedral_type.expression == opls_expr
        assert (
            typed_ethane.dihedrals[0].dihedral_type.name
            == "OPLSTorsionPotential"
        )

    def test_K_to_kcal(self):
        input_value = 1 * u.Kelvin / u.nm**2
        new_value = convert_kelvin_to_energy_units(
            input_value,
            "kcal/mol",
        )

        assert new_value == u.unyt_quantity(
            0.0019872041457050975, "kcal/(mol*nm**2)"
        )

    def test_kcal_per_mol_to_kJ_per_mol(self):
        input_value = 2 * u.kcal / u.mol * u.gram**2
        new_value = convert_kelvin_to_energy_units(
            input_value,
            "kJ/mol",
        )

        assert new_value == u.unyt_quantity(2, "kcal/mol*g**2")

    def test_input_not_unyt_units(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_input_unyt value is a <class 'float'>, "
            r"not a <class 'unyt.unit_object.Unit'>.",
        ):
            input_value = 2.0
            convert_kelvin_to_energy_units(
                input_value,
                "kJ/mol",
            )

    def test_kcal_per_mol_to_float_output(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str value is a <class 'float'>, "
            r"not a <class 'str'>.",
        ):
            input_value = 2 * u.kcal / u.mol * u.gram**2
            convert_kelvin_to_energy_units(
                input_value,
                1.0,
            )

    def test_output_units_in_K(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str can not be in K energy units.",
        ):
            input_value = 2 * u.kcal / u.mol * u.gram**2
            convert_kelvin_to_energy_units(
                input_value,
                "K",
            )

    def test_kcal_per_mol_to_string_m(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str value must be in units of energy/mol, "
            r"\(length\)\**2\*\(mass\)/\(time\)\**2, but not in K energy units.",
        ):
            input_value = 2 * u.kcal / u.mol * u.gram**2
            convert_kelvin_to_energy_units(
                input_value,
                "m",
            )

    def test_conversion_for_topology_dihedrals(self, typed_ethane):
        expected_units_dim = dict(
            zip(["c0", "c1", "c2", "c3", "c4", "c5"], ["energy"] * 6)
        )
        base_units = u.UnitSystem("atomic", "Å", "mp", "fs", "nK", "rad")
        base_units["energy"] = "kcal/mol"
        potentials = _convert_potential_types(
            typed_ethane, "dihedral", expected_units_dim, base_units
        )
        assert (
            str(
                potentials.iterator[0]
                .dihedral_type.parameters["c0"]
                .units.dimensions
            )
            == "(length)**2*(mass)/(time)**2"
        )
        assert potentials.iterator[0].dihedral_type.parameters[
            "c0"
        ].units == u.Unit("kcal/mol")

    def test_conversion_for_topology_angles(self, typed_ethane):
        expected_units_dim = dict(k="energy/angle**2", theta_eq="angle")
        base_units = u.UnitSystem("atomic", "Å", "mp", "fs", "nK", "rad")
        base_units["energy"] = "kcal/mol"
        potentials = _convert_potential_types(
            typed_ethane, "angle", expected_units_dim, base_units
        )
        assert (
            str(
                potentials.iterator[0]
                .angle_type.parameters["k"]
                .units.dimensions
            )
            == "(length)**2*(mass)/((angle)**2*(time)**2)"
        )
        assert potentials.iterator[0].angle_type.parameters[
            "k"
        ].units == u.Unit("kcal/mol/rad**2")

    def test_conversion_for_topology_bonds(self, typed_ethane):
        expected_units_dim = dict(k="energy/length**2", r_eq="length")
        base_units = u.UnitSystem("atomic", "Å", "mp", "fs", "nK", "rad")
        base_units["energy"] = "kcal/mol"
        potentials = _convert_potential_types(
            typed_ethane, "bond", expected_units_dim, base_units
        )
        assert (
            str(
                potentials.iterator[0]
                .bond_type.parameters["k"]
                .units.dimensions
            )
            == "(mass)/(time)**2"
        )
        assert potentials.iterator[0].bond_type.parameters["k"].units == u.Unit(
            "kcal/mol/angstrom**2"
        )

    def test_conversion_for_topology_sites(self, typed_ethane):
        expected_units_dim = dict(sigma="length", epsilon="energy")
        base_units = u.UnitSystem("atomic", "Å", "mp", "fs", "nK", "rad")
        base_units["energy"] = "kcal/mol"
        potentials = _convert_potential_types(
            typed_ethane, "atom", expected_units_dim, base_units
        )
        assert (
            str(
                potentials.iterator[0]
                .atom_type.parameters["epsilon"]
                .units.dimensions
            )
            == "(length)**2*(mass)/(time)**2"
        )
        assert potentials.iterator[0].atom_type.parameters[
            "epsilon"
        ].units == u.Unit("kcal/mol")

    def test_lammps_dimensions_to_energy(self):
        from gmso.utils.units import LAMMPS_UnitSystems

        units = u.Unit("kg")
        outdims = LAMMPS_UnitSystems._dimensions_to_energy(units.dimensions)
        assert outdims == units.dimensions == u.dimensions.mass
        units = u.Unit("J")
        outdims = LAMMPS_UnitSystems._dimensions_to_energy(units.dimensions)
        assert outdims == sympy.Symbol("(energy)")
        assert (
            units.dimensions
            == u.dimensions.length**2 * u.dimensions.mass / u.dimensions.time**2
        )
        units = u.Unit("kcal/nm")
        outdims = LAMMPS_UnitSystems._dimensions_to_energy(units.dimensions)
        assert outdims == sympy.Symbol("(energy)") / u.dimensions.length
        assert (
            units.dimensions
            == u.dimensions.length * u.dimensions.mass / u.dimensions.time**2
        )

    def test_lammps_conversion_parameters_base_units(self):
        from gmso.utils.units import LAMMPS_UnitSystems

        parameter = 100 * u.Unit("kcal/mol*fs/Å")
        base_unyts = LAMMPS_UnitSystems(
            "real"
        )  # "lammps_real", "Å", "amu", "fs", "K", "rad"
        paramStr = base_unyts.convert_parameter(
            parameter, conversion_factorDict=None
        )
        assert paramStr == "100.000"
        parameter = 100 * u.Unit("K*fs/amu/nm")
        paramStr = base_unyts.convert_parameter(
            parameter, conversion_factorDict=None
        )
        assert paramStr == "10.000"
        parameter = 100 * u.Unit("km*g*ms*kJ*degree")
        base_unyts = LAMMPS_UnitSystems(
            "si"
        )  # "lammps_si", "m", "kg", "s", "K", "rad",
        paramStr = base_unyts.convert_parameter(
            parameter, conversion_factorDict=None, n_decimals=6
        )
        assert paramStr == str(round(100 * np.pi / 180, 6))
        parameter = 1 * u.Unit("g*kJ*Coulomb*m*degree")
        base_unyts = LAMMPS_UnitSystems(
            "si"
        )  # "lammps_si", "m", "kg", "s", "K", "rad"
        paramStr = base_unyts.convert_parameter(
            parameter, conversion_factorDict=None, n_decimals=6
        )
        assert np.isclose(float(paramStr), np.pi / 180, 1e-3)

    def test_lammps_conversion_parameters_lj(self):
        from gmso.utils.units import LAMMPS_UnitSystems

        parameter = 1 * u.Unit("g*kJ*Coulomb*m*degree")
        conversion_factorDict = {
            "mass": 3 * u.Unit("g"),
            "energy": 3 * u.Unit("kJ"),
            "charge": 3 * u.Unit("Coulomb"),
            "length": 3 * u.Unit("m"),
        }
        base_unyts = LAMMPS_UnitSystems("lj")
        paramStr = base_unyts.convert_parameter(
            parameter,
            conversion_factorDict=conversion_factorDict,
            n_decimals=6,
        )
        assert np.isclose(float(paramStr), 1 / 3**4, atol=1e-6)
