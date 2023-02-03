from copy import deepcopy

import pytest
import sympy
import unyt as u
from mbuild.tests.base_test import BaseTest
from unyt.testing import assert_allclose_units

from gmso.utils.conversions import convert_kelvin_to_energy_units


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
            "0.5 * k0 + 0.5 * k1 * (1 + cos(phi)) + 0.5 * k2 * (1 - cos(2*phi)) + \
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
