import pytest
import sympy
import unyt as u
from unyt.testing import assert_allclose_units
from copy import deepcopy
from mbuild.tests.base_test import BaseTest
from gmso.utils.conversions import check_convert_kelvin_to_energy_units


class TestKelvinToEnergy(BaseTest):
    def test_K_to_kcal(self):
        input_value = 1 * u.Kelvin / u.nm ** 2
        new_value = check_convert_kelvin_to_energy_units(
            input_value,
            'kcal/mol',
        )

        assert new_value == u.unyt_quantity(0.0019872041457050975, "kcal/(mol*nm**2)")

    def test_kcal_per_mol_to_kJ_per_mol(self):
        input_value = 2 * u.kcal / u.mol * u.gram**2
        new_value = check_convert_kelvin_to_energy_units(
            input_value,
            'kJ/mol',
        )

        assert new_value == u.unyt_quantity(2, "kcal/mol*g**2")

    def test_input_not_unyt_units(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_input_unyt value is a <class 'float'>, "
                  r"not a <class 'unyt.unit_object.Unit'>."
        ):
            input_value = 2.0
            check_convert_kelvin_to_energy_units(
                input_value,
                'kJ/mol',
            )

    def test_kcal_per_mol_to_float_output(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str value is a <class 'float'>, "
                  r"not a <class 'str'>."
        ):
            input_value = 2 * u.kcal / u.mol * u.gram ** 2
            check_convert_kelvin_to_energy_units(
                input_value,
                1.0,
            )

    def test_output_units_in_K(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str can not be in K energy units."
        ):
            input_value = 2 * u.kcal / u.mol * u.gram ** 2
            check_convert_kelvin_to_energy_units(
                input_value,
                'K',
            )

    def test_kcal_per_mol_to_string_m(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: The entered energy_output_unyt_units_str value must be in units of energy/mol, "
                  r"\(length\)\**2\*\(mass\)/\(time\)\**2, but not in K energy units."
        ):
            input_value = 2 * u.kcal / u.mol * u.gram ** 2
            check_convert_kelvin_to_energy_units(
                input_value,
                'm',
            )
