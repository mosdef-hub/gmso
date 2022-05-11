import re

import pytest
import sympy
import unyt as u

from gmso.exceptions import MissingParameterError, UnknownParameterError
from gmso.lib.potential_templates import PotentialTemplate
from gmso.tests.base_test import BaseTest


class TestTemplate(BaseTest):
    def test_potential_template(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={"a": "energy", "b": "length"},
        )

        assert template.expression == sympy.sympify("a*x+b")
        assert template.expected_parameters_dimensions == {
            "a": u.dimensions.energy,
            "b": u.dimensions.length,
        }

        assert (
            template.expression.free_symbols - template.independent_variables
            is not None
        )

    def test_template_set_expression(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={"a": "length", "b": "length"},
        )
        with pytest.raises(NotImplementedError):
            template.set_expression(expression="a*y+b")

    def test_parameterization_non_dict_expected_dimensions(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={"a": "length", "b": "length"},
        )

        with pytest.raises(TypeError):
            template.assert_can_parameterize_with(object())

    def test_parameterization_unknown_dimension(self):
        with pytest.raises(
            AttributeError,
            match="^module 'unyt.dimensions' has no attribute 'missing'$",
        ):
            invalid_dimension_template = PotentialTemplate(
                expression="a*x+c",
                independent_variables="x",
                expected_parameters_dimensions={"a": "missing", "b": "length"},
            )

    def test_unknown_missing_parameters(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={"a": "energy", "b": "length"},
        )

        with pytest.raises(
            UnknownParameterError,
            match=re.escape(
                "Parameter c is not one of the expected parameters ['a', 'b']"
            ),
        ):
            template.assert_can_parameterize_with(
                {
                    "c": 1.0 * u.kcal / u.mol,
                    "d": 2.0 * u.meter,
                    "a": 1.0 * u.kcal / u.mol,
                    "b": 2.0 * u.meter,
                }
            )

        with pytest.raises(
            MissingParameterError,
            match=re.escape(
                "Parameter 'b' missing from the provided parameters ['a']"
            ),
        ):
            template.assert_can_parameterize_with(
                {
                    "a": 1.0 * u.kcal / u.mol,
                }
            )

    def test_complex_dimensions(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={
                "a": "energy**2/mass**3/length",
                "b": "length**-2/mass**-2",
            },
        )

        template.assert_can_parameterize_with(
            {
                "a": 25 * (u.kcal / u.mol) ** 2 / u.kg**3 / u.meter,
                "b": 50 * (u.gram * u.gram) / (u.nm * u.nm),
            }
        )

    def test_non_unyt_error(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={
                "a": "dimensionless",
                "b": "length",
            },
        )

        with pytest.raises(ValueError):
            template.assert_can_parameterize_with({"a": 1.0, "b": 2.0})

    def test_dimensionless_errors(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
            expected_parameters_dimensions={
                "a": "dimensionless",
                "b": "length",
            },
        )

        with pytest.raises(AssertionError):
            template.assert_can_parameterize_with(
                {"a": 1.0 * u.nm, "b": 2.0 * u.dimensionless}
            )

        with pytest.raises(AssertionError):
            template.assert_can_parameterize_with(
                {"a": 1.0 * u.dimensionless, "b": 2.0 * u.dimensionless}
            )
