import pytest
import sympy
import unyt as u

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
        )
        assert template.expected_parameters_dimensions == {}
        with pytest.raises(NotImplementedError):
            template.set_expression(expression="a*y+b")

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
