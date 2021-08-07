import pytest
import sympy

from gmso.lib.potential_templates import PotentialTemplate
from gmso.tests.base_test import BaseTest


class TestTemplate(BaseTest):
    def test_potential_template(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
        )

        assert template.expression == sympy.sympify("a*x+b")

        assert (
            template.expression.free_symbols - template.independent_variables
            is not None
        )

    def test_template_set_expression(self):
        template = PotentialTemplate(
            expression="a*x+b",
            independent_variables={"x"},
        )
        with pytest.raises(NotImplementedError):
            template.set_expression(expression="a*y+b")
