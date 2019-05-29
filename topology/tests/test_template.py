import sympy

from topology.lib.potential_templates import PotentialTemplate
from topology.tests.base_test import BaseTest


class TestTempalte(BaseTest):
    def test_potential_template(self):
        template = PotentialTemplate(
            expression='a*x+b',
            independent_variables={'x'},
        )
    
        assert template.expression == sympy.sympify('a*x+b')
        assert template.template
    
        assert template.expression.free_symbols - template.independent_variables is not None
