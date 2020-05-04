import sympy

import pytest

from gmso.lib.potential_templates import PotentialTemplate
from gmso.tests.base_test import BaseTest


class TestTemplate(BaseTest):

    @pytest.fixture()
    def template(self):
        template = PotentialTemplate(
            expression='a*x+b',
            independent_variables={'x'},
        )
        return template

    def test_potential_template(self, template):
    
        assert template.expression == sympy.sympify('a*x+b')

        assert template.expression.free_symbols - template.independent_variables is not None

    def test_potential_template_name_change(self, template):
        with pytest.raises(AttributeError) as e:
            template.name = 'new_name'
            assert 'Properties for a potential template cannot be changed' in str(e.value)

    def test_potential_template_expression_change(self, template):
        with pytest.raises(AttributeError) as e:
            template.expression = 'new_expression'
            assert 'Properties for a potential template cannot be changed' in str(e.value)

    def test_template_independent_variables_change(self, template):
        with pytest.raises(AttributeError) as e:
            template.independent_variables = 'new_variables'
            assert 'Properties for a potential template cannot be changed' in str(e.value)

    def test_call_to_set_expression(self, template):
        with pytest.raises(NotImplementedError):
            template.set_expression()
