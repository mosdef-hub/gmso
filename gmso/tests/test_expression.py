import sympy
import unyt as u

from gmso.utils.expression import _PotentialExpression
from gmso.tests.base_test import BaseTest


class TestExpression(BaseTest):

    def test_expression(self):
        expression = _PotentialExpression(
            expression='a*x+b',
            independent_variables='x',
            parameters={
                'a': 1.0 * u.dimensionless,
                'b': 2.0 * u.dimensionless
            }
        )

        assert expression.expression == sympy.sympify('a*x+b')
        assert 'a' in expression.parameters.keys()
        assert 'b' in expression.parameters.keys()
        assert expression.parameters['a'] == 1.0 * u.dimensionless
        assert expression.parameters['b'] == 2.0 * u.dimensionless

