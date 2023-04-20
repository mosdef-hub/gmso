import pytest
import sympy
import unyt as u

from gmso.tests.base_test import BaseTest
from gmso.utils.expression import PotentialExpression, _are_equal_parameters


class TestExpression(BaseTest):
    def test_expression(self):
        expression = PotentialExpression(
            expression="a*x+b",
            independent_variables="x",
            parameters={"a": 1.0 * u.dimensionless, "b": 2.0 * u.dimensionless},
        )

        assert expression.expression == sympy.sympify("a*x+b")
        assert "a" in expression.parameters.keys()
        assert "b" in expression.parameters.keys()
        assert expression.parameters["a"] == 1.0 * u.dimensionless
        assert expression.parameters["b"] == 2.0 * u.dimensionless

    def test_expression_multiple_indep_vars(self):
        expression = PotentialExpression(
            expression="a^2+2*a*b+b^2+2*theta*phi",
            independent_variables={"theta", "phi"},
            parameters={"a": 2.0 * u.nm, "b": 2.0 * u.rad},
        )
        theta = sympy.Symbol("theta")
        phi = sympy.Symbol("phi")
        assert theta in expression.independent_variables
        assert phi in expression.independent_variables
        assert theta in expression.expression.free_symbols
        assert phi in expression.expression.free_symbols

    def test_invalid_expression(self):
        with pytest.raises(ValueError) as e:
            expression = PotentialExpression(
                expression="a*x+b",
                independent_variables="x",
                parameters={"sigma": 1.0 * u.nm, "phi": 1.0 * u.rad},
            )
            assert (
                "ValueError: Missing necessary dependencies to "
                "evaluate potential expression. Missing symbols: {b, a}" in e
            )

    def test_invalid_indep_vars(self):
        with pytest.raises(ValueError) as e:
            expression = PotentialExpression(
                expression="a*x+b", independent_variables="j", parameters=None
            )
            assert (
                "symbol j is not in expression's free symbols Cannot "
                "use an independent variable that doesn't exist in the "
                "expression's free symbols {x, a, b}" in e
            )

    def test_non_parametric_expression(self):
        expression = PotentialExpression(
            expression="a^2+2*a*b+b^2",
            independent_variables="a",
            parameters=None,
        )
        assert expression.is_parametric is False
        with pytest.raises(AttributeError) as e:
            assert expression.parameters
            assert (
                "Object of type _PotentialExpression "
                "has no attribute parameters" in e
            )

    def test_set_indep_variables(self):
        expression = PotentialExpression(
            expression="a^2+2*a*b+b^2",
            independent_variables="a",
            parameters=None,
        )
        expression.independent_variables = {"b"}
        assert sympy.Symbol("b") in expression.independent_variables
        assert sympy.Symbol("a") not in expression.independent_variables

    def test_set_indep_variables_invalid(self):
        expression = PotentialExpression(
            expression="a^2+2*a*b+b^2",
            independent_variables="a",
            parameters=None,
        )

        with pytest.raises(ValueError) as e:
            expression.independent_variables = "y"

        assert expression.independent_variables == {sympy.Symbol("a")}

    def test_set_expression(self):
        expression = PotentialExpression(
            "a^x + b^y + c^z",
            independent_variables={"x", "y", "z"},
            parameters={"a": 2.6 * u.nm, "b": 2.7 * u.nm, "c": 22.8 * u.hertz},
        )
        expression.expression = "a^(2*x) + b^(2*y) + c^(2*z)"
        assert sympy.Symbol("x") in expression.independent_variables

    def test_set_expression_invalid(self):
        expression = PotentialExpression(
            "a^x + b^y + c^z",
            independent_variables={"x", "y", "z"},
            parameters={"a": 2.6 * u.nm, "b": 2.7 * u.nm, "c": 22.8 * u.hertz},
        )
        with pytest.raises(ValueError) as e:
            expression.expression = "2 * theta^2 + 3 * phi^2"

        assert sympy.sympify("a^x + b^y + c^z") == expression.expression

    def test_set_parameters(self):
        expression = PotentialExpression(
            "a^x + b^y + c^z",
            independent_variables={"x", "y", "z"},
            parameters={"a": 2.6 * u.nm, "b": 2.7 * u.nm, "c": 22.8 * u.hertz},
        )
        expression.parameters = {
            "a": 2.7 * u.nm,
            "b": 2.8 * u.nm,
            "c": 220.0 * u.hertz,
        }
        assert expression.parameters["a"] == u.unyt_quantity(2.7, units="nm")
        assert expression.parameters["b"] == u.unyt_quantity(2.8, units="nm")
        assert expression.parameters["c"] == u.unyt_quantity(
            220.0, units="hertz"
        )

    def test_set_parameters_extra(self):
        expression = PotentialExpression(
            "a^x + b^y + c^z",
            independent_variables={"x", "y", "z"},
            parameters={"a": 2.6 * u.nm, "b": 2.7 * u.nm, "c": 22.8 * u.hertz},
        )

        expression.parameters = {
            "a": 2.7 * u.nm,
            "b": 2.8 * u.nm,
            "c": 220.0 * u.hertz,
            "d": 229.0 * u.hertz,
        }

        assert expression.parameters["a"] == u.unyt_quantity(2.7, units="nm")
        assert expression.parameters["b"] == u.unyt_quantity(2.8, units="nm")
        assert expression.parameters["c"] == u.unyt_quantity(
            220.0, units="hertz"
        )
        assert "d" not in expression.parameters

    def test_set_parameters_invalid(self):
        expression = PotentialExpression(
            "a^x + b^y + c^z",
            independent_variables={"x", "y", "z"},
            parameters={"a": 2.6 * u.nm, "b": 2.7 * u.nm, "c": 22.8 * u.hertz},
        )
        with pytest.raises(ValueError):
            expression.parameters = {
                "l": 2.7 * u.nm,
                "m": 2.8 * u.nm,
                "n": 220.0 * u.hertz,
            }

        assert expression.parameters["a"] == u.unyt_quantity(2.6, units="nm")
        assert expression.parameters["b"] == u.unyt_quantity(2.7, units="nm")
        assert expression.parameters["c"] == u.unyt_quantity(
            22.8, units="hertz"
        )
        assert "l" not in expression.parameters

    def test_expression_equality(self):
        expression_1 = PotentialExpression(
            expression="exp(2)+exp(4)+2*phi", independent_variables={"phi"}
        )

        expression_2 = PotentialExpression(
            expression="exp(4) + exp(2) + phi*2", independent_variables={"phi"}
        )

        expression_3 = PotentialExpression(
            expression="exp(4) + exp(2) + phi * 8",
            independent_variables={"phi"},
        )

        assert expression_1.expression == expression_2.expression
        assert expression_3 != expression_2
        assert expression_1 != expression_3

    def test_parametric_equality(self):
        expression_1 = PotentialExpression(
            expression="e^2+e^4+2*phi",
            independent_variables={"phi"},
            parameters={"e": 2.2400 * u.dimensionless},
        )

        expression_2 = PotentialExpression(
            expression="e^4 + e^2 + phi*2",
            independent_variables={"phi"},
            parameters={"e": 2.2400 * u.dimensionless},
        )

        expression_3 = PotentialExpression(
            expression="e^4 + e^2 + phi * 8",
            independent_variables={"phi"},
            parameters={"e": 2.2400 * u.dimensionless},
        )

        assert expression_1.expression == expression_2.expression
        assert expression_3 != expression_2
        assert expression_1 != expression_3

    def test_clone(self):
        expr = PotentialExpression(
            expression="a^2+2*a*b+b^2+2*theta*phi",
            independent_variables={"theta", "phi"},
            parameters={"a": 2.0 * u.nm, "b": 2.0 * u.rad},
        )

        expr_clone = expr.clone()

        assert expr_clone.expression == expr.expression
        assert id(expr_clone.expression) != id(expr.expression)

        assert expr_clone.parameters == expr.parameters
        assert id(expr_clone.parameters) != id(expr.parameters)

        assert expr_clone.independent_variables == expr.independent_variables
        assert id(expr_clone.independent_variables) != id(
            expr.independent_variables
        )

        assert expr == expr_clone

    def test_from_non_parametric(self):
        non_parametric = PotentialExpression(
            expression="x**2+z*x*y+y**2", independent_variables={"z"}
        )

        parametric = PotentialExpression.from_non_parametric(
            non_parametric,
            parameters={"x": 2.9 * u.dimensionless, "y": 10000 * u.m},
            valid=True,
        )

        assert parametric.expression == non_parametric.expression
        assert id(parametric.expression) != id(non_parametric.expression)
        assert (
            parametric.independent_variables
            == non_parametric.independent_variables
        )
        parametric.independent_variables.add("X")
        assert (
            parametric.independent_variables
            != non_parametric.independent_variables
        )

    def test_from_non_parametric_errors(self):
        with pytest.raises(
            TypeError,
            match="Expected <object object at .*> to be of type "
            "<class 'gmso.utils.expression.PotentialExpression'> "
            "but found <class 'object'>.",
        ):
            parametric = PotentialExpression.from_non_parametric(
                non_parametric=object(), parameters={}
            )

        non_parametric = PotentialExpression(
            expression="x**2+z*x*y+y**2", independent_variables={"z"}
        )

        parametric = PotentialExpression.from_non_parametric(
            non_parametric,
            parameters={"x": 2.9 * u.dimensionless, "y": 10000 * u.m},
            valid=False,
        )

        with pytest.raises(
            ValueError,
            match="Cannot create a parametric expression from a parametric expression.",
        ):
            PotentialExpression.from_non_parametric(parametric, parameters={})

        with pytest.raises(
            ValueError,
            match="Missing necessary dependencies to evaluate potential expression. Missing symbols: {y}",
        ):
            parametric = PotentialExpression.from_non_parametric(
                non_parametric,
                parameters={"x": 2.9 * u.dimensionless, "z": 10000 * u.m},
                valid=False,
            )

        parametric = PotentialExpression.from_non_parametric(
            non_parametric,
            parameters={"x": 2.9 * u.dimensionless, "z": 10000 * u.m},
            valid=True,
        )

    def test_clone_with_unyt_arrays(self):
        expression = PotentialExpression(
            expression="x**2 + y**2 + 2*x*y*theta",
            independent_variables="theta",
            parameters={
                "x": [2.0, 4.5] * u.nm,
                "y": [3.4, 4.5] * u.kcal / u.mol,
            },
        )

        expression_clone = expression.clone()
        assert expression_clone == expression

    def test_expression_equality_different_params(self):
        expr1 = PotentialExpression(
            independent_variables="r",
            parameters={"a": 2.0 * u.nm, "b": 3.0 * u.nm},
            expression="a+r*b",
        )

        expr2 = PotentialExpression(
            independent_variables="r",
            parameters={"c": 2.0 * u.nm, "d": 3.0 * u.nm},
            expression="c+r*d",
        )

        assert expr1 != expr2

    def test_expression_equality_same_params_different_values(self):
        expr1 = PotentialExpression(
            independent_variables="r",
            parameters={"a": 2.0 * u.nm, "b": 3.0 * u.nm},
            expression="a+r*b",
        )

        expr2 = PotentialExpression(
            independent_variables="r",
            parameters={"a": 2.0 * u.nm, "b": 3.5 * u.nm},
            expression="a+r*b",
        )

        assert expr1 != expr2

    def test_are_equal_parameters(self):
        u1 = {"a": 2.0 * u.nm, "b": 3.5 * u.nm}
        u2 = {"c": 2.0 * u.nm, "d": 3.5 * u.nm}
        assert _are_equal_parameters(u1, u2) is False
