"""Manage Potential functional expressions and variables."""
import warnings

import sympy
import unyt as u

from gmso.utils.decorators import register_pydantic_json
from gmso.utils.misc import unyt_to_hashable

__all__ = ["PotentialExpression"]


@register_pydantic_json(method="json")
class PotentialExpression:
    """A general Expression class with parameters.

    This class is used by `gmso.core.potential.Potential` class and its
    descendants to manage expression and independent variables for the
    Potential class and meant to be used within the scope of aforementioned
    classes. This class can be used in two forms:
        * parametric: With expression, independent_variables and parameters
        * non-parametric: With expression and independent_variables but no parameters

    Non-Parametric form of the class is used by PotentialTemplate class while
    parametric form is used by `gmso.core.AtomType`, `gmso.core.BondType` and so forth.

    Parameters
    ----------
    expression: str or sympy.Expr
        The string expression

    independent_variables: str or sympy.Symbol or list or set thereof
        The independent variables in the expression of the potential.

    parameters: dict, default=None
        A dictionary of parameter whose key is a string and values are parameters
    """

    __slots__ = (
        "_parameters",
        "_expression",
        "_independent_variables",
        "_is_parametric",
    )

    def __init__(self, expression, independent_variables, parameters=None):
        self._expression = self._validate_expression(expression)
        self._independent_variables = self._validate_independent_variables(
            independent_variables
        )
        self._is_parametric = False

        if parameters is not None:
            self._is_parametric = True
            self._parameters = self._validate_parameters(parameters)
            self._verify_validity(
                self._expression, self._independent_variables, self._parameters
            )
        else:
            self._verify_validity(
                self._expression, self.independent_variables, None
            )

    @property
    def expression(self):
        """Return the expression of the functional form."""
        return self._expression

    @expression.setter
    def expression(self, expr):
        """Set the functional form's expression."""
        self.set(expression=expr)

    @property
    def parameters(self):
        """Return the parameters of the expression, if applicable."""
        if not self._is_parametric:
            raise AttributeError(
                f"Object of type {self.__class__.__name__} has no attribute parameters"
            )
        return self._parameters

    @parameters.setter
    def parameters(self, new_params):
        """Set the parameters of the expression, if applicable."""
        if not self._is_parametric:
            raise AttributeError(
                f"Object of type {self.__class__.__name__} has no attribute parameters"
            )
        self.set(parameters=new_params)

    @property
    def independent_variables(self):
        """Return the independent variables involved in the expression."""
        return self._independent_variables

    @independent_variables.setter
    def independent_variables(self, indep_vars):
        """Set the independent variables for the potential expression."""
        self.set(independent_variables=indep_vars)

    @property
    def is_parametric(self):
        """Is the expression parametric."""
        return self._is_parametric

    def set(self, expression=None, parameters=None, independent_variables=None):
        """Set the expression, parameters, and independent variables for this potential.

        Parameters
        ----------
        expression: sympy.Expression or string
            The mathematical expression corresponding to the potential
            If None, the expression remains unchanged
        parameters: dict
            {parameter: value} in the expression
            If None, the parameters remain unchanged
        independent_variables: str or sympy.Symbol or list or set thereof

        Notes
        -----
        Be aware of the symbols used in the `expression` and `parameters`.
        If unnecessary parameters are supplied, an error is thrown.
        If only a subset of the parameters are supplied, they are updated
        while the non-passed parameters default to the existing values
        """
        if expression is not None:
            expression = self._validate_expression(expression)
        else:
            expression = self._expression

        if independent_variables is not None:
            independent_variables = self._validate_independent_variables(
                independent_variables
            )
        else:
            independent_variables = self.independent_variables

        if self._is_parametric:
            if parameters is not None:
                parameters = self._validate_parameters(parameters)
                total_free_symbols = self._expression.free_symbols.union(
                    expression.free_symbols
                )
                for key in list(parameters.keys()):
                    if sympy.Symbol(key) not in total_free_symbols:
                        parameters.pop(key)
                if len(parameters) == 0:
                    raise ValueError(
                        f"`parameters` argument includes no "
                        f"variables found in expression. Expected "
                        f"at least one of {self._parameters.keys()}"
                    )
                for key in self._parameters:
                    if key not in parameters:
                        parameters[key] = self._parameters[key]
            else:
                parameters = self._parameters

            self._verify_validity(expression, independent_variables, parameters)

            self._parameters.update(parameters)
        else:
            self._verify_validity(expression, independent_variables)

        self._expression = expression
        self._independent_variables = independent_variables

        if self._is_parametric:
            for key in list(self._parameters.keys()):
                if sympy.Symbol(key) not in self.expression.free_symbols:
                    self._parameters.pop(key)

    def __hash__(self):
        """Return hash of the potential expression."""
        if self._is_parametric:
            return hash(
                tuple(
                    (
                        self.expression,
                        tuple(self.independent_variables),
                        tuple(self.parameters.keys()),
                        tuple(
                            unyt_to_hashable(val)
                            for val in self.parameters.values()
                        ),
                    )
                )
            )
        else:
            return hash(
                tuple((self.expression, tuple(self.independent_variables)))
            )

    def __eq__(self, other):
        """Determine if two expressions are equivalent."""
        return hash(self) == hash(other)

    def __repr__(self):
        """Representation of the potential expression."""
        descr = list(f"<PotentialExpression, ")
        descr.append(f"expression: {self.expression}, ")
        descr.append(
            f"{len(self.independent_variables)} independent variables>"
        )

        return "".join(descr)

    @staticmethod
    def _validate_expression(expression):
        """Check to see that an expression is a valid sympy expression."""
        if expression is None or isinstance(expression, sympy.Expr):
            pass
        elif isinstance(expression, str):
            expression = sympy.sympify(expression)
        else:
            raise ValueError(
                "Please enter a string, sympy expression or None for expression"
            )

        return expression

    @staticmethod
    def _validate_parameters(parameters):
        """Check to see that parameters is a valid dictionary with units."""
        if not isinstance(parameters, dict):
            raise ValueError("Please enter a dictionary for parameters")
        for key, val in parameters.items():
            if isinstance(val, list):
                for params in val:
                    if not isinstance(params, u.unyt_array):
                        raise ValueError(
                            "Parameter value {} lacks a unyt".format(val)
                        )
            else:
                if not isinstance(val, u.unyt_array):
                    raise ValueError(
                        "Parameter value {} lacks a unyt".format(val)
                    )
            if not isinstance(key, str):
                raise ValueError("Parameter key {} is not a str".format(key))

        return parameters

    @staticmethod
    def _validate_independent_variables(indep_vars):
        """Check to see that independent_variables is a set of valid sympy symbols."""
        if isinstance(indep_vars, str):
            indep_vars = {sympy.symbols(indep_vars)}
        elif isinstance(indep_vars, sympy.Symbol):
            indep_vars = {indep_vars}
        elif isinstance(indep_vars, (list, set)):
            if all([isinstance(val, sympy.Symbol) for val in indep_vars]):
                pass
            elif all([isinstance(val, str) for val in indep_vars]):
                indep_vars = set([sympy.symbols(val) for val in indep_vars])
            else:
                raise ValueError(
                    "`independent_variables` argument was a list "
                    "or set of mixed variables. Please enter a "
                    "list or set of either only strings or only "
                    "sympy symbols"
                )
        else:
            raise ValueError(
                "Please enter a string, sympy expression, "
                "list or set thereof for independent_variables"
            )

        return indep_vars

    @staticmethod
    def json(potential_expression):
        """Convert the provided potential expression to a json serializable dictionary."""
        if not isinstance(potential_expression, PotentialExpression):
            raise TypeError(
                f"{potential_expression} is not of type _PotentialExpression"
            )
        else:
            json_dict = {
                "expression": str(potential_expression.expression),
                "independent_variables": list(
                    str(idep)
                    for idep in potential_expression.independent_variables
                ),
            }
            if potential_expression.is_parametric:
                json_dict["parameters"] = potential_expression.parameters

        return json_dict

    @staticmethod
    def _verify_validity(
        expression, independent_variables_symbols, parameters=None
    ):
        """Verify whether or not the parameters, independent_variables and expression are consistent."""
        for sym in independent_variables_symbols:
            if sym not in expression.free_symbols:
                raise ValueError(
                    f"symbol {sym} is not in expression's free symbols. "
                    f"Cannot use an independent variable which doesn't "
                    f"exist in the expression's free symbols {expression.free_symbols}"
                )
        if parameters is not None:
            parameter_symbols = sympy.symbols(set(parameters.keys()))
            used_symbols = parameter_symbols.union(
                independent_variables_symbols
            )
            unused_symbols = expression.free_symbols - used_symbols
            if len(unused_symbols) > 0:
                warnings.warn(
                    f"You supplied parameters with "
                    f"unused symbols {unused_symbols}"
                )

            if used_symbols != expression.free_symbols:
                symbols = sympy.symbols(set(parameters.keys()))
                if symbols != expression.free_symbols:
                    missing_syms = (
                        expression.free_symbols
                        - symbols
                        - independent_variables_symbols
                    )
                    if missing_syms:
                        raise ValueError(
                            f"Missing necessary dependencies to evaluate "
                            f"potential expression. Missing symbols: {missing_syms}"
                        )

                    extra_syms = symbols ^ expression.free_symbols
                    warnings.warn(
                        f"Potential expression and parameter symbols do not agree, "
                        f"extraneous symbols: {extra_syms}"
                    )
