"""Manage Potential functional expressions and variables."""

import warnings
from copy import deepcopy
from functools import lru_cache
from typing import Dict

import sympy
import unyt as u

__all__ = ["PotentialExpression"]


def _are_equal_parameters(u1, u2):
    """Compare two parameters of unyt quantities/arrays.

    This method compares two dictionaries (`u1` and `u2`) of
    `unyt_quantities` and returns True if:
        * u1 and u2 have the exact same key set
        * for each key, the value in u1 and u2 have the same unyt quantity
    """
    if u1.keys() != u2.keys():
        return False
    else:
        for k, v in u1.items():
            if not u.allclose_units(v, u2[k]):
                return False

        return True


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

    verify_validity: bool, default=True
        If true verify validity of the expression, parameters and independent variables
    """

    __slots__ = (
        "_parameters",
        "_expression",
        "_independent_variables",
        "_is_parametric",
    )

    def __init__(
        self,
        expression,
        independent_variables,
        parameters=None,
        verify_validity=True,
    ):
        self._expression = (
            self._validate_expression(expression)
            if verify_validity
            else expression
        )
        self._independent_variables = (
            self._validate_independent_variables(independent_variables)
            if verify_validity
            else independent_variables
        )
        self._is_parametric = False

        if parameters is not None:
            self._is_parametric = True
            self._parameters = (
                self._validate_parameters(parameters)
                if verify_validity
                else parameters
            )

        if verify_validity:
            self._verify_validity(
                self._expression,
                frozenset(self._independent_variables),
                frozenset(self._parameters) if self._is_parametric else None,
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

            self._verify_validity(
                expression,
                frozenset(independent_variables),
                frozenset(parameters),
            )

            self._parameters.update(parameters)
        else:
            self._verify_validity(expression, frozenset(independent_variables))

        self._expression = expression
        self._independent_variables = independent_variables

        if self._is_parametric:
            for key in list(self._parameters.keys()):
                if sympy.Symbol(key) not in self.expression.free_symbols:
                    self._parameters.pop(key)

    def __repr__(self):
        """Representation of the potential expression."""
        descr = list(f"<PotentialExpression, ")
        descr.append(f"expression: {self.expression}, ")
        descr.append(
            f"{len(self.independent_variables)} independent variables>"
        )

        return "".join(descr)

    @staticmethod
    @lru_cache(maxsize=128)
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

    def __eq__(self, other):
        """Equality checks for two expressions."""
        if other is self:
            return True
        if not isinstance(other, PotentialExpression):
            return False
        if not self.is_parametric:
            return (
                self.expression == other.expression
                and self.independent_variables == other.independent_variables
            )
        else:
            return (
                self.expression == other.expression
                and self.independent_variables == other.independent_variables
                and _are_equal_parameters(self.parameters, other.parameters)
            )

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

    def clone(self, fast_copy=False):
        """Return a clone of this potential expression, faster alternative to deepcopying."""
        if not fast_copy:
            return PotentialExpression(
                deepcopy(self._expression),
                deepcopy(self._independent_variables),
                (
                    {
                        k: (
                            u.unyt_quantity(v.value, v.units)
                            if v.value.shape == ()
                            else u.unyt_array(v.value, v.units)
                        )
                        for k, v in self._parameters.items()
                    }
                    if self._is_parametric
                    else None
                ),
                verify_validity=False,
            )
        elif fast_copy:
            return PotentialExpression(
                self._expression,
                self._independent_variables,
                (
                    {
                        k: (
                            u.unyt_quantity(v.value, v.units)
                            if v.value.shape == ()
                            else u.unyt_array(v.value, v.units)
                        )
                        for k, v in self._parameters.items()
                    }
                    if self._is_parametric
                    else None
                ),
                verify_validity=False,
            )

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
    @lru_cache(maxsize=128)
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
            parameter_symbols = sympy.symbols(parameters)
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
                symbols = sympy.symbols(parameters)
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

    @classmethod
    def from_non_parametric(
        cls,
        non_parametric: "PotentialExpression",
        parameters: Dict[str, u.unyt_array],
        valid: bool = False,
    ) -> "PotentialExpression":
        """Create a parametric expression from a non-parametric one.

        Parameters
        ----------
        non_parametric: PotentialExpression
            The non-parametric potential expression to create the parametric one from

        parameters: dict
            The dictionary of parameters for the newly created parametric expression.

        valid: bool, default=False
            Whether to validate expression/independent_variables and with the parameters.

        Notes
        -----
        When `valid=True`, the validation checks (on whether or not the expression/independent_variables
        match with the provided parameters is not performed. Use with caution.

        Returns
        -------
        PotentialExpression
            The parametric potential expression from the provided parameters
        """
        if not isinstance(non_parametric, cls):
            raise TypeError(
                f"Expected {non_parametric} to be of type {cls} but found "
                f"{type(non_parametric)}."
            )

        if non_parametric.is_parametric:
            raise ValueError(
                "Cannot create a parametric expression from a parametric "
                "expression."
            )

        else:
            return cls(
                expression=deepcopy(non_parametric.expression),
                parameters=parameters,
                independent_variables=deepcopy(
                    non_parametric.independent_variables
                ),
                verify_validity=not valid,
            )
