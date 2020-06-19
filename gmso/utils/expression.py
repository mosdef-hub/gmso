import warnings

import sympy

import unyt as u


class _PotentialExpression:
    """A general Expression class with parameters

    Parameters
    ----------
    expression: str or sympy.Expr
        The string expression

    parameters: dict
        A dictionary of parameter whose key is a string and values are parameters

    independent_variables: str or sympy.Symbol or list or set thereof
        The independent variables in the expression of the potential.
    """
    def __init__(self,
                 expression,
                 independent_variables,
                 parameters=None):
        self._expression = self._validate_expression(expression)
        self._independent_variables = self._validate_independent_variables(independent_variables)
        self._is_parametric = False
        if parameters is not None:
            self._is_parametric = True
            self._parameters = self._validate_parameters(parameters)

    @property
    def expression(self):
        return self._expression

    @property
    def parameters(self):
        if not self._is_parametric:
            raise AttributeError(
                f'{self} of type {self.__class__.__name__} has no attribute parameters'
            )
        return self._parameters

    @property
    def independent_variables(self):
        return self._independent_variables

    def set(self,
            expression=None,
            parameters=None,
            independent_variables=None
            ):
        if expression is not None:
            expression = self._validate_expression(expression)
        else:
            expression = self._expression

        if self._is_parametric:
            if parameters is not None:
                parameters = self._validate_parameters(parameters)
            else:
                parameters = self._parameters

            if independent_variables is not None:
                independent_variables = self._validate_independent_variables(independent_variables)
            else:
                independent_variables = self.independent_variables

        self._verify_validity(expression,
                              independent_variables,
                              parameters)

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.expression,
                    tuple(self.independent_variables),
                    tuple(self.parameters.keys())
                )
            )
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __repr__(self):
        descr = list('<PotentialExpression')
        descr.append(' expression, '.format(self.expression))
        descr.append('{:d} independent variables, '.format(self.independent_variables))

        return ''.join(descr)

    @staticmethod
    def _validate_expression(expression):
        """Check to see that an expression is a valid sympy expression"""
        if expression is None or isinstance(expression, sympy.Expr):
            pass
        elif isinstance(expression, str):
            expression = sympy.sympify(expression)
        else:
            raise ValueError(
                'Please enter a string, sympy expression or None for expression'
            )

        return expression

    @staticmethod
    def _validate_parameters(parameters):
        """Check to see that parameters is a valid dictionary with units"""
        if not isinstance(parameters, dict):
            raise ValueError("Please enter dictionary for parameters")
        for key, val in parameters.items():
            if isinstance(val, list):
                for params in val:
                    if not isinstance(params, u.unyt_array):
                        raise ValueError('Parameter value {} lacks a unyt'.format(val))
            else:
                if not isinstance(val, u.unyt_array):
                    raise ValueError('Parameter value {} lacks a unyt'.format(val))
            if not isinstance(key, str):
                raise ValueError('Parameter key {} is not a str'.format(key))

        return parameters

    @staticmethod
    def _validate_independent_variables(indep_vars):
        """Check to see that independent_variables is a set of valid sympy symbols"""
        if isinstance(indep_vars, str):
            indep_vars = {sympy.symbols(indep_vars)}
        elif isinstance(indep_vars, sympy.symbol.Symbol):
            indep_vars = {indep_vars}
        elif isinstance(indep_vars, (list, set)):
            if all([isinstance(val, sympy.symbol.Symbol) for val in indep_vars]):
                pass
            elif all([isinstance(val, str) for val in indep_vars]):
                indep_vars = set([sympy.symbols(val) for val in indep_vars])
            else:
                raise ValueError(
                    '`independent_variables` argument was a list '
                    'or set of mixed variables. Please enter a '
                    'list or set of either only strings or only '
                    'sympy symbols'
                )
        else:
            raise ValueError(
                'Please enter a string, sympy expression, '
                'or list or set thereof for independent_variables'
            )

        return indep_vars

    @staticmethod
    def _verify_validity(expression,
                         independent_variables_symbols,
                         parameters):
        parameter_symbols = sympy.symbols(parameters.keys())
        used_symbols = parameter_symbols.union(independent_variables_symbols)
        unused_symbols = expression.free_symbols - used_symbols
        if len(unused_symbols) > 0:
            warnings.warn(
                f'You supplied parameters with '
                f'unued symbols {unused_symbols}'
            )

        if used_symbols != expression.free_symbols:
            symbols = sympy.symbols(set(parameters.keys()))
            if symbols != expression.free_symbols:
                missing_syms = expression.free_symbols - symbols - independent_variables_symbols
                if missing_syms:
                    raise ValueError(
                        f'Missing necessary dependencies to evaluate '
                        f'potential expression. Missing symbols: {missing_syms}'
                    )

                extra_syms = symbols ^ expression.free_symbols
                warnings.warn(
                    f'Potential expression and parameter symbols do not agree, '
                    f'extraneous symbols: {extra_syms}'
                )

