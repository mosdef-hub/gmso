from abc import ABC, abstractmethod

import sympy


class AbstractPotential(ABC):
    """An abstract potential class

    Potential stores a general interaction between components of a chemical
    topology that can be specified by a mathematical expression. The functional
    form of the potential is stored as a `sympy` expression and the parameters
    are stored explicitly. This class is agnostic to the instantiation of the
    potential, which can be e.g. a non-bonded potential, a bonded potential, an
    angle potential, a dihedral potential, etc. and is designed to be inherited
    by classes that represent these potentials.

    Parameters
    ----------
    name : str, default="Potential"
        The name of the potential.
    expression : str or sympy.Expr, default='a*x+b'
        The mathematical expression describing the functional form of the
        potential.
    independent_variables : str or sympy.Symbol or list or set thereof
        The independent variables in the expression of the potential.
    """
    __slots__ = (
        '_name',
        '_expression',
        '_independent_variables'
    )

    @abstractmethod
    def __init__(self,
                 name='Potential',
                 expression=None,
                 independent_variables=None):
        if expression is None:
            assert independent_variables is None,\
                'Cannot specify an independent variable when expression is not provided.'
            expression, independent_variables = self._default_expression()

        if independent_variables is None:
            assert expression is None,\
                'Cannot specify an expression when independent variable is not provided.'
            expression, independent_variables = self._default_expression()

        self._name = name
        self._expression = self._validate_expression(expression)
        self._independent_variables = self._validate_independent_variables(independent_variables)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def expression(self):
        return self._expression

    @expression.setter
    def expression(self, expression):
        self._expression = self._validate_expression(expression)

    @property
    def independent_variables(self):
        return self._independent_variables

    @independent_variables.setter
    def independent_variables(self, independent_variables):
        self._independent_variables = self._validate_independent_variables(independent_variables)

    @abstractmethod
    def set_expression(self, **kwargs):
        raise NotImplementedError

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
                raise ValueError('`independent_variables` argument was a list '
                                 'or set of mixed variables. Please enter a '
                                 'list or set of either only strings or only '
                                 'sympy symbols')
        else:
            raise ValueError("Please enter a string, sympy expression, "
                             "or list or set thereof for independent_variables")

        return indep_vars

    @staticmethod
    def _validate_expression(expression):
        """Check to see that an expression is a valid sympy expression"""
        if expression is None or isinstance(expression, sympy.Expr):
            pass
        elif isinstance(expression, str):
            expression = sympy.sympify(expression)
        else:
            raise ValueError("Please enter a string, sympy expression, "
                             "or None for expression")

        return expression

    @staticmethod
    def _default_expression():
        expression = sympy.sympify('a*x+b')
        independent_variables = {'x'}
        return expression, independent_variables
