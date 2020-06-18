from abc import abstractmethod
from typing import Union, Iterable, TypeVar

import sympy

from pydantic import Field, validator, root_validator

from gmso.abc.gmso_base import GMSOBase

IndependentVarType = Union[Union[str, sympy.core.Symbol],
                           Iterable[Union[str, sympy.core.Symbol]]]
PotentialT = TypeVar('PotentialT', bound='Potential')


class Potential(GMSOBase):
    __base_doc__ = """An abstract potential class.

    Potential stores a general interaction between components of a chemical
    topology that can be specified by a mathematical expression. The functional
    form of the potential is stored as a `sympy` expression and the parameters
    are stored explicitly. This class is agnostic to the instantiation of the
    potential, which can be e.g. a non-bonded potential, a bonded potential, an
    angle potential, a dihedral potential, etc. and is designed to be inherited
    by classes that represent these potentials.
    """
    name_: str = Field(
        '',
        description='The name of the potential. Defaults to class name'
    )

    expression_: Union[str, sympy.Expr] = Field(
        'a*x+b',
        description='The mathematical expression describing the functional form of the potential.'
    )

    independent_variables_: IndependentVarType = Field(
        {sympy.Symbol('x')},
        description='The independent variables in the expression of the potential.'
    )

    @property
    def name(self):
        return self.__dict__.get('name_')

    @property
    def expression(self):
        return self.__dict__.get('expression_')

    @property
    def independent_variables(self):
        return self.__dict__.get('independent_variables_')

    @validator('expression_')
    def _validate_expression(cls, expression):
        if isinstance(expression, sympy.Expr):
            pass
        elif isinstance(expression, str):
            return sympy.sympify(expression)
        else:
            raise ValueError(
                'Please enter a string or sympy expression'
            )
        return expression

    @abstractmethod
    def set_expression(self, **kwargs):
        raise NotImplementedError

    @validator('independent_variables_')
    def _validate_independent_variables(cls, indep_vars):
        if isinstance(indep_vars, str):
            indep_vars = {sympy.symbols(indep_vars)}
        elif isinstance(indep_vars, sympy.Symbol):
            indep_vars = {indep_vars}
        elif isinstance(indep_vars, Iterable):
            if all([isinstance(val, sympy.Symbol) for val in indep_vars]):
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

    @root_validator(pre=True)
    def inject_name(cls, values):
        print(values)
        if not values.get('name'):
            values['name'] = cls.__name__
        return values

    class Config:
        fields = {
            'name_': 'name',
            'expression_': 'expression',
            'independent_variables_': 'independent_variables'
        }

        alias_to_fields = {
            'name': 'name_',
            'expression': 'expression_',
            'independent_variables': 'independent_variables_'
        }

        validate_assignment = True
