from typing import Any
from abc import abstractmethod

from pydantic import Field

from gmso.abc.gmso_base import GMSOBase
from gmso.utils.expression import _PotentialExpression


class AbstractPotential(GMSOBase):
    __base_doc__ = """An abstract potential class.

    AbstractPotential stores a general interaction between components of a chemical
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

    potential_expression_: _PotentialExpression = Field(
        _PotentialExpression(expression='a*x+b', independent_variables={'x'}),
        description='The mathematical expression for the potential'
    )

    def __init__(self,
                 name='Potential',
                 expression='a*x+b',
                 independent_variables=None,
                 potential_expression=None,
                 **kwargs):
        if potential_expression is None:
            if expression is None:
                expression = 'a*x+b'

            if independent_variables is None:
                independent_variables = {'x'}

            potential_expression = _PotentialExpression(
                expression=expression,
                independent_variables=independent_variables,
                parameters=None
            )

        super().__init__(
            name=name,
            potential_expression=potential_expression,
            **kwargs
        )

    @property
    def name(self):
        return self.__dict__.get('name_')

    @property
    def independent_variables(self):
        """Optional[Union[set, str]]\n\tThe independent variables in the `Potential`'s expression"""
        return self.potential_expression_.independent_variables

    @property
    def expression(self):
        """Optional[Union[str, sympy.Expr]]\n\tThe mathematical expression of the functional form of the potential"""
        return self.potential_expression_.expression

    @property
    def potential_expression(self):
        return self.__dict__.get('potential_expression_')

    @abstractmethod
    def set_expression(self):
        raise NotImplementedError

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    self.potential_expression
                )
            )
        )

    def __repr__(self):
        desc = "<{} {}, id {}>".format(self.__class__.__name__, self.name, id(self))
        return desc

    def __setattr__(self, key: Any, value: Any) -> None:
        if key == 'expression':
            self.potential_expression_.expression = value
        elif key == 'independent_variables':
            self.potential_expression_.independent_variables = value
        elif key == 'set_ref_':
            return
        else:
            super().__setattr__(key, value)

    class Config:
        fields = {
            'name_': 'name',
            'potential_expression_': 'potential_expression'
        }

        alias_to_fields = {
            'name': 'name_',
            'potential_expression': 'potential_expression_'
        }
