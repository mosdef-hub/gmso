from typing import Optional, Any

import unyt as u
from pydantic import Field, validator

from gmso.abc.abstract_potential import AbstractPotential
from gmso.utils.expression import _PotentialExpression
from gmso.utils.decorators import confirm_dict_existence
from gmso.exceptions import GMSOError


class ParametricPotential(AbstractPotential):
    __base_doc__ = """A parametric potential class.

    Potential stores a general interaction between components of a chemical
    topology that can be specified by a mathematical expression. The functional
    form of the potential is stored as a `sympy` expression and the parameters
    are stored explicitly. This class is agnostic to the instantiation of the
    potential, which can be e.g. a non-bonded potential, a bonded potential, an
    angle potential, a dihedral potential, etc. and is designed to be inherited
    by classes that represent these potentials.
    """

    # FIXME: Use proper forward referencing??
    topology_: Optional[Any] = Field(
        None,
        description="the topology of which this potential is a part of"
    )

    set_ref_: Optional[str] = Field(
        None,
        description='The string name of the bookkeeping set in gmso.Topology class. '
                    'This is used to track property based hashed object\'s '
                    'changes so that a dictionary/set can keep track of them'
    )

    def __init__(self,
                 name="ParametricPotential",
                 expression='a*x+b',
                 parameters=None,
                 independent_variables=None,
                 topology=None,
                 **kwargs
                 ):

        if expression is None:
            expression = 'a*x+b'

        if parameters is None:
            parameters = {
                'a': 1.0 * u.dimensionless,
                'b': 1.0 * u.dimensionless
            }

        if independent_variables is None:
            independent_variables = {'x'}

        _potential_expression = _PotentialExpression(
            expression=expression,
            independent_variables=independent_variables,
            parameters=parameters
        )

        super().__init__(
            name=name,
            potential_expression=_potential_expression,
            topology=topology,
            **kwargs
        )

    @property
    def parameters(self):
        """Optional[dict]\n\tThe parameters of the potential and their values, as unyt quantities"""
        return self.potential_expression_.parameters

    @property
    def topology(self):
        return self.__dict__.get('topology_')

    @property
    def set_ref(self):
        return self.__dict__.get('set_ref_')

    @validator('topology_')
    def is_valid_topology(cls, value):
        if value is None:
            return None
        else:
            from gmso.core.topology import Topology
            if not isinstance(value, Topology):
                raise TypeError(f'{type(value).__name__} is not of type Topology')
        return value

    @confirm_dict_existence
    def __setattr__(self, key: Any, value: Any) -> None:
        if key == 'parameters':
            self.potential_expression_.parameters = value
        else:
            super().__setattr__(key, value)

    @confirm_dict_existence
    def set_expression(self, expression=None, parameters=None, independent_variables=None):
        """Set the expression, parameters, and independent variables for this potential.

        Parameters
        ----------
        expression: sympy.Expression or string
            The mathematical expression corresponding to the potential
            If None, the expression remains unchanged
        parameters: dict
            {parameter: value} in the expression
            If None, the parameters remain unchanged

        Notes
        -----
        Be aware of the symbols used in the `expression` and `parameters`.
        If unnecessary parameters are supplied, an error is thrown.
        If only a subset of the parameters are supplied, they are updated
        while the non-passed parameters default to the existing values
       """
        self.potential_expression_.set(
            expression=expression,
            independent_variables=independent_variables,
            parameters=parameters
        )

    @classmethod
    def from_template(cls, potential_template, parameters, topology=None):
        """Create a potential object from the potential_template

        Parameters
        ----------
        potential_template : gmso.lib.potential_templates.PotentialTemplate,
                            The potential template object
        parameters : dict,
                    The parameters of the potential object to create
        topology : gmso.Topology, default=None
                   The topology to which the created potential object belongs to

        Returns
        -------
        gmso.ParametricPotential
            The potential object created

        Raises
        ------
        GMSOError
            If potential_template is not of instance PotentialTemplate
        """
        from gmso.lib.potential_templates import PotentialTemplate
        if not isinstance(potential_template, PotentialTemplate):
            raise GMSOError(f'Object {type(potential_template)} is not an instance of PotentialTemplate.')

        return cls(name=potential_template.name,
                   expression=potential_template.expression,
                   independent_variables=potential_template.independent_variables,
                   parameters=parameters,
                   topology=topology)

    class Config:
        fields = {
            'topology_': 'topology',
            'set_ref_': 'set_ref'
        }
        alias_to_fields = {
            'topology': 'topology_'
        }
        validate_assignment = True
