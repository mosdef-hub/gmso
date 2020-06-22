import warnings

import sympy
import unyt as u

from gmso.utils.misc import unyt_to_hashable
from gmso.utils.expression import _PotentialExpression
from gmso.utils.decorators import confirm_dict_existence
from gmso.exceptions import GMSOError


class Potential(object):
    """An abstract potential class.

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
    parameters : dict {str: unyt.unyt_quantity},
            default={'a': 1.0*u.dimensionless, 'b': 1.0*u.dimensionless}
        The parameters of the potential and their values, as unyt quantities.
        The keys are names of the variables included in `expression` and values
        are the numerical values of these parameters recorded as instances of
        `unyt.unyt_quantity`, which combine both value and unit information.
    independent_variables : str or sympy.Symbol or list or set thereof
        The independent variables in the expression of the potential.
    topology: gmso.core.Topology, the topology of which this potential is a part of, default=None
    set_ref: (str), the string name of the bookkeeping set in topology class.

    """

    def __init__(self,
                 name="Potential",
                 expression='a*x+b',
                 parameters=None,
                 independent_variables=None,
                 template=False,
                 topology=None
                 ):
        self._name = name

        if expression is None:
            expression = 'a*x+b'

        if parameters is None:
            parameters = {
                'a': 1.0 * u.dimensionless,
                'b': 1.0 * u.dimensionless
            }

        if independent_variables is None:
            independent_variables = {'x'}

        if template:
            self._template = template
            parameters = None

        self._potential_expression = _PotentialExpression(
            expression=expression,
            independent_variables=independent_variables,
            parameters=parameters
        )

        if topology is not None:
            self._topology = topology
        else:
            self._topology = None

    @property
    def name(self):
        return self._name

    @name.setter
    @confirm_dict_existence
    def name(self, val):
        self._name = val

    @property
    def parameters(self):
        return self._potential_expression.parameters

    @parameters.setter
    @confirm_dict_existence
    def parameters(self, newparams):
        self._potential_expression.parameters = newparams

    @property
    def independent_variables(self):
        return self._potential_expression.independent_variables

    @independent_variables.setter
    @confirm_dict_existence
    def independent_variables(self, indep_vars):
        self._potential_expression.independent_variables = indep_vars

    @property
    def template(self):
        return self._template

    @property
    def expression(self):
        return self._potential_expression.expression

    @expression.setter
    @confirm_dict_existence
    def expression(self, expression):
        self._potential_expression.expression = expression

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, top):
        self._topology = top

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
        self._potential_expression.set(
            expression=expression,
            independent_variables=independent_variables,
            parameters=parameters
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    self._potential_expression
                )
            )
        )

    def __repr__(self):
        desc = "<Potential {}, id {}>".format(self._name, id(self))
        return desc

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
        gmso.Potential
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
                   topology=topology,
                   template=False)
