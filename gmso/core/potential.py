import warnings

import sympy
import unyt as u

from gmso.abc.abstract_potential import AbstractPotential
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.decorators import confirm_dict_existence
from gmso.exceptions import GMSOError


class ParametricPotential(AbstractPotential):
    """A parametric potential class

    In addition to the attributes from the parent class, a ParametricPotential class
    has parameters, which are values for the dependent variables in its expression.
    Generally, a parametric potential is stored in topology which is the container for
    the potential.

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
    """

    __slots__ = (
        '_parameters',
        '_topology',
        '_dict_ref',
    )

    def __init__(self,
                 name="Potential",
                 expression=None,
                 parameters=None,
                 independent_variables=None,
                 topology=None,
                 set_ref=None
                 ):
        super(ParametricPotential, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables
        )

        if parameters is None:
            parameters = {'a': 1.0 * u.dimensionless,
                          'b': 1.0 * u.dimensionless}

        self._parameters = _validate_parameters(parameters)
        self._validate_expression_parameters()

        if topology is not None:
            self._topology = topology
            self._dict_ref = set_ref
        else:
            self._topology = None
            self._dict_ref = None

    @AbstractPotential.name.setter
    @confirm_dict_existence
    def name(self, val):
        self._name = val

    @AbstractPotential.independent_variables.setter
    @confirm_dict_existence
    def independent_variables(self, indep_vars):
        self._independent_variables = self._validate_independent_variables(indep_vars)

    @AbstractPotential.expression.setter
    @confirm_dict_existence
    def expression(self, expression):
        self._expression = self._validate_expression(expression)
        self._validate_expression_parameters()

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    @confirm_dict_existence
    def parameters(self, newparams):
        newparams = _validate_parameters(newparams)

        self._parameters.update(newparams)
        self._validate_expression_parameters()

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, top):
        self._topology = top

    @confirm_dict_existence
    def set_expression(self,
                       expression=None,
                       parameters=None,
                       independent_variables=None):
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
        if expression is not None:
            self._expression = self._validate_expression(expression)

        if parameters is None:
            parameters = self._parameters
        else:
            parameters = _validate_parameters(parameters)
            if not set(self._parameters).intersection(set(parameters)):
                if expression is None:
                    raise ValueError('`parameters` argument includes no '
                                     'variables found in expression. Expected '
                                     'at least one of {}'.format(
                                        self._parameters.keys()))
            self._parameters.update(parameters)

        if independent_variables is not None:
            self._independent_variables = self._validate_independent_variables(independent_variables)

        if not set(parameters.keys()).isdisjoint(self._expression.free_symbols):
            raise ValueError('Mismatch between parameters and expression symbols')

        self._validate_expression_parameters()

    def _validate_expression_parameters(self):
        # Check for unused symbols
        parameter_symbols = sympy.symbols(set(self._parameters.keys()))
        independent_variable_symbols = self._independent_variables
        used_symbols = parameter_symbols.union(independent_variable_symbols)
        unused_symbols = self._expression.free_symbols - used_symbols
        if len(unused_symbols) > 0:
            warnings.warn('You supplied parameters with '
                          'unused symbols {}'.format(unused_symbols))

        if used_symbols != self._expression.free_symbols:
            symbols = sympy.symbols(set(self.parameters.keys()))
            if symbols != self._expression.free_symbols:
                missing_syms = self._expression.free_symbols - symbols - self._independent_variables
                if missing_syms:
                    raise ValueError("Missing necessary parameters to evaluate "
                                     "potential expression. Missing symbols: {}"
                                     "".format(missing_syms))
                extra_syms = symbols ^ self.expression.free_symbols
                warnings.warn("Potential expression and parameter"
                              " symbols do not agree,"
                              " extraneous symbols:"
                              " {}".format(extra_syms))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    self.expression,
                    tuple(self.independent_variables),
                    tuple(self.parameters.keys()),
                    tuple(unyt_to_hashable(val) for val in self.parameters.values())
                )
            )
        )

    def __repr__(self):
        desc = "<ParametricPotential {}, id {}>".format(self._name, id(self))
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


