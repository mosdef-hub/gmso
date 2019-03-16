import warnings
import numpy as np
import sympy
import unyt as u

from topology.testing.utils import allclose
from topology.core.potential import Potential


class AtomType(Potential):
    """An atom type."""

    def __init__(self,
                 name='AtomType',
                 mass=0.0 * u.gram / u.mol,
                 charge=0.0 * u.elementary_charge,
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 parameters={
                    'sigma': 0.3 * u.nm,
                    'epsilon': 0.3 * u.Unit('kJ')},
                 independent_variables = {'r'}):

        super(AtomType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables)
        self._mass = _validate_mass(mass)
        self._charge = _validate_charge(charge)

        self._validate_expression_parameters()

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, val):
        self._charge = _validate_charge(val)

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, val):
        self._mass = _validate_mass(val)



    def _validate_expression_parameters(self):
        # Check for unused symbols
        parameter_symbols = sympy.symbols(set(self._parameters.keys()))
        independent_variable_symbols = self._independent_variables
        used_symbols = parameter_symbols.union(independent_variable_symbols)
        unused_symbols = self.expression.free_symbols - used_symbols
        if len(unused_symbols) > 0:
            warnings.warn('You supplied parameters with '
                          'unused symbols {}'.format(unused_symbols))

        if used_symbols != self.expression.free_symbols:
            symbols = sympy.symbols(set(self.parameters.keys()))
            if symbols != self.expression.free_symbols:
                missing_syms = self.expression.free_symbols - symbols - self._independent_variables
                if missing_syms:
                    raise ValueError("Missing necessary parameters to evaluate "
                                     "NB function. Missing symbols: {}"
                                     "".format(missing_syms))
                extra_syms = symbols ^ self.expression.free_symbols
                warnings.warn("NB function and parameter"
                              " symbols do not agree,"
                              " extraneous symbols:"
                              " {}".format(extra_syms))

    def __eq__(self, other):
        name_match = (self.name == other.name)
        charge_match = allclose(
            self.charge,
            other.charge,
            atol=1e-6 * u.elementary_charge,
            rtol=1e-5 * u.elementary_charge)
        mass_match = allclose(
            self.mass,
            other.mass,
            atol=1e-6 * u.gram / u.mol,
            rtol=1e-5 * u.gram / u.mol)
        parameter_match = (self.parameters == other.parameters)
        expression_match = (self.expression == other.expression)

        return all([
            name_match, charge_match, mass_match, parameter_match,
            expression_match
        ])

    def __repr__(self):
        desc = "<AtomType {}, id {}>".format(self._name, id(self))
        return desc


def _validate_charge(charge):
    if not isinstance(charge, u.unyt_array):
        warnings.warn("Charges are assumed to be elementary charge")
        charge *= u.elementary_charge
    elif charge.units.dimensions != u.elementary_charge.units.dimensions:
        warnings.warn("Charges are assumed to be elementary charge")
        charge = charge.value * u.elementary_charge
    else:
        pass

    return charge


def _validate_mass(mass):
    if not isinstance(mass, u.unyt_array):
        warnings.warn("Masses are assumed to be g/mol")
        mass *= u.gram / u.mol
    elif mass.units.dimensions != (u.gram / u.mol).units.dimensions:
        warnings.warn("Masses are assumed to be g/mol")
        mass = mass.value * u.gram / u.mol
    else:
        pass

    return mass
