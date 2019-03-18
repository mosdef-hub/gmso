import warnings
import unyt as u

from topology.testing.utils import allclose
from topology.core.potential import Potential


class AtomType(Potential):
    """An atom type, inheriting from the Potential class.

    AtomType represents an atom type and includes the functional form describing its interactions and,
    optionally, other properties such as mass and charge.
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
    mass : unyt.unyt_quantity, optional, default=0.0 * unyt.elementary_charge
        The mass of the atom type.
    charge : unyt.unyt_quantity, optional, default=0.0 * unyt.g / u.mol
        The charge of the atom type.
    expression : str or sympy.Expr,
                 default='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
        The mathematical expression describing the functional form of the
        potential describing this atom type, i.e. a Lennard-Jones potential.
    parameters : dict of str : unyt.unyt_quantity pairs,
        default={'sigma': 0.3 * u.nm, 'epsilon': 0.3 * u.Unit('kJ')},
        The parameters of the potential describing this atom type and their
        values, as unyt quantities.
    independent_variables : str or sympy.Symbol or a list or set thereof
        The independent variables in the expression of the potential
        describing this atom type.

    """

    def __init__(self,
                 name='AtomType',
                 mass=0.0 * u.gram / u.mol,
                 charge=0.0 * u.elementary_charge,
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 parameters={
                    'sigma': 0.3 * u.nm,
                    'epsilon': 0.3 * u.Unit('kJ')},
                 independent_variables={'r'}):

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
