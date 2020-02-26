import warnings
import unyt as u

from gmso.core.potential import Potential
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.decorators import confirm_dict_existence
from gmso.utils._constants import ATOM_TYPE_DICT


class AtomType(Potential):
    """An atom type, inheriting from the Potential class.

    AtomType represents an atom type and includes the functional form
    describing its interactions and, optionally, other properties such as mass
    and charge.  This class inhereits from Potential, which stores the
    non-bonded interaction between atoms or sites. The functional form of the
    potential is stored as a `sympy` expression and the parameters, with units,
    are stored explicitly.

    Parameters
    ----------
    name : str, default="AtomType"
        The name of the potential.
    mass : unyt.unyt_quantity, optional, default=0.0 * unyt.g / u.mol
        The mass of the atom type.
    charge : unyt.unyt_quantity, optional, default=0.0 * unyt.elementary_charge
        The charge of the atom type.
    expression : str or sympy.Expr,
                 default='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
        The mathematical expression describing the functional form of the
        potential describing this atom type, i.e. a Lennard-Jones potential.
        The default value is a 12-6 Lennard-Jones potential.
    parameters : dict of str : unyt.unyt_quantity pairs,
        default={'sigma': 0.3 * u.nm, 'epsilon': 0.3 * u.Unit('kJ')},
        The parameters of the potential describing this atom type and their
        values, as unyt quantities.
    independent_variables : str, sympy.Symbol, or list-like of str, sympy.Symbol
        The independent variables of the functional form previously described.
    atomclass : str, default=''
        The class of the atomtype
    doi : str
        Digital Object Identifier of publication where this atom type was specified
    desc : str
        Simple description of the atom type
    overrides : set of str
        Set of other atom types that this atom type overrides
    definition : str
        SMARTS string defining this atom type
    topology: gmso.core.Topology, the topology of which this atom_type is a part of, default=None
    set_ref: (str), the string name of the bookkeeping set in topology class.

    """

    def __init__(self,
                 name='AtomType',
                 mass=0.0 * u.gram / u.mol,
                 charge=0.0 * u.elementary_charge,
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 parameters=None,
                 independent_variables=None,
                 atomclass='', doi='', overrides=None, definition='',
                 description='',
                 topology=None):
        if parameters is None:
            parameters = {'sigma': 0.3 * u.nm,
                          'epsilon': 0.3 * u.Unit('kJ')}
        if independent_variables is None:
            independent_variables = {'r'}

        if overrides is None:
            overrides = set()

        super(AtomType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            topology=topology)
        self._mass = _validate_mass(mass)
        self._charge = _validate_charge(charge)
        self._atomclass = _validate_str(atomclass)
        self._doi = _validate_str(doi)
        self._overrides = _validate_set(overrides)
        self._description = _validate_str(description)
        self._definition = _validate_str(definition)
        self._set_ref = ATOM_TYPE_DICT
        self._validate_expression_parameters()

    @property
    def set_ref(self):
        return self._set_ref

    @property
    def charge(self):
        return self._charge

    @charge.setter
    @confirm_dict_existence
    def charge(self, val):
        self._charge = _validate_charge(val)

    @property
    def mass(self):
        return self._mass

    @mass.setter
    @confirm_dict_existence
    def mass(self, val):
        self._mass = _validate_mass(val)

    @property
    def atomclass(self):
        return self._atomclass

    @atomclass.setter
    @confirm_dict_existence
    def atomclass(self, val):
        self._atomclass = val

    @property
    def doi(self):
        return self._doi

    @doi.setter
    @confirm_dict_existence
    def doi(self, doi):
        self._doi = _validate_str(doi)

    @property
    def overrides(self):
        return self._overrides

    @overrides.setter
    @confirm_dict_existence
    def overrides(self, overrides):
        self._overrides = _validate_set(overrides)

    @property
    def description(self):
        return self._description

    @description.setter
    @confirm_dict_existence
    def description(self, description):
        self._description = _validate_str(description)

    @property
    def definition(self):
        return self._definition

    @definition.setter
    @confirm_dict_existence
    def definition(self, definition):
        self._definition = _validate_str(definition)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    unyt_to_hashable(self.mass),
                    unyt_to_hashable(self.charge),
                    self.expression,
                    tuple(self.independent_variables),
                    tuple(self.parameters.keys()),
                    tuple(unyt_to_hashable(val) for val in self.parameters.values())
                )
            )
        )

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


def _validate_str(val):
    if not isinstance(val, str):
        raise ValueError("Passed value {} is not a string".format(val))
    return val


def _validate_set(val):
    if not isinstance(val, set):
        raise ValueError("Passed value {} is not a set".format(val))
    if not all([isinstance(char, str) for char in val]):
        raise ValueError("Passed overrides of non-string to overrides")
    return val
