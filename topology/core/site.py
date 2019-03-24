import warnings

import numpy as np
import unyt as u

from topology.core.atom_type import AtomType
from topology.testing.utils import allclose


class Site(object):
    """A general site."""

    def __init__(self,
                 name='Site',
                 position=None,
                 charge=None,
                 mass=None,
                 element=None,
                 atom_type=None):
        if name is not None:
            self.name = str(name)
        if position is None:
            self.position = u.nm * np.zeros(3)
        else:
            self.position = _validate_position(position)

        self._element = element
        self._atom_type = _validate_atom_type(atom_type)
        self._charge = _validate_charge(charge)
        self._mass = _validate_mass(mass)
        self._connections = list()

    def add_connection(self, connection):
        connection = _validate_connection(connection)
        self._connections.append(connection)

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, element):
        self._element = element

    @property
    def connections(self):
        return self._connections

    @property
    def n_connections(self):
        return len(self._connections)

    @property
    def charge(self):
        if self._charge is not None:
            return self._charge
        elif self.atom_type is not None:
            return self.atom_type.charge
        else:
            return None

    @charge.setter
    def charge(self, charge):
        self._charge = _validate_charge(charge)

    @property
    def mass(self):
        if self._mass is not None:
            return self._mass
        elif self.atom_type is not None:
            return self.atom_type.mass
        else:
            return None

    @mass.setter
    def mass(self, mass):
        self._mass = _validate_mass(mass)

    @property
    def atom_type(self):
        return self._atom_type

    @atom_type.setter
    def atom_type(self, val):
        val = _validate_atom_type(val)
        self._atom_type = val

    def __eq__(self, other):
        if not allclose(self.position, other.position):
            return False
        if not allclose(self.charge, other.charge, atol=1e-22):
            return False
        if self.atom_type != other.atom_type:
            return False

        return True

    def __hash__(self):
        return id(self)

    def __repr__(self):
        return "<Site {}, id {}>".format(self.name, id(self))

def _validate_position(position):
    if not isinstance(position, u.unyt_array):
        warnings.warn('Positions are assumed to be in nm')
        position *= u.nm

    input_unit = position.units

    position = np.asarray(position, dtype=float, order='C')
    position = np.reshape(position, newshape=(3, ), order='C')

    position *= input_unit
    position.convert_to_units(u.nm)

    return position

def _validate_charge(charge):
    if charge is None:
        return None
    elif not isinstance(charge, u.unyt_array):
        warnings.warn("Charges are assumed to be elementary charge")
        charge *= u.elementary_charge
    elif charge.units.dimensions != u.elementary_charge.units.dimensions:
        warnings.warn("Charges are assumed to be elementary charge")
        charge = charge.value * u.elementary_charge
    else:
        pass

    return charge

def _validate_mass(mass):
    if mass is None:
        return None
    elif not isinstance(mass, u.unyt_array):
        warnings.warn("Masses are assumed to be g/mol")
        mass *= u.gram/u.mol
    elif mass.units.dimensions != (u.gram/u.mol).units.dimensions:
        warnings.warn("Masses are assumed to be g/mol")
        mass = mass.value * u.gram/u.mol
    else:
        pass

    return mass

def _validate_atom_type(val):
    if val is None:
        return None
    elif not isinstance(val, AtomType):
        raise ValueError("Passed value {} is not an AtomType".format(val))
    else:
        return val

def _validate_connection(connection):
    from topology.core.connection import Connection
    if not isinstance(connection, Connection):
        raise ValueError("Passed value {} is not a Connection".format(connection))
    return connection
