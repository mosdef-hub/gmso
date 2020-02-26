import warnings
from functools import reduce


import numpy as np
import unyt as u
from boltons.setutils import IndexedSet

from topology.core.atom_type import AtomType
from topology.exceptions import TopologyError


class Site(object):
    """An interaction site object in the topology hierarchy.

    Site is the object that represents any general interaction site in a molecular simulation.
    Sites have been designed to be as general as possible, making no assumptions about representing atoms or beads, or having mass or charge.
    That is, a Site can represent an atom in an atomistic system, a bead in a coarse-grained system, and much more.

    Parameters
    ----------
    name : str, optional, default='Site'
       Name of the site
    position : unyt array or numpy array or list, optional, default=None
       The position of the site in Cartesian space.
       If a unyt array is not passed, units are assumed to be in 'nm'.
    charge : unyt quantity or float, optional, default=None
       The charge of the site.
       Unyt quantities are converted to units of elementary charge, float values are assumed to be in units of elementary charge.  
       If no value is passed, site attempts to grab a charge from site.atom_type.
    mass : unyt quantity or float, optional, default=None
       The mass of the site.  
       Unyt quantities are converted to units of g/mol, float values are assumed to be in units of g/mol.  
       If no value is passed, site attempts to grab a mass from site.atom_type.
    element : 'Element' object, optional, default=None
       The element of the site represented by the `Element` object.  
       See `element.py` for more information.
    atom_type : 'AtomType' object, optional, default=None
       The atom type of the site containing functional forms, interaction parameters, and other properties such as mass and charge.  
       See `atom_type.py` for more information.

    Attributes
    ----------
    connections : IndexedSet
       Set that contains connection information for the site
    n_connections : int
       Number of connections for the site

    """
    def __init__(self,
                 name='Site',
                 position=None,
                 charge=None,
                 mass=None,
                 element=None,
                 atom_type=None):
        if name is not None:
            self.name = str(name)
        if name is None:  # Guardrail against checking deliberate None
            self.name = 'Site'
        if position is None:
            self.position = u.nm * np.zeros(3)
        else:
            self.position = _validate_position(position)

        self._element = element
        self._atom_type = _validate_atom_type(atom_type)
        self._charge = _validate_charge(charge)
        self._mass = _validate_mass(mass)
        self._connections = IndexedSet()

    def add_connection(self, connection):
        connection = _validate_connection(self, connection)
        if connection:
            self._connections.add(connection)

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


def _validate_connection(site, connection):
    if not isinstance(site, Site):
        raise ValueError("Passed value {} is not a site".format(site))
    from topology.core.connection import Connection
    if not isinstance(connection, Connection):
        raise ValueError("Passed value {} is not a Connection".format(connection))
    if site not in connection.connection_members:
        raise TopologyError("Error: Site not in connection members. Cannot add the connection.")
    if connection in site.connections:
        return None
    return connection
