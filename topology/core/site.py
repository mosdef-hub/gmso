import warnings

import numpy as np
import unyt as u


class Site(object):
    """A general site."""

    def __init__(self,
                 name,
                 position=None,
                 charge=None,
                 element=None,
                 atom_type=None):
        self.name = str(name)
        if position is None:
            self.position = u.nm * np.zeros(3)
        else:
            self.position = _validate_position(position)
        if element:
            self.element = element

        self.atom_type = atom_type
        self._charge = charge
        self._connections = list()

    def add_connection(self, other_site):
        self._connections.append(other_site)

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
        else:
            return self.atom_type.charge

    @charge.setter
    def charge(self, charge):
        self._charge = float(charge)


def _validate_position(position):
    if not isinstance(position, u.unyt_array):
        warnings.warn('Positions are assumed to be in nm')
        position *= u.nm

    input_unit = position.units

    position = np.asarray(position, dtype=float, order='C')
    np.reshape(position, newshape=(3, ), order='C')

    position *= input_unit
    position.convert_to_units(u.nm)

    return position
