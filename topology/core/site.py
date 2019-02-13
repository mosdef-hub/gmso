import numpy as np


class Site(object):
    """A general site."""
    def __init__(self, name, position=None, charge=None, 
            element=None, atom_type=None):
        self.name = name
        if position is None:
            self.position = np.zeros(3)
        else:
            self.position = position
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
        self._charge = charge

