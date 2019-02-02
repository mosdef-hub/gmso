import numpy as np


class Site(object):
    """A general site."""
    def __init__(self, name, position=None, element=None, atom_type=None):
        self.name = name
        if position is None:
            self.position = np.zeros(3)
        else:
            self.position = position
        if element:
            self.element = element
        if atom_type:
            self.atom_type = atom_type
        self._connections = list()

    def add_connection(self, other_site):
        self._connections.append(other_site)

    @property
    def connections(self):
        return self._connections
