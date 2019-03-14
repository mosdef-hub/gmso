import numpy as np
import unyt as u

from topology.core.connection import Connection


class Topology(object):
    """A topology.

    Parameters
    ----------
    name : str, optional
        A name for the Topology.
    """
    def __init__(self, name="Topology", box=None):
        if name is not None:
            self._name = name
        self._box = box
        self._site_list = list()
        self._connection_list = list()

        self._atom_types = list()
        self._connection_types = list()

        self._atom_type_functionals = list()
        self._connection_type_functionals = list()

    def add_site(self, site):
        self._site_list.append(site)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = str(name)

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, box):
        self._box = box

    @property
    def n_sites(self):
        return len(self._site_list)

    @property
    def atom_types(self):
        return self._atom_types

    @property
    def connection_types(self):
        return self._connection_types

    @property
    def atom_type_functionals(self):
        return [atype.nb_function for atype in self.atom_types]

    @property
    def connection_type_functionals(self):
        return [ctype.potential_function for ctype in self.connection_types]

    def positions(self):
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self.site_list):
            xyz[i, :] = site.position
        return xyz

    @property
    def site_list(self):
        return self._site_list

    @property
    def connection_list(self):
        return self._connection_list

    @property
    def n_connections(self):
        return len(self._connection_list)

    def add_connection(self, connection):
        self._connection_list.append(connection)

    def update_top(self):
        """ Update the entire topology's unique types

        Notes
        -----
        Will update: atom_types, connnection_types
        """
        self.update_atom_types()
        self.update_connection_list()
        self.update_connection_types()
        
    def update_atom_types(self):
        self._atom_types = list(set(
            site.atom_type for site in self.site_list))

    def update_connection_list(self):
        for site in self.site_list:
            for neighbor in site.connections:
                temp_connection = Connection(site, neighbor, update=False)
                if temp_connection not in self.connection_list:
                    self.add_connection(Connection(site, neighbor, update=True))

    def update_connection_types(self):
        self._connection_types = list(set(
            connection.connection_type for connection in self.connection_list))

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('{:d} connections, '.format(self.n_connections))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)
