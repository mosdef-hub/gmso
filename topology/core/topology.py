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
            self.name = name
        if box:
            self.box = box
        self._site_list = list()
        self._connection_list = list()

    def add_site(self, site):
        self._site_list.append(site)

    @property
    def n_sites(self):
        return len(self._site_list)

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

    def update_connection_list(self):
        for site in self.site_list:
            for neighbor in site.connections:
                temp_connection = Connection(site, neighbor, update=False)
                if temp_connection not in self.connection_list:
                    self.add_connection(Connection(site, neighbor, update=True))

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('{:d} connectiosn, '.format(self.n_connections))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)
