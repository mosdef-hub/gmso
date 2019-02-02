import numpy as np

from topology.core.connection import Connection


class Topology(object):
    """A topology.

    Parameters
    ----------
    name : str, optional
        A name for the Topology.
    """
    def __init__(self, name=None, box=None):
        if name:
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
        xyz = np.empty(shape=(self.n_sites, 3))
        for i, site in enumerate(self.site_list):
            xyz[i, :] = site.position
        return xyz

    @property
    def site_list(self):
        return self._site_list

    @property
    def connection_list(self):
        return self._connection_list

    def add_connection(self, connection):
        self._connection_list.append(connection)

    def check_connection_list(self):
        for site in self._site_list:
            for neighbor in site.connections:
                if len(self._connection_list) > 0:
                    for connect in self._connection_list:
                        if (site, neighbor) == (connect.site1, connect.site2):
                            continue
                        if (site, neighbor) == (connect.site2, connect.site1):
                            continue
                        new_connect = Connection(site1=site, site2=neighbor)
                        self.add_connection(new_connect)
                else:
                    new_connect = Connection(site1=site, site2=neighbor)
                    self.add_connection(new_connect)
