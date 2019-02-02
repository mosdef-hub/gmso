import numpy as np


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
        self.site_list = list()

    def add_site(self, site):
        self.site_list.append(site)

    @property
    def n_sites(self):
        return len(self.site_list)

    def positions(self):
        xyz = np.empty(shape=(self.n_sites, 3))
        for i, site in enumerate(self.site_list):
            xyz[i, :] = site.position
        return xyz
