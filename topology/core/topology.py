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
