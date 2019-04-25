import warnings

from boltons.setutils import IndexedSet


SubTopology(object):
    """A sub-topology."""
    def __init__(self, name="Topology"):
        if name is not None:
            self._name = name
        self._sites = IndexedSet()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = str(name)

    def add_site(self, site):
        if site in self.sites:
            warnings.warn("Redundantly adding Site {}".format(site))
        self._sites.add(site)
