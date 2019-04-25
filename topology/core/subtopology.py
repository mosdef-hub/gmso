import warnings

from boltons.setutils import IndexedSet

from topology.core.topology import Topology


class SubTopology(object):
    """A sub-topology."""
    def __init__(self, name="Topology", parent=None):
        if name is not None:
            self._name = name
        if parent is not None:
            self._parent = _validate_parent(parent)
        self._sites = IndexedSet()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = str(name)

    @property
    def sites(self):
        return self._sites

    @property
    def n_sites(self):
        return len(self.sites)

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        self._parent = _validate_parent(parent)

    def add_site(self, site):
        if site in self.sites:
            warnings.warn("Redundantly adding Site {}".format(site))
        self._sites.add(site)

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')
        descr.append('{:d} sites, '.format(self.n_sites))
        descr.append('id: {}>'.format(id(self)))


def _validate_parent(parent):
    if isinstance(parent, Topology):
        return parent
    else:
        raise TypeError('Argument {} is not type Topology'.format(parent))
