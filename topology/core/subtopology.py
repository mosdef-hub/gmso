import warnings

from boltons.setutils import IndexedSet

from topology.core.topology import Topology
from topology.core.site import Site


class SubTopology(object):
    """A sub-topology."""
    def __init__(self, name="Sub-Topology", parent=None):
        if name is not None:
            self._name = str(name)
        if parent is None:
            self._parent = parent
        else:
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
        warnings.warn(
            'Setting a parent is potentially dangerous. Consider using '
            'Topology.add_subtopology instead'
        )
        if parent is None:
            raise NotImplementedError(
                'Setting parents to None is not yet supported'
            )
        self._parent = _validate_parent(parent)

    def add_site(self, site):
        site = _validate_site_addability(site)
        if site in self.sites:
            warnings.warn("Redundantly adding Site {}".format(site))
        self._sites.add(site)
        if self.parent:
            self.parent.add_site(site)

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

def _validate_site_addability(site):
    """Ensure a site is a site and not already a part of a top/subtop"""
    if not isinstance(site, Site):
        raise TypeError('Argument {} is not a Site. See topology/core/site.py')
    # TODO: Some sort of a check on site.parent
    return site
