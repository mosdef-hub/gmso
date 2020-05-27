import warnings

from boltons.setutils import IndexedSet

from gmso.core.topology import Topology
from gmso.core.atom import Atom


class SubTopology(object):
    """A sub-topology i.e. topology within a topology

    This class provides a hierarchical topological representation to
    the topology as it imperative with many chemical structures to have
    separation of layers/ boundaries. A sub-topology can be added to a
    gmso.Topology object which will be the parent of the sub-topology.

    Parameters
    ----------
    name : str, optional, default='Sub-Topology'
        Name of the sub-topology
    parent : gmso.Topology, optional, default=None
        The parent topology of this SubTopology

    Attributes
    ----------
    sites : IndexedSet of gmso.Site objects
        Collection of sites within this sub-topology
    n_sites : int
        Number of sites withing this sub-topology
    """

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
        """Add a site to this sub-topology

        This method adds a site to the sub-topology.
        If the sub-topology has a parent, the site will
        also be added to the parent topology.

        Parameters
        ----------
        site : gmso.Atom
            The site to be added to this sub-topology

        Raises
        ------
        TypeError
            If the parameter site is not of type topology.Site
        """

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

        return ''.join(descr)

def _validate_parent(parent):
    if isinstance(parent, Topology):
        return parent
    else:
        raise TypeError('Argument {} is not type Topology'.format(parent))

def _validate_site_addability(site):
    """Ensure a site is a site and not already a part of a top/subtop"""
    if not isinstance(site, Atom):
        raise TypeError('Argument {} is not a Site. See gmso/core/site.py')
    # TODO: Some sort of a check on site.parent
    return site
