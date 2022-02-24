"""A topology within a topology."""
import warnings
from copy import deepcopy

from boltons.setutils import IndexedSet

from gmso.core.angle import Angle
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.dihedral import Dihedral
from gmso.core.improper import Improper
from gmso.core.topology import Topology


class SubTopology(object):
    """A sub-topology i.e. topology within a topology.

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
        """Return the name of the sub-topology."""
        return self._name

    @name.setter
    def name(self, name):
        """Set the name of the sub-topology."""
        self._name = str(name)

    @property
    def sites(self):
        """Return the sites associated with the sub-topology."""
        return self._sites

    @property
    def n_sites(self):
        """Return the number of sites associated with the sub-topology."""
        return len(self.sites)

    @property
    def parent(self):
        """Return the parent of the sub-topology."""
        return self._parent

    @parent.setter
    def parent(self, parent):
        """Set the parent of the sub-topology."""
        warnings.warn(
            "Setting a parent is potentially dangerous. Consider using "
            "Topology.add_subtopology instead"
        )
        if parent is None:
            raise NotImplementedError(
                "Setting parents to None is not yet supported"
            )
        self._parent = _validate_parent(parent)

    def add_site(self, site, update_types=True):
        """Add a site to this sub-topology.

        This method adds a site to the sub-topology.
        If the sub-topology has a parent, the site will
        also be added to the parent topology. If the
        update_types parameter is set to true (default
        behavior), this method will also check if there
        is an gmso.AtomType associated with the site and
        it to the sub-topology's AtomTypes collection.

        Parameters
        ----------
        site : gmso.Atom
            The site to be added to this sub-topology
        update_types : (bool), default=True
            If true, add this site's atom type to the sub-topology's set of AtomTypes

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
            self.parent.add_site(site, update_types=update_types)

    def __repr__(self):
        """Return a formatted representation of the sub-topology."""
        return (
            f"<SubTopology {self.name},\n "
            f"{self.n_sites} sites,\n "
            f"id: {id(self)}>"
        )

    def __str__(self):
        """Return a string representation of the sub-topology."""
        return (
            f"<SubTopology {self.name}, "
            f"{self.n_sites} sites, "
            f"id: {id(self)}>"
        )

    def json_dict(self):
        """Return a json serializable dictionary of this subtopology."""
        subtop_dict = {"name": self.name, "atoms": []}

        for site in self._sites:
            subtop_dict["atoms"].append(self.parent.get_index(site))

        return subtop_dict

    def to_top(self):
        """Return a Topology formed by sites of this Sub-Topology."""

        top = Topology(name=self.name)
        connection_type_attr_map = {
            Bond: "bond_type",
            Angle: "angle_type",
            Dihedral: "dihedral_type",
            Improper: "improper_type",
        }
        copy_excludes = {"topology_", "set_ref_"}
        sites_map = {}
        for site in self._sites:
            site_copy = site.copy(deep=True, exclude={"atom_type_"})
            if site.atom_type:
                atom_type_copy = site.atom_type.copy(
                    deep=True, exclude=copy_excludes
                )
                site_copy.atom_type = atom_type_copy
            sites_map[id(site)] = site_copy
            top.add_site(site, update_types=False)

        # create new connections
        for connection in self._parent._connections:
            new_members = []
            for member in connection.connection_members:
                if id(member) in sites_map:
                    new_members.append(sites_map[id(member)])

            if len(new_members) == 0:
                break

            if len(new_members) < len(connection.connection_members):
                raise Exception(  # ToDo: raise a better error
                    "One or more sites in this Topology are connected to other sites"
                    "which are not part of this subTopology."
                )
            ConnClass = type(connection)
            connection_copy = ConnClass(connection_members=new_members)
            conn_type = getattr(connection, connection_type_attr_map[ConnClass])
            if conn_type:
                setattr(
                    connection_copy,
                    connection_type_attr_map[ConnClass],
                    conn_type.copy(deep=True, exclude=copy_excludes),
                )

            top.add_connection(connection_copy, update_types=False)

        top.combining_rule = self._parent.combining_rule
        top.scaling_factors = deepcopy(self._parent.scaling_factors)
        top.update_topology()

        atom_type_names = set(atom_type.name for atom_type in top.atom_types)

        # create new pairpotential types
        for ptype in self._parent.pairpotential_types:
            at1_name = ptype.member_types[0]
            at2_name = ptype.member_types[1]
            if at1_name in atom_type_names and at2_name in atom_type_names:
                top.add_pairpotentialtype(
                    ptype.copy(deep=True, exclude=copy_excludes)
                )

        return top

    @classmethod
    def from_sites(cls, sites, parent, update_types=False):
        subtop = cls(parent=parent)
        for site in sites:
            subtop.add_site(site, update_types=update_types)
        return subtop


def _validate_parent(parent):
    """Ensure the parent is a topology."""
    if isinstance(parent, Topology):
        return parent
    else:
        raise TypeError("Argument {} is not type Topology".format(parent))


def _validate_site_addability(site):
    """Ensure a site is a site and not already a part of a top/subtop."""
    if not isinstance(site, Atom):
        raise TypeError("Argument {} is not a Site. See gmso/core/atom.py")
    # TODO: Some sort of a check on site.parent
    return site
