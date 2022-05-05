from collections import defaultdict

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType

potential_attribute_map = {
    Atom: "atom_type",
    Bond: "bond_type",
    Angle: "angle_type",
    Dihedral: "dihedral_type",
    Improper: "improper_type",
}


def get_identifier(potential):
    """Get identifier for a topology potential."""
    if isinstance(potential, AtomType):
        return potential.name
    if isinstance(potential, (BondType, AngleType, DihedralType, ImproperType)):
        return potential.member_types or potential.member_classes


def unique_potentials(potential_types):
    """Filter unique potentials based on pre-defiend identifiers."""
    visited = defaultdict(set)
    for potential_type in potential_types:
        identifier = get_identifier(potential_type)
        if identifier not in visited[type(potential_type)]:
            visited[type(potential_type)].add(identifier)

            yield potential_type


class PotentialFilters:
    UNIQUE_NAME_CLASS = "unique_name_class"


potential_filters = {PotentialFilters.UNIQUE_NAME_CLASS: unique_potentials}


class TopologyPotentialView:
    """A potential view based on different filters for a topology's potentials.

    Parameters
    ----------
    iterator: typing.Iterator, required=True
        An iterator of either topology sites or connections from which to extract potentials from

    Parameters
    ----------
    filter_by: str or function, default=None
        If provided, filter the collected potentials by some function
        see, default_filters for names of the default potential filters

    Examples
    --------
    To use a TopologyPotentialView, a Topology must have few sites with AtomTypes or BondTypes.

    >>> from gmso.core.topology import Topology
    >>> from gmso.core.atom import Atom
    >>> from gmso.core.atom_type import AtomType
    >>> top = Topology(name="ViewTopology")
    >>> sites = [Atom(name=f"Atom_{j}") for j in range(10)]
    >>> atom_type1 = AtomType(name='atom_type1')
    >>> atom_type2 = AtomType(name='atom_type2')
    >>> for site in sites:
    ...     site.atom_type = atom_type1 if int(site.name[-1]) % 2 == 0 else atom_type2
    ...     top.add_site(site)
    >>> top.update_topology()
    >>> for atom_type in top.atom_types:
    ...     print(atom_type.name)
    atom_type1
    atom_type2
    >>> top.get_index(atom_type2)
    1

    Notes
    -----
    The implementation of this class is inspired from networkx.classes.reportviews.NodeView by extending
    the idea of a view with filteration capabilities. See the source for NodeView for further details

    https://github.com/networkx/networkx/blob/12c1a00cd116701a763f7c57c230b8739d2ed085/networkx/classes/reportviews.py#L115-L279
    """

    attribute = None

    def __init__(self, iterator, filter_by=None):
        self.iterator = iterator
        self.filter_by = filter_by

    def __iter__(self):
        yield from self.yield_view()

    def index(self, item):
        for j, potential in enumerate(self.yield_view()):
            if potential is item:
                return j

    def _collect_potentials(self):
        """Collect potentials from the iterator"""
        visited = set()
        for item in self.iterator:
            potential = getattr(
                item, potential_attribute_map[type(item)]
            )  # Since this use is internal, KeyErrors N/A
            if potential and id(potential) not in visited:
                visited.add(id(potential))
                yield potential

    def yield_view(self):
        """Yield a view of the potentials of the iterator provided.

        Yields
        ------
        gmso.core.ParametricPotential
            An instance of gmso.core.ParametricPotential from the attributes of Sites/Connections
            in the iterator.
        """
        if not self.filter_by:
            yield from self._collect_potentials()

        else:
            if isinstance(self.filter_by, str):
                potential_filter = potential_filters[self.filter_by]
            else:
                potential_filter = self.filter_by

            yield from potential_filter(self._collect_potentials())

    def __call__(self, filter_by=None):
        """The call method, for turning property decorators into functions"""
        if filter_by == self.filter_by:
            return self

        return TopologyPotentialView(
            iterator=self.iterator, filter_by=filter_by
        )

    def __repr__(self):
        name = self.__class__.__name__
        return f"<{name}({tuple(self)})>"

    def __len__(self):
        return len(
            list(self.yield_view())
        )  # This will be costly? But How frequent?
