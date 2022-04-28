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


default_filters = {PotentialFilters.UNIQUE_NAME_CLASS: unique_potentials}


class TopologyPotentialView:
    """A potential view based on different filters for a topology's potentials.

    Parameters
    ----------
    iterator: typing.Iterator, required=True
        An iterator of either topology sites or connections from which to extract potentials from

    Notes
    -----
    The implementation of this class has been inspired by networkx.classes.reportviews.NodeView by extending
    the idea of a view with filteration capabilities. See the source for NodeView for further details

    https://github.com/networkx/networkx/blob/12c1a00cd116701a763f7c57c230b8739d2ed085/networkx/classes/reportviews.py#L115-L279
    """

    def __init__(self, iterator):
        self.iterator = iterator

    def __iter__(self):
        yield from self.yield_view()

    def index(self, item):
        for j, potential in enumerate(self.yield_view()):
            if potential is item:
                return j

    def yield_view(self, filter_by=None, repeat=False):
        """Yield a view of the potentials of the iterator provided.

        Parameters
        ----------
        filter_by: str or func, default=None
            If provided, filter the collected potentials by some function
            see, default_filters for names of the default potential filters

        repeat: bool, default=False
            If true, repeat same objects in the iterator

        Yields
        ------
        gmso.core.ParametricPotential
            An instance of gmso.core.ParametricPotential from the attributes of Sites/Connections
            in the iterator.

        Notes
        -----
        It is allowed that multiple sites/connections point to the same parametric potential
        by design. By default this function will just return a single object even if they are
        referenced by multiple sites or connections, if so isn't desired, set the argument
        repeat=True.
        """
        visited = set()
        if not filter_by:
            for item in self.iterator:
                potential = getattr(
                    item, potential_attribute_map[type(item)]
                )  # Since this use is internal, KeyErrors N/A
                if repeat and potential:
                    yield potential
                else:
                    if potential and id(potential) not in visited:
                        visited.add(id(potential))
                        yield potential
        else:
            collected_potentials = filter(
                lambda p: p is not None,
                (
                    getattr(item, potential_attribute_map[type(item)])
                    for item in self.iterator
                ),
            )

            for item in filter_by(collected_potentials):
                yield item

    def __call__(self, filter_by=None):
        """The call method, for turning property decorators into functions"""
        if filter_by in default_filters:
            filter_by = default_filters[filter_by]

        yield from self.yield_view(filter_by)

    def __len__(self):
        return len(
            list(self.yield_view())
        )  # This will be costly? But How frequent?
