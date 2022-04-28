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
    if isinstance(potential, AtomType):
        return potential.name
    if isinstance(potential, (BondType, AngleType, DihedralType, ImproperType)):
        return potential.member_types or potential.member_classes


def unique_potentials(potential_types):
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
    """A potential view based on different filters for a topology potentials.

    Parameters
    ----------
    iterator: typing.Iterator, required=True
        An iterator of either topology sites or connections from which to extract potentials from

    """

    def __init__(self, iterator):
        self.iterator = iterator

    def __iter__(self):
        yield from self.yield_view()

    def index(self, item):
        for j, potential in enumerate(self.yield_view()):
            if potential is item:
                return j

    def yield_view(self, filter_by=None):
        visited = set()
        if not filter_by:
            for item in self.iterator:
                potential = getattr(
                    item, potential_attribute_map[type(item)]
                )  # Since this use is internal, KeyErrors N/A
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

        if filter_by in default_filters:
            filter_by = default_filters[filter_by]

        yield from self.yield_view(filter_by)

    def __len__(self):
        return len(
            list(self.yield_view())
        )  # This will be costly? But How frequent?
