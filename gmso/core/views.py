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


def filter_unique_potentials(potential_types):
    visited = defaultdict(set)
    for potential_type in potential_types:
        if get_identifier(potential_type) not in visited[type(potential_type)]:
            visited[type(potential_type)].add(
                potential_type.atom_types or potential_type.atom_classes
            )

            yield potential_type


class TopologyPotentialView:
    """A potential view based on different filters for a topology potentials.

    iterator: typing.Iterator, required=True
        An iterator of either topology sites or connections from which to extract potentials from

    attribute: str, required=True
        The attribute information to extract i.e atom_type, bond_type, angle_type etc...
    """

    def __init__(self, iterator):
        self.iterator = iterator

    def __iter__(self):
        yield from self.yield_view(False, None)

    def index(self, item):
        for j, potential in enumerate(self.yield_view(False)):
            if potential is item:
                return j

    def yield_view(self, unique=False, checker=None):
        if not unique:
            for item in self.iterator:
                yield getattr(
                    item, potential_attribute_map[type(item)]
                )  # Since this use is internal, KeyErrors N/A
        else:
            if checker is None:
                checker = filter_unique_potentials

            collected_potentials = (
                getattr(item, potential_attribute_map[type(item)])
                for item in self.iterator
            )

            for item in checker(collected_potentials):
                yield item

    def __call__(self, unique=False, checker=None):
        yield from self.yield_view(unique, checker)

    def __len__(self):
        return len(self.iterator)
