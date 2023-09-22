import uuid
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
from gmso.utils.sorting import sort_by_types

__all__ = ["TopologyPotentialView", "PotentialFilters"]

potential_attribute_map = {
    Atom: "atom_type",
    Bond: "bond_type",
    Angle: "angle_type",
    Dihedral: "dihedral_type",
    Improper: "improper_type",
}


class MissingFilterError(KeyError):
    """Error to be raised when there's a missing builtin filter."""


def get_name_or_class(potential):
    """Get identifier for a topology potential based on name or membertype/class."""
    if isinstance(potential, AtomType):
        return potential.name
    elif isinstance(
        potential, (BondType, AngleType, DihedralType, ImproperType)
    ):
        return potential.member_types or potential.member_classes


def get_parameters(potential):
    """Return hashable version of parameters for a potential."""
    return (
        tuple(potential.get_parameters().keys()),
        tuple(map(lambda x: x.to_value(), potential.get_parameters().values())),
    )


def filtered_potentials(potential_types, identifier):
    """Filter and return unique potentials based on pre-defined identifier function."""
    visited = defaultdict(set)

    for potential_type in potential_types:
        potential_id = identifier(potential_type)
        if potential_id not in visited[type(potential_type)]:
            visited[type(potential_type)].add(potential_id)

            yield potential_type


class PotentialFilters:
    UNIQUE_NAME_CLASS = "unique_name_class"
    UNIQUE_SORTED_NAMES = "unique_sorted_names"
    UNIQUE_EXPRESSION = "unique_expression"
    UNIQUE_PARAMETERS = "unique_parameters"
    UNIQUE_ID = "unique_id"
    REPEAT_DUPLICATES = "repeat_duplicates"

    @staticmethod
    def all():
        return set(
            f"{PotentialFilters.__name__}.{k}"
            for k, v in PotentialFilters.__dict__.items()
            if not k.startswith("__") and not callable(v)
        )


potential_identifiers = {
    PotentialFilters.UNIQUE_NAME_CLASS: get_name_or_class,
    PotentialFilters.UNIQUE_SORTED_NAMES: sort_by_types,
    PotentialFilters.UNIQUE_EXPRESSION: lambda p: str(p.expression),
    PotentialFilters.UNIQUE_PARAMETERS: get_parameters,
    PotentialFilters.UNIQUE_ID: lambda p: id(p),
    PotentialFilters.REPEAT_DUPLICATES: lambda _: str(uuid.uuid4()),
}


class TopologyPotentialView:
    """A potential view based on different filters for a topology's potentials.

    Parameters
    ----------
    iterator: typing.Iterator, required=True
        An iterator of either topology sites or connections from which to extract potentials from

    Parameters
    ----------
    filter_by: str or function, default=PotentialFilters.UNIQUE_ID
        If provided, filter the collected potentials by some unique identifier
        of a potential.

    Notes
    -----
    If `filter_by` is provided and is a custom function, the collected potentials are
    filtered on the basis of the return value of the `filter_by` function. A single potential
    is passed to the filter_by function and it should return an identifier (should be hashable)
    for that potential thus describing its uniqueness in the context of which the filter is being
    used. Some simple examples are given below.

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

    def __init__(self, iterator, filter_by=PotentialFilters.UNIQUE_ID):
        self.iterator = iterator
        self.filter_by = filter_by

    def __iter__(self):
        yield from self.yield_view()

    def index(self, item):
        for j, potential in enumerate(self.yield_view()):
            if potential is item:
                return j
        return None

    def equality_index(self, item):
        for j, potential in enumerate(self.yield_view()):
            if potential == item:
                return j
        return None

    def _collect_potentials(self):
        """Collect potentials from the iterator"""
        for item in self.iterator:
            potential = getattr(
                item, potential_attribute_map[type(item)]
            )  # Since this use is internal, KeyErrors N/A
            if potential:
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
                try:
                    identifier_func = potential_identifiers[self.filter_by]
                except KeyError:
                    raise MissingFilterError(
                        f"Potential filter {self.filter_by} is not among the built-in"
                        f"filters. Please use one of the builtin filters or define a custom "
                        f"filter callable. Builtin Filters are \n {PotentialFilters.all()}"
                    )
            else:
                identifier_func = self.filter_by

            yield from filtered_potentials(
                self._collect_potentials(), identifier=identifier_func
            )

    def __call__(self, filter_by=PotentialFilters.UNIQUE_ID):
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
