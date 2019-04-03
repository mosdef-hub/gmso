import warnings

from topology.core.connection import Connection
from topology.core.bond_type import BondType
from topology.exceptions import TopologyError


class Bond(Connection):
    """A 2-partner connection between sites.

    Partners
    --------
    connection_members: list of topology.Site
        Should be length 2
    connection_type : topology.BondType

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None):
        connection_members = _validate_two_partners(connection_members)
        connection_type = _validate_bondtype(connection_type)

        super(Bond, self).__init__(connection_members=connection_members,
                connection_type=connection_type)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __hash__(self):
        if self.connection_type:
            return hash(
                tuple(
                    (
                        self.connection_type,
                        tuple(self.connection_members),
                    )
                )
            )
        return hash(tuple(self.connection_members))

def _validate_two_partners(connection_members):
    """Ensure 2 partners are involved in Bond"""
    if len(connection_members) != 2:
        raise TopologyError("Trying to create a Bond "
                "with {} connection members". format(len(connection_members)))

    return connection_members


def _validate_bondtype(contype):
    """Ensure connection_type is a BondType """
    if contype is None:
        warnings.warn("Non-parametrized Bond detected")
    elif not isinstance(contype, BondType):
        raise TopologyError("Supplied non-BondType {}".format(contype))
    return contype
