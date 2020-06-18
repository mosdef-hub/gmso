import warnings

from gmso.core.connection import Connection
from gmso.core.bond_type import BondType
from gmso.exceptions import GMSOError


class Bond(Connection):
    """A 2-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 2 members in its connection_members.
    The connection_type in this class corresponds to gmso.BondType.

    Parameters
    ---------
    connection_members: list of gmso.Atom
        2 sites of a bond.
    connection_type : gmso.BondType, optional, default=None
        BondType of this bond.
    name : str, optional, default="Bond"
        Name of the bond.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods.
    Additional _validate methods are presented.
    """

    def __init__(self, connection_members=None, connection_type=None, name="Bond"):
        connection_members = _validate_two_partners(connection_members)
        connection_type = _validate_bondtype(connection_type)

        super(Bond, self).__init__(connection_members=connection_members,
                connection_type=connection_type, name=name)


def _validate_two_partners(connection_members):
    """Ensure 2 partners are involved in Bond"""
    assert connection_members is not None, "connection_members is not given"
    if len(connection_members) != 2:
        raise GMSOError("Trying to create a Bond "
                "with {} connection members". format(len(connection_members)))
    return connection_members


def _validate_bondtype(contype):
    """Ensure connection_type is a BondType """
    if contype is None:
        warnings.warn("Non-parametrized Bond detected")
    elif not isinstance(contype, BondType):
        raise GMSOError("Supplied non-BondType {}".format(contype))
    return contype
