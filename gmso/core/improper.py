import warnings

from gmso.core.connection import Connection
from gmso.core.improper_type import ImproperType
from gmso.exceptions import GMSOError


class Improper(Connection):
    """A 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.ImproperType
    The connectivity of an improper is:

                   m2
                   |
                   m1
                  / \
                 m3  m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Parameters
    --------
    connection_members: list of gmso.Atom
        4 sites of a improper. Central site first, then three
        sites connected to the central site.
    connection_type : gmso.ImproperType, optional, default=None
        ImproperType of this improper.
    name : str, optional, default=Improper
        Name of the improper.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Additional _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None, name="Improper"):
        connection_members = _validate_four_partners(connection_members)
        connection_type = _validate_impropertype(connection_type)

        super(Improper, self).__init__(connection_members=connection_members,
                connection_type=connection_type, name=name)


def _validate_four_partners(connection_members):
    """Ensure 4 partners are involved in Improper"""
    assert connection_members is not None, "connection_members is not given"
    if len(connection_members) != 4:
        raise GMSOError("Trying to create an Improper "
                "with {} connection members". format(len(connection_members)))
    return connection_members


def _validate_impropertype(contype):
    """Ensure connection_type is a ImproperType """
    if contype is None:
        warnings.warn("Non-parametrized Improper detected")
    elif not isinstance(contype, ImproperType):
        raise GMSOError("Supplied non-ImproperType {}".format(contype))
    return contype
