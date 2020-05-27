import warnings

from gmso.core.connection import Connection
from gmso.core.angle_type import AngleType
from gmso.exceptions import GMSOError


class Angle(Connection):
    """A 3-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 3 members in its connection members.
    The connection_type in this class corresponds to gmso.AngleType.

    Parameters
    ----------
    connection_members: list of gmso.Atom
        3 sites of an angle.
    connection_type : gmso.AngleType, optional, default=None
        AngleType of this angle.
    name : str, optional, default="Angle"
        Name of the angle. 

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Additional _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None, name="Angle"):
        connection_members = _validate_three_partners(connection_members)
        connection_type = _validate_angletype(connection_type)

        super(Angle, self).__init__(connection_members=connection_members,
                connection_type=connection_type, name=name)


def _validate_three_partners(connection_members):
    """Ensure 3 partners are involved in Angle"""
    assert connection_members is not None, "connection_members is not given"
    if len(connection_members) != 3:
        raise GMSOError("Trying to create an Angle "
                "with {} connection members". format(len(connection_members)))

    return connection_members


def _validate_angletype(contype):
    """Ensure connection_type is a AngleType """
    if contype is None:
        warnings.warn("Non-parametrized Angle detected")
    elif not isinstance(contype, AngleType):
        raise GMSOError("Supplied non-AngleType {}".format(contype))
    return contype
