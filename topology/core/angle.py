import warnings

from topology.core.connection import Connection
from topology.core.angle_type import AngleType
from topology.exceptions import TopologyError

class Angle(Connection):
    """A 3-partner connection between sites.
    
    Partners
    --------
    connection_members: list of topology.Site
        Should be length 3
    connection_type : topology.AngleType

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None):
        connection_members = _validate_three_partners(connection_members)
        connection_type = _validate_angletype(connection_type)

        super(Angle, self).__init__(connection_members=connection_members, 
                connection_type=connection_type)


def _validate_three_partners(connection_members):
    """Ensure 3 partners are involved in Angle"""
    if len(connection_members) != 3:
        raise TopologyError("Trying to create an Angle " 
                "with {} connection members". format(len(connection_members)))
    
    return connection_members


def _validate_angletype(contype):
    """Ensure connection_type is a AngleType """
    if contype is None:
        warnings.warn("Non-parametrized Angle detected")
    elif not isinstance(contype, AngleType):
        raise TopologyError("Supplied non-AngleType {}".format(contype))
    return contype
