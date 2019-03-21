import warnings

from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.bond_type import BondType
from topology.exceptions import TopologyError

class Bond(Connection):
    """A 2-partner connection between sites.
    
    Partners
    --------
    connected_members: list of topology.Site
        Should be length 2
    connection_type : topology.BondType

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, connected_members=[], connection_type=None):
        connected_members = _validate_two_partners(connected_members)
        connection_type = _validate_bondtype(connection_type)

        super(Bond, self).__init__(connected_members=connected_members, 
                connection_type=connection_type)


def _validate_two_partners(connected_members):
    """Ensure 2 partners are involved in Bond"""
    if len(connected_members) != 2:
        raise TopologyError("Trying to create a Bond " 
                "with {} bond partners". format(len(connected_members)))
    
    return connected_members


def _validate_bondtype(contype):
    """Ensure connection_type is a BondType """
    if contype is None:
        warnings.warn("Non-parametrized Bond detected")
    elif not isinstance(contype, BondType):
        raise TopologyError("Supplied non-BondType {}".format(contype))
    return contype
