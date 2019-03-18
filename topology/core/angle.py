import warnings

from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.angle_type import AngleType
from topology.exceptions import TopologyError

class Angle(Connection):
    """A 3-partner connection between sites.
    
    Partners
    --------
    bond_partners: list of topology.Site
        Should be length 3
    connection_type : topology.AngleType

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, bond_partners=[], connection_type=None):
        bond_partners = _validate_three_partners(bond_partners)
        connection_type = _validate_angletype(connection_type)

        super(Angle, self).__init__(bond_partners=bond_partners, 
                connection_type=connection_type)


def _validate_three_partners(bond_partners):
    """Ensure 2 partners are involved in Bond"""
    if len(bond_partners) != 3:
        raise TopologyError("Trying to create an Angle " 
                "with {} bond partners". format(len(bond_partners)))
    
    return bond_partners


def _validate_angletype(ctype):
    """Ensure connection_type is a AngleType """
    if ctype is None:
        warnings.warn("Non-parametrized Angle detected")
    elif not isinstance(ctype, AngleType):
        raise TopologyError("Supplied non-AngleType {}".format(ctype))
    return ctype
