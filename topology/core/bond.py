from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.bond_type import BondType
from topology.exceptions import TopologyError

class Bond(Connection):
    """A 2-partner connection between sites.
    
    Partners
    --------
    bond_partners: list of topology.Site
        Should be length 2
    connection_type : topology.BondType

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, bond_partners=[], connection_type=None):
        bond_partners = _validate_two_partners(bond_partners)
        connection_type = _validate_bondtype(connection_type)

        super(Bond, self).__init__(bond_partners=bond_partners, 
                connection_type=connection_type)


def _validate_two_partners(bond_partners):
    """Ensure 2 partners are involved in Bond"""
    if len(bond_partners) != 2:
        raise TopologyError("Trying to create a Bond " 
                "with {} bond partners". format(len(bond_partners)))
    
    return bond_partners


def _validate_bondtype(ctype):
    """Ensure connection_type is a BondType """
    if not isinstance(ctype, BondType):
        raise TopologyError("Supplied non-BondType {}".format(c_type))
    return ctype
