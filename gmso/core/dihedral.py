import warnings

from gmso.core.connection import Connection
from gmso.core.dihedral_type import DihedralType
from gmso.exceptions import GMSOError


class Dihedral(Connection):
    """A 4-partner connection between sites.

    Partners
    --------
    connection_members: list of gmso.Site
        Should be length 4
    connection_type : gmso.DihedralType
    name : name of the dihedral
        inherits the name attribute from Connection

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Addiitonal _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None, name="Dihedral"):
        connection_members = _validate_four_partners(connection_members)
        connection_type = _validate_dihedraltype(connection_type)

        super(Dihedral, self).__init__(connection_members=connection_members,
                connection_type=connection_type, name=name)


def _validate_four_partners(connection_members):
    """Ensure 4 partners are involved in Dihedral"""
    if len(connection_members) != 4:
        raise GMSOError("Trying to create an Dihedral "
                "with {} connection members". format(len(connection_members)))

    return connection_members


def _validate_dihedraltype(contype):
    """Ensure connection_type is a DihedralType """
    if contype is None:
        warnings.warn("Non-parametrized Dihedral detected")
    elif not isinstance(contype, DihedralType):
        raise GMSOError("Supplied non-DihedralType {}".format(contype))
    return contype
