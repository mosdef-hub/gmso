import warnings

from gmso.core.connection import Connection
from gmso.core.dihedral_type import DihedralType
from gmso.exceptions import GMSOError


class Dihedral(Connection):
    """A 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.DihedralType.
    The connectivity of a dihedral is:

           m1–m2–m3–m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Parameters
    --------
    connection_members: list of gmso.Atom
        4 sites of a dihedral.
    connection_type : gmso.DihedralType, optional, default=None
        DihedralType of this dihedral.
    name : str, optional, default=Dihedral
        Name of the dihedral.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Additional _validate methods are presented
    """

    def __init__(self, connection_members=[], connection_type=None, name="Dihedral"):
        connection_members = _validate_four_partners(connection_members)
        connection_type = _validate_dihedraltype(connection_type)

        super(Dihedral, self).__init__(connection_members=connection_members,
                connection_type=connection_type, name=name)


def _validate_four_partners(connection_members):
    """Ensure 4 partners are involved in Dihedral"""
    assert connection_members is not None, "connection_members is not given"
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
