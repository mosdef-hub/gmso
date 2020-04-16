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
    connection_members: list of gmso.Site
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

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples

        Returns
        _______
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        _____
        For an dihedral:
            i, j, k, l == l, k, j, i
        where i, j, k, and l are the connection members.
        """
        return frozenset([
                tuple(self.connection_members),
                tuple(reversed(self.connection_members))
                ])

    def _equivalent_member_hash(self):
        """Get a unique hash representing the connection

        Returns
        _______
        int
            A unique hash to represent the connection members

        Notes
        _____
        For an dihedral:
            i, j, k, l == l, k, j, i
        where i, j, k, and l are the connection members.
        Here i and j are interchangeable, j and k are interchangeable,
        and k and l are interchangeble, as long as each are adjacent to
        one another.
        """

        return hash(tuple([frozenset([
            frozenset([self.connection_members[0],
                       self.connection_members[1]]),
            frozenset([self.connection_members[1],
                       self.connection_members[2]]),
            frozenset([self.connection_members[2],
                       self.connection_members[3]])
            ])]))

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
