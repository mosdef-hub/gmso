from typing import Optional, Tuple

from pydantic import Field

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.dihedral_type import DihedralType


class Dihedral(Connection):
    __base_doc__ = """A 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.DihedralType.
    The connectivity of a dihedral is:
        m1–m2–m3–m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods

    Additional _validate methods are presented
    """
    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ..., description="The 4 atoms involved in the dihedral."
    )

    dihedral_type_: Optional[DihedralType] = Field(
        default=None, description="DihedralType of this dihedral."
    )

    @property
    def dihedral_type(self):
        return self.__dict__.get("dihedral_type_")

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get("dihedral_type_")

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples

        Returns
        _______
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        _____
        For a dihedral:

            i, j, k, l == l, k, j, i

        where i, j, k, and l are the connection members.
        """
        return frozenset(
            [self.connection_members, tuple(reversed(self.connection_members))]
        )

    def _equivalent_members_hash(self):
        """Returns a unique hash representing the connection
        Returns
        _______
        int
            A unique hash to represent the connection members
        Notes
        _____
        For a dihedral:
            i, j, k, l == l, k, j, i
        where i, j, k, and l are the connection members.
        Here i and j are interchangeable, j and k are interchangeable,
        and k and l are interchangeble, as long as each are adjacent to
        one another.
        """

        return hash(
            frozenset(
                [
                    frozenset(
                        [self.connection_members[0], self.connection_members[1]]
                    ),
                    frozenset(
                        [self.connection_members[1], self.connection_members[2]]
                    ),
                    frozenset(
                        [self.connection_members[2], self.connection_members[3]]
                    ),
                ]
            )
        )

    def __setattr__(self, key, value):
        if key == "connection_type":
            super(Dihedral, self).__setattr__("dihedral_type", value)
        else:
            super(Dihedral, self).__setattr__(key, value)

    class Config:
        fields = {
            "dihedral_type_": "dihedral_type",
            "connection_members_": "connection_members",
        }
        alias_to_fields = {
            "dihedral_type": "dihedral_type_",
            "connection_members": "connection_members_",
        }
