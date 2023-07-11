from typing import Callable, ClassVar, Optional, Tuple

from pydantic import ConfigDict, Field

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
    __members_creator__: ClassVar[Callable] = Atom.parse_obj

    connection_members: Tuple[Atom, Atom, Atom, Atom] = Field(
        ..., description="The 4 atoms involved in the dihedral."
    )

    dihedral_type: Optional[DihedralType] = Field(
        default=None, description="DihedralType of this dihedral."
    )

    restraint: Optional[dict] = Field(
        default=None,
        description="""
        Restraint for this dihedral, must be a dict with the following keys:
        'k' (unit of energy/(mol * angle**2)), 'phi_eq' (unit of angle), 'delta_phi' (unit of angle).
        Refer to https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
        for more information.
        """,
    )

    @property
    def dihedral_type(self):
        return self.__dict__.get("dihedral_type")

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get("dihedral_type")

    @property
    def restraint(self):
        """Return the restraint of this dihedral."""
        return self.__dict__.get("restraint")

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

    def __setattr__(self, key, value):
        if key == "connection_type":
            super(Dihedral, self).__setattr__("dihedral_type", value)
        else:
            super(Dihedral, self).__setattr__(key, value)
