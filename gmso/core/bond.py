"""Module for 2-partner connections between sites."""
from typing import Callable, ClassVar, Optional, Tuple

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.bond_type import BondType

try:
    from pydantic.v1 import Field
except ImportError:
    from pydantic import Field


class Bond(Connection):
    __base_doc__ = """A 2-partner connection between sites.

    This is a subclass of the gmso.abc.Connection superclass.
    This class has strictly 2 members in its connection_members.
    The connection_type in this class corresponds to gmso.BondType.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods.
    Additional _validate methods are presented.
    """
    __members_creator__: ClassVar[Callable] = Atom.parse_obj

    connection_members_: Tuple[Atom, Atom] = Field(
        ..., description="The 2 atoms involved in the bond."
    )

    bond_type_: Optional[BondType] = Field(
        default=None, description="BondType of this bond."
    )
    restraint_: Optional[dict] = Field(
        default=None,
        description="""
        Restraint for this bond, must be a dict with the following keys:
        'b0' (unit of length), 'kb' (unit of energy/(mol * length**2)).
        Refer to https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
        for more information.
        """,
    )

    @property
    def bond_type(self):
        """Return parameters of the potential type."""
        return self.__dict__.get("bond_type_")

    @property
    def connection_type(self):
        """Return parameters of the potential type."""
        # ToDo: Deprecate this?
        return self.__dict__.get("bond_type_")

    @property
    def restraint(self):
        """Return the restraint of this bond."""
        return self.__dict__.get("restraint_")

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples.

        Returns
        -------
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        -----
        For a bond:

            i, j == j, i

        where i and j are the connection members.
        """
        return frozenset(
            [self.connection_members, tuple(reversed(self.connection_members))]
        )

    def __setattr__(self, key, value):
        """Handle attribute assignment."""
        if key == "connection_type":
            super(Bond, self).__setattr__("bond_type", value)
        else:
            super(Bond, self).__setattr__(key, value)

    class Config:
        """Pydantic configuration for Bond."""

        fields = {
            "bond_type_": "bond_type",
            "connection_members_": "connection_members",
            "restraint_": "restraint",
        }
        alias_to_fields = {
            "bond_type": "bond_type_",
            "connection_members": "connection_members_",
            "restraint": "restraint_",
        }
