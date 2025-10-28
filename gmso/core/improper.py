"""Support for improper style connections (4-member connection)."""

from typing import Callable, ClassVar, Optional, Tuple

from pydantic import ConfigDict, Field, model_validator

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.improper_type import ImproperType


class Improper(Connection):
    """A 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.ImproperType
    The connectivity of an improper is:
             m2
        m3 - m1 - m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Notes
    -----
    Inherits some methods from Connection:
    __eq__, __repr__, _validate methods

    Additional _validate methods are presented.
    """

    __members_creator__: ClassVar[Callable] = Atom.model_validate

    connectivity: ClassVar[Tuple[Tuple[int]]] = ((0, 1), (0, 2), (0, 3))

    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ...,
        description="The 4 atoms of this improper. Central atom first, "
        "then the three atoms connected to the central site.",
        alias="connection_members",
    )

    bonds_: Tuple[Bond, Bond, Bond] = Field(
        default=None,
        description="""
        List of connection bonds.
        Ordered to align with connection_members, such that self.bonds_[0] is
        the bond between (self.connection_members[0], self.connection_members[1]).
        """,
        alias="bonds",
    )

    improper_type_: Optional[ImproperType] = Field(
        default=None,
        description="ImproperType of this improper.",
        alias="improper_type",
    )

    bond_orders_: Optional[Tuple[str, str]] = Field(
        default=None,
        description="""
        List of connection members bond orders.
        """,
        alias="bond_orders",
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **Connection.model_config["alias_to_fields"],
            **{
                "improper_type": "improper_type_",
                "bond_orders": "bond_orders_",
            },
        )
    )

    @property
    def improper_type(self):
        """Return Potential object for this connection if it exists."""
        return self.__dict__.get("improper_type_")

    @property
    def connection_type(self):
        """Return Potential object for this connection if it exists."""
        # ToDo: Deprecate this?
        return self.__dict__.get("improper_type_")

    @property
    def bonds(self):
        """Return the bonds that makeup this Improper.

        Connectivity is ((0,1), (0,2), (0,3))
        """
        return self.__dict__.get("bonds_")

    @property
    def bonds_orders(self):
        """Return the bond_order strings of this improper."""
        return [b.bond_order for b in self.bonds]

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples.

        Returns
        -------
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        -----
        For an improper:

            i, j, k, l == i, k, j, l

        where i, j, k, and l are the connection members.
        """
        equiv_members = [
            self.connection_members[0],
            self.connection_members[2],
            self.connection_members[1],
            self.connection_members[3],
        ]

        return frozenset([self.connection_members, tuple(equiv_members)])

    def __setattr__(self, key, value):
        """Set attribute override to support connection_type key."""
        if key == "connection_type":
            super(Improper, self).__setattr__("improper_type_", value)
        else:
            super(Improper, self).__setattr__(key, value)

    @model_validator(mode="before")
    @classmethod
    def set_dependent_value_default(cls, data):
        """Automatically set bonds for this improper if connection_members is defined."""
        if "bonds" not in data and "connection_members" in data:
            atoms = data["connection_members"]
            data["bonds"] = (
                Bond(connection_members=(atoms[0], atoms[1])),
                Bond(connection_members=(atoms[0], atoms[2])),
                Bond(connection_members=(atoms[0], atoms[3])),
            )
        return data
